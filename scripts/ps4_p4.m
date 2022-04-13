clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get data
load('../data/radarmeasdata_missle.mat')

% Constants / given data
sigma_rhoa = 10; %m
sigma_rhob = 30; %m

global l_A l_B g
l_A = 4.1e5;
l_B = 4.4e5;
g = 9.80665;

% Initialize z_prime, R matrices
global n_z
n_z = 2; % number of measurements
zhist = zeros(n_z*length(thist),1);
R = zeros(n_z*length(thist),n_z*length(thist));

% Build z, R matrices and find cholesky factorization:

for ii=1:length(thist)
    
    zhist(n_z*(ii-1)+1) = rhoahist(ii);
    zhist(n_z*(ii-1)+2) = rhobhist(ii);
    
    R(n_z*(ii-1)+1,n_z*(ii-1)+1) = sigma_rhoa;
    R(n_z*(ii-1)+2,n_z*(ii-1)+2) = sigma_rhob;
    
end

R_a = chol(R);
%z = z = R_a'*zhist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute initial guess, x_g
delta_l = l_B - l_A;

% Get first and last measurements
rhoa_i = rhoahist(1);
rhob_i = rhobhist(1);
rhoa_f = rhoahist(length(rhoahist))
rhob_f = rhobhist(length(rhobhist))
delta_t = thist(length(thist)) - thist(1);

%Compute initial azimuth angle and position relative to site A
theta_i = acos((rhoa_i^2 + delta_l^2 - rhob_i^2)/(2*rhoa_i*delta_l));
% Note - round this to pi
theta_i = pi;
y_1_i = l_A + (rhoa_i^2 + delta_l^2 - rhob_i^2)/(2*delta_l);
y_2_i = rhoa_i*sin(theta_i);

%Compute final azimuth angle and position relative to site A
theta_f = acos((rhoa_f^2 + delta_l^2 - rhob_f^2)/(2*rhoa_f*delta_l));
y_1_f = l_A + (rhoa_f^2 + delta_l^2 - rhob_f^2)/(2*delta_l);
y_2_f = rhoa_f*sin(theta_f);

% Compute initial velocity estimate from 
v_i = (y_1_f - y_1_i)/delta_t;

x_g = [y_1_i;
       v_i;
       y_2_i; %real(y_1_f - y_1_i)/delta_t;
       ((y_2_f - y_2_i) + .5*g*delta_t^2)/delta_t]
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Solve for x to minimize J
   
   
% Compute cost at x_g
[x_hat, P_xx, J_pct_chg] = ps4p4NLLS(x_g, zhist, thist, R_a)
num_calc = 1;
% Continue to iterate until change in J is small (less than 1%)
while J_pct_chg > .01
    [x_hat, P_xx, J_pct_chg] = ps4p4NLLS(x_hat, zhist, thist, R_a);
    num_calc = num_calc+1;
end

x_hat
P_xx
J_pct_chg
num_calc

% Plot, compare initial and final output trajectory
figure(1)
t = 0:1:120;
y_1_initial = x_g(1) + x_g(2).*t;
y_2_initial = x_g(3) + x_g(4).*t - .5*g.*t.*t;
y_1_final = x_hat(1) + x_hat(2).*t;
y_2_final = x_hat(3) + x_hat(4).*t - .5*g.*t.*t;
plot(y_1_initial,y_2_initial,'b--',y_1_final,y_2_final,'b-')

% Want z and h matrices, initial and final (z should be the same)
[z_0, h_0, H_0] = ps4p4matrices(zhist, x_g, thist, R_a);
[z_f, h_f, H_f] = ps4p4matrices(zhist, x_hat, thist, R_a);

% Split the normalized z and h vectors. A are odd indices, B are even
% indices
num_A = 1;
num_B = 1;
for ii=1:length(z_0)
    
    if mod(ii,2)==1 %If this is an odd index, add to A vector
        z_A(num_A) = z_0(ii);
        h_0_A(num_A) = h_0(ii);
        h_f_A(num_A) = h_f(ii);
        
        num_A = num_A+1;
    end 
    if mod(ii,2)==0 %If this is an even index, 
        z_B(num_B) = z_0(ii);
        h_0_B(num_B) = h_0(ii);
        h_f_B(num_B) = h_f(ii);
        
        num_B = num_B+1;
    end

end
figure(2)
plot(z_A,z_B,'bo',h_0_A,h_0_B,'b--',h_f_A,h_f_B,'b-')

% Gauss-Newton nonlinear least squares solver
function [x_hat, P_xx, J_pct_chg] = ps4p4NLLS(x_g, zhist, thist, R_a)

    % Generate normalized matrices
    [z, h, H] = ps4p4matrices(zhist, x_g, thist, R_a);
    
    % Calculate d_x_hat
    d_x_hat = inv(H'*H)*H'*(z - h);
    
    % Evaluate cost function J at alpha = 1
    alpha=1;
    [J_0] = cost(z, h)
    [z, h_alpha, H_alpha] = ps4p4matrices(zhist, (x_g + alpha*d_x_hat), thist, R_a); % get h_x matrix at alpha = 1
    [J_alpha] = cost(z, h_alpha)
    
    % Iterate to make sure alpha value reduces cost, J
    while (J_alpha > J_0)
        alpha = alpha/2;
        [z, h_alpha, H_alpha] = ps4p4matrices(zhist, (x_g + alpha*d_x_hat), thist, R_a); % get h_x matrix at alpha = 1
        [J_alpha] = cost(z, h_alpha)
    end %while
    
    J_pct_chg = (abs(J_alpha - J_0)/J_0);
    
    % Compute x_hat and associated P_xx
    x_hat = x_g + alpha*d_x_hat;
    P_xx = inv(H_alpha'*H_alpha);

end


function [z, h_x, H_x] = ps4p4matrices(z_prime, x, t, Ra)

global l_A l_B g n_z

% Normalize z
z = Ra'*z_prime;

% For x_g, construct h(x) vector and H(x) matrix
h_prime_x = zeros(n_z*length(t),1);
H_prime_x = zeros(n_z*length(t),length(x));

for ii = 1:length(t)
    delta_t = t(ii) - t(1);
    delta_y_1_A = l_A - x(1) -x(2)*delta_t;
    delta_y_1_B = l_B - x(1) -x(2)*delta_t;
    delta_y_2 = x(3) + delta_t*x(4) - .5*g*delta_t^2;
    
    h_prime_x(n_z*(ii-1)+1) = sqrt(delta_y_1_A^2 + delta_y_2^2);
    h_prime_x(n_z*(ii-1)+2) = sqrt(delta_y_1_B^2 + delta_y_2^2);
    
    H_prime_x(n_z*(ii-1)+1,1) = -delta_y_1_A/h_prime_x(n_z*(ii-1)+1);
    H_prime_x(n_z*(ii-1)+1,2) = -delta_t*delta_y_1_A/h_prime_x(n_z*(ii-1)+1);
    H_prime_x(n_z*(ii-1)+1,3) = delta_y_2/h_prime_x(n_z*(ii-1)+1);
    H_prime_x(n_z*(ii-1)+1,4) = delta_t*delta_y_2/h_prime_x(n_z*(ii-1)+1);
    
    H_prime_x(n_z*(ii-1)+2,1) = -delta_y_1_B/h_prime_x(n_z*(ii-1)+2);
    H_prime_x(n_z*(ii-1)+2,2) = -delta_t*delta_y_1_B/h_prime_x(n_z*(ii-1)+2);
    H_prime_x(n_z*(ii-1)+2,3) = delta_y_2/h_prime_x(n_z*(ii-1)+2);
    H_prime_x(n_z*(ii-1)+2,4) = delta_t*delta_y_2/h_prime_x(n_z*(ii-1)+2);
    
end
% Now normalize h_prime_x = h_x and H_prime_x = H_x
h_x = Ra'*h_prime_x; 
H_x = Ra'*H_prime_x; 

end %matrices function

function J = cost(z, h_x)
% Compute cost and return it
J = norm((z - h_x).^2);
end %cost
   
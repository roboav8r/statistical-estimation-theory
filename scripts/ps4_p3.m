% Clear workspace
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Accept inputs and parameters
% Time series
thist =[0; 0.1000; 0.2000; 0.3000; 0.4000; 0.5000;
0.6000; 0.7000; 0.8000; 0.9000; 1.0000];

% Measurement seriess
zhist = [7.7969; 1.4177; -3.0970; -7.6810; -9.8749; -6.1828;
-0.8212; 4.5074; 8.2259; 9.5369; 6.2827];


% Build noise covariance matrix R and find Cholesky factorization R_a
R = zeros(length(zhist),length(zhist));
for ii = 1:length(zhist)
    for jj = 1:length(zhist)
        if ii==jj
            R(ii,jj) = 1;
        elseif abs(ii-jj)==1
            R(ii,jj) = 0.5;
        end
    end
end

R_a = chol(R);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Solve the system
% Give initial guess of x, x_0:
x_0 = [10, 2*pi, 0]'

% Propagate forward one step. Compute new estimate x_hat, covariance P_xx,
% and the percentage change in cost function.
[x_hat, P_xx, J_pct_chg] = ps4p3NLLS(x_0, zhist, thist, R_a, .01);

% Continue propagating until the change in cost function is below a
% tolerance of 1%
while J_pct_chg > .01
    [x_hat, P_xx, J_pct_chg] = ps4p3NLLS(x_hat, zhist, thist, R_a, .01);
end

% Output the final estimate to the terminal
x_hat
P_xx

% Plot output
t = 0:.01:1;
y_0 = x_0(1)*cos(x_0(2)*t + x_0(3));
y_f = x_hat(1)*cos(x_hat(2)*t + x_hat(3));
plot(thist,zhist,'bo',t,y_0,'b--',t,y_f,'b-')
legend('measurement z_{t}','Initial signal estimate','Final signal estimate')
title('Gauss-Newton estimate of sinusoidal signal parameters')


% Gauss-Newton 
function [x_hat, P_xx, J_pct_chg] = ps4p3NLLS(x_g, zhist, thist, R_a, tolerance)

    % Normalize z
    z = R_a'*zhist;
    
    % For x_g, construct h(x) vector
    h_prime_x = zeros(length(zhist),1);
    for ii = 1:length(zhist)
        h_prime_x(ii) = x_g(1)*cos(x_g(2)*thist(ii) + x_g(3));
    end
    % Now normalize h_prime_x = h_x
    h_x = R_a'*h_prime_x; 
    
    % Generate H matrix (linearized dz/dx matrix)
    H = ps4p3H(x_g, z, thist, R_a);
    
    % Calculate d_x_hat / "step" toward new solution
    d_x_hat = inv(H'*H)*H'*(z - h_x);
    
    % Evaluate cost function J at alpha = 1 (full step toward new solution) 
    % to determine if this results in a least-squares error reduction
    
    % Compute current cost
    [J_0,h_x_0] = cost(z, x_g, thist, R_a);
    % Estimate new cost at one full step toward new solution
    alpha=1;
    [J_alpha, h_x_alpha] = cost(z, (x_g + alpha*d_x_hat), thist, R_a);
    
    % If the step does not result in 
    while (J_alpha - J_0)/ J_0 > tolerance
        alpha = alpha/2;
        [J_alpha, h_x_alpha] = cost(z, (x_g + alpha*d_x_hat), thist, R_a);
    end %while
    J_pct_chg = (abs(J_alpha - J_0)/J_0);
    
    % Compute new x_hat (parameter estimate)
    x_hat = x_g + alpha*d_x_hat;
    
    % Compute covariance P_xx for new estimate x_hat
    H_new = ps4p3H(x_hat, z, thist, R_a);
    P_xx = inv(H_new'*H_new);

end % ps4p3NLLS function 

% Compute H(x), linearized dz/dx matrix
function H_x = ps4p3H(x, z, t, R_a)

    H_prime_x = zeros(length(z),length(x));
    for ii=1:length(z)
        H_prime_x(ii,1) = cos(x(2)*t(ii) + x(3));
        H_prime_x(ii,2) = -t(ii)*x(1)*sin(x(2)*t(ii) + x(3));
        H_prime_x(ii,3) = -x(1)*sin(x(2)*t(ii) + x(3));
    end
    
    % Normalize
    H_x = R_a'*H_prime_x;

end %H function


% The Least-Squares cost function
function [J, h_x] = cost(z, x, t, Ra)
% For x_g, construct h(x) vector
h_prime_x = zeros(length(z),1);
for ii = 1:length(z)
    h_prime_x(ii) = x(1)*cos(x(2)*t(ii) + x(3));
end
% Now normalize h_prime_x = h_x
h_x = Ra'*h_prime_x; 

% Compute cost and return it
J = norm((z - h_x).^2);

end %cost function
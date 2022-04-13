clear all
clc

format long

% parameters of the problem
nx=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part A, i
r_bar = 76;
theta_bar = -3*pi/180;
sigma_r_2 = 1^2;
sigma_theta_2 = (pi/180)^2;
[x_lin_i, P_lin_i] = lin_transform(r_bar, theta_bar, r_bar, theta_bar, sigma_r_2, sigma_theta_2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part A, ii
sigma_theta_2 = (15*pi/180)^2;
[x_lin_ii, P_lin_ii] = lin_transform(r_bar, theta_bar, r_bar, theta_bar, sigma_r_2, sigma_theta_2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part B, i
global alpha beta kappa
alpha = 1e-3;
beta=2;
kappa=0;

r_bar = 76;
theta_bar = -3*pi/180;
sigma_r_2 = 1^2;
sigma_theta_2 = (pi/180)^2;

[x_unscent_i, P_unscent_i] = unscent_transform(r_bar, theta_bar, sigma_r_2, sigma_theta_2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part B, ii
sigma_theta_2 = (15*pi/180)^2;

[x_unscent_ii, P_unscent_ii] = unscent_transform(r_bar, theta_bar, sigma_r_2, sigma_theta_2)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part C, i
n_rand = 1000000;
z_rand = zeros(nx,n_rand);

z_bar = [76; -3*pi/180];
P_zz=[1 0;
    0 (pi/180)^2];

% Generate matrix of random z values
for ii=1:n_rand
    z_rand_i(:,ii) = make_rand_x(z_bar,P_zz);
end %for

% Compute x based on these z values
x_rand(1,:) = z_rand_i(1,:).*cos(z_rand_i(2,:));
x_rand(2,:) = z_rand_i(1,:).*sin(z_rand_i(2,:));

% Compute mean and covariance of actual sample data
x_actual_i = mean(x_rand')'
P_actual_i = cov(x_rand')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part C, ii
z_bar = [76; -3*pi/180];
P_zz=[1 0;
    0 (15*pi/180)^2];

% Generate matrix of random z values
for ii=1:n_rand
    z_rand_ii(:,ii) = make_rand_x(z_bar,P_zz);
end %for

% Compute x based on these z values
x_rand_ii(1,:) = z_rand_ii(1,:).*cos(z_rand_ii(2,:));
x_rand_ii(2,:) = z_rand_ii(1,:).*sin(z_rand_ii(2,:));

% Compute mean and covariance of actual sample data
x_actual_ii = mean(x_rand_ii')'
P_actual_ii = cov(x_rand_ii')


% Compute unscented mean and covariance

function [x_lin, P_lin] = lin_transform(r, theta, r_bar, theta_bar, sigma_r_2, sigma_theta_2)

    % compute jacobian about z_bar
    A = [cos(theta_bar) -r_bar*sin(theta_bar);
        sin(theta_bar) r_bar*cos(theta_bar)];

    % compute delta_z
    delta_z = [r - r_bar;
        theta - theta_bar];
    
    % Compute expected value per derivation
    x_lin = [r_bar*cos(theta_bar); 
        r_bar*sin(theta_bar)] + A*delta_z; 
    
    % Compute covariance estimate per derivation
    P_zz=[sigma_r_2 0;
        0 sigma_theta_2];
    P_lin = A*P_zz*A';

end %lin_transform


function [x_unscent, P_unscent] = unscent_transform(r_bar, theta_bar, sigma_r_2, sigma_theta_2)
    global alpha beta kappa

    % Generate gaussian covariance matrix; decompose into S matrix
    P_zz=[sigma_r_2 0;
        0 sigma_theta_2];
    S = chol(P_zz);
    
    % Compute sigma points
    nx=2;
    lambda = alpha^2*(nx + kappa) -nx;
    
    sigma_pts = zeros(nx, 2*nx+1);
    sigma_pts(:,1) = [r_bar, theta_bar]';
    
    for ii=1:nx
        sigma_pts(:,ii+1) = [r_bar, theta_bar]' + sqrt(nx + lambda)*S(:,ii);
    end
    
    for ii=(nx+1):2*nx
        sigma_pts(:,ii+1) = [r_bar, theta_bar]' - sqrt(nx + lambda)*S(:,ii-nx);
    end
    
    % Pass sigma points through nonlinear transform
    x_xformed(1,:) = sigma_pts(1,:).*cos(sigma_pts(2,:));
    x_xformed(2,:) = sigma_pts(1,:).*sin(sigma_pts(2,:));
    
    % Generate weight vector for mean
    W_m = ones(2*nx+1,1)./(2*(nx+lambda));
    W_m(1) = lambda/(nx+lambda); %W_m_0
    
    % Generate weight vector for covariance
    W_c = W_m;
    W_c(1) = W_m(1) + 1 - alpha^2 + beta; %W_c_0
    
    % apply weighting, then sum
    x_weighted = x_xformed.*W_m'; %Create matrix where rows=x, col = weighted sigma points
    x_unscent=sum(x_weighted')'; %then sum them
    
    P_xx = zeros(nx);
    for ii=1:(2*nx+1)
        delta_x = x_xformed(:,ii) - x_unscent;
        P_xx = P_xx + W_c(ii)*(delta_x*delta_x');
    end
    
    P_unscent = P_xx;
        

end
format long

clear all
clc

% Load data
global Fk Gammak Hk P0 Qk Rk thist xhat0 zhist
run('../data/kf_example03a')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KALMAN FILTER, DATA SET A
% Initial propagation step
[x_bar, P_bar] = kf_prop(xhat0, P0);

% Initial update step
[x_hat(:,:,1),P(:,:,1)] = kf_update(x_bar, P_bar, zhist(1));

for k = 2:length(zhist) 
    [x_bar, P_bar] = kf_prop(x_hat(:,:,k-1), P(:,:,k-1));
    [x_hat(:,:,k),P(:,:,k)] = kf_update(x_bar, P_bar, zhist(k));
end %for loop

disp('Kalman filter results - data set a')
% Print output
x_hat(:,:,length(zhist))
P(:,:,length(zhist))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SRIF, DATA SET A
disp('SRIF - data set a')
% Generate initial R_vv, R_xx matrices
R_xx_0 = [inv(chol(P0))]';
R_vv_0 = [inv(chol(Qk))]';

z_x_0 = R_xx_0*xhat0;

% Initial propagation step
[z_bar_x, R_bar_xx] = srif_prop(R_vv_0, R_xx_0,z_x_0);

% Initial update step
[z_x(:,:,1), R_xx(:,:,1)] = srif_update(z_bar_x, R_bar_xx, zhist(1));

% Remaining updates
for k = 2:length(zhist)
    [z_bar_x, R_bar_xx] = srif_prop(R_vv_0, R_xx(:,:,k-1),z_x(:,:,k-1));
    [z_x(:,:,k), R_xx(:,:,k)] = srif_update(z_bar_x, R_bar_xx, zhist(k));
end

% Print output
inv(R_xx(:,:,length(zhist)))*z_x(:,:,length(zhist))
inv(R_xx(:,:,length(zhist)))*inv(R_xx(:,:,length(zhist)))'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load data set B
% Load data
clear all
global Fk Gammak Hk P0 Qk Rk thist xhat0 zhist
run('../data/kf_example03b')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KALMAN FILTER, DATA SET B
% Initial propagation step
[x_bar, P_bar] = kf_prop(xhat0, P0);

% Initial update step
[x_hat(:,:,1),P(:,:,1)] = kf_update(x_bar, P_bar, zhist(1));

for k = 2:length(zhist) 
    [x_bar, P_bar] = kf_prop(x_hat(:,:,k-1), P(:,:,k-1));
    [x_hat(:,:,k),P(:,:,k)] = kf_update(x_bar, P_bar, zhist(k));
end %for loop

disp('Kalman filter results - data set b')
% Print output
x_hat(:,:,length(zhist))
P(:,:,length(zhist))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SRIF, DATA SET A
disp('SRIF - data set b')
% Generate initial R_vv, R_xx matrices
R_xx_0 = [inv(chol(P0))]';
R_vv_0 = [inv(chol(Qk))]';

z_x_0 = R_xx_0*xhat0;

% Initial propagation step
[z_bar_x, R_bar_xx] = srif_prop(R_vv_0, R_xx_0,z_x_0);

% Initial update step
[z_x(:,:,1), R_xx(:,:,1)] = srif_update(z_bar_x, R_bar_xx, zhist(1));

% Remaining updates
for k = 2:length(zhist)
    [z_bar_x, R_bar_xx] = srif_prop(R_vv_0, R_xx(:,:,k-1),z_x(:,:,k-1));
    [z_x(:,:,k), R_xx(:,:,k)] = srif_update(z_bar_x, R_bar_xx, zhist(k));
end

% Print output
inv(R_xx(:,:,length(zhist)))*z_x(:,:,length(zhist))
inv(R_xx(:,:,length(zhist)))*inv(R_xx(:,:,length(zhist)))'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS
function [x_bar, P_bar] = kf_prop(xhat, P)
    global Fk Gammak Qk
    x_bar = Fk*xhat;
    P_bar = Fk*P*Fk' + Gammak*Qk*Gammak';
end % prop function

function [x_hat,P] = kf_update(x_bar, P_bar, z)
    global Hk Rk
    
    nu = z - Hk*x_bar;
    S = Hk*P_bar*Hk' + Rk;
    W = P_bar*Hk'*inv(S);
    epsilon_nu = nu*inv(S)*nu;
    
    x_hat = x_bar + W*nu;
    P = P_bar - W*S*W';
end %kf

function [z_bar_x, R_bar_xx] = srif_prop(R_vv, R_xx, z_x)
    global Fk Gammak
    n_z = 1;
    n_v = length(R_vv);
    n_x = length(R_xx);
    
    % QR factorize
    [Q_a, R_a] = qr([[R_vv zeros(n_v,n_x)];
        [-R_xx*inv(Fk)*Gammak R_xx*inv(Fk)]]);
    
    % Perform orthonormal transformation
    z_matrix = Q_a'*[zeros(n_v,n_z); z_x];
    
    % Extract & Return R_bar_xx and z_bar_x
    z_bar_x = z_matrix((n_v+1:n_v+n_x),:);
    R_bar_xx = R_a((n_v+1:n_v+n_x),(n_v+1:n_v+n_x));

end % prop function

function [z_x, R_xx] = srif_update(z_bar_x, R_bar_xx, z)
    global Hk Rk
    n_z = 1;
    %n_v = length(R_vv);
    n_x = length(R_bar_xx);
    
    % Compute cholesky factorization of R
    R_a = chol(Rk);
    
    % Transform z and H
    z_a = inv(R_a)'*z;
    H_a = inv(R_a)'*Hk;
    
    % Compute next QR factorization
    [Q_b, R_b] = qr([R_bar_xx; H_a]);
    
    % transform z matrix
    z_x_r = Q_b'*[z_bar_x; z_a];
    
    % Extract z_x and R_xx
    z_x = z_x_r(1:n_x,:);
    R_xx = R_b(1:n_x,:);
    
    
end %kf
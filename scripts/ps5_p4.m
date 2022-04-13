clear all
clc

% Load data
global Fk Hk Qk Rk Gammak
run('../data/kf_example02a')

% Initialize loop and variables
k = 0;
x_hat = zeros(2,1,length(zhist));
P = zeros(2,2,length(zhist));

% Initial propagation step:
x_bar = Fk*xhat0
P_bar = Fk*P0*Fk' + Gammak*Qk*Gammak'

% Initial correction step:
[x_hat(:,:,1),P(:,:,1)] = kf(x_bar, P_bar, zhist(1));

for k = 1:length(zhist) 
    x_bar(:,:,k) = Fk*x_hat(:,:,k);
    P_bar(:,:,k) = Fk*P(:,:,k)*Fk' + Gammak*Qk*Gammak';
    
    % Write to output vector
    xhat1(k) = x_hat(1,1,k);
    xhat2(k) = x_hat(2,1,k);
    P1(k) = sqrt(P(1,1,k));
    P2(k) = sqrt(P(1,1,k));
    
    [x_hat(:,:,k+1),P(:,:,k+1)] = kf(x_bar(:,:,k), P_bar(:,:,k), zhist(k));
end

subplot(2,1,1)
plot(thist,xhat1, thist,P1)

subplot(2,1,2)
plot(thist,xhat2, thist,P2)

x_hat(:,:,50)
P(:,:,50)

% Problem 4
sys = ss(Fk, Gammak, Hk, 0, -1);
[km, L, Pss, Wss] = kalman(sys, Qk, Rk)

% error transition matrix
ETM = (eye(2) - Wss*Hk)*Fk
disp('Eigenvalues of the error transition matrix:')
eig(ETM)
function [x_hat,P] = kf(x_bar, P_bar, z)
    global Hk Rk
    
    nu = z - Hk*x_bar;
    S = Hk*P_bar*Hk' + Rk;
    W = P_bar*Hk'*inv(S);
    
    x_hat = x_bar + W*nu;
    P = P_bar - W*S*W';
end
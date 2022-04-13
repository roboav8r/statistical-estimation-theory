clear all
clc

% Load data
global zhist thist Fk Hk Rk Gammak xhat0 P0 
run('../data/kf_example02b')

% Compute distribution boundaries
alpha = .01; % for 99%
N = 1; %only one iteration
r1 = chi2inv(alpha/2,N*length(zhist))/N
r2 = chi2inv(1-alpha/2,N*length(zhist))/N

% First gamma value
[x_hat_1, P_1] = ps5p5(40, 1);

%Second gamma value
[x_hat_2, P_2] = ps5p5(.4, 2);

%Third gamma value
[x_hat_3, P_3] = ps5p5(.004, 3);

% Convert x_hat from an nx x 1 x nz array into a nz x nx matrix
for jj = 1:size(x_hat_1,3)
   xhat1(jj, :) = x_hat_1(:,:,jj)';
end

for jj = 1:size(x_hat_2,3)
   xhat2(jj, :) = x_hat_2(:,:,jj)';
end

for jj = 1:size(x_hat_3,3)
   xhat3(jj, :) = x_hat_3(:,:,jj)';
end

% Compare differences between best filter (2) with 1 and 3
x_tilde_21 = xhat2-xhat1;
x_tilde_23 = xhat2-xhat3;

err_21 = x_tilde_21'*x_tilde_21
err_23 = x_tilde_23'*x_tilde_23
sqrt(P_2(1,1,51))
sqrt(P_2(2,2,51))


function [x_hat, P] = ps5p5(Qk, index)

    global zhist thist Fk Hk Rk Gammak xhat0 P0 
    
    % Initialize loop and variables
    k = 0;
    x_hat = zeros(2,1,length(zhist));
    P = zeros(2,2,length(zhist));
    
    % Initial propagation step:
    x_bar = Fk*xhat0;
    P_bar = Fk*P0*Fk' + Gammak*Qk*Gammak';
    
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
        
        [x_hat(:,:,k+1),P(:,:,k+1),epsilon_nu(k)] = kf(x_bar(:,:,k), P_bar(:,:,k), zhist(k));
    end %for loop
    
    figure(index)
    subplot(2,1,1)
    plot(thist,xhat1, thist,P1)
    subplot(2,1,2)
    plot(thist,xhat2, thist,P2)
    
    %x_hat_f = x_hat(:,:,50)
    epsilon_nu_bar = mean(epsilon_nu)
    P_f = P(:,:,50)
    %x_RMS = rms(x_hat(2,1,10:50))
end %ps5p5

% Update step for Kalman Filter
function [x_hat,P,epsilon_nu] = kf(x_bar, P_bar, z)
    global Hk Rk
    
    % Compute innovation, nu
    nu = z - Hk*x_bar;

    % Compute Kalman gain, W
    S = Hk*P_bar*Hk' + Rk;
    W = P_bar*Hk'*inv(S);

    % Compute epsilon_nu
    epsilon_nu = nu*inv(S)*nu;
    
    % Compute posterior state estimate and covariance
    x_hat = x_bar + W*nu;
    P = P_bar - W*S*W';

end %kf
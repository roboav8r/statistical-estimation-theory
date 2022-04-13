clear all
clc

% Load data and variables using kf_example script
global zhist thist Fk Hk Rk Gammak xhat0 P0 
run('../data/kf_example02a')

% First gamma value
ps5p3(Qk, 1)

function ps5p3(Qk, index)

    global zhist thist Fk Hk Rk Gammak Qk xhat0 P0 
    
    % Initialize loop and variables
    k = 0;
    x_hat = zeros(2,1,length(zhist));
    P = zeros(2,2,length(zhist));
    
    % Initial propagation step:
    x_bar = Fk*xhat0;
    P_bar = Fk*P0*Fk' + Gammak*Qk*Gammak';
    
    % Initial update step:
    [x_hat(:,:,1),P(:,:,1)] = kf(x_bar, P_bar, zhist(1));
    
    % Remaining predict/update steps for the rest of the data
    for k = 1:length(zhist) 
        % Predict/propagate x and covariance, P (bar means prediction)
        x_bar(:,:,k) = Fk*x_hat(:,:,k);
        P_bar(:,:,k) = Fk*P(:,:,k)*Fk' + Gammak*Qk*Gammak';
        
        % Write to output vector
        xhat1(k) = x_hat(1,1,k);
        xhat2(k) = x_hat(2,1,k);
        P1(k) = sqrt(P(1,1,k));
        P2(k) = sqrt(P(1,1,k));
        
        % Update step using kf function
        [x_hat(:,:,k+1),P(:,:,k+1)] = kf(x_bar(:,:,k), P_bar(:,:,k), zhist(k));
    end %for loop
    
    figure(index)
    subplot(2,1,1)
    plot(thist,xhat1, thist,P1)
    legend('x^{\^}_1(t)','\sigma_1(t)')
    xlabel('Index, k')
    title('Kalman Filter series estimation of LTI state variables')
    subplot(2,1,2)
    plot(thist,xhat2, thist,P2)
    legend('x^{\^}_2(t)','\sigma_2(t)')
    xlabel('Index, k')

    x_hat(:,:,50)
    P(:,:,50)

end %ps5p5

function [x_hat,P] = kf(x_bar, P_bar, z)
    global Hk Rk
    
    % Compute innovation, nu
    nu = z - Hk*x_bar;

    % Compute prior covariance, S and Kalman gain, W
    S = Hk*P_bar*Hk' + Rk;
    W = P_bar*Hk'*inv(S);
    
    % 
    x_hat = x_bar + W*nu;
    P = P_bar - W*S*W';
end %kf
%% John Duncan
% ASE 396P-6: Statistical Estimation Theory
% Final Exam
% Problem 1

% final_p1.m
% Problem set 5, Problem 6: Truth model simulation / KF Consistency

% Clear & configure workspace
clear all
clc
format long

% Add functions folder to file path
addpath ../functions/

% Load system and measurement parameters
global Fk Hk Rk Qk Gammak xhat0 P0 
run('../data/kf_example02a')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 50 monte carlo runs
disp('50 Monte Carlo Runs')
n_monte_carlo_runs = 50;
for ii = 1:n_monte_carlo_runs
    
    % From truth model, generate zhist, thist - length 35
    [xhist,zhist] = mcltisim(Fk,Gammak,Hk,Qk,Rk,Fk*xhat0,P0,35);
    
    % Compute KF state estimate
    [x_hat, P] = ps5p6(zhist, xhat0, P0);

    % Compute errors for x(10) and x(35) for this run & store in ii-th matrix row
    x_tilde_10(ii,:) = (xhist(11,:) - x_hat(11,:));
    x_tilde_35(ii,:) = (xhist(36,:) - x_hat(36,:));
    
    % Get P(10) and P(35) for this run and store in ii-th element of an array
    P_10(:,:,ii) = P(:,:,11);
    P_35(:,:,ii) = P(:,:,36);
        
end

% Compute expected values of x~(10) and x~(35) for these Monte Carlo runs
exp_x_tilde_10 = mean(x_tilde_10)'
exp_x_tilde_35 = mean(x_tilde_35)'

% Compute Exp[x~(10)x~(10)'] and Exp[x~(35)x~(35)'] for these Monte Carlo runs
% Compute running average, then divide by n_monte_carlo_runs
x_tilde_cov_10_sum = [0 0; 0 0];
P_10_sum = [0 0; 0 0];

x_tilde_cov_35_sum = [0 0; 0 0];
P_35_sum = [0 0; 0 0];

for n = 1:n_monte_carlo_runs
    x_tilde_cov_10(:,:,n) = x_tilde_10(n,:)'*x_tilde_10(n,:);
    x_tilde_cov_10_sum = x_tilde_cov_10_sum + x_tilde_10(n,:)'*x_tilde_10(n,:);
    P_10_sum = P_10_sum + P_10(:,:,n);
    
    x_tilde_cov_35(:,:,n) = x_tilde_35(n,:)'*x_tilde_35(n,:);
    x_tilde_cov_35_sum = x_tilde_cov_35_sum + x_tilde_35(n,:)'*x_tilde_35(n,:);
    P_35_sum = P_35_sum + P_35(:,:,n);
    
    % Compute moving average of (1,1) and (2,2) elements of P and
    % exp(x~x~')
    x_tilde_cov_10_bar_1(n) = x_tilde_cov_10_sum(1,1)/n;
    x_tilde_cov_10_bar_2(n) = x_tilde_cov_10_sum(2,2)/n;
    P_10_bar_1(n) = P_10_sum(1,1)/n;
    P_10_bar_2(n) = P_10_sum(2,2)/n;
    
    x_tilde_cov_35_bar_1(n) = x_tilde_cov_35_sum(1,1)/n;
    x_tilde_cov_35_bar_2(n) = x_tilde_cov_35_sum(2,2)/n;
    P_35_bar_1(n) = P_35_sum(1,1)/n;
    P_35_bar_2(n) = P_35_sum(2,2)/n;
    
end
exp_x_tilde_cov_10 = x_tilde_cov_10_sum/n_monte_carlo_runs^2
exp_P_10 = P_10_sum/n_monte_carlo_runs

exp_x_tilde_cov_35 = x_tilde_cov_35_sum/n_monte_carlo_runs^2
exp_P_35 = P_35_sum/n_monte_carlo_runs

% Plot output
n_monte_carlo_vec = 1:1:n_monte_carlo_runs;

figure(1) % x(10), 50 monte carlo runs
subplot(4,1,1)
plot(n_monte_carlo_vec, x_tilde_10(:,1),'b-',[1 50], [exp_x_tilde_10(1) exp_x_tilde_10(1)],'b--')
legend('x~(1)', 'Expected Value (mean)')
title('x_{1}^{~} and expected value')
subplot(4,1,2)
plot(n_monte_carlo_vec, x_tilde_10(:,2),'b-',[1 50], [exp_x_tilde_10(2) exp_x_tilde_10(2)],'b--')
title('x_{2}^{~} and expected value')
legend('x~(2)', 'Expected Value (mean)')
subplot(4,1,3)
plot(n_monte_carlo_vec, P_10_bar_1,'b--',n_monte_carlo_vec, x_tilde_cov_10_bar_1,'b-')
title('P_{11}(10) and x_{1}^{~}^{T}x_{1}^{~}')
legend('Expected Value of P_{11}','Weighted average of x~(1)*x~(1)^{T}')
subplot(4,1,4)
plot(n_monte_carlo_vec,P_10_bar_2,'b--',n_monte_carlo_vec, x_tilde_cov_10_bar_2,'b-')
title('P_{22}(10) and x_{2}^{~}^{T}x_{2}^{~}')
legend('Expected Value of P_{22}','Weighted average of x~(2)*x~(2)^{T}')
sgtitle('x(10), 50 monte carlo runs','t')

figure(2) % x(35), 50 monte carlo runs
subplot(4,1,1)
plot(n_monte_carlo_vec, x_tilde_35(:,1),'b-',[1 50], [exp_x_tilde_35(1) exp_x_tilde_35(1)],'b--')
title('x_{1}^{~} and expected value')
legend('x~(1)', 'Expected Value (mean)')
subplot(4,1,2)
plot(n_monte_carlo_vec, x_tilde_35(:,2),'b-',[1 50], [exp_x_tilde_35(2) exp_x_tilde_35(2)],'b--')
title('x_{2}^{~} and expected value')
legend('x~(2)', 'Expected Value (mean)')
subplot(4,1,3)
plot(n_monte_carlo_vec, P_35_bar_1,'b--',n_monte_carlo_vec, x_tilde_cov_35_bar_1,'b-')
title('P_{11}(35) and x_{1}^{~}^{T}x_{1}^{~}')
legend('Expected Value of P_{11}','Weighted average of x~(1)*x~(1)^{T}')
subplot(4,1,4)
plot(n_monte_carlo_vec,P_35_bar_2,'b--',n_monte_carlo_vec, x_tilde_cov_35_bar_2,'b-')
title('P_{22}(35) and x_{2}^{~}^{T}x_{2}^{~}')
legend('Expected Value of P_{22}','Weighted average of x~(2)*x~(2)^{T}')
sgtitle('x(35), 50 monte carlo runs','t')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1000 monte carlo runs

% 1000 monte carlo runs
disp('1000 Monte Carlo Runs')
n_monte_carlo_runs = 1000;
for ii = 1:n_monte_carlo_runs
    
    %Generate zhist, thist - length 35
    [xhist,zhist] = mcltisim(Fk,Gammak,Hk,Qk,Rk,Fk*xhat0,P0,35);
    [x_hat, P] = ps5p6(zhist, xhat0, P0);

    % Compute errors for x(10) and x(35) for this run & store in matrix row
    x_tilde_10(ii,:) = (xhist(11,:) - x_hat(11,:));
    x_tilde_35(ii,:) = (xhist(36,:) - x_hat(36,:));
    
    % Get P_10 and P_35 for this run and store in ii-th element of an array
    P_10(:,:,ii) = P(:,:,11);
    P_35(:,:,ii) = P(:,:,36);    
   
end

% Compute expected values of x~(10) and x~(35) for these Monte Carlo runs
exp_x_tilde_10 = mean(x_tilde_10)'
exp_x_tilde_35 = mean(x_tilde_35)'

% Compute Exp[x~(10)x~(10)'] and Exp[x~(35)x~(35)'] for these Monte Carlo runs
% Compute running average, then divide by n_monte_carlo_runs
x_tilde_cov_10_sum = [0 0; 0 0];
P_10_sum = [0 0; 0 0];

x_tilde_cov_35_sum = [0 0; 0 0];
P_35_sum = [0 0; 0 0];

for n = 1:n_monte_carlo_runs
    x_tilde_cov_10(:,:,n) = x_tilde_10(n,:)'*x_tilde_10(n,:);
    x_tilde_cov_10_sum = x_tilde_cov_10_sum + x_tilde_10(n,:)'*x_tilde_10(n,:);
    P_10_sum = P_10_sum + P_10(:,:,n);
    
    x_tilde_cov_35(:,:,n) = x_tilde_35(n,:)'*x_tilde_35(n,:);
    x_tilde_cov_35_sum = x_tilde_cov_35_sum + x_tilde_35(n,:)'*x_tilde_35(n,:);
    P_35_sum = P_35_sum + P_35(:,:,n);
    
    % Compute moving average of (1,1) and (2,2) elements of P and
    % exp(x~x~')
    x_tilde_cov_10_bar_1(n) = x_tilde_cov_10_sum(1,1)/n;
    x_tilde_cov_10_bar_2(n) = x_tilde_cov_10_sum(2,2)/n;
    P_10_bar_1(n) = P_10_sum(1,1)/n;
    P_10_bar_2(n) = P_10_sum(2,2)/n;
    
    x_tilde_cov_35_bar_1(n) = x_tilde_cov_35_sum(1,1)/n;
    x_tilde_cov_35_bar_2(n) = x_tilde_cov_35_sum(2,2)/n;
    P_35_bar_1(n) = P_35_sum(1,1)/n;
    P_35_bar_2(n) = P_35_sum(2,2)/n;
    
end
exp_x_tilde_cov_10 = x_tilde_cov_10_sum/n_monte_carlo_runs
exp_P_10 = P_10_sum/n_monte_carlo_runs

exp_x_tilde_cov_35 = x_tilde_cov_35_sum/n_monte_carlo_runs
exp_P_35 = P_35_sum/n_monte_carlo_runs

% Plot output
n_monte_carlo_vec = 1:1:n_monte_carlo_runs;

figure(3) % x(10), 1000 monte carlo runs
subplot(4,1,1)
plot(n_monte_carlo_vec, x_tilde_10(:,1),'b-',[1 50], [exp_x_tilde_10(1) exp_x_tilde_10(1)],'b--')
legend('x~(1)', 'Expected Value (mean)')
title('x_{1}^{~} and expected value')
subplot(4,1,2)
plot(n_monte_carlo_vec, x_tilde_10(:,2),'b-',[1 50], [exp_x_tilde_10(2) exp_x_tilde_10(2)],'b--')
legend('x~(2)', 'Expected Value (mean)')
title('x_{2}^{~} and expected value')
subplot(4,1,3)
plot(n_monte_carlo_vec, P_10_bar_1,'b--',n_monte_carlo_vec, x_tilde_cov_10_bar_1,'b-')
legend('Expected Value of P_{11}','Weighted average of x~(1)*x~(1)^{T}')
title('P_{11}(10) and x_{1}^{~}^{T}x_{1}^{~}')
subplot(4,1,4)
plot(n_monte_carlo_vec,P_10_bar_2,'b--',n_monte_carlo_vec, x_tilde_cov_10_bar_2,'b-')
legend('Expected Value of P_{22}','Weighted average of x~(2)*x~(2)^{T}')
title('P_{22}(10) and x_{2}^{~}^{T}x_{2}^{~}')
sgtitle('x(10), 50 monte carlo runs','t')

figure(4) % x(35), 1000 monte carlo runs
subplot(4,1,1)
plot(n_monte_carlo_vec, x_tilde_35(:,1),'b-',[1 50], [exp_x_tilde_35(1) exp_x_tilde_35(1)],'b--')
legend('x~(2)', 'Expected Value (mean)')
title('x_{1}^{~} and expected value')
subplot(4,1,2)
plot(n_monte_carlo_vec, x_tilde_35(:,2),'b-',[1 50], [exp_x_tilde_35(2) exp_x_tilde_35(2)],'b--')
legend('x~(2)', 'Expected Value (mean)')
title('x_{2}^{~} and expected value')
subplot(4,1,3)
plot(n_monte_carlo_vec, P_35_bar_1,'b--',n_monte_carlo_vec, x_tilde_cov_35_bar_1,'b-')
legend('Expected Value of P_{11}','Weighted average of x~(1)*x~(1)^{T}')
title('P_{11}(35) and x_{1}^{~}^{T}x_{1}^{~}')
subplot(4,1,4)
plot(n_monte_carlo_vec,P_35_bar_2,'b--',n_monte_carlo_vec, x_tilde_cov_35_bar_2,'b-')
legend('Expected Value of P_{22}','Weighted average of x~(2)*x~(2)^{T}')
title('P_{22}(35) and x_{2}^{~}^{T}x_{2}^{~}')
sgtitle('x(35), 50 monte carlo runs','t')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS

function [x_hat, P] = ps5p6(zhist, xhat0, P0)

    global Fk Qk Gammak Qk
    
    % Initialize loop and variables
    x_hat = zeros(length(zhist)+1,2);
    P = zeros(2,2,length(zhist)+1);
    
    % Initialize state estimate
    x_hat(1,:) = xhat0';
    P(:,:,1) = P0;
    
    for k = 1:length(zhist)
        % Propagation / a priori estimates
        x_bar = Fk*x_hat(k,:)';
        P_bar = Fk*P(:,:,k)*Fk' + Gammak*Qk*Gammak';
        
        % Measurement update / a posteriori estimates
        [x_hat_kp1,P(:,:,k+1)] = meas_update(x_bar, P_bar, zhist(k));
        
        x_hat(k+1,:) = x_hat_kp1';
        
    end %for loop
end %ps5p6

function [x_hat,P] = meas_update(x_bar, P_bar, z)
    global Hk Rk
    
    nu = z - Hk*x_bar;
    S = Hk*P_bar*Hk' + Rk;
    W = P_bar*Hk'*inv(S);
    
    x_hat = x_bar + W*nu;
    P = P_bar - W*S*W';
end %kf
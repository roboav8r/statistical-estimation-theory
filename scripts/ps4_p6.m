%% John Duncan
% ASE 396P 
% Problem set 4
% Problem 6

% ps4_p6.m
% Calling file for problem set 4, problem 6

% Clear workspace
clear all
clc

% Part A: Compute x_hat using full Gauss-Newton step size of dx_hat (alpha=1)
disp('Initial guess x_g = [1.5], alpha = 1')
x_g = [1.5];
[x_g_vec_1,J_vec_1] = newton_4_6(x_g, 0, 4, 1)


% Part B: Compute x_hat using partial G-N step size of .1*dx_hat (alpha=.1)
disp('Initial guess x_g = [1.5], alpha = .1')
x_g = [1.5];
[x_g_vec_2,J_vec_2] = newton_4_6(x_g, 0, 4, .75)


% Compute values of the cost function
x = -6:.01:6;
J = abs(0 - atan(x)).^2;

% Plot x_hat vs cost J for both step sizes
plot(x,J,'-',x_g_vec_1,J_vec_1,'o-',x_g_vec_2,J_vec_2,'x-')
title('Effects of step size on Gauss-Newton solver convergence')
legend('Cost Function J','Alpha=1 step size (diverges)','Alpha=.1 step size (converges)')
xlabel('x')
ylabel('Cost J(x) = atan(x)')

function [x_g_vec, J_vec] = newton_4_6(x_g_0, z, n, alpha)
    
    x_g_vec = zeros(1,n);
    x_g_vec(1) = x_g_0;

    %x_g_n = x_g_0; %Initialize loop guess x_g_n
    J_vec(1) = abs(z - atan(x_g_0))^2
    
    for ii=1:n-1
        
        %Compute H
        H = 1/(x_g_vec(ii)^2 + 1);
        
        % Compute V
        V = H'*H;
        
        % Evaluate f_x_hat_g
        f_x_hat_g = H'*(z-atan(x_g_vec(ii)));
        
        % Compute d_x_hat
        d_x_hat = inv(V)*f_x_hat_g;
        
        % Determine alpha
        %alpha=.1;
        
%        f_x = [(x_g_n(1) + x_g_n(2) + x_g_n(1)*x_g_n(2) + 5);
%               (x_g_n(1)^2 + 2*x_g_n(2) - x_g_n(2)^2 -2)]; %Compute f(x_g)
%           
%        df_dx = [(1+x_g_n(2)) (1+x_g_n(1));
%                2*x_g_n(1)    2 - 2*x_g_n(2)]; %compute Jacobian
%        
%        %Compute new guess
        %n = ii
        x_g_vec(ii + 1) = x_g_vec(ii)  + alpha*d_x_hat;
        J_vec(ii+1) = abs(z - atan(x_g_vec(ii+1)))^2;
%        
        %Evaluate cost function J[z - h(x)]
%        f_x = [(x_g_n(1) + x_g_n(2) + x_g_n(1)*x_g_n(2) + 5);
%               (x_g_n(1)^2 + 2*x_g_n(2) - x_g_n(2)^2 -2)];
%        
%        norm_f_x = norm(f_x)
  
    end %For loop
    
    %Assign & return final values
    %x_star = x_g_n;
    %J_final = abs(z - atan(x_star))^2
%    norm_f_x_g = norm_f_x;

end %function
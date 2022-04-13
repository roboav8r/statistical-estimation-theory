%% John Duncan
% ASE 396P 
% Problem set 4
% Problem 1

% ps4_p1.m
% Calling file for problem set 4, problem 1

% Clear workspace
clear all
clc

% Part A
disp('Initial guess x_g = [4; -4]')
x_g = [4; -4];
[x_g_final, norm_f_x_g] = newton_4_1(x_g, 6)

% Part B
disp('Initial guess x_g = [6; 0]')
x_g = [6; 0];
[x_g_final, norm_f_x_g] = newton_4_1(x_g, 6)

% Part C
disp('Initial guess x_g = [-5; 5]')
x_g = [-5; 5];
[x_g_final, norm_f_x_g] = newton_4_1(x_g, 6)


function [x_g_final, norm_f_x_g] = newton_4_1(x_g_0, n)

   x_g_n = x_g_0; %Initialize loop guess x_g_n
   for ii=1:n
       f_x = [(x_g_n(1) + x_g_n(2) + x_g_n(1)*x_g_n(2) + 5);
              (x_g_n(1)^2 + 2*x_g_n(2) - x_g_n(2)^2 -2)]; %Compute f(x_g)
          
       df_dx = [(1+x_g_n(2)) (1+x_g_n(1));
               2*x_g_n(1)    2 - 2*x_g_n(2)]; %compute Jacobian
       
       %Compute new guess
       n = ii
       x_g_n = x_g_n  - inv(df_dx)*f_x
       
       %Evaluate new norm |f(x_g)|
       f_x = [(x_g_n(1) + x_g_n(2) + x_g_n(1)*x_g_n(2) + 5);
              (x_g_n(1)^2 + 2*x_g_n(2) - x_g_n(2)^2 -2)];
       
       norm_f_x = norm(f_x)
 
   end %For loop
   
   %Assign & return final values
   x_g_final = x_g_n;
   norm_f_x_g = norm_f_x;
   
   
end %function
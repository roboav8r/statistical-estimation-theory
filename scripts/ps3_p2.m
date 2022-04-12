%% John Duncan
% ASE 396P 
% Exam 1
% Problem 8

% ps3_p2.m
% Calling file for problem set 3, problem 2 / Exam problem #8

% Clear workspace
clear all
clc

% Given information

zprime = [34392.7080226075
95049.241211972
-165957.584445193
12813.3430731501
33744.0080015867
-107847.178062601];

Hprime = [537.6671395461 -303.514415613978 2118.36542708511
1833.88501459509 239.837126577055 5422.13861846274
-2258.84686100365 2504.87785780803 -10638.8090086238
862.173320368121 1938.60592091941 -33.0051183925756
318.765239858981 -944.920858109565 2337.75095208492
-1307.68829630527 2124.4464264323 -7112.64499263293];

R = [4.48969254194419e-06 -5.3639754880881e-08 -2.82108367550151e-07 1.15527795296165e-07 -2.49965420713856e-06 6.36392589368105e-07
-5.3639754880881e-08 6.07007287029812e-06 -2.62065228561023e-06 3.81520661439239e-06 4.6837147858155e-07 8.24799673475808e-07
-2.82108367550151e-07 -2.62065228561023e-06 6.68658685773165e-06 -8.8120668398492e-07 2.99242554994465e-07 2.08470773900216e-08
1.15527795296165e-07 3.81520661439239e-06 -8.8120668398492e-07 6.07371651695821e-06 2.63707372141155e-06 3.55951797531251e-07
-2.49965420713856e-06 4.6837147858155e-07 2.99242554994465e-07 2.63707372141155e-06 4.55596815512913e-06 -1.40588403974735e-06
6.36392589368105e-07 8.24799673475808e-07 2.08470773900216e-08 3.55951797531251e-07 -1.40588403974735e-06 2.43477996714807e-06];

%%% Main routine
disp('Normal Least Squares Solution and Covariance:')
[x_hat_LS, P_LS] = least_squares(Hprime,zprime,R)

disp('Square Root Information Based Solution and Covariance:')
[x_hat_SRIBLS,Rotilde] = sribls(Hprime,zprime,R)
P_SRIBLS=inv(Rotilde)*inv(Rotilde')

%%%% FUNCTIONS
function [x_hat_LS, P_LS] = least_squares(Hprime,zprime,R)

    %%% Header removed for brevity %%%
    %Function to compute the solution to the normal equations (least squares)
    P_LS = inv(Hprime'*inv(R)*Hprime); % Compute covariance
    x_hat_LS = P_LS*Hprime'*inv(R)*zprime;

end %least_squares

function [xhat,Rotilde] = sribls(Hprime,zprime,R)
    % sribls : Square-root-information-based least squares routine.
    
    %%% Header removed for brevity %%%
    
    % Transform measurement model
    R_a = chol(R); %Compute Cholesky factorization
    z = inv(R_a')*zprime;
    H = inv(R_a')*Hprime;
    
    % Perform QR decomposition of H
    [Qtilde,Rtilde] = qr(H);
    ztilde = Qtilde'*z;
    
    % Extract Rotilde matrix from Rtilde; extract zotilde from ztilde 
    % (the first n rows, where n = number of columns in h)
    [m,n]=size(H);
    Rotilde = Rtilde(1:n,:);
    zotilde = ztilde(1:n,:);
    
    % Compute x_hat_SRIBLS
    xhat=inv(Rotilde)*zotilde;

end % sribls

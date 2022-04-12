%% John Duncan
% ASE 396P 
% Exam 1
% Problem 3

% ps2_p5.m
% Matlab script for problem set 2, problem 5 / Exam problem #3

% Clear workspace
clear all
clc

%Initial probability of false alarm:
P_F = .01;

disp('Part C: sigma ratio 4')
P_D = chi2detector(1, 4, P_F, 1) %ratio of 1/4 with given P_F and chi2 degree=1

disp('Part C: sigma ratio 6')
P_D = chi2detector(1, 6, P_F, 1) %ratio of 1/6 with given P_F and chi2 degree=1

disp('Part C: sigma ratio 8')
P_D = chi2detector(1, 8, P_F, 1) %ratio of 1/8 with given P_F and chi2 degree=1

disp('Part C: sigma ratio 10')
P_D = chi2detector(1, 10, P_F, 1) %ratio of 1/10 with given P_F and chi2 degree=1

disp('Part D: multiple measurements with sigma ratio 4')
n=1; %initialize count
while P_D < .99 %Keep running loop until P_D exceeds 99%
    P_D = chi2detector(1, 4, .001, n); %ratio of 1/4 with given P_F and chi2 degree=n
    n = n+1; %Increase n
end
n %Print n

function PD = chi2detector(sigma0, sigma1, PF, n)

    %%% Inputs:
    % sigma0, sigma1: deviations of Brahms / Stravinsky, respectively
    % PF: probability of false alarm
    % n: degree of chi-squared variable
    
    %%% Output:
    % PD: probability of detection
    
    nu_star = chi2inv((1-PF), n);
    PD = 1 - chi2cdf(nu_star*(sigma0^2/sigma1^2), n);

end %chi2detector


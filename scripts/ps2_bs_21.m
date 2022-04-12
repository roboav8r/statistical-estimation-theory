
%% John Duncan
% ASE 396P 
% Exam 1
% Problem 5

% ps2_bs_21.m
% Calling file for Bar-Shalom problem 2-1 / Exam problem #5

% Clear workspace
clear all
clc

disp('Part A')
p_1 = .5;
p_2 = 1-p_1;
sigma=1;
z=1.5;
[x_MAP, x_MMSE] = bs21(p_1, p_2, sigma, z)

disp('Part B')
p_1 = .5;
p_2 = 1-p_1;
sigma=1;
z=3;
[x_MAP, x_MMSE] = bs21(p_1, p_2, sigma, z)


disp('Part C')
p_1 = .3;
p_2 = 1-p_1;
sigma=1;
z=1.5;
[x_MAP, x_MMSE] = bs21(p_1, p_2, sigma, z)


disp('Part D')
p_1 = .5;
p_2 = 1-p_1;
sigma=.1;
z=1.8;
[x_MAP, x_MMSE] = bs21(p_1, p_2, sigma, z)





function [x_MAP, x_MMSE] = bs21(p1, p2, sigma, z)

%%% Compute x_MAP
%Evaluate both numerator A Posteriori components in a vector
x = [1 2]';
x_AP_vec = [p1*exp(-1/2*((z - x(1))/sigma)^2) p2*exp(-1/2*((z - x(2))/sigma)^2)]';

% Find maximum value of x_MAP vector, then get the index and select the
% associated x_MAP value
if x_AP_vec(1)==x_AP_vec(2); % If there are two maxima
    x_MAP = x; % Return vector / both values as maximum value
else % Choose actual maximum
    [MAP, MAP_index] = max(x_AP_vec);
    x_MAP = x(MAP_index);
end

% Compute MSE as derived in notes
%MSE = [(p1+4*p2) - 2.*x_MAP.*(p1+2*p2) + x_MAP.^2] ./ (sqrt(2*pi)*sigma)

%%% Compute x_MMSE
x_MMSE = [(z-1)*p1*exp(-1/2*((z - 1)/sigma)^2);
          (z-2)*p2*exp(-1/2*((z - 2)/sigma)^2)]./ (sqrt(2*pi)*sigma);
      
x_MMSE=sum(x_MMSE);

% Compute variances associated with x_MMSE.
% NOTE: Compares x_MMSE to both potential values of x and returns a vector
% of variance values.
%var = (x-x_MMSE).^2

end
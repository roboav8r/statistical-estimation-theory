% Clear workspace
clear all
clc

%%% Part 1: 100,000 points uniformly sampled between [0,100]
% Generate random data
n_points_1 = 100000; % Number of data points for part 1
range_1 = 100; % Maximum range for part 1
% Create random vector with specified size and range
r_1 = range_1*rand(1,n_points_1);

% Plot Data
figure(1)
hist_1 = histogram(r_1,range_1); % Convert data into a histogram
b = bar([1:range_1],hist_1.Values./n_points_1); % Normalize by number of occurrences
title('Problem 1.8.1: 100,000 Uniformly distributed points in range [0,100]')
xlabel('Number')
ylabel('Probability')

%%% Second part
% Generate random data (average of 2 points)
n_points_2 = 100000;
range_2 = 100;
r_2 = zeros(1,n_points_2);
for n_2=1:n_points_2
    % Take average of 2 random numbers in range and assign it to r_2
    r_2(n_2) = mean(range_2*rand(1,2));
end % for

% Plot Data
figure(2)
hist_2 = histogram(r_2,range_2);
b = bar([1:range_2],hist_2.Values./n_points_2);
title('Problem 1.8.2: 100,000 points, average of 2 uniformly sampled points in range [0,100]')
xlabel('Number')
ylabel('Probability')

%%% Third part
% Generate random data (averaged)
n_points_3 = 100000;
range_3 = 100;
r_3 = zeros(1,n_points_3);
for n_3=1:n_points_3
    % Take average of 100 random numbers in range and assign it to r_3
    r_3(n_3) = mean(range_3*rand(1,100));
end % for

% Plot Data
figure(3)
hist_3 = histogram(r_3,range_3);
b = bar([1:range_3],hist_3.Values./n_points_3);
title('Problem 1.8.3: 100,000 points, average of 100 uniformly sampled points in range [0,100]')
xlabel('Number')
ylabel('Probability')

%%% Fourth part
% Generate random data (averaged)
n_points_4 = 100000;
range_4 = 100;
r_4 = zeros(1,n_points_4);
for n_4=1:n_points_4
    % Take average of 1000 random numbers in range and assign it to r_4
    r_4(n_4) = mean(range_4*rand(1,1000));
end % for

% Plot Data
figure(4)
hist_4 = histogram(r_4,range_4);
b = bar([1:range_4],hist_4.Values./n_points_4);
title('Problem 1.8.4: 100,000 points, average of 1000 uniformly sampled points in range [0,100]')
xlabel('Number')
ylabel('Probability')
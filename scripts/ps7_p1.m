clear all
clc

% Use long precision
format long

% Add functions to path
addpath ../functions/

% Define constants
r = 6378000; %m
mu = 3.986005e14;

% Define initial conditions (HEO with 12-hour period, AKA the poor man's
% GEO)
x0 = [(39700 + 6378)*1000;
    0;
    0;
    0;
    1500;
    0];

% Define propagator parameters
T = 60; %seconds; 1 minute


% Initial propagation step
uk_0 = [0 0 0]';
vk_0 = [0 0 0]';

[x_1,Fk_1,GAMMAk_1]=propagateOrbit(0,T,x0,uk_0,vk_0,mu);
x_k(1,:) = x_1';

for ii=2:720 %minutes in 12 hours
    [x_loop,Fk_loop,GAMMAk_loop]=propagateOrbit(T*(ii-1),T,x_k(ii-1,:)',uk_0,vk_0,mu);
    x_k(ii,:) = x_loop';
end

figure(1)
hold on
% Plot X vs y altitude
plot(x_k(:,1),x_k(:,2))
%axis([-5e7 5e7 -5e7 5e7])
%axis square

% plot earth radius
th = 0:pi/50:2*pi;
x_earth = r * cos(th);
y_earth = r * sin(th);
plot(x_earth, y_earth);
title('Propagated Orbit around Earth')

hold off

% Compute % magnitude of difference between initial and final radii
r_initial = x_k(1,1:3);
r_final = x_k(length(x_k),1:3);

delta_r = r_final - r_initial;
pct_delta_r = norm(delta_r)*100/norm(r_initial)


% Validate F(k)
uk_0 = [0 0 0]';
vk_0 = [0 0 0]';

dt = 60; %seconds
dx = .1; %x_hat = x0 + (1 + dx); in this case, 10%

[x_1,Fk_1,GAMMAk_1]=propagateOrbit(0,dt,x0,uk_0,vk_0,mu); %propagate x0

x_hat_0 = x0.*[1+dx]; % Create x_hat by multiplying x0 by 10%, then propagate
[x_hat_1,Fk_hat_1,GAMMAk_hat_1]=propagateOrbit(0,dt,x_hat_0,uk_0,vk_0,mu);

dx_1_actual = x_1 - x_hat_1
dx_1_linear = Fk_1*(x_1 - x_hat_1)

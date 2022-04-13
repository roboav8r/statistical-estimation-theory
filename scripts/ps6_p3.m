format long

clear all
clc

% Load data
global Fk Gammak Hk P0 Qk Rk thist xhat0 zhist
run('../data/kf_example03a')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KALMAN FILTER, DATA SET A
% Initial propagation step
[x_bar, P_bar] = kf_prop(xhat0, P0);

% Initial update step
[x_hat(:,:,1),P(:,:,1)] = kf_update(x_bar, P_bar, zhist(1));

for k = 2:length(zhist) 
    [x_bar, P_bar] = kf_prop(x_hat(:,:,k-1), P(:,:,k-1));
    [x_hat(:,:,k),P(:,:,k)] = kf_update(x_bar, P_bar, zhist(k));
end %for loop

disp('Kalman filter results - data set a')
% Print output
x_hat(:,:,length(zhist))
P(:,:,length(zhist))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SRIF, DATA SET A
disp('SRIF - data set a')
% Generate initial R_vv, R_xx matrices
R_xx_0 = [inv(chol(P0))]';
R_vv_0 = [inv(chol(Qk))]';

z_x_0 = R_xx_0*xhat0;

% Initial propagation step
%[z_bar_x_0, R_bar_xx_0] = srif_prop(R_vv_0, R_xx_0,z_x_0);
[z_bar_x(:,:,1), z_bar_v_0, R_bar_xx(:,:,1), R_bar_vx(:,:,1), R_bar_vv_0] = srif_prop(R_vv_0, R_xx_0, z_x_0);

% Initial update step
[z_x(:,:,1), R_xx(:,:,1)] = srif_update(z_bar_x(:,:,1), R_bar_xx(:,:,1), zhist(1));

% Remaining updates
for k = 2:length(zhist)
    [z_bar_x(:,:,k), z_bar_v(:,:,k-1), R_bar_xx(:,:,k), R_bar_vx(:,:,k), R_bar_vv(:,:,k-1)] = srif_prop(R_vv_0, R_xx(:,:,k-1), z_x(:,:,k-1));
    %[z_bar_x, R_bar_xx] = srif_prop(R_vv_0, R_xx(:,:,k-1),z_x(:,:,k-1));
    [z_x(:,:,k), R_xx(:,:,k)] = srif_update(z_bar_x(:,:,k), R_bar_xx(:,:,k), zhist(k));
end

% Print output
inv(R_xx(:,:,length(zhist)))*z_x(:,:,length(zhist))
inv(R_xx(:,:,length(zhist)))*inv(R_xx(:,:,length(zhist)))'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KALMAN FILTER SMOOTHING
N = length(zhist);

x_KF_smooth(N,:) = x_hat(:,:,N)';

for kk = N-1:-1:1 % Step backwards from N-1
    
    % Get P_bar at index kk
    [x_bar, P_bar] = kf_prop(x_hat(:,:,kk), P(:,:,kk));
    
    % Compute smoother gain
    C_KF = P(:,:,kk)*Fk'*inv(P_bar);
    
    % Compute x_KF_smooth
    x_KF_smooth(kk,:) = (x_hat(:,:,kk) + C_KF*(x_KF_smooth(kk+1,:)' - x_bar))';
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SRIF SMOOTHING

% Initialize smoothed data with final SRIF filter value
z_x_smooth(N,:) = z_x(:,:,N)';
R_xx_smooth(:,:,N) = R_xx(:,:,N);

x_SRIF_smooth(N,:) = (inv(R_xx_smooth(:,:,N))*z_x_smooth(N,:)')';
P_SRIF_smooth(:,:,N) = inv(R_xx_smooth(:,:,N)*inv(R_xx_smooth(:,:,N)))';

nv = length(R_vv_0);
nx = length(R_xx_0);

for kk = N-1:-1:1 % Step backwards from N-1
    
    % Build block matrix to be factorized
    R_block = [(R_bar_vv(:,:,kk) + R_bar_vx(:,:,kk+1)*Gammak) (R_bar_vx(:,:,kk+1)*Fk);
        (R_xx_smooth(:,:,kk+1)*Gammak) (R_xx_smooth(:,:,kk+1)*Fk)];
    
    % QR transform block matrix
    [Q_a, R_a] = qr(R_block);
    T_a = Q_a';
    
    % Extract R_xx_smooth from R_a matrix
    R_xx_smooth(:,:,kk) = R_a((nv+1:nv+nx),(nv+1:nv+nx));
    
    % Generate and transform z matrix
    z_vector = [z_bar_v(:,:,kk);
        z_x_smooth(kk+1,:)'];
    z_transformed = T_a*z_vector;
    
    % Extract z_x_smooth from transformed z matrix
    z_x_smooth(kk,:) = z_transformed(nv+1:nv+nx)'
    
    % Compute x_SRIF_smooth and P_SRIF_smooth
    x_SRIF_smooth(kk,:) = (inv(R_xx_smooth(:,:,kk))*z_x_smooth(kk,:)')';
    P_SRIF_smooth(:,:,kk) = inv(R_xx_smooth(:,:,kk)*inv(R_xx_smooth(:,:,kk)))';
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate data for plotting

% KF data - forward pass
for ii = 1:length(zhist)
    x_KF(ii,:) = x_hat(:,:,ii)';
end

% SRIF data - forward pass
for ii = 1:length(zhist)
    x_SRIF(ii,:) = (inv(R_xx(:,:,ii))*z_x(:,:,ii))';
end

% Plot
figure(1)
subplot(3,1,1)
plot(thist,x_KF(:,1),'b-',thist,x_KF_smooth(:,1),'r-')
title('Kalman Filter & Smoother')
ylabel('x_1')
subplot(3,1,2)
plot(thist,x_KF(:,2),'b-',thist,x_KF_smooth(:,2),'r-')
ylabel('x_2')
subplot(3,1,3)
plot(thist,x_KF(:,3),'b-',thist,x_KF_smooth(:,3),'r-')
ylabel('x_3')
legend('Filtered data','Smoothed Data')

figure(2)
subplot(3,1,1)
plot(thist,x_SRIF(:,1),'b-',thist,x_SRIF_smooth(:,1),'r-')
title('SRI Filter & Smoother')
ylabel('x_1')
subplot(3,1,2)
plot(thist,x_SRIF(:,2),'b-',thist,x_SRIF_smooth(:,2),'r-')
ylabel('x_2')
subplot(3,1,3)
plot(thist,x_SRIF(:,3),'b-',thist,x_SRIF_smooth(:,3),'r-')
ylabel('x_3')
legend('Filtered data','Smoothed Data')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS
function [x_bar, P_bar] = kf_prop(xhat, P)
    global Fk Gammak Qk
    x_bar = Fk*xhat;
    P_bar = Fk*P*Fk' + Gammak*Qk*Gammak';
end % prop function

function [x_hat,P] = kf_update(x_bar, P_bar, z)
    global Hk Rk
    
    nu = z - Hk*x_bar;
    S = Hk*P_bar*Hk' + Rk;
    W = P_bar*Hk'*inv(S);
    epsilon_nu = nu*inv(S)*nu;
    
    x_hat = x_bar + W*nu;
    P = P_bar - W*S*W';
end %kf

function [z_bar_x, z_bar_v, R_bar_xx, R_bar_vx, R_bar_vv] = srif_prop(R_vv, R_xx, z_x)
    global Fk Gammak
    n_z = 1;
    n_v = length(R_vv);
    n_x = length(R_xx);
    
    % QR factorize
    [Q_a, R_a] = qr([[R_vv zeros(n_v,n_x)];
        [-R_xx*inv(Fk)*Gammak R_xx*inv(Fk)]]);
    
    % Perform orthonormal transformation
    z_matrix = Q_a'*[zeros(n_v,n_z); z_x];
    
    % Extract & Return R_bar_xx and z_bar_x
    z_bar_x = z_matrix((n_v+1:n_v+n_x),:); %k+1
    R_bar_xx = R_a((n_v+1:n_v+n_x),(n_v+1:n_v+n_x)); %k+1
    
    % And additional matrices for smoothing
    z_bar_v = z_matrix((1:n_v),:);%k
    R_bar_vv = R_a((1:n_v),(1:n_v)); %k
    R_bar_vx = R_a((1:n_v),(n_v+1:n_v+n_x)); %k+1

end % prop function

function [z_x, R_xx] = srif_update(z_bar_x, R_bar_xx, z)
    global Hk Rk
    n_z = 1;
    %n_v = length(R_vv);
    n_x = length(R_bar_xx);
    
    % Compute cholesky factorization of R
    R_a = chol(Rk);
    
    % Transform z and H
    z_a = inv(R_a)'*z;
    H_a = inv(R_a)'*Hk;
    
    % Compute next QR factorization
    [Q_b, R_b] = qr([R_bar_xx; H_a]);
    
    % transform z matrix
    z_x_r = Q_b'*[z_bar_x; z_a];
    
    % Extract z_x and R_xx
    z_x = z_x_r(1:n_x,:);
    R_xx = R_b(1:n_x,:);
    
    
end %kf
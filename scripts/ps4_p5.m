kmax=50;

global f sigma
f = .9; % Value less than one
sigma = .25;
P_0 = 1;
P(1) = P_0*f^2 + (1-f^2)*sigma;

for k = 2:kmax
    P(k) = P(k-1)*f^2 + (1-f^2)*sigma;
end

plot(P)
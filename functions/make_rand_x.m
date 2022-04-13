function rand_x = make_rand_x(x_bar,P_xx)

n = length(x_bar);

z = randn(n,1);

A = chol(P_xx);

rand_x = A*z + x_bar;

end
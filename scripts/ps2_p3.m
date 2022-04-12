% Setup
format long

% Define matrices and vectors
A = [1 1; 1 1e-8]
B = [1 1; 1 1e-3]

y_A = [1; 1-3e-8]
y_B = [1; 1-3e-3]

% Initial observation: Compute condition numbers
cond_A = cond(A)
cond_B = cond(B)
cond_AA = cond(A'*A)
cond_BB = cond(B'*B)

% 3b Left division operator
x_A_leftdivision = A\y_A
x_B_leftdivision = B\y_B

% 3c Least Squares solution
x_A_leastsquares = (A'*A)\(A'*y_A)
x_B_leastsquares = (B'*B)\(B'*y_B)

% 3d Least Squares inversion
x_A_leastsquaresinv=inv(A'*A)*(A'*y_A)
x_B_leastsquaresinv=inv(B'*B)*(B'*y_B)


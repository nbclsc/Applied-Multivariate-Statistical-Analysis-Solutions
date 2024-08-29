C = [1 -1; 1 1]';
X = [6 10 8; 9 6 3]';
X = X*C';
n = height(X);
p = width(X);
xbar = mean(X)';
S = cov(X);
mu0 = [9 5]';
mu0 = C*mu0;

% Compute T^2.
T2 = n*(xbar - mu0)'/S*(xbar - mu0);

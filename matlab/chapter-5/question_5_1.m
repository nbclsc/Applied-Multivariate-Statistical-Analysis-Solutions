% (a.)
X = [2 8 6 8; 12 9 9 10]';
n = height(X);
p = width(X);
xbar = mean(X)';
S = cov(X);
mu0 = [7 11]';

% Compute T^2.
T2 = n*(xbar - mu0)'/S*(xbar - mu0);

% (c.)
an_alpha = 0.05;
critical_value = (((n-1)*p)/(n-p))*icdf('F', 1 - an_alpha, 2, 2);

% This returns 0 (false), so we fail to reject H0, that mu = mu0.
T2 > critical_value
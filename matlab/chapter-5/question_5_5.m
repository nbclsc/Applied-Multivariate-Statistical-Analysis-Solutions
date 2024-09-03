closed = readmatrix(fullfile(data_folder, 'Table4.1.xlsx'));
open = readmatrix(fullfile(data_folder, 'Table4.5.xlsx'));
X = [closed(:,2).^0.25 open(:,2).^0.25];

n = height(X);
p = width(X);

xbar = mean(X)';
S = cov(X);
Sinv = inv(S);

mu0 = [0.55 0.60]';
an_alpha = 0.05;

T2 = n*(xbar - mu0)'*Sinv*(xbar - mu0);

const = (p*(n-1))/(n-p);
f_crit = icdf('F', 1-an_alpha, p, n-p);

const*f_crit

% Returns 0, so fail to reject H0 (the null hypothesis).
T2 > const*f_crit
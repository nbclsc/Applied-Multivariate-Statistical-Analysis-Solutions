t0501 = readmatrix(fullfile(data_folder, 'Table5.1.xlsx'));
X = t0501(:,2:end);

n = height(X);
p = width(X);

% H0: mu = m0
% H1: mu not equal m0

mu0 = [4 50 10]';
an_alpha = 0.10;

xbar = mean(X)';
S = cov(X);

% Compute the Hotelling's T2 value.
T2 = n*(xbar - mu0)'/S*(xbar - mu0);

% Compute the critical value from the reference distribution.
critical_value = (n - 1)*p/(n - p)*icdf('F', 1-an_alpha, p, n - p);

% Returns 1, so reject H0 in facor of H1.
T2 > critical_value

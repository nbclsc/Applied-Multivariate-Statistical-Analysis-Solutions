X = [2 8 6 8; 12 9 9 10]';
mu0 = [7 11]';
[n, p] = size(X);

xbar = mean(X)';

S = zeros(p, p);
for i = 1:n
    S = S + (X(i,:)' - xbar)*(X(i,:)' - xbar)';
end

S0 = zeros(p, p);
for i = 1:n
    S0 = S0 + (X(i,:)' - mu0)*(X(i,:)' - mu0)';
end

% (a)

T2 = (n-1)*det(S0)/det(S) - (n-1);
T2

% (b)

% (5-13)
lmbda = (det(S)/det(S0))^(n/2);
lmbda

wilks_lmbda = det(S)/det(S0);
wilks_lmbda


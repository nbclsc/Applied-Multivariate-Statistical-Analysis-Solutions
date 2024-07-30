% Exercise 4.4 (a)
mu = [2 -3 1]';
Sigma = [1 1 1; 1 3 2; 1 2 2];
b = [3 -2 1]';
b'*mu
b'*Sigma*b

% Exercise 4.4 (b)
A = [0 1 0; -1 1/3 0];
A*Sigma*A'

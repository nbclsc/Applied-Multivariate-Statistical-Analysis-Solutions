output_path = 'C:\Users\nbclsc\Desktop\applied-multivariate-statistics\solutions\matlab\chapter-4\';
sort(diag((Xex4_12 - ones(10,1)*xbar')*Sinv*(Xex4_12 - ones(10,1)*xbar')'))

V = [2 -0.8*sqrt(2); -0.8*sqrt(2) 1]
inv(V)

% Exercise 4.2 (c) 
mu = [0;2];
Sigma = [2 sqrt(2)/2; sqrt(2)/2 1];
[V,D] = eig(Sigma);
c = sqrt(chi2inv(0.5,2));
MyPlotEllipse(V,mu,c*sqrt(D(1,1)),c*sqrt(D(2,2)),[-2 2],[0 4],output_path,'sol4.2')

% Exercise 4.3 (d)
Sigma = [1 -2 0; -2 5 0; 0 0 2]
A = [0.5 0.5 0; 0 0 1];
A*Sigma*A'

% X11 = [1 -2; -2 5];
% b = [0.5 0.5]';
% b'*X11*b

% Exercise 4.3 (e)
A = [-5/2 1 -1; 0 1 0];
A*Sigma*A'

% Exercise 4.4 (b)
mu = [2 -3 1]';
Sigma = [1 1 1; 1 3 2; 1 2 2];
b = [3 -2 1]';
b'*mu
b'*Sigma*b

% Exercise 4.4 (b)
A = [0 1 0; -1 1/3 0];
A*Sigma*A'

% Exercise 4.8
x = [-3:.1:3];
y = normpdf(x,0,1);


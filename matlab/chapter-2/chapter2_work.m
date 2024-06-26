output_path = 'C:\Users\nbclsc\Desktop\applied-multivariate-statistics\solutions\matlab\chapter-2\';

% Exercsie 2.1(a)
x = [5,1,3]'; y = [-1,3,1]';

starts = zeros(2,3);  % Starts at the origin.
ends = [x'; y'];  % Ends at hte point.

% quiver3 args are x,y,z,u,v,w. x,y,z are are start positions and u,v,w are
% end positions.
a = quiver3(starts(:,1), starts(:,2), starts(:,3), ends(:,1), ends(:,2), ends(:,3));
axis equal
saveas(a, append(output_path, 'sol2.1a.png'), 'png')


% Exercsie 2.1(c)
starts = zeros(2,3);  % Starts at the origin.
ends = [(x-mean(x))'; (y-mean(y))'];  % Ends at the point. Subtract the mean values.

% quiver3 args are x,y,z,u,v,w. x,y,z are the start positions and u,v,w are
% the end positions.
a = quiver3(starts(:,1), starts(:,2), starts(:,3), ends(:,1), ends(:,2), ends(:,3));
axis equal
saveas(a, append(output_path,'sol2.1c.png'), 'png')


% Exercise 2.18
A = [4 -sqrt(2); -sqrt(2) 3];
[V,D] = eig(A);
rref(A - D(1,1)*eye(width(A)))
rref(A - D(2,2)*eye(width(A)))
c = 1;
MyPlotEllipse(V,[0;0],c/sqrt(D(1,1)),c/sqrt(D(2,2)),[-c c],[-c c],output_path,'sol2.18.c1')
c = 4;
MyPlotEllipse(V,[0;0],c/sqrt(D(1,1)),c/sqrt(D(2,2)),[-c c],[-c c],output_path,'sol2.18.c4')

% Exercise 2.20
A = [2 1; 1 3]
eig(A)
[V,D] = eig(A);
rref(A - D(1,1)*eye(width(A)))
rref(A - D(2,2)*eye(width(A)))
a = 5 - sqrt(5);
b = 5 + sqrt(5);
A12 = V*sqrt(D)*V';

% Exercise 2.21
A = [1 1; 2 -2; 2 2]
A'*A
A*A'
[V,D] = eig(A*A');

% Exercise 2.22
A = [4 8 8; 3 6 -9]
A'*A
A*A'
[V,D] = eig(A*A');

% Exercise 2.24
A = [4 0 0; 0 9 0; 0 0 1]
[P,D] = eig(A)
[P,D] = eig(inv(A))

% Exercise 2.32
popsig = [4 -1 0.5 -0.5 0; -1 3 1 -1 0; 0.5 1 6 1 -1; -0.5 -1 1 4 0; 0 0 -1 0 2]
popmu = [2 4 -1 3 0]'
A = [1 -1; 1 1]
B = [1 1 1; 1 1 -2]
x1 = diag([1 1 0 0 0])  % Use the identity matrix instead of slicing things up.
x2 = diag([0 0 1 1 1])
x1*popmu % (a)
[A zeros(2,3)]*(x1*popmu) % (b)
x1*popsig*x1 % (c)
[A zeros(2,3)]*(x1*popsig*x1)*[A zeros(2,3)]' % (d)
x2*popmu % (e)
[zeros(2,2) B]*(x2*popmu) % (f)
x2*popsig*x2 % (g)
[zeros(2,2) B]*(x2*popsig*x2)*[zeros(2,2) B]' % (h)
x1*popsig*x2 % (i)
[A zeros(2,3)]*(x1*popsig*x2)*[zeros(2,2) B]' % (j)

% Exercise 2.33
A = [2 -1 0; 1 1 3]
B = [1 2; 1 -1]
x1 = diag([1 1 1 0 0])  % Use the identity matrix instead of slicing things up.
x2 = diag([0 0 0 1 1])
x1*popmu % (a)
[A zeros(2,3)]*(x1*popmu) % (b)
x1*popsig*x1 % (c)
[A zeros(2,2)]*(x1*popsig*x1)*[A zeros(2,2)]' % (d)
x2*popmu % (e)
[zeros(2,3) B]*(x2*popmu) % (f)
x2*popsig*x2 % (g)
[zeros(2,3) B]*(x2*popsig*x2)*[zeros(2,3) B]' % (h)
x1*popsig*x2 % (i)
[A zeros(2,2)]*(x1*popsig*x2)*[zeros(2,3) B]' % (j)

% Exercise 2.34
b = [2 -1 4 0]'
d = [-1 3 -2 1]'
(b'*d) <= (b'*b)*(d'*d)

% Exercise 2.35
b = [-4 3]'
d = [1 1]'
B = [2 -2; -2 5]
(b'*d)^2 <= (b'*B*b)*(d'*inv(B)*d)

% Exercise 2.35
B = [4 3; 3 4]
[P,D] = eig(B)
max(diag(D))
min(diag(D))


% Exercise 2.38
A = [13 -4 2; -4 13 -2; 2 -2 10]
roots([1 -36 405 -1458])
[P,D] = eig(A)

% Exercise 2.41
popmu = [3 2 -2 0]'
popsig = 3*eye(4)
A = [1 -1 0 0; 1 1 -2 0; 1 1 1 -3]
A*popmu
A*popsig*A'

% Exercise 2.42
popmu = [3 2 -2 0]'
popsig = 2*eye(4) + ones(4)
A = [1 -1 0 0; 1 1 -2 0; 1 1 1 -3]
A*popmu
A*popsig*A'

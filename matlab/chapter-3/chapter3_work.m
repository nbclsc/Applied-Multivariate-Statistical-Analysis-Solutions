output_path = 'C:\Users\nbclsc\Desktop\applied-multivariate-statistics\solutions\matlab\chapter-3\';
% Exercise 3.1 (a)
X = [9 1; 5 3; 1 2];
mean_pt = mean(X);
hold on
    % Plot data.
    scatter(X(:,1),X(:,2),'blue','filled')
    % Plot mean.
    plot(mean_pt(1),mean_pt(2),'d')
    % Text for mean.
    text(mean_pt(1)-0.5,mean_pt(2)-0.2, ...
        join(["$$\bar{x}=\left[\begin{array}{c}",mean_pt(1),"\\",mean_pt(2),"\end{array}\right]$$"],''), ...
        'interpreter','latex')
    % Text for data.
    for r = 1:height(X)
        anno = join(["$$x_",r," = \left[\begin{array}{c}",X(r,1),"\\",X(r,2),"\end{array}\right]$$"],'');
        text(X(r,1)-0.5,X(r,2)-0.2,anno,'interpreter','latex');
    end
    xlim([0 10])
    ylim([0 4])
    saveas(gcf,'sol3.1a.png')
hold off

% Exercise 3.1 (b)
% Compute the deviation vectors.
d1 = X(:,1) - mean_pt(1)*ones([3,1]);
d2 = X(:,2) - mean_pt(2)*ones([3,1]);
% Combine the data with the deviation vectors. First two rows are data,
% second two are the deviation vectors.
D = [X d1 d2]';
start = zeros(size(D));

% Plot the y_1 and y_2 vectors and the d_1, d_2 deviation vectors.
quiver3(start(:,1), start(:,2), start(:,3), D(:,1), D(:,2), D(:,3));
% Text for data.
for r = 1:height(D)
    if r < 3
        % Labels for the data, y_1 and y_2.
        anno = join(["$$\textbf{y}_",r," = \left[\begin{array}{c}",D(r,1),"\\",D(r,2),"\\",D(r,3),"\end{array}\right]$$"],'');
        text(D(r,1)-0.5,D(r,2)-0.2,D(r,3)-0.2,anno,'interpreter','latex');
    else
        % Labels for the deviation vectors, d_1 and d_2.
        anno = join(["$$\textbf{d}_",r-2," = \left[\begin{array}{c}",D(r,1),"\\",D(r,2),"\\",D(r,3),"\end{array}\right]$$"],'');
        text(D(r,1)-0.5,D(r,2)-0.2,D(r,3)-0.2,anno,'interpreter','latex');
    end
end

% Exercise 3.1 (c)
acosd(d1'*d2/(norm(d1)*norm(d2)))


% Exercise 3.2 (a)
X = [3 4; 6 -2; 3 1];
mean_pt = mean(X);
hold on
    % Plot data.
    scatter(X(:,1),X(:,2),'blue','filled')
    % Plot mean.
    plot(mean_pt(1),mean_pt(2),'d')
    % Text for mean.
    text(mean_pt(1)+0.5,mean_pt(2)+0.2, ...
        join(["$$\bar{x}=\left[\begin{array}{c}",mean_pt(1),"\\",mean_pt(2),"\end{array}\right]$$"],''), ...
        'interpreter','latex')
    % Text for data.
    for r = 1:height(X)
        anno = join(["$$x_",r," = \left[\begin{array}{c}",X(r,1),"\\",X(r,2),"\end{array}\right]$$"],'');
        text(X(r,1)-1.5,X(r,2)-0.5,anno,'interpreter','latex');
    end
    xlim([0 10])
    ylim([-5 5])
    saveas(gcf,'sol3.2a.png')
hold off

% Exercise 3.2 (b)
% Compute the deviation vectors.
d1 = X(:,1) - mean_pt(1)*ones([3,1]);
d2 = X(:,2) - mean_pt(2)*ones([3,1]);
% Combine the data with the deviation vectors. First two rows are data,
% second two are the deviation vectors.
D = [X d1 d2]';
start = zeros(size(D));

% Plot the y_1 and y_2 vectors and the d_1, d_2 deviation vectors.
quiver3(start(:,1), start(:,2), start(:,3), D(:,1), D(:,2), D(:,3));
% Text for data.
for r = 1:height(D)
    if r < 3
        % Labels for the data, y_1 and y_2.
        anno = join(["$$\textbf{y}_",r," = \left[\begin{array}{c}",D(r,1),"\\",D(r,2),"\\",D(r,3),"\end{array}\right]$$"],'');
        text(D(r,1)-0.5,D(r,2)-0.2,D(r,3)-0.2,anno,'interpreter','latex');
    else
        % Labels for the deviation vectors, d_1 and d_2.
        anno = join(["$$\textbf{d}_",r-2," = \left[\begin{array}{c}",D(r,1),"\\",D(r,2),"\\",D(r,3),"\end{array}\right]$$"],'');
        text(D(r,1)-0.5,D(r,2)-0.2,D(r,3)-0.2,anno,'interpreter','latex');
    end
end

% Exercise 3.2 (c)
acosd(d1'*d2/(norm(d1)*norm(d2)))

% Exercise 3.4 (a)
% Divide original units by 1,000,000 so we're in units of millions.
y1 = [3497900 2485475 1782875 1725450 1645575 1469800]'/1000000;
a1 = ones(height(y1),1);
% The projection of y1 onto a1.
x1bar1 = ((y1'*a1)/(norm(a1)*norm(a1)))*a1;

% Exercise 3.4 (b)
d1 = y1 - x1bar1;

% Exercise 3.4 (c)
clear clf
hold on
start = zeros(1,2);
quiver(start(:,1), start(:,2), 1, 1, 3);
quiver(start(:,1), start(:,2), 1, 0, 3);
quiver(3, 0, 0, 1, 3);
xlim([-1,5])
ylim([-2,4])
% Use the DRAWBRACE created by Pål Næverlid Sævik
drawbrace([3.2, -0], [3.2, -3],10,1,'Color','k') % Draws a curly brace for deviance vector.
% Text for the norm of deviance vector d1, norm(d1).
anno = join(["$$\left\|\textbf{d}_1\right\|"," = \left\|\textbf{y}_1 - \bar{x}_1 \textbf{1}_6 \right\| = 1.7167$$"],'');
text(3.2+0.4,1.5+1.5,anno,'Rotation',270,'interpreter','latex');

drawbrace([0, 0.2], [3, 0.2],10,1,'Color','k') % Draws a curly brace for projection vector
% % Text for the norm of the projection vector, norm(x1bar1).
anno = join(["$$\left\|\rm{Proj}_{\textbf{1}_6}\textbf{y}_1 \right\| = \left\|\bar{x}_1\textbf{1}_6\right\|"," = 5.1468$$"],'');
text(0.2,-0.7,anno,'interpreter','latex');

drawbrace([-0.1, 0.2], [2.9, 3.2],10,0,'Color','k') % Draws a curly brace for y vector.
% Text for the norm of y2, norm(y1).
anno = join(["$$\left\|\textbf{y}_1 \right\| = "," 5.4256$$"],'');
text(0.5,1.4,anno,'Rotation',40,'interpreter','latex');

% Include text for angle theta.
text(0.3,0.2,'\theta')

% Include the right-angle symbol on plot.
plot([3,2.7],[0.3,0.3])
plot([2.7,2.7],[0,0.3])
set(gca,'xtick',[],'ytick',[]);
hold off
saveas(gcf,'sol3.4c.png')

% Exercise 3.4 (d)
y2 = [0.623 0.593 0.512 0.500 0.463 0.395]';
% The projection of y2 onto a1.
x2bar1 = ((y2'*a1)/(norm(a1)*norm(a1)))*a1;

% Compute the deviance vector for y2.
d2 = y2 - x2bar1;

clear clf
hold on
start = zeros(1,2);
quiver(start(:,1), start(:,2), 1, 1, 3);
quiver(start(:,1), start(:,2), 1, 0, 3);
quiver(3, 0, 0, 1, 3);
xlim([-1,5])
ylim([-2,4])
% Use the DRAWBRACE created by Pål Næverlid Sævik
drawbrace([3.2, -0], [3.2, -3],10,1,'Color','k') % Draws a curly brace for deviance vector.
% Text for the norm of deviance vector d2, norm(d2).
anno = join(["$$\left\|\textbf{d}_2\right\|"," = \left\|\textbf{y}_2 - \bar{x}_2 \textbf{1}_6 \right\| = 0.1873$$"],'');
text(3.2+0.4,1.5+1.5,anno,'Rotation',270,'interpreter','latex');

drawbrace([0, 0.2], [3, 0.2],10,1,'Color','k') % Draws a curly brace for projection vector.
% Text for the norm of the projection vector, norm(x2bar1).
anno = join(["$$\left\|\rm{Proj}_{\textbf{1}_6}\textbf{y}_2 \right\| = \left\|\bar{x}_2\textbf{1}_6\right\|"," = 1.2599$$"],'');
text(0.2,-0.7,anno,'interpreter','latex');

drawbrace([-0.1, 0.2], [2.9, 3.2],10,0,'Color','k') % Draws a curly brace for y vector.
% Text for the norm of y2, norm(y2).
anno = join(["$$\left\|\textbf{y}_2 \right\| = "," 1.2737$$"],'');
text(0.5,1.4,anno,'Rotation',40,'interpreter','latex');

% Include text for angle theta.
text(0.3,0.2,'\theta')

% Include the right-angle symbol on plot.
plot([3,2.7],[0.3,0.3])
plot([2.7,2.7],[0,0.3])
set(gca,'xtick',[],'ytick',[]);
hold off
saveas(gcf,'sol3.4d.png')

% Exercise 3.4 (e)
% Combine the deviation vectors.
D = [d1 d2]';
start = zeros(size(D));

clear clf
% Plot the d_1, d_2 deviation vectors.
quiver3(start(:,1), start(:,2), start(:,3), D(:,1), D(:,2), D(:,3));
% Text for data.
for r = 1:height(D)
    % Labels for the deviation vectors, d_1 and d_2.
    anno = join(["$$\textbf{d}_",r," = \left[\begin{array}{c}",D(r,1),"\\",D(r,2),"\\",D(r,3),"\end{array}\right]$$"],'');
    text(D(r,1)-0.05,D(r,2)-0.05,D(r,3)-0.05,anno,'interpreter','latex');
end

% Compute the angle between d1 and d2.
acosd(d1'*d2/(norm(d1)*norm(d2)))
% Same answer using (2-36) on page 72, R = V^{-1/2} \Sigma V^{-1/2}.
acosd(diag(sqrt(inv(diag(diag(cov(y1,y2)))))*cov(y1,y2)*sqrt(inv(diag(diag(cov(y1,y2))))),1))

% Exercise 3.5 (a)
X = [9 1; 5 3; 1 2];
det(cov(X))

% Exercise 3.5 (b)
X = [3 4; 6 -2; 3 1];
det(cov(X))

% Exercise 3.6 (a)
X = [-1 3 -2; 2 4 2; 5 2 3];
xbar = mean(X)';
D = X - ones(height(X),1)*xbar';

% Exercise 3.6 (b)
covX = D'*D*(1/2);

% Exercise 3.6 (c)
trace(covX)

% Exercise 3.7
% Assuming that the mean vector is 0 for these 3 matrices.
xbar = [0;0];
S = [5 4; 4 5];
[V,D] = eig(S);
c = 1;
MyPlotEllipse(V,xbar,c/sqrt(D(1,1)),c/sqrt(D(2,2)),[-c c],[-c c],output_path,'sol3.7.1')

S = [5 -4; -4 5];
[V,D] = eig(S);
c = 1;
MyPlotEllipse(V,xbar,c/sqrt(D(1,1)),c/sqrt(D(2,2)),[-c c],[-c c],output_path,'sol3.7.2')

S = [3 0; 0 3];
[V,D] = eig(S);
c = 1;
MyPlotEllipse(V,xbar,c/sqrt(D(1,1)),c/sqrt(D(2,2)),[-c c],[-c c],output_path,'sol3.7.3')

% Exercise 3.8 (a)
S1 = eye(3);
S2 = 1.5*eye(3)-0.5*ones(3,3);

trace(S1)
trace(S2)

% Exercise 3.8 (b)
det(S1)
det(S2)

% Exercise 3.9 (a)
X = [12 17 29;
     18 20 38;
     14 16 30;
     20 18 38;
     16 19 35];
D = X - mean(X);
a = [1 1 -1]';
all(X*a == zeros(height(X),1))
all(D*a == zeros(height(X),1))

% Exercise 3.9 (b)
% S computed 3 different ways.
S1 = (1/(height(X)-1))*D'*D;
S2 = (1/(height(X)-1))*X'*(eye(height(X)) - (1/height(X))*ones(height(X),height(X)))*X;
S3 = (1/(height(X)-1))*X'*D;
% Make sure all computations of S are equal.
isequal(round(S1,2),round(S2,2))
isequal(S1,S3)
isequal(round(S2,2),round(S3,2))

det(S1)
S1*a

% Exercise 3.10 (a)
X = [3 1 0;
     6 4 6;
     4 2 2;
     7 0 3;
     5 3 4];
rank(X)
D = X - mean(X);
all(D*a == zeros(height(X),1))

% Exercise 3.10 (b)
S1 = (1/(height(X)-1))*D'*D;
S2 = (1/(height(X)-1))*X'*(eye(height(X)) - (1/height(X))*ones(height(X),height(X)))*X;
S3 = (1/(height(X)-1))*X'*D;
% Make sure all computations of S are equal.
isequal(round(S1,2),round(S2,2))
isequal(S1,S3)
isequal(round(S2,2),round(S3,2))

det(S1)

% Exercise 3.11 (a)
S = [252.04 -68.43; -68.43 123.67];
Dn12 = sqrt(inv(diag(diag(S))));
R = Dn12*S*Dn12
D12 = sqrt(diag(diag(S)));
D12*R*D12

% Exercise 3.14 (a)
X = [9 5 1; 1 3 2]';

c = [-1 2]';
Y1 = X*c;
EY1 = mean(Y1);
SY1 = (1/(height(Y1)-1))*(Y1 - ones(height(Y1),1)*EY1)'*(Y1 - ones(height(Y1),1)*EY1);

b = [2 3]';
Y2 = X*b;
EY2 = mean(Y2);
SY2 = (1/(height(Y2)-1))*(Y2 - ones(height(Y2),1)*EY2)'*(Y2 - ones(height(Y2),1)*EY2);

SY12 = (1/(height(Y1)-1))*(Y1 - ones(height(Y1),1)*EY1)'*(Y2 - ones(height(Y2),1)*EY2);

% Exercise 3.14 (b)
EX = (1/height(X))*ones(1,height(X))*X;
SX = (1/(height(X)-1))*(X - ones(height(X),1)*EX)'*(X - ones(height(X),1)*EX);

EcX = EX*c;
ScX = c'*SX*c;
EY1 == round(EcX)
SY1 == ScX

EbX = EX*b;
SbX = b'*SX*b;
EY2 == EbX
SY2 == SbX

ScXb = c'*SX*b;
SY12 == ScXb

% Exercise 3.15 (a)
X = [1 6 8; 4 2 3; 3 6 3]';

b = [1 1 1]';
Y1 = X*b;
EY1 = mean(Y1);
SY1 = (1/(height(Y1)-1))*(Y1 - ones(height(Y1),1)*EY1)'*(Y1 - ones(height(Y1),1)*EY1);

c = [1 2 -3]';
Y2 = X*c;
EY2 = mean(Y2);
SY2 = (1/(height(Y2)-1))*(Y2 - ones(height(Y2),1)*EY2)'*(Y2 - ones(height(Y2),1)*EY2);

SY12 = (1/(height(Y1)-1))*(Y1 - ones(height(Y1),1)*EY1)'*(Y2 - ones(height(Y2),1)*EY2);

% Exercise 3.15 (b)
EX = (1/height(X))*ones(1,height(X))*X;
SX = (1/(height(X)-1))*(X - ones(height(X),1)*EX)'*(X - ones(height(X),1)*EX);

EbX = EX*b;
SbX = b'*SX*b;
EY1 == EbX
SY1 == SbX

EcX = EX*c;
ScX = c'*SX*c;
EY2 == round(EcX)
SY2 == ScX

SbXc = b'*SX*c;
SY12 == SbXc

% Exercise 3.18 (a)
xbar = [0.766 0.508 0.438 0.161]';
S = [0.856 0.635 0.173 0.096;
     0.635 0.568 0.128 0.067;
     0.173 0.127 0.171 0.039;
     0.096 0.067 0.039 0.043];
b = [1 1 1 1]';
b'*xbar
b'*S*b

% Exercise 3.18 (b)
c = [1 -1 0 0]';
c'*xbar
c'*S*c
c'*S*b

% Exercise 3.19
S = S(1:3,1:3)
det(S)
D12 = inv(sqrt(diag(diag(S))));
R = D12*S*D12;
prod(diag(S))* det(R)
round(det(S),10) == round(prod(diag(S))* det(R),10)

% Exercise 3.20 (a)
X = [12.5 13.7;
    14.5 16.5;
    8.0 17.4;
    9.0 11.0;
    19.5 23.6;
    8.0 13.2;
    9.0 32.1;
    7.0 12.3;
    7.0 11.8;
    9.0 24.4;
    6.5 18.2;
    10.5 22.0;
    10.0 32.5;
    4.5 18.7;
    7.0 15.8;
    8.5 15.6;
    6.5 12.0;
    8.0 12.8;
    3.5 26.1;
    8.0 14.5;
    17.5 42.3;
    10.5 17.5;
    12.0 21.8;
    6.0 10.4;
    13.0 25.6];
Xbar = mean(X)';
S = (1/(height(X)-1))*(X - ones(height(X),1)*Xbar')'*(X - ones(height(X),1)*Xbar');
b = [1 -1]';

bXbar = b'*Xbar;
bS = b'*S*b;

% Exercise 3.20 (b)
Xt = X*b;
Xtbar = mean(Xt)';
XtS = (1/(height(Xt)-1))*(Xt - ones(height(Xt),1)*Xtbar')'*(Xt - ones(height(Xt),1)*Xtbar');
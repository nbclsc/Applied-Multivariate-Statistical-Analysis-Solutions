% Exercise 4.2 (c)
mu = [0;2];
Sigma = [2 sqrt(2)/2; sqrt(2)/2 1];
[V,D] = eig(Sigma);
c = sqrt(chi2inv(0.5,2));
MyPlotEllipse(V,mu,c*sqrt(D(1,1)),c*sqrt(D(2,2)),[-2 2],[0 4],output_path,'sol4.2')
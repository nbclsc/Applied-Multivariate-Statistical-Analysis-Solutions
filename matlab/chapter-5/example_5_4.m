closed = readmatrix(fullfile(data_folder, 'Table4.1.xlsx'));
open = readmatrix(fullfile(data_folder, 'Table4.5.xlsx'));
X = [closed(:,2).^0.25 open(:,2).^0.25];

n = height(X);
p = width(X);

xbar = mean(X)';
S = cov(X);

a_mu = [0.562 0.589]';
an_alpha = 0.05;

% Comput the test statistic for the H0 value of mu.
test_value = n*( xbar - a_mu)'/S*(xbar - a_mu);

% Compute the critical value from the reference distribution.
critical_value = (n-1)*p/(n-p)*icdf('F',1-an_alpha,p,n-p);

% Returns 1, so test value is within the confidence region.
test_value <= critical_valu

[e, lmbda] = eig(S);

% Notice the n in the denominator.
xleng = sqrt(lmbda(1,1))*sqrt(((n-1)*p)/(n*(n-p))*icdf('F',1-an_alpha,p,n-p));
yleng = sqrt(lmbda(2,2))*sqrt(((n-1)*p)/(n*(n-p))*icdf('F',1-an_alpha,p,n-p));

ci95 = xbar + [-1 1].*sqrt(diag(S))*sqrt(((n-1)*p)/(n*(n-p))*icdf('F',1-an_alpha,p,n-p));


% Plot the ellipse seen in figure 5.2.
my_plot_ellipse(e, xbar, xleng, yleng, [0.51 0.62] , [0.5 0.68])
line([ci95(1,1), ci95(1,1)], [0, ci95(2,2)], 'Color', 'r', 'LineStyle', '--'); % Vertical dashed line for lower CI for x1
line([ci95(1,2), ci95(1,2)], [0, ci95(2,2)], 'Color', 'r', 'LineStyle', '--'); % Vertical dashed line for upper CI for x1

line([0, ci95(1,2)], [ci95(2,1), ci95(2,1)], 'Color', 'r', 'LineStyle', '--'); % Horizontal dashed line for lower CI for x2
line([0, ci95(1,2)], [ci95(2,2), ci95(2,2)], 'Color', 'r', 'LineStyle', '--'); % Horizontal dashed line for upper CI for x2

% Label the x1 (closed) CI projections.
text(ci95(1,1), ci95(2,1)-0.03, '0.516', 'Rotation', 90)
text(ci95(1,2), ci95(2,1)-0.03, '0.612', 'Rotation', 90)

% Label the x2 (open) CI projections.
text(ci95(1,1)+0.01, ci95(2,2), '0.651')
text(ci95(1,1)+0.01, ci95(2,1), '0.555')

title('Simultaneous T^2-intervals for component means as shadows', 'of the confidence ellipse  on the axis - Microwave Data')
xlabel('x_{1}: Open Radiation')
ylabel('x_{2}: Closed Radiation')
saveas(gcf, fullfile('.', 'example5.4.png'), 'png')

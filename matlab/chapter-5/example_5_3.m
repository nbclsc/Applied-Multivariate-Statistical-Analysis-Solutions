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

% Plot the ellipse seen in figure 5.1.
my_plot_ellipse(e, xbar, xleng, yleng, [0.51 0.62] , [0.5 0.68])
line([xbar(1), xbar(1)], [0, xbar(2)], 'Color', 'r', 'LineStyle', '--'); % Vertical dashed line
line([0.51, xbar(1)], [xbar(2), xbar(2)], 'Color', 'r', 'LineStyle', '--'); % Horizontal dashed line
text(xbar(1), 0.50+0.005, '$\bar{x}_{1}$', 'Interpreter', 'latex')
text(0.51+0.003, xbar(2), '$\bar{x}_{2}$', 'Interpreter', 'latex')
title('95 % CI for \mu using Microwave Data')
xlabel('x_{1}')
ylabel('x_{2}')
saveas(gcf, fullfile('.', 'example5.3.png'), 'png')

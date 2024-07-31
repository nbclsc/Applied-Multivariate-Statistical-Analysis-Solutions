clep_cqt = readmatrix(fullfile(data_folder, 'Table5.2.xlsx'));
X = clep_cqt(:, 2:end);

n = height(X);
p = width(X);

xbar = mean(X)';
S = cov(X);

an_alpha = 0.05;
f_val = icdf('F', 1-an_alpha, p, n-p);

const = (((n - 1)*p)/(n - p))*f_val;

% Compute the 95% CI.
ci95 = xbar + [-1 1].*sqrt(const*(diag(S))/n);

% Using a linear combination. Compute the 95% CI.
a = [0 1 -1]';
a'*xbar + [-1 1].*sqrt(const*(a'*S*a/n))

X12 = X(:, 1:2);
X13 = X(:, [1 3]);
X23 = X(:, 2:3);

xbar12 = mean(X12)';
S12 = cov(X12);

xbar13 = mean(X13)';
S13 = cov(X13);

xbar23 = mean(X23)';
S23 = cov(X23);

[e12, lmbda12] = eig(S12);
[e13, lmbda13] = eig(S13);
[e23, lmbda23] = eig(S23);


xleng12 = sqrt(lmbda12(1,1))*sqrt(const/n);
yleng12 = sqrt(lmbda12(2,2))*sqrt(const/n);

xleng13 = sqrt(lmbda13(1,1))*sqrt(const/n);
yleng13 = sqrt(lmbda13(2,2))*sqrt(const/n);

xleng23 = sqrt(lmbda23(1,1))*sqrt(const/n);
yleng23 = sqrt(lmbda23(2,2))*sqrt(const/n);

% Plot the ellipse seen in figure 5.3.

% X1 and X2
my_plot_ellipse(e12, xbar12, xleng12, yleng12, [500 560] , [50 60])
line([ci95(1,1), ci95(1,1)], [0, ci95(2,2)], 'Color', 'r', 'LineStyle', '--'); % Vertical dashed line for lower CI for x1
line([ci95(1,2), ci95(1,2)], [0, ci95(2,2)], 'Color', 'r', 'LineStyle', '--'); % Vertical dashed line for upper CI for x1

line([0, ci95(1,2)], [ci95(2,1), ci95(2,1)], 'Color', 'r', 'LineStyle', '--'); % Horizontal dashed line for lower CI for x2
line([0, ci95(1,2)], [ci95(2,2), ci95(2,2)], 'Color', 'r', 'LineStyle', '--'); % Horizontal dashed line for upper CI for x2

% Label the x1 (closed) CI projections.
text(ci95(1,1), ci95(2,1)-0.10, sprintf('%6.2f', ci95(1,1)), 'Rotation', 90)
text(ci95(1,2), ci95(2,1)-0.10, sprintf('%6.2f', ci95(1,2)), 'Rotation', 90)

% Label the x2 (open) CI projections.
text(ci95(1,1)+10, ci95(2,2), sprintf('%6.2f', ci95(2,1)))
text(ci95(1,1)+10, ci95(2,1), sprintf('%6.2f', ci95(2,2)))

title('Simultaneous T^2-intervals for component means as shadows', 'of the confidence ellipse  on the axis: \mu_{1} & \mu_{2}')
xlabel('\mu_{1}: Social science and history')
ylabel('\mu_{2}: Verbal')
saveas(gcf, fullfile('.', 'example5.5a.png'), 'png')

% X1 and X3
my_plot_ellipse(e13, xbar13, xleng13, yleng13, [500 560] , [22 28])
line([ci95(1,1), ci95(1,1)], [0, ci95(3,2)], 'Color', 'r', 'LineStyle', '--'); % Vertical dashed line for lower CI for x1
line([ci95(1,2), ci95(1,2)], [0, ci95(3,2)], 'Color', 'r', 'LineStyle', '--'); % Vertical dashed line for upper CI for x1

line([0, ci95(1,2)], [ci95(3,1), ci95(3,1)], 'Color', 'r', 'LineStyle', '--'); % Horizontal dashed line for lower CI for x2
line([0, ci95(1,2)], [ci95(3,2), ci95(3,2)], 'Color', 'r', 'LineStyle', '--'); % Horizontal dashed line for upper CI for x2

% Label the x1 (closed) CI projections.
text(ci95(1,1), ci95(3,1)-1, sprintf('%6.2f', ci95(1,1)), 'Rotation', 90)
text(ci95(1,2), ci95(3,1)-1, sprintf('%6.2f', ci95(1,2)), 'Rotation', 90)

% Label the x2 (open) CI projections.
text(ci95(1,1)+10, ci95(3,2), sprintf('%6.2f', ci95(3,1)))
text(ci95(1,1)+10, ci95(3,1), sprintf('%6.2f', ci95(3,2)))

title('Simultaneous T^2-intervals for component means as shadows', 'of the confidence ellipse  on the axis: \mu_{1} & \mu_{3}')
xlabel('\mu_{1}: Social science and history')
ylabel('\mu_{3}: Science')
saveas(gcf, fullfile('.', 'example5.5b.png'), 'png')


% X2 and X3
my_plot_ellipse(e23, xbar23, xleng23, yleng23, [50 60] , [22 28])
line([ci95(2,1), ci95(2,1)], [0, ci95(3,2)], 'Color', 'r', 'LineStyle', '--'); % Vertical dashed line for lower CI for x1
line([ci95(2,2), ci95(2,2)], [0, ci95(3,2)], 'Color', 'r', 'LineStyle', '--'); % Vertical dashed line for upper CI for x1

line([0, ci95(2,2)], [ci95(3,1), ci95(3,1)], 'Color', 'r', 'LineStyle', '--'); % Horizontal dashed line for lower CI for x2
line([0, ci95(2,2)], [ci95(3,2), ci95(3,2)], 'Color', 'r', 'LineStyle', '--'); % Horizontal dashed line for upper CI for x2

% Label the x1 (closed) CI projections.
text(ci95(2,1), ci95(3,1)-1, sprintf('%6.2f', ci95(2,1)), 'Rotation', 90)
text(ci95(2,2), ci95(3,1)-1, sprintf('%6.2f', ci95(2,2)), 'Rotation', 90)

% Label the x2 (open) CI projections.
text(ci95(1,1)+10, ci95(3,2), sprintf('%6.2f', ci95(3,1)))
text(ci95(1,1)+10, ci95(3,1), sprintf('%6.2f', ci95(3,2)))

title('Simultaneous T^2-intervals for component means as shadows', 'of the confidence ellipse  on the axis: \mu_{2} & \mu_{3}')
xlabel('\mu_{2}: Verbal')
ylabel('\mu_{3}: Science')
saveas(gcf, fullfile('.', 'example5.5c.png'), 'png')



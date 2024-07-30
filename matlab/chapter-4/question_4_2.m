% Exercise 4.2 (c)
mu = [0; 2];
Sigma = [2 sqrt(2)/2; sqrt(2)/2 1];
[V,D] = eig(Sigma);
c = sqrt(chi2inv(0.5,2));

my_plot_ellipse(V,mu,c*sqrt(D(1,1)),c*sqrt(D(2,2)),[-2 2],[0 4])
line([mu(1), mu(1)], [0, mu(2)], 'Color', 'r', 'LineStyle', '--'); % Vertical dashed line
line([-2, mu(1)], [mu(2), mu(2)], 'Color', 'r', 'LineStyle', '--'); % Horizontal dashed line
text(mu(1), 0+0.10, '$\bar{x}_{1}$', 'Interpreter', 'latex')
text(-2+0.10, mu(2), '$\bar{x}_{2}$', 'Interpreter', 'latex')
title('50% Probability Contour')
xlabel('x_{1}')
ylabel('x_{2}')
saveas(gcf, fullfile('.','sol4.2.png'), 'png')
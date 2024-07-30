% Exercise 4.23 (a)
x = sort([-0.6 3.1 25.3 -16.8 -7.1 -6.2 25.2 22.6 26.0]);
% x = sort([-1 -.1 .16 .41 .62 .8 1.26 1.54 1.71 2.3])
prob_x = (linspace(1, length(x), length(x)) - 0.5) / length(x)
quant_x = icdf('Normal', prob_x, 0, 1)

p = scatter(quant_x, x)
title('Q-Q Plot')
xlabel('q_{(j)}')
ylabel('x_{(j)}')
saveas(p, append('.\', 'sol4.23a', '.png'), 'png')

% Exercise 4.23 (b)
% Using the sample covariance matrix:
xq_cov = cov([x', quant_x'])
xq_cov(1,2)/sqrt(xq_cov(1,1)*xq_cov(2,2))
% Compute "by-hand", using the vectors. Note that mean(quant_x) = 0
(x - mean(x))*(quant_x - mean(quant_x))'/sqrt((x - mean(x))*(x - mean(x))'*(quant_x - mean(quant_x))*(quant_x - mean(quant_x))')

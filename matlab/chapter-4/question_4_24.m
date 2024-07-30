% Exercise 4.24 (a)
data_4_24 = readtable(fullfile(data_folder, 'Exercise1.4.xlsx'));
sales = data_4_24.X1_sales;
x1 = sort(sales);
prob_x1 = (linspace(1, length(x1), length(x1)) - 0.50) / length(x1);
quant_x1 = icdf('Normal',prob_x1, 0, 1)';

p =scatter(quant_x1, x1)
title('Q-Q Plot Sales')
xlabel('q_{(j)}')
ylabel('x_{(j)}')
saveas(p, append('.\', 'sol4.24a_sales', '.png'), 'png')

profits = data_4_24.Xz_profits;
x2 = sort(profits);
prob_x2 = (linspace(1, length(x2), length(x2)) - 0.50) / length(x2);
quant_x2 = icdf('Normal', prob_x2, 0, 1)';  % Transpose to make the output a column vector.

p = scatter(quant_x2, x2)
title('Q-Q Plot Profits')
xlabel('q_{(j)}')
ylabel('x_{(j)}')
saveas(p, append('.\', 'sol4.24a_profits', '.png'), 'png')

% Exercise 4.24 (b)

x1_quant_x1_ss = (x1 - mean(x1))'*(quant_x1 - mean(quant_x1));
x1_ss = ((x1 - mean(x1))'*(x1 - mean(x1)));
quant_x1_ss = ((quant_x1 - mean(quant_x1))'*(quant_x1 - mean(quant_x1)));
sales_rq = x1_quant_x1_ss / sqrt(x1_ss * quant_x1_ss);

x2_quant_x2_ss = (x2 - mean(x2))'*(quant_x2 - mean(quant_x2));
x2_ss = ((x2 - mean(x2))'*(x2 - mean(x2)));
quant_x2_ss = ((quant_x2 - mean(quant_x2))'*(quant_x2 - mean(quant_x2)));
profits_rq = x2_quant_x2_ss / sqrt(x2_ss * quant_x2_ss);

% For profits (x2), profits_rq is the same answer as
% corr(x2, quant_x2)
% or same as off-diagonal of
% cov(x2, quant_x2)/sqrt(var(x2)*var(quant_x2))
% or same as off-diagonal of
% corrcoef(x2, quant_x2)

[samp, crit] = ppcc_simulation(length(x2), 10000000, 0.1);
crit < sales_rq
crit < profits_rq
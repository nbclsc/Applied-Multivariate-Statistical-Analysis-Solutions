% Exercise 4.25
% To be consistent with the book each observation is a column vector.
data_4_25 = readtable(fullfile(data_folder, 'Exercise1.4.xlsx'));
data_4_25 = [data_4_25.X1_sales data_4_25.Xz_profits data_4_25.x3_assets]';
x_bar = mean(data_4_25')';
S = cov(data_4_25');
d2 = sort(sum((data_4_25 - x_bar)'/S.*(data_4_25 - x_bar)', 2));
% Another way, sort(diag((data_4_25 - x_bar)'/S*(data_4_25 - x_bar)));
prob = (linspace(1, length(d2), length(d2)) - 0.50) / length(d2);
quant_d2 = icdf('chi2', prob, height(data_4_25));  % Same as given in book.

p = scatter(quant_d2, d2);
title('\chi^2 Plot Company data (Exercise 1.4)')
xlabel('$q_{c, 3}(\frac{j-0.5}{10})$','Interpreter','latex')
ylabel('d_{(j)}^{2}')
saveas(p, append('.\', 'sol4.25', '.png'), 'png')

[linspace(1,10,10)' d2 quant_d2']

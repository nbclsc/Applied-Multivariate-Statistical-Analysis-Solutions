% Exercise 4.26 (a)
% To be consistent with the book each observation is a column vector.
data_4_26 = readmatrix(fullfile(data_folder, 'Exercise1.2.xlsx'))';
x_bar = mean(data_4_26')';
S = cov(data_4_26');
d2 = sum((data_4_26 - x_bar)'/S.*(data_4_26 - x_bar)', 2);
% Another way, diag((data_4_26 - x_bar)'/S*(data_4_26 - x_bar));

% Exercise 4.26 (b)
sum(d2 < icdf('chi2',0.50,2)) / width(data_4_26)

% Exercise 4.26 (c)
d2 = sort(d2);
prob = (linspace(1, length(d2), length(d2)) - 0.5) / length(d2);
quant_d2 = icdf('chi2', prob, height(data_4_26));

p = scatter(quant_d2, d2)
title('\chi^2 Plot Used Car Selling Price')
xlabel('$q_{c, 2}(\frac{j - 0.5}{10})$', 'Interpreter', 'latex')
ylabel('d_{(j)}^{2}')
saveas(p, append('.\', 'sol4.26c', '.png'), 'png')

[linspace(1,10,10)' d2 quant_d2']

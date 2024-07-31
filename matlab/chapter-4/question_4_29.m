% Exercise 4.29 (a)
data_4_29 = readmatrix(fullfile(data_folder, 'Table1.5.xlsx'));
X = data_4_29(:,5:6)';
x_bar = mean(X')';
S = cov(X');
d2 = sum(((X - x_bar)'/S).*(X - x_bar)', 2);

d2 = ((X - x_bar)'/S*(X - x_bar);
% For d2, same answer as diag((X - x_bar)'/S*(X - x_bar))

% Exercise 4.29 (b)
sum(d2 <= icdf('chi2', 0.50, 2)) / length(d2)

% Exercise 4.29 (c)
d2 = sort(d2);
prob = (linspace(1, length(d2), length(d2)) - 0.5) / length(d2);
quant_d2 = icdf('chi2', prob, height(x_bar));

p = scatter(quant_d2, d2);
title('\chi^2 Plot NO_2 and O_3')
xlabel('$q_{c,2}(\frac{j - 0.5}{42})$', 'Interpreter', 'latex')
ylabel('d_{(j)}^{2}')
saveas(p, append('.\', 'sol4.29c', '.png'), 'png')

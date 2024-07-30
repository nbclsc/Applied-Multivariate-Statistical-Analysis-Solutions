% Exercise 4.28
data_4_28 = readmatrix(fullfile(data_folder, 'Table1.5.xlsx'));
x1 = sort(data_4_28(:,2));
prob = (linspace(1, length(x1), length(x1)) - 0.50) / length(x1);
quant_x1 = icdf('Normal', prob, 0, 1);

p = scatter(quant_x1, x1);
title('Q-Q Plot Solar Radiation')
xlabel('q_{(j)}')
ylabel('x_{(j)}')
saveas(p, append('.\', 'sol4.28', '.png'), 'png')

corr(quant_x1',x1)  % 0.9693

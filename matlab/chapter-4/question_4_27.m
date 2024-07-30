% Exercise 4.27
data_4_27 = readmatrix(fullfile(data_folder, 'Table4.1.xlsx'));
x1 = sort(log(data_4_27(:,2)));
% x1 = sort(data_4_27(:,2).^(0.25));
prob = (linspace(1, length(x1), length(x1)) - 0.5) / length(x1);
quant_x1 = icdf('Normal', prob, 0, 1);

p = scatter(quant_x1, x1);
title('Q-Q Plot Oven Data Log Radiation')
xlabel('q_{(j)}')
ylabel('x_{(j)}')
saveas(p, append('.\', 'sol4.27', '.png'), 'png')

corr(x1, quant_x1')  % 0.9706

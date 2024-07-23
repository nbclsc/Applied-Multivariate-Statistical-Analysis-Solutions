% Exercise 4.39
% Psychological Profile Data
data_4_39 = readmatrix(fullfile(data_folder, 'Table4.6.xlsx'));
X = data_4_39(:,1:5);

% Exercise 4.39 (a)

% Compute simulation of critical values for Q-Q correlation test.
n = height(X);
[~, crit_4_39] = ppcc_simulation(n, 1000000, [0.01 0.05 0.10]);

% x1: Independence Score.
x1 = X(:,1);
x1_col = "Independence Score";

% This returns 1, so YrHgt considered not normally distributed.
abs(sum((x1 >= (mean(x1) - 1*var(x1))) & (x1 <= (mean(x1) + 1*var(x1)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x1, r_Q_x1] = corr_q_q(x1);

% 0 for 0.01, but not 0.05 and 0.10.
r_Q_x1 < crit_4_39

simple_qq_plot(qq_x1(:,1), qq_x1(:,3), append(x1_col, " (raw)"));
saveas(gcf, append(".\", "sol4.39.qq.1", ".png"))

clear x1 x1_col r_Q_x1 qq_x1

% x2: Support Score.
x2 = X(:,2);
x2_col = "Support Score";

% This returns 1, so YrHgt considered not normally distributed.
abs(sum((x2 >= (mean(x2) - 1*var(x2))) & (x2 <= (mean(x2) + 1*var(x2)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x2, r_Q_x2] = corr_q_q(x2);

% 0 for 0.01, but not 0.05 and 0.10.
r_Q_x2 < crit_4_39

simple_qq_plot(qq_x2(:,1), qq_x2(:,3), append(x2_col, " (raw)"));
saveas(gcf, append(".\", "sol4.39.qq.2", ".png"))

clear x2 x2_col r_Q_x2 qq_x2

% x3: Benevolence Score.
x3 = X(:,3);
x3_col = "Benevolence Score";

% This returns 1, so YrHgt considered not normally distributed.
abs(sum((x3 >= (mean(x3) - 1*var(x3))) & (x3 <= (mean(x3) + 1*var(x3)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x3, r_Q_x3] = corr_q_q(x3);

% 0 for all. Considered normal at 0.01, 0.05, and 0.10 levels.
r_Q_x3 < crit_4_39

simple_qq_plot(qq_x3(:,1), qq_x3(:,3), append(x3_col, " (raw)"));
saveas(gcf, append(".\", "sol4.39.qq.3", ".png"))

clear x3 x3_col r_Q_x3 qq_x3

% x4: Conformity Score.
x4 = X(:,4);
x4_col = "Conformity Score";

% This returns 1, so YrHgt considered not normally distributed.
abs(sum((x4 >= (mean(x4) - 1*var(x4))) & (x4 <= (mean(x4) + 1*var(x4)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x4, r_Q_x4] = corr_q_q(x4);

% 0 for all. Considered normal at 0.01, 0.05, and 0.10 levels.
r_Q_x4 < crit_4_39

simple_qq_plot(qq_x4(:,1), qq_x4(:,3), append(x4_col, " (raw)"));
saveas(gcf, append(".\", "sol4.39.qq.4", ".png"))

clear x4 x4_col r_Q_x4 qq_x4

% x5: Leadership Score.
x5 = X(:,5);
x5_col = "Leadership Score";

% This returns 1, so YrHgt considered not normally distributed.
abs(sum((x5 >= (mean(x5) - 1*var(x5))) & (x5 <= (mean(x5) + 1*var(x5)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x5, r_Q_x5] = corr_q_q(x5);

% 
r_Q_x5 < crit_4_39

simple_qq_plot(qq_x5(:,1), qq_x5(:,3), append(x5_col, " (raw)"));
saveas(gcf, append(".\", "sol4.39.qq.5", ".png"))

clear x5 x5_col r_Q_x5 qq_x5

% Exercise 4.39 (b)

% Check multivariate normal.
output = chi2_q_q(X);

scatter(output(:,3), output(:,1))
title('\chi^{2} Plot Psychological Profile Data (Raw)')
xlabel('$q_{c,7}\left(\frac{j - 0.5}{130}\right)$', 'Interpreter', 'latex')
ylabel('d_{(j)}^{2}')
saveas(gcf, append(".\", "sol4.39.b", ".png"))

% Exercise 4.39 (c)

% x1: Independence Score.
x1 = X(:,1);
x1_col = "Independence Score";

% Perform power transformation.
power_x1 = box_cox_power_transform(x1, x1_col, ...
    -5, 5, append(".\", "sol4.39.power.1", ".png"));

% Use power transformation suggested.
x1_tr = x1.^0.50;

[qq_x1_tr, r_Q_x1_tr] = corr_q_q(x1_tr);

% Transformation worked!
r_Q_x1_tr < crit_4_39

simple_qq_plot(qq_x1_tr(:,1), qq_x1_tr(:,3), append(x1_col, " (transformed)"));
saveas(gcf, append(".\", "sol4.39.qq.tr.1", ".png"))

sum(abs((x1 - mean(x1))/sqrt(var(x1))) > 3)
sum(abs((x1_tr - mean(x1_tr))/sqrt(var(x1_tr))) > 3)

clear x1 x1_tr x1_col power_x1 qq_x1_tr r_Q_x1_tr

% x2: Support Score.
x2 = X(:,2);
x2_col = "Support Score";

power_x2 = box_cox_power_transform(x2, x2_col, ...
    -5, 5, append(".\", "sol4.39.power.2", ".png"));

% Use power transformation suggested.
x2_tr = x2.^1.4;

[qq_x2_tr, r_Q_x2_tr] = corr_q_q(x2_tr);

% Transformation worked!
r_Q_x2_tr < crit_4_39

simple_qq_plot(qq_x2_tr(:,1), qq_x2_tr(:,3), append(x2_col, " (transformed)"));
saveas(gcf, append(".\", "sol4.39.qq.tr.2", ".png"))

sum(abs((x2 - mean(x2))/sqrt(var(x2))) > 3)
sum(abs((x2_tr - mean(x2_tr))/sqrt(var(x2_tr))) > 3)

clear x2 x2_tr x2_col power_x2 qq_x2_tr r_Q_x2_tr

% x5: Leadership Score.
x5 = X(:,5);
x5_col = "Leadership Score";

power_x5 = box_cox_power_transform(x5, x5_col, ...
    -5, 5, append(".\", "sol4.39.power.5", ".png"));

% Use power transformation suggested.
x5_tr = x5.^0.40;

[qq_x5_tr, r_Q_x5_tr] = corr_q_q(x5_tr);

% Transformation worked!
r_Q_x5_tr < crit_4_39

simple_qq_plot(qq_x5_tr(:,1), qq_x5_tr(:,3), append(x5_col, " (transformed)"));
saveas(gcf, append(".\", "sol4.39.qq.tr.5", ".png"))

sum(abs((x5 - mean(x5))/sqrt(var(x5))) > 3)
sum(abs((x5_tr - mean(x5_tr))/sqrt(var(x5_tr))) > 3)

clear x5 x5_tr x5_col power_x5 qq_x5_tr r_Q_x5_tr

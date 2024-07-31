% Exercise 4.32
% Radiotherapy Data
data_4_32 = readmatrix(fullfile(data_folder, 'Table1.7.xlsx'));
% Source: Data courtesy of Mrs. Annette Tealey, R.N. Values of x2 and x3
% less than 1.0 are due to errors in the data-collection process. Rows
% containing values of x2 and x3 less than 1.0 may be omitted.

% Compute simulation of critical values for Q-Q correlation test.
n = height(data_4_32);
[~, crit_4_32] = ppcc_simulation(n, 1000000, [0.01 0.05 0.10]);

% x1: Symptoms
x1 = data_4_32(:,1);
x1_col = "Symptoms";

% This returns 1, so symptoms considered not normally distributed.
abs(sum((x1 >= (mean(x1) - 1*var(x1))) & (x1 <= (mean(x1) + 1*var(x1)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x1, r_Q_x1] = corr_q_q(x1);

% Okay at the 0.01-level, but not 0.05 or 0.10.
r_Q_x1 < crit_4_32

simple_qq_plot(qq_x1(:,1), qq_x1(:,3), append(x1_col, " (raw)"));

% Perform power transformation.
power_x1 = box_cox_power_transform(x1, x1_col, ...
    -5, 5, append(".\", "sol4.32.power.1", ".png"));

% Use power transformation suggested.
x1_tr = x1.^0.5;

[qq_x1_tr, r_Q_x1_tr] = corr_q_q(x1_tr);

% Transformation worked!
r_Q_x1_tr < crit_4_32

simple_qq_plot(qq_x1_tr(:,1), qq_x1_tr(:,3), append(x1_col, " (transformed)"));

% Saved as sol4.32.qq.1.png
subplot(1,2,1);
simple_qq_plot(qq_x1(:,1), qq_x1(:,3), append(x1_col, " (raw)"))
subplot(1,2,2);
simple_qq_plot(qq_x1_tr(:,1), qq_x1_tr(:,3), append(x1_col, " (transformed)"));

clear x1 x1_tr x1_col r_Q_x1 r_Q_x1_tr power_x1 qq_x1 qq_x1_tr

% x2: Activity
x2 = data_4_32(data_4_32(:,2) >= 1.0,2);
n_x2 = length(x2);
x2_col = "Activity";
[~, crit_4_32_x2] = ppcc_simulation(n_x2, 1000000, [0.01 0.05 0.10]);

% This returns 1, so activity considered not normally distributed.
abs(sum((x2 >= (mean(x2) - 1*var(x2))) & (x2 <= (mean(x2) + 1*var(x2)))) / n_x2 - 0.683) > 1.396/sqrt(n_x2)

% Look at the QQ plot of the raw data.
[qq_x2, r_Q_x2] = corr_q_q(x2);

% True at all levels, so not normal at 0.01, 0.05, or 0.10 levels.
r_Q_x2 < crit_4_32_x2

% Can see tons of values at 1, so causing the problem.
simple_qq_plot(qq_x2(:,1), qq_x2(:,3), append(x2_col, " (raw)"));

% Perform power transformation.
power_x2 = box_cox_power_transform(x2, x2_col, ...
    -5, 5, append(".\", "sol4.32.power.2", ".png"));

% Use power transformation suggested.
x2_tr = x2.^-0.50;

[qq_x2_tr, r_Q_x2_tr] = corr_q_q(x2_tr);

% Transformation improved some, but so many values of 1 are causing issues.
r_Q_x2_tr < crit_4_32_x2

simple_qq_plot(qq_x2_tr(:,1), qq_x2_tr(:,3), append(x2_col, " (transformed)"));

% Saved as sol4.32.qq.2.png
subplot(1,2,1);
simple_qq_plot(qq_x2(:,1), qq_x2(:,3), append(x2_col, " (raw)"))
subplot(1,2,2);
simple_qq_plot(qq_x2_tr(:,1), qq_x2_tr(:,3), append(x2_col, " (transformed)"));

clear x2 x2_tr x2_col r_Q_x2 r_Q_x2_tr power_x2 qq_x2 qq_x2_tr n_x2 crit_4_32_x2

% x3: Sleep
x3 = data_4_32((data_4_32(:,3) >= 1.0),3);
n_x3 = length(x3);
x3_col = "Sleep";
[~, crit_4_32_x3] = ppcc_simulation(n_x3, 1000000, [0.01 0.05 0.10]);

% This returns 1, so sleep is not normally distributed.
abs(sum((x3 >= (mean(x3) - 1*var(x3))) & (x3 <= (mean(x3) + 1*var(x3)))) / n_x3 - 0.683) > 1.396/sqrt(n_x3)

% Look at the QQ plot of the raw data.
[qq_x3, r_Q_x3] = corr_q_q(x3);

% False at all levels, so normal at 0.01, 0.05, or 0.10 levels.
r_Q_x3 < crit_4_32_x3

% Save the qq plot on the raw data.
simple_qq_plot(qq_x3(:,1), qq_x3(:,3), append(x3_col, " (raw)"));
saveas(gcf, append(".\", "sol4.32.qq.3", ".png"))

% Perform power transformation, but don't really need it.
power_x3 = box_cox_power_transform(x3, x3_col, ...
    -5, 5, append(".\", "sol4.32.power.3", ".png"));

% Use power transformation suggested.
x3_tr = x3.^power_x3;

[qq_x3_tr, r_Q_x3_tr] = corr_q_q(x3_tr);

% Transformation improved some.
r_Q_x3_tr < crit_4_32_x3

simple_qq_plot(qq_x3_tr(:,1), qq_x3_tr(:,3), append(x3_col, " (transformed)"));

clear x3 x3_tr x3_col r_Q_x3 r_Q_x3_tr power_x3 qq_x3 qq_x3_tr n_x3 crit_4_32_x3

% x4: Eat
x4 = data_4_32(:,4);
x4_col = "Eat";

% This returns 1, so eat considered not normally distributed.
abs(sum((x4 >= (mean(x4) - 1*var(x4))) & (x4 <= (mean(x4) + 1*var(x4)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x4, r_Q_x4] = corr_q_q(x4);

% True at all levels, so not normal at 0.01, 0.05, or 0.10 levels.
r_Q_x4 < crit_4_32

% QQ plot of the raw data.
simple_qq_plot(qq_x4(:,1), qq_x4(:,3), append(x4_col, " (raw)"));

% Perform power transformation.
power_x4 = box_cox_power_transform(x4, x4_col, ...
    -5, 5, append(".\", "sol4.32.power.4", ".png"));

% Use power transformation suggested.
x4_tr = x4.^0.25;

[qq_x4_tr, r_Q_x4_tr] = corr_q_q(x4_tr);

% Transformation improved some, but not at 0.05 and 0.1 levels.
r_Q_x4_tr < crit_4_32

simple_qq_plot(qq_x4_tr(:,1), qq_x4_tr(:,3), append(x4_col, " (transformed)"));

% Saved as sol4.32.qq.4.png
subplot(1,2,1);
simple_qq_plot(qq_x4(:,1), qq_x4(:,3), append(x4_col, " (raw)"))
subplot(1,2,2);
simple_qq_plot(qq_x4_tr(:,1), qq_x4_tr(:,3), append(x4_col, " (transformed)"));

min((x4 - mean(x4))/sqrt(var(x4)))

clear x4 x4_tr x4_col r_Q_x4 r_Q_x4_tr power_x4 qq_x4 qq_x4_tr

% x5: Appetite
x5 = data_4_32(:,5);
x5_col = "Appetite";

% This returns 0, so appetite considered normally distributed.
abs(sum((x5 >= (mean(x5) - 1*var(x5))) & (x5 <= (mean(x5) + 1*var(x5)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x5, r_Q_x5] = corr_q_q(x5);

% True at all levels, so considered normal at 0.01, 0.05, or 0.10 levels.
r_Q_x5 < crit_4_32

% Save the qq plot on the raw data.
simple_qq_plot(qq_x5(:,1), qq_x5(:,3), append(x5_col, " (raw)"));
saveas(gcf, append(".\", "sol4.32.qq.5", ".png"))

% Perform power transformation.
power_x5 = box_cox_power_transform(x5, x5_col, ...
    -5, 5, append(".\", "sol4.32.power.5", ".png"));

% Use power transformation suggested.
x5_tr = x5.^power_x5;

[qq_x5_tr, r_Q_x5_tr] = corr_q_q(x5_tr);

% Transformation improved some, but don'treally need it.
r_Q_x5_tr < crit_4_32

simple_qq_plot(qq_x5_tr(:,1), qq_x5_tr(:,3), append(x5_col, " (transformed)"));

clear x5 x5_tr x5_col r_Q_x5 r_Q_x5_tr power_x5 qq_x5 qq_x5_tr

% x6: Skin Reaction
x6 = data_4_32(:,6);
x6_col = "Skin Reaction";

% This returns 0, so symptoms considered normally distributed.
abs(sum((x6 >= (mean(x6) - 1*var(x6))) & (x6 <= (mean(x6) + 1*var(x6)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x6, r_Q_x6] = corr_q_q(x6);

% True at all levels, so not normal at 0.01, 0.05, or 0.10 levels.
r_Q_x6 < crit_4_32

% Data is discrete, so I'll keep the raw data, no transformation.
simple_qq_plot(qq_x6(:,1), qq_x6(:,3), append(x6_col, " (raw)"));
saveas(gcf, append(".\", "sol4.32.qq.6", ".png"))

% Just out of curiosity. Perform power transformation.
power_x6 = box_cox_power_transform(x6, x6_col, ...
    -5, 5, append(".\", "sol4.32.power.6", ".png"));

% Use power transformation suggested.
x6_tr = x6.^power_x6;

[qq_x6_tr, r_Q_x6_tr] = corr_q_q(x6_tr);

% Transformation made things worse.
r_Q_x6_tr < crit_4_32

simple_qq_plot(qq_x6_tr(:,1), qq_x6_tr(:,3), append(x6_col, " (transformed)"));

clear x6 x6_tr x6_col r_Q_x6 r_Q_x6_tr power_x6 qq_x6 qq_x6_tr
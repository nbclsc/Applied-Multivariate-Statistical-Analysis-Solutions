% Exercise 4.33
% Stiffness Data
data_4_33 = readmatrix(fullfile(data_folder, 'Table4.3.xlsx'));

% Compute simulation of critical values for Q-Q correlation test.
n = height(data_4_33);
[~, crit_4_33] = ppcc_simulation(n, 1000000, [0.01 0.05 0.10]);

% x1: Sending shockwave down the board.
x1 = data_4_33(:,2);
x1_col = "Shockwave";

% This returns 1, so shockwave considered not normally distributed.
abs(sum((x1 >= (mean(x1) - 1*var(x1))) & (x1 <= (mean(x1) + 1*var(x1)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x1, r_Q_x1] = corr_q_q(x1);

% Okay at the 0.01-level, but not 0.05 or 0.10.
r_Q_x1 < crit_4_33

simple_qq_plot(qq_x1(:,1), qq_x1(:,3), append(x1_col, " (raw)"));

% Perform power transformation.
power_x1 = box_cox_power_transform(x1, x1_col, ...
    append(".\", "sol4.33.power.1", ".png"));

% Use power transformation suggested.
x1_tr = x1.^-0.50;

[qq_x1_tr, r_Q_x1_tr] = corr_q_q(x1_tr);

% Transformation worked!
r_Q_x1_tr < crit_4_33

simple_qq_plot(qq_x1_tr(:,1), qq_x1_tr(:,3), append(x1_col, " (transformed)"));

% Saved as sol4.33.qq.1.png
subplot(1,2,1);
simple_qq_plot(qq_x1(:,1), qq_x1(:,3), append(x1_col, " (raw)"))
subplot(1,2,2);
simple_qq_plot(qq_x1_tr(:,1), qq_x1_tr(:,3), append(x1_col, " (transformed)"));

clear x1 x1_tr x1_col r_Q_x1 r_Q_x1_tr power_x1 qq_x1 qq_x1_tr gcf

% x2: Sending shockwave down the board.
x2 = data_4_33(:,3);
x2_col = "Vibration";

% This returns 1, so vibration considered not normally distributed.
abs(sum((x2 >= (mean(x2) - 1*var(x2))) & (x2 <= (mean(x2) + 1*var(x2)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x2, r_Q_x2] = corr_q_q(x2);

% Okay at the 0.01-level, but not 0.05 or 0.10.
r_Q_x2 < crit_4_33

simple_qq_plot(qq_x2(:,1), qq_x2(:,3), append(x2_col, " (raw)"));

% Perform power transformation.
power_x2 = box_cox_power_transform(x2, x2_col, ...
    -5, 5, append(".\", "sol4.33.power.2", ".png"));

% Use power transformation suggested.
x2_tr = x2.^-0.75;

[qq_x2_tr, r_Q_x2_tr] = corr_q_q(x2_tr);

% Transformation worked!
r_Q_x2_tr < crit_4_33

simple_qq_plot(qq_x2_tr(:,1), qq_x2_tr(:,3), append(x2_col, " (transformed)"));

% Saved as sol4.33.qq.2.png
subplot(1,2,1);
simple_qq_plot(qq_x2(:,1), qq_x2(:,3), append(x2_col, " (raw)"))
subplot(1,2,2);
simple_qq_plot(qq_x2_tr(:,1), qq_x2_tr(:,3), append(x2_col, " (transformed)"));

clear x2 x2_tr x2_col r_Q_x2 r_Q_x2_tr power_x2 qq_x2 qq_x2_tr gcf

% x3: Static test 1.
x3 = data_4_33(:,4);
x3_col = "Static Test 1";

% This returns 1, so static test 1 considered not normally distributed.
abs(sum((x3 >= (mean(x3) - 1*var(x3))) & (x3 <= (mean(x3) + 1*var(x3)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x3, r_Q_x3] = corr_q_q(x3);

% Okay at the 0.01-level, but not 0.05 or 0.10.
r_Q_x3 < crit_4_33

simple_qq_plot(qq_x3(:,1), qq_x3(:,3), append(x3_col, " (raw)"));

% Perform power transformation.
power_x3 = box_cox_power_transform(x3, x3_col, ...
    -5, 5, append(".\", "sol4.33.power.3", ".png"));

% Use power transformation suggested.
x3_tr = x3.^-0.75;

[qq_x3_tr, r_Q_x3_tr] = corr_q_q(x3_tr);

% Transformation worked!
r_Q_x3_tr < crit_4_33

simple_qq_plot(qq_x3_tr(:,1), qq_x3_tr(:,3), append(x3_col, " (transformed)"));

% Saved as sol4.33.qq.3.png
subplot(1,2,1);
simple_qq_plot(qq_x3(:,1), qq_x3(:,3), append(x3_col, " (raw)"))
subplot(1,2,2);
simple_qq_plot(qq_x3_tr(:,1), qq_x3_tr(:,3), append(x3_col, " (transformed)"));

clear x3 x3_tr x3_col r_Q_x3 r_Q_x3_tr power_x3 qq_x3 qq_x3_tr gcf

% x4: Static test 2.
x4 = data_4_33(:,5);
x4_col = "Static Test 2";

% This returns 1, so static test 2 considered not normally distributed.
abs(sum((x4 >= (mean(x4) - 1*var(x4))) & (x4 <= (mean(x4) + 1*var(x4)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x4, r_Q_x4] = corr_q_q(x4);

% Okay at the 0.01-level, but not 0.05 or 0.10.
r_Q_x4 < crit_4_33

simple_qq_plot(qq_x4(:,1), qq_x4(:,3), append(x4_col, " (raw)"));
saveas(gcf, append(".\", "sol4.33.qq.4", ".png"))

% Perform power transformation.
power_x4 = box_cox_power_transform(x4, x4_col, ...
    -5, 5, append(".\", "sol4.33.power.4", ".png"));

% Use power transformation suggested.
x4_tr = x4.^power_x4;

[qq_x4_tr, r_Q_x4_tr] = corr_q_q(x4_tr);

% Transformation worked!
r_Q_x4_tr < crit_4_33

simple_qq_plot(qq_x4_tr(:,1), qq_x4_tr(:,3), append(x4_col, " (transformed)"));

clear x4 x4_tr x4_col r_Q_x4 r_Q_x4_tr power_x4 qq_x4 qq_x4_tr gcf

% Checking bivariate normality.
chi50 = icdf('Chi2',0.5, 2);

% x1 and x2
bivar_output = zeros(6,2);

X12 = [data_4_33(:,2) data_4_33(:,3)];
d2_12 = diag((X12 - mean(X12))*inv(cov(X12))*(X12 - mean(X12))');
bivar_output(1,1) = sum(d2_12 <= chi50)/n;
X12tr = [data_4_33(:,2).^-0.50 data_4_33(:,3).^-0.75];
d2_12tr = diag((X12tr - mean(X12tr))*inv(cov(X12tr))*(X12tr - mean(X12tr))');
bivar_output(1,2) = sum(d2_12tr <= chi50)/n;

% x1 and x3
X13 = [data_4_33(:,2) data_4_33(:,4)];
d2_13 = diag((X13 - mean(X13))*inv(cov(X13))*(X13 - mean(X13))');
bivar_output(2,1) = sum(d2_13 <= chi50)/n;
X13tr = [data_4_33(:,2).^-0.50 data_4_33(:,4).^-0.75];
d2_13tr = diag((X13tr - mean(X13tr))*inv(cov(X13tr))*(X13tr - mean(X13tr))');
bivar_output(2,2) = sum(d2_13tr <= chi50)/n;

% x1 and x4
X14 = [data_4_33(:,2) data_4_33(:,5)];
d2_14 = diag((X14 - mean(X14))*inv(cov(X14))*(X14 - mean(X14))');
bivar_output(3,1) = sum(d2_14 <= chi50)/n;
X14tr = [data_4_33(:,2).^-0.50 data_4_33(:,5)];
d2_14tr = diag((X14tr - mean(X14tr))*inv(cov(X14tr))*(X14tr - mean(X14tr))');
bivar_output(3,2) = sum(d2_14tr <= chi50)/n;

% x2 and x3
X23 = [data_4_33(:,3) data_4_33(:,4)];
d2_23 = diag((X23 - mean(X23))*inv(cov(X23))*(X23 - mean(X23))');
bivar_output(4,1) = sum(d2_23 <= chi50)/n;
X23tr = [data_4_33(:,3).^-0.75 data_4_33(:,4).^-0.75];
d2_23tr = diag((X23tr - mean(X23tr))*inv(cov(X23tr))*(X23tr - mean(X23tr))');
bivar_output(4,2) = sum(d2_23tr <= chi50)/n;

% x2 and x4
X24 = [data_4_33(:,3) data_4_33(:,5)];
d2_24 = diag((X24 - mean(X24))*inv(cov(X24))*(X24 - mean(X24))');
bivar_output(5,1) = sum(d2_24 <= chi50)/n;
X24tr = [data_4_33(:,3).^-0.75 data_4_33(:,5)];
d2_24tr = diag((X24tr - mean(X24tr))*inv(cov(X24tr))*(X24tr - mean(X24tr))');
bivar_output(5,2) = sum(d2_24tr <= chi50)/n;

% x3 and x4
X34 = [data_4_33(:,4) data_4_33(:,5)];
d2_34 = diag((X34 - mean(X34))*inv(cov(X34))*(X34 - mean(X34))');
bivar_output(6,1) = sum(d2_34 <= chi50)/n;
X34tr = [data_4_33(:,4).^-0.75 data_4_33(:,5)];
d2_34tr = diag((X34tr - mean(X34tr))*inv(cov(X34tr))*(X34tr - mean(X34tr))');
bivar_output(6,2) = sum(d2_34tr <= chi50)/n;
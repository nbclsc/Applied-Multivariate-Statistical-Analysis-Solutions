% Exercise 4.35
% Paper Quality Data
data_4_35 = readmatrix(fullfile(data_folder, 'Table1.2.xlsx'));
% Delete outlier and old paper obs mentioned on page 20.
data_4_35 = data_4_35(ismember(data_4_35(:,1), [16:21 25 34 38:41])==0,:);

% Compute simulation of critical values for Q-Q correlation test.
n = height(data_4_35);
[~, crit_4_35] = ppcc_simulation(n, 1000000, [0.01 0.05 0.10]);

% x1: Paper Density.
x1 = data_4_35(:,2);
x1_col = "Paper Density";

% This returns 1, so density considered not normally distributed.
abs(sum((x1 >= (mean(x1) - 1*var(x1))) & (x1 <= (mean(x1) + 1*var(x1)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x1, r_Q_x1] = corr_q_q(x1);

% Okay at all three levels.
r_Q_x1 < crit_4_35

simple_qq_plot(qq_x1(:,1), qq_x1(:,3), append(x1_col, " (raw)"));
saveas(gcf, append(".\", "sol4.35.qq.1", ".png"))

% Perform power transformation.
power_x1 = box_cox_power_transform(x1, x1_col, ...
    append(".\", "sol4.35.power.1", ".png"));

% Use power transformation suggested.
x1_tr = x1.^power_x1;

[qq_x1_tr, r_Q_x1_tr] = corr_q_q(x1_tr);

% Transformation only worked at the 0.01-level.
r_Q_x1_tr < crit_4_35

simple_qq_plot(qq_x1_tr(:,1), qq_x1_tr(:,3), append(x1_col, " (transformed)"));

sum(abs((x1 - mean(x1))/sqrt(var(x1))) > 3)
sum(abs((x1_tr - mean(x1_tr))/sqrt(var(x1_tr))) > 3)

clear x1 x1_tr x1_col r_Q_x1 r_Q_x1_tr power_x1 qq_x1 qq_x1_tr gcf

% x2: Paper Strength Machine Direction.
x2 = data_4_35(:,3);
x2_col = "Paper Strength Machine Direction";

% This returns 1, so machine direction considered not normally distributed.
abs(sum((x2 >= (mean(x2) - 1*var(x2))) & (x2 <= (mean(x2) + 1*var(x2)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x2, r_Q_x2] = corr_q_q(x2);

% Okay at all three levels.
r_Q_x2 < crit_4_35

simple_qq_plot(qq_x2(:,1), qq_x2(:,3), append(x2_col, " (raw)"));
saveas(gcf, append(".\", "sol4.35.qq.2", ".png"))

% Perform power transformation.
power_x2 = box_cox_power_transform(x2, x2_col, ...
    append(".\", "sol4.35.power.2", ".png"));

% Use power transformation suggested.
x2_tr = x2.^power_x2;

[qq_x2_tr, r_Q_x2_tr] = corr_q_q(x2_tr);

% Transformation worked!
r_Q_x2_tr < crit_4_35

simple_qq_plot(qq_x2_tr(:,1), qq_x2_tr(:,3), append(x2_col, " (transformed)"));

sum(abs((x2 - mean(x2))/sqrt(var(x2))) > 3)
sum(abs((x2_tr - mean(x2_tr))/sqrt(var(x2_tr))) > 3)

clear x2 x2_tr x2_col r_Q_x2 r_Q_x2_tr power_x2 qq_x2 qq_x2_tr gcf

% x3: Paper Strength Cross Direction.
x3 = data_4_35(:,4);
x3_col = "Paper Strength Cross Direction";

% This returns 1, so cross direction considered not normally distributed.
abs(sum((x3 >= (mean(x3) - 1*var(x3))) & (x3 <= (mean(x3) + 1*var(x3)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x3, r_Q_x3] = corr_q_q(x3);

% Okay at all 3 levels.
r_Q_x3 < crit_4_35

simple_qq_plot(qq_x3(:,1), qq_x3(:,3), append(x3_col, " (raw)"));
saveas(gcf, append(".\", "sol4.35.qq.3", ".png"))

% Perform power transformation.
power_x3 = box_cox_power_transform(x3, x3_col, ...
    append(".\", "sol4.35.power.3", ".png"));

% Use power transformation suggested.
x3_tr = x3.^power_x3;

[qq_x3_tr, r_Q_x3_tr] = corr_q_q(x3_tr);

% Transformation worked!
r_Q_x3_tr < crit_4_35

simple_qq_plot(qq_x3_tr(:,1), qq_x3_tr(:,3), append(x3_col, " (transformed)"));

sum(abs((x3 - mean(x3))/sqrt(var(x3))) > 3)
sum(abs((x3_tr - mean(x3_tr))/sqrt(var(x3_tr))) > 3)

clear x3 x3_tr x3_col r_Q_x3 r_Q_x3_tr power_x3 qq_x3 qq_x3_tr gcf

% Checking bivariate normality.
chi50 = icdf('Chi2',0.5, 2);
bivar_output = zeros(3,1);

% x1 and x2
X12 = [data_4_35(:,2) data_4_35(:,3)];
d2_12 = diag((X12 - mean(X12))*inv(cov(X12))*(X12 - mean(X12))');
bivar_output(1,1) = sum(d2_12 <= chi50)/n;

% x1 and x3
X13 = [data_4_35(:,2) data_4_35(:,4)];
d2_13 = diag((X13 - mean(X13))*inv(cov(X13))*(X13 - mean(X13))');
bivar_output(2,1) = sum(d2_13 <= chi50)/n;

% x2 and x3
X23 = [data_4_35(:,3) data_4_35(:,4)];
d2_23 = diag((X23 - mean(X23))*inv(cov(X23))*(X23 - mean(X23))');
bivar_output(3,1) = sum(d2_23 <= chi50)/n;

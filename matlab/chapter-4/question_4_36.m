% Exercise 4.36
% National Track Records for Women Data
data_4_36 = readmatrix(fullfile(data_folder, 'Table1.9.xlsx'));

% Compute simulation of critical values for Q-Q correlation test.
n = height(data_4_36);
[~, crit_4_36] = ppcc_simulation(n, 1000000, [0.01 0.05 0.10]);

% x1: 100m.
x1 = data_4_36(:,2);
x1_col = "100m";

% This returns 1, so 100m considered not normally distributed.
abs(sum((x1 >= (mean(x1) - 1*var(x1))) & (x1 <= (mean(x1) + 1*var(x1)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x1, r_Q_x1] = corr_q_q(x1);

% Okay at all three levels.
r_Q_x1 < crit_4_36

simple_qq_plot(qq_x1(:,1), qq_x1(:,3), append(x1_col, " (raw)"));
saveas(gcf, append(".\", "sol4.36.qq.1", ".png"))

% Perform power transformation.
power_x1 = box_cox_power_transform(x1, x1_col, ...
    append(".\", "sol4.36.power.1", ".png"));

% Use power transformation suggested.
x1_tr = x1.^power_x1;

[qq_x1_tr, r_Q_x1_tr] = corr_q_q(x1_tr);

% Transformation worked!
r_Q_x1_tr < crit_4_36

simple_qq_plot(qq_x1_tr(:,1), qq_x1_tr(:,3), append(x1_col, " (transformed)"));

sum(abs((x1 - mean(x1))/sqrt(var(x1))) > 3)
sum(abs((x1_tr - mean(x1_tr))/sqrt(var(x1_tr))) > 3)

clear x1 x1_tr x1_col r_Q_x1 r_Q_x1_tr power_x1 qq_x1 qq_x1_tr gcf

% x2: 200m.
x2 = data_4_36(:,3);
x2_col = "200m";

% This returns 0 for once, so 200m considered normally distributed.
abs(sum((x2 >= (mean(x2) - 1*var(x2))) & (x2 <= (mean(x2) + 1*var(x2)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x2, r_Q_x2] = corr_q_q(x2);

% Okay at the 0.01-level, but not 0.05 or 0.10.
r_Q_x2 < crit_4_36

simple_qq_plot(qq_x2(:,1), qq_x2(:,3), append(x2_col, " (raw)"));

% Perform power transformation.
power_x2 = box_cox_power_transform(x2, x2_col, ...
    -10,2, append(".\", "sol4.36.power.2", ".png"));

% Use power transformation suggested.
x2_tr = x2.^-6;

[qq_x2_tr, r_Q_x2_tr] = corr_q_q(x2_tr);

% Transformation worked!
r_Q_x2_tr < crit_4_36

simple_qq_plot(qq_x2_tr(:,1), qq_x2_tr(:,3), append(x2_col, " (transformed)"));

% Saved as sol4.36.qq.2.png
subplot(1,2,1);
simple_qq_plot(qq_x2(:,1), qq_x2(:,3), append(x2_col, " (raw)"))
subplot(1,2,2);
simple_qq_plot(qq_x2_tr(:,1), qq_x2_tr(:,3), append(x2_col, " (transformed)"));

sum(abs((x2 - mean(x2))/sqrt(var(x2))) > 3)
sum(abs((x2_tr - mean(x2_tr))/sqrt(var(x2_tr))) > 3)

clear x2 x2_tr x2_col r_Q_x2 r_Q_x2_tr power_x2 qq_x2 qq_x2_tr gcf

% x3: 400m.
x3 = data_4_36(:,4);
x3_col = "400m";

% This returns 1, so 400m considered not normally distributed.
abs(sum((x3 >= (mean(x3) - 1*var(x3))) & (x3 <= (mean(x3) + 1*var(x3)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x3, r_Q_x3] = corr_q_q(x3);

% Okay at the 0.01-level, but not 0.05 or 0.10.
r_Q_x3 < crit_4_36

simple_qq_plot(qq_x3(:,1), qq_x3(:,3), append(x3_col, " (raw)"));

% Perform power transformation.
power_x3 = box_cox_power_transform(x3, x3_col, ...
    -7, 2, append(".\", "sol4.36.power.3", ".png"));

% Use power transformation suggested.
x3_tr = x3.^-4.75;

[qq_x3_tr, r_Q_x3_tr] = corr_q_q(x3_tr);

% Transformation worked!
r_Q_x3_tr < crit_4_36

simple_qq_plot(qq_x3_tr(:,1), qq_x3_tr(:,3), append(x3_col, " (transformed)"));

% Saved as sol4.36.qq.3.png
subplot(1,2,1);
simple_qq_plot(qq_x3(:,1), qq_x3(:,3), append(x3_col, " (raw)"))
subplot(1,2,2);
simple_qq_plot(qq_x3_tr(:,1), qq_x3_tr(:,3), append(x3_col, " (transformed)"));

sum(abs((x3 - mean(x3))/sqrt(var(x3))) > 3)
sum(abs((x3_tr - mean(x3_tr))/sqrt(var(x3_tr))) > 3)

clear x3 x3_tr x3_col r_Q_x3 r_Q_x3_tr power_x3 qq_x3 qq_x3_tr gcf

% x4: 800m.
x4 = data_4_36(:,5);
x4_col = "800m";

% This returns 1, so 800m considered not normally distributed.
abs(sum((x4 >= (mean(x4) - 1*var(x4))) & (x4 <= (mean(x4) + 1*var(x4)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x4, r_Q_x4] = corr_q_q(x4);

% True for all, so not considered normal distributed at 0.01, 0.05 or 0.10.
r_Q_x4 < crit_4_36

simple_qq_plot(qq_x4(:,1), qq_x4(:,3), append(x4_col, " (raw)"));

% Perform power transformation.
power_x4 = box_cox_power_transform(x4, x4_col, ...
    -10, 2, append(".\", "sol4.36.power.4", ".png"));

% Use power transformation suggested.
x4_tr = x4.^-9.5;

[qq_x4_tr, r_Q_x4_tr] = corr_q_q(x4_tr);

% Transformation worked!
r_Q_x4_tr < crit_4_36

simple_qq_plot(qq_x4_tr(:,1), qq_x4_tr(:,3), append(x4_col, " (transformed)"));

% Saved as sol4.36.qq.4.png
subplot(1,2,1);
simple_qq_plot(qq_x4(:,1), qq_x4(:,3), append(x4_col, " (raw)"))
subplot(1,2,2);
simple_qq_plot(qq_x4_tr(:,1), qq_x4_tr(:,3), append(x4_col, " (transformed)"));

sum(abs((x4 - mean(x4))/sqrt(var(x4))) > 3)
sum(abs((x4_tr - mean(x4_tr))/sqrt(var(x4_tr))) > 3)

clear x4 x4_tr x4_col r_Q_x4 r_Q_x4_tr power_x4 qq_x4 qq_x4_tr gcf

% x5: 1500m.
x5 = data_4_36(:,6);
x5_col = "1500m";

% This returns 1, so 1500m considered not normally distributed.
abs(sum((x5 >= (mean(x5) - 1*var(x5))) & (x5 <= (mean(x5) + 1*var(x5)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x5, r_Q_x5] = corr_q_q(x5);

% True for all, so not considered normal distributed at 0.01, 0.05 or 0.10.
r_Q_x5 < crit_4_36

simple_qq_plot(qq_x5(:,1), qq_x5(:,3), append(x5_col, " (raw)"));

% Perform power transformation.
power_x5 = box_cox_power_transform(x5, x5_col, ...
    -10, 2, append(".\", "sol4.36.power.5", ".png"));

% Use power transformation suggested.
x5_tr = x5.^-8;

[qq_x5_tr, r_Q_x5_tr] = corr_q_q(x5_tr);

% Transformation worked!
r_Q_x5_tr < crit_4_36

simple_qq_plot(qq_x5_tr(:,1), qq_x5_tr(:,3), append(x5_col, " (transformed)"));

% Saved as sol4.36.qq.5.png
subplot(1,2,1);
simple_qq_plot(qq_x5(:,1), qq_x5(:,3), append(x5_col, " (raw)"))
subplot(1,2,2);
simple_qq_plot(qq_x5_tr(:,1), qq_x5_tr(:,3), append(x5_col, " (transformed)"));

sum(abs((x5 - mean(x5))/sqrt(var(x5))) > 3)
sum(abs((x5_tr - mean(x5_tr))/sqrt(var(x5_tr))) > 3)

clear x5 x5_tr x5_col r_Q_x5 r_Q_x5_tr power_x5 qq_x5 qq_x5_tr gcf

% x6: 3000m.
x6 = data_4_36(:,7);
x6_col = "3000m";

% This actually returns 10, so 3000m considered to be normally distributed.
abs(sum((x6 >= (mean(x6) - 1*var(x6))) & (x6 <= (mean(x6) + 1*var(x6)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x6, r_Q_x6] = corr_q_q(x6);

% True for all, so not considered normal distributed at 0.01, 0.05 or 0.10.
r_Q_x6 < crit_4_36

simple_qq_plot(qq_x6(:,1), qq_x6(:,3), append(x6_col, " (raw)"));

% Perform power transformation.
power_x6 = box_cox_power_transform(x6, x6_col, ...
    -10, 2, append(".\", "sol4.36.power.6", ".png"));

% Use power transformation suggested.
x6_tr = x6.^-7.2;

[qq_x6_tr, r_Q_x6_tr] = corr_q_q(x6_tr);

% Transformation worked!
r_Q_x6_tr < crit_4_36

simple_qq_plot(qq_x6_tr(:,1), qq_x6_tr(:,3), append(x6_col, " (transformed)"));

% Saved as sol4.36.qq.6.png
subplot(1,2,1);
simple_qq_plot(qq_x6(:,1), qq_x6(:,3), append(x6_col, " (raw)"))
subplot(1,2,2);
simple_qq_plot(qq_x6_tr(:,1), qq_x6_tr(:,3), append(x6_col, " (transformed)"));

sum(abs((x6 - mean(x6))/sqrt(var(x6))) > 3)
sum(abs((x6_tr - mean(x6_tr))/sqrt(var(x6_tr))) > 3)

clear x6 x6_tr x6_col r_Q_x6 r_Q_x6_tr power_x6 qq_x6 qq_x6_tr gcf

% x7: Marathon.
x7 = data_4_36(:,8);
x7_col = "Marathon";

% This returns 1, so marathon time considered not normally distributed.
abs(sum((x7 >= (mean(x7) - 1*var(x7))) & (x7 <= (mean(x7) + 1*var(x7)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x7, r_Q_x7] = corr_q_q(x7);

% Okay at the 0.01-level, but not 0.05 or 0.10.
r_Q_x7 < crit_4_36

simple_qq_plot(qq_x7(:,1), qq_x7(:,3), append(x7_col, " (raw)"));

% Perform power transformation.
power_x7 = box_cox_power_transform(x7, x7_col, ...
    -7, 2, append(".\", "sol4.36.power.7", ".png"));

% Use power transformation suggested.
x7_tr = x7.^-7;

[qq_x7_tr, r_Q_x7_tr] = corr_q_q(x7_tr);

% Transformation worked!
r_Q_x7_tr < crit_4_36

simple_qq_plot(qq_x7_tr(:,1), qq_x7_tr(:,3), append(x7_col, " (transformed)"));

% Saved as sol4.36.qq.7.png
subplot(1,2,1);
simple_qq_plot(qq_x7(:,1), qq_x7(:,3), append(x7_col, " (raw)"))
subplot(1,2,2);
simple_qq_plot(qq_x7_tr(:,1), qq_x7_tr(:,3), append(x7_col, " (transformed)"));

sum(abs((x7 - mean(x7))/sqrt(var(x7))) > 3)
sum(abs((x7_tr - mean(x7_tr))/sqrt(var(x7_tr))) > 3)

clear x7 x7_tr x7_col r_Q_x7 r_Q_x7_tr power_x7 qq_x7 qq_x7_tr gcf

% Check multivariate normal.
X = data_4_36(:,2:8);
chi50 = icdf('Chi2',0.5, width(X));

output = chi2_q_q(X);

col_tr = [1 -6 -4.75 -9.5 -8 -7.2 -7];
X_tr = zeros(height(X), width(X));
for i=1:width(X)
    X_tr(:,i) = X(:,i).^col_tr(i);
end

output_tr = chi2_q_q(X_tr);

subplot(1,2,1)
scatter(output(:,3), output(:,1))
title('\chi^{2} Plot National Track Records for Women (Raw)')
xlabel('$q_{c,7}\left(\frac{j - 0.5}{54}\right)$', 'Interpreter', 'latex')
ylabel('d_{(j)}^{2}')
subplot(1,2,2)
scatter(output_tr (:,3), output_tr (:,1))
title('\chi^{2} Plot National Track Records for Women (Transformed)')
xlabel('$q_{c,7}\left(\frac{j - 0.5}{54}\right)$', 'Interpreter', 'latex')
ylabel('d_{(j)}^{2}')
saveas(p, append('.\', 'sol4.36chi2', '.png'), 'png')

sum(output(:,1) < icdf('Chi2', 0.5, width(X)))/n
sum(output_tr(:,1) < icdf('Chi2', 0.5, width(X)))/n

% Find the 3 largest statistical distances.
asd = diag((X_tr - mean(X_tr))/cov(X_tr)*(X_tr - mean(X_tr))');
[a, b] = maxk(asd, 3)

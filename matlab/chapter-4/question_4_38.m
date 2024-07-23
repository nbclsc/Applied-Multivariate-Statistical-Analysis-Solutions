% Exercise 4.38
% Bulls Data
data_4_38 = readmatrix(fullfile(data_folder, 'Table1.10.xlsx'));
X = data_4_38(:,[3:5 7:9]);

% (3) YrHgt = Yearling height at shoulder (inches)
% (4) FtFrBody = Fat free body (pounds)
% (5) PrctFFB = Percent fat-free body
% (7) BkFat = Back fat (inches)
% (8) SaleHt = Sale height at shoulder (inches)
% (9) Sale Wt = Sale weight

% Compute simulation of critical values for Q-Q correlation test.
n = height(X);
[~, crit_4_38] = ppcc_simulation(n, 1000000, [0.01 0.05 0.10]);

% x1: Yearling Height at Shoulder (inches).
x1 = X(:,1);
x1_col = "Yearling Height at Shoulder (inches)";

% This returns 1, so YrHgt considered not normally distributed.
abs(sum((x1 >= (mean(x1) - 1*var(x1))) & (x1 <= (mean(x1) + 1*var(x1)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x1, r_Q_x1] = corr_q_q(x1);

% Okay at the 0.01-level, but not 0.05 or 0.10.
r_Q_x1 < crit_4_38

simple_qq_plot(qq_x1(:,1), qq_x1(:,3), append(x1_col, " (raw)"));
saveas(gcf, append(".\", "sol4.38.qq.1", ".png"))

% Perform power transformation.
power_x1 = box_cox_power_transform(x1, x1_col, ...
    -5, 5, append(".\", "sol4.38.power.1", ".png"));

% Use power transformation suggested.
x1_tr = x1.^-2.5;

[qq_x1_tr, r_Q_x1_tr] = corr_q_q(x1_tr);

% Transformation worked!
r_Q_x1_tr < crit_4_38

simple_qq_plot(qq_x1_tr(:,1), qq_x1_tr(:,3), append(x1_col, " (transformed)"));

sum(abs((x1 - mean(x1))/sqrt(var(x1))) > 3)
sum(abs((x1_tr - mean(x1_tr))/sqrt(var(x1_tr))) > 3)

clear x1 x1_tr x1_col r_Q_x1 r_Q_x1_tr power_x1 qq_x1 qq_x1_tr gcf

% x2: Fat Free Body (pounds).
x2 = X(:,2);
x2_col = "Fat Free Body (pounds)";

% This returns 1, so FtFrBody considered not normally distributed.
abs(sum((x2 >= (mean(x2) - 1*var(x2))) & (x2 <= (mean(x2) + 1*var(x2)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x2, r_Q_x2] = corr_q_q(x2);

% True for all, so not considered normal at any of the 3 levels.
r_Q_x2 < crit_4_38

simple_qq_plot(qq_x2(:,1), qq_x2(:,3), append(x2_col, " (raw)"));

% Perform power transformation.
power_x2 = box_cox_power_transform(x2, x2_col, ...
    -5, 5, append(".\", "sol4.38.power.2", ".png"));

% Use power transformation suggested.
x2_tr = x2.^-2.2;

[qq_x2_tr, r_Q_x2_tr] = corr_q_q(x2_tr);

% Transformation worked!
r_Q_x2_tr < crit_4_38

simple_qq_plot(qq_x2_tr(:,1), qq_x2_tr(:,3), append(x2_col, " (transformed)"));

% Saved as sol4.38.qq.2.png
subplot(1,2,1);
simple_qq_plot(qq_x2(:,1), qq_x2(:,3), append(x2_col, " (raw)"))
subplot(1,2,2);
simple_qq_plot(qq_x2_tr(:,1), qq_x2_tr(:,3), append(x2_col, " (transformed)"));

sum(abs((x2 - mean(x2))/sqrt(var(x2))) > 3)
sum(abs((x2_tr - mean(x2_tr))/sqrt(var(x2_tr))) > 3)

clear x2 x2_tr x2_col r_Q_x2 r_Q_x2_tr power_x2 qq_x2 qq_x2_tr gcf

% x3: Percent Fat-Free Body.
x3 = X(:,3);
x3_col = "Percent Fat-Free Body";

% This returns 1, so PrctFFB considered not normally distributed.
abs(sum((x3 >= (mean(x3) - 1*var(x3))) & (x3 <= (mean(x3) + 1*var(x3)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x3, r_Q_x3] = corr_q_q(x3);

% Okay at the 0.01-level and 0.05, but not 0.10.
r_Q_x3 < crit_4_38

simple_qq_plot(qq_x3(:,1), qq_x3(:,3), append(x3_col, " (raw)"));

% Perform power transformation.
power_x3 = box_cox_power_transform(x3, x3_col, ...
    -5, 5, append(".\", "sol4.38.power.3", ".png"));

% Use power transformation suggested.
x3_tr = x3.^-3.5;

[qq_x3_tr, r_Q_x3_tr] = corr_q_q(x3_tr);

% Transformation worked!
r_Q_x3_tr < crit_4_38

simple_qq_plot(qq_x3_tr(:,1), qq_x3_tr(:,3), append(x3_col, " (transformed)"));

% Saved as sol4.38.qq.3.png
subplot(1,2,1);
simple_qq_plot(qq_x3(:,1), qq_x3(:,3), append(x3_col, " (raw)"))
subplot(1,2,2);
simple_qq_plot(qq_x3_tr(:,1), qq_x3_tr(:,3), append(x3_col, " (transformed)"));

sum(abs((x3 - mean(x3))/sqrt(var(x3))) > 3)
sum(abs((x3_tr - mean(x3_tr))/sqrt(var(x3_tr))) > 3)

clear x3 x3_tr x3_col r_Q_x3 r_Q_x3_tr power_x3 qq_x3 qq_x3_tr gcf

% x4: Back Fat (inches).
x4 = X(:,4);
x4_col = "Back Fat (inches)";

% This returns 1, so BkFat considered not normally distributed.
abs(sum((x4 >= (mean(x4) - 1*var(x4))) & (x4 <= (mean(x4) + 1*var(x4)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x4, r_Q_x4] = corr_q_q(x4);

% True for all values, so not considered normally distributed at any level.
r_Q_x4 < crit_4_38

simple_qq_plot(qq_x4(:,1), qq_x4(:,3), append(x4_col, " (raw)"));

% Perform power transformation.
power_x4 = box_cox_power_transform(x4, x4_col, ...
    -5, 5, append(".\", "sol4.38.power.4", ".png"));

% Use power transformation suggested.
x4_tr = x4.^-0.30;

[qq_x4_tr, r_Q_x4_tr] = corr_q_q(x4_tr);

% Transformation worked!
r_Q_x4_tr < crit_4_38

simple_qq_plot(qq_x4_tr(:,1), qq_x4_tr(:,3), append(x4_col, " (transformed)"));

% Saved as sol4.38.qq.4.png
subplot(1,2,1);
simple_qq_plot(qq_x4(:,1), qq_x4(:,3), append(x4_col, " (raw)"))
subplot(1,2,2);
simple_qq_plot(qq_x4_tr(:,1), qq_x4_tr(:,3), append(x4_col, " (transformed)"));

sum(abs((x4 - mean(x4))/sqrt(var(x4))) > 3)
sum(abs((x4_tr - mean(x4_tr))/sqrt(var(x4_tr))) > 3)

clear x4 x4_tr x4_col r_Q_x4 r_Q_x4_tr power_x4 qq_x4 qq_x4_tr gcf

% x5: Sale Height at Shoulder (inches).
x5 = X(:,5);
x5_col = "Sale Height at Shoulder (inches)";

% This returns 1, so SaleHt considered not normally distributed.
abs(sum((x5 >= (mean(x5) - 1*var(x5))) & (x5 <= (mean(x5) + 1*var(x5)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x5, r_Q_x5] = corr_q_q(x5);

% 0 for all, so considered normally distributed at 0.01, 0.05, and 0.10 levels.
r_Q_x5 < crit_4_38

simple_qq_plot(qq_x5(:,1), qq_x5(:,3), append(x5_col, " (raw)"));
saveas(gcf, append(".\", "sol4.38.qq.5", ".png"))

% Perform power transformation.
power_x5 = box_cox_power_transform(x5, x5_col, ...
    -5, 5, append(".\", "sol4.38.power.5", ".png"));

% Use power transformation suggested.
x5_tr = log(x5);

[qq_x5_tr, r_Q_x5_tr] = corr_q_q(x5_tr);

% Transformation worked!
r_Q_x5_tr < crit_4_38

simple_qq_plot(qq_x5_tr(:,1), qq_x5_tr(:,3), append(x5_col, " (transformed)"));

sum(abs((x5 - mean(x5))/sqrt(var(x5))) > 3)
sum(abs((x5_tr - mean(x5_tr))/sqrt(var(x5_tr))) > 3)

clear x5 x5_tr x5_col r_Q_x5 r_Q_x5_tr power_x5 qq_x5 qq_x5_tr gcf

% x6: Sale Weight (pounds).
x6 = X(:,6);
x6_col = "Sale Weight (pounds)";

% This returns 1, so SaleWt considered not normally distributed.
abs(sum((x6 >= (mean(x6) - 1*var(x6))) & (x6 <= (mean(x6) + 1*var(x6)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x6, r_Q_x6] = corr_q_q(x6);

% 0 for all, so considered normally distributed at 0.01, 0.05, and 0.10 levels.
r_Q_x6 < crit_4_38

simple_qq_plot(qq_x6(:,1), qq_x6(:,3), append(x6_col, " (raw)"));
saveas(gcf, append(".\", "sol4.38.qq.6", ".png"))

% Perform power transformation.
power_x6 = box_cox_power_transform(x6, x6_col, ...
    -5, 5, append(".\", "sol4.38.power.6", ".png"));

% Use power transformation suggested.
x6_tr = x6.^-0.50;

[qq_x6_tr, r_Q_x6_tr] = corr_q_q(x6_tr);

% Transformation worked!
r_Q_x6_tr < crit_4_38

simple_qq_plot(qq_x6_tr(:,1), qq_x6_tr(:,3), append(x6_col, " (transformed)"));

sum(abs((x6 - mean(x6))/sqrt(var(x6))) > 3)
sum(abs((x6_tr - mean(x6_tr))/sqrt(var(x6_tr))) > 3)

clear x6 x6_tr x6_col r_Q_x6 r_Q_x6_tr power_x6 qq_x6 qq_x6_tr gcf

% Check multivariate normal.
chi50 = icdf('Chi2',0.5, width(X));

output = chi2_q_q(X);

col_tr = [1 -2.2 -3.5 -0.3 0 1];
X_tr = zeros(height(X), width(X));
for i=1:width(X)
    if col_tr(i) ~= 0
        X_tr(:,i) = X(:,i).^col_tr(i);
    else
        X_tr(:,i) = log(X(:,i));
    end
end

output_tr = chi2_q_q(X_tr);

subplot(1,2,1)
scatter(output(:,3), output(:,1))
title('\chi^{2} Plot Bull Data (Raw)')
xlabel('$q_{c,7}\left(\frac{j - 0.5}{54}\right)$', 'Interpreter', 'latex')
ylabel('d_{(j)}^{2}')
subplot(1,2,2)
scatter(output_tr (:,3), output_tr (:,1))
title('\chi^{2} Plot Bull Data (Transformed)')
xlabel('$q_{c,7}\left(\frac{j - 0.5}{54}\right)$', 'Interpreter', 'latex')
ylabel('d_{(j)}^{2}')
saveas(gcf, append('.\', 'sol4.38chi2', '.png'), 'png')

sum(output(:,1) < icdf('Chi2', 0.5, width(X)))/n
sum(output_tr(:,1) < icdf('Chi2', 0.5, width(X)))/n

% Find the 3 largest statistical distances.
asd = diag((X_tr - mean(X_tr))/cov(X_tr)*(X_tr - mean(X_tr))');
[a, b] = maxk(asd, 3)


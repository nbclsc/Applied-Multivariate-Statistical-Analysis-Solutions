% Exercise 4.34
% Bone Material Content Data
data_4_34 = readmatrix(fullfile(data_folder, 'Table1.8.xlsx'));

% Compute simulation of critical values for Q-Q correlation test.
n = height(data_4_34);
[~, crit_4_34] = ppcc_simulation(n, 1000000, [0.01 0.05 0.10]);

% x1: Dominant Radius.
x1 = data_4_34(:,2);
x1_col = "Dominant Radius";

% This returns 1, so dominant radius considered not normally distributed.
abs(sum((x1 >= (mean(x1) - 1*var(x1))) & (x1 <= (mean(x1) + 1*var(x1)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x1, r_Q_x1] = corr_q_q(x1);

% Okay at the 0.01-level, but not 0.05 or 0.10.
r_Q_x1 < crit_4_34

simple_qq_plot(qq_x1(:,1), qq_x1(:,3), append(x1_col, " (raw)"));

% Perform power transformation.
power_x1 = box_cox_power_transform(x1, x1_col, ...
    -5, 5, append(".\", "sol4.34.power.1", ".png"));

% Use power transformation suggested.
x1_tr = x1.^power_x1;

[qq_x1_tr, r_Q_x1_tr] = corr_q_q(x1_tr);

% Transformation worked!
r_Q_x1_tr < crit_4_34

simple_qq_plot(qq_x1_tr(:,1), qq_x1_tr(:,3), append(x1_col, " (transformed)"));

% Saved as sol4.34.qq.1.png
subplot(1,2,1);
simple_qq_plot(qq_x1(:,1), qq_x1(:,3), append(x1_col, " (raw)"))
subplot(1,2,2);
simple_qq_plot(qq_x1_tr(:,1), qq_x1_tr(:,3), append(x1_col, " (transformed)"));

sum(abs((x1 - mean(x1))/sqrt(var(x1))) > 3)
sum(abs((x1_tr - mean(x1_tr))/sqrt(var(x1_tr))) > 3)

clear x1 x1_tr x1_col r_Q_x1 r_Q_x1_tr power_x1 qq_x1 qq_x1_tr gcf

% x2: Radius.
x2 = data_4_34(:,3);
x2_col = "Radius";

% This returns 1, so radius considered not normally distributed.
abs(sum((x2 >= (mean(x2) - 1*var(x2))) & (x2 <= (mean(x2) + 1*var(x2)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x2, r_Q_x2] = corr_q_q(x2);

% Okay at the 0.01-level, but not 0.05 or 0.10.
r_Q_x2 < crit_4_34

simple_qq_plot(qq_x2(:,1), qq_x2(:,3), append(x2_col, " (raw)"));
saveas(gcf, append(".\", "sol4.34.qq.2", ".png"))

% Perform power transformation.
power_x2 = box_cox_power_transform(x2, x2_col, ...
    -5, 5, append(".\", "sol4.34.power.2", ".png"));

% Use power transformation suggested.
x2_tr = x2.^power_x2;

[qq_x2_tr, r_Q_x2_tr] = corr_q_q(x2_tr);

% Transformation was an omprovment on data already considered normal.
r_Q_x2_tr < crit_4_34

simple_qq_plot(qq_x2_tr(:,1), qq_x2_tr(:,3), append(x2_col, " (transformed)"));

sum(abs((x2 - mean(x2))/sqrt(var(x2))) > 3)
sum(abs((x2_tr - mean(x2_tr))/sqrt(var(x2_tr))) > 3)

clear x2 x2_tr x2_col r_Q_x2 r_Q_x2_tr power_x2 qq_x2 qq_x2_tr gcf

% x3: Dominant Humerus.
x3 = data_4_34(:,4);
x3_col = "Dominant Humerus";

% This returns 1, so dominant hunerus considered not normally distributed.
abs(sum((x3 >= (mean(x3) - 1*var(x3))) & (x3 <= (mean(x3) + 1*var(x3)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x3, r_Q_x3] = corr_q_q(x3);

% Okay at the 0.01-level, but not 0.05 or 0.10.
r_Q_x3 < crit_4_34

simple_qq_plot(qq_x3(:,1), qq_x3(:,3), append(x3_col, " (raw)"));
saveas(gcf, append(".\", "sol4.34.qq.3", ".png"))

% Perform power transformation.
power_x3 = box_cox_power_transform(x3, x3_col, ...
    -5, 5, append(".\", "sol4.34.power.3", ".png"));

% Use power transformation suggested.
x3_tr = x3.^power_x3;

[qq_x3_tr, r_Q_x3_tr] = corr_q_q(x3_tr);

% Transformation worked!
r_Q_x3_tr < crit_4_34

simple_qq_plot(qq_x3_tr(:,1), qq_x3_tr(:,3), append(x3_col, " (transformed)"));

sum(abs((x3 - mean(x3))/sqrt(var(x3))) > 3)
sum(abs((x3_tr - mean(x3_tr))/sqrt(var(x3_tr))) > 3)

clear x3 x3_tr x3_col r_Q_x3 r_Q_x3_tr power_x3 qq_x3 qq_x3_tr gcf

% x4: Humerus.
x4 = data_4_34(:,5);
x4_col = "Humerus";

% This returns 1, so humerus considered not normally distributed.
abs(sum((x4 >= (mean(x4) - 1*var(x4))) & (x4 <= (mean(x4) + 1*var(x4)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x4, r_Q_x4] = corr_q_q(x4);

% Okay at the 0.01-level, but not 0.05 or 0.10.
r_Q_x4 < crit_4_34

simple_qq_plot(qq_x4(:,1), qq_x4(:,3), append(x4_col, " (raw)"));
saveas(gcf, append(".\", "sol4.34.qq.4", ".png"))

% Perform power transformation.
power_x4 = box_cox_power_transform(x4, x4_col, ...
    -5, 5, append(".\", "sol4.34.power.4", ".png"));

% Use power transformation suggested.
x4_tr = x4.^power_x4;

[qq_x4_tr, r_Q_x4_tr] = corr_q_q(x4_tr);

% Transformation worked!
r_Q_x4_tr < crit_4_34

simple_qq_plot(qq_x4_tr(:,1), qq_x4_tr(:,3), append(x4_col, " (transformed)"));

sum(abs((x4 - mean(x4))/sqrt(var(x4))) > 3)
sum(abs((x4_tr - mean(x4_tr))/sqrt(var(x4_tr))) > 3)

clear x4 x4_tr x4_col r_Q_x4 r_Q_x4_tr power_x4 qq_x4 qq_x4_tr gcf

% x5: Dominant Ulna.
x5 = data_4_34(:,6);
x5_col = "Dominant Ulna";

% This returns 1, so dominant ulna considered not normally distributed.
abs(sum((x5 >= (mean(x5) - 1*var(x5))) & (x5 <= (mean(x5) + 1*var(x5)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x5, r_Q_x5] = corr_q_q(x5);

% Okay at the 0.01-level, but not 0.05 or 0.10.
r_Q_x5 < crit_4_34

simple_qq_plot(qq_x5(:,1), qq_x5(:,3), append(x5_col, " (raw)"));
saveas(gcf, append(".\", "sol4.34.qq.5", ".png"))

% Perform power transformation.
power_x5 = box_cox_power_transform(x5, x5_col, ...
    -5, 5, append(".\", "sol4.34.power.5", ".png"));

% Use power transformation suggested.
x5_tr = x5.^power_x5;

[qq_x5_tr, r_Q_x5_tr] = corr_q_q(x5_tr);

% Transformation worked!
r_Q_x5_tr < crit_4_34

simple_qq_plot(qq_x5_tr(:,1), qq_x5_tr(:,3), append(x5_col, " (transformed)"));

sum(abs((x5 - mean(x5))/sqrt(var(x5))) > 3)
sum(abs((x5_tr - mean(x5_tr))/sqrt(var(x5_tr))) > 3)

clear x5 x5_tr x5_col r_Q_x5 r_Q_x5_tr power_x5 qq_x5 qq_x5_tr gcf

% x6: Ulna.
x6 = data_4_34(:,7);
x6_col = "Ulna";

% This returns 1, so ulna considered not normally distributed.
abs(sum((x6 >= (mean(x6) - 1*var(x6))) & (x6 <= (mean(x6) + 1*var(x6)))) / n - 0.683) > 1.396/sqrt(n)

% Look at the QQ plot of the raw data.
[qq_x6, r_Q_x6] = corr_q_q(x6);

% Okay at the 0.01-level, but not 0.05 or 0.10.
r_Q_x6 < crit_4_34

simple_qq_plot(qq_x6(:,1), qq_x6(:,3), append(x6_col, " (raw)"));
saveas(gcf, append(".\", "sol4.34.qq.6", ".png"))

% Perform power transformation.
power_x6 = box_cox_power_transform(x6, x6_col, ...
    -5, 5, append(".\", "sol4.34.power.6", ".png"));

% Use power transformation suggested.
x6_tr = x6.^power_x6;

[qq_x6_tr, r_Q_x6_tr] = corr_q_q(x6_tr);

% Transformation worked!
r_Q_x6_tr < crit_4_34

simple_qq_plot(qq_x6_tr(:,1), qq_x6_tr(:,3), append(x6_col, " (transformed)"));

sum(abs((x6 - mean(x6))/sqrt(var(x6))) > 3)
sum(abs((x6_tr - mean(x6_tr))/sqrt(var(x6_tr))) > 3)

clear x6 x6_tr x6_col r_Q_x6 r_Q_x6_tr power_x6 qq_x6 qq_x6_tr gcf

% Checking bivariate normality.
chi50 = icdf('Chi2',0.5, 2);

% Store proportions. Column 1 is raw data. Column 2 is transformed data.
bivar_output = zeros(nchoosek(width(data_4_34)-1,2), 2);

% List of transformations for each variable. Only x1 had a transformation.
col_tr = [2.2144 1 1 1 1 1];

% Compute statistical distance for pairs of variables and compute
% proportion of data within Chi2 CDF value where 50% of the data exist 
% on 2 degrees of freedom.
row_number = 0;
for i=1:(width(data_4_34)-1)
    incr = i + 1;
    for j=incr:(width(data_4_34)-1)
        row_number = row_number + 1;

        % Raw data.
        X = [data_4_34(:,i+1) data_4_34(:,j+1)];
        d2 = diag((X - mean(X))/cov(X)*(X - mean(X))');
        bivar_output(row_number,1) = sum(d2 <= chi50)/n;

        % Transformed data.
        X_tr = [data_4_34(:,i+1).^col_tr(i) data_4_34(:,j+1).^col_tr(j)];
        d2_tr = diag((X_tr - mean(X_tr))/cov(X_tr)*(X_tr - mean(X_tr))');
        bivar_output(row_number,2) = sum(d2_tr <= chi50)/n;

        clear X d2 X_tr d2_tr
    end
end

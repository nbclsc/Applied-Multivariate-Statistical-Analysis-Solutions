% Exercise 4.41
% Snow Data

% In northern climates, roads must be cleared of snow quickly following a
% storm. One measure of storm severity is
% x1 = its duration in hours
% while the effectiveness of snow removal can be quantified by
% x2 = the number of hours crews, men, and machine, spend to clear snow.
% Here are the results for 25 incidents in Wisconsin.

data_4_41 = readtable(fullfile(data_folder, 'Table3.2.xlsx'));
X = table2array(data_4_41);

% Exercise 4.41 (a)

scatter(X(:,1), X(:,2))
xlabel("Snow Duration (hours)")
ylabel("Work (hours)")
title("Snow Data")
saveas(gcf, append('.\', 'sol4.41.a', '.png'), 'png')

% Exercise 4.41 (b)

% Compute simulation of critical values for Q-Q correlation test.
n = height(X);
[~, crit_4_41] = ppcc_simulation(n, 1000000, [0.01 0.05 0.10]);

x1 = X(:,1);
x1_col = "Snow Duration (hours)";

% This returns 1, so duration is considered normally distributed by this metric.
abs(sum((x1 >= (mean(x1) - 1*var(x1))) & (x1 <= (mean(x1) + 1*var(x1)))) / n - 0.683) > 1.396/sqrt(n)

% Perform power transformation.
power_x1 = box_cox_power_transform(x1, x1_col, ...
    -5, 5, append(".\", "sol4.41.power.1", ".png"));

% Use power transformation suggested.
x1_tr = log(x1);

[qq_x1_tr, r_Q_x1_tr] = corr_q_q(x1_tr);

% Transformation worked!
r_Q_x1_tr < crit_4_41

simple_qq_plot(qq_x1_tr(:,1), qq_x1_tr(:,3), append(x1_col, " (transformed)"));
saveas(gcf, append(".\", "sol4.41.qq.tr.1", ".png"))

sum(abs((x1 - mean(x1))/sqrt(var(x1))) > 3)
sum(abs((x1_tr - mean(x1_tr))/sqrt(var(x1_tr))) > 3)

clear x1 x1_tr x1_col power_x1 qq_x1_tr r_Q_x1_tr 

% Exercise 4.41 (c)

x2 = X(:,2);
x2_col = "Work (hours)";

% This returns 1, so work is not normally distributed by this metric.
abs(sum((x2 >= (mean(x2) - 1*var(x2))) & (x2 <= (mean(x2) + 1*var(x2)))) / n - 0.683) > 1.396/sqrt(n)

% Perform power transformation.
power_x2 = box_cox_power_transform(x2, x2_col, ...
    -5, 5, append(".\", "sol4.41.power.2", ".png"));

% Use power transformation suggested.
x2_tr = x2.^-0.70;

[qq_x2_tr, r_Q_x2_tr] = corr_q_q(x2_tr);

% Transformation worked!
r_Q_x2_tr < crit_4_41

simple_qq_plot(qq_x2_tr(:,1), qq_x2_tr(:,3), append(x2_col, " (transformed)"));
saveas(gcf, append(".\", "sol4.41.qq.tr.2", ".png"))

sum(abs((x2 - mean(x2))/sqrt(var(x2))) > 3)
sum(abs((x2_tr - mean(x2_tr))/sqrt(var(x2_tr))) > 3)

clear x2 x2_tr x2_col power_x2 qq_x2_tr r_Q_x2_tr

% Exercise 4.40 (d)

% Set up 2D lambda values.
lambda1 = linspace(-1, 1, 150);
lambda2 = lambda1;
% Meshgrid 1st argument is x-axis (cols) and 2nd argument is y-axis (rows).
[l1, l2] = meshgrid(lambda1, lambda2);
% [l1, l2] = ngrid(lambda1, lambda2);

% Transformed observations.
x1_lambda = (X(:,1).^(lambda1) - 1) ./ lambda1; % 25 x 150
x2_lambda = (X(:,2).^(lambda2) - 1) ./ lambda2; % 25 x 150
if sum(lambda1==0)
    x1_lambda(:, lambda1==0) = log(X(:,1));
end
if sum(lambda2==0)
    x2_lambda(:, lambda2==0) = log(X(:,2));
end

% Compute values for function to maximize equation (4-40).
% To be consistent with meshgrid, rows are X2 and columns are X1.
l_lambda = zeros(length(lambda2), length(lambda1));
n = height(X);
part1 = sum(log(X(:,1)));  % Variable X1 = Size (acres)
part2 = sum(log(X(:,2)));  % Variable X2 = Visitors (millions)
% Index i is rows (X2) and j is columns (X1).
for i=1:length(lambda2)
    for j=1:length(lambda1)
        l_lambda(i,j) = -(n/2)*log(det(cov(x2_lambda(:,i), x1_lambda(:,j)))) + (lambda1(j) - 1) * part1 + (lambda2(i) - 1) * part2;
    end
end

% Find the maximum by stacking columns (linear).
[a, max_idx] = max(l_lambda(:));
% We used l1 on x-axis for meshgrid, which corresponds to the columns and
% l2 on the y-axis corresponding to the rows. The same for matrix l_lambda,
% the rows correspond to l2 and columns to l1, that's why l2_max is first,
% since it's for rows and l1_max is second for columns in the output of
% applying ind2sub below.
[l2_max, l1_max] = ind2sub(size(l_lambda), max_idx);

% Plot the surface.
contour(l1, l2, l_lambda, 'ShowText', 'on')
title('Contour Plot of $\ell\left( \lambda_{1}, \lambda_{2} \right)$ for Snow Data', 'Interpreter', 'latex')
xlabel('\lambda_{1}')
ylabel('\lambda_{2}')
hold on
plot(lambda1(l1_max), lambda2(l2_max), 'o')
text(lambda1(l1_max), lambda2(l2_max), ...
    sprintf('(%4.4f, %4.4f)', lambda1(l1_max), lambda2(l2_max)), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
hold off
saveas(gcf, append('.\', 'sol4.41.d.contour', '.png'), 'png')

surf(l1, l2, l_lambda)
colormap default
hold on
plot3(lambda1(l1_max), lambda2(l2_max), a, ...
    'bo', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'b')
text(lambda1(l1_max), lambda2(l2_max), a+2, ...
    sprintf('(%5.3f, %5.3f)', lambda1(l1_max), lambda2(l2_max)), ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right')
hold off


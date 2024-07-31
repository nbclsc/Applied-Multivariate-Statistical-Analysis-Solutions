% Exercise 4.30
% Used Car Data

% Exercise 4.30 (a)
data_4_30 = readmatrix(fullfile(data_folder, 'Exercise1.2.xlsx'));

lambda = linspace(-1, 2, 100);

x1 = data_4_30(:, 1);
x1_lambda = (x1.^(lambda) - 1) ./ lambda;
if sum(lambda==0)
    x1_lambda(:, lambda==0) = log(x1);
end
x1_lambda_bar = (1/length(x1))*sum(x1_lambda);
l1_lambda = -(length(x1)/2)*log((1/length(x1))*sum((x1_lambda - x1_lambda_bar).^2)) + (lambda - 1)*sum(log(x1));
[max_lambda1, argmax_lambda1] = max(l1_lambda);
[min_lambda1, ~] = min(l1_lambda);
p = plot(lambda, l1_lambda);
line([lambda(argmax_lambda1) lambda(argmax_lambda1)], [min_lambda1 max_lambda1], 'Color', 'black', 'LineStyle', '--')
hold on
plot(lambda(argmax_lambda1), max_lambda1, 'o-','MarkerFaceColor','black','MarkerEdgeColor','black')
text(lambda(argmax_lambda1)+0.1, max_lambda1-2.5, sprintf('$\\hat{\\lambda}_{1} = %.4f$', lambda(argmax_lambda1)), 'Interpreter', 'latex')
title('Box Cox Plot Used Car Vehicle Age')
xlabel('$\lambda$', 'Interpreter', 'latex')
ylabel('$\ell(\lambda)$', 'Interpreter', 'latex');
hold off
saveas(p, append('.\', 'sol4.30a', '.png'), 'png')
clear x1 x1_lambda x1_lambda_bar l1_lambda max_lambda1 argmax_lambda1 min_lambda1 argmin_lambda1

% Exercise 4.30 (b)
lambda = linspace(-1, 2, 100);
x2 = data_4_30(:, 2);
x2_lambda = (x2.^(lambda) - 1) ./ lambda;
if sum(lambda==0)
    x2_lambda(:, lambda==0) = log(x2);
end
x2_lambda_bar = (1/length(x2))*sum(x2_lambda);
l2_lambda = -(length(x2)/2)*log((1/length(x2))*sum((x2_lambda - x2_lambda_bar).^2)) + (lambda -1)*sum(log(x2));
[max_lambda2, argmax_lambda2] = max(l2_lambda);
[min_lambda2, ~] = min(l2_lambda);
p = plot(lambda, l2_lambda);
line([lambda(argmax_lambda2) lambda(argmax_lambda2)], [min_lambda2 max_lambda2], 'Color', 'black', 'LineStyle', '--')
hold on
plot(lambda(argmax_lambda2), max_lambda2, 'o-','MarkerFaceColor','black','MarkerEdgeColor','black')
text(lambda(argmax_lambda2)+0.1, max_lambda2-2.5, sprintf('$\\hat{\\lambda}_{2} = %.4f$', lambda(argmax_lambda2)), 'Interpreter', 'latex')
title('Box Cox Plot Used Car Selling Price')
xlabel('$\lambda$', 'Interpreter', 'latex')
ylabel('$\ell(\lambda)$', 'Interpreter', 'latex');
hold off
saveas(p, append('.\', 'sol4.30b', '.png'), 'png')
clear x2 x2_lambda x2_lambda_bar l2_lambda max_lambda2 argmax_lambda2 min_lambda2 argmin_lambda2

% Exercise 4.30 (c)
X = data_4_30;

% Set up 2D lambda values. Using the same lambda for each variable.
N = 100;
% lambda1 = linspace(-2, 3, N);
% lambda2 = lambda1;

lambda1 = linspace(0, 2, N);
lambda2 = lambda1;
% Meshgrid 1st argument is x-axis (cols) and 2nd argument is y-axis (rows).
[l1, l2] = meshgrid(lambda1, lambda2);

% Transformed observations.
x1_lambda = (X(:,1).^(lambda1) - 1) ./ lambda1; % 10 x N
x2_lambda = (X(:,2).^(lambda2) - 1) ./ lambda2; % 10 x N
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
title('Contour Plot of $\ell\left( \lambda_{1}, \lambda_{2} \right)$ for Used Car Data', 'Interpreter', 'latex')
xlabel('\lambda_{1}')
ylabel('\lambda_{2}')
hold on
plot(lambda1(l1_max), lambda2(l2_max), 'o')
text(lambda1(l1_max), lambda2(l2_max), ...
    sprintf('(%4.4f, %4.4f)', lambda1(l1_max), lambda2(l2_max)), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
hold off
saveas(gcf, append('.\', 'sol4.30c', '.png'), 'png')

surf(l1, l2, l_lambda)
colormap default
hold on
plot3(lambda1(l1_max), lambda2(l2_max), a, ...
    'bo', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'b')
text(lambda1(l1_max), lambda2(l2_max), a+5, ...
    sprintf('(%5.3f, %5.3f)', lambda1(l1_max), lambda2(l2_max)), ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right')
hold off

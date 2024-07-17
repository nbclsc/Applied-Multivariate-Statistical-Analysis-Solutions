% Load and join data.
t1 = readtable(fullfile(data_folder, 'Table4.1.xlsx'), 'ReadVariableNames', false);
t1.Properties.VariableNames = {'OvenNo', 'RadiationClosed'};
t2 = readtable(fullfile(data_folder, 'Table4.5.xlsx'), 'ReadVariableNames', false);
t2.Properties.VariableNames = {'OvenNo', 'RadiationOpen'};
T = join(t1, t2, 'Keys', 'OvenNo');
X = table2array(T(:,2:3));
clear t1 t2 T

% Set up 2D lambda values.
lambda1 = linspace(0, 0.5, 200);
lambda2 = lambda1;
% Meshgrid 1st argument is x-axis (cols) and 2nd argument is y-axis (rows).
[l1, l2] = meshgrid(lambda1, lambda2);

% Transformed observations.
x1_lambda = (X(:,1).^(lambda1) - 1) ./ lambda1; % 42 x 200
x2_lambda = (X(:,2).^(lambda2) - 1) ./ lambda2; % 42 x 200
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
part1 = sum(log(X(:,1)));
part2 = sum(log(X(:,2)));
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
xlabel('\lambda_{1}')
ylabel('\lambda_{2}')
hold on
plot(lambda1(l1_max), lambda2(l2_max), 'o')
hold off

surf(l1, l2, l_lambda)
colormap default
hold on
plot3(lambda1(l1_max), lambda2(l2_max), a, ...
    'bo', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'b')
text(lambda1(l1_max), lambda2(l2_max), a+2, ...
    sprintf('(%5.3f, %5.3f)', lambda1(l1_max), lambda2(l2_max)), ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right')
hold off
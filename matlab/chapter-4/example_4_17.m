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

% Compute values for function to maximize.
l_lambda = zeros(length(lambda1), length(lambda2));
for i=1:length(lambda1)
    for j=1:length(lambda2)
        l_lambda(i,j) = -(height(X)/2)*log(det(cov(x1_lambda(:,i), x2_lambda(:,j)))) + (lambda1(i) - 1) * sum(log(X(:,1))) + (lambda2(j) - 1) * sum(log(X(:,2)));
    end
end

% Find the maximum.
[a, max_idx] = max(l_lambda(:));
% We used l1 on x-axis for meshgrid, which corresponds to the columns and
% l2 on the y-axis corresponding to the rows. In the matrix l_lambda, the
% rows correspond to l1 and columns to l2, that's why things are switched
% in the output of applying ind2sub below.
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
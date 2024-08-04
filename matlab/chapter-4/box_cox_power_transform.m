function [best_lambda] = box_cox_power_transform(x, col_label, lam_min, lam_max, output_file)
    % Univariate Box-Cox power transform to find the optimal power transformation attempting to make the input data normally distributed.
    
    % If there are 0 values need to add small amount so don't get Inf.
    if sum(x == 0) > 0
        x = x + 0.0000001;
    end

    lambda = linspace(lam_min, lam_max, 500);

    x_lambda = (x.^(lambda) - 1) ./ lambda;
    if sum(lambda==0)
        x_lambda(:, lambda==0) = log(x);
    end

    x_lambda_bar = (1/length(x))*sum(x_lambda);
    l_lambda = -(length(x)/2)*log((1/length(x))*sum((x_lambda - x_lambda_bar).^2)) + (lambda - 1)*sum(log(x));
    [max_lambda, argmax_lambda] = max(l_lambda);
    [min_lambda, ~] = min(l_lambda);
    best_lambda = lambda(argmax_lambda);
    
    p = plot(lambda, l_lambda);
    line([best_lambda best_lambda], [min_lambda max_lambda], 'Color', 'black', 'LineStyle', '--')
    hold on
    plot(best_lambda, max_lambda, 'o-','MarkerFaceColor','black','MarkerEdgeColor','black')
    text(best_lambda+0.1, max_lambda-2.5, sprintf('$\\hat{\\lambda}_{1} = %.4f$', best_lambda), 'Interpreter', 'latex')
    title(sprintf('Box-Cox Power Transformation Plot %s', col_label))
    xlabel('$\lambda$', 'Interpreter', 'latex')
    ylabel('$\ell(\lambda)$', 'Interpreter', 'latex');
    hold off
    
    saveas(p, output_file, 'png')
end
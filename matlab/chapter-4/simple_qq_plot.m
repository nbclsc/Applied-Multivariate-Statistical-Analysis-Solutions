function [p] = simple_qq_plot(xj, qj, col_name)
    % Create a Q-Q plot based on input data.
    title_text = sprintf('Q-Q Plot %s', col_name);
    p = scatter(qj, xj);
    title(title_text)
    xlabel("q_{(j)}")
    ylabel("x_(j)")
end
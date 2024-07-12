function [output, r_q] = corr_q_q(x)
    % Function to compute QQ correlation from (4-31).
    x_s = sort(x);
    N = length(x_s);
    prob_x = (linspace(1, N, N) - 0.5) / N;
    quant_x = icdf('Normal', prob_x, 0, 1);
    output = [x_s prob_x' quant_x'];

    % Compute the correlation coeeficient between the data and the
    % quantiles.
    r_q = corr(x_s, quant_x');
end
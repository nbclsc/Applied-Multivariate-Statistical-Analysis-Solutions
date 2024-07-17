function [output] = chi2_q_q(X)
    % Compute the values needed to construct a chi-square plot.
    % Detailed on page 184. For a chi square plot, plot column 1 vs
    % column 3.

    % The statistical distance values based on the input data (4-32).
    d2 = diag((X - mean(X))/cov(X)*(X - mean(X))');

    % Sort distances.
    d2_s = sort(d2);
    
    % Compute the probability values.
    N = length(d2_s);
    prob_d2 = (linspace(1, N, N) - 0.5) / N;
    
    % Compute the chi-squared quantiles from the probabilities.
    quant_d2 = icdf('Chi2', prob_d2, width(X));
    
    output = [d2_s prob_d2' quant_d2'];
end
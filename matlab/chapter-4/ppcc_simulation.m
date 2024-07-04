function [ppcc_values, critical_values] = ppcc_simulation(n, num_simulations, quant_levels)
    % Function to perform the PPCC simulation.
    
    % Set a 'default' number of simulations.
    if num_simulations == 0
        num_simulations = 100000;
    end

    % # Generate a random sample from a normal distribution.
    ppcc_values = zeros(num_simulations, 1);

    for i = 1:num_simulations

        % Start with a sample of size n from N(0,1).
        sample_data = normrnd(0, 1, n, 1);

        % Sort the sample data.
        ordered_data = sort(sample_data);

        % Compute the probability values for the sample.
        prob = (linspace(1, n, n) - 0.50) / n;

        % Theoretical quantiles from a normal distribution.
        theoretical_quantiles = icdf('Normal', prob, 0, 1)';

        % Compute the correlation coefficient.
        correlation_coefficient = corr(ordered_data, theoretical_quantiles);

        % Store the PPCC value.
        ppcc_values(i) = correlation_coefficient;
    end

    % Determine critical values for common significance levels.
    % Levels from the book in Table 4.2: [0.01, 0.05, 0.10]
    critical_values = quantile(ppcc_values, quant_levels);
end
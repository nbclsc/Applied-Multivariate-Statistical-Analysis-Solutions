ppcc_simulation <- function(n, num_simulations = 10000) {
  # Function to perform the PPCC simulation.
  
  # Initialize a vector to store PPCC values.
  ppcc_values <- numeric(num_simulations)
  
  for (i in 1:num_simulations) {
    # Generate a random sample from a normal distribution.
    sample_data <- rnorm(n)
    
    # Sort the sample data.
    ordered_data <- sort(sample_data)
    
    # Theoretical quantiles from a normal distribution.
    theoretical_quantiles <- qnorm(ppoints(n))
    
    # Compute the correlation coefficient.
    correlation_coefficient <- cor(ordered_data, theoretical_quantiles)
    
    # Store the PPCC value.
    ppcc_values[i] <- correlation_coefficient
  }
  
  # Determine critical values for common significance levels.
  critical_values <- quantile(ppcc_values, 0.1)
  
  # Return the distribution of PPCC values and critical values.
  list(
    ppcc_values = ppcc_values,
    critical_values = critical_values
  )
}

# Load the funtion that performs the simulation.
source("./ppcc_simulation.R")

# Perform the simulation for a sample size of 9.
n <- 9
num_simulations <- 1000000
result <- ppcc_simulation(n, num_simulations)

# Print critical value.
print(result$critical_values)

# Plot the distribution of PPCC values.
hist(result$ppcc_values,
     main = paste("Distribution of Probability Plot Correlation Coefficient values (n =", n, ")"),
     xlab = "Probability Plot Correlation Coefficient Simulated values",
     breaks = 50)
abline(v = result$critical_values, col = "red", lty = 2)
text(result$critical_values-.015,80000,latex2exp::TeX("$\\alpha = $0.10"))
# text(result$critical_values+.02,80000,expression(paste(alpha,round(result$critical_values,3))))
# text(result$critical_values-.02,80000,latex2exp::TeX(sprintf("$\\alpha =%3.4f$",0.10)))




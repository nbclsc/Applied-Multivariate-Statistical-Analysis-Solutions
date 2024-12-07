# Make sure the group (treatment) label column is a factor variable.
# Otherwise the output will not match.
df <- readxl::read_excel("../../data/Exercise6.8.xlsx")
    |> dplyr::mutate(Treat = as.factor(Treat))

# Use the car package for MANOVA. Personally, I'd rather use R for this
# than Python. Python results using statsmodels.multivariate.manova don't
# match the book (output here).
fit <- car::Manova(lm(cbind(resp1, resp2) ~ Treat, data = df))

# Output should match the corresponding Python notebook for Exercise 6.8 in
# parts (b) and (c).
summary(fit)

# fit$SSP  # B matrix
# fit$SSPE  # W matrix
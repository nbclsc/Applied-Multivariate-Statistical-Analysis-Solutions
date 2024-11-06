# Load the data file for Table 6.4.
# film_df <- readxl::read_excel("./Applied-Multivariate-Statistical-Analysis-Solutions/data/Table6.4.xlsx")
film_df <- readxl::read_excel("../../data/Table6.4.xlsx")

# Rename the columns to something shorter.
names(film_df) <- c('factor1','factor2','x1','x2','x3')
film_df <- film_df |> dplyr::mutate(dplyr::across(c(factor1, factor2), factor))

# Fit the MANOVA model.
fit <- manova(cbind(x1, x2, x3) ~ factor1 * factor2, data=film_df)
results <- summary(fit, test="Wilks")

print("The SSP matrices found in the ANOVA table above Panel 6.1 on page 319.")
results$SS["factor1"]
results$SS["factor2"]
results$SS["factor1:factor2"]
results$SS["Residuals"]

print("")
print("")
print("The ANOVA tables broken out for each variable.")
print("Found in the first part of the SAS output on pages 319-321.")
print("These are diagonal values from the matrices above.")
summary(aov(x1 ~ factor1 * factor2, data = film_df))
summary(aov(x2 ~ factor1 * factor2, data = film_df))
summary(aov(x3 ~ factor1 * factor2, data = film_df))

print("")
print("")
print("Wilk's lambda values.")
results

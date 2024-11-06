# Load the data file for Table 6.4.
# "./Applied-Multivariate-Statistical-Analysis-Solutions/data/Table6.4.xlsx"
film_df = readxl::read_excel("../../data/Table6.4.xlsx")

# Rename the columns to something shorter.
names(film_df) = c('factor1','factor2','x1','x2','x3')

# Fit the MANOVA model.
fit <- manova(cbind(x1, x2, x3) ~ factor1 * factor2, data=film_df)
results <- summary(manova(fit), test="Wilks")
results
# file.exists()

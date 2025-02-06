# Make sure the treatment column are factor variables.
# Otherwise the output will not match.
df <- readxl::read_excel("../../data/Exercise6.14.xlsx") |>
    dplyr::mutate(Factor1 = as.factor(Factor1),
                  Factor2 = as.factor(Factor2))
df <- readxl::read_excel("./data/Exercise6.14.xlsx") %>%
    dplyr::mutate(Factor1 = as.factor(Factor1),
                  Factor2 = as.factor(Factor2))


# Use the car package for MANOVA. Personally, I'd rather use R for this
# than Python. Python results using statsmodels.multivariate.manova don't
# match the book (output here).
fit_q_6_14 <- car::Manova(lm(cbind(x1, x2) ~ Factor1*Factor2, data = df))

summary(fit_q_6_14)

# fit$SSP  # B matrices
# fit$SSPE  # W matrix

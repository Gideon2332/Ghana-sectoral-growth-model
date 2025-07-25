# ---------------------------------------------
# 1. CLEAR ENVIRONMENT & LOAD LIBRARIES
# ---------------------------------------------
rm(list = ls())
cat("\014")

# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tseries)
library(urca)
library(knitr)
library(kableExtra)
library(dynamac)
library(vars)
library(lmtest)
library(strucchange)
library(moments)
library(psych)
library(cointReg)

# ---------------------------------------------
# 2. LOAD DATA
# ---------------------------------------------
data <- read.csv("C:/Users/PC/OneDrive/Desktop/New folder/Cletus_New/New folder/data_cleaned.csv")

# Create additional variables if needed
data$log_edu <- log(data$edu)
data$edu_PC <- data$log_edu * data$PoliticalCycle

# ---------------------------------------------
# 3. DESCRIPTIVE STATISTICS & HISTOGRAMS
# ---------------------------------------------
summary_vars <- data %>%
  select(Agriculture_GDP, Industry_GDP, Services_GDP,
         log_agric, log_fdi, log_gcf,
         FDI_percent_GDP, Capital_Formation, PoliticalCycle)

describe(summary_vars)
apply(summary_vars, 2, skewness)
apply(summary_vars, 2, kurtosis)

# Histograms
hist(data$log_agric, breaks = 10, main = "Histogram of log_agric", col = "lightblue")
hist(data$log_fdi, breaks = 10, main = "Histogram of log_fdi", col = "lightgreen")

# ---------------------------------------------
# 4. STATIONARITY TESTS
# ---------------------------------------------
vars_to_test <- data %>%
  select(log_agric, log_fdi, log_gcf, Agriculture_GDP, Industry_GDP, Services_GDP)

stationarity_tests <- function(x) {
  adf <- adf.test(x)$p.value
  pp_test <- ur.pp(x, type = "Z-tau", model = "constant")
  pp_stat <- pp_test@teststat[1]
  pp_pval <- pp_test@cval[1, "5pct"]
  pp_result <- ifelse(pp_stat < pp_pval, "< 5% (Stationary)", "> 5% (Non-stationary)")
  return(c(ADF = round(adf, 4),
           ADF_result = ifelse(adf < 0.05, "Stationary", "Non-stationary"),
           PP_stat = round(pp_stat, 3),
           PP_result = pp_result))
}

results_df <- as.data.frame(t(sapply(vars_to_test, stationarity_tests)))
results_df <- tibble::rownames_to_column(results_df, "Variable")

kable(results_df, caption = "Stationarity Test Results (ADF and PP)") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)

# ---------------------------------------------
# 5. BASELINE ARDL MODEL (AGRIC)
# ---------------------------------------------
model_agric <- dynardl(
  log_agric ~ log_fdi + log_gcf + PoliticalCycle,
  data = data,
  lags = list(log_agric = 1, log_fdi = 1, log_gcf = 1, PoliticalCycle = 1),
  ec = TRUE,
  simulate = FALSE
)
summary(model_agric)

# ---------------------------------------------
# 6. DISAGGREGATED ARDL WITH INTERACTIONS
# ---------------------------------------------
model_disagg <- dynardl(
  log_agric ~ log_edu + PoliticalCycle + edu_PC + log_fdi,
  data = data,
  lags = list(log_agric = 1, log_edu = 1, PoliticalCycle = 1, edu_PC = 1, log_fdi = 1),
  ec = TRUE
)
summary(model_disagg)

# ---------------------------------------------
# 7. ROBUSTNESS: FMOLS & DOLS
# ---------------------------------------------
y <- data$log_agric
X <- data[, c("log_edu", "PoliticalCycle", "edu_PC", "log_fdi")]
X <- as.matrix(na.omit(X))
y <- y[complete.cases(X)]

# FMOLS
fmols_model <- cointRegFM(y = y, x = X)
fmols_table <- data.frame(
  Variable = colnames(X),
  Estimate = fmols_model$theta,
  `Std. Error` = fmols_model$sd.theta,
  `t value` = fmols_model$t.theta,
  `p value` = fmols_model$p.theta
)
kable(fmols_table, caption = "FMOLS Estimation Results", digits = 4) %>%
  kable_styling()

# DOLS
dols_model <- cointRegD(y = y, x = X, p = 1, q = 1, bandwidth = 4)
dols_table <- data.frame(
  Variable = colnames(X),
  Estimate = dols_model$theta,
  `Std. Error` = dols_model$sd.theta,
  `t value` = dols_model$t.theta,
  `p value` = dols_model$p.theta
)
kable(dols_table, caption = "DOLS Estimation Results", digits = 4) %>%
  kable_styling()

# ---------------------------------------------
# 8. DIAGNOSTICS & STABILITY
# ---------------------------------------------
bgtest(model_disagg$model)
bptest(model_disagg$model)
resettest(model_disagg$model)
jarque.bera.test(residuals(model_disagg$model))

cusum_test <- efp(residuals(model_disagg$model) ~ 1, type = "OLS-CUSUM")
plot(cusum_test, main = "CUSUM Test - Disaggregated ARDL")

# ---------------------------------------------
# 9. SUBSAMPLE ANALYSIS (PRE- & POST-OIL)
# ---------------------------------------------
pre_oil <- subset(data, Year <= 2006)
post_oil <- subset(data, Year >= 2007)

model_pre <- dynardl(log_agric ~ log_edu + PoliticalCycle + edu_PC + log_fdi,
                     data = pre_oil,
                     lags = list(log_agric = 1, log_edu = 1, PoliticalCycle = 1, edu_PC = 1, log_fdi = 1))
model_post <- dynardl(log_agric ~ log_edu + PoliticalCycle + edu_PC + log_fdi,
                      data = post_oil,
                      lags = list(log_agric = 1, log_edu = 1, PoliticalCycle = 1, edu_PC = 1, log_fdi = 1))

summary(model_pre)
summary(model_post)

# ---------------------------------------------
# 10. SHOCK ANALYSIS (COVID + 2014 OIL CRASH)
# ---------------------------------------------
data$covid_dummy <- ifelse(data$Year %in% c(2020, 2021), 1, 0)
data$shock2014 <- ifelse(data$Year == 2014, 1, 0)

model_shocks <- dynardl(
  log_agric ~ log_edu + PoliticalCycle + edu_PC + log_fdi + covid_dummy + shock2014,
  data = data,
  lags = list(log_agric = 1, log_edu = 1, PoliticalCycle = 1, edu_PC = 1, log_fdi = 1),
  levels = c("covid_dummy", "shock2014"),
  model = "levels"
)
summary(model_shocks)

# ---------------------------------------------
# END
# ---------------------------------------------

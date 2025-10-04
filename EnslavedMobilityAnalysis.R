###################################################

# Enslaved Mobility Analysis (251003)

###################################################


library(dplyr)
library(tidyr)
library(vegan)
library(ppcor)
library(TraMineR)
library(tibble)
library(MASS)
library(pscl)
library(stringr)
library(cluster)
library(ggplot2)
library(lmtest)
library(broom)
library(readr)

# 1. Data

edges <- read.csv("https://raw.githubusercontent.com/JoonHwang-psu/EnslavedMobilityNetwork/refs/heads/main/WeightedEdgelist_JH.csv", 
                  stringsAsFactors = FALSE)

# 2. overall mobility diversity and entropy

# Count unique locations per person (overall mobility diversity)
location_count <- edges %>%
  group_by(Source) %>%
  summarise(unique_locations = n_distinct(Target)) %>%
  ungroup()

# Compute overall mobility entropy
wide_matrix <- edges %>%
  group_by(Source, Target) %>%
  summarise(weight = sum(Weight), .groups = "drop") %>%
  pivot_wider(names_from = Target, values_from = weight, values_fill = 0) %>%
  column_to_rownames("Source") %>%
  as.matrix()

entropy_values <- diversity(wide_matrix, index = "shannon")

# Convert entropy to a dataframe
entropy_df <- data.frame(
  Source = names(entropy_values),
  entropy = entropy_values
)

# Merge counts and entropy
mobility_stats <- location_count %>%
  left_join(entropy_df, by = "Source")

# Sum escape attempts
escape_totals <- edges %>%
  filter(Layer == "escape_attempt") %>%
  group_by(Source) %>%
  summarise(total_escape_attempts = sum(Weight), .groups = "drop")

mobility_stats <- mobility_stats %>%
  left_join(escape_totals, by = "Source") %>%
  mutate(total_escape_attempts = ifelse(is.na(total_escape_attempts), 0, total_escape_attempts))

# Get gender per person
gender_df <- edges %>%
  group_by(Source) %>%
  summarise(gender = first(gender))

# Merge with mobility_stats
mobility_stats <- mobility_stats %>%
  left_join(gender_df, by = "Source")

mobility_stats$gender <- as.factor(mobility_stats$gender)

# Create binary escape variable
mobility_stats$any_escape <- ifelse(mobility_stats$total_escape_attempts > 0, 1, 0)

sum(mobility_stats$any_escape)

# 2.1. mobility diversity analysis

# 2.1.1. logistic regression
logit_model <- glm(any_escape ~ unique_locations + gender, family = binomial, data = mobility_stats)
summary(logit_model)


# 2.1.2. negative binomial (NB) / Poisson regression
nb_model <- glm.nb(total_escape_attempts ~ unique_locations + gender, data = mobility_stats)
summary(nb_model)

poisson_model <- glm(total_escape_attempts ~ unique_locations + gender, 
                     data = mobility_stats, family = poisson)
summary(poisson_model)

# Poisson vs. NB
AIC(nb_model, poisson_model)
dispersion <- sum(residuals(poisson_model, type="pearson")^2) / df.residual(poisson_model)
print(dispersion) 
lrtest(poisson_model, nb_model)

# 2.1.3. zero-inflated negative binomial (ZINB) regression
zinb_model <- zeroinfl(total_escape_attempts ~ unique_locations + gender | unique_locations + gender, 
                       data = mobility_stats, dist = "negbin")
summary(zinb_model)

# ZINB vs. NB/Poisson
vuong(zinb_model, nb_model)
vuong(zinb_model, poisson_model)

# 2.1.4. calculating OR/IRR 

# For logistic regression (OR)
logit_or <- exp(coef(logit_model))
logit_ci <- exp(confint(logit_model))  # 95% CI for OR

# For Poisson regression (IRR)
poisson_irr <- exp(coef(poisson_model))
poisson_ci <- exp(confint(poisson_model))

# For Negative Binomial regression (IRR)
nb_irr <- exp(coef(nb_model))
nb_ci <- exp(confint(nb_model))

# For ZINB model
zinb_summary <- summary(zinb_model)

# Count part IRR
zinb_count_coefs <- zinb_summary$coefficients$count
zinb_count_irr <- exp(zinb_count_coefs[, "Estimate"])
zinb_count_ci <- exp(cbind(
  zinb_count_coefs[, "Estimate"] - 1.96 * zinb_count_coefs[, "Std. Error"],
  zinb_count_coefs[, "Estimate"] + 1.96 * zinb_count_coefs[, "Std. Error"]
))
colnames(zinb_count_ci) <- c("2.5 %", "97.5 %")

# Zero-inflation part OR
zinb_zero_coefs <- zinb_summary$coefficients$zero
zinb_zero_or <- exp(zinb_zero_coefs[, "Estimate"])
zinb_zero_ci <- exp(cbind(
  zinb_zero_coefs[, "Estimate"] - 1.96 * zinb_zero_coefs[, "Std. Error"],
  zinb_zero_coefs[, "Estimate"] + 1.96 * zinb_zero_coefs[, "Std. Error"]
))
colnames(zinb_zero_ci) <- c("2.5 %", "97.5 %")


cat("\n--- Logistic Regression (Odds Ratios) ---\n")
print(round(cbind(OR = logit_or, logit_ci), 3))

cat("\n--- Poisson Regression (Incidence Rate Ratios) ---\n")
print(round(cbind(IRR = poisson_irr, poisson_ci), 3))

cat("\n--- Negative Binomial Regression (Incidence Rate Ratios) ---\n")
print(round(cbind(IRR = nb_irr, nb_ci), 3))

cat("\n--- ZINB Model (Count part: IRR) ---\n")
print(round(cbind(IRR = zinb_count_irr, zinb_count_ci), 3))

cat("\n--- ZINB Model (Zero-inflation part: OR) ---\n")
print(round(cbind(OR = zinb_zero_or, zinb_zero_ci), 3))

# 2.2. mobility entropy

# 2.2.1. logistic regression
logit_model <- glm(any_escape ~ entropy + gender, 
                     family = binomial, data = mobility_stats)
summary(logit_model)

# 2.2.2. NB/Poisson regression
nb_model <- glm.nb(total_escape_attempts ~ entropy + gender, data = mobility_stats)
summary(nb_model)

poisson_model <- glm(total_escape_attempts ~ entropy + gender, 
                     data = mobility_stats, family = poisson)
summary(poisson_model)

# NB vs. Poisson
AIC(nb_model, poisson_model)
dispersion <- sum(residuals(poisson_model, type="pearson")^2) / df.residual(poisson_model)
print(dispersion) 
lrtest(poisson_model, nb_model)

# 2.2.3. ZINB regression
zinb_model <- zeroinfl(total_escape_attempts ~ gender + entropy | gender + entropy, 
                         data = mobility_stats, dist = "negbin")
summary(zinb_model)

# ZINB vs NB/Poisson
vuong(zinb_model, nb_model)
vuong(zinb_model, poisson_model)

# 2.2.4. calculating OR/IRR

# For logistic regression (OR)
logit_or <- exp(coef(logit_model))
logit_ci <- exp(confint(logit_model))  # 95% CI for OR

# For Poisson regression (IRR)
poisson_irr <- exp(coef(poisson_model))
poisson_ci <- exp(confint(poisson_model))

# For Negative Binomial regression (IRR)
nb_irr <- exp(coef(nb_model))
nb_ci <- exp(confint(nb_model))

# For ZINB model
zinb_summary <- summary(zinb_model)

# Count part IRR
zinb_count_coefs <- zinb_summary$coefficients$count
zinb_count_irr <- exp(zinb_count_coefs[, "Estimate"])
zinb_count_ci <- exp(cbind(
  zinb_count_coefs[, "Estimate"] - 1.96 * zinb_count_coefs[, "Std. Error"],
  zinb_count_coefs[, "Estimate"] + 1.96 * zinb_count_coefs[, "Std. Error"]
))
colnames(zinb_count_ci) <- c("2.5 %", "97.5 %")

# Zero-inflation part OR
zinb_zero_coefs <- zinb_summary$coefficients$zero
zinb_zero_or <- exp(zinb_zero_coefs[, "Estimate"])
zinb_zero_ci <- exp(cbind(
  zinb_zero_coefs[, "Estimate"] - 1.96 * zinb_zero_coefs[, "Std. Error"],
  zinb_zero_coefs[, "Estimate"] + 1.96 * zinb_zero_coefs[, "Std. Error"]
))
colnames(zinb_zero_ci) <- c("2.5 %", "97.5 %")


cat("\n--- Logistic Regression (Odds Ratios) ---\n")
print(round(cbind(OR = logit_or, logit_ci), 3))

cat("\n--- Poisson Regression (Incidence Rate Ratios) ---\n")
print(round(cbind(IRR = poisson_irr, poisson_ci), 3))

cat("\n--- Negative Binomial Regression (Incidence Rate Ratios) ---\n")
print(round(cbind(IRR = nb_irr, nb_ci), 3))

cat("\n--- ZINB Model (Count part: IRR) ---\n")
print(round(cbind(IRR = zinb_count_irr, zinb_count_ci), 3))

cat("\n--- ZINB Model (Zero-inflation part: OR) ---\n")
print(round(cbind(OR = zinb_zero_or, zinb_zero_ci), 3))

# 3. place-specific mobility

# calculating place-specific mobility
unique(edges$Target)
unique(edges$Layer)

edges <- edges %>%
  mutate(
    Target = str_trim(Target),                       # remove leading/trailing whitespace
    Target = str_to_lower(Target),                   # convert to lowercase
    Target = str_replace_all(Target, " ", "_"),      # replace spaces with underscores
    Target = str_replace_all(Target, ",", "")        # remove commas
  )

location_visits <- edges %>%
  filter(Layer != "escape_attempt") %>%
  group_by(Source, Target) %>%
  summarise(total_visits = sum(Weight), .groups = "drop") %>%
  pivot_wider(names_from = Target, values_from = total_visits, values_fill = 0)

mobility_stats <- left_join(mobility_stats, location_visits, by = "Source")

names(mobility_stats)

# 3.1. logistic regression
logit_model <- glm(any_escape ~ white_bluff + savannah + sassafras_field + gender,
                   family = binomial, data = mobility_stats)
summary(logit_model)

# 3.2. NB / Poisson regression
nb_model <- glm.nb(total_escape_attempts ~ white_bluff + savannah + sassafras_field + gender,
                   data = mobility_stats)
summary(nb_model)

poisson_model <- glm(total_escape_attempts ~ white_bluff + savannah + sassafras_field + gender,
                     data = mobility_stats, family = poisson)
summary(poisson_model)

# NB vs. Poisson
AIC(nb_model, poisson_model)
dispersion <- sum(residuals(poisson_model, type="pearson")^2) / df.residual(poisson_model)
print(dispersion) 
lrtest(poisson_model, nb_model)

# 3.3. ZINB regression
zinb_model <- zeroinfl(total_escape_attempts ~ white_bluff + savannah + sassafras_field + gender |
                        white_bluff + savannah + sassafras_field + gender,
                       data = mobility_stats, dist = "negbin")
summary(zinb_model)

# ZINB vs. NB/Poisson
vuong(zinb_model, nb_model)
vuong(zinb_model, poisson_model)

# calculating OR/IRR

# For logistic regression (OR)
logit_or <- exp(coef(logit_model))
logit_ci <- exp(confint(logit_model))  # 95% CI for OR

# For Poisson regression (IRR)
poisson_irr <- exp(coef(poisson_model))
poisson_ci <- exp(confint(poisson_model))

# For Negative Binomial regression (IRR)
nb_irr <- exp(coef(nb_model))
nb_ci <- exp(confint(nb_model))

# For ZINB model
zinb_summary <- summary(zinb_model)

# Count part IRR
zinb_count_coefs <- zinb_summary$coefficients$count
zinb_count_irr <- exp(zinb_count_coefs[, "Estimate"])
zinb_count_ci <- exp(cbind(
  zinb_count_coefs[, "Estimate"] - 1.96 * zinb_count_coefs[, "Std. Error"],
  zinb_count_coefs[, "Estimate"] + 1.96 * zinb_count_coefs[, "Std. Error"]
))
colnames(zinb_count_ci) <- c("2.5 %", "97.5 %")

# Zero-inflation part OR
zinb_zero_coefs <- zinb_summary$coefficients$zero
zinb_zero_or <- exp(zinb_zero_coefs[, "Estimate"])
zinb_zero_ci <- exp(cbind(
  zinb_zero_coefs[, "Estimate"] - 1.96 * zinb_zero_coefs[, "Std. Error"],
  zinb_zero_coefs[, "Estimate"] + 1.96 * zinb_zero_coefs[, "Std. Error"]
))
colnames(zinb_zero_ci) <- c("2.5 %", "97.5 %")


cat("\n--- Logistic Regression (Odds Ratios) ---\n")
print(round(cbind(OR = logit_or, logit_ci), 3))

cat("\n--- Poisson Regression (Incidence Rate Ratios) ---\n")
print(round(cbind(IRR = poisson_irr, poisson_ci), 3))

cat("\n--- Negative Binomial Regression (Incidence Rate Ratios) ---\n")
print(round(cbind(IRR = nb_irr, nb_ci), 3))

cat("\n--- ZINB Model (Count part: IRR) ---\n")
print(round(cbind(IRR = zinb_count_irr, zinb_count_ci), 3))

cat("\n--- ZINB Model (Zero-inflation part: OR) ---\n")
print(round(cbind(OR = zinb_zero_or, zinb_zero_ci), 3))

# 4. task-specific mobility

# calculating task-specific mobility
task_visits <- edges %>%
  filter(Layer != "escape_attempt") %>%
  group_by(Source, Layer) %>%
  summarise(total_visits = sum(Weight), .groups = "drop") %>%
  pivot_wider(names_from = Layer, values_from = total_visits, values_fill = 0)

mobility_stats <- left_join(mobility_stats, task_visits, by = "Source")

cor(mobility_stats[, c("sanctioned_travel", "unscheduled_time", "agri_labor", "healthcare", "maintenance")])
summary(mobility_stats[, c("sanctioned_travel", "unscheduled_time", "agri_labor", "healthcare", "maintenance")])

mobility_stats <- mobility_stats %>%
  mutate(any_unscheduled_time = as.factor(unscheduled_time > 0))

table(mobility_stats$any_unscheduled_time)

# 4.1. logistic regression
logit_model <- glm(any_escape ~ any_unscheduled_time + sanctioned_travel + agri_labor + healthcare + gender,
                   family = binomial, data = mobility_stats)
summary(logit_model)

# 4.2. NB / Poisson regression
nb_model <- glm.nb(total_escape_attempts ~ any_unscheduled_time + sanctioned_travel + agri_labor + healthcare + gender,
                   data = mobility_stats)
summary(nb_model)

poisson_model <- glm(total_escape_attempts ~ any_unscheduled_time + sanctioned_travel + agri_labor + healthcare + gender, 
                     data = mobility_stats, family = poisson)
summary(poisson_model)

# NB vs. Poisson
AIC(nb_model, poisson_model)
dispersion <- sum(residuals(poisson_model, type="pearson")^2) / df.residual(poisson_model)
print(dispersion) 
lrtest(poisson_model, nb_model)

# 4.3. ZINB regression
zinb_model <- zeroinfl(total_escape_attempts ~ any_unscheduled_time + sanctioned_travel + agri_labor + healthcare + gender |
                         any_unscheduled_time + sanctioned_travel + agri_labor + healthcare + gender,
                       data = mobility_stats, dist = "negbin")
summary(zinb_model)

# ZINB vs. NB / Poisson
vuong(zinb_model, nb_model)
vuong(zinb_model, poisson_model)

# 4.4. calculating OR/IRR

# For logistic regression (OR)
logit_or <- exp(coef(logit_model))
logit_ci <- exp(confint(logit_model))  # 95% CI for OR

# For Poisson regression (IRR)
poisson_irr <- exp(coef(poisson_model))
poisson_ci <- exp(confint(poisson_model))

# For Negative Binomial regression (IRR)
nb_irr <- exp(coef(nb_model))
nb_ci <- exp(confint(nb_model))

# For ZINB model
zinb_summary <- summary(zinb_model)

# Count part IRR
zinb_count_coefs <- zinb_summary$coefficients$count
zinb_count_irr <- exp(zinb_count_coefs[, "Estimate"])
zinb_count_ci <- exp(cbind(
  zinb_count_coefs[, "Estimate"] - 1.96 * zinb_count_coefs[, "Std. Error"],
  zinb_count_coefs[, "Estimate"] + 1.96 * zinb_count_coefs[, "Std. Error"]
))
colnames(zinb_count_ci) <- c("2.5 %", "97.5 %")

# Zero-inflation part OR
zinb_zero_coefs <- zinb_summary$coefficients$zero
zinb_zero_or <- exp(zinb_zero_coefs[, "Estimate"])
zinb_zero_ci <- exp(cbind(
  zinb_zero_coefs[, "Estimate"] - 1.96 * zinb_zero_coefs[, "Std. Error"],
  zinb_zero_coefs[, "Estimate"] + 1.96 * zinb_zero_coefs[, "Std. Error"]
))
colnames(zinb_zero_ci) <- c("2.5 %", "97.5 %")


cat("\n--- Logistic Regression (Odds Ratios) ---\n")
print(round(cbind(OR = logit_or, logit_ci), 3))

cat("\n--- Poisson Regression (Incidence Rate Ratios) ---\n")
print(round(cbind(IRR = poisson_irr, poisson_ci), 3))

cat("\n--- Negative Binomial Regression (Incidence Rate Ratios) ---\n")
print(round(cbind(IRR = nb_irr, nb_ci), 3))

cat("\n--- ZINB Model (Count part: IRR) ---\n")
print(round(cbind(IRR = zinb_count_irr, zinb_count_ci), 3))

cat("\n--- ZINB Model (Zero-inflation part: OR) ---\n")
print(round(cbind(OR = zinb_zero_or, zinb_zero_ci), 3))


# 5. place-task-specific models

edges <- edges %>%
  mutate(task_place = paste(Layer, Target, sep = "_"))

# Filter to exclude escape_attempt (optional, based on your focus)
edges_filtered <- edges %>%
  filter(Layer != "escape_attempt")

# Summarize total visits by person and task-place
task_place_matrix <- edges_filtered %>%
  group_by(Source, task_place) %>%
  summarise(total = sum(Weight), .groups = "drop") %>%
  pivot_wider(names_from = task_place, values_from = total, values_fill = 0)

mobility_stats <- mobility_stats %>%
  left_join(task_place_matrix, by = "Source")


# 5.1. Savanna

# 5.1.1. logistic regression
logit_model <- glm(any_escape ~ healthcare_savannah + sanctioned_travel_savannah + unscheduled_time_savannah + gender,
                   family = binomial, data = mobility_stats)
summary(logit_model)

# 5.1.2. NB/Poisson regression
nb_model <- glm.nb(total_escape_attempts ~ healthcare_savannah + sanctioned_travel_savannah + unscheduled_time_savannah + gender,
                   data = mobility_stats)
summary(nb_model)

poisson_model <- glm(total_escape_attempts ~ healthcare_savannah + sanctioned_travel_savannah + unscheduled_time_savannah + gender, 
                     data = mobility_stats, family = poisson)
summary(poisson_model)

# NB vs. Poisson
AIC(nb_model, poisson_model)
dispersion <- sum(residuals(poisson_model, type="pearson")^2) / df.residual(poisson_model)
print(dispersion) 
lrtest(poisson_model, nb_model)

# 5.1.3. ZINB regression
zinb_model <- zeroinfl(total_escape_attempts ~ gender + unscheduled_time_savannah + sanctioned_travel_savannah + healthcare_savannah|
                         gender + unscheduled_time_savannah + sanctioned_travel_savannah + healthcare_savannah,
                       data = mobility_stats, dist = "negbin")
summary(zinb_model)

# ZINB vs. NB/Poisson
vuong(zinb_model, nb_model)
vuong(zinb_model, poisson_model)

# 5.1.4. calculating OR/IRR 

# For logistic regression (OR)
logit_or <- exp(coef(logit_model))
logit_ci <- exp(confint(logit_model))  # 95% CI for OR

# For Poisson regression (IRR)
poisson_irr <- exp(coef(poisson_model))
poisson_ci <- exp(confint(poisson_model))

# For Negative Binomial regression (IRR)
nb_irr <- exp(coef(nb_model))
nb_ci <- exp(confint(nb_model))

nb_ci_profile <- exp(confint(nb_model))
print(round(nb_ci_profile, 3))

# For ZINB model
zinb_summary <- summary(zinb_model)

# Count part IRR
zinb_count_coefs <- zinb_summary$coefficients$count
zinb_count_irr <- exp(zinb_count_coefs[, "Estimate"])
zinb_count_ci <- exp(cbind(
  zinb_count_coefs[, "Estimate"] - 1.96 * zinb_count_coefs[, "Std. Error"],
  zinb_count_coefs[, "Estimate"] + 1.96 * zinb_count_coefs[, "Std. Error"]
))
colnames(zinb_count_ci) <- c("2.5 %", "97.5 %")

# Zero-inflation part OR
zinb_zero_coefs <- zinb_summary$coefficients$zero
zinb_zero_or <- exp(zinb_zero_coefs[, "Estimate"])
zinb_zero_ci <- exp(cbind(
  zinb_zero_coefs[, "Estimate"] - 1.96 * zinb_zero_coefs[, "Std. Error"],
  zinb_zero_coefs[, "Estimate"] + 1.96 * zinb_zero_coefs[, "Std. Error"]
))
colnames(zinb_zero_ci) <- c("2.5 %", "97.5 %")


cat("\n--- Logistic Regression (Odds Ratios) ---\n")
print(round(cbind(OR = logit_or, logit_ci), 3))

cat("\n--- Poisson Regression (Incidence Rate Ratios) ---\n")
print(round(cbind(IRR = poisson_irr, poisson_ci), 3))

cat("\n--- Negative Binomial Regression (Incidence Rate Ratios) ---\n")
print(round(cbind(IRR = nb_irr, nb_ci), 3))

cat("\n--- ZINB Model (Count part: IRR) ---\n")
print(round(cbind(IRR = zinb_count_irr, zinb_count_ci), 3))

cat("\n--- ZINB Model (Zero-inflation part: OR) ---\n")
print(round(cbind(OR = zinb_zero_or, zinb_zero_ci), 3))

# 5.2. White Bluff

# 5.2.1. logistic regression
logit_model <- glm(any_escape ~ healthcare_white_bluff + sanctioned_travel_white_bluff + unscheduled_time_white_bluff + gender,
                   family = binomial, data = mobility_stats)
summary(logit_model)

# 5.2.2. NB / Poisson regression
nb_model <- glm.nb(total_escape_attempts ~ healthcare_white_bluff + sanctioned_travel_white_bluff + unscheduled_time_white_bluff + gender,
                   data = mobility_stats)
summary(nb_model)

poisson_model <- glm(total_escape_attempts ~ healthcare_white_bluff + sanctioned_travel_white_bluff + unscheduled_time_white_bluff + gender, 
                     data = mobility_stats, family = poisson)
summary(poisson_model)

# NB vs. Poisson
AIC(nb_model, poisson_model)
dispersion <- sum(residuals(poisson_model, type="pearson")^2) / df.residual(poisson_model)
print(dispersion) 
lrtest(poisson_model, nb_model)

# 5.2.3. ZINB regression
zinb_model <- zeroinfl(total_escape_attempts ~ healthcare_white_bluff + sanctioned_travel_white_bluff + unscheduled_time_white_bluff + gender |
                         healthcare_white_bluff + sanctioned_travel_white_bluff + unscheduled_time_white_bluff + gender,
                       data = mobility_stats, dist = "negbin")
summary(zinb_model)

vuong(zinb_model, nb_model)

vuong(zinb_model, poisson_model)

# 5.2.4. calculating OR/IRR 

# For logistic regression (OR)
logit_or <- exp(coef(logit_model))
logit_ci <- exp(confint(logit_model))  # 95% CI for OR

# For Poisson regression (IRR)
poisson_irr <- exp(coef(poisson_model))
poisson_ci <- exp(confint(poisson_model))

# For Negative Binomial regression (IRR)
nb_irr <- exp(coef(nb_model))
nb_ci <- exp(confint(nb_model))

nb_ci_profile <- exp(confint(nb_model))
print(round(nb_ci_profile, 3))

# For ZINB model
zinb_summary <- summary(zinb_model)

# Count part IRR
zinb_count_coefs <- zinb_summary$coefficients$count
zinb_count_irr <- exp(zinb_count_coefs[, "Estimate"])
zinb_count_ci <- exp(cbind(
  zinb_count_coefs[, "Estimate"] - 1.96 * zinb_count_coefs[, "Std. Error"],
  zinb_count_coefs[, "Estimate"] + 1.96 * zinb_count_coefs[, "Std. Error"]
))
colnames(zinb_count_ci) <- c("2.5 %", "97.5 %")

# Zero-inflation part OR
zinb_zero_coefs <- zinb_summary$coefficients$zero
zinb_zero_or <- exp(zinb_zero_coefs[, "Estimate"])
zinb_zero_ci <- exp(cbind(
  zinb_zero_coefs[, "Estimate"] - 1.96 * zinb_zero_coefs[, "Std. Error"],
  zinb_zero_coefs[, "Estimate"] + 1.96 * zinb_zero_coefs[, "Std. Error"]
))
colnames(zinb_zero_ci) <- c("2.5 %", "97.5 %")


cat("\n--- Logistic Regression (Odds Ratios) ---\n")
print(round(cbind(OR = logit_or, logit_ci), 3))

cat("\n--- Poisson Regression (Incidence Rate Ratios) ---\n")
print(round(cbind(IRR = poisson_irr, poisson_ci), 3))

cat("\n--- Negative Binomial Regression (Incidence Rate Ratios) ---\n")
print(round(cbind(IRR = nb_irr, nb_ci), 3))

cat("\n--- ZINB Model (Count part: IRR) ---\n")
print(round(cbind(IRR = zinb_count_irr, zinb_count_ci), 3))

cat("\n--- ZINB Model (Zero-inflation part: OR) ---\n")
print(round(cbind(OR = zinb_zero_or, zinb_zero_ci), 3))

# 5.3. Plantation core

# 5.3.1. logistic regression
logit_model <- glm(any_escape ~ healthcare_plantation_core + sanctioned_travel_plantation_core + unscheduled_time_plantation_core + gender,
                   family = binomial, data = mobility_stats)
summary(logit_model)

# 5.3.2. NB / Poisson regression
nb_model <- glm.nb(total_escape_attempts ~ healthcare_plantation_core + sanctioned_travel_plantation_core + unscheduled_time_plantation_core + gender,
                   data = mobility_stats)
summary(nb_model)

poisson_model <- glm(total_escape_attempts ~ healthcare_plantation_core + sanctioned_travel_plantation_core + unscheduled_time_plantation_core + gender, 
                     data = mobility_stats, family = poisson)
summary(poisson_model)

# NB vs. Poisson
AIC(nb_model, poisson_model)
dispersion <- sum(residuals(poisson_model, type="pearson")^2) / df.residual(poisson_model)
print(dispersion) 
lrtest(poisson_model, nb_model)

# 5.3.3. ZINB regression
zinb_model <- zeroinfl(total_escape_attempts ~ healthcare_plantation_core + sanctioned_travel_plantation_core + unscheduled_time_plantation_core + gender |
                         healthcare_plantation_core + sanctioned_travel_plantation_core + unscheduled_time_plantation_core + gender,
                       data = mobility_stats, dist = "negbin")
summary(zinb_model)

# ZINB vs. NB / Poisson
vuong(zinb_model, nb_model)
vuong(zinb_model, poisson_model)

# 5.3.4. calculating OR / IRR

# For logistic regression (OR)
logit_or <- exp(coef(logit_model))
logit_ci <- exp(confint(logit_model))  # 95% CI for OR

# For Poisson regression (IRR)
poisson_irr <- exp(coef(poisson_model))
poisson_ci <- exp(confint(poisson_model))

# For Negative Binomial regression (IRR)
nb_irr <- exp(coef(nb_model))
nb_ci <- exp(confint(nb_model))

nb_ci_profile <- exp(confint(nb_model))
print(round(nb_ci_profile, 3))

# For ZINB model
zinb_summary <- summary(zinb_model)

# Count part IRR
zinb_count_coefs <- zinb_summary$coefficients$count
zinb_count_irr <- exp(zinb_count_coefs[, "Estimate"])
zinb_count_ci <- exp(cbind(
  zinb_count_coefs[, "Estimate"] - 1.96 * zinb_count_coefs[, "Std. Error"],
  zinb_count_coefs[, "Estimate"] + 1.96 * zinb_count_coefs[, "Std. Error"]
))
colnames(zinb_count_ci) <- c("2.5 %", "97.5 %")

# Zero-inflation part OR
zinb_zero_coefs <- zinb_summary$coefficients$zero
zinb_zero_or <- exp(zinb_zero_coefs[, "Estimate"])
zinb_zero_ci <- exp(cbind(
  zinb_zero_coefs[, "Estimate"] - 1.96 * zinb_zero_coefs[, "Std. Error"],
  zinb_zero_coefs[, "Estimate"] + 1.96 * zinb_zero_coefs[, "Std. Error"]
))
colnames(zinb_zero_ci) <- c("2.5 %", "97.5 %")

cat("\n--- Logistic Regression (Odds Ratios) ---\n")
print(round(cbind(OR = logit_or, logit_ci), 3))

cat("\n--- Poisson Regression (Incidence Rate Ratios) ---\n")
print(round(cbind(IRR = poisson_irr, poisson_ci), 3))

cat("\n--- Negative Binomial Regression (Incidence Rate Ratios) ---\n")
print(round(cbind(IRR = nb_irr, nb_ci), 3))

cat("\n--- ZINB Model (Count part: IRR) ---\n")
print(round(cbind(IRR = zinb_count_irr, zinb_count_ci), 3))

cat("\n--- ZINB Model (Zero-inflation part: OR) ---\n")
print(round(cbind(OR = zinb_zero_or, zinb_zero_ci), 3))


# 6. Task-specific mobility diversity

individual_task_places <- edges %>%
  filter(Layer != "escape_attempt") %>%
  group_by(Source, Layer) %>%
  summarise(unique_places = n_distinct(Target), .groups = "drop") %>%
  mutate(Layer = paste0("places_visited_", Layer)) %>%
  pivot_wider(names_from = Layer, values_from = unique_places, values_fill = 0)

mobility_stats <- left_join(mobility_stats, individual_task_places, by = "Source")

summary(mobility_stats[, c("places_visited_unscheduled_time", "places_visited_sanctioned_travel", "places_visited_agri_labor", "places_visited_healthcare")])

# 6.1. logistic regression
logit_model <- glm(any_escape ~ places_visited_unscheduled_time + places_visited_sanctioned_travel + places_visited_agri_labor + places_visited_healthcare + gender,
                   family = binomial, data = mobility_stats)
summary(logit_model)

# 6.2. NB / Poisson regression
nb_model <- glm.nb(total_escape_attempts ~ places_visited_unscheduled_time + places_visited_sanctioned_travel + places_visited_agri_labor + places_visited_healthcare + gender,
                   data = mobility_stats)
summary(nb_model)

poisson_model <- glm(total_escape_attempts ~ places_visited_unscheduled_time + places_visited_sanctioned_travel + places_visited_agri_labor + places_visited_healthcare + gender, 
                     data = mobility_stats, family = poisson)
summary(poisson_model)

# NB vs. Poisson
AIC(nb_model, poisson_model)
dispersion <- sum(residuals(poisson_model, type="pearson")^2) / df.residual(poisson_model)
print(dispersion) 
lrtest(poisson_model, nb_model)

# 6.3. ZINB regression
zinb_model <- zeroinfl(total_escape_attempts ~ places_visited_unscheduled_time + places_visited_sanctioned_travel + places_visited_agri_labor + places_visited_healthcare + gender |
                         places_visited_unscheduled_time + places_visited_sanctioned_travel + places_visited_agri_labor + places_visited_healthcare + gender,
                       data = mobility_stats, dist = "negbin")
summary(zinb_model)

# ZINB vs. NB / Poisson
vuong(zinb_model, nb_model)
vuong(zinb_model, poisson_model)

# 6.4. calculating OR / IRR

# For logistic regression (OR)
logit_or <- exp(coef(logit_model))
logit_ci <- exp(confint(logit_model))  # 95% CI for OR

# For Poisson regression (IRR)
poisson_irr <- exp(coef(poisson_model))
poisson_ci <- exp(confint(poisson_model))

# For Negative Binomial regression (IRR)
nb_irr <- exp(coef(nb_model))
nb_ci <- exp(confint(nb_model))

nb_ci_profile <- exp(confint(nb_model))
print(round(nb_ci_profile, 3))

# For ZINB model
zinb_summary <- summary(zinb_model)

# Count part IRR
zinb_count_coefs <- zinb_summary$coefficients$count
zinb_count_irr <- exp(zinb_count_coefs[, "Estimate"])
zinb_count_ci <- exp(cbind(
  zinb_count_coefs[, "Estimate"] - 1.96 * zinb_count_coefs[, "Std. Error"],
  zinb_count_coefs[, "Estimate"] + 1.96 * zinb_count_coefs[, "Std. Error"]
))
colnames(zinb_count_ci) <- c("2.5 %", "97.5 %")

# Zero-inflation part OR
zinb_zero_coefs <- zinb_summary$coefficients$zero
zinb_zero_or <- exp(zinb_zero_coefs[, "Estimate"])
zinb_zero_ci <- exp(cbind(
  zinb_zero_coefs[, "Estimate"] - 1.96 * zinb_zero_coefs[, "Std. Error"],
  zinb_zero_coefs[, "Estimate"] + 1.96 * zinb_zero_coefs[, "Std. Error"]
))
colnames(zinb_zero_ci) <- c("2.5 %", "97.5 %")

cat("\n--- Logistic Regression (Odds Ratios) ---\n")
print(round(cbind(OR = logit_or, logit_ci), 3))

cat("\n--- Poisson Regression (Incidence Rate Ratios) ---\n")
print(round(cbind(IRR = poisson_irr, poisson_ci), 3))

cat("\n--- Negative Binomial Regression (Incidence Rate Ratios) ---\n")
print(round(cbind(IRR = nb_irr, nb_ci), 3))

cat("\n--- ZINB Model (Count part: IRR) ---\n")
print(round(cbind(IRR = zinb_count_irr, zinb_count_ci), 3))

cat("\n--- ZINB Model (Zero-inflation part: OR) ---\n")
print(round(cbind(OR = zinb_zero_or, zinb_zero_ci), 3))


# 7. Task-specific mobility entropy

entropy_list <- list()

for (task in unique(edges$Layer)) {
  task_edges <- edges %>%
    filter(Layer == task) %>%
    group_by(Source, Target) %>%
    summarise(weight = sum(Weight), .groups = "drop") %>%
    pivot_wider(names_from = Target, values_from = weight, values_fill = 0)
  
  if (ncol(task_edges) <= 2) {
    cat("Skipping task:", task, "- only one location\n")
    next
  }
  
  entropy_matrix <- task_edges %>%
    column_to_rownames("Source") %>%
    as.matrix()
  
  entropy_values <- diversity(entropy_matrix, index = "shannon")
  
  entropy_df <- data.frame(Source = names(entropy_values),
                           entropy = entropy_values,
                           task = task)
  entropy_list[[task]] <- entropy_df
}

task_entropy_df <- do.call(rbind, entropy_list)

entropy_wide <- task_entropy_df %>%
  pivot_wider(
    names_from = task,
    values_from = entropy,
    names_prefix = "entropy_"
  )

mobility_stats <- mobility_stats %>%
  left_join(entropy_wide, by = "Source")

mobility_stats <- mobility_stats %>%
  mutate(across(starts_with("entropy_"), ~replace_na(., 0)))

# 7.1. logistic regression
logit_model <- glm(any_escape ~ entropy_unscheduled_time + entropy_sanctioned_travel + entropy_agri_labor + entropy_healthcare + gender,
                   family = binomial, data = mobility_stats)
summary(logit_model)

# 7.2. NB / Poisson regression
nb_model <- glm.nb(total_escape_attempts ~ entropy_unscheduled_time + entropy_sanctioned_travel + entropy_agri_labor + entropy_healthcare + gender,
                   data = mobility_stats)
summary(nb_model)

poisson_model <- glm(total_escape_attempts ~ entropy_unscheduled_time + entropy_sanctioned_travel + entropy_agri_labor + entropy_healthcare + gender, 
                     data = mobility_stats, family = poisson)
summary(poisson_model)

# NB vs. Poisson
AIC(nb_model, poisson_model)
dispersion <- sum(residuals(poisson_model, type="pearson")^2) / df.residual(poisson_model)
print(dispersion) 
lrtest(poisson_model, nb_model)

# 7.3. ZINB regression
zinb_model <- zeroinfl(total_escape_attempts ~ entropy_unscheduled_time + entropy_sanctioned_travel + entropy_agri_labor + entropy_healthcare + gender |
                         entropy_unscheduled_time + entropy_sanctioned_travel + entropy_agri_labor + entropy_healthcare + gender,
                       data = mobility_stats, dist = "negbin")
summary(zinb_model)

# ZINB vs. NB / Poisson
vuong(zinb_model, nb_model)
vuong(zinb_model, poisson_model)

# 7.4. calculating OR / IRR

# For logistic regression (OR)
logit_or <- exp(coef(logit_model))
logit_ci <- exp(confint(logit_model))  # 95% CI for OR

# For Poisson regression (IRR)
poisson_irr <- exp(coef(poisson_model))
poisson_ci <- exp(confint(poisson_model))

# For Negative Binomial regression (IRR)
nb_irr <- exp(coef(nb_model))
nb_ci <- exp(confint(nb_model))

nb_ci_profile <- exp(confint(nb_model))
print(round(nb_ci_profile, 3))

# For ZINB model
zinb_summary <- summary(zinb_model)

# Count part IRR
zinb_count_coefs <- zinb_summary$coefficients$count
zinb_count_irr <- exp(zinb_count_coefs[, "Estimate"])
zinb_count_ci <- exp(cbind(
  zinb_count_coefs[, "Estimate"] - 1.96 * zinb_count_coefs[, "Std. Error"],
  zinb_count_coefs[, "Estimate"] + 1.96 * zinb_count_coefs[, "Std. Error"]
))
colnames(zinb_count_ci) <- c("2.5 %", "97.5 %")

# Zero-inflation part OR
zinb_zero_coefs <- zinb_summary$coefficients$zero
zinb_zero_or <- exp(zinb_zero_coefs[, "Estimate"])
zinb_zero_ci <- exp(cbind(
  zinb_zero_coefs[, "Estimate"] - 1.96 * zinb_zero_coefs[, "Std. Error"],
  zinb_zero_coefs[, "Estimate"] + 1.96 * zinb_zero_coefs[, "Std. Error"]
))
colnames(zinb_zero_ci) <- c("2.5 %", "97.5 %")

cat("\n--- Logistic Regression (Odds Ratios) ---\n")
print(round(cbind(OR = logit_or, logit_ci), 3))

cat("\n--- Poisson Regression (Incidence Rate Ratios) ---\n")
print(round(cbind(IRR = poisson_irr, poisson_ci), 3))

cat("\n--- Negative Binomial Regression (Incidence Rate Ratios) ---\n")
print(round(cbind(IRR = nb_irr, nb_ci), 3))

cat("\n--- ZINB Model (Count part: IRR) ---\n")
print(round(cbind(IRR = zinb_count_irr, zinb_count_ci), 3))

cat("\n--- ZINB Model (Zero-inflation part: OR) ---\n")
print(round(cbind(OR = zinb_zero_or, zinb_zero_ci), 3))


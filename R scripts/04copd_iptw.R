# Author: Marleen Bokern
# Date: 03/2023
# Purpose: generate propensity scores and IPTW weights for COPD cohort. Fit Cox models for all three outcomes.

packages <- c("tidyverse", "MetBrewer", "arrow", "mice", "survey", "cobalt", "ggplot2", "WeightIt", "gtsummary", "flextable", "openxlsx", "tools")
installed_packages <- packages %in% rownames(installed.packages())

if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

lapply(packages, library, character.only = TRUE)

#########################################################
setwd(Datadir_copd)

#read in parquet dataset
df <- read_parquet("copd_wave1_60d.parquet")

palette <- met.brewer("Cassatt2")

#keep only people who are in one of the treatment groups, remove all unnecessary variables
subset_df <- df[!is.na(df$treatgroup),]

subset_df <- subset_df[, c("patid", "age_index", "gender", "death_date", "imd", "bmicat", "eth", "diabetes_present", "hypertension_present", "cvd_present", "allcancers_present", "past_asthma_present", "kidney_present", "immunosuppression_present", "flu_vacc_present", "pneumo_vacc_present", "exacerb_present", "smokstatus","timeinstudy1", "timeinstudy2", "timeinstudy3", "timeinstudy_death_any", "timeout1", "timeout2", "timeout3", "timeout_death_any", "pos_covid_test_present", "covid_hes_present", "covid_death_present", "any_death_present", "treat", "treatgroup", "missing_ons")]

#check missingness pattern
#md.pattern(subset_df, rotate.names = TRUE)

#check number of missings per variable in model
missing_counts <- colSums(is.na(subset_df[, c("age_index", "gender", "imd", "bmicat","eth", "diabetes_present", "hypertension_present", "cvd_present", "allcancers_present", "past_asthma_present", "kidney_present", "immunosuppression_present", "flu_vacc_present", "pneumo_vacc_present", "exacerb_present", "smokstatus")]))
missing_counts 

#set reference levels for ethnicity and bmi
subset_df$eth <- relevel(subset_df$eth, ref = "White")
subset_df$bmicat <- relevel(subset_df$bmicat, ref = "Normal (18.5-24.9)")

# Define the propensity score formula
ps_formula <- treatgroup ~ age_index + gender + eth + imd + bmicat + diabetes_present + hypertension_present + cvd_present + allcancers_present + past_asthma_present + kidney_present + immunosuppression_present + flu_vacc_present + pneumo_vacc_present + exacerb_present + smokstatus

# Generating propensity scores as specified in Austin (2011)
ps1 <- glm(formula = ps_formula,
                    data = subset_df,
                    family = "binomial") #for logistic regression

subset_df$ps <- ps1$fitted.values

# Create a summary of the model
model_summary <- summary(ps1)$coefficients
model_summary <- as.data.frame(model_summary)

model_summary$OR <- exp(model_summary[, "Estimate"])

z_value <- qnorm(0.975) # 0.975 for a two-sided 95% confidence interval
model_summary$lower_ci <- model_summary$OR - z_value * model_summary[, "Std. Error"]
model_summary$upper_ci <- model_summary$OR + z_value * model_summary[, "Std. Error"]

model_summary <- model_summary[, -3]
model_summary$`p-value` <- ifelse(model_summary$`Pr(>|z|)` < 0.001, "< 0.001", sprintf("%.4f", model_summary$`Pr(>|z|)`))
model_summary <- model_summary[, -3]

colnames(model_summary)[colnames(model_summary) == 'lower_ci'] <- 'Lower CI'
colnames(model_summary)[colnames(model_summary) == 'upper_ci'] <- 'Upper CI'

model_rows <- c("Intercept", "Age at baseline", "Gender: Female", "Ethnicity: South Asian", "Ethnicity: Black", "Ethnicity: Mixed", "Ethnicity: Unknown", "IMD 2", "IMD 3", "IMD 4", "IMD 5", "IMD: Missing",  "BMI: Underweight (< 18.5)", "BMI: Overweight (25-29.9)", "BMI: Obese (>=30)", "Diabetes", "Hypertension", "CVD", "Cancer", "Past asthma", "Chronic Kidney Disease", "Immunosuppression", "Influenza Vaccine (past 12 months)", "Pneumococcal Vaccine (past 5 years)", "COPD exacerbation (past 12 months)", "Former smoking")

#make new variable with row names
model_summary$Variable <- model_rows

#move Variable to first column
model_summary <- model_summary[,c(7,1:6)]

file_path <- file.path(Graphdir, "iptw_diagnostics", "ps_model_results.xlsx")

# Extract the coefficients from the summary and write them to Excel
write.xlsx(model_summary, file_path, rowNames = FALSE)
rm("model_summary")

# Create a table of the propensity score model
tbl_ps <- ps1 %>%
  tbl_regression(exponentiate = TRUE,
                 label = list(
                   age_index = "Age",
                   gender = "Gender",
                   eth = "Ethnicity",
                   imd = "Index of Multiple Deprivation",
                   bmicat = "BMI",
                   diabetes_present = "Diabetes",
                   hypertension_present = "Hypertension",
                   cvd_present = "Cardiovascular disease",
                   allcancers_present = "Cancer",
                   past_asthma_present = "Past asthma",
                   kidney_present = "Kidney impairment",
                   immunosuppression_present = "Immunosuppression",
                   smokstatus = "Smoking",
                   flu_vacc_present = "Influenza vaccine",
                   pneumo_vacc_present = "Pneumococcal vaccine",
                   exacerb_present = "Any exacerbation in past 12 months")) %>%
  add_n(location = "level")

file_path <- paste0(Tables, "copd_ps_model_death.docx")

# Check if the file exists
if (file.exists(file_path)) {
  # If the file exists, delete it
  file.remove(file_path)
}

tbl_ps %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = file_path, align = "left")

#save min and max propensity scores in each treatment group, save limits of area of common support
ps_trim <- subset_df %>% 
  dplyr::select(treat, ps) %>% 
  group_by(treat) %>% 
  summarise(min = min(ps), max = max(ps)) %>% 
  ungroup() %>% 
  summarise(min = max(min), max = min(max))

#exclude people outside of common support
subset_df <- subset_df %>% 
  filter(ps >= ps_trim$min & ps <= ps_trim$max)



# IPTW --------------------------------------------------------------------
#generate IPTW weights

# generate ATE weights
ate_stabilised <- weightit(formula = ps_formula,
                           data = subset_df,
                           method = "glm",
                           estimand = "ATE",
                           stabilize = TRUE)

ate_unstabilised <- weightit(formula = ps_formula,
                             data = subset_df,
                             method = "glm",
                             estimand = "ATE",
                             stabilize = FALSE)

# generate ATT weights
att_unstabilised <- weightit(formula = ps_formula,
                             data = subset_df,
                             method = "glm",
                             estimand = "ATT",
                             stabilize = FALSE,
                             focal = "1")


# add ATE and ATT weights to dataset
subset_df$ate_weight_stab <- ate_stabilised[[1]]
subset_df$ate_weight_unstab <- ate_unstabilised[[1]]
subset_df$att_weight_unstab <- att_unstabilised[[1]]

#save as parquet
write_parquet(subset_df, "copd_wave1_60d_iptw.parquet")

weight_df <- subset_df %>% dplyr::select(patid, ate_weight_stab, ate_weight_unstab, att_weight_unstab)

write_parquet(weight_df, "copd_wave1_60d_iptw_weights.parquet")

# IPTW Diagnostics --------------------------------------------------------
# Density plot of propensity scores
plot <- ggplot(subset_df, aes(x = ps, color = factor(treatgroup))) +
   geom_density(size = 0.8) +
  labs(x = "Propensity Score",
       y = "Density",
       title = "Propensity score distribution (unweighted)",
       color = NULL) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c(palette[4], palette[9]), labels = c("LABA/LAMA", "ICS/LABA")) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(fill = c(palette[4], palette[9]))))

file_path <- file.path(Graphdir, "iptw_diagnostics", "ps_dist_unweighted.png")
ggsave(file_path, plot, width = 8, height = 4)

# Function to generate density plot of weighted propensity scores
generate_density_plot <- function(data, weight, title, file_name) {
  plot <- ggplot(data, aes(x = !!rlang::sym(weight), color = factor(treatgroup))) + 
    geom_density(aes(weight = !!rlang::sym(weight)), size = 0.8) +
    labs(x = "Weight",
         y = "Density",
         title = title,
         color = NULL) +
    scale_x_continuous(limits = c(0, 3)) +
    scale_color_manual(values = c(palette[4], palette[9]), labels = c("LABA/LAMA", "ICS/LABA")) +
    theme_minimal() +
    guides(color = guide_legend(override.aes = list(fill = c(palette[4], palette[9]))))
  
  file_path <- file.path(Graphdir, "iptw_diagnostics", file_name)
  ggsave(file_path, plot, width = 8, height = 4)
  
  return(plot)
}

# Generate density plots
generate_density_plot(subset_df, "ate_weight_unstab", "Distribution of unstabilised ATE Weights by Treatment Group", "dist_ate_unstab.png")
generate_density_plot(subset_df, "ate_weight_stab", "Distribution of stabilised ATE Weights Propensity Scores by Treatment Group", "dist_ate_stab.png")
generate_density_plot(subset_df, "att_weight_unstab", "Distribution of ATT Weights by Treatment Group", "dist_att_unstab.png")

# Function to generate density plot of weighted propensity scores
generate_density_plot <- function(data, weight, title, file_name) {
  plot <- ggplot(data, aes(x = ps, color = factor(treatgroup))) + 
    geom_density(aes(weight = !!rlang::sym(weight)), size = 0.8) +
    labs(x = "Propensity Score",
         y = "Density",
         title = title,
         color = NULL) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_color_manual(values = c(palette[4], palette[9]), labels = c("LABA/LAMA", "ICS/LABA")) +
    theme_minimal() +
    guides(color = guide_legend(override.aes = list(fill = c(palette[4], palette[9]))))
  
  file_path <- file.path(Graphdir, "iptw_diagnostics", file_name)
  ggsave(file_path, plot, width = 8, height = 4)
  
  return(plot)
}

# Generate density plots
generate_density_plot(subset_df, "ate_weight_unstab", "Distribution of unstabilised ATE-Weighted Propensity Scores by Treatment Group", "ps_dist_ate_unstab.png")
generate_density_plot(subset_df, "ate_weight_stab", "Distribution of stabilised ATE-Weighted Propensity Scores by Treatment Group", "ps_dist_ate_stab.png")
generate_density_plot(subset_df, "att_weight_unstab", "Distribution of unstabilised ATT-Weighted Propensity Scores by Treatment Group", "ps_dist_att_unstab.png")


#summarise mean and median ate and att weights by treatment group, with IQR and range
weight_summary <- subset_df %>% 
  group_by(treat) %>% 
  summarise(mean_att_weight_unstab = mean(att_weight_unstab, na.rm = TRUE),
            median_att_weight_unstab = median(att_weight_unstab, na.rm = TRUE),
            p25_att_weight_unstab = quantile(att_weight_unstab, 0.25, na.rm = TRUE),
            p75_att_weight_unstab = quantile(att_weight_unstab, 0.75, na.rm = TRUE),
            min_att_weight_unstab = min(att_weight_unstab, na.rm = TRUE),
            max_att_weight_unstab = max(att_weight_unstab, na.rm = TRUE),
            mean_ate_weight_stab = mean(ate_weight_stab, na.rm = TRUE),
            median_ate_weight_stab = median(ate_weight_stab, na.rm = TRUE),
            p25_ate_weight_stab = quantile(ate_weight_stab, 0.25, na.rm = TRUE),
            p75_ate_weight_stab = quantile(ate_weight_stab, 0.75, na.rm = TRUE),
            min_ate_weight_stab = min(ate_weight_stab, na.rm = TRUE),
            max_ate_weight_stab = max(ate_weight_stab, na.rm = TRUE),
            mean_ate_weight_unstab = mean(ate_weight_unstab, na.rm = TRUE),
            median_ate_weight_unstab = median(ate_weight_unstab, na.rm = TRUE),
            p25_ate_weight_unstab = quantile(ate_weight_unstab, 0.25, na.rm = TRUE),
            p75_ate_weight_unstab = quantile(ate_weight_unstab, 0.75, na.rm = TRUE),
            min_ate_weight_unstab = min(ate_weight_unstab, na.rm = TRUE),
            max_ate_weight_unstab = max(ate_weight_unstab, na.rm = TRUE))

#export this to excel
file_path <- file.path(Graphdir,  "iptw_diagnostics", "weights_summary.xlsx")
write.xlsx(weight_summary, file_path, rowNames = FALSE)

covariates <- c("exacerb_present",  "pneumo_vacc_present", "flu_vacc_present", "immunosuppression_present", "kidney_present", "past_asthma_present", "allcancers_present", "cvd_present", "hypertension_present", "diabetes_present", "smokstatus", "imd", "eth", "bmicat", "gender", "age_index")

#assess SMDs
data_ate_unstab <- subset_df[c(covariates, "treatgroup", "ate_weight_unstab")]
data_ate_stab <- subset_df[c(covariates, "treatgroup", "ate_weight_stab")]
data_att_unstab <- subset_df[c(covariates, "treatgroup", "att_weight_unstab")]

#generate bal.tab tables
ate_weighted_table_unstab <- bal.tab(data_ate_unstab, treat = data_ate_unstab$treatgroup, weights = data_ate_unstab$ate_weight_unstab)[[1]]
ate_weighted_table_stab <- bal.tab(data_ate_stab, treat = data_ate_stab$treatgroup, weights = data_ate_stab$ate_weight_stab)[[1]]
att_weighted_table_unstab <- bal.tab(data_att_unstab, treat = data_att_unstab$treatgroup, weights = data_att_unstab$att_weight_unstab)[[1]]
bal_unweighted_table <- bal.tab(data_ate_unstab, treat = data_ate_unstab$treatgroup)[[1]]

#remove the treatgroup and weight rows from the table
ate_weighted_table_unstab <- ate_weighted_table_unstab[!rownames(ate_weighted_table_unstab) %in% c("treatgroup", "ate_weight_unstab"), ]
ate_weighted_table_stab <- ate_weighted_table_stab[!rownames(ate_weighted_table_stab) %in% c("treatgroup", "ate_weight_stab"), ]
att_weighted_table_unstab <- att_weighted_table_unstab[!rownames(att_weighted_table_unstab) %in% c("treatgroup", "att_weight_unstab"), ]
bal_unweighted_table <- bal_unweighted_table[!rownames(bal_unweighted_table) %in% c("treatgroup", "ate_weight_unstab"), ]

#check that rownames are the same
all.equal(rownames(ate_weighted_table_unstab), rownames(ate_weighted_table_stab), rownames(att_weighted_table_unstab), rownames(bal_unweighted_table))

# generate a dataframe with the SMDs, unweighted, ATT and ATE weighted
smd_data <- data.frame(variable = rownames(ate_weighted_table_unstab),
                       unweighted = bal_unweighted_table$Diff.Un,
                       "ATT (unstabilised)" = att_weighted_table_unstab$Diff.Adj,
                       "ATE (unstabilised)" = ate_weighted_table_unstab$Diff.Adj,
                       "ATE (stabilised)" = ate_weighted_table_stab$Diff.Adj)

#export this to excel
file_path <- file.path(Graphdir,  "iptw_diagnostics", "smd_table.xlsx")
write.xlsx(smd_data, file_path, rowNames = FALSE)

new_names <- c("variable", "unweighted", "att_unstab", "ate_unstab", "ate_stab")
colnames(smd_data) <- new_names

#generate absolute SMDs
smd_data$abs_smd_unweighted <- abs(smd_data$unweighted)
smd_data$abs_smd_att_unstab <- abs(smd_data$att_unstab)
smd_data$abs_smd_ate_unstab <- abs(smd_data$ate_unstab)
smd_data$abs_smd_ate_stab <- abs(smd_data$ate_stab)


# Create plot of SMDs
ggplot(smd_data, aes(x = abs_smd_ate_stab, y = variable)) +
  geom_point(aes(shape = "ATE", fill = "ATE"), color = palette[9], size = 3) +
  geom_point(data = smd_data, aes(x = abs_smd_att_unstab, shape = "ATT", fill = "ATT"), color = palette[4], size = 2) +
  geom_point(data = smd_data, aes(x = abs_smd_unweighted, shape = "unweighted", fill = "unweighted"), color = palette[7], size = 2) +
  scale_shape_manual(values = c("ATE" = 23, "ATT" = 21, "unweighted" = 24)) +
  scale_fill_manual(values = c("ATE" = palette[9], "ATT" = palette[4], "unweighted" = palette[7])) +
  labs(title = "Comparison of Absolute SMDs",
       x = "Absolute Standardized Mean Difference (SMD)",
       y = "Variable") +
  scale_y_discrete(labels = c("smokstatus_ex-smoker" = "Former smoking", "pneumo_vacc_present_Yes" = "Pneumococcal vaccine (past 5 years)", "kidney_present_Yes" = "Chronic kidney disease", "immunosuppression_present_Yes" = "Immunosuppression", "imd_Missing" = "Missing IMD", "imd_5" = "IMD 5", "imd_4" = "IMD 4", "imd_3" = "IMD 3", "imd_2" = "IMD 2", "imd_1" = "IMD 1", "hypertension_present_Yes" = "Hypertension", "gender_Female" = "Female gender", "flu_vacc_present_Yes" = "Influenza vaccine (past year)", "exacerb_present_Yes" = "COPD exacerbation (past 12 months)", "diabetes_present_Yes" = "Diabetes", "cvd_present_Yes" = "Cardiovascular disease", "eth_White" = "Ethnicity: White", "eth_Unknown" = "Ethnicity: Unknown", "eth_South Asian" = "Ethnicity: South Asian", "eth_Mixed" = "Ethnicity: Mixed", "eth_Black" = "Ethnicity: Black", "bmicat_Underweight (<18.5)" = "BMI: Underweight (< 18.5)", "bmicat_Overweight (25-29.9)" = "BMI: Overweight (25-29.9)", "bmicat_Obese (>=30)" = "BMI: Obese (>=30)", "bmicat_Normal (18.5-24.9)" = "BMI: Normal (18.5-24.9)", "past_asthma_present_Yes" = "Past asthma", "allcancers_present_Yes" = "Cancer", "age_index" = "Age at index")) + 
  theme_minimal() +
  guides(fill = guide_legend(override.aes = list(color = c(palette[9], palette[4], palette[7])))) +
  labs(fill = NULL, shape = NULL)

file_path <- file.path(Graphdir, "iptw_diagnostics", "SMD_plot.png")
ggsave(file_path, width = 8, height = 4)

# Specify the weight variables
weight_vars <- c("ate_weight_stab", "ate_weight_unstab", "att_weight_unstab")
#recode smokstatus to Current smoking, Former smoking
subset_df$smokstatus <- recode(subset_df$smokstatus, "current smoker" = "Current smoking", "ex-smoker" = "Former smoking")

# Loop through each weight variable
for (weight_var in weight_vars) {
  # Subset the data and select the columns
  tab_data <- subset_df %>% dplyr::select(c("age_index", "gender", "bmicat", "eth", "imd", "diabetes_present", "hypertension_present","cvd_present", "allcancers_present", "past_asthma_present", "kidney_present", "immunosuppression_present", "flu_vacc_present",  "pneumo_vacc_present",  "exacerb_present", "smokstatus", "treat", weight_var))

  # Create the survey design
  design <- svydesign(ids = ~1, weights = as.formula(paste0("~", weight_var)), data = tab_data)
  
  # Generate the table
  tab <- tbl_svysummary(data = design,
                        by = treat,
                        label = list(age_index ~ "Age",
                                     gender ~ "Gender",
                                     eth ~ "Ethnicity",
                                     bmicat ~ "BMI",
                                     imd ~ "Index of Multiple Deprivation",
                                     diabetes_present ~ "Diabetes",
                                     hypertension_present ~ "Hypertension",
                                     cvd_present ~ "Cardiovascular disease",
                                     allcancers_present ~ "Cancer",
                                     past_asthma_present ~ "Past asthma",
                                     kidney_present ~ "Kidney impairment",
                                     immunosuppression_present ~ "Immunosuppression",
                                     smokstatus ~ "Smoking",
                                     flu_vacc_present ~ "Influenza vaccine",
                                     pneumo_vacc_present ~ "Pneumococcal vaccine",
                                     exacerb_present ~ "Any exacerbation in past 12 months",
                                     !!sym(weight_var) ~ paste0("Weight (", weight_var, ")")),
                        percent = "column",
                        digits = all_continuous() ~ 2,
                        missing = "ifany",
                        missing_text = "Missing",
                        statistic = list(
                          all_continuous() ~ "{mean} ({sd})",
                          age_index ~ c("{mean} ({sd})", "{median}  \n ({p25}-{p75})"),
                          all_categorical() ~ "{n} ({p}%)"),
                        type = list(
                          c(age_index) ~ "continuous2"))  %>%
    modify_header(label ~ "", all_stat_cols() ~ "**{level}**  \n N = {n}")  %>%
    modify_caption("Patient Characteristics") %>%
    modify_column_alignment(columns = c(stat_1, stat_2), align = "right") %>%
    italicize_levels()

  # Export the table to Word
  file_path <- file.path(Graphdir, "iptw_diagnostics", paste0("baseline_table_", weight_var, ".docx"))

  # Check if the file exists
  if (file.exists(file_path)) {
    # If the file exists, delete it
    file.remove(file_path)
  }

  tab %>%
    as_flex_table() %>%
    flextable::save_as_docx(path = file_path, align = "left")
}

# Create an empty summary dataframe
iptw_dist <- data.frame(
  Treat = character(),
  Mean = numeric(),
  SD = numeric(),
  Median = numeric(),
  Q25 = numeric(),
  Q75 = numeric(),
  min = numeric(),
  max = numeric())

# Specify the weight variables
weight_vars <- c("ate_weight_stab", "ate_weight_unstab", "att_weight_unstab")

# Create a function to plot IPTW weight distributions
plot_weight_distribution <- function(weights) {
  ggplot(data.frame(Weights = weights), aes(x = Weights)) +
    geom_histogram(binwidth = 0.1, fill = palette[3], color = "black") +
    labs(x = "IPTW Weights", y = "Frequency", title = "Distribution of IPTW Weights")
}

# Loop through each weight variable
for (weight_var in weight_vars) {
  
  # Loop through each treatment level
  for (treat_level in unique(subset_df$treat)) {
    weights_subset <- unlist(subset_df[subset_df$treat == treat_level, weight_var])
    
    # Calculate statistics
    mean_weight <- mean(weights_subset, na.rm = TRUE)
    sd_weight <- sd(weights_subset, na.rm = TRUE)
    median_weight <- median(weights_subset, na.rm = TRUE)
    q25 <- quantile(weights_subset, 0.25, na.rm = TRUE)
    q75 <- quantile(weights_subset, 0.75, na.rm = TRUE)
    min_weight <- min(weights_subset, na.rm = TRUE)
    max_weight <- max(weights_subset, na.rm = TRUE)
    
    # Add statistics to the summary dataframe
    iptw_dist <- rbind(iptw_dist, data.frame(
      Treat = treat_level,
      Weight_Var = weight_var,
      Mean = mean_weight,
      SD = sd_weight,
      Median = median_weight,
      Q25 = q25,
      Q75 = q75,
      min = min_weight,
      max = max_weight))
    
    # Print statistics
    cat("Treatment:", treat_level, "\n")
    cat("Weight Variable:", weight_var, "\n")
    cat("Mean:", mean_weight, "\n")
    cat("Standard Deviation:", sd_weight, "\n")
    cat("Median:", median_weight, "\n")
    cat("IQR (25th percentile):", q25, "\n")
    cat("IQR (75th percentile):", q75, "\n\n")
    
    # Call the plot function with the correct argument name
    print(plot_weight_distribution(weights_subset))
    
    # Remove temporary variables
    rm(list = c("weights_subset", "mean_weight", "sd_weight", "median_weight", "q25", "q75", "min_weight", "max_weight"))
  }
  
  # Calculate overall statistics
  weights <- unlist(subset_df[, weight_var])
  
  # Calculate statistics
  mean_weight <- mean(weights, na.rm = TRUE)
  sd_weight <- sd(weights, na.rm = TRUE)
  median_weight <- median(weights, na.rm = TRUE)
  q25 <- quantile(weights, 0.25, na.rm = TRUE)
  q75 <- quantile(weights, 0.75, na.rm = TRUE)
  min_weight <- min(weights, na.rm = TRUE)
  max_weight <- max(weights, na.rm = TRUE)
  
  # Add statistics to the summary dataframe
  iptw_dist <- rbind(iptw_dist, data.frame(
    Treat = "Overall",
    Weight_Var = weight_var,
    Mean = mean_weight,
    SD = sd_weight,
    Median = median_weight,
    Q25 = q25,
    Q75 = q75,
    min = min_weight,
    max = max_weight))
  
  # Print statistics
  cat("Treatment: Overall\n")
  cat("Weight Variable:", weight_var, "\n")
  cat("Mean:", mean_weight, "\n")
  cat("Standard Deviation:", sd_weight, "\n")
  cat("Median:", median_weight, "\n")
  cat("IQR (25th percentile):", q25, "\n")
  cat("IQR (75th percentile):", q75, "\n\n")
  
  # Call the plot function with the correct argument name
  #print(plot_weight_distribution(weights))
  
  # Remove temporary variables
  rm(list = c("weights", "mean_weight", "sd_weight", "median_weight", "q25", "q75", "min_weight", "max_weight"))
}
# Print the summary statistics
#print(iptw_dist)

# check ps model ----------------------------------------------------------

# #check of PS model: #restrict population to those with PS ~0.7 to see why the dip in propensity score occurs.
# 
# #restrict population to those with PS between 0.65 and 0.75
# subset_df <- subset_df %>% 
#   filter(ps >= 0.68 & ps <= 0.73)
# 
# tab1 <- tbl_summary(subset_df %>% dplyr::select(age_index, gender, bmicat, eth, smok,  diabetes_present, hypertension_present, cvd_present, allcancers_present, asthma_present, kidney_present, immunosuppression_present, flu_vacc_present, pneumo_vacc_present,  exacerb_present, pos_covid_test_present, covid_hes_present, covid_death_present, treat),
#                     by = treat,
#                     label = list(age_index ~ "Age",
#                                  gender ~ "Gender",
#                                  eth ~ "Ethnicity",
#                                  bmicat ~ "BMI",
#                                  diabetes_present ~ "Diabetes",
#                                  hypertension_present ~ "Hypertension",
#                                  cvd_present ~ "Cardiovascular disease",
#                                  allcancers_present ~ "Cancer",
#                                  asthma_present ~ "Past asthma",
#                                  kidney_present ~ "Kidney impairment",
#                                  immunosuppression_present ~ "Immunosuppression",
#                                  smok ~ "Smoking",
#                                  flu_vacc_present ~ "Influenza vaccine",
#                                  pneumo_vacc_present ~ "Pneumococcal vaccine",
#                                  pos_covid_test_present ~ "Positive COVID-19 test",
#                                  covid_hes_present ~ "COVID-19 hospitalisation",
#                                  covid_death_present ~ "COVID-19 death",
#                                  exacerb_present ~ "Any exacerbation in past 12 months"),
#                     percent = "column",
#                     digits = all_continuous() ~ 2,
#                     missing = "ifany",
#                     missing_text = "Missing",
#                     statistic = list(
#                       all_continuous() ~ "{mean} ({sd})",
#                       age_index ~ c("{mean} ({sd})", "{median}  \n ({p25}-{p75})"),
#                       all_categorical() ~ "{n} ({p}%)"),
#                       type = list(
#                       c(age_index) ~ "continuous2"))  %>% 
#   add_p() %>%
#   modify_header(label ~ "", all_stat_cols() ~ "**{level}**  \n N = {n}")  %>%
#   # modify_caption("Patient Characteristics") %>%
#   modify_column_alignment(columns = c(stat_1, stat_2), align = "right") %>% 
#   italicize_levels()
# 
# tab1
# #export to word
# tab1 %>%
#   as_flex_table() %>%
#   flextable::save_as_docx(path = paste0(Tables, "copd_baseline_w1_60d_ps0.7.docx"), align = "left")


####check SMDs manually
#smd_result <- tidy_smd(.df = subset_df, .vars = diabetes_present, .group = treat, .wts = c("ate_weight_stab", "ate_weight_unstab"))

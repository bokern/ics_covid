# Author: Marleen Bokern
# Date: 03/2023
# Purpose: generate propensity scores and IPTW weights for COPD cohort. Fit Cox models for all three outcomes.

packages <- c("tidyverse", "MetBrewer", "arrow", "survey", "ggplot2", "gtsummary", "survival", "survminer", "openxlsx")
installed_packages <- packages %in% rownames(installed.packages())

if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

lapply(packages, library, character.only = TRUE)

setwd(Datadir_copd)

#read in parquet dataset
df <- read_parquet("SA_copd_wave1_60d_iptw_no_triple.parquet")
subset_df <- df[!is.na(df$treatgroup),]

palette <- met.brewer("Cassatt2")

# COX AND LOGISTIC REGRESSION --------------------------------------------------------

#create object to save results
estimates <- data.frame()
subset_df$treat <- factor(subset_df$treat)
subset_df$treat <- relevel(subset_df$treat, ref = "LABA/LAMA")

# Define the survival outcome variables
outcomes <- c("covid_hes_present", "covid_death_present", "any_death_present")
time_to_outcome_vars <- c("timeinstudy2", "timeinstudy3", "timeinstudy_death_any")
weight_vars <- c("ate_weight_stab")

km_weighted <- function(weight_var, outcome, time_to_outcome_var) {
  km_surv <<- Surv(subset_df[[time_to_outcome_var]], subset_df[[outcome]]) 
  
  km_curve <- survfit(km_surv ~ subset_df$treat, data = subset_df, weights = as.numeric(subset_df[[weight_var]]))
  
  # Round the number of events to 2 decimal places
  km_curve$n.event <- round(km_curve$n.event, 2)
  
  ggsurvplot <- ggsurvplot(km_curve, data = subset_df, conf.int = T, censor = F, ylim = c(0.95, 1), xlab = "Time in days", risk.table = "absolute", risk.table.title = "Number at risk",  cumevents = TRUE, fontsize = 13, tables.height = 0.15, legend.labs = c("ICS", "LABA/LAMA"), legend.title = "", palette = c(palette[9], palette[4]), xlim = c(0, 183)) 
  ggsurvplot$plot <- ggsurvplot$plot + scale_x_continuous(breaks = c(0, 50, 100, 150, 183))
  ggsurvplot$table$theme$axis.text.y$colour <- "black"
  ggsurvplot$table$theme$axis.text.y$size <- 40
  ggsurvplot$table$theme$axis.text.x$size <- 40
  ggsurvplot$table$labels$x <- ""
  ggsurvplot$table$theme$plot.title$size <- 44
  ggsurvplot$cumevents$theme$axis.text.y$colour <- "black"
  ggsurvplot$cumevents$theme$axis.text.y$size <- 40
  ggsurvplot$cumevents$theme$axis.text.x$size <- 40
  ggsurvplot$cumevents$labels$x <- ""
  ggsurvplot$cumevents$theme$plot.title$size <- 44
  ggsurvplot$plot$theme$axis.title.x$size <- 44
  ggsurvplot$plot$theme$axis.title.y$size <- 44
  ggsurvplot$plot$theme$axis.text.x$size <- 44
  ggsurvplot$plot$theme$axis.text.y$size <- 44
  ggsurvplot$plot$theme$legend.text$size <- 44
  ggsurvplot$plot$theme$legend.key.size <- unit(4, "lines")
  
  # Define the file path
  file_path <- file.path(Graphdir, "cox_regression", paste0("km_", outcome, "_", weight_var, "_no_triple.png"))
  
  # Save the plot
  png(file_path, width = 2000, height = 1500)
  print(ggsurvplot)
  dev.off()
  
  return(ggsurvplot)
}

# Apply the function to each weight variable, outcome, and time to outcome variable
plots_list <- expand.grid(weight_vars, outcomes, time_to_outcome_vars)
plots_list$Var4 <- factor(plots_list$Var2, levels = outcomes)
plots_list$Var5 <- factor(plots_list$Var3, levels = time_to_outcome_vars)
plots_list$Var4 <- as.numeric(plots_list$Var4)
plots_list$Var5 <- as.numeric(plots_list$Var5)

# keep if the content of Var2 has the same index in the outcomes list as the content of Var3 in the time_to_outcome_vars list
plots_list <- plots_list[plots_list[,4] == plots_list[,5],]

# Loop through each row of plots_list
for (i in 1:nrow(plots_list)) {
  print(i)
  # Extract variables for the current row
  weight_var <- plots_list[i, 1]
  outcome <- plots_list[i, 2]
  time_to_outcome_var <- plots_list[i, 3]
  
  # Apply km_weighted function to the current combination of variables
  km_weighted(weight_var, outcome, time_to_outcome_var)
}

# Loop through each outcome and fit the unadjusted and IPTW models
for (j in seq_along(outcomes)){
  outcome_event <- outcomes[j]
  print(j)
  
  # Determine the position of the current outcome in the 'outcomes' vector
  outcome_position <- which(outcomes == outcome_event)
  
  # Use the 'switch' function to select the appropriate time-in-study variable
  timeinstudy_var <- switch(outcome_position,
                            as.name("timeinstudy2"),
                            as.name("timeinstudy3"),
                            as.name("timeinstudy_death_any"))
  print(timeinstudy_var)
  
  outcome_label <- switch(outcome_position,
                          as.name("COVID-19 hospitalisation"),
                          as.name("COVID-19 death"),
                          as.name("all-cause mortality"))
  
  print(outcome_label)
  
  stopifnot(outcome_label != "COVID-19 hospitalisation" || timeinstudy_var == "timeinstudy2")
  stopifnot(outcome_label != "COVID-19 death" || timeinstudy_var == "timeinstudy3")
  stopifnot(outcome_label != "all-cause mortality" || timeinstudy_var == "timeinstudy_death_any")
  
  # Formula for unadjusted model using the selected time-in-study variable
  formula_unadj_cox <- as.formula(paste("Surv(", timeinstudy_var, ", ", outcome_event, ") ~ treat"))
  
  # Fit the unadjusted model
  model_unadj_cox <- coxph(formula_unadj_cox, 
                           data = subset_df)
  
  summary_unadj_cox <- summary(model_unadj_cox)
  
  #plot schoenfeld residuals
  schoenfeld_resid <- cox.zph(model_unadj_cox)
  file_path <- file.path(Graphdir, "cox_regression", paste0("schoenfeld_resid_",outcome_event,"_unadj_no_triple.png"))
  png(file_path)
  
  plot(schoenfeld_resid, main = paste("Schoenfeld residuals for", outcome_label, "(unadjusted)"))
  dev.off()
  
  #formally test proportional hazards assumption
  res_p <- cox.zph(model_unadj_cox)$table[, "p"]
  
  # Extract coefficients and confidence intervals for unadjusted model
  summary_unadj_cox <- summary(model_unadj_cox)
  
  coef_unadj_cox <- model_unadj_cox$coefficients
  
  hr_unadj_cox <- coef_unadj_cox %>% exp()
  
  #extract the standard error of the coefficient
  se_unadj_cox <- summary_unadj_cox$coef["treatICS", "se(coef)"]
  
  #calculate the CI using the normal SE
  ci_unadj_cox <- (coef_unadj_cox + c(-1, 1) * qnorm(0.975) * se_unadj_cox) %>% exp()
  
  # Fit logistic regression model without timeinstudy_var
  model_unadj_log <- glm(get(outcome_event) ~ treat, data = subset_df, family = "binomial")
  
  # Summary of the logistic regression model
  summary(model_unadj_log)
  
  # Calculate the residuals
  residuals <- residuals(model_unadj_log)
  fitted_values <- fitted.values(model_unadj_log)
  file_path_resid <- file.path(Graphdir, "cox_regression", paste0("SA_log_resid_", outcome_event, "_unadj_no_triple_jittered.png"))
  png(file_path_resid)
  
  # Create a plot of the jittered residuals versus the jittered fitted values
  plot(
    jitter(fitted_values, factor = 0.3),
    jitter(residuals, factor = 10),
    main = paste("Log residuals for", outcome_label, "(unadjusted)"),
    cex.main = 0.9,
    xlab = "Fitted values",
    ylab = "Residuals",
    pch = 16,  # Use solid circles for points
    cex = 0.5  # Reduce point size for better visibility
  )
  
  dev.off()
  
  # Extract coefficients, odds ratios, standard errors, and confidence intervals as shown in the previous example
  coef_unadj_log <- coef(model_unadj_log)[2]
  or_unadj_log <- exp(coef_unadj_log)
  se_unadj_log <- summary(model_unadj_log)$coef[, "Std. Error"][2]
  z_value <- qnorm(0.975) # 95% confidence interval
  ci_unadj_log <- c(
  coef_unadj_log - z_value * se_unadj_log,
  coef_unadj_log + z_value * se_unadj_log)
  or_unadj_log_ci <- exp(ci_unadj_log)
  
  # Save results in the 'estimates' dataframe
  result_position <- j  # Starting position for the outcome
  
  estimates[result_position, "outcome_event"] <- gsub("_present", "", outcome_event)
  estimates[result_position, paste0("coef_unadj_cox")] <- coef_unadj_cox
  estimates[result_position, paste0("hr_unadj_cox")] <- hr_unadj_cox
  estimates[result_position, paste0("se_unadj_cox")] <- se_unadj_cox
  estimates[result_position, paste0("ci_lower_unadj_cox")] <- ci_unadj_cox[1]
  estimates[result_position, paste0("ci_upper_unadj_cox")] <- ci_unadj_cox[2]
  estimates[result_position, paste0("res_p")] <- res_p[1]
  estimates[result_position, paste0("coef_unadj_log")] <- coef_unadj_log
  estimates[result_position, paste0("or_unadj_log")] <- or_unadj_log
  estimates[result_position, paste0("se_unadj_log")] <- se_unadj_log
  estimates[result_position, paste0("ci_lower_unadj_log")] <- or_unadj_log_ci[1]
  estimates[result_position, paste0("ci_upper_unadj_log")] <- or_unadj_log_ci[2]
  
  #adjusted models 
  for (w in seq_along(weight_vars)) {
    print(w)
    
    weight_name <- weight_vars[w]
    print(weight_name)
    
    # Determine the position of the current outcome in the 'outcomes' vector
    weight_position <- which(weight_vars == weight_name)
    
    # Use the 'switch' function to select the appropriate time-in-study variable
    weight_var <- switch(weight_position,
                         as.name("ate_weight_stab"), 
                         as.name("att_weight_stab"), 
                         as.name("ate_weight_unstab"), 
                         as.name("att_weight_unstab"))
    print(weight_var)
    
    weight_label <- switch(weight_position,
                           as.name("stabilised ATE weights"),
                           as.name("stabilised ATT weights"), 
                           as.name("unstabilised ATE weights"), 
                           as.name("unstabilised ATT weights"))
    
    print(weight_label)
    
    # Formula for IPTW model using the selected time-in-study variable
    formula_iptw_cox <- as.formula(paste("Surv(", timeinstudy_var, ", ", outcome_event, ") ~ treat"))
    # Fit the IPTW model
    model_iptw_cox <- coxph(formula_iptw_cox, data = subset_df, weights = (subset_df[[weight_name]]))
    
    # Plot schoenfeld residuals for IPTW model
    schoenfeld_resid_iptw <- cox.zph(model_iptw_cox, transform = 'identity')
    file_path_iptw <- file.path(Graphdir, "cox_regression", paste0("schoenfeld_resid_", outcome_event, "_IPTW_", weight_name, "_no_triple.png"))
    png(file_path_iptw)
    plot(schoenfeld_resid_iptw, main = paste("Schoenfeld residuals for", outcome_label, "using", weight_label),
         cex.main = 0.9)
    dev.off()
    res_p_iptw <- cox.zph(model_iptw_cox)$table[, "p"]
    
    # Extract coefficients and confidence intervals for IPTW model
    coef_iptw_cox <- coef(model_iptw_cox)
    hr_iptw_cox <- exp(coef_iptw_cox)
    se_iptw_cox <- summary(model_iptw_cox)$coef["treatICS", "se(coef)"]
    se_iptw_robust_cox <- summary(model_iptw_cox)$coef["treatICS", "robust se"]
    ci_iptw_cox <- exp(coef_iptw_cox + c(-1, 1) * qnorm(0.975) * se_iptw_robust_cox)
    
    # Fit logistic regression model without timeinstudy_var
    model_iptw_log <- glm(get(outcome_event) ~ treat, data = subset_df, family = "quasibinomial", weights = (subset_df[[weight_name]]))
    
    # Summary of the logistic regression model
    summary(model_iptw_log)
    
    # Extract coefficients, odds ratios, standard errors, and confidence intervals as shown in the previous example
    coef_iptw_log <- coef(model_iptw_log)[2]
    or_iptw_log <- exp(coef_iptw_log)
    se_iptw_log <- summary(model_iptw_log)$coef[, "Std. Error"][2]
    z_value <- qnorm(0.975)
    ci_iptw_log <- c(
      coef_iptw_log - z_value * se_iptw_log,
      coef_iptw_log + z_value * se_iptw_log)
    or_iptw_log_ci <- exp(ci_iptw_log)
    
    # Calculate the residuals
    residuals <- residuals(model_iptw_log, type = "response")
    fitted_values <- predict(model_iptw_log, type = "response")
    file_path_resid <- file.path(Graphdir, "cox_regression", paste0("SA_log_resid_", outcome_event, "_IPTW_", weight_name, "_no_triple_jittered.png"))
    
    png(file_path_resid)
    # Create a plot of the residuals versus the fitted values
    plot(
      jitter(fitted_values, factor = 0.3),
      jitter(residuals, factor = 10),
      main = paste("Log residuals for", outcome_label, "using", weight_label),
      cex.main = 0.9,
      xlab = "Fitted values",
      ylab = "Residuals",
      pch = 16,  # Use solid circles for points
      cex = 0.5  # Reduce point size for better visibility
    )
    
    dev.off()
    
    
    # Save estimates to dataframe
    estimates[result_position, paste0("coef_iptw_cox_", weight_name)] <- coef_iptw_cox
    estimates[result_position, paste0("hr_iptw_cox_", weight_name)] <- hr_iptw_cox
    estimates[result_position, paste0("se_normal_iptw_cox_", weight_name)] <- se_iptw_cox
    estimates[result_position, paste0("se_robust_iptw_cox_", weight_name)] <- se_iptw_robust_cox
    estimates[result_position, paste0("ci_lower_iptw_cox_", weight_name)] <- ci_iptw_cox[1]
    estimates[result_position, paste0("ci_upper_iptw_cox_", weight_name)] <- ci_iptw_cox[2]
    estimates[result_position, paste0("res_p_iptw_cox_", weight_name)] <- res_p_iptw[1]
    estimates[result_position, paste0("coef_iptw_log_", weight_name)] <- coef_iptw_log
    estimates[result_position, paste0("or_iptw_log_", weight_name)] <- or_iptw_log
    estimates[result_position, paste0("se_iptw_normal_log_", weight_name)] <- se_iptw_log
    estimates[result_position, paste0("ci_lower_iptw_log_", weight_name)] <- or_iptw_log_ci[1]
    estimates[result_position, paste0("ci_upper_iptw_log_", weight_name)] <- or_iptw_log_ci[2]
    
    rm(list = c("coef_iptw_cox", "hr_iptw_cox", "se_iptw_cox", "se_iptw_robust_cox", "ci_iptw_cox"))
  }
  
  rm(list = c("coef_unadj_cox", "hr_unadj_cox", "se_unadj_cox", "ci_unadj_cox", "model_unadj_cox"))
}

#estimates <- t(estimates)
estimates <- as.data.frame(estimates)

# Write the modified data frame to Parquet
arrow::write_parquet(estimates, file.path(Tables, "QBA", "SA_cox_log_regression_estimates_no_triple.parquet"))


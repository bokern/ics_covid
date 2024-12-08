############################
# Project: Outcome Misclassification adjustment -- RECORD LEVEL
# Programmer name: Richard MacLehose, adapted by Marleen Bokern
# Date Started:  1/16/23
#############################  
# Updated versions can be found at: https://sites.google.com/site/biasanalysis/short-code 
# code to run record-level PBA for deaths in COPD dataset. Runs main analysis, analysis excluding triple therapy users, analysis using 6-month-post-prescription exposure definition, and analysis including people with asthma in a function.
#############################

start_time_total <- Sys.time()
packages <- c("tidyverse", "sandwich", "ggplot2", "MASS", "MetBrewer", "survival", "survminer", "survey", "purrr", "future", "furrr", "arrow", "survival", "survminer", "openxlsx", "data.table")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE, verbose = FALSE))
}))

palette <- met.brewer("Cassatt2")

#NUMBER OF PBA SIMULATIONS:
samp = 10000
set.seed(123)

setwd(Datadir_copd)

death_record_qba <- function(inputfile, samples, output_ext) {
  
  # Load analytical file and select and format relevant variables
  D <- read_parquet(inputfile)
  D <- D %>% dplyr::select(patid, treat, treatgroup, covid_death_present, any_death_present, timeinstudy3, ate_weight_stab)

  D$covid_death_present <- as.numeric(D$covid_death_present)

  D <- D[!is.na(D$treatgroup),]

  D$treatgroup_num <- ifelse(D$treat == "ICS", 1,
                             ifelse(D$treat == "LABA/LAMA", 0, NA))

  D$treat <- factor(D$treat)
  D$treat <- relevel(D$treat, ref = "LABA/LAMA")

  df <- D %>% filter(any_death_present == 1)

  total_exp <- sum(D$treat == "ICS")
  total_unexp <- sum(D$treat == "LABA/LAMA")

  death_exp <- sum(df$treat == "ICS")
  death_unexp <- sum(df$treat == "LABA/LAMA")

  unique_patid <- D[D$any_death_present == 1, ]
  unique_patid <- unique_patid %>% dplyr::select(patid, covid_death_present, treat) %>% unique()

  #count unique patids for each combination of covid_hes_present and treat
  a <- sum(unique_patid$covid_death_present == 1 & unique_patid$treat == "ICS")
  b <- sum(unique_patid$covid_death_present == 1 & unique_patid$treat == "LABA/LAMA")
  c <- sum(unique_patid$covid_death_present == 0 & unique_patid$treat == "ICS")
  d <- sum(unique_patid$covid_death_present == 0 & unique_patid$treat == "LABA/LAMA")

  #number of obs
  n=dim(D)[1]

  #include id in dataset
  D$id=1:n

  survey_design_iptw <- svydesign(ids = ~1, weights = D$ate_weight_stab, data = D)

  # Create the design object with the survey design and sampling weights
  survey_design_unadj <- update(svydesign(ids = ~1, data = D), d = D$covid_death_present, treat = D$treat)

  # Fit the logistic regression model using svyglm
  model_unadj_LR <- svyglm(covid_death_present ~ treat, design = survey_design_unadj, family = quasibinomial())
  model_ate_LR <- svyglm(covid_death_present ~ treat, design = survey_design_iptw, family = quasibinomial())

  or_ate_conventional <- exp(model_ate_LR$coefficients[2])
  or_unadj_conventional <- exp(model_unadj_LR$coefficients[2])

  # Fit a Cox proportional hazards model without IPTW weighting
  cox_model_unadj <- coxph(Surv(timeinstudy3, covid_death_present) ~ treat, data = D)
  #summary(cox_model_unadj)

  # Fit a Cox proportional hazards model with IPTW weighting
  cox_model_iptw <- coxph(Surv(timeinstudy3, covid_death_present) ~ treat, weights = D$ate_weight_stab, data = D)
  #summary(cox_model_iptw)

  hr_ate_conventional <- exp(cox_model_iptw$coefficients)
  hr_unadj_conventional <- exp(cox_model_unadj$coefficients)

# # Sampling ----------------------------------------------------------------

  #Bias parameters
  se1.a <- 39.6
  se1.b <- 15.4
  se0.a <- se1.a
  se0.b <- se1.b

  #bias parameters for sensitivity (beta distribution)
  sp1.a <- 91.2
  sp1.b <- 3.86
  sp0.a <- sp1.a
  sp0.b <- sp1.b

  # Sample Sens/Spec for non-differential misclassification
  se1 = rbeta(samp, se1.a, se1.b)
  se0 = se1

  sp1 = rbeta(samp, sp1.a, sp1.b)
  sp0 = sp1

  # NOTE: if you have differential misclassification, you would specify a correlation rho1 between se1 and sp1, and then sample from a bivariate normal distribution
  #sample from the quantile funciton of the beta distribution

  # Correlation
  # rho1 = .80
  # V = matrix(c(1, rho1, rho1, 1), ncol = 2)
  # Z = mvrnorm(I, c(0, 0), V)
  # U = pnorm(Z)
  # se1 = qbeta(U[, 1], se1.a, se1.b)
  # sp1 = qbeta(U[, 2], sp1.a, sp1.b)

  # #aggregate sampled bias parameters
  imp.d=data.frame(iter=1:samp,se_E1=se1,se_E0=se0,sp_E1=sp1,sp_E0=sp0)

  #Create A0,B0,C0,D0
  imp.d$A0=(a-(a+c)*(1-imp.d$sp_E1))/(imp.d$se_E1+imp.d$sp_E1-1)
  imp.d$C0=(a+c)-imp.d$A0
  imp.d$B0=(b-(b+d)*(1-imp.d$sp_E0))/(imp.d$se_E0+imp.d$sp_E0-1)
  imp.d$D0=(b+d)-imp.d$B0

  #Sample pr(D1|e=1) pr(D1|e=0)
  flag1 = (imp.d$A0<=0)|(imp.d$B0<=0)|(imp.d$C0<=0)|(imp.d$D0<=0)
  imp.d = imp.d[flag1=="FALSE",]
  impossible <- sum(flag1)
  I = samp-sum(flag1)
  imp.d$pde1 = rbeta(I,imp.d$A0,imp.d$C0)
  imp.d$pde0 = rbeta(I,imp.d$B0,imp.d$D0)

  #Compute predictive values for imputation
  imp.d$PPV_e1 = imp.d$se_E1*imp.d$pde1/(imp.d$se_E1*imp.d$pde1+(1-imp.d$pde1)*(1-imp.d$sp_E1))
  imp.d$PPV_e0 = imp.d$se_E0*imp.d$pde0/(imp.d$se_E0*imp.d$pde0+(1-imp.d$pde0)*(1-imp.d$sp_E0))
  imp.d$NPV_e1 = (imp.d$sp_E1*(1-imp.d$pde1))/((imp.d$sp_E1*(1-imp.d$pde1))+(1-imp.d$se_E1)*(imp.d$pde1))
  imp.d$NPV_e0 = (imp.d$sp_E0*(1-imp.d$pde0))/((imp.d$sp_E0*(1-imp.d$pde0))+(1-imp.d$se_E0)*(imp.d$pde0))

  DT <- setDT(D)
  #save non-dead patients as separate object
  DT_no_death <- DT[any_death_present == 0]

  DT_no_death[, p := 0]
  DT_no_death[, d := 0]

  #remove those patients from DT
  DT <- DT[any_death_present == 1]

  #setup output file to store results
  output <- data.frame()

  start_time <- Sys.time()

  #Loop through I iterations and impute new exposures at each step
  # pba_deaths_record <- function(iter, se_E1, se_E0, sp_E1, sp_E0, A0, C0, B0, D0,
  #                              pde1, pde0, PPV_e1, PPV_e0, NPV_e1, NPV_e0, DT_main) {
  # 
  #   #DT <- copy(DT_main) #copy of DT to stay safe, not sure if needed
  #   stopifnot(sum(DT$any_death_present) == nrow(DT))
  # 
  #   DT[, p := treatgroup_num * covid_death_present * PPV_e1 + #probability of outcome in treated group
  #        treatgroup_num * (1 - covid_death_present) * (1 - NPV_e1) + #probability of no outcome in treated group
  #        (1 - treatgroup_num) * covid_death_present * PPV_e0 + #probability of outcome in control group
  #        (1 - treatgroup_num) * (1 - covid_death_present) * (1 - NPV_e0)] #probability of no outcome in control group
  # 
  #   DT[, d := rbinom(.N, 1, p)]
  # 
  #   missing_cols <- setdiff(names(DT), names(DT_no_death))
  # 
  #   # Fill missing columns in DT_no_hosp with 0
  #   for (col in missing_cols) {
  #     DT_no_hosp[[col]] <- 0
  #   }
  # 
  #   # Append the rows from DT_no_hosp to DT
  #   DT <- rbindlist(list(DT, DT_no_death), use.names = TRUE)
  # 
  #   # Check if the number of rows in DT is less than or equal to the number of unique patid
  #   stopifnot(nrow(DT) == uniqueN(DT$patid))
  # 
  #   #count where d=1 and treat=ICS
  #   exp1d1 <- sum(DT$d ==  1 & DT$treat == "ICS", na.rm = TRUE)
  #   exp1d0 <- sum(DT$d ==  0 & DT$treat == "ICS", na.rm = TRUE)
  #   exp0d1 <- sum(DT$d ==  1 & DT$treat == "LABA/LAMA", na.rm = TRUE)
  #   exp0d0 <- sum(DT$d ==  0 & DT$treat == "LABA/LAMA", na.rm = TRUE)
  # 
  #   D_iptw <- DT[!is.na(ate_weight_stab)]
  # 
  #   invisible({
  # 
  #     # Create a new survey design for each iteration
  #     design_unadj_i <- svydesign(ids = ~1, data = DT, weights = ~1);
  #     design_ate_i <- svydesign(ids = ~1, data = D_iptw, weights = ~ate_weight_stab);
  # 
  #     z <- rnorm(1);
  #     # Unadjusted logistic model
  #     model_unadj_LR <- svyglm(d ~ treat, design = design_unadj_i, family = binomial());
  #     coef_unadj_LR <- summary(model_unadj_LR)$coeff[2, 1];
  #     se_unadj_LR <- SE(model_unadj_LR)[2];
  # 
  #     # Weighted logistic model
  #     model_ate_LR <- svyglm(d ~ treat, design = design_ate_i, family = quasibinomial());
  #     coef_ate_LR <- summary(model_ate_LR)$coeff[2, 1];
  #     se_ate_LR <- SE(model_ate_LR)[2];
  # 
  #     # Unadjusted Cox model
  #     model_unadj_cox <- svycoxph(Surv(timeinstudy3, d) ~ treat, design = design_unadj_i);
  #     coef_unadj_cox <- summary(model_unadj_cox)$coefficients[1, "coef"];
  #     se_unadj_cox <- sqrt(vcov(model_unadj_cox)["treatICS", "treatICS"]);
  # 
  #     # Weighted Cox model
  #     model_ate_cox <- svycoxph(Surv(timeinstudy3, d) ~ treat, design = design_ate_i);
  #     coef_ate_cox <- summary(model_ate_cox)$coeff[1];
  #     se_ate_cox <- sqrt(diag(vcov(model_ate_cox))["treatICS"]);
  #   })
  # 
  #   output <- tibble(
  #     a1 = exp1d1,
  #     b1 = exp0d1,
  #     c1 = exp1d0,
  #     d1 = exp0d0,
  #     or_unadj = exp(coef_unadj_LR),
  #     or_tot = exp(coef_unadj_LR + z * se_unadj_LR),
  #     ci_lower_LR = exp(coef_unadj_LR - z * se_unadj_LR),
  #     ci_upper_LR = exp(coef_unadj_LR + z * se_unadj_LR),
  #     or_adj_ate = exp(coef_ate_LR),
  #     or_tot_ate = exp(coef_ate_LR + z * se_ate_LR),
  #     ci_lower_ate_LR = exp(coef_ate_LR - z * se_ate_LR),
  #     ci_upper_ate_LR = exp(coef_ate_LR + z * se_ate_LR),
  #     hr_unadj = exp(coef_unadj_cox),
  #     hr_tot = exp(coef_unadj_cox + z * se_unadj_cox),
  #     ci_lower_cox = exp(coef_unadj_cox + z * se_unadj_cox),
  #     ci_upper_cox = exp(coef_unadj_cox - z * se_unadj_cox),
  #     hr_adj_ate = exp(coef_ate_cox),
  #     hr_tot_ate = exp(coef_ate_cox + z * se_ate_cox),
  #     ci_lower_ate_cox = exp(coef_ate_cox - z * se_ate_cox),
  #     ci_upper_ate_cox = exp(coef_ate_cox + z * se_ate_cox),
  #     se = se_E1,
  #     sp = sp_E1,
  #     PPV_e1 = PPV_e1,
  #     NPV_e1 = NPV_e1,
  #     PPV_e0 = PPV_e0,
  #     NPV_e0 = NPV_e0,
  #     PrevD_exp = pde1,
  #     PrevD_unexp = pde0
  #   )
  # 
  #   return(output)
  # }

  # # Use future_map to apply this function in parallel
  # plan(multisession)
  # 
  # # Combine the results into a single data frame
  # output_df <- pmap_dfr(.l=imp.d[1:I,], .f=pba_deaths_record, DT_main = DT, .progress = TRUE)
  # output_df <- as.data.frame(output_df)
  # 
  # write_parquet(output_df, paste0("qba_death_record_full_sample", output_ext, ".parquet"))
   output_df <- read_parquet(paste0("qba_death_record_full_sample", output_ext, ".parquet"))
  # end_time <- Sys.time()
  # run_time <- end_time - start_time
  # print(run_time)
  # #extract the unit from runtime
  # run_time_unit <- attr(run_time, "units")

  result_vars <- c("or_tot", "hr_tot", "or_tot_ate", "hr_tot_ate")

  # Initialize an empty list to store results
  res <- list()

  # Loop through result_vars
  for (var in result_vars) {
    # Calculate percentiles
    percentiles <- quantile(output_df[[var]], probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

    # Calculate I-length(output_df$or_unadj) and I
    impossible_values <- samp-length(output_df[[var]])
    simulations <- length(output_df[[var]])

    # Combine percentiles, impossible_values, and simulations into a single vector
    res[[var]] <- c(percentiles, impossible_values, simulations)
  }

  # Create a data frame qba_df from the results stored in res
  qba_df <- do.call(cbind, res)
  qba_df <- as.data.frame(qba_df)
  new_row_names <- c("p2.5", "p50", "p97.5", "impossible_values", "sims")
  row.names(qba_df) <- new_row_names

  sims_values <- qba_df["sims", ]

  # Check if all 'sims' values are equal across columns
  stopifnot(all(sims_values == sims_values[, 1]))

  # Add the suffix to each variable name
  qba_df <- rename_with(qba_df, ~paste0(., "_death_record"))
  write_parquet(qba_df, file.path(Tables, "QBA", paste0("qba_death_record_results", output_ext, ".parquet")))

  # add stats to excel file to save the number of iterations, date, time, outcome, summary level vs record level
  # import runtime.xlsx
  # runtime_df <- read.xlsx(file.path(Tables, "QBA", "runtime.xlsx"))
  # 
  # time_stats <- data.frame(
  #   outcome = "death",
  #   analysis = "record level",
  #   sa = output_ext,
  #   niter = format(I, scientific = FALSE),
  #   date = format(Sys.Date(), "%Y-%m-%d"),
  #   time = format(Sys.time(), "%H:%M:%S"),  # Format Sys.time() to display only hour, minute, second
  #   run_time = run_time,
  #   run_time_unit = run_time_unit
  # )
  # 
  # #add time stats to runtime.df
  # runtime_df <- rbind(runtime_df, time_stats)
  # #add this line to runtime.xlsx
  # write.xlsx(runtime_df, file.path(Tables, "QBA", "runtime.xlsx"))
  
  # Plotting results --------------------------------------------------------
  sims <- qba_df["sims", 1]
  median_value_or <- median(output_df$or_tot_ate, na.rm = TRUE)
  pal <- c("IPT-weighted OR (no QBA)" = palette[4], "Median OR after QBA" = palette[10])
  pal2 <- c("IPT-weighted OR (no QBA)" = "dashed", "Median OR after QBA" = "dotted")
  plot <- ggplot(data = output_df, aes(x = or_tot_ate)) +
    geom_histogram(binwidth = 0.002, fill = palette[6]) +
    geom_vline(aes(xintercept = or_ate_conventional, color = "IPT-weighted OR (no QBA)", linetype = "IPT-weighted OR (no QBA)"), size = 0.9) +
    geom_vline(aes(xintercept = median_value_or, color = "Median OR after QBA", linetype = "Median OR after QBA"), size = 0.9) +
    xlab("Odds Ratio") +
    ylab("Frequency") +
    ggtitle(paste0("Odds ratios adjusted for outcome misclassification, IPT-weighted (n = ", sims, ")")) +
    xlim(0, 3) +
    theme_classic() +
    scale_color_manual(values = pal) +
    scale_linetype_manual(values = pal2) +
    labs(color = "", linetype = "") +
    theme(plot.background = element_rect(fill = "white"),
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12), 
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.text = element_text(size = 11),
          legend.key.size = unit(1, "cm"))
  file_path <- file.path(Graphdir, "QBA", "copd_death_w1",  paste0("adjusted_OR_ipt_RL_pba_death", output_ext, ".png"))
  ggsave(file_path, plot, width = 8, height = 4)
  
  median_value_or <- median(output_df$or_tot, na.rm = TRUE)
  pal <- c("unadjusted OR (no QBA)" = palette[4], "Median OR after QBA" = palette[10])
  pal2 <- c("unadjusted OR (no QBA)" = "dashed", "Median OR after QBA" = "dotted")
  plot <- ggplot(data = output_df, aes(x = or_tot)) +
    geom_histogram(binwidth = 0.002, fill = palette[6]) +
    geom_vline(aes(xintercept = or_unadj_conventional, color = "unadjusted OR (no QBA)", linetype = "unadjusted OR (no QBA)"), size = 0.9) +
    geom_vline(aes(xintercept = median_value_or, color = "Median OR after QBA", linetype = "Median OR after QBA"), size = 0.9) +
    xlab("Odds Ratio") +
    ylab("Frequency") +
    ggtitle(paste0("Odds ratios adjusted for outcome misclassification, unadjusted (n = ", sims, ")")) +
    xlim(0, 3) +
    theme_classic() +
    scale_color_manual(values = pal) +
    scale_linetype_manual(values = pal2) +
    labs(color = "", linetype = "") +
    theme(plot.background = element_rect(fill = "white"),
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12), 
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.text = element_text(size = 11),
          legend.key.size = unit(1, "cm"))
  file_path <- file.path(Graphdir, "QBA", "copd_death_w1", paste0("unadj_OR_RL_pba_death", output_ext, ".png"))
  ggsave(file_path, plot, width = 8, height = 4)
  
  
  median_value_hr <- median(output_df$hr_tot_ate, na.rm = TRUE)
  pal <- c("IPT-weighted HR (no QBA)" = palette[4], "Median HR after QBA" = palette[10])
  pal2 <- c("IPT-weighted HR (no QBA)" = "dashed", "Median HR after QBA" = "dotted")
  plot <- ggplot(data = output_df, aes(x = hr_tot_ate)) +
    geom_histogram(binwidth = 0.002, fill = palette[6]) +
    geom_vline(aes(xintercept = hr_ate_conventional, color = "IPT-weighted HR (no QBA)", linetype = "IPT-weighted HR (no QBA)"), size = 0.9) +
    geom_vline(aes(xintercept = median_value_hr, color = "Median HR after QBA", linetype = "Median HR after QBA"), size = 0.9) +
    xlab("Hazard Ratio") +
    ylab("Frequency") +
    ggtitle(paste0("Hazard ratios adjusted for outcome misclassification, IPT-weighted (n = ", sims, ")")) +
    xlim(0, 3) +
    theme_classic() +
    scale_color_manual(values = pal) +
    scale_linetype_manual(values = pal2) +
    labs(color = "", linetype = "") +
    theme(plot.background = element_rect(fill = "white"),
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12), 
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.text = element_text(size = 11),
          legend.key.size = unit(1, "cm"))
  file_path <- file.path(Graphdir, "QBA", "copd_death_w1",  paste0("adjusted_HR_ipt_RL_pba_death", output_ext, ".png"))
  ggsave(file_path, plot, width = 8, height = 4)
  
  median_value_hr <- median(output_df$hr_tot, na.rm = TRUE)
  pal <- c("unadjusted HR (no QBA)" = palette[4], "Median HR after QBA" = palette[10])
  pal2 <- c("unadjusted HR (no QBA)" = "dashed", "Median HR after QBA" = "dotted")
  plot <- ggplot(data = output_df, aes(x = hr_tot)) +
    geom_histogram(binwidth = 0.002, fill = palette[6]) +
    geom_vline(aes(xintercept = hr_unadj_conventional, col = "unadjusted HR (no QBA)", linetype = "unadjusted HR (no QBA)"), size = 0.9) +
    geom_vline(aes(xintercept = median_value_hr, col = "Median HR after QBA", linetype = "Median HR after QBA"), size = 0.9) +
    xlab("Hazard Ratio") +
    ylab("Frequency") +
    ggtitle(paste0("Hazard ratios adjusted for outcome misclassification, unadjusted (n = ", sims, ")")) +
    xlim(0, 3) +
    theme_classic() +
    scale_color_manual(values = pal) +
    scale_linetype_manual(values = pal2) +
    labs(color = "", linetype = "") +
    theme(plot.background = element_rect(fill = "white"),
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12), 
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.text = element_text(size = 11),
          legend.key.size = unit(1, "cm"))
  file_path <- file.path(Graphdir, "QBA", "copd_death_w1", paste0("unadj_hr_RL_pba_death", output_ext, ".png"))
  ggsave(file_path, plot, width = 8, height = 4)
  
    plot <- ggplot(data = output_df, aes(x = se)) +
    geom_histogram(binwidth = 0.002, fill = palette[2]) +
    xlab("Sampled sensitivity") +
    ylab("Frequency") +
    ggtitle(paste0("Sampled sensitivity (alpha = ", se1.a, ", beta = ", se1.b, ", n = ", samp, ")")) +
    xlim(0, 1) +
    theme_classic()+
    theme(plot.background = element_rect(fill = "white"),
            axis.text.x = element_text(size = 12), 
            axis.text.y = element_text(size = 12), 
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14)) +
    scale_color_manual(values = palette[6], labels = c("Sensitivity")) +
    labs(color = "")
  file_path <- file.path(Graphdir, "QBA", "copd_death_w1", paste0("sampled_se_RL_pba_death", output_ext, ".png"))
  ggsave(file_path, plot, width = 8, height = 4)
  
  plot <- ggplot(data = output_df, aes(x = sp)) +
    geom_histogram(binwidth = 0.002, fill = palette[2]) +
    xlab("Sampled specificity") +
    ylab("Frequency") +
    ggtitle(paste0("Sampled specificity (alpha = ", sp1.a, ", beta = ", sp1.b, ", n = ", samp, ")")) +
    xlim(0, 1) +
    theme_classic() +
    theme(plot.background = element_rect(fill = "white"),
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12), 
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    scale_color_manual(values = palette[6], labels = c("Specificity")) +
    labs(color = "")
  file_path <- file.path(Graphdir, "QBA", "copd_death_w1",  paste0("sampled_sp_RL_pba_death", output_ext, ".png"))
  ggsave(file_path, plot, width = 8, height = 4)
  
  plot_ppv <- ggplot(data = output_df, aes(x = PPV_e1, color = "ICS")) +
    geom_density(adjust = 1, aes(color = "ICS"), fill = NA, size = 1) +
    geom_density(aes(x = PPV_e0, color = "LABA/LAMA"), adjust = 1, fill = NA, size = 1) +
    xlab("Sampled PPV") +
    ylab("Density") +
    ggtitle("Sampled PPV by treatment group") +
    xlim(0, 1) +
    theme_classic() +
    theme(plot.background = element_rect(fill = "white"),
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12), 
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.text = element_text(size = 11)) +
    scale_color_manual(values = c("ICS" = palette[9], "LABA/LAMA" = palette[4]), labels = c("ICS", "LABA/LAMA")) +
    guides(color = guide_legend(title = NULL, override.aes = list(fill = c(palette[9], palette[4]))))
  file_path <- file.path(Graphdir, "QBA", "copd_death_w1", paste0("sampled_PPV_RL_pba_death", output_ext, ".png"))
  ggsave(file_path, plot_ppv, width = 8, height = 4)
  
  plot_npv <- ggplot(data = output_df, aes(x = NPV_e1, color = "ICS")) +
    geom_density(adjust = 1, aes(color = "ICS"), fill = NA, size = 1) +
    geom_density(aes(x = NPV_e0, color = "LABA/LAMA"), adjust = 1, fill = NA, size = 1) +
    xlab("Sampled NPV") +
    ylab("Density") +
    ggtitle("Sampled NPV by treatment group") +
    xlim(0, 1) +
    theme_classic() +
    theme(plot.background = element_rect(fill = "white"),
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12), 
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.text = element_text(size = 11)) +
    scale_color_manual(values = c("ICS" = palette[9], "LABA/LAMA" = palette[4]), labels = c("ICS", "LABA/LAMA")) +
    guides(color = guide_legend(title = NULL, override.aes = list(fill = c(palette[9], palette[4]))))
  file_path <- file.path(Graphdir, "QBA", "copd_death_w1", paste0("sampled_NPV_RL_pba_death", output_ext, ".png"))
  ggsave(file_path, plot_npv, width = 8, height = 4)
  
  # Extract the legend from one of the plots
  legend_plot <- plot_ppv +
    scale_color_manual(values = c("ICS" = palette[9], "LABA/LAMA" = palette[4])) +
    scale_fill_manual(values = c("ICS" = palette[9], "LABA/LAMA" = palette[4])) +
    theme(
      legend.text = element_text(size = 18),
      legend.box = "horizontal",
      legend.spacing.x = unit(0.6, "cm"),
      legend.key.size = unit(0.6, "cm"),
      plot.margin = margin(t = -20, r = 0, b = 0, l = 0)
    ) +
    guides(
      color = guide_legend(
        override.aes = list(
          fill = c(palette[9], palette[4]),
          size = 5,
          shape = 15
        ), 
        title = NULL,
        nrow = 1,
        keywidth = unit(0.8, "cm"),
        keyheight = unit(0.8, "cm")
      )
    )
  
  shared_legend <- get_legend(legend_plot)
  
  # Combine the PPV and NPV plots without legends
  plot_ppv <- plot_ppv + theme(legend.position = "none")
  plot_npv <- plot_npv + theme(legend.position = "none") 
  #remove y axis label from the NPV plot
  plot_npv <- plot_npv + theme(axis.title.y = element_blank())
  combined_plot <- cowplot::plot_grid(
    plot_ppv, plot_npv, ncol = 2, labels = c("A", "B")
  )
  
  # Add the legend below the combined plot
  final_plot <- cowplot::plot_grid(
    combined_plot, shared_legend, ncol = 1, rel_heights = c(0.9, 0.1)
  )
  
  # Save the final plot
  file_path <- file.path(Graphdir, "QBA", "copd_death_w1", paste0("sampled_PPV_NPV_RL_pba_death", output_ext, ".png"))
  ggsave(file_path, final_plot, width = 10, height = 4)
  
  plot <- ggplot(data = output_df, aes(x = PrevD_exp, color = "ICS")) +
    geom_density(adjust = 1, fill = NA, size = 1) +
    geom_density(aes(x = PrevD_unexp, color = "LABA/LAMA"), adjust = 1, fill = NA, size = 1) +
    xlab("Sampled Outcome Prevalence") +
    ylab("Density") +
    ggtitle("Sampled Outcome Prevalence by treatment group") +
    xlim(0, 1) +
    theme_classic() +
    theme(plot.background = element_rect(fill = "white"),
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12), 
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.text = element_text(size = 11)) +
    scale_color_manual(values = c("ICS" = palette[9], "LABA/LAMA" = palette[4]), labels = c("ICS", "LABA/LAMA")) +
    guides(color = guide_legend(title = NULL, override.aes = list(fill = c(palette[9], palette[4]))))
  file_path <- file.path(Graphdir, "QBA", "copd_death_w1", paste0("sampled_Prev_RL_pba_death", output_ext, ".png"))
  ggsave(file_path, plot, width = 8, height = 4)

}

#list inputfiles
inputfile <- c("copd_wave1_60d_iptw.parquet", "SA_copd_wave1_60d_iptw_no_triple.parquet", "SA_copd_wave1_6m_iptw.parquet", "SA_copd_wave1_60d_iptw_all.parquet")
samples <- c("qba_death_record_full_sample.parquet", "SA_qba_death_record_full_sample_no_triple.parquet", "SA_qba_death_record_full_sample_6m.parquet", "SA_qba_death_record_full_sample_all.parquet")
output_ext <- c("", "_no_triple",  "_6m", "_all")

#run funciton for 6m and all
# inputfile <- c("SA_copd_wave1_6m_iptw.parquet", "SA_copd_wave1_60d_iptw_all.parquet")
# samples <- c("SA_qba_death_record_full_sample_6m.parquet", "SA_qba_death_record_full_sample_all.parquet")
# output_ext <- c("_6m", "_all")


#run the function
mapply(death_record_qba, inputfile, samples, output_ext)

end_time_total <- Sys.time()
print(end_time_total - start_time_total)
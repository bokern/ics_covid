# Author: Marleen Bokern
# Date: 10/2023
# Purpose: application of summary level PBA to deaths. Code adapted from : https://academic.oup.com/ije/article/52/5/1624/7152433

###SET UP
packages <- c("tidyverse", "episensr", "arrow", "ggplot2", "dplyr", "sandwich", "MASS", "MetBrewer", "openxlsx", "arrow")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE, verbose = FALSE))
}))

palette <- met.brewer("Cassatt2")

setwd(Datadir_copd)
set.seed(123)

death_summary_qba <- function(inputfile, conventional_result, samples, output_ext) {
  
  D <- read_parquet(inputfile)
  
  D$treat <- factor(D$treat)
  D$treat <- relevel(D$treat, ref = "LABA/LAMA")
  
    D <- D %>%
    group_by(treat) %>%
    summarise(
      any_death_present = sum(any_death_present > 0, na.rm = TRUE),
      covid_death_present = sum(covid_death_present > 0, na.rm = TRUE),
      n = n_distinct(patid)
    ) %>%
    ungroup()
  
  I <- 100000
  
  estimates <- read_parquet(file.path(Tables, "QBA", conventional_result))
# Bias analysis -----------------------------------------------------------

#the QBA steps are conducted only in those who died of any cause (any_death_present == 1)

#Inputs for this macro are as follow:
#Data
#a=number of exposed cases
#b=number of unexposed cases
#c=number of exposed controls
#d=number of unexposed controls

# from analysis without QBA, ATE (stab) weighted Cox model:

total_exp <- sum(D$treat == "ICS")
total_unexp <- sum(D$treat == "LABA/LAMA")

death_exp <- sum(D$treat == "ICS" & D$any_death_present == 1)
death_unexp <- sum(D$treat == "LABA/LAMA" & D$any_death_present == 1)

a1 <- sum(D$treat == "ICS" & D$covid_death_present == 1)
b1 <- sum(D$treat == "LABA/LAMA" & D$covid_death_present == 1)
c1 <- death_exp - a1 #non-covid deaths among exposed
d1 <- death_unexp - b1 #non-covid deaths among unexposed

#bias parameters for sensitivity (beta distribution)
#se1.a=alpha parameter for sens among cases
#se1.b=beta parameter for sens among cases
#se0.a=alpha parameter for sens among controls
#se0.b=beta parameter for sens among controls

#this is the set up for non-differential misclassification. For differential misclassification, se and sp for exposed and unexposed should be different. Correlation should be specified for differential misclassification
se1.a <- 39.6
se1.b <- 15.4
se0.a <- se1.a
se0.b <- se1.b

#bias parameters for sensitivity (beta distribution)
#se1.a=alpha parameter for sens among cases
#se1.b=beta parameter for sens among cases
#se0.a=alpha parameter for sens among controls
#se0.b=beta parameter for sens among controls

sp1.a <- 91.2
sp1.b <- 3.86
sp0.a <- sp1.a
sp0.b <- sp1.b

#calculate unadjusted RR for COVID death among the entire cohort
unadj.rr_covid <- (a1 / (total_exp)) / (b1 / (total_unexp))
print(unadj.rr_covid)

unadj.or_covid <- (a1 / (total_exp - a1)) / (b1 / (total_unexp - b1))
print(unadj.or_covid)

#as a sense check: calculate unadjusted RR for all-cause death among the entire cohort
unadj.rr_all_cause_death <- (death_exp / total_exp) / (death_unexp / total_unexp)
print(unadj.rr_all_cause_death)

# Function for summary level correction -----------------------------------

pba.summary.exp.rr <- function(a1, b1, c1, d1, total_exp, total_unexp, death_exp, death_unexp, se1.a, se1.b, se0.a, se0.b, sp1.a, sp1.b, sp0.a, sp0.b, SIMS){
  
  #among all deaths
  n_case <- a1 + b1 #total covid deaths
  n_ctrl <- c1 + d1 #total non-covid deaths
  n_exp  <- a1 + c1 #total deaths among ICS
  n_unexp <- b1 + d1 #total deaths among LABA/LAMA
  
  print(c(n_case, n_ctrl, n_exp, n_unexp))
  
  # draw sensitivities from beta distribution
  se <- rbeta(I, se1.a, se1.b)
  sp <- rbeta(I, sp1.a, sp1.b)
  
  # Create a data frame with only the first 1000 draws for both variables
  data_for_plot <- data.frame(Sensitivity = se[1:2000], Specificity = sp[1:2000])
  
  # Create the scatter plot
  plot <- ggplot(data_for_plot, aes(x = Sensitivity, y = Specificity)) +
    geom_point(alpha = 0.2) + 
    theme_minimal()
  
  # Save the plot
  ggsave(filename = file.path(Graphdir, "QBA", "copd_death_w1", "sens_spec_scatter_death_summary.png"),
         plot = plot,
         width = 6,
         height = 6)
 
  # calculate bias-adjusted cell frequencies: only among deaths
  ac1 <- round((a1 - death_exp*(1-sp))/(se - (1-sp))) #bias-adjusted cases, exposed 
  cc1 <- round(death_exp - ac1) #bias-adjusted cases, unexposed
  bc1 <- round((b1 - death_unexp*(1-sp))/(se - (1-sp))) #bias-adjusted controls, exposed
  dc1 <- round(death_unexp - bc1) #bias-adjusted controls, unexposed
  
  #calculate prevalence of exposure in cases and controls, accounting for sampling error
  PrevD_exp <- rbeta(I,ac1,cc1)
  PrevD_unexp <- rbeta(I,bc1,dc1)
  
  #calculate PPV and NPV of exposure classification in cases and controls
  #these must be calculated separately for exposed and unexposed
  PPV_exp <- (se*PrevD_exp)/((se*PrevD_exp)+(1-sp)*(1-PrevD_exp))
  PPV_unexp <- (se*PrevD_unexp)/((se*PrevD_unexp)+(1-sp)*(1-PrevD_unexp))
  NPV_exp <- (sp*(1-PrevD_exp))/((1-se)*PrevD_exp+sp*(1-PrevD_exp))
  NPV_unexp <- (sp*(1-PrevD_unexp))/((1-se)*PrevD_unexp+sp*(1-PrevD_unexp))
  
  #calculate the expected number of cases among exp and unexp
  #this incorporates error from the misclassification process by using binomial trials
  ab <- rbinom(I,ac1,PPV_exp) + rbinom(I,cc1,1-NPV_exp) ## exposed with outcome
  cb <- n_exp-ab ##exposed without outcome
  bb <- rbinom(I,bc1,PPV_unexp) + rbinom(I,dc1,1-NPV_unexp) ###unexposed with outcome
  db <- n_unexp-bb ###unexposed without outcome

  # add back in the people who survived in the unexposed group
  ac <- ab   
  bc <- bb   
  cc <- total_exp - ac  # add back in the people who survived in the exposed group
  dc <- total_unexp - bc # add back in the people who survived in the unexposed group
  
  flag <- (ac<=0)|(bc<=0)|(cc<=0)|(dc<=0)
  
  #calculate bias adjusted RR with second source of uncertainty
  # this incorporates the uncertainty stemming from the misclassification process
  rr_bb <- (ac/(ac+cc))/(bc/(bc+dc))
  or_bb <- (ac*dc)/(bc*cc)
  
  # calculate bias-adjusted risk ratios, third source of uncertainty, bias-adjusted standard error
  se_bb <- sqrt(1/ab+1/bb-1/(ab+cb)-1/(bb+db))
  se_bb_or <- sqrt(1/ac+1/bc+1/cc+1/dc)
  #draw from normal distribution to incorporate standard error
  z <- rnorm(I)
  rr_bb_cb <- exp(log(rr_bb) - (z*se_bb))
  or_bb_cb <- exp(log(or_bb) - (z*se_bb_or))
  
  # bind vectors as dataframe
  bound <- data.frame(se, sp, ac1, bc1, cc1, dc1, ac, bc, cc, dc, ab, bb, cb, db, PrevD_exp, PrevD_unexp, PPV_exp, PPV_unexp, NPV_exp, NPV_unexp, rr_bb_cb, or_bb_cb)
  
  # housekeeping - retain rows with positive corrected cell freqs
  summary_pba <<- bound %>%
    filter(ac > 0) %>%
    filter(bc > 0) %>%
    filter(cc > 0) %>%
    filter(dc > 0) %>%
    filter(ab > 0) %>%
    filter(bb > 0) %>%
    filter(cb > 0) %>%
    filter(db > 0)
  
  rr <- c(quantile(summary_pba$rr_bb_cb,c(0.025,0.5,0.975)), I-length(summary_pba$rr_bb_cb), length(summary_pba$rr_bb_cb))
  or <- c(quantile(summary_pba$rr_bb_cb,c(0.025,0.5,0.975)), I-length(summary_pba$or_bb_cb), length(summary_pba$or_bb_cb))
  names(rr) <- c("2.5th %tile","50th %tile", "97.5th %tile","Impossible values", "Simulations")
  names(or) <- c("2.5th %tile","50th %tile", "97.5th %tile","Impossible values", "Simulations")
  return(list(rr = rr, or = or,  data = bound))
  
}

start_time <- Sys.time()

qba_result <- pba.summary.exp.rr(a1 = a1,
                                 b1 = b1,
                                 c1 = c1,
                                 d1 = d1,
                                 total_exp = total_exp,
                                 total_unexp = total_unexp,
                                 death_exp = death_exp,
                                 death_unexp = death_unexp,
                                 se1.a = se1.a,
                                 se1.b = se1.b,
                                 se0.a = se0.a,
                                 se0.b = se0.b,
                                 sp1.a = sp1.a,
                                 sp1.b = sp1.b,
                                 sp0.a = sp0.a,
                                 sp0.b = sp0.b,
                                 SIMS = I
                                 )

write_parquet(qba_result$data, file.path(Tables, "QBA", samples))

# Access the main result
print(round(qba_result$rr, 3))
print(round(qba_result$or, 3))

# Create dataframes for rr and or
rr_df <- data.frame(rr = qba_result$rr)
or_df <- data.frame(or = qba_result$or)

# Combine rr_df and or_df into qba_df
qba_df <- cbind(rr_df, or_df)
#rename columns
colnames(qba_df) <- c("rr_death_summary", "or_death_summary")
#rename rownames
rownames(qba_df) <- c("2.5th %tile","50th %tile", "97.5th %tile", "Impossible values", "Simulations")

write_parquet(qba_df, file.path(Tables, "QBA", paste0("qba_death_summary_results", output_ext, ".parquet")))

end_time <- Sys.time()
run_time <- end_time - start_time
print(run_time)

#extract the unit from runtime
run_time_unit <- attr(run_time, "units")

# add stats to excel file to save the number of iterations, date, time, outcome, summary level vs record level
#import runtime.xlsx
runtime_df <- read.xlsx(file.path(Tables, "QBA", "runtime.xlsx"))

time_stats <- data.frame(
  outcome = "death",
  analysis = "summary level_func",
  sa = output_ext,
  niter = format(I, scientific = FALSE),
  date = format(Sys.Date(), "%Y-%m-%d"),
  time = format(Sys.time(), "%H:%M:%S"),  # Format Sys.time() to display only hour, minute, second
  run_time = run_time,
  run_time_unit = run_time_unit
)

#add time stats to runtime.df
runtime_df <- rbind(runtime_df, time_stats)
#add this line to runtime.xlsx
write.xlsx(runtime_df, file.path(Tables, "QBA", "runtime.xlsx"))


# Plotting key metrics ----------------------------------------------------
median_value_rr <- median(qba_result$data$rr_bb_cb, na.rm = TRUE)
pal <- c("unadjusted RR" = palette[4], "Median RR after QBA" = palette[10])
pal2 <- c("unadjusted RR" = "dashed", "Median RR after QBA" = "dotted")
plot <- ggplot(data = data.frame(x = qba_result$data$rr_bb_cb), aes(x = x)) +
  geom_histogram(binwidth = 0.001, fill = palette[6]) +
  geom_vline(aes(xintercept = unadj.rr_covid, color = "unadjusted RR", linetype = "unadjusted RR"), size = 0.9) +
  geom_vline(aes(xintercept = median_value_rr, color = "Median RR after QBA", linetype = "Median RR after QBA"), size = 0.9) +
  xlab("Risk Ratio") +
  ylab("Frequency") +
  ggtitle(paste0("Risk ratios adjusted for outcome misclassification (n = ", format(length(summary_pba$rr_bb_cb), scientific = FALSE), ")"))+
  xlim(0, 5) +
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
file_path <- file.path(Graphdir, "QBA", "copd_death_w1", paste0("adjusted_RR_SL", output_ext, ".png"))
ggsave(file_path, plot, width = 8, height = 4)

median_value_or <- median(qba_result$data$or_bb_cb, na.rm = TRUE)
pal <- c("unadjusted OR" = palette[4], "Median OR after QBA" = palette[10])
pal2 <- c("unadjusted OR" = "dashed", "Median OR after QBA" = "dotted")
plot <- ggplot(data = data.frame(x = qba_result$data$or_bb_cb), aes(x = x)) +
  geom_histogram(binwidth = 0.001, fill = palette[6]) +
  geom_vline(aes(xintercept = unadj.or_covid, color = "unadjusted OR", linetype = "unadjusted OR"), size = 0.9) +
  geom_vline(aes(xintercept = median_value_or, color = "Median OR after QBA", linetype = "Median OR after QBA"), size = 0.9) +
  xlab("Odds Ratio") +
  ylab("Frequency") +
  ggtitle(paste0("Odds ratios adjusted for outcome misclassification (n = ", format(length(summary_pba$or_bb_cb), scientific = FALSE), ")"))+
  xlim(0, 5) +
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
file_path <- file.path(Graphdir, "QBA", "copd_death_w1", paste0("adjusted_OR_SL", output_ext, ".png"))
ggsave(file_path, plot, width = 8, height = 4)

plot <- ggplot(data = data.frame(x = qba_result$data$se), aes(x = x)) +
  geom_histogram(binwidth = 0.001, fill = palette[2]) +
  xlab("Sampled Sensitivity") +
  ylab("Frequency") +
  ggtitle(paste0("Sampled sensitivity (alpha = ", se1.a, ", beta = ", se1.b, ", n = ", format(I, scientific = FALSE),  ")")) +
  xlim(0, 1) +
  theme_classic()+
  theme(plot.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
file_path <- file.path(Graphdir, "QBA", "copd_death_w1", paste0("sampled_se_SL_pba", output_ext, ".png"))
ggsave(file_path, plot, width = 8, height = 4)

plot <- ggplot(data = data.frame(x = qba_result$data$sp), aes(x = x)) +
  geom_histogram(binwidth = 0.001, fill = palette[2]) +
  xlab("Sampled Specificity") +
  ylab("Frequency") +
  ggtitle(paste0("Sampled specificity (alpha = ", sp1.a, ", beta = ", sp1.b, ", n = ", format(I, scientific = FALSE),  ")")) +
  xlim(0, 1) +
  theme_classic() +
  theme(plot.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
file_path <- file.path(Graphdir, "QBA", "copd_death_w1", paste0("sampled_sp_SL_pba", output_ext, ".png"))
ggsave(file_path, plot, width = 8, height = 4)

plot_ppv <- ggplot(data = qba_result$data, aes(x = PPV_exp, color = "ICS")) +
  geom_density(adjust = 1, fill = NA, size = 1) +
  geom_density(aes(x = PPV_unexp, color = "LABA/LAMA"), adjust = 1, fill = NA, size = 1) +
  xlab("Sampled PPV") +
  ylab("Density") +
  ggtitle("Sampled PPV by treatment group") +
  xlim(0, 1) +
  theme_classic()+
  theme(plot.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 11)) +
  scale_color_manual(values = c("ICS" = palette[9], "LABA/LAMA" = palette[4]), labels = c("ICS", "LABA/LAMA")) +
  guides(color = guide_legend(title = NULL, override.aes = list(fill = c(palette[9], palette[4]))))
file_path <- file.path(Graphdir, "QBA", "copd_death_w1", paste0("sampled_PPV_SL_pba_death", output_ext, ".png"))
ggsave(file_path, plot_ppv, width = 8, height = 4)

plot_npv <- ggplot(data = qba_result$data, aes(x = NPV_exp, color = "ICS")) +
  geom_density(adjust = 1, fill = NA, size = 1) +
  geom_density(aes(x = NPV_unexp, color = "LABA/LAMA"), adjust = 1, fill = NA, size = 1) +
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
file_path <- file.path(Graphdir, "QBA", "copd_death_w1", paste0("sampled_NPV_SL_pba_death", output_ext, ".png"))
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
file_path <- file.path(Graphdir, "QBA", "copd_death_w1", paste0("sampled_PPV_NPV_SL_pba_death", output_ext, ".png"))
ggsave(file_path, final_plot, width = 10, height = 4)


plot <- ggplot(data = qba_result$data, aes(x = PrevD_exp, color = "ICS")) +
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
file_path <- file.path(Graphdir, "QBA", "copd_death_w1", paste0("sampled_Prev_SL_pba_death", output_ext, ".png"))
ggsave(file_path, plot, width = 8, height = 4)

}

#list inputfiles
inputfile <- c("copd_wave1_60d_iptw.parquet", "SA_copd_wave1_60d_iptw_no_triple.parquet", "SA_copd_wave1_6m_iptw.parquet", "SA_copd_wave1_60d_iptw_all.parquet")
conventional_result <- c("cox_log_regression_estimates.parquet", "SA_cox_log_regression_estimates_no_triple.parquet", "SA_cox_log_regression_estimates_6m.parquet", "SA_cox_log_regression_estimates_all.parquet")
samples <- c("qba_death_summary_full_sample.parquet", "SA_qba_death_summary_full_sample_no_triple.parquet", "SA_qba_death_summary_full_sample_6m.parquet", "SA_qba_death_summary_full_sample_all.parquet")
output_ext <- c("", "_no_triple",  "_6m", "_all")
# 
# inputfile <- c("SA_copd_wave1_60d_iptw_no_triple.parquet")
# conventional_result <- c("SA_cox_log_regression_estimates_no_triple.parquet")
# samples <- c("qba_death_summary_full_sample.parquet")
# output_ext <- c("")
#run the function

mapply(death_summary_qba, inputfile, conventional_result, samples, output_ext)

# Author: Marleen Bokern
# Date: 12/1623
# Purpose: Combine all effect estimates from 

packages <- c("tidyverse", "ggplot2", "gt", "arrow", "patchwork", "forestplot", "cowplot")
installed_packages <- packages %in% rownames(installed.packages())

if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

lapply(packages, library, character.only = TRUE)

setwd(Datadir_copd)

estimates_conventional <- read_parquet(file.path(Tables, "QBA", "cox_log_regression_estimates.parquet")) %>% dplyr::select(-contains("coef_"), -contains("se_"), -contains("res_p"), -contains("or"), -contains("att_weight"), -contains("log")) %>% filter(!stringr::str_detect(outcome_event, "any_"))

hosp_record <- read_parquet(file.path(Tables, "QBA", "qba_hosp_record_results.parquet"))
death_record <- read_parquet(file.path(Tables, "QBA", "qba_death_record_results.parquet"))
#remove or_ columns
hosp_record <- hosp_record %>% select(-contains("or_"))
death_record <- death_record %>% select(-contains("or_"))

results_table <- estimates_conventional %>% 
  pivot_longer(!outcome_event, names_to = "weight", values_to = "value") %>%
  mutate(weightlabel = case_when(
    stringr::str_detect(weight, "ate_weight_stab") ~ "IPTW",
    stringr::str_detect(weight, "unadj") ~ "Unweighted",
    TRUE ~ weight )) %>%
  mutate(weight = case_when(
    stringr::str_detect(weight, "hr") ~ "point_estimate",
    stringr::str_detect(weight, "ci_lower") ~ "ci_lower",
    stringr::str_detect(weight, "ci_upper") ~ "ci_upper",
    TRUE ~ weight ))  %>%
  pivot_wider(values_from = value, names_from = weight)

# Combine all effect estimates and keep rownames
estimates <- bind_cols(hosp_record, death_record)

#remove the 4th and 5th row
estimates <- estimates[-c(4),]

#transpose the data frame
estimates <- t(estimates)
#convert estimates to df
estimates <- as.data.frame(estimates)

#display v4 not in scientific notation
estimates$V4 <- format(estimates$V4, scientific = FALSE)

#rename the columns v1, v2, v3
estimates <- estimates %>%
  rename(
    point_estimate = V2,
    ci_lower = V1,
    ci_upper = V3, 
    simulations = V4
  )

#rbind the new row for death as row 8
final_estimates <- estimates

#add rownames as a column
final_estimates <- final_estimates %>% 
  rownames_to_column("estimate")

#add conventional estimates
final_estimates <- final_estimates %>% 
  bind_rows(results_table)

final_estimates <- final_estimates %>% 
  mutate(reg = case_when(
    stringr::str_detect(tolower(as.character(estimate)), "hr_") ~ "HR",
    TRUE ~ estimate)) %>% 
  mutate(outcome = case_when(
    stringr::str_detect(tolower(as.character(estimate)), "death") ~ "Death",
    stringr::str_detect(tolower(as.character(outcome_event)), "death") ~ "Death",
    stringr::str_detect(tolower(as.character(estimate)), "hosp") ~ "Hosp",
    stringr::str_detect(tolower(as.character(outcome_event)), "hes") ~ "Hosp",
    TRUE ~ estimate)) %>%
  mutate(qba = case_when(
    stringr::str_detect(tolower(as.character(estimate)), "record") ~ "Record-level PBA",
    TRUE ~ "Conventional")) %>% 
  mutate(estimate_type = case_when(
    stringr::str_detect(tolower(as.character(estimate)), "ate_") ~ "IPTW",
    stringr::str_detect(tolower(as.character(weightlabel)), "iptw") ~ "IPTW",
    TRUE ~ "Unweighted"))

forest_table <- final_estimates %>%
  mutate(across(
    c(point_estimate, ci_lower, ci_upper),
    as.numeric)) %>%
  mutate(across(
    c(point_estimate, ci_lower, ci_upper), 
    ~ ifelse(is.na(.x), "", sprintf("%.2f", round(.x, 2))))) %>% 
  mutate(estimate_lab = paste0(point_estimate, 
                               ifelse(ci_lower == "" | ci_upper == "",
                                      "",
                                      paste0(" (", ci_lower, " - ", ci_upper, ")"))))

forest_table <- forest_table %>%
  mutate(outcome = fct_recode(outcome,
                              "COVID-19 hospitalisation" = "Hosp",
                              "COVID-19 death" = "Death"
  ))

# Split out hospitalisations and deaths
forest_table_hosp <- forest_table %>%
  filter(outcome == "COVID-19 hospitalisation")

forest_table_hosp <- forest_table_hosp %>%
  arrange(outcome, estimate_type, qba) %>% 
  mutate(
    point_estimate = as.numeric(point_estimate),
    lower_CI = as.numeric(ci_lower),
    upper_CI = as.numeric(ci_upper)) %>%
  mutate(order = case_when(
    qba == "Conventional" & estimate_type == "Unweighted" ~ 1,
    qba == "Conventional" & estimate_type == "IPTW" ~ 2,
    qba == "Record-level PBA" & estimate_type == "Unweighted" ~ 3,
    qba == "Record-level PBA" & estimate_type == "IPTW" ~ 4
  )) %>% #create orderno that is reverse of order
  mutate(orderno = (nrow(forest_table_hosp) + 1) - order)


#sort by orderno
forest_table_hosp <- forest_table_hosp[order(forest_table_hosp$orderno),]

# Generate forest plot for hospitalizations
forest_plot_hosp <- ggplot(forest_table_hosp, aes(x = point_estimate, y = orderno)) +
  geom_point(aes(x = point_estimate), shape = 16, size = 5) +
  geom_linerange(aes(xmin = lower_CI, xmax = upper_CI)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Effect estimate", y = "") +
  theme_classic() +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20)) +
  scale_x_continuous(trans = 'log10', breaks = c(1, 2, 4, 8)) +
  coord_cartesian(xlim = c(0.5, 8))

est_hosp <- ggplot(forest_table_hosp, aes(y = orderno)) +
  geom_text(aes(x = 0, label = qba), hjust = 0, size = 22/.pt, fontface = "bold") +
  geom_text(aes(x = 1.2, label = estimate_type), hjust = 0, size = 22/.pt) +
  geom_text(aes(x = 2.2, label = estimate_lab), hjust = 0, size = 22/.pt) +
  theme_void() +
  coord_cartesian(xlim = c(0, 3))

# Split out deaths
forest_table_death <- forest_table %>%
  filter(outcome == "COVID-19 death")

forest_table_death <- forest_table_death %>%
  arrange(outcome, estimate_type, qba) %>% 
  mutate(
    point_estimate = as.numeric(point_estimate),
    lower_CI = as.numeric(ci_lower),
    upper_CI = as.numeric(ci_upper)) %>%
  mutate(order = case_when(
    qba == "Conventional" & estimate_type == "Unweighted" ~ 1,
    qba == "Conventional" & estimate_type == "IPTW" ~ 2,
    qba == "Record-level PBA" & estimate_type == "Unweighted" ~ 3,
    qba == "Record-level PBA" & estimate_type == "IPTW" ~ 4
  )) %>% 
  mutate(orderno = (nrow(forest_table_hosp) + 1) - order)

# Generate forest plot for death
forest_plot_death <- ggplot(forest_table_death, aes(x = point_estimate, y = orderno)) +
  geom_point(aes(x = point_estimate), shape = 16, size = 5) +
  geom_linerange(aes(xmin = lower_CI, xmax = upper_CI)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Effect estimate", y = "") +
  theme_classic() +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20)) +
  scale_x_continuous(trans = 'log10', breaks = c(1, 2, 4, 8))+
  coord_cartesian(xlim = c(0.5, 8))

est_death <- ggplot(forest_table_death, aes(y = orderno)) +
  #geom_text(aes(x = 0, label = qba), hjust = 0, size = 22/.pt) +
  #geom_text(aes(x = 0, label = estimate_type), hjust = 0, size = 22/.pt) +
  geom_text(aes(x = 0, label = estimate_lab), hjust = 0, size = 22/.pt) +
  theme_void() +
  coord_cartesian(xlim = c(0, 3))

# define titles
title <- ggdraw() +
  draw_label("COVID-19 hospitalisation", fontface = 'bold', size = 20, x = 0.45, hjust = 0.5, vjust = 1) + 
  draw_label("COVID-19 death", fontface = 'bold', size = 20, x = 0.85, hjust = 0.5, vjust = 1)

print(title)

layout <- c(
  area(t = 0, l = 0, b = 30, r = 180), # Reduced right boundary to increase spacing
  area(t = 0, l = 180, b = 30, r = 260), # Increased left boundary for more spacing
  area(t = 0, l = 280, b = 30, r = 340), # Further increased left boundary for spacing
  area(t = 0, l = 330, b = 30, r = 410) # Final area with increased left boundary for spacing
)

# Assuming est_hosp, forest_plot_hosp, est_death, and forest_plot_death are defined elsewhere
final_forest <- title /(est_hosp + forest_plot_hosp + est_death + forest_plot_death + plot_layout(design = layout)) +
  plot_layout(heights = c(0.2, 1))

# Save plot for hospitalizations
file_path_hosp <- file.path(Graphdir, "QBA", "forest_plot_full_main_cox.png")
ggsave(filename = file_path_hosp, plot = final_forest, width = 22, height = 7, units = "in", dpi = 1000)

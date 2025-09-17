# Author: Marleen Bokern
# Date: 10/2023
# Purpose: Master script for analysis of hospitalisations and deaths (wave 1) in COPD population
#suppress warnings

#set filepaths using filepaths.R
setwd(Github_folder)

#import and format analytic file
source("01copd_import_analytic_file_func.R")
setwd(Github_folder)
source("01copd_import_analytic_file_60d_hosp.R")

#make baseline tables
setwd(Github_folder)
source("02copd_baseline_table_w1_func.R")
setwd(Github_folder)
source("02copd_baseline_table_w1_no_triple.R")
setwd(Github_folder)
source("02copd_baseline_table_w1_all.R")

setwd(Github_folder)
suppressWarnings(source("04copd_iptw.R"))
setwd(Github_folder)
suppressWarnings(source("04copd_iptw_all.R"))
setwd(Github_folder)
suppressWarnings(source("04copd_iptw_no_triple.R"))
setwd(Github_folder)
suppressWarnings(source("04copd_iptw_6m.R"))
setwd(Github_folder)
source("05copd_cox_hosp_death.R")
setwd(Github_folder)
source("05copd_cox_hosp_death_6m.R")
setwd(Github_folder)
source("05copd_cox_hosp_death_all.R")
setwd(Github_folder)
source("05copd_cox_hosp_death_no_triple.R")
setwd(Github_folder)
source("05copd_cox_hosp_death_all_missing_covid.R")


#Death summary level QBA

setwd(Github_folder)
source("copd_pba_deaths_w1_summary.R")
setwd(Github_folder)
source("SA_copd_pba_deaths_w1_summary_no_triple.R")
setwd(Github_folder)
source("SA_copd_pba_deaths_w1_summary_all_missing_covid.R")


setwd(Github_folder)
source("copd_pba_deaths_w1_summary_func.R")
setwd(Github_folder)
source("SA_copd_pba_deaths_w1_summary_all_missing_covid.R")
setwd(Github_folder)
source("SA_copd_pba_deaths_w1_summary_no_missing_covid.R")

 
setwd(Github_folder)
source("copd_pba_deaths_w1_record_func.R")##run
setwd(Github_folder)
source("SA_copd_pba_deaths_w1_record_parallel_all_missing_covid.R")##run


setwd(Github_folder)
source("copd_pba_hosp_w1_summary.R")
setwd(Github_folder)
source("copd_pba_hosp_w1_summary_dist2_func.R")


setwd(Github_folder)
source("copd_pba_hosp_w1_record_parallel.R")
setwd(Github_folder)
source("SA_copd_pba_hosp_w1_record_parallel_no_triple.R")

setwd(Github_folder)
source("copd_pba_forest_plot_func.R")
setwd(Github_folder)
source("copd_pba_forest_plot_func_main_paper.R")
setwd(Github_folder)
source("SA_copd_pba_forest_plot_no_triple.R")
setwd(Github_folder)
source("copd_pba_forest_plot_no_missing_covid.R")
setwd(Github_folder)
source("copd_pba_forest_plot_all_missing_covid.R")




# Author: Marleen Bokern
# Date: 03/2023
# Purpose: basic data formatting for analytic file, generate binary treatment group variable

packages <- c("tidyverse", "dplyr", "arrow", "openxlsx")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

lapply(packages, library, character.only = TRUE)
#######################################################

setwd(Datadir_copd)

#read in csv dataset
df <- read.csv("copd_analytic_file_w1_60d.csv", colClasses = c("patid" = "character"), header = T)

#count baseline ics and baseline control
sum(df$baseline_ics == 1, na.rm = TRUE)
sum(df$baseline_control == 1, na.rm = TRUE)

#count pos_covid_test_date not empty where baseline_ics is 1 or baseline_control is 1
length(unique(df$patid[(df$baseline_ics == 1 | df$baseline_control == 1) & df$pos_covid_test_date != ""]))
length(unique(df$patid[(df$baseline_ics == 1 | df$baseline_control == 1) & df$covid_hes_date != ""]))
length(unique(df$patid[(df$baseline_ics == 1 | df$baseline_control == 1)]))


file_path <- paste0(Projectdir_stata, "patients_missing_ons.csv")
df_missing_ons <- read.csv(file_path, colClasses = c("patid" = "character"), header = T)

# Create missing_ons variable indicating presence in df_missing_ons
df$missing_ons <- ifelse(df$patid %in% df_missing_ons$patid, 1, 0)

df$death_date_cprd <- df$deathdate
df <- subset(df, select = -deathdate)

# Replace values in the 'death' column based on conditions considering empty strings as missing values
df$death_date <- ifelse(df$death_date_ons != "", df$death_date_ons,
                        ifelse(df$death_date_ons == "" & df$death_date_cprd != "", df$death_date_cprd, NA))

# Reformatting variables --------------------------------------------------

#create binary treatment group variable. Set to 1 if ICS group, set to 0 if LABA/LAMA
df$treatgroup <- ifelse(df$baseline_ics == 1, 1, NA)
df$treatgroup[df$baseline_control == 1] <- 0
df$treatgroup <- as.factor(df$treatgroup)

df_baseline <- df %>% dplyr::select(c("patid", "treatgroup", "baseline_ics", "baseline_control", "baseline_triple"))
sum(df_baseline$baseline_ics == 1, na.rm = TRUE)
sum(df_baseline$baseline_control == 1, na.rm = TRUE)
sum(df_baseline$baseline_triple == 1, na.rm = TRUE)
sum(df_baseline$baseline_ics == 1 & df_baseline$baseline_control == 1)

write.xlsx(df_baseline, "copd_baseline_exposure.xlsx")

#convert covariates coded as date variables to dates 
date_cols <- c("regstart", "death_date_cprd", "death_date", "regend", "dob", "do35bday", "lcd", "enddate", "diabetes_date", "hypertension_date", "past_asthma_date", "cvd_date", "allcancers_date", "kidney_date", "immunosuppression_date", "smokdate", "flu_vacc_date", "pneumo_vacc_date", "pos_covid_test_date", "covid_hes_date", "covid_death_date", "death_date_ons")
df[date_cols] <- lapply(df[date_cols], as.Date, format = "%d%b%Y")

df <- subset(df, death_date_ons >= as.Date("2020-03-01") | is.na(death_date_ons)) #this excludes a couple of people who had an ons_death date before 01 mar 2020

#create binary variables for each covariate
df$diabetes_present <- factor(ifelse(!is.na(df$diabetes_date), "Yes", "No"))
df$hypertension_present <- factor(ifelse(!is.na(df$hypertension_date), "Yes", "No"))
df$cvd_present <- factor(ifelse(!is.na(df$cvd_date), "Yes", "No"))
df$allcancers_present <- factor(ifelse(!is.na(df$allcancers_date), "Yes", "No"))
df$past_asthma_present <- factor(ifelse(!is.na(df$past_asthma_date), "Yes", "No"))
df$covid_present <- factor(ifelse(!is.na(df$pos_covid_test_date), "Yes", "No"))
df$kidney_present <- factor(ifelse(!is.na(df$kidney_date), "Yes", "No"))
df$immunosuppression_present <- factor(ifelse(!is.na(df$immunosuppression_date), "Yes", "No"))
df$flu_vacc_present <- factor(ifelse(!is.na(df$flu_vacc_date), "Yes", "No"))
df$pneumo_vacc_present <- factor(ifelse(!is.na(df$pneumo_vacc_date), "Yes", "No"))

df$exacerbations <- ifelse(is.na(df$exacerbations), 0, df$exacerbations)
df$exacerb_present <- factor(ifelse(df$exacerbations != 0, "Yes", "No"))

#recode smoking status
df$smok <- recode(df$smokstatus,
                  "current smoker" = "Current smoking",
                  "current/ex-smoker" = "Current/Former smoking",
                  "ex-smoker" = "Former smoking",
                  "nonspecified - depends on quantity" = "Unclear",
                  .default = "Unclear")

# Convert "smok" to a factor
df$smok <- factor(df$smok)

# Convert "imd" to a factor
df$imd <- factor(df$imd)
df$imd <- fct_explicit_na(df$imd, "Missing")

# Define the value labels for treatgroup, gender and smoking
labels_treat <- c("LABA/LAMA", "ICS")
labels_gender <- c("Male", "Female")
#labels_smoking <- c("Current smoking", "Current/Former smoking", "Former smoking", "Unclear")
labels_ics_ever <- c("ICS")
labels_control_ever <- c("LABA/LAMA")

# Convert the variable to a factor with the specified labels
df$treat <- ifelse(df$treatgroup == 0, "LABA/LAMA",
                   ifelse(df$treatgroup == 1, "ICS", NA))
df$gender <- factor(df$gender, levels = c(1,2), labels = labels_gender)
#df$smok <- factor(df$smokstatus, labels = labels_smoking)
df$ics_ever <- factor(df$ics_ever, labels = labels_ics_ever)
df$control_ever <- factor(df$control_ever, labels = labels_control_ever)



# Recode ethnicity 
table(df$eth5)

# create a new variable 'eth' in the data frame
df$eth <- NA_integer_

# replace eth values based on eth5 values
df$eth[df$eth5 == "0. White"] <- 0
df$eth[df$eth5 == "1. South Asian"] <- 1
df$eth[df$eth5 == "2. Black"] <- 2
df$eth[df$eth5 == "4. Mixed"] <- 3
df$eth[df$eth5 == "5. Not Stated" | df$eth5 == "3. Other" | df$eth5 == ""] <- 4

#This also works
# df$ethn <- ifelse(df$eth5 == "", NA, df$eth5)
# df$ethn <- factor(df$ethn)
# df$ethn <- ifelse(df$eth5 == 6, NA, df$ethn)

labels_eth <- c("White", "South Asian", "Black", "Mixed", "Unknown")
df$eth <- factor(df$eth, levels = c(0,1,2,3,4), labels = labels_eth)

# recode BMI as categorical
df$bmicat <- NA

# Categorize BMI values
df$bmicat[df$bmi < 18.5] <- 1
df$bmicat[df$bmi >= 18.5 & df$bmi < 25] <- 2
df$bmicat[df$bmi >= 25 & df$bmi < 30] <- 3
df$bmicat[df$bmi >= 30] <- 4

df$bmicat[is.na(df$bmicat)] <- 2

# Assign category labels
df$bmicat <- factor(df$bmicat, levels = c(1, 2, 3, 4),
                    labels = c("Underweight (<18.5)", "Normal (18.5-24.9)",
                               "Overweight (25-29.9)", "Obese (>=30)"))

# Create a less granular categorization
df$obese4cat <- NA
df$obese4cat[df$bmicat %in% c("Underweight (<18.5)", "Normal (18.5-24.9)", "Overweight (25-29.9)")] <- 1
df$obese4cat[df$bmicat == "Obese I (30-34.9)"] <- 2
df$obese4cat[df$bmicat == "Obese II (35-39.9)"] <- 3
df$obese4cat[df$bmicat == "Obese III (40+)"] <- 4

# Assign category labels
df$obese4cat <- factor(df$obese4cat, levels = c(1, 2, 3, 4),
                       labels = c("No record of obesity", "Obese I (30-34.9)",
                                  "Obese II (35-39.9)", "Obese III (40+)"))

# code end of follow ups
#create new variables denoting the length of time until the outcome of interest. 
#generate a timeout variable for all cause deaths as a sense check.
df$timeout1 <- df$pos_covid_test_date
df$timeout2 <- df$covid_hes_date
df$timeout3 <- df$covid_death_date
df$timeout_death_any <- df$death_date
df$timeout_hes_any <- df$hes_date

#WORK WITH OUTCOME DATES
timeout_variables <- c("timeout1", "timeout2", "timeout3", "timeout_death_any")

for (var_name in timeout_variables) {
  
  df[[var_name]][is.na(df[[var_name]])] <- df$enddate[is.na(df[[var_name]])]
  
  #replace the timeout variables with the end date if the timeout is after the end date
  df[[var_name]][df[[var_name]] > df$enddate] <- df$enddate[df[[var_name]] > df$enddate]
  
  # Replace values greater than "2020-08-31" with "2020-08-31"
  df[[var_name]][df[[var_name]] > as.Date("2020-08-31")] <- as.Date("2020-08-31")
  
  # Calculate the sum of values that are equal to "2020-08-31"
  sum_result <- sum(df[[var_name]] == as.Date("2020-08-31"), na.rm = TRUE)
  
  cat("Sum for", var_name, ":", sum_result, "\n")
} 

# Generate a variable that is 1 if outcome is present, and 0 if outcome date is empty
df$pos_covid_test_present <- factor(ifelse(!is.na(df$pos_covid_test_date), "Yes", "No"))
df$pos_covid_test_present <- ifelse(is.na(df$pos_covid_test_date), 0, 1)

df$covid_hes_present <- factor(ifelse(!is.na(df$covid_hes_date), "Yes", "No"))
df$covid_hes_present <- ifelse(is.na(df$covid_hes_date), 0, 1)

df$covid_death_present <- factor(ifelse(!is.na(df$covid_death_date), "Yes", "No"))
df$covid_death_present <- ifelse(is.na(df$covid_death_date), 0, 1)

df$any_death_present <- ifelse(df$death_date < as.Date("2020-09-01"), 1, 0)
df$any_death_present <- ifelse(df$covid_death_present == 1, 1, df$any_death_present)
df$any_death_present <- ifelse(is.na(df$any_death_present), 0, df$any_death_present)

#assert that there are more all cause deaths than covid deaths
sum(df$any_death_present == 1 & df$death_date_ons < as.Date("2020-09-01"), na.rm = TRUE) > sum(df$covid_death_present == 1 & df$death_date_ons < as.Date("2020-09-01"), na.rm = TRUE)

# Set the time origin to 01 Mar 2020
df$time_origin <- as.Date("2020-03-01")

# Convert the exit time variable to a survival time object
df$time_origin = as.numeric(df$time_origin, "%d%b%Y")

df[timeout_variables] <- lapply(df[timeout_variables], as.numeric, format = "%d%b%Y")

df$timeinstudy1 <- df$timeout1 - df$time_origin
df$timeinstudy2 <- df$timeout2 - df$time_origin
df$timeinstudy3 <- df$timeout3 - df$time_origin
df$timeinstudy_death_any <- df$timeout_death_any - df$time_origin

df <- df %>% dplyr::select(-c("yob", "mob", "day", "dob", "yo35bday", "do35bday", "lcd", "diabetes_date", "hypertension_date", "past_asthma_date", "cvd_date", "allcancers_date", "kidney_date", "immunosuppression_date", "smokdate", "flu_vacc_date", "pneumo_vacc_date", "dobmi"))

#if these columns exist, remove them

# Filter the dataframe based on conditions
missing_ons_laba <- df %>%
  filter(missing_ons == 1 & treat == "LABA/LAMA") %>%
  pull(patid) %>%
  unique()

# Filter the dataframe based on conditions
missing_ons_ics <- df %>%
  filter(missing_ons == 1 & treat == "ICS") %>%
  pull(patid) %>%
  unique()

# Filter the dataframe based on conditions
surv_laba <- df %>%
  filter(missing_ons == 0 & any_death_present == 0  & treat == "LABA/LAMA") %>%
  summarise(unique_patids = n_distinct(patid))

# View the count of unique patids
surv_laba$unique_patids

# Filter the dataframe based on conditions
surv_ics <- df %>%
  filter(missing_ons == 0 & any_death_present == 0  & treat == "ICS") %>%
  summarise(unique_patids = n_distinct(patid))

# View the count of unique patids
surv_ics$unique_patids

length(unique((df$patid[df$treat == "ICS"])))
length(unique((df$patid[df$treat == "LABA/LAMA"])))

df <- df %>%
  mutate(
    hos_pre_death = ifelse(
      !is.na(covid_hes_date) &
        !is.na(death_date) &
        covid_hes_date < death_date, 1, 0))

#export as parquet file
arrow::write_parquet(df, "copd_wave1_60d.parquet")
length(unique((df$patid[df$treat == "ICS"])))
length(unique((df$patid[df$treat == "LABA/LAMA"])))

length(unique((df$patid[df$treat == "ICS" & df$covid_death_present == 1])))
length(unique((df$patid[df$treat == "LABA/LAMA" & df$covid_death_present == 1])))

# #drop variables aside from regend, treatgroup, treat, patid, and timeinstudy1, timeinstudy2, timeinstudy3, timeinstudy_death_any, pos_covid_test_present, covid_hes_present, covid_death_present, any_death_present, hos_pre_death, death_date, death_date_cprd, death_date_ons, death, baseline_ics, baseline_control, baseline_triple
# 
# df <- df %>% dplyr::select(c("patid", "enddate", "regend", "treatgroup", "treat", "timeinstudy1", "timeinstudy2", "timeinstudy3", "timeinstudy_death_any", "pos_covid_test_present", "covid_hes_present", "covid_death_present", "any_death_present", "hos_pre_death", "death_date", "death_date_cprd", "death_date_ons", "death", "baseline_ics", "baseline_control", "baseline_triple"))
# 
# df2 <- read_parquet("copd_wave1_60d_hosp.parquet")
# 
# #drop variables aside from regend, treatgroup, treat, patid, and timeinstudy1, timeinstudy2, timeinstudy3, timeinstudy_death_any, pos_covid_test_present, covid_hes_present, covid_death_present, any_death_present, hos_pre_death, death_date, death_date_cprd, death_date_ons, death, baseline_ics, baseline_control, baseline_triple
# 
# df2 <- df2 %>% dplyr::select(c("patid", "enddate", "regend", "treatgroup", "treat", "timeinstudy1", "timeinstudy2", "timeinstudy3", "timeinstudy_death_any", "pos_covid_test_present", "covid_hes_present", "covid_death_present", "any_death_present", "hos_pre_death", "death_date", "death_date_cprd", "death_date_ons", "baseline_ics", "baseline_control"))
# 
# df <- merge(df, df2[, c("patid", "treatgroup")], by = "patid", all.x = TRUE)


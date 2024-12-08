/*=========================================================================
DO FILE NAME:		    10cr_copd_all_drugs.do

AUTHOR:					Marleen Bokern

VERSION:				v1

DATE VERSION CREATED: 	11/2022

DATASETS CREATED:       all_drugs_copd.dta
						
DESCRIPTION OF FILE:	appends all inhaled medication files

*=========================================================================*/
/************************************************************************
**************************************************************************
CREATE DATASET CONTAINING ALL INHALED MED PRESCRIPTIONS
- append all inhaled medication drug issue files
************************************************************************
*************************************************************************/
clear all 
capture log close
log using $Logdir/10cr_copd_all_drugs.log , replace

cd $Datadir_copd

***append datasets containing individual prescriptions
use "${file_stub}_budesonide_single_clean.dta"
gen budesonide_single = 1

append using "${file_stub}_budesonide_laba_clean.dta", gen(budesonide_laba)
append using "${file_stub}_laba_single_clean.dta", gen(laba_single)
append using "${file_stub}_laba_lama_clean.dta", gen(laba_lama)
append using "${file_stub}_lama_single_clean.dta", gen(lama_single)
append using "${file_stub}_triple_budesonide_clean.dta", gen(triple_budesonide)
sort patid issuedate

***generate new variable for class of inhaler
gen class = 0 if budesonide_single == 1
replace class = 1 if budesonide_laba == 1
replace class = 2 if laba_single == 1
replace class = 3 if laba_lama == 1
replace class = 4 if lama_single == 1
replace class = 5 if triple_budesonide == 1


label define class 0 "budesonide_single" 1 "budesonide_laba" 2 "laba_single" 3 "laba_lama" 4 "lama_single" 5 "triple_budesonide"
label values class class
drop budesonide_single budesonide_laba laba_single lama_single laba_lama triple_budesonide
drop termfromemis productname

compress
save "10cr_copd_all_drugs_budesonide.dta", replace

/*********************************************************
1. INITIATIONS
-initiations are defined as a new episode where the end of the last episode is at least 365 days after the end of the last exposure end date
**********************************************************/
sort patid issuedate

**drop all unnecessary variables, only keep patid, end of exposure date and issue date
keep patid issuedate end_exp class

****create counter for each prescription by date and class
sort patid issuedate
bysort patid class (issuedate): gen episode = _n

***reshape data
***end_exp = issuedate + exposure length
rename (issuedate end_exp) (date1 date2)
reshape long date, i(patid class episode) j(startend) 

***Encode the start and end for ranking the order
gen startend2 = 0 if startend == 1
replace startend2 = 1 if startend == 2

***sort by patid and class, date (and within date, start and end), subract number of preceding end dates from number of preceding starts
by patid class (date startend2), sort: gen int in_proc = sum(startend == 1) - sum(startend == 2)
replace in_proc = 1 if in_proc > 1

by patid class (date): gen block_num = 1 if in_proc == 1 & in_proc[_n - 1] != 1 //**1s now denote start of new block
by patid class (date): replace block_num = sum(block_num) //***count number of blocks per patient

***assert that a patient's first date is a start and the last is an end
by patid class block_num (date), sort: assert startend == 1 if _n == 1
by patid class block_num (date): assert startend == 2 if _n == _N

***keep only start and end of blocks
by patid class block_num (date): keep if _n == 1 | _n == _N

**************************************************************************************
drop episode in_proc startend2 
reshape wide date, i(patid class block_num) j(startend)

bysort patid class (date1): gen episode = _n
keep patid episode date1 date2 class

***merge in first reg date 
merge m:1 patid using "copd_patient_included_all.dta", keepusing(regstart enddate) keep(match) nogen

***tag people whose regstart is less than 12 months before start of episode
gen lag = date1 - regstart
gen excl = 1 if lag < 365

****generate variable init denoting initiation
sort patid class date1
bysort patid class (date1): gen init = 1 if (episode == 1 & excl ==.) | date1 - date2[_n - 1] > 365 & excl ==.

***generate variable to denote month of episode start
g dm = mofd(date1)
format dm %tm

bysort patid class (date1): assert date2 < date1[_n + 1]

***assess initation over time 
histogram dm if init == 1 & date1 > td(01dec2017), discrete frequency title(Initiations of any inhaled medication per month)
graph export $Graphdir/initiation_any_budesonide.jpg, replace

compress
save "copd_treatment_episodes_initiation_budesonide", replace
/************************************************************************
*************************************************************************
ASSESS DISCONTINUATIONS, TAKING INTO ACCOUNT 60 DAY ALLOWABLE GAP
*************************************************************************
*************************************************************************/
*60 days gap
reshape long date, i(patid episode class) j(startend)
format date %td

gen startend2 = 0 if startend == 1
replace startend2 = 1 if startend == 2

*** flag startdates within 60d of previous date (sorted within patient on the date, starts before ends if same date)
by patid class (date startend2), sort: gen no_gap = 1 if startend == 1 & (date - date[_n - 1] <= 60) 

*** flag end dates that are followed by start date within 60d
replace no_gap = 1 if startend == 2 & no_gap[_n + 1] == 1

*flag if episode belongs to the same treatment block as other episodes
egen no_gap_max = max(no_gap), by (patid class episode)

keep if no_gap ==. //*** flags starts and ends of episodes
**keep if (no_gap_max == 1 & no_gap ==.) | (no_gap==. & no_gap_max ==.)

drop no_gap no_gap_max episode startend2 dm lag excl init

*change the episode no as rx for reshaping the wide form
egen rx = seq(), f(1) b(2)

reshape wide date, i(patid class rx) j(startend)

g dm = mofd(date1)
format dm %tm

replace date2 = date2 + 60
bysort patid class (date1): assert date2 < date1[_n + 1]

***flag if registration end or death is before calculated discontinuation date
gen flag = 1 if enddate < date2

g dm_disc = mofd(date2) if flag != 1
format dm_disc %tm

****assess discontinuation over time 
histogram dm_disc if dm_disc < tm(2021m4) & dm_disc > tm(2017m4) & flag != 1, xlabel(#10) discrete frequency title(Discontinuations per month defined as 60 days after end of exposure)
graph export $Graphdir/discontinuation60_any_budesonide.jpg, replace

compress
save "copd_treatment_episodes_60d_budesonide", replace

/************************************************************************
*************************************************************************
ASSESS DISCONTINUATIONS, DEFINED AS 6M AFTER ISSUEDATE
*************************************************************************
*************************************************************************/

clear all
use "10cr_copd_all_drugs_budesonide.dta", replace

***merge in registration start and end dates
drop substancestrength dosageid
merge m:1 patid using "copd_patient_included_all.dta", keepusing(regstart enddate) keep(match) nogen

assert issuedate <= enddate

*** get discontinuation date = 6m post last prescription
gen disc_date6m = issuedate + 183
format disc_date6m %td

***flag prescriptions where the registration end is before the calculated discontinuation date
gen flag = 1 if enddate < disc_date6m

**drop all unnecessary variables, only keep patid, class, end of exposure date and issue date
keep patid issuedate disc_date6m class flag

****create counter for each prescription
sort patid issuedate
bysort patid class (issuedate): gen episode = _n

rename (issuedate disc_date6m) (date1 date2)
reshape long date, i(patid episode class flag) j(startend) //when using reshape, each combination of the variables in i() should be unique

*Encode the start and end for ranking the order
gen startend2 = 0 if startend == 1
replace startend2 = 1 if startend == 2

***sort by patid, date (and within date, start and end), subract number of preceding end dates from number of preceding starts
by patid class (date startend2 episode), sort: gen int in_proc = sum(startend == 1) - sum(startend == 2)
replace in_proc = 1 if in_proc > 1

***generate new variable "block_num" that indicates start of a new block
bysort patid class (date episode): gen block_num = 1 if in_proc == 1 & in_proc[_n - 1] != 1 //**1s now denote start of new block

***count number of blocks per patient
bysort patid class (date episode): replace block_num = sum(block_num)

***assert that a patient's first date is a start and the last is an end
by patid class block_num (date episode), sort: assert startend == 1 if _n == 1
by patid class block_num (date episode): assert startend == 2 if _n == _N

***keep only start and end of blocks
sort patid class block_num date episode
bysort patid class block_num (date episode): keep if _n == 1 | _n == _N

**************************************************************************************
drop in_proc episode startend2

bysort patid: replace flag = 1 if flag[_n + 1] == 1 & startend == 1
reshape wide date, i(patid flag block_num class) j(startend)

bysort patid class (date1): gen episode = _n 

***merge in registration start and end dates
merge m:1 patid using "copd_patient_included_all.dta", keepusing(regstart enddate) keep(match) nogen

bysort patid class (episode date1): assert date2 < date1[_n + 1] if date2 !=.

***generate variable for month and year of episode start
g dm = mofd(date1)
format dm %tm

***replace discontinuation as missing if registration end or death is before calculated discontinuation date
gen disc_date = date2
format disc_date %td
replace disc_date =. if date2 > enddate

g dm_disc = mofd(disc_date) if flag != 1
format dm_disc %tm

histogram dm_disc if dm_disc < tm(2021m4) & dm_disc > tm(2017m4) & flag != 1, xlabel(#10) discrete frequency title(Discontinuations per month defined as 6m after issue date)
graph export $Graphdir/discontinuation6m_any_budesonide.jpg, replace

compress
save "copd_treatment_episodes_6m_budesonide", replace

log close 
clear all

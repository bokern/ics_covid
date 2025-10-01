/*=========================================================================
DO FILE NAME:		    11cr_copd_treatment_combinations60d.do

AUTHOR:					Marleen Bokern

VERSION:				v1

DATE VERSION CREATED: 	11/2022

DATASETS CREATED:       ${file_stub}_treatment_ep_for_combination60d
						${file_stub}_treatment_eps_w1_full60d
						${file_stub}_treatment_baseline_w160d
						
DESCRIPTION OF FILE:	uses code from Kate Mansfield to generate treatment episodes

*=========================================================================*/

clear all
capture log close
log using "$Logdir\11cr_copd_treatment_combinations60d.log", replace

***use file with treatment episodes according to 60d discontinuation definition
use "$Datadir_copd\copd_treatment_episodes_60d"
unique patid if class == 3

expand 2, gen(dup)

*dup = 0 has starts, dup = 1 has ends
gen date = date1 if dup == 0
replace date = date2 if dup == 1
format (date) %td

gen onoff = 1 if dup == 0
replace onoff = 0 if dup == 1

*flag denotes an exposure end date that occurs after the registration end date
replace flag =. if dup == 0

drop rx date1 date2 dm_disc enddate regstart dm

compress

/************************************************************************
generate on and off indicators for each drug class
*************************************************************************/

label list
gen ics_single = 1 if class == 0
gen ics_laba = 1 if class == 1
gen laba_single = 1 if class == 2
gen laba_lama = 1 if class == 3
gen lama_single = 1 if class == 4
gen triple_inhaler = 1 if class == 5

* generate a state var for each different thing that can happen (start/stop)
local expClasses ics_single ics_laba laba_single lama_single laba_lama triple_inhaler
di "`expClasses'"
foreach a in `expClasses' {
	generate on_`a'= 0
	generate off_`a'= 0
	replace on_`a' = 1 if `a' == 1 & onoff == 1
	replace off_`a' = 1 if `a' == 1 & onoff == 0
}

* collapse on patid and date taking the max value of each indicator
collapse (max) on_* off_* flag, by(patid date)

*drop rows that happen after the end date
merge m:1 patid using "$Datadir_copd\copd_Patient_included_all.dta", keepusing(enddate) keep(match) nogen
assert (on_ics_single & on_ics_laba & on_laba_single & on_lama_single & on_laba_lama & on_triple_inhaler) == 0 if flag == 1
assert (off_ics_single == 1 | off_ics_laba == 1 | off_laba_single == 1 | off_lama_single == 1 | off_laba_lama == 1 | off_triple_inhaler == 1) if flag == 1

assert date > enddate if flag == 1

drop if flag == 1
drop enddate flag

/************************************************************************
update indicator variables
*************************************************************************/
* generate current status variable
* and identify base state of vars based on switching on drug
foreach b in `expClasses' {
	gen current_`b'= 0
	recode current_`b' 0 = 1 if on_`b' == 1
}

sort patid date

* update current status for each drug based on what happened before
foreach c in `expClasses' {
	replace current_`c' = current_`c'[_n - 1] if 	/// set the `drug' flag to the same as the record 
									/// (i.e. previous date) before 
		patid == patid[_n - 1] &			/// IF this is same patient
		current_`c'[_n - 1] != 0 &			/// the record before isn't 0 - we don't update 
									/// a switched on record (don't change a 1 to a 0)		
		off_`c'!= 1						// 	`drug' not switched off on this date
	rename current_`c' `c'
} /*end foreach c in `expClasses' */

gen triple = 1 if (ics_laba == 1 & lama_single == 1)|(ics_single == 1 & laba_lama == 1)|(ics_single == 1 & laba_single == 1 & lama_single == 1)|(ics_laba == 1 & laba_lama == 1)
replace triple = 1 if triple_inhaler == 1
replace triple = 0 if triple ==.
drop if date > td(01mar2021) //drop changes that happen after the end of follow-up

gen ics_group = 1 if ics_single == 1 & laba_single == 1 | ics_laba == 1 | triple == 1
replace ics_group = 0 if ics_group ==.
gen control_group = 1 if laba_single == 1 & lama_single == 1 | laba_lama == 1 
replace control_group = 0  if ics_group == 1 | control_group ==.

gen ics_only = 1 if ics_single == 1 & ics_group == 0 & control_group == 0 & lama_single == 0 & laba_single == 0
replace ics_only = 0 if ics_only ==.
gen laba_only = 1 if laba_single == 1 & ics_group == 0 & control_group == 0 & ics_single == 0 & lama_single == 0
replace laba_only = 0 if laba_only ==.
gen lama_only = 1 if lama_single == 1 & ics_group == 0 & control_group == 0 &  ics_single == 0 & laba_single == 0
replace lama_only = 0 if lama_only ==.
gen ics_lama = 1 if ics_single == 1 & ics_group == 0 & control_group == 0 & lama_single == 1 & laba_single == 0
replace ics_lama = 0 if ics_lama ==.
gen no_med = 1 if ics_single == 0 & ics_laba == 0 & laba_single == 0 & lama_single == 0 & laba_lama == 0 & triple == 0
replace no_med = 0 if no_med ==.

unique patid 
return list

***check every combination has been accounted for
by patid date: gen sum = sum(ics_group + control_group + ics_only + laba_only + lama_only + ics_lama + no_med)
assert sum == 1
drop sum

compress
save "$Datadir_copd\\${file_stub}_treatment_eps_full60d", replace

************************************************************************************************************************************************************
*** find baseline exposure --> keep last observation before index date
drop if date > td(01mar2020)
bysort patid: gen last = _N
bysort patid (date): gen ep = _n
keep if ep == last
gen baseline_ics = 1 if ics_group == 1
replace baseline_ics = 0 if baseline_ics ==.
gen baseline_control = 1 if control_group == 1
replace baseline_control = 0 if baseline_control ==.
gen baseline_triple = 1 if triple == 1
replace baseline_triple = 0 if baseline_triple ==.

assert ics_group == 1 if triple == 1
assert baseline_ics == 1 if baseline_triple == 1

unique patid if triple == 1
unique patid if (ics_single == 1 & laba_single == 1 | ics_laba == 1) & triple == 0
unique patid if ics_group == 1
unique patid if control_group == 1
unique patid if ics_only == 1
unique patid if laba_only == 1
unique patid if lama_only == 1
unique patid if ics_lama == 1
unique patid if no_med == 1
drop last ep on_ics_single on_ics_laba on_laba_single on_lama_single on_laba_lama off_ics_single off_ics_laba off_laba_single off_lama_single off_laba_lama date

compress
save "$Datadir_copd\\${file_stub}_treatment_baseline_60d", replace

clear all

***count people who change treatment groups during follow up
use "$Datadir_copd\\${file_stub}_treatment_eps_full60d"
merge m:1 patid using "$Datadir_copd\\${file_stub}_treatment_baseline_60d", keepusing(baseline_control baseline_ics baseline_triple)

sort patid date
unique patid if date > td(01mar2020) & date < td(01sep2020) //anyone with a treatment change during follow-up
unique patid if date > td(01mar2020) & date < td(01sep2020) & (baseline_ics == 1 & ics_group == 0) // people in ics group at baseline who leave ics group
unique patid if date > td(01mar2020) & date < td(01sep2020) & (baseline_control == 1 & control_group == 0) // people in control group at baseline who leave control group
unique patid if date > td(01mar2020) & date < td(01sep2020) & (baseline_ics == 0 & ics_group == 1) // people not in ics group at baseline who enter ics group
unique patid if date > td(01mar2020) & date < td(01sep2020) & (baseline_control == 0 & control_group == 1) // people not in control group at baseline who enter control group
unique patid if date > td(01mar2020) & date < td(01sep2020) & (baseline_ics == 1 & control_group == 1) // people in ics group at baseline who enter control group
unique patid if date > td(01mar2020) & date < td(01sep2020) & (baseline_control == 1 & ics_group == 1) // people in control group at baseline who enter ics group

sort patid date
by patid: gen ics_ever = 1 if baseline_ics == 1 | (ics_group == 1) & date > td(01mar2020)
by patid: ereplace ics_ever = max(ics_ever)
by patid: gen control_ever = 1 if baseline_control == 1 | (control_group == 1) & date > td(01mar2020)
by patid: ereplace control_ever = max(control_ever)
by patid: gen triple_ever = 1 if baseline_triple == 1 | (triple == 1) & date > td(01mar2020)
by patid: ereplace triple_ever = max(triple_ever)

compress
save "$Datadir_copd\\${file_stub}_treatment_eps_60d", replace
log close 
clear all
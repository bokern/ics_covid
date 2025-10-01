/*=========================================================================
DO FILE NAME:		    09cr_copd_treatment_ep.do

AUTHOR:					Marleen Bokern

VERSION:				v1

DATE VERSION CREATED: 	12/2022

DATASETS CREATED:       "${file_stub}_`drugclass'_clean"
						
DESCRIPTION OF FILE:	loops through and cleans drug issue files. deals with duplicates, gaps, overlaps

*=========================================================================*/
/************************************************************
HOUSEKEEPING
*************************************************************/
clear all 
capture log close
log using "$Logdir\09cr_copd_treatment_ep.log", replace
cd "$Datadir_copd"

putexcel set "$Datadir_copd/medication_summary.xlsx", replace
putexcel A2 = "n_pat"
putexcel A3 = "n_rx"

/************************************************************
SET UP LOOPING THROUGH ALL DRUG CLASSES
*************************************************************/
local drug ics_single ics_laba laba_lama laba_single lama_single triple_therapy
*local drug budesonide_single budesonide_laba 
*local drug triple_budesonide

di `"`drug'"'

local ncol = 2

di "`col'"
foreach drugclass in `drug'{

local col: word `ncol' of `c(ALPHA)'

putexcel `col'1 = `"`drugclass'"'

cd "$Datadir_copd"

noi di "Cleaning drug issue file for `drugclass'"
noi di "Import drug issue file for `drugclass'"

***use file containing all prescriptions for each drug class (drug issue files merged with code list)
use "${file_stub}_`drugclass'_checked", clear

***merge variable comes from merge with dosage lookup file
capture drop _merge

noi di "Merge with included patient file, keep matches"
merge m:1 patid using "${file_stub}_Patient_included_all.dta", keepusing(dob enddate)

***people in using only are included in the cohort but have no prescription. People in master only dont meet inclusion criteria
keep if _merge == 3
drop _merge

noi di "Drop if issuedate is after enddate"
count if issuedate > enddate
drop if issuedate > enddate

/********************************************************************
DEAL WITH DATES
*********************************************************************/
noi di "Drop prescriptions before March 2017 and after April 2021"
***drop prescriptions before mar 2017 & after may 2021
count if issuedate < td(01mar2017)
drop if issuedate < td(01mar2017) 
count if issuedate >= td(01may2021)
drop if issuedate >= td(01may2021)

**check for weird issue and enterdates
noi di "Check for missing issuedates and issuedates before dob'"
count if missing(issuedate) //***0
count if issuedate < dob
drop if issuedate < dob 

***replace implausible quantities as missing
noi di "Replace quantities < 10 and > 1000 as missing"
gen quantity2 = quantity
replace quantity2 =. if quantity2 < 10 |quantity2 > 1000

***replace values > 100 and < 7 in duration variable as missing
noi di "Replace durations < 7 and > 100 as missing"
gen length = duration
// Count implausible durations before replacement
count if length < 7 | length > 100
local n_implausible_length = r(N)
local total_before_length = r(N) + r(N) // Total observations before dropping
count
local total_obs = r(N)
local pct_implausible_length = (`n_implausible_length' / `total_obs') * 100

replace length =. if length < 7 //***replaces all the 0s
replace length =. if length > 100
tab length 

unique patid
return list
putexcel `col'2 = (r(unique))
putexcel `col'3 = (r(N))

/********************************************************************
DEAL WITH DUPLICATES
*********************************************************************/

***find duplicates in terms of patid, prodcodeid, issuedate and duplicates in terms of all variables
noi di "Tag duplicates in terms of patid, issuedate and prodcodeid, and duplicates in terms of all variables"
duplicates tag patid issuedate prodcodeid, generate(duplicate)

noi di "Count duplicates in terms of patid, issuedate and prodcodeid"
tab duplicate

// Calculate and save the number of duplicates for this drug class
count if duplicate > 0
local n_duplicates = r(N)
noi di "Number of duplicate records for `drugclass': `n_duplicates'"
putexcel A19 = "n_duplicates" // Added row for number of duplicates
putexcel `col'19 = (`n_duplicates') // Added line to save number of duplicates


duplicates tag, gen(duplicate_all)

noi di "Count duplicates in terms of all variables"
tab duplicate_all
count if duplicate != duplicate_all

***flag those where there are differences in other variables (other than patid, date, product)
gen dup_flag = 1 if duplicate != duplicate_all

noi di "Generate counter = _N and counter2 = _n"
sort patid issuedate prodcodeid
by patid issuedate prodcodeid : gen counter = _N
by patid issuedate prodcodeid : gen counter2 = _n

***fill in missing lengths, quantities and doses from duplicates in terms of pid, date, product
noi di "Fill in missing lengths from duplicates in terms of pid, date, product"
sort patid issuedate prodcode
gen length2 = length
by patid issuedate prodcode: replace length2 = length2[_n - 1] if length2 ==. & dup_flag == 1 & dup_flag[_n - 1] == 1
by patid issuedate prodcode: replace length2 = length2[_n + 1] if length2 ==. & dup_flag == 1 & dup_flag[_n + 1] == 1

noi di "fill in missing quantities from duplicates in terms of pid, date, product"
sort patid issuedate prodcode
by patid issuedate prodcode: replace quantity2 = quantity2[_n - 1] if quantity2 ==. & dup_flag == 1 & dup_flag[_n - 1] == 1
by patid issuedate prodcode: replace quantity2 = quantity2[_n + 1] if quantity2 ==. & dup_flag == 1 & dup_flag[_n + 1] == 1

noi di "fill in missing daily doses from duplicates in terms of pid, date, product"
sort patid issuedate prodcode
gen daily_dose2 = daily_dose
by patid issuedate prodcode: replace daily_dose2 = daily_dose2[_n - 1] if daily_dose2 ==. & dup_flag == 1 & dup_flag[_n - 1] == 1
by patid issuedate prodcode: replace daily_dose2 = daily_dose2[_n + 1] if daily_dose2 ==. & dup_flag == 1 & dup_flag[_n + 1] == 1

**generate exposure duration from quantity and daily dose. Replace implausible values as missing
noi di "Generate exposure duration from quantity and daily dose. Replace implausible values as missing"
gen days = quantity2 / daily_dose2
replace days = round(days,1)
order days, after(duration)
tab days
replace days =. if days > 100
replace days =. if days < 7

***combine length and quantity of "duplicate" prescriptions 
noi di "Combine length and quantity of duplicate prescriptions "
bysort patid issuedate prodcodeid : egen sum_quant = sum(quantity2) if duplicate >= 1
bysort patid issuedate prodcodeid : egen sum_length = sum(length2) if duplicate >= 1
bysort patid issuedate prodcodeid : egen sum_days = sum(days) if duplicate >= 1
order (sum_quant sum_length sum_days), after(quantity)

**keep "last" observations where duplicate
noi di "keep last observations where duplicate"
keep if counter == counter2
drop counter2 dup_flag duplicate_all duplicate

***insert quantity and length to sum_quant or sum_length variables where these are missing
noi di "insert quantity and length to sum_quant or sum_length variables where these are missing"
replace sum_quant = quantity2 if sum_quant ==.
replace sum_quant =. if sum_quant < 10
replace sum_length = length2 if sum_length ==.
replace sum_length =. if sum_length < 7 //*when adding lengths in row above, missing values get added together as 0
replace sum_days =. if sum_days == 0
replace sum_days = days if sum_days ==.
****sum_quant and sum_length are now the "correct" lengths and quantities

tab sum_quant
tab sum_length
tab sum_days

noi di "count if duration == 0 & days ==. "
count if duration == 0 & days ==. 

gen diff_duration = days - sum_length
tab diff_duration

noi di "count if diff_duration != 0 "
count if diff_duration != 0 
order diff_duration, after(days)

***
local metric sum_days
tabstat `metric', stat(mean sd median p25 p75) save
return list
matrix list r(StatTotal)
matrix stats = r(StatTotal)
putexcel `col'4 = matrix(stats) 
putexcel A4 = "duration_mean"
putexcel A5 = "duration_sd"
putexcel A6 = "duration_median"
putexcel A7 = "duration_p25"
putexcel A8 = "duration_p75"

local metric sum_length
tabstat `metric', stat(mean sd median p25 p75) save
return list
matrix list r(StatTotal)
matrix stats = r(StatTotal)
putexcel `col'9 = matrix(stats) 
putexcel A9 = "`metric'_mean"
putexcel A10 = "`metric'_sd"
putexcel A11 = "`metric'_median"
putexcel A12 = "`metric'_p25"
putexcel A13 = "`metric'_p75"

local metric sum_quant
tabstat `metric', stat(mean sd median p25 p75) save
return list
matrix list r(StatTotal)
matrix stats = r(StatTotal)
putexcel `col'14 = matrix(stats) 
putexcel A14 = "`metric'_mean"
putexcel A15 = "`metric'_sd"
putexcel A16 = "`metric'_median"
putexcel A17 = "`metric'_p25"
putexcel A18 = "`metric'_p75"
***

by drugsubstancename, sort : tabstat sum_quant duration days daily_dose2 sum_length, stat(n mean sd median p25 p75 min max)

/********************************************************************
DEAL WITH GAPS BETWEEN PRESCRIPTIONS
*********************************************************************/

**generate prescription counter by patid
noi di "generate prescription counter by patid"
sort patid issuedate
bysort patid : gen prescription = _n

*generate variable of gaps between prescription issue dates
noi di "generate variable of gaps between prescription issue dates"
by patid (issuedate): gen pres_gap = issuedate - issuedate[_n - 1]
order pres_gap, after(issuedate)
replace pres_gap =. if prescription == 1 //*** should be 0

count if pres_gap ==. 
assert pres_gap ==. if prescription == 1
assert pres_gap !=. if prescription != 1
hist pres_gap if pres_gap < 200, percent xtitle(Days between presciptions) //***0s are people with 2 different products on same day
graph export $Graphdir/`drugclass'_pres_gap.jpg, replace

tabstat duration days pres_gap sum_length, stat(n median p25 p75 mean sd min max) 

by drugsubstancename, sort : tabstat duration days pres_gap length, stat(n median p25 p75 mean sd min max)
twoway (hist sum_length if sum_length < 500, percent lcolor(gs12) fcolor(gs12)) (hist days if days < 500, percent fcolor(red) lcolor(red)), legend(off) xtitle("Duration (red: quantity/daily dose)")
graph export $Graphdir/`drugclass'_length.jpg, replace

g dm = mofd(issuedate)
format dm %tm

histogram dm [fweight = counter], discrete percent ylabel(#6) xtitle(Calendar month) xlabel(#10) fcolor(ltblue) lcolor(none)
graph export $Graphdir/`drugclass'_hist.jpg, replace
order dosageid, last

***use duration variable where available (length variable = duration but has values <7 and >100 as missing). use calculated quantity/daily dose where duration not available. where neither are available, use median for each dru
noi di "Generate exp_dur as duration where plausible, quantity/DD or median length for each drug"
gen exp_dur = sum_length
replace exp_dur = sum_days if exp_dur ==.
// Count prescriptions before median imputation
count
local total_prescriptions = r(N)
count if exp_dur ==.
local n_missing_before_median = r(N)


egen median_dur = median(sum_length), by(drugsubstancename)
replace exp_dur = (median_dur) if exp_dur ==.

// Count prescriptions where median was imputed
local n_median_imputed = `n_missing_before_median'
local pct_median_imputed = (`n_median_imputed' / `total_prescriptions') * 100

drop median_dur

assert exp_dur !=.

// Export implausible lengths and median imputation statistics to Excel
putexcel A20 = "n_implausible_length"
putexcel `col'20 = (`n_implausible_length')
putexcel A21 = "pct_implausible_length" 
putexcel `col'21 = (`pct_implausible_length')
putexcel A22 = "n_median_imputed"
putexcel `col'22 = (`n_median_imputed')
putexcel A23 = "pct_median_imputed"
putexcel `col'23 = (`pct_median_imputed')


*generate exposure end dates
noi di "Generate exposure end dates"
gen end_pres = issuedate + exp_dur
format (end_pres) %td
order exp_dur end_pres, after(issuedate)

assert end_pres > issuedate

**generate variable to denote the number of days of overlap
noi di "Generate variable to denote the number of days of overlap"
sort patid prescription
gen overlap_len = end_pres[_n - 1] - issuedate
replace overlap_len =. if prescription == 1
order overlap_len, after(end_pres)
assert overlap_len ==. if prescription == 1
assert overlap_len !=. if prescription != 1

*****variable with number of overlap days to add to each prescription
noi di "Generate variable with number of overlap days to add to each prescription"
gen sum_overlap1 = overlap_len
sort patid prescription

noi di "adds up all previous overlaps/gaps in each patient, limit sum_overlap1 so cannot be negative"
bysort patid (prescription) : replace sum_overlap1 = sum(overlap_len) //***this just adds up all previous overlaps/gaps in each patient --> can become negative (which is impossible irl)
bysort patid (prescription) : replace sum_overlap1 = 0 if sum_overlap1 < 0 //****limits sum_overlap1 so cannot be negative
order sum_overlap1, after(overlap_len)

gen sum_overlap = cond(_n == 1 | sum_overlap1 == 0, sum_overlap, .) //***sum_overlap is 0 if its first prescription or if there was a gap/no overlap 
order sum_overlap, after(overlap_len)
replace sum_overlap = sum_overlap[_n - 1] + overlap_len if missing(sum_overlap) //***adds overlap_len to previously accrued sum_overlap

***cap maximum allowed overlap at 90 days
replace sum_overlap = 90 if sum_overlap > 90

***generate end_exp which is the end date of the prescription plus any previously accrued overlap
gen end_exp = end_pres + sum_overlap

noi di "Generate episode denoting if prescription belongs to same episode as previous Rx"
sort patid prescription
gen episode = 1 if end_exp[_n - 1] >= issuedate & prescription != 1
order episode, after(overlap_len)

drop sum_overlap1

***/
compress
save "${file_stub}_`drugclass'_clean", replace
local ncol = `ncol' + 1

}

log close 
clear all
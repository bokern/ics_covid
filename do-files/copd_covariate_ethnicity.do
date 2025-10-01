/*=========================================================================
DO FILE NAME:			copd_covariate_ethnicity.do

AUTHOR:					Marleen Bokern, adpated from Emily Herrett

VERSION:				v1

DATE VERSION CREATED: 	11/2022

DATASETS CREATED:       copd_Observation_ethnicity.dta
						
DESCRIPTION OF FILE:	ethnicity codelist applied

*=========================================================================*/

/*******************************************************************************
>> HOUSEKEEPING
*******************************************************************************/
clear all

*****run global file

****export cl_ethnicity as txt file
use "$Codelistsdir/cl_ethnicity.dta"
export delim "$Codelistsdir/cl_ethnicity.txt"

clear all


capture log close
log using $Logdir/covariate_ethnicity_copd.log, replace

cd "$Datadir_copd"

****ethnicity codelist was originally from Angel, adapted using code from Emily Herrett
****loop thorugh all observation files to generate dataset of all ethnicity codes

foreach file of numlist 1/$no_Observation {
	noi di "Merging ethnicity observations, File `file'"
    use "$Copd_aurum_extract\\${file_stub}_Extract_Observation_`file'", clear
	keep patid medcodeid eventdate 
	merge m:1 patid using "${file_stub}_Patid_list_included_all.dta", keep(match) nogen
    merge m:1 medcodeid using "$Codelistsdir/cl_ethnicity.dta", keep(match) nogen
	drop observations
	if `file' == 1 {
		save "${file_stub}_Observation_ethnicity.dta", replace
	}
	if `file' > 1 {
		append using "${file_stub}_Observation_ethnicity.dta"
	}
	compress
	save "${file_stub}_Observation_ethnicity.dta", replace
}

use "${file_stub}_Observation_ethnicity.dta"
unique patid

drop if eventdate > td(30apr2021)

duplicates drop //**0

*drop duplicates obs where the same ethnicity was recorded on the same system date
duplicates drop patid eventdate medcodeid, force //**9,521 


**turn sysdate format into years
gen eventyear = year(eventdate)
sum eventyear //1996-2014

**tag people with the same ethnicity recorded in the same year
duplicates tag patid eventyear medcodeid, gen(duplicate)
tab duplicate if duplicate > 0

lab var duplicate "duplicate ethnicity recorded in the same year"

**gen indicators for obs of ethnicity per patient
sort patid eventdate
bysort patid: gen count = [_n]
bysort patid: gen total = [_N]
sum count total

gen totalobs = total

**ADD UP ETHNICITIES
bysort patid eth5: gen eth5count = [_N]
tab eth5count

bysort patid eth16: gen eth16count = [_N]
tab eth16count

*generate variables for each of the 5 categories indicating how often each category was recorded for each patient
gen white5count = eth5count if eth5 == 0
gen sa5count = eth5count if eth5 == 1
gen black5count = eth5count if eth5 == 2
gen other5count = eth5count if eth5 == 3
gen mixed5count = eth5count if eth5 == 4
gen notstated5count = eth5count if eth5 == 5

tab white5count if count == 1
tab white5count if count == total

gen british16count = eth16count if eth16 == 1
gen irish16count = eth16count if eth16 == 2
gen otherwhite16count = eth16count if eth16 == 3
gen whitecarib16count = eth16count if eth16 == 4
gen whiteaf16count = eth16count if eth16 == 5
gen whiteasian16count = eth16count if eth16 == 6
gen othermixed16count = eth16count if eth16 == 7
gen indian16count = eth16count if eth16 == 8
gen pak16count = eth16count if eth16 == 9
gen bangla16count = eth16count if eth16 == 10
gen otherasian16count = eth16count if eth16 == 11
gen carib16count = eth16count if eth16 == 12
gen african16count = eth16count if eth16 == 13
gen otherblack16count = eth16count if eth16 == 14
gen chinese16count = eth16count if eth16 == 15
gen other16count = eth16count if eth16 == 16
gen notstated16count = eth16count if eth16 == 17


**MAKE COUNTS CONSTANT WITHIN PATIENT 
*if a patient has 7 rows total, 5 with white and 2 with sa: previously white5count would be 5 in rows where eth5 ==white and missing where sa. After this loop all 7 rows of this patient say white5count==5 and sa5count ==2

local p  "white sa black other mixed notstated"
foreach i of local p {
sort patid count
replace `i'5count = `i'5count[_n - 1] if `i'5count[_n] ==. & `i'5count[_n - 1]!=. & patid[_n] == patid[_n - 1] & totalobs > 1
gsort patid -count
replace `i'5count = `i'5count[_n - 1] if `i'5count[_n] ==. & `i'5count[_n - 1]!=. & patid[_n] == patid[_n - 1] & totalobs > 1
}

local p "british irish otherwhite whitecarib whiteaf whiteasian othermixed indian pak bangla otherasian carib african otherblack chinese other notstated"
foreach i of local p {
sort patid count
replace `i'16count = `i'16count[_n - 1] if `i'16count[_n] ==. & `i'16count[_n - 1] !=. & patid[_n] == patid[_n - 1] & totalobs > 1
gsort patid -count
replace `i'16count = `i'16count[_n - 1] if `i'16count[_n] ==. & `i'16count[_n - 1] !=. & patid[_n] == patid[_n - 1] & totalobs > 1
}


**DUMMY FOR WHETHER ONLY ETHNICITY IS NOT STATED
gen notstatedonly = 0
replace notstatedonly = 1 if white5count ==. & sa5count ==. & black5count ==. & other5count ==. & mixed5count ==. & notstated5count !=. 

**MAKE CONSTANT WITHIN EACH PATIENT : no changes made
sort patid count
replace notstatedonly = notstatedonly[_n - 1] if notstatedonly[_n] == 0 & notstatedonly[_n - 1] != 0 & patid[_n] == patid[_n - 1] & totalobs > 1
gsort patid -count
replace notstatedonly = notstatedonly[_n - 1] if notstatedonly[_n] == 0 & notstatedonly[_n - 1] != 0 & patid[_n] == patid[_n - 1] & totalobs > 1
sort patid count

gen enter = 1 if count == 1
gen exit = 1 if count == total
tab enter exit
tab notstatedonly if enter == 1, missing
tab notstatedonly if exit == 1, missing


**MOST COMMON ETHNICITY EXCLUDING NOT STATED
egen eth5max = rowmax(white5count sa5count black5count other5count mixed5count)
tab eth5max


egen eth16max = rowmax(british16count irish16count otherwhite16count whitecarib16count whiteaf16count whiteasian16count othermixed16count indian16count pak16count bangla16count otherasian16count carib16count african16count otherblack16count chinese16count other16count) 
tab eth16max

**ETH 5
gen mostcommoneth5 = eth5 if eth5max == totalobs // give most common eth a value if a person only had 1 ethnicity recorded in all observations
replace mostcommoneth5 = eth5 if totalobs == 1   //makes mostcomoneth==eth if a person only has 1 observation
replace mostcommoneth5 = 0 if eth5max == white5count & eth5max !=.
replace mostcommoneth5 = 1 if eth5max == sa5count & eth5max !=.
replace mostcommoneth5 = 2 if eth5max == black5count & eth5max !=.
replace mostcommoneth5 = 3 if eth5max == other5count & eth5max !=.
replace mostcommoneth5 = 4 if eth5max == mixed5count & eth5max !=.
replace mostcommoneth5 = 5 if notstatedonly == 1
label values mostcommoneth5 eth5

tab mostcommoneth5 if enter == 1
tab mostcommoneth5 if exit == 1

**PEOPLE WITH 2 ETHNICITIES THAT ARE EQUALLY MOST COMMON- exclude not stated
**this creates a dummy which is equal to 1 if 2 ethnicities are equally common, but only 1 has been coded as being "mostcommoneth5"
gen equallycommon5 = 0
replace equallycommon5 = 1 if eth5max == white5count & mostcommoneth != 0 & totalobs != 1 & notstatedonly == 0
replace equallycommon5 = 1 if eth5max == sa5count & mostcommoneth != 1 & totalobs != 1 & notstatedonly == 0
replace equallycommon5 = 1 if eth5max == black5count & mostcommoneth != 2 & totalobs != 1 & notstatedonly == 0
replace equallycommon5 = 1 if eth5max == other5count & mostcommoneth != 3 & totalobs != 1 & notstatedonly == 0
replace equallycommon5 = 1 if eth5max == mixed5count & mostcommoneth != 4 & totalobs != 1 & notstatedonly == 0
tab equallycommon5


**update mostcommoneth to separate those with equally common ethnicities
replace mostcommoneth5 = 6 if equallycommon5 == 1 & notstatedonly == 0


**ETH 16
gen  mostcommoneth16 = eth16 if eth16max == totalobs // give most common eth a value if a person only had 1 ethnicity recorded in all observations
replace mostcommoneth16 = eth16 if totalobs == 1   //makes mostcomoneth==eth if a person only has 1 observation
replace mostcommoneth16 = 1 if eth16max == british16count & eth16max !=.
replace mostcommoneth16 = 2 if eth16max == irish16count & eth16max !=.
replace mostcommoneth16 = 3 if eth16max == otherwhite16count & eth16max !=.
replace mostcommoneth16 = 4 if eth16max == whitecarib16count & eth16max !=.
replace mostcommoneth16 = 5 if eth16max == whiteaf16count & eth16max !=.
replace mostcommoneth16 = 6 if eth16max == whiteasian16count & eth16max !=.
replace mostcommoneth16 = 7 if eth16max == othermixed16count & eth16max !=.
replace mostcommoneth16 = 8 if eth16max == indian16count & eth16max !=.
replace mostcommoneth16 = 9 if eth16max == pak16count & eth16max !=.
replace mostcommoneth16 = 10 if eth16max == bangla16count & eth16max !=.
replace mostcommoneth16 = 11 if eth16max == otherasian16count & eth16max !=.
replace mostcommoneth16 = 12 if eth16max == carib16count & eth16max !=.
replace mostcommoneth16 = 13 if eth16max == african16count & eth16max !=.
replace mostcommoneth16 = 14 if eth16max == otherblack16count & eth16max !=.
replace mostcommoneth16 = 15 if eth16max == chinese16count & eth16max !=.
replace mostcommoneth16 = 16 if eth16max == other16count & eth16max !=.
replace mostcommoneth16 = 17 if notstatedonly == 1
label values mostcommoneth16 eth16

tab mostcommoneth16 if enter == 1
tab mostcommoneth16 if exit == 1


**PEOPLE WITH 2 ETHNICITIES THAT ARE EQUALLY MOST COMMON- exclude not stated
**this creates a dummy which is equal to 1 if 2 ethnicities are equally common, but only 1 has been coded as being "mostcommoneth16"
gen equallycommon16 = 0
replace equallycommon16 = 1 if eth16max == british16count & mostcommoneth16 != 1 & totalobs != 1 & notstatedonly == 0
replace equallycommon16 = 1 if eth16max == irish16count & mostcommoneth16 != 2 & totalobs != 1 & notstatedonly == 0
replace equallycommon16 = 1 if eth16max == otherwhite16count & mostcommoneth16 != 3 & totalobs != 1 & notstatedonly == 0
replace equallycommon16 = 1 if eth16max == whitecarib16count & mostcommoneth16 != 4 & totalobs != 1 & notstatedonly == 0
replace equallycommon16 = 1 if eth16max == whiteaf16count & mostcommoneth16 != 5 & totalobs != 1 & notstatedonly == 0
replace equallycommon16 = 1 if eth16max == whiteasian16count & mostcommoneth16 != 6 & totalobs != 1 & notstatedonly == 0
replace equallycommon16 = 1 if eth16max == othermixed16count & mostcommoneth16 != 7 & totalobs != 1 & notstatedonly == 0
replace equallycommon16 = 1 if eth16max == indian16count & mostcommoneth16 != 8 & totalobs != 1 & notstatedonly == 0
replace equallycommon16 = 1 if eth16max == pak16count & mostcommoneth16 != 9 & totalobs != 1 & notstatedonly == 0
replace equallycommon16 = 1 if eth16max == bangla16count & mostcommoneth16 != 10 & totalobs != 1 & notstatedonly == 0
replace equallycommon16 = 1 if eth16max == otherasian16count & mostcommoneth16 != 11 & totalobs != 1 & notstatedonly == 0
replace equallycommon16 = 1 if eth16max == carib16count & mostcommoneth16 != 12 & totalobs != 1 & notstatedonly == 0
replace equallycommon16 = 1 if eth16max == african16count & mostcommoneth16 != 13 & totalobs != 1 & notstatedonly == 0
replace equallycommon16 = 1 if eth16max == otherblack16count & mostcommoneth16 != 14 & totalobs != 1 & notstatedonly == 0
replace equallycommon16 = 1 if eth16max == chinese16count & mostcommoneth16 != 15 & totalobs != 1 & notstatedonly == 0
replace equallycommon16 = 1 if eth16max == other16count & mostcommoneth16 != 16 & totalobs != 1 & notstatedonly == 0
tab equallycommon16

**update mostcommoneth to separate those with equally common ethnicities
replace mostcommoneth16 = 18 if equallycommon16 == 1 & notstatedonly == 0


label define eth5 6"equally common", modify
label define eth16 18"equally common", modify
compress

**ETHNICITY
**patients with valid ethnicity ever recorded
gen anyethever = 0
replace anyethever = 1 if medcodeid != ""
tab anyethever //100% as these are only patients with ethncity codes-need to make this variable 0 when attached to denominator population


**patients with valid ethnicity ever recorded
gen validethever = 0
replace validethever = 1 if eth16 != 17
replace validethever = 0 if medcodeid == "" //**0 changes
tab validethever


*make ethever constant within patients : 0 changes
sort patid count
replace anyethever=anyethever[_n - 1] if anyethever[_n] == 0 & anyethever[_n - 1] == 1 & patid[_n] == patid[_n - 1]
replace validethever=validethever[_n - 1] if validethever[_n] == 0 & validethever[_n - 1] == 1 & patid[_n] == patid[_n - 1]
gsort patid -count
replace anyethever = anyethever[_n - 1] if anyethever[_n] == 0 & anyethever[_n - 1] == 1 & patid[_n] == patid[_n - 1]
replace validethever = validethever[_n - 1] if validethever[_n] == 0 & validethever[_n - 1] == 1 & patid[_n] == patid[_n - 1]
sort patid count
replace anyethever = 0 if anyethever ==.
replace validethever = 0 if validethever ==.


**count of valid ethnicities recorded
gen validethcount = 1 if eth16 != 17 & eth16 !=.
replace validethcount = 0 if eth16 == 17 | eth16 ==.

sort patid count
replace validethcount = validethcount[_n] + validethcount[_n - 1] if patid[_n] == patid[_n - 1]

**total number of ethnicities recorded (including multiple recordings of the same ethnicity)
gen totalvalideth = validethcount if exit == 1
tab totalvalideth

*make totaleth constant for each patient
sort patid count
replace totalvalideth = totalvalideth[_n - 1] if totalvalideth[_n] ==. & totalvalideth[_n - 1] !=. & patid[_n] == patid[_n - 1]
gsort patid -count
replace totalvalideth = totalvalideth[_n - 1] if totalvalideth[_n] ==. & totalvalideth[_n - 1] !=. & patid[_n] == patid[_n - 1]
sort patid count

*dummy for multiple ethnicities excluding not stated
gen morethanoneeth = 0 if totalvalideth <= 1
replace morethanoneeth = 1 if totalvalideth > 1
tab morethanoneeth

**GEN VARIABLE FOR THE FIRST YEAR ETH WAS RECORDED
sort patid count
gen firstyear = year(eventdate) if count == 1
destring firstyear, replace
sum firstyear //1905-2019

replace firstyear = firstyear[_n - 1] if patid[_n] == patid[_n - 1] & [_n] != 1 

*are ethnicities matching
*do not give unknown ethnicity a unique counter
sort patid eth16 count
gen uniqueeth = 0
replace uniqueeth = 1 if eth16[_n] != eth16[_n - 1] & patid[_n] == patid[_n - 1] & eth16[_n] != 17 & eth16[_n - 1] != 17
replace uniqueeth = 1 if patid[_n] > patid[_n - 1]
replace uniqueeth = 0 if eth16 == 17
tab uniqueeth

*count unique ethnicities 
**EXCLUDE UNKNOWN
sort patid count
replace uniqueeth = uniqueeth[_n] + uniqueeth[_n - 1] if patid[_n] == patid[_n - 1] & uniqueeth[_n] !=. & uniqueeth[_n - 1] !=. & count != 1 
gsort patid -count
sum uniqueeth //people have up to a maximum of 9 different ethnicities

sort patid eth16 uniqueeth
replace uniqueeth=uniqueeth[_n - 1] if uniqueeth[_n]==. & uniqueeth[_n - 1]!=. & patid[_n]==patid[_n - 1] //**0 changes

*ethsum gives the number of different ethnic groups recorded
*totaluniqueeth gives the number of ethnicities recorded per patient- excluding duplicates

*count of unique ethnicities per patient
sort patid uniqueeth
bysort patid: gen totaluniqueeth = uniqueeth[_N] 
tab totaluniqueeth if enter == 1
tab totaluniqueeth if exit == 1

**dummy for yes no to having multiple unique ethnicities
**UNKNOWN Ethnicity is excluded from all counts
gen sameeth = 1 if totaluniqueeth == 1
replace sameeth = 0 if totaluniqueeth > 1
tab sameeth if enter == 1
tab sameeth if exit == 1

**indicator for whether all of the ethnicities fall under the same high level group
**ethnicity which is unknown is ignored
*ie/ if a person has 3 ethnicities, white, british and unknown, then they are considered to have matching eth5
sort patid count
gen eth5same = 1 if enter == 1
replace eth5same = 1 if eth5[_n] == eth5[_n - 1] & patid[_n] == patid[_n - 1]
replace eth5same = 1 if eth5 == 5
tab eth5same, missing

*if any eth5same values are missing- then replace as 0 and make constant
replace eth5same = 0 if eth5same ==.
replace eth5same = eth5same[_n - 1] if eth5same[_n - 1] == 0 & patid[_n] == patid[_n - 1]
gsort patid -count
replace eth5same = eth5same[_n - 1] if eth5same[_n - 1] == 0 & patid[_n] == patid[_n - 1]
sort patid count

tab eth5same if enter == 1
tab eth5same if exit == 1


sort patid count
gen eth16same = 1 if enter == 1
replace eth16same = 1 if eth16[_n] == eth16[_n - 1] & patid[_n] == patid[_n - 1]
replace eth16same = 1 if eth16 == 17
tab eth16same, missing

*if any eth5same values are missing- then replace as 0 and make constant
replace eth16same = 0 if eth16same ==.
replace eth16same = eth16same[_n - 1] if eth16same[_n - 1] == 0 & patid[_n] == patid[_n - 1]
gsort patid -count
replace eth16same = eth16same[_n - 1] if eth16same[_n - 1] == 0 & patid[_n] == patid[_n - 1]
sort patid count

tab eth16same if enter == 1
tab eth16same if exit == 1


*FIXED VARIABLE FOR LATEST ETH
gsort patid -eventdate
gen latesteth = medcodeid 
gen latestdesc = term 
gen latesteth16 = eth16 
gen latesteth5 = eth5

replace latesteth =latesteth[_n - 1] if patid[_n] == patid[_n - 1] & [_n] != 1 
replace latestdesc = latestdesc[_n - 1] if patid[_n] == patid[_n - 1] & [_n] != 1 
replace latesteth16 = latesteth16[_n - 1] if patid[_n] == patid[_n - 1] & [_n] != 1 
replace latesteth5 = latesteth5[_n - 1] if patid[_n] == patid[_n - 1] & [_n] != 1 

*label values latesteth16 eth16
label values latesteth5 eth5
compress
save "${file_stub}_Observation_all_eth.dta", replace

**MAKE PATIENT LEVEL FILE
**drop duplicate patids
use "${file_stub}_Patient_included_all.dta", clear

merge 1:m patid using "${file_stub}_Observation_all_eth.dta", gen(merge2)

codebook patid //
duplicates drop patid, force
codebook patid //

keep patid latesteth16 latesteth5 mostcommon*

**GEN ONE VARIABLE FOR ETHNICITY
gen eth5 = mostcommoneth5 //main ethnicity is most common in CPRD
label values eth5 eth5
tab eth5, missing

*remove equally common group
replace eth5=latesteth5 if eth5 >= 5 & latesteth5 !=.  //replace ethnicity with latest eth5 if mostcommoneth5 is not stated/equal/missing
tab eth5, missing

gen eth16 = mostcommoneth16 
label values eth16 eth16
replace eth16 = latesteth16 if eth16 >= 17 & latesteth16 < 17  //replace ethnicity with latest eth16 if mostcommoneth5 is not stated/equal/missing
tab eth16, missing
tab eth16,m

label define eth16 17"Unknown", modify

compress

keep patid eth5 eth16  
compress 
save "${file_stub}_Patient_ethnicity", replace
log close

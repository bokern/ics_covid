
/*=========================================================================
DO FILE NAME:		    07cr_copd_product.do

AUTHOR:					Marleen Bokern

VERSION:				v1

DATE VERSION CREATED: 	12/2022

DATASETS CREATED:       copd_DrugIssue_ics_single
						copd_DrugIssue_ics_laba
						copd_DrugIssue_laba_lama
						copd_DrugIssue_laba_single
						copd_DrugIssue_lama_single
						copd_DrugIssue_triple_therapy
						
DESCRIPTION OF FILE:	Generate list of prescriptions for each drug type, merged with patient list with excl/incl criteria applied and dosage lookup file

*=========================================================================*/
/**************************************************************************************************************
HOUSEKEEPING
***************************************************************************************************************/

clear all

capture log close 
log using $Logdir\07cr_copd_product.log, replace

****get dosages lookup file as dta
cd "$Datadir_copd"

import delimited "$Denominator\CPRD Aurum\Lookups\2022_05\common_dosages.txt", clear
save "$Mainfolder\common_dosages", replace

glob common_dosages "$Mainfolder\common_dosages.dta"

/*******************************************************************************************************************
merge each drug type drug issue file with included patient list and dosage lookup file
*******************************************************************************************************************/

local drugissue ics_single ics_laba laba_lama laba_single lama_single triple_therapy budesonide_single budesonide_laba 
local drugissue triple_budesonide
disp `"`drugissue'"'

foreach drug in `drugissue' {
	noi di "File `drug'"
    use "${file_stub}_DrugIssue_`drug'.dta", clear
	merge m:1 patid using "${file_stub}_Patient_included_all", keep(match) nogen
	merge m:1 dosageid using "$common_dosages", keep(match master) nogen
	drop dose_duration dose_max_average change_dose dose_unit dose_interval dose_frequency choice_of_dose dose_number enterdate day yob mob yo35bday do35bday lcd copd_date
	cap drop drugissues
	compress
	save "${file_stub}_`drug'.dta", replace
}

log close
clear all
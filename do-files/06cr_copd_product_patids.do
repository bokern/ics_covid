/*=========================================================================
DO FILE NAME:		    06cr_copd_product_patids.do

AUTHOR:					Marleen Bokern

VERSION:				v1

DATE VERSION CREATED: 	06/2023

DATASETS CREATED:       copd_DrugIssue_ics_single
						copd_DrugIssue_ics_laba
						copd_DrugIssue_laba_lama
						copd_DrugIssue_laba_single
						copd_DrugIssue_lama_single
						copd_DrugIssue_triple_therapy
						
DESCRIPTION OF FILE:	Generates list of drug issues for each drug type. 

*=========================================================================*/
/************************************************************************
HOUSEKEEPING
************************************************************************/

clear all
capture log close 
log using $Logdir\06cr_copd_product_patids.log, replace

cd "$Datadir_copd"

// Specify file names
glob file_stub 			= 	"copd"
glob file_Patient 		= 	"${file_stub}_Extract_Patient_"
glob file_Practice 		=	"${file_stub}_Extract_Practice_"
glob file_Staff			=	"${file_stub}_Extract_Staff_"
glob file_Consultation 	= 	"${file_stub}_Extract_Consultation_"
glob file_Observation	= 	"${file_stub}_Extract_Observation_"
glob file_Referral 		= 	"${file_stub}_Extract_Referral_"
glob file_Problem 		= 	"${file_stub}_Extract_Problem_"
glob file_DrugIssue		= 	"${file_stub}_Extract_DrugIssue_"

// Specify number of different files
glob no_Patient = 1
glob no_Practice = 1
glob no_Staff = 1
glob no_Consultation = 12
glob no_Observation = 48
glob no_Referral = 1
glob no_Problem = 1
glob no_DrugIssue = 49

/***************************************************************************************
MANAGE CODELISTS - use if necessary
****************************************************************************************/

/***check before that codelists are unique
cd $Codelistsdir\Product_codelists

local products : dir "" files "*.dta"
disp `"`products'"'


foreach x in `products'{
	use `x'
	duplicates report
	duplicates drop
	save `x', replace
}
clear all

****************/

cd $Datadir_copd

/*****************************************************************************************
Inhalers - 
TRIPLE THERAPY has been extracted in inclusion/ exclusion criteria file
******************************************************************************************/

local product ics_single ics_laba laba_lama laba_single lama_single
local product budesonide_single budesonide_laba 
local product triple_budesonide
disp `"`product'"'

foreach drug of local product {
	foreach file of numlist 1/$no_DrugIssue {
		noi di "Merging `drug', File `file'"
		use "$Copd_aurum_extract\\${file_stub}_Extract_DrugIssue_`file'", clear
		drop if issuedate > td(30apr2021)
		drop pracid estnhscost
		merge m:1 prodcodeid using "$Codelistsdir\Product_codelists\cl_`drug'.dta", keep(match) nogen
		cap drop drugissues
		if `file' == 1{
			save "${file_stub}_DrugIssue_`drug'.dta", replace
		}
		if `file' > 1{
			append using "${file_stub}_DrugIssue_`drug'.dta"
		}
		compress
		save "${file_stub}_DrugIssue_`drug'.dta", replace
	}
}

/*****************************************************************************************************
CHEMOTHERAPY - IMMUNOSUPPRESSION COVARIATE
******************************************************************************************************/
clear all

local product immunosuppressants
disp `"`product'"'

import delim "$Codelistsdir/Comorbidities/immunosuppressants_aurum_mar20.csv", stringcols(1)
collapse (max) biological cancer_chemo steroidinj methotrexate azathioprine mercaptopurine other_dmard other_immunosupp oralsteroid, by(prodcodeid)
save "$Codelistsdir/Comorbidities/cl_immunosuppressants", replace

foreach file of numlist 1/$no_DrugIssue {
	noi di "Merging immunosuppressants Drug Issue, File `file'"
    use "$Copd_aurum_extract\\${file_stub}_Extract_DrugIssue_`file'", clear
	drop if issuedate > td(30apr2021)
	drop pracid estnhscost
    merge m:1 prodcodeid using "$Codelistsdir/Comorbidities/cl_`product'", keep(match) nogen
	cap drop drugissues
	if `file' == 1{
		save "${file_stub}_DrugIssue_`product'.dta", replace
	}
	if `file' > 1{
		append using "${file_stub}_DrugIssue_`product'.dta"
	}
	compress
	save "${file_stub}_DrugIssue_`product'.dta", replace
}

clear all
log close

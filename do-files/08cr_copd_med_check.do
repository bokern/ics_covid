/*=========================================================================
DO FILE NAME:		    08cr_copd_med_check.do

AUTHOR:					Marleen Bokern

VERSION:				v1

DATE VERSION CREATED: 	01/2023

DATASETS CREATED:       copd_{drugclass}_checked.dta
						
DESCRIPTION OF FILE:	checks drug issue files for each drug class for missing drug substance names and imputes from term variable where necessary, drops unnecessary variables 
*=========================================================================*/
/****************************************************
HOUSEKEEPING
*****************************************************/
*run globals
clear all
capture log close
log using "$Logdir\08cr_copd_med_check.log", replace
cd "$Datadir_copd"
ssc install ereplace
ssc install moss

/****************************************************************************
ICS SINGLE
*****************************************************************************/

global drugclass ics_single
use "${file_stub}_${drugclass}", clear //this file is the DrugIssue file of single ICS inhalers (drug issue files merged with ICS single codelist, matches kept)

sort drugsubstancename

***replace missing drugsubstancenames using info from term variable 
count if missing(drugsubstancename) //*368

replace drugsubstancename = "ciclesonide" if regex(termfromemis, "ciclesonide") |regex(termfromemis, "alvesco") & drugsubstancename==""
replace drugsubstancename = "mometasone furoate" if regex(termfromemis, "mometasone furoate")|regex(termfromemis, "asmanex") & drugsubstancename==""
replace drugsubstancename = "fluticasone propionate" if regex(termfromemis, "fluticasone propionate")|regex(termfromemis, "flixotide") & drugsubstancename==""
replace drugsubstancename = "beclometasone dipropionate" if (regex(termfromemis, "beclomethasone")|regex(termfromemis, "becotide")|regex(termfromemis, "becloforte")|regex(termfromemis, "asmabec")|regex(termfromemis, "betamethasone")) & drugsubstancename==""
replace drugsubstancename = "budesonide" if regex(termfromemis, "budesonide") | regex(termfromemis, "pulmicort") & drugsubstancename==""
count if missing(drugsubstancename) //*0
assert drugsubstancename != ""

compress
save "${file_stub}_${drugclass}_checked", replace
clear all

/****************************************************************************
ICS LABA
*****************************************************************************/
global drugclass ics_laba
use "${file_stub}_${drugclass}", clear

sort drugsubstancename
count if missing(drugsubstancename) //*0
assert drugsubstancename != ""
***no missing drugsubstancenames

compress
save "${file_stub}_${drugclass}_checked", replace
clear all

/****************************************************************************
LABA LAMA
*****************************************************************************/
global drugclass laba_lama
use "${file_stub}_${drugclass}", clear

sort drugsubstancename
count if missing(drugsubstancename) //*0
assert drugsubstancename != ""
***no missing drugsubstancenames

compress
save "${file_stub}_${drugclass}_checked", replace
clear all

/****************************************************************************
LABA SINGLE
*****************************************************************************/
global drugclass laba_single
use "${file_stub}_${drugclass}", clear

sort drugsubstancename
count if missing(drugsubstancename) //*0
replace drugsubstancename = "formoterol" if regex(termfromemis, "formoterol")  & drugsubstancename == ""
replace drugsubstancename = "salmeterol xinafoate" if regex(termfromemis, "salmeterol xinafoate")  & drugsubstancename == ""

assert drugsubstancename !=""

compress
save "${file_stub}_${drugclass}_checked", replace
clear all

/****************************************************************************
LAMA SINGLE
*****************************************************************************/
global drugclass lama_single
use "${file_stub}_${drugclass}", clear

sort drugsubstancename
count if missing(drugsubstancename) //*0
***no missing drugsubstancenames

assert drugsubstancename != ""
compress
save "${file_stub}_${drugclass}_checked", replace
clear all

/****************************************************************************
TRIPLE THERAPY
*****************************************************************************/
global drugclass triple_therapy
use "${file_stub}_${drugclass}", clear

sort drugsubstancename
count if missing(drugsubstancename) //*0
***no missing drugsubstancenames

assert drugsubstancename != ""
compress
save "${file_stub}_${drugclass}_checked", replace

clear all
/****************************************************************************
BUDESONIDE SINGLE
*****************************************************************************/
global drugclass budesonide_single
use "${file_stub}_${drugclass}", clear


sort drugsubstancename
count if missing(drugsubstancename)
replace drugsubstancename = "budesonide" if regex(termfromemis, "budesonide") | regex(termfromemis, "pulmicort") & drugsubstancename==""
count if missing(drugsubstancename)

***no missing drugsubstancenames

assert drugsubstancename != ""
compress
save "${file_stub}_${drugclass}_checked", replace

clear all
/****************************************************************************
BUDESONIDE LABA
*****************************************************************************/
global drugclass budesonide_laba
use "${file_stub}_${drugclass}", clear

sort drugsubstancename
count if missing(drugsubstancename)
replace drugsubstancename = "budesonide" if regex(termfromemis, "budesonide") | regex(termfromemis, "pulmicort") & drugsubstancename==""
count if missing(drugsubstancename)

***no missing drugsubstancenames

assert drugsubstancename != ""
compress
save "${file_stub}_${drugclass}_checked", replace

clear all
/****************************************************************************
BUDESONIDE TRIPLE
*****************************************************************************/
global drugclass triple_budesonide
use "${file_stub}_${drugclass}", clear

sort drugsubstancename
count if missing(drugsubstancename)
*replace drugsubstancename = "budesonide" if regex(termfromemis, "budesonide") | regex(termfromemis, "pulmicort") & drugsubstancename==""
count if missing(drugsubstancename)

***no missing drugsubstancenames

assert drugsubstancename != ""
compress
save "${file_stub}_${drugclass}_checked", replace

clear all
log close 

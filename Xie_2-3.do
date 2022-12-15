/*******************************************************************************
	Replication paper 2: Auffhammer and Kellogg (2011)
	Clearing the Air? The Effects of Gasoline Content Regulation on Air Quality.
	
	5 June 2022
	Jack (Quan Cheng) Xie
	Topics in Econometrics
	Uppsala University--Masters in Economics
	
*******************************************************************************/

cls
set more off

// Dependencies for attachments
* ssc install outreg2
* ssc install tabout
* ssc install estout
* ssc install texsave

*** Preamble ***

// Base directory containing replication package
// https://www.openicpsr.org/openicpsr/project/112465/version/V1/view
global base "C:\Users\qcx20\Documents\School\4-Topics_in_Econometrics\Replication"

global support "${base}\support"
global attachments "${base}\attachments"
capture mkdir "$support"
capture mkdir "$attachments"

cd "${base}\112465-V1\AER20090377_DataAndPrograms"


capture log close
log using "${support}\Xie_2.log", replace


/*******************************************************************************
	Part A1. Difference-in-differences Data
	Cleaning according to section II
*******************************************************************************/

use AER20090377_FinalData, clear

*** Drop monitor-days where nine hours from 9:00-21:00 unobserved ***

// Removal (i): drop if less than 9 hours of observation between 9 AM and 9 PM
drop if valid < 9

// only data between 1989 and 2003
keep if (1989 <= year) & (year <= 2003)


*** Drop monitor-years where 25% unobserved summer days***

gen Season = 1

// winter 12-2, spring 3-5, summer 6-8, fall 9-11
forvalues i = 1/3 {
	quiet replace Season = `i'+1 if (`i'*3 <= month) & (month < `i'*3 + 3)
}

tab month Season

// Removal (ii): keep if 75% days observed per Season
// group Dec with next year
gen year2 = year
replace year2 = year + 1 if month==12

// drop by Season (if less than 75% of summer days)
bysort panelid year2 Season: egen SeasonDays = count(day)

gen temp = .
replace temp = SeasonDays if Season == 3
drop SeasonDays year2

bysort panelid year: egen summer_days = max(temp)
drop temp

keep if summer_days >= 0.75 * (30 + 31 + 31) & summer_days !=.
drop summer_days

tab year

*** Weather shocks W_{it} ***

gen DOW = dow(Date)
gen DOY = doy(Date)

// time weather interactions
xi i.DOW*TempMax, prefix(_DM)
drop _DMDOW_*

xi i.DOW*TempMin, prefix(_Dm)
drop _DmDOW_*

xi i.DOW*Rain, prefix(_Dr)
drop _DrDOW_*

xi i.DOW*Snow, prefix(_Ds)
drop _DsDOW_*


// temperature function terms
foreach var of varlist(Temp* Rain Snow){
	forvalues i = 1/3 {
		gen `var'`i' = `var'^`i'
	}
}

// Rain and Snow cubics not used
drop Rain3 Snow3

// interactions
gen TempMaxMin = TempMax * TempMin
gen RainTempMax = Rain * TempMax

// lagged temperature
sort panelid Date
by panelid: gen tracker = _n

foreach var of varlist(TempMax TempMin){
	gen `var'L1 = TempMax[_n-1]
	replace `var'L1 = . if tracker==1
}
gen TempMaxMaxL1 = TempMax * TempMaxL1
gen TempMaxMinL1 = TempMax * TempMinL1
drop tracker

// time temp interactions
foreach var of varlist(TempMax1-TempMaxMinL1) {
	gen DOY`var' = DOY * `var'
}


*** Region and treatment trends Trend_{rct} ***

local region1 9 23 25 33 34 36 42 44 50
local region2 17 18 19 20 26 27 29 31 38 39 46 55
local region3 1 5 10 11 12 13 21 22 24 28 37 40 45 47 48 51 54
local region4 2 4 6 8 15 16 30 32 35 41 49 53 56

gen region=.
forvalues i = 1/4 {
	foreach n of local region`i' {
		replace region = `i' if state_code == `n'
	}
}

// create interactions of region with year, DOW, and DOY
xi i.year*i.region, prefix(_RY)
xi i.DOW*i.region, prefix(_RW)
xi i.region*DOY, prefix(_RD)


// turn RFG off if CARB
replace treat_rfg = 0 if treat_CARB~=0

// turn off RVP in non-summer months
replace treat_rvpI = 0 if Season != 3
replace treat_rvpII = 0 if Season != 3

// recreate county treatments as indicators for if county ever treated
drop RVPCty RFGCty CARBCty

bysort fips: egen RVP_Cty = max(treat_rvpII)
bysort fips: egen RFG_Cty = max(treat_rfg)
bysort fips: egen CARB_Cty = max(treat_CARB)

// indicators for if treated both RVP and RFG
gen RVP_RFG_Cty = 0
gen CARB_RFG_Cty = 0

// variables for multiple treatments
replace CARB_RFG_Cty = 1 if RFG_Cty==1 & CARB_Cty==1
replace RVP_RFG_Cty = 1 if RFG_Cty==1 & RVP_Cty==1
replace CARB_Cty = 0 if CARB_RFG_Cty==1

// reset individual treatment to 0 if multiple
replace RFG_Cty = 0 if (RVP_RFG_Cty==1  | CARB_RFG_Cty==1)
replace RVP_Cty = 0 if RVP_RFG_Cty==1

// CARB > RVP
replace RVP_RFG_Cty = 0 if CARB_RFG_Cty==1


*** Add daily time trend ***

// standardize date to annual increments
gen DateS = Date / 365	
gen DateS2 = DateS^2

// treatment types (including combinations)
local treats RVP RFG RVP_RFG CARB CARB_RFG

// region-type treatment time trends
foreach treat of local treats {
	forvalues i = 1/4 {
		
		// 1st and 2nd order treatment type-region time trends
		gen Trend`treat'`i' = 0
		gen QTrend`treat'`i' = 0
		
		replace Trend`treat'`i' = DateS if `treat'_Cty==1 & region==`i'
		replace QTrend`treat'`i' = DateS2 if `treat'_Cty==1 & region==`i'
	}
}


*** More data trimming ***

// Removal (iii): Drop if no weather or income variables

// drop if zeros for ozone
drop if ozone_max==0 | epa_8hr==0

// drop if no weather
foreach var of varlist(TempMax1-TempMaxMinL1) {
	drop if `var' == .
}

// drop if no income
merge m:1 state_code county_code year using AER20090377_IncomeData.dta, ///
	generate(_merge_income)
	
keep if _merge_income==3

// drop if no income
drop if income == .

// rescale to billions
replace income = income / 1000000000


// drop variables not needed
drop state-regtype
drop sulfur-rfgtype
drop _merge*


// take logs of daily ozone
gen lozone_max = log(ozone_max)
gen lozone_8hr = log(epa_8hr)

rename epa_8hr ozone_8hr

// cluster by state-year
sort state_code year
egen StateYear = group(state_code year)

tab year

save "${support}\DD_data.dta", replace

/*******************************************************************************
	Part A2. Difference-in-differences results
	
	Replicate results with errors and robustness checks
	
	Robustness specification:
	* Neighbor merge
	* Do not mean difference log ozone and treatment
	* Add constant to regression
	
*******************************************************************************/

use "${support}\DD_data.dta", replace

*** Replicate with merge error ***

// Removal (iv): drop if neighbors treated with more stringent regulation
// this was an error in author's code
merge m:1 fips using "AER20090377_NeighborData.dta", ///
	generate(_merge_stricter_neighbor)

keep if _merge_stricter_neighbor == 1
// drop if treated_neighbor == 1 | _merge_stricter_neighbor == 2
drop treated_neighbor _merge_stricter_neighbor

// Removal (v): keep summer only
keep if Season == 3


// set panel data--allow monitor FEs
xtset panelid Date, daily

// mean-differenced variables (mean diff treatment?)
local mean_diff lozone_* treat* Temp* Rain* Snow*  DOY* Trend* QTrend* DateS  _D* _R* income

// non-mean-differenced variables
local non_mean_diff fips month year StateYear state_code urban

// keep relevant variables
keep `mean_diff' `non_mean_diff' panelid

// take panel mean differences
sort panelid

// replication without mean-differenced
save "${support}\DD_repnd_full.dta", replace

foreach var of varlist(`mean_diff') {
	
	by panelid: egen mean_`var' = mean(`var')
	
	gen `var'_diff = `var' - mean_`var'
	
	drop `var' mean_`var'
	
	rename `var'_diff `var'
}

// replication data
save "${support}\DD_rep_full.dta", replace

*** Robustness check data ***

use "${support}\DD_data.dta", replace

// drop if neighbors treated with more stringent regulation
// this was an error in author's code
merge m:1 fips using "AER20090377_NeighborData.dta", ///
	generate(_merge_stricter_neighbor)
	
// keep if _merge_stricter_neighbor == 1
drop if treated_neighbor == 1 | _merge_stricter_neighbor == 2
drop treated_neighbor _merge_stricter_neighbor
	
// keep summer only
keep if Season == 3

// set panel data--allow monitor FEs
xtset panelid Date, daily

// mean-differenced variables (mean diff treatment?)
local mean_diff Temp* Rain* Snow* DOY* Trend* QTrend* DateS _D* _R* income

// non-mean-differenced variables
local non_mean_diff lozone_* treat* fips month year StateYear state_code urban

// keep relevant variables
keep `mean_diff' `non_mean_diff' panelid

// take panel mean differences
sort panelid

// correct merge data
save "${support}\DD_cmerge_full.dta", replace


*** Save restricted sample for monitors observed through all years ***

local datasets rep repnd cmerge

foreach data of local datasets {
	use "${support}\DD_`data'_full.dta", clear	
		
	bysort panelid year: gen NumYear = _n == 1
	by panelid: replace NumYear = sum(NumYear)
	by panelid: replace NumYear = NumYear[_N]
	tab NumYear

	keep if NumYear == 15

	save "${support}\DD_`data'_15years.dta", replace
}


*** Table 1 data summary ***

// tabout documentation: http://fmwww.bc.edu/RePEc/bocode/t/tabout.html
// https://www.stata.com/meeting/oceania16/slides/watson-oceania16.pdf
* tabout year using `basedir'\test.tex, replace style(tex)

// https://stackoverflow.com/questions/34296678/stata-esttab-twoway-tabulate-labels

est clear

use "${support}\DD_rep_full.dta", clear

// observations per year
eststo: estpost tab year

// counties (fips) per year
keep year fips
duplicates drop

eststo: estpost tab year

// monitors per year
use "${support}\DD_rep_full.dta", clear
keep year panelid
duplicates drop

eststo: estpost tab year


// urban and rural counties per year
use "${support}\DD_rep_full.dta", clear
keep year panelid urban
duplicates drop

eststo: estpost tab year if urban == 1

eststo: estpost tab year if urban == 3

// treated monitors per year
use "${support}\DD_rep_full.dta", clear

collapse(max) treat*, by(state_code fips year)

foreach var of varlist(treat_rvpI treat_rvpII treat_rfg treat_CARB) {
	
	// treated (mean is zero)
	eststo: estpost tab year if `var' > 0
}

noisily esttab using "${attachments}\table1.tex", replace ///
	unstack cell(b(fmt(%15.0fc))) noobs label booktabs f ///
    mtitles("Observations" "Counties" "Total monitors" "Urban" "Rural" ///
	"RVPI\textsuperscript{a}" "RVPII\textsuperscript{b}" "RFG95" "CARB") collabels(none) sfmt(%9.0fc) ///
	nonumbers mgroups("Daily" "Counts of active monitors" "Counts of treatments", ///
	pattern(1 0 1 0 0 1 0 0 0) prefix(\multicolumn{@span}{c}{) suffix(}) ///
	span erepeat(\cmidrule(lr){@span}))
	


*** Run DD Regressions ***

// treatment and full control variables
local treats treat_rvpI treat_rvpII treat_rfg treat_CARB
local controls income Temp* Rain* Snow*  DOY* _D* _R* Trend* QTrend*

// samples and datasets
local samples full 15years
local datasets rep repnd cmerge

// regression options
local reg_options noconstant vce(cluster StateYear)

est clear

quietly {

	foreach data of local datasets {
		
		if "`data'" == "cmerge" {
			local CMv "Yes"
		}
		else {
			local CMv "No"
		}
		
		if "`data'" == "rep" {
			local MDv "Yes"
		}
		else {
			local MDv "No"
		}
		
		
		foreach sample of local samples {
			
			noisily disp "`data' `sample' sample regressions"
			use "${support}\DD_`data'_`sample'.dta", clear
			
			// monitor
			eststo: regress lozone_max `treats' _RY*, `reg_options'
			estadd local  MFE  "Yes"
			estadd local  RYE  "Yes"
			estadd local  Ctrls  "No"
			estadd local  MD	`MDv'
			estadd local  CM	`CMv'
			
			// with time region FEs and weather controls
			eststo: regress lozone_max `treats' `controls', `reg_options'
			estadd local  MFE  "Yes"
			estadd local  RYE  "Yes"
			estadd local  Ctrls  "Yes"
			estadd local  MD	`MDv'
			estadd local  CM	`CMv'

		}
	}
}

// label variables
label var treat_rvpI "RVPI (9.5, 10.5 psi)"
label var treat_rvpII "RVPII (7.8 psi)"
label var treat_rfg  "Federal RFG"
label var treat_CARB "CARB Gasoline"
label var income "Income"

// attachments regression
noisily esttab using "${attachments}\table2.tex", replace f   ///
	keep(treat* income) ///
	nocons se(3) b(3) label booktabs nomtitle collabels(none) compress ///
	stats(MFE RYE Ctrls MD CM N, label("Monitor FEs" "Region-Year FEs" ///
	"Full controls" "Mean Difference" "Corrected Merge" "Observations") fmt(%9.0fc)) ///
	mgroups("Full sample" "All year monitors" "Full sample" "All year monitors" ///
	"Full sample" "All year monitors", ///
	pattern(1 0 1 0 1 0 1 0 1 0 1 0) prefix(\multicolumn{@span}{c}{) suffix(}) ///
	span erepeat(\cmidrule(lr){@span}))


	
/*******************************************************************************
	Part A3. Difference-in-differences results extension
*******************************************************************************/

use "${support}\DD_data.dta", clear

// drop if neighbors treated with more stringent regulation
merge m:1 fips using "AER20090377_NeighborData.dta", ///
	generate(_merge_stricter_neighbor)

// keep if _merge_stricter_neighbor == 1
drop if treated_neighbor == 1 | _merge_stricter_neighbor == 2
drop treated_neighbor _merge_stricter_neighbor

// gen week = week(Date)

sort panelid year month Date
by panelid year month: egen month_start = min(Date)
by panelid year month: egen month_end 	= max(Date)
by panelid year month: egen month_days 	= count(Date)

collapse	///
	(max) treat* month_start month_end month_days income ///
	(mean) ozone* TempMax TempMin Rain Snow, ///
	by(panelid state_code county_code fips *_Cty site_id region urban year Season month StateYear)


*** Weather shocks W_{it} ***
xi i.year*i.region, prefix(_RY)
xi i.month*i.region, prefix(_RM)

foreach var of varlist(TempMax* Rain Snow) {
	xi i.month*`var', prefix(_M)
	drop _Mmonth_*
 }


// temperature function terms
foreach var of varlist(TempMax* TempMin Rain Snow){
	forvalues i = 1/3 {
		gen `var'`i' = `var'^`i'
	}
}

// Rain and Snow cubics not used
drop Rain3 Snow3

// interactions
gen TempMaxMin = TempMax * TempMin
gen RainTempMax = Rain * TempMax

// lagged temperature
sort panelid year month
by panelid: gen tracker = _n

foreach var of varlist(TempMax TempMin){
	gen `var'L1 = TempMax[_n-1]
	replace `var'L1 = . if tracker==1
}
gen TempMaxMaxL1 = TempMax * TempMaxL1
gen TempMaxMinL1 = TempMax * TempMinL1
drop tracker

// time temp interactions
foreach var of varlist(TempMax1-TempMaxMinL1) {
	gen Month`var' = month * `var'
}


*** Add monthly time trend ***

// create time trend
sort year month
egen YearMonth = group(year month)

// annualize
gen YearMonthS = YearMonth / 12
gen YearMonthS2 = YearMonth^2

// treatment types (including combinations)
local treats RVP RFG RVP_RFG CARB CARB_RFG

// region-type treatment time trends
foreach treat of local treats {
	forvalues i = 1/4 {
		
		// 1st and 2nd order treatment type-region time trends
		gen Trend`treat'`i' = 0
		gen QTrend`treat'`i' = 0
		
		replace Trend`treat'`i' = YearMonthS if `treat'_Cty==1 & region==`i'
		replace QTrend`treat'`i' = YearMonthS2 if `treat'_Cty==1 & region==`i'
	}
}



*** staggered treatment time dummies ***

// keep only summer months
// replace treat_rvpI = 0 if Season != 3
// replace treat_rvpII = 0 if Season != 3
keep if Season == 3

sort panelid year month

// month index
by panelid: gen t_month = _n

foreach var of varlist(treat*) {

	disp "`var'"
	// treatment type
	local treat = subinstr("`var'", "treat_", "", .)
	
	disp "`treat'"
	
	// code =. for non-treated
	gen temp_`treat' = .
	replace temp_`treat' = 1 if `var' == 1
	
	// first (min) month of treatment
	by panelid: egen start`treat' = min(t_month * temp_`treat')
	drop temp_`treat'

	// treatment starts at t_`treat'==0
	gen t_`treat' =  t_month - start`treat'
	
	// sort panelid
	// by panelid: egen ever_`treat' = max(`var')
	// by panelid: egen always_`treat' = min(`var')
	
	// categorical variables
	xi i.t_`treat', prefix("_E")

	// rename categoroies
	foreach var of varlist(_E*) {

		// label name of variable category
		local labelname : variable label `var'
		
		// negative _n
		local labelname = subinstr("`labelname'", "==-", "_n", .)
		
		// positive _p
		local labelname = subinstr("`labelname'", "==", "_p", .)
		
		// rename with label
		rename `var' `labelname'
		
	}
}


// drop if zeros for ozone
drop if ozone_max==0 | ozone_8hr==0

// drop if no weather
foreach var of varlist(TempMax1-TempMaxMinL1) {
	drop if `var' == .
}


gen lozone_max = log(ozone_max)
gen lozone_8hr = log(ozone_8hr)

// set panel data
xtset panelid YearMonth

// variables to standardize (not include treat)
local mean_diff lozone_* treat* Temp* Rain* Snow* Month* Trend* QTrend* _M* _R* income

// variables to not standardized 
local non_mean_diff fips month year StateYear ///
	region state_code urban panelid t_*

// keep relevant variables
keep `mean_diff' `non_mean_diff'

// take panel mean differences
sort panelid


foreach var of varlist(`std_vars') {
	
	// mean and standard deviation
	by panelid: egen mean_`var' = mean(`var')
	
	// standardize with mean 0 and variance 1
	gen `var'_diff = `var' - mean_`var'
	
	drop `var' mean_`var'
	
	rename `var'_diff `var'
}

// save staggered regression data
save "${support}\stag_data.dta", replace


*** Run regressions ***

use "${support}\stag_data.dta", clear

// control variables
local controls income Temp* Rain* Snow* Month* _M* _R* Trend* QTrend*

// regression parameters
local reg_options nocons vce(cluster StateYear)

local windows 3 6 12

local treats rvpII rfg CARB

est clear

// coefficients
local coeffs
forvalues i = -12/12 {
	if `i' < 0 {
		local n = -`i'
		gen t_n`n' = .
		local coeffs `coeffs' t_n`n'
		label var t_n`n' "$\alpha_{\hat t=`i'}$"
	}
	else {
		gen t_p`i' = .
		local coeffs `coeffs' t_p`i'
		label var t_p`i' "$\alpha_{t=`i'}$"
	}

}

// save results

foreach treat of local treats {

	// first regression
	local first 1
	
	forvalues i = -12/12 {
		if `i' < 0 {
			local n = -`i'
			replace t_n`n' = .
			replace t_n`n' = t_`treat'_n`n'
		}
		else {
			replace t_p`i' = .
			replace t_p`i' = t_`treat'_p`i'
		}
	}
		

	foreach w of local windows {

		
		// save type
		if `first'==1 {
			local save replace
			local first 0
		}
		else {
			local save append
		}
		
		// treatment-time indicators
		local indicators t_n`w'-t_p`w'
		
		eststo: regress lozone_max `indicators'	///
			_RY*, 								///
			`reg_options'	///
			
		noisily {
			disp "`treat' window `w'"
			testparm t_n*
			local fstat = round(r(F), 0.0001)
			local pval = round(r(p), 0.0001)
		}
		
		estadd local  WD  "$\pm`w'$"
		estadd local  MFE  "Yes"
		estadd local  RYE  "Yes"
		estadd local  Ctrls  "No"
		estadd local  MD  "Yes"
		estadd local  CM  "Yes"
		estadd local  Fstat "`fstat'"
		estadd local  Pval 	"`pval'"
	
		// save control results
		gen Beta_`treat'_`w' = .
		gen SE_`treat'_`w' = .
		
		replace Beta_`treat'_`w' = _b[t_p0] if t_`treat' == 0
		replace SE_`treat'_`w' = _se[t_p0] if t_`treat' == 0
		
		forvalues i = 1/`w' {
			
			// positive lags
			local pvar t_p`i'
			replace Beta_`treat'_`w' = _b[`pvar'] if t_`treat' == `i'
			replace SE_`treat'_`w' = _se[`pvar'] if t_`treat' == `i'
			
			// negative leads
			local nvar t_n`i'
			replace Beta_`treat'_`w' = _b[`nvar'] if t_`treat' == -`i'
			replace SE_`treat'_`w' = _se[`nvar'] if t_`treat' == -`i'
		}
		
		eststo: regress lozone_max `indicators'	///
			`controls',						 	///
			`reg_options'
		
		// test lead coefficients are jointly zero
		// test lead (negative t) coefficients are jointly zero
		noisily {
			disp "`treat' window `w'"
			testparm t_n*
			local fstat = round(r(F), 0.0001)
			local pval = round(r(p), 0.0001)
		}
		
		estadd local  WD  "$\pm`w'$"
		estadd local  MFE  "Yes"
		estadd local  RYE  "Yes"
		estadd local  Ctrls  "Yes"
		estadd local  MD  "Yes"
		estadd local  CM  "Yes"
		estadd local  Fstat "`fstat'"
		estadd local  Pval 	"`pval'"
			
		// save control results
		gen Beta_`treat'_`w'C = .
		gen SE_`treat'_`w'C = .
		
		replace Beta_`treat'_`w'C = _b[t_p0] if t_`treat' == 0
		replace SE_`treat'_`w'C = _se[t_p0] if t_`treat' == 0
		
		forvalues i = 1/`w' {
			
			local pvar t_p`i'
			local nvar t_n`i'
			
			// positive lags
			replace Beta_`treat'_`w'C = _b[`pvar'] if t_`treat' == `i'
			replace SE_`treat'_`w'C = _se[`pvar'] if t_`treat' == `i'
			
			// negative leads
			replace Beta_`treat'_`w'C = _b[`nvar'] if t_`treat' == -`i'
			replace SE_`treat'_`w'C = _se[`nvar'] if t_`treat' == -`i'
			
		}
	}

}

label var income `County income'

noisily esttab using "${attachments}\tableA1.tex", replace f   ///
	keep(`coeffs' income) ///
	nocons se(3) b(3) label booktabs nomtitle collabels(none) compress ///
	stats(WD MFE RYE Ctrls MD CM Fstat Pval N, label("Window" "Monitor FEs" ///
	"Region-Year FEs" "Full controls" "Mean differenced" "Corrected merge" ///
	"F-stat" "P-value" "Observations") ///
	fmt(%9.0fc)) substitute(\_ _) order(`coeffs') ///
	mgroups("RVP Phase II" "RFG Federal" "CARB Gasoline", ///
	pattern(1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0) ///
	prefix(\multicolumn{@span}{c}{) suffix(}) ///
	span erepeat(\cmidrule(lr){@span}))


// save data for plotting
keep t_rvpI t_rvpII t_rfg t_CARB Beta* SE*
duplicates drop
save "${support}\stag_estimates.dta", replace


*** Plot ***
local treats rvpII rfg CARB
local windows "3 6 12"
local i = 1
quietly{
	foreach treat of local treats {
		
		local first_treat 1
		
		foreach w of local windows {
				
			if `first_treat'==1 {
				local title = upper("`treat'")
				local first_treat 0
			}
			else {
				local title ""
			}
			
			use "${support}\stag_estimates.dta", clear
			local Beta Beta_`treat'_`w'C
			local SE SE_`treat'_`w'C
			local t t_`treat'
			
			keep `t' `Beta' `SE'
			duplicates drop
			
			sort `t'

			gen lo = `Beta' - 1.96 * `SE'
			gen hi = `Beta' + 1.96 * `SE'

			local xrange (-`w' <= t_`treat' & t_`treat' <= `w')
			
			local main_fmt ///
				msymbol(triangle_hollow) mlcolor(dknavy) msize(medlarge)	///
				mfcolor(gs1) lcolor(dknavy) lwidth(medthick)				///
				connect(l) sort // msize(medium)
			
			local error_fmt lcolor(red) lwidth(thin)

			local graph_fmt													///
				yline(0, lpattern(solid) lcolor(black) lstyle(foreground))	///
				xline(0, lpattern(dash) lcolor(gs6) lstyle(foreground))		///
				legend(off)													///
				title(`title') ytitle(`ytitle') xtitle("")					///
				graphregion(color(white)) bgcolor(white)

			capture graph drop
			twoway ///
				line hi `t' if `xrange', `error_fmt'	 ||	///
				line lo `t' if `xrange', `error_fmt'	 ||	///
				scatter `Beta' `t' if `xrange', `main_fmt'		///
				`graph_fmt'
			
			graph save "${support}\\`treat'`w'", replace
			
		}
	}
}


/// attachments event combined study graph

local gphs "${support}\rvpII3.gph" "${support}\rfg3.gph" 	"${support}\CARB3.gph"	///
		"${support}\rvpII6.gph"	"${support}\rfg6.gph" 	"${support}\CARB6.gph"	///
		"${support}\rvpII12.gph"	"${support}\rfg12.gph"	"${support}\CARB12.gph"
		

graph combine "`gphs'",							///
	rows(3) cols(3) ycommon 					///
	l1title("Period estimates")					///
	b2title("Periods post-treatment (summer months)")			///
	graphregion(color(white))					///
	saving("${support}\fig1", replace)

graph export "${attachments}\fig1.png", replace

/*******************************************************************************
	Part B1. Regression discontinuity data cleaning
*******************************************************************************/

global treatments "RVP RFG CARB RVP_RFG"

foreach treatment of global treatments {
	
	disp "RD cleaning: `treatment'"
	
	use AER20090377_FinalData, clear
	
	// filter treatment
	if "`treatment'" == "RVP" {
		keep if CARBCty==0 & RVPCty==1 & RFGCty==0
	}
	else if "`treatment'" == "RFG" {
		keep if CARBCty==0 & RVPCty==0 & RFGCty==1
	}
	else if "`treatment'" == "CARB" {
		keep if CARBCty==1
	}
	else if "`treatment'" == "RVP_RFG" {
		keep if CARBCty==0 & RVPCty==1 & RFGCty==1
	}

	drop if year > 2003
	drop if valid < 9

	// drop monitors with less than 75% of all observations over date range
	// create Season variable
	gen Season = 1

	// winter 12-2, spring 3-5, summer 6-8, fall 9-11
	forvalues i = 1/3 {
		quiet replace Season = `i'+1 if (`i'*3 <= month) & (month < `i'*3 + 3)
	}

	tab month Season

	// group Dec with next year
	// keep if 75% days observed per Season
	gen year2 = year
	replace year2 = year + 1 if month==12
	
	bysort panelid year2 Season: egen NumDays = count(day)
	keep if NumDays>=0.75 * 30 * 3
	quiet drop NumDays year2

	// drop monitors with less than 75% of all observations over date range
	sum Date
	gen Range = `r(max)' - `r(min)'
	replace Range = 0.75 * Range
	bysort fips site_id: egen NumDays = count(day)
	keep if NumDays>=Range		
	drop NumDays Range
	
	
	// filter monitors for observation dates
	
	if inlist("`treatment'", "RVP", "CARB", "RVP_RFG") {
		
		// 75% three years before RVP in 1992
		
		disp "d(01jun1989)<=Date & Date<d(01jun1992)"
		
		bysort fips site_id: egen NumDays = count(day) ///
			if d(01jun1989)<=Date & Date<d(01jun1992)
			
		replace NumDays = 0 if NumDays==.
		bysort fips site_id: egen NumDays2 = max(NumDays)
		keep if NumDays2 >= 0.75 * 3 * 365
		drop NumDays NumDays2
		
	}
	if inlist("`treatment'", "RFG", "RVP_RFG") {

		// 75% three years between RFG in 1995 and RVP in 1992
	
		disp "d(01jun1992)<=Date & Date<=d(01jun1995)"
		
		bysort fips site_id: egen NumDays = count(day) ///
			if d(01jun1992)<=Date & Date<=d(01jun1995)
			
		replace NumDays = 0 if NumDays==.
		bysort fips site_id: egen NumDays2 = max(NumDays)
		keep if NumDays2 >= 0.75 * 3 * 365
		drop NumDays NumDays2
		
	}
	
	if "`treatment'" == "CARB" {
		
		// 75% four years between start of RVP and CARB
		
		disp "Date>=d(01jun1992) & Date<=d(01jun1995)"
		
		bysort fips site_id: egen NumDays = count(day) ///
			if d(01jun1992)<=Date & Date<d(01jun1996)
			
		replace NumDays = 0 if NumDays==.
		bysort fips site_id: egen NumDays2 = max(NumDays)
		keep if NumDays2 >= 0.75 * 4 * 365
		drop NumDays NumDays2

	}

	
	*** Weather variables ***
	
	// time weather interactions
	gen DOW = dow(Date)
	gen DOY = doy(Date)

	local weather_vars "TempMax TempMin Rain Snow"
	local prefixes1 "_DM _Dm _Dr _Ds"
	local prefixes2 "_MT _Mt _MR _MS"


	// iterate polynomial weather variables
	forval i = 1/4 {
		
		local var `: word `i' of `weather_vars''
		local pfx1 `: word `i' of `prefixes1''
		
		xi i.DOW*`var', prefix(`pfx1')
		drop `pfx1'DOW_*
		
		// iterate polynomials
		forval p = 1/3 {
			
			local pfx2 `: word `i' of `prefixes2''`p'
			
			// polynomial
			gen `var'`p' = `var'^`p'
			
			// generate interaction
			xi i.Season*`var'`p', prefix(`pfx2')
			drop `pfx2'Season*
		}
	}
	// Rain and Snow cubics not used
	drop *Rain3* *Snow3*

	// dropping these creates unestimatable regressions
	// day of week cross season effects
	xi i.DOW*i.Season, prefix(_DSe)

	// month effects
	xi i.month, prefix(_M)
	
	// interactions
	gen TempMaxMin = TempMax * TempMin
	gen RainTempMax = Rain * TempMax

	// lagged temperature
	sort panelid Date
	by panelid: gen tracker = _n

	foreach var of varlist(TempMax TempMin){
		gen `var'L1 = TempMax[_n-1]
		replace `var'L1 = . if tracker==1
	}

	// lagged interactions
	gen TempMaxMaxL1 = TempMax * TempMaxL1
	gen TempMaxMinL1 = TempMax * TempMinL1
	drop tracker


	// remaining interactions
	local weather_vars "TempMaxMin TempMaxL1 TempMinL1 TempMaxMaxL1 TempMaxMinL1 RainTempMax"
	local prefixes "_MTt _MTL _MtL _MTM _MTm _MRT"

	local n_vars : word count `weather_vars'

	forval i = 1/`n_vars' {
		
		local var `: word `i' of `weather_vars''
		local pfx `: word `i' of `prefixes''

		// generate interaction
		xi i.Season*`var', prefix(`pfx')
		drop `pfx'Season*
	}
	
	// log dependent variables
	gen lozone_max = log(ozone_max)
	gen lozone_8hr = log(epa_8hr)

	// drop stuff we don't need
	drop county-regtype
	drop sulfur-rfgtype

	
	
	*** Time Variables ***
	// annualize date
	gen DateS = Date / 365
	egen MaxDate = max(DateS)
	egen MinDate = min(DateS)

	/*	Chebychev polynomials
		https://en.wikipedia.org/wiki/Chebyshev_polynomials#Recurrence_definition
		
		T_0 = 1
		T_1 = Z
		T_{n}(Z) = 2 * Z * T_{n-1}(Z) - T_{n-2}(Z)
	*/

	// 0th and 1st order
	gen Z = 2 * (DateS - MinDate) / (MaxDate - MinDate) - 1
	gen Time0 = 1
	gen Time1 = Z

	// 2nd to 8th order
	forval n = 2/8 {
		
		local n_1 = `n' - 1
		local n_2 = `n' - 2
		
		disp "`n' `n_1' `n_2'"
		gen Time`n' = 2 * Z * Time`n_1' - Time`n_2'
	}

	drop Z Time0
	
	save "${support}\RD_data_`treatment'.dta", replace
}

/*******************************************************************************
	Part B2. Regression discontinuity replication
*******************************************************************************/

// treatment variables with 30 day linear phase-in
local treatvar_RVP TreatRVPII
local treatvar_RFG TreatRFG
local treatvar_RVP_RFG TreatRVPII TreatRFG
local treatvar_CARB TreatRVPII TreatRFG TreatCARB

// global treatments "RVP RFG CARB RVP_RFG"
quietly {
	foreach treatment of global treatments {
		
		local treatvars `treatvar_`treatment''

		
		disp "RD analysis: `treatment' (`treatvars') "
		
		use "${support}\RD_data_`treatment'.dta", clear
		
		// clustering variables
		quiet gen YearMonth = year * 100 + month
		quiet gen YearSeason = year * 100 + Season
		
		
		// regression variables
		foreach var of varlist(`treatvars'  Time1-Time8) {
			quiet gen Beta_`var' = .
			quiet gen SE_`var' = .
			quiet gen T_`var' = .
			quiet gen P_`var' = .
			
			local kind = subinstr("`var'", "Treat", "", .)
			label var Beta_`var' "RD `kind' effect"
			label var SE_`var' "`kind' standard error"
			label var T_`var' "`kind' t-score"
			label var P_`var' "`kind' p-value"
		}

		// unestimatable regression: model overfit?
		if "`treatment'"== "CARB" {
			drop if county_code==71 & site_id==4001
		}
		
		
		// RD estimates by monitor
		egen monitor = group(fips site_id)
		gen PredLogO3 = .

		sum monitor
		local total_monitors = `r(max)'
		
		forvalues i = 1/`total_monitors' {
			
			// RD estimation
			quiet regress lozone_max `treatvars' 			///
				Time1-Time8 _DMDOWXTempM_1-_MRTSeaXRainT_4	///
				if monitor==`i', vce(cluster YearSeason)
				
			// predict log ozone
			quiet predict tempvar
			quiet replace PredLogO3 = tempvar if monitor==`i'
			drop tempvar
		
			
			// save regression results
			foreach var of varlist(`treatvars'  Time1-Time8  ) {
				quiet replace Beta_`var' = _b[`var'] if monitor == `i'
				quiet replace SE_`var' = _se[`var'] if monitor == `i'
				quiet replace T_`var' = Beta_`var' /  SE_`var' if monitor == `i'

				// calculate p-value
				// https://www.statalist.org/forums/forum/general-stata-discussion/general/1308562-capture-p-values-from-regression-results
				quiet replace P_`var' = 2*ttail(e(df_r), abs(T_`var'))
		 
			}
		}
		

		// save estimates
		keep state state_code county_code site_id fips Beta* SE* T_* P_*
		duplicates drop
		
		// add county name
		sort state_code county_code site_id
		merge state_code county_code using AER20090377_CountyList.dta, uniqusing
		keep if _merge==3
		drop _merge
		
		keep state state_code county_desc county_code site_id Beta* SE* T_* P_*
		sort state_code county_code site_id	

		label var state "State"
		label var county_desc "County"
		label var site_id "Monitor ID"

		drop if state == ""
		save "${support}\tableA_`treatment'.dta", replace
		
		// summarize 
		noisily {
			sum Beta_Treat*
		}
		
		// round
		foreach var of varlist(Beta* SE* T_* P_*) {
			
			replace `var' = . if `var' == 0
			replace `var' = round(`var', 0.001)
		}

		// add stars
		foreach beta of varlist(Beta_Treat*){

			local suffix = subinstr("`beta'", "Beta", "", .)
			
			
			replace P`suffix' = . if P`suffix' == 0
			
			gen Beta`suffix'_str = string(Beta`suffix')
			
			replace Beta`suffix'_str = Beta`suffix'_str + "*" if P`suffix' < 0.01
			replace Beta`suffix'_str = Beta`suffix'_str + "*" if P`suffix' < 0.05
			replace Beta`suffix'_str = Beta`suffix'_str + "*" if P`suffix' < 0.1
			
			
			local labelname : variable label Beta`suffix'
			label var Beta`suffix'_str "`labelname'"
			
			drop Beta`suffix'
			rename Beta`suffix'_str Beta`suffix'
		}
		
		// add index
		gen index = _n
		label var index "Index"
		
		// attachments
		noisily texsave index state county_desc site_id  ///
			Beta_Treat* SE_Treat* T_Treat* P_Treat* 	///
			using "${attachments}\tableA_`treatment'.tex",	///
			replace varlabels
	}
}

*** Create kernel density plots ***

global treatments "RVP RFG CARB"

foreach treatment of global treatments {

	use "${support}\tableA_`treatment'.dta", clear

	sum Beta_Treat`treatment'
	local m = round(r(mean), 0.0001)
	local x = `m'+0.1
	
	hist Beta_Treat`treatment', text(8 `x' "mean = `m'") 			///
		xline(`m', lpattern(dash) lcolor(gs6) lstyle(foreground)) 	///
		addplot(kdensity Beta_Treat`treatment' , bwidth(0.05) lwidth(thick)) ///
		title("") legend(off)										///
		graphregion(color(white)) bgcolor(white)					///
		saving("${support}\RD_hist_`treatment'.gph", replace)
}

local gphs "${support}\RD_hist_RVP.gph" "${support}\RD_hist_RFG.gph" "${support}\RD_hist_CARB.gph"
	
graph combine "`gphs'",						///
	rows(3) xcommon ycommon 				///
	graphregion(color(white)) 				///
	saving("${support}/fig2.gph", replace)
	
graph export "${attachments}/fig2.png", replace



/*******************************************************************************
	Part B3. Regression discontinuity robustness checks
	
	Regression Discontinuity in Time: Considerations for Empirical Applications
	https://www.nber.org/papers/w23602
*******************************************************************************/

global treatments "RVP RFG CARB"

// treatment variables with 30 day linear phase-in
local treatvar_RVP TreatRVPII
local treatvar_RFG TreatRFG
local treatvar_CARB TreatCARB

// binary treatment
// local treatvar_RVP treat_rvpII
// local treatvar_RFG treat_rfg
// local treatvar_CARB treat_CARB

// RD window: 6 months
global RD_window = 3/12
	
local RDcontrols _DMDOWXTempM_1-_MRTSeaXRainT_4

est clear

// quiet {
foreach treat of global treatments {
	
	noisily disp "RD robustness: `treat' (`treatvar_`treat'') "
	
	use "${support}\RD_data_`treat'.dta", clear
	

	gen YearMonth = year * 100 + month
	gen YearSeason = year * 100 + Season
	
	sort panelid Date
	
	gen treatvar = `treatvar_`treat''
	
	
	// get first (min) month of treatment
	gen temp_treat = .
	replace temp_treat = 1 if treatvar == 1
	sort panelid Date
	by panelid: egen start_treat = min(Date * temp_treat)
	drop temp_treat

	
	// treatment starts at t_treat==0, annualized
	gen t_treat =  (Date - start_treat) / 365
	gen txtreat = t_treat * treatvar
	
	// quadratic
	gen t_treat2 = t_treat^2
	gen txtreat2 = t_treat2 * treatvar
	
	// cubic
	gen t_treat3 = t_treat^3
	gen txtreat3 = t_treat3 * treatvar
	
	// combine polynomials
	local run_var 	t_treat txtreat
	local run_var2 `run_var'2 t_treat2 txtreat2
	local run_var3 `run_var'3 t_treat3 `txtreat3
	
	// linear trend
	eststo: reg lozone_max treatvar t_treat 		///
		`RDcontrols'								///
		if abs(t_treat) <= $RD_window, 				///
		vce(cluster YearSeason)
	
	estadd local Order "1"
	estadd local Slope "No"
	estadd local Ctrls "Yes"
	
	predict tempvar
	gen Pred_lozone_0 = tempvar
	drop tempvar
		
	
	// linear with differential slopes
	eststo: reg lozone_max treatvar `run_var' 				///
		`RDcontrols'										///
		if abs(t_treat) <= $RD_window, ///
		vce(cluster YearSeason)
	
	estadd local Order "1"
	estadd local Slope "Yes"
	estadd local Ctrls "Yes"
	
	predict tempvar
	gen Pred_lozone_1 = tempvar
	drop tempvar
	
	// quadratic with differential slopes
	eststo: reg lozone_max treatvar `run_var' `run_var'2	///
		`RDcontrols'										///
		if abs(t_treat) <= $RD_window, ///
		vce(cluster YearSeason)
		
	estadd local Order "2"
	estadd local Slope "Yes"
	estadd local Ctrls "Yes"
	
	predict tempvar
	gen Pred_lozone_2 = tempvar
	drop tempvar
	
	// cubic with differential slopes
	eststo: reg lozone_max treatvar `run_var' `run_var'2 `run_var'3	///
		`RDcontrols'												///
		if abs(t_treat) <= $RD_window, ///
		vce(cluster YearSeason)
		
	estadd local Order "3"
	estadd local Slope "Yes"
	estadd local Ctrls "Yes"
	
	predict tempvar
	gen Pred_lozone_3 = tempvar
	drop tempvar
		
	save "${support}\RD_plot_`treat'.dta", replace 
	
}
//}

label var treatvar "Treatment Estimate"

// attachments regression
noisily esttab using "${attachments}\table3.tex", replace f   ///
	keep(treatvar) ///
	nocons se(3) b(3) label booktabs nomtitle collabels(none) compress ///
	stats(Order Slope Ctrls N, ///
	label("Polynomial Order" "Differential Slope" "Full controls" "Observations") ///
	fmt(%9.0fc)) ///
	mgroups("RVP Phase II" "RFG Federal" "CARB Gasoline", ///
	pattern(1 0 0 0 1 0 0 0 1 0 0 0) ///
	prefix(\multicolumn{@span}{c}{) suffix(}) ///
	span erepeat(\cmidrule(lr){@span}))


/* Plot of grouped RD predictions
global treatments "RVP RFG CARB"
local first_row = 0

foreach treat of global treatments{

	local first_treat = 1
	
	use "${support}\RD_plot_`treat'.dta", clear

	local xrange abs(t_treat) <= $RD_window
	
	forvalues i = 0/3 {
		
		if `first_treat'==1 {
			local title "`treat'"
			local first_treat = 0
		}
		else {
			local title ""
		}
		
		twoway ///
			scatter lozone_max t_treat if `xrange', msize(tiny) color(gs12) || 	///
			line Pred_lozone_`i' t_treat if `xrange', sort lwidth(medthick)		///
			xline(0, lpattern(dash) lcolor(gs6) lstyle(foreground))				///
			legend(off) xtitle("") title(`title') 								///
			graphregion(color(white)) bgcolor(white)							///
			saving("${support}\\`treat'_RD`i'", replace)
	}
}

local gphs ///
	"${support}\RVP_RD0.gph" "${support}\RFG_RD0.gph" "${support}\CARB_RD0.gph"	///
	"${support}\RVP_RD1.gph" "${support}\RFG_RD1.gph" "${support}\CARB_RD1.gph"	///
	"${support}\RVP_RD2.gph" "${support}\RFG_RD2.gph" "${support}\CARB_RD2.gph"	///
	"${support}\RVP_RD3.gph" "${support}\RFG_RD3.gph" "${support}\CARB_RD3.gph"
		
graph combine "`gphs'",							///
	rows(4) cols(3) ycommon 					///
	l2title("Log ozone")						///
	b2title("Periods post-treatment (years)")	///
	graphregion(color(white))					///
	saving("${attachments}\fig3", replace)

graph export "${attachments}\fig3.png", replace

*/

log close

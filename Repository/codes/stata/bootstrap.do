** this dofile bootstraps the estimation of the trade elasticities and policy parameters
* and produces table A1: Tab_basic.tex



clear
set more off
set matsize 10000
set type double


cd "Repository/estimation"


capture log close




************prep the regression dataset with AGGREGATED sectors ********************
*** tariffs from WITS-TRAINS
use ../rawdata/alltariffs_wiodsample_agg, clear
drop if year<2000
g sec=1
g tau_n=tariff/100
replace tau_n=1+tau_n
label var tau_n "1 plus tariffs, mean"


merge 1:1 iso1 iso2 sec year using  ../rawdata/schengen_new_agg_notau

***services sectors with missing tariffs
tab sec if _m==2
replace sec=2 if sec==.
replace tau_n=1 if sec==2 

drop _m

sort iso1 iso2 sec year

egen id1=group(iso1)
egen id2=group(iso2)
egen idn=group(id1 id2)


*** ln tariffs
g ltau_n=ln(tau_n)

replace rta_c=0 if bothEEA_c==1

g nEEA = bothEU
replace nEEA = 1 if bothEEA==1

egen eu1=max(bothEU), by(iso1 year)
egen eu2=max(bothEU), by(iso2 year)
egen ea1=max(bothEEA), by(iso1 year)
egen ea2=max(bothEEA), by(iso2 year)	
replace ea1=1 if iso1=="CHE"
replace ea2=1 if iso2=="CHE"

replace nEEA=1 if eu1==1 & ea2==1
replace nEEA=1 if eu2==1 & ea1==1


replace nEEA=0 if iso1==iso2
replace rta_c=0 if iso1==iso2
replace eupta=0 if iso1==iso2

replace rta_c=0 if nEEA==1 & (iso1=="CHE" | iso2=="CHE")
replace rta_c=0 if nEEA==1 & (iso1=="NOR" | iso2=="NOR")


drop bothEEA
rename nEEA  bothEEA



keep iso1 iso2 imp_int  intst ltau* bothEU botheuro rta_c eukor year sec id1 id2 idn eupta
rename bothEU bothEEA

save eucol_nnew_dataset_IH, replace







************ prep the regression data with 50 SECTORS********************
*** tariffs from WITS-TRAINS
use  ../rawdata/alltariffs_wiodsample, clear
drop if year<2000
drop if sec>22
g tau_n=tariff/100
replace tau_n=1+tau_n
label var tau_n "1 plus tariffs, simple mean"
sort iso1 iso2 sec year
rename sec_new sector
merge 1:1 iso1 iso2 sec year using  ../rawdata/schengen_new_notau

***services sectors with missing tariffs, internal trade with missing tariffs
tab sec if _m==2

replace tau_n=1 if sec>22
drop _m
replace tau=1 if bothEU==1

sort iso1 iso2 sec year
duplicates drop

egen id1=group(iso1)
egen id2=group(iso2)
egen idn=group(id1 id2)


*** ln tariffs
g ltau_n=ln(tau_n)
drop tariff tau_n

*undo rta_c
replace rta_c=0 if bothEEA_c==1

g nEEA = bothEU
replace nEEA = 1 if bothEEA==1

egen eu1=max(bothEU), by(iso1 year)
egen eu2=max(bothEU), by(iso2 year)
egen ea1=max(bothEEA), by(iso1 year)
egen ea2=max(bothEEA), by(iso2 year)	
replace ea1=1 if iso1=="CHE"
replace ea2=1 if iso2=="CHE"

replace nEEA=1 if eu1==1 & ea2==1
replace nEEA=1 if eu2==1 & ea1==1


replace nEEA=0 if iso1==iso2
replace rta_c=0 if iso1==iso2
replace eupta=0 if iso1==iso2

replace rta_c=0 if nEEA==1 & (iso1=="CHE" | iso2=="CHE")
replace rta_c=0 if nEEA==1 & (iso1=="NOR" | iso2=="NOR")

drop bothEEA
rename nEEA  bothEEA

keep iso1 iso2 imp_int  intst ltau* bothEU botheuro rta_c eukor year sec id1 id2 idn eupta
rename bothEU bothEEA

save eucol_nnew_dataset_sec_IH, replace





global coeff _b[ltau_n] _b[bothEEA] _b[intst] _b[eukor] _b[rta_c] _b[botheuro] _b[eupta]
global acoeff  _b[bothEEA] _b[intst] _b[eukor] _b[rta_c] _b[botheuro] _b[eupta]

forvalues s=1(1)22 {
capture log close
log using bt_sec`s', replace
use eucol_nnew_dataset_sec_IH, clear
keep if sec==`s'
g Aimp = imp_int*exp(ltau_n)
di "`s'"

bootstrap $coeff, seed(1986) cluster(idn) reps(100) size(1400) saving(draws`s', replace): ppmlhdfe imp_int ltau_n  bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn)
bootstrap $acoeff, seed(1986) cluster(idn) reps(1000) size(1400) saving(adraws`s', replace): ppmlhdfe Aimp bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn)
bootstrap $coeff, seed(1986) cluster(idn) reps(1000) size(1400) saving(draws`s'_t, replace): ppmlhdfe imp_int ltau_n bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn##c.year)
bootstrap $acoeff, seed(1986) cluster(idn) reps(1000) size(1400) saving(adraws`s'_t, replace): ppmlhdfe Aimp bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn##c.year)

ppmlhdfe imp_int ltau_n  bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn) cluster(idn)
ppmlhdfe Aimp  bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn) cluster(idn)
ppmlhdfe imp_int ltau_n bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn##c.year) cluster(idn)
ppmlhdfe Aimp  bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn##c.year) cluster(idn)
}




forvalues s=23(1)26 {
capture log close
log using bt_sec`s', replace
use eucol_nnew_dataset_sec_IH, clear
keep if sec==`s'
g Aimp = imp_int*exp(ltau_n)
di "`s'"
bootstrap $acoeff, seed(1986) cluster(idn) reps(1000) size(1400)  saving(adraws`s', replace): ppmlhdfe Aimp bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn)
bootstrap $acoeff, seed(1986) cluster(idn) reps(1000) size(1400) saving(adraws`s'_t, replace): ppmlhdfe Aimp bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn##c.year)

ppmlhdfe imp_int   bothEEA intst eukor rta_c  botheuro eupta , a(id1#year id2#year idn) cluster(idn)
ppmlhdfe Aimp  bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn##c.year) cluster(idn)
}





*aggregate regression to back out services elasticity


use eucol_nnew_dataset_IH, clear
keep if sec==1
ppmlhdfe imp_int ltau_n  bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn) cluster(idn)
di _b[ltau_n]+2.026
ppmlhdfe imp_int ltau_n  bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn##c.year) cluster(idn)
di _b[ltau_n]+2.026

g impA = imp_int*exp(ltau)
ppmlhdfe impA ltau_n  bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn) cluster(idn)
di _b[ltau_n]+2.026


** table A1: Aggregate results
use eucol_nnew_dataset_IH, clear
capture log close
log using log_aggregate, replace
estimates clear
keep if sec==1
eststo: ppmlhdfe imp_int  bothEEA , a(id1#year id2#year idn) cluster(idn)
*eststo: ppmlhdfe imp_int ltau_n  bothEEA, a(id1#year id2#year idn) cluster(idn)
eststo: ppmlhdfe imp_int ltau_n  bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn) cluster(idn)

use eucol_nnew_dataset_IH, clear
 keep if sec==2
eststo: ppmlhdfe imp_int  bothEEA , a(id1#year id2#year idn) cluster(idn)
eststo: ppmlhdfe imp_int  bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn) cluster(idn)
 
 
estout *, cells(b(star fmt(3)) se(par fmt(2))) drop() starlevels(* 0.1 ** 0.05 *** 0.01) legend label varlabels(_cons Constant) stats(N N_g ll chi2, fmt(%9.0g %9.0g %9.3f %9.3f) label (Observations CountryPairs Loglikelihood Chi2 )) style(fixed)


esttab * using ../tables/Tab_basic.tex, coeflabels(bothEEA "Both Single Market" intst "Schengen" eukor "EU-KOR PTA" botheuro "Both Euro" eupta "EU PTAs" rta_c "Other RTAs" ltau_n "Tariffs") drop(_cons) replace cells(b(star fmt(3)) se(par fmt(2)))  starlevels(* 0.1 ** 0.05 *** 0.01) noobs order(bothEEA botheuro intst eukor eupta rta_c ltau_n) fragment nomtitle 


 
log close


**** bootstrap aggregate results

*Manufacturing
global coeff _b[ltau_n] _b[bothEEA] _b[intst] _b[eukor] _b[rta_c] _b[botheuro] _b[eupta]
global acoeff  _b[bothEEA] _b[intst] _b[eukor] _b[rta_c] _b[botheuro] _b[eupta]

capture log close
log using bt_secManu, replace
use eucol_nnew_dataset_IH, clear
keep if sec==1
g Aimp = imp_int*exp(ltau_n)
di "`s'"

bootstrap $coeff, seed(1986) cluster(idn) reps(1000) saving(drawsManu, replace): ppmlhdfe imp_int ltau_n  bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn)
bootstrap $acoeff, seed(1986) cluster(idn) reps(1000) saving(adrawsManu, replace): ppmlhdfe Aimp bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn)

ppmlhdfe imp_int ltau_n  bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn) cluster(idn)
ppmlhdfe Aimp  bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn) cluster(idn)

capture log close
log using bt_secSERV, replace
use eucol_nnew_dataset_IH, clear
keep if sec==2
g Aimp = imp_int*exp(ltau_n)
di "`s'"
bootstrap $acoeff, seed(1986) cluster(idn) reps(1000) saving(adrawSERV, replace): ppmlhdfe Aimp bothEEA intst eukor rta_c  botheuro eupta, a(id1#year id2#year idn)

ppmlhdfe imp_int   bothEEA intst eukor rta_c  botheuro eupta , a(id1#year id2#year idn) cluster(idn)










*** process output main bootstrap

forvalues i=1(1)22 {

u draws`i', clear
rename _bs_* BS*_`i'
g line = _n
save temp, replace

u adraws`i', clear
rename _bs_* ABS*
forvalues l=1(1)6 {
	local j=`l'+1
	g ABS`j'_`i' =ABS`l'
	drop ABS`l'
	}

g line = _n
merge 1:1 line using temp
drop _merge

egen dropme=rowmiss(BS* ABS*)
drop if dropme!=0
drop dropme

forvalues k=2(1)6 {
replace BS`k'=ABS`k' if BS1>-1
}

drop ABS*
replace BS1=-1 if BS1>-1
drop line
g line=_n
save ddraws`i', replace
}



forvalues i=23(1)50 {
u adraws`i', clear
capture drop _b_cons
capture rename _b_* _bs_*
local j=2
foreach var of varlist _bs_* {
	rename `var' BS`j'_`i'
local j=`j'+1	
}
	
g BS1_`i'=-1.6189
g line=_n
save ddraws`i', replace
}


u ddraws1, clear
	forvalues i=2(1)50 {
	di `i'
		merge 1:1 line using ddraws`i'
		drop _merge
	}
		keep if line<1001
	
export delim BS1* using ../simulation/input/BDraws/sectoral/bltau-sec-nnew.txt, replace delim(tab)
export delim BS2* using ../simulation/input/BDraws/sectoral/bbothEEA-sec-nnew.txt, replace delim(tab)
export delim BS3* using ../simulation/input/BDraws/sectoral/bschengen-sec-nnew.txt, replace delim(tab)
export delim BS4* using ../simulation/input/BDraws/sectoral/beukor-sec-nnew.txt, replace delim(tab)
export delim BS5* using ../simulation/input/BDraws/sectoral/bbothrta-sec-nnew.txt, replace delim(tab)
export delim BS6* using ../simulation/input/BDraws/sectoral/bbotheuro-sec-nnew.txt, replace delim(tab)

rename BS1* bltau*
rename BS2* bbothEEA*
rename BS3* bschengen*
rename BS4* beukor*
rename BS5* bbothrta*
rename BS6* bbotheuro*
rename BS7* beupta*

rename line draw
order draw
save coef_final_adj, replace


******************** make dataset with unadjusted coefficients


forvalues i=1(1)22 {
u draws`i', clear
rename _bs_* BS*_`i'
egen dropme=rowmiss(BS*)
drop if dropme!=0
drop dropme
g line = _n
save temp`i', replace
}

u temp1, clear
forvalues i=2(1)22 {
merge 1:1 line using temp`i'
drop _mer
erase temp`i'.dta
}
erase temp1.dta

forvalues i=23(1)50 {
merge 1:1 line using ddraws`i', keepusing(BS2* BS3* BS4* BS5* BS6* BS7*)
drop _mer
}
rename line draw
order draw
drop if draw>1000

rename BS1* bltau*
rename BS2* bbothEEA*
rename BS3* bschengen*
rename BS4* beukor*
rename BS5* bbothrta*
rename BS6* bbotheuro*
rename BS7* beupta*
save coef_final_org, replace






***************************************************************************
**************    prepare treatment coefficients    ***********************
******                 SECTORAL ESTIMATIONS                   *************
*****			means of 1000 cleaned sectoral draws			***********
***************************************************************************



import delim ../simulation/input/BDraws/sectoral/bltau-sec-nnew.txt, varnames(1) clear
g line=_n
save bltau-sec-nnewB, replace
drop line
collapse (mean) _all,
g draw=_n
reshape long bs1_, i(draw) j(sec)
drop draw
outsheet b* using "../simulation/input/bltau_smdB.txt", replace nonames delimiter(",")



local vars "schengen botheuro botheea eukor bothrta"

foreach x of local vars {

import delim ../simulation/input/BDraws/sectoral/b`x'-sec-nnew.txt, varnames(1) clear

collapse (mean) _all
g draw=_n
*rename bs_* b*
reshape long bs, i(draw) j(sec) string
replace sec=substr(sec,3,.)
destring sec, replace
drop draw
outsheet b* using "../simulation/input/b`x'_smdB.txt", replace nonames delimiter(",")
}




******************************************************************************************
* the same with time trends **************************************************************
******************************************************************************************



forvalues i=1(1)22 {
u draws`i'_t, clear
rename _bs_* BS*_`i'
g line = _n
save temp, replace

u adraws`i'_t, clear
rename _bs_* ABS*
forvalues l=1(1)6 {
	local j=`l'+1
	g ABS`j'_`i' =ABS`l'
	drop ABS`l'
	}

g line = _n
merge 1:1 line using temp
drop _merge

egen dropme=rowmiss(BS* ABS*)
drop if dropme!=0
drop dropme

forvalues k=2(1)6 {
replace BS`k'=ABS`k' if BS1>-1
}

drop ABS*
replace BS1=-1 if BS1>-1
drop line
g line=_n
save ddraws`i'_t, replace
}





forvalues i=23(1)50 {
u adraws`i'_t, clear
capture drop *cons
capture rename _b_* _bs_*
local j=2
foreach var of varlist _bs_* {
	rename `var' BS`j'_`i'
local j=`j'+1	
}
g BS1_`i'=-1.6189
	
egen dropme=rowmiss(BS*)
drop if dropme!=0
drop dropme
g line=_n
save ddraws`i'_t, replace
}


u ddraws1_t, clear
	forvalues i=2(1)50 {
	di `i'
		merge 1:1 line using ddraws`i'_t
		drop _merge
	}
		keep if line<1001
	
export delim BS1* using ../simulation/input/BDraws/sectoral-trend/bltau-sec-nnew.txt, replace delim(tab)
export delim BS2* using ../simulation/input/BDraws/sectoral-trend/bbothEEA-sec-nnew.txt, replace delim(tab)
export delim BS3* using ../simulation/input/BDraws/sectoral-trend/bschengen-sec-nnew.txt, replace delim(tab)
export delim BS4* using ../simulation/input/BDraws/sectoral-trend/beukor-sec-nnew.txt, replace delim(tab)
export delim BS5* using ../simulation/input/BDraws/sectoral-trend/bbothrta-sec-nnew.txt, replace delim(tab)
export delim BS6* using ../simulation/input/BDraws/sectoral-trend/bbotheuro-sec-nnew.txt, replace delim(tab)



rename BS1* bltau*
rename BS2* bbothEEA*
rename BS3* bschengen*
rename BS4* beukor*
rename BS5* bbothrta*
rename BS6* bbotheuro*
rename BS7* beupta*

rename line draw
order draw
save coef_trend_final_adj, replace



******************** make dataset with unadjusted coefficients


forvalues i=1(1)22 {
u draws`i'_t, clear
rename _bs_* BS*_`i'
*egen dropme=rowmiss(BS*)
*drop if dropme!=0
*drop dropme
g line = _n
save temp`i', replace
}

u temp1, clear
forvalues i=2(1)22 {
merge 1:1 line using temp`i'
drop _mer
erase temp`i'.dta
}
erase temp1.dta

forvalues i=23(1)50 {
merge 1:1 line using ddraws`i'_t, keepusing(BS2* BS3* BS4* BS5* BS6* BS7*)
drop _mer
}
rename line draw
order draw
drop if draw>1000

rename BS1* bltau*
rename BS2* bbothEEA*
rename BS3* bschengen*
rename BS4* beukor*
rename BS5* bbothrta*
rename BS6* bbotheuro*
rename BS7* beupta*
save coef_trend_final_org, replace





***************************************************************************
**************    prepare treatment coefficients    ***********************
******                 SECTORAL ESTIMATIONS                   *************
*****			means of 1000 cleaned sectoral draws			***********
***************************************************************************



import delim ../simulation/input/BDraws/sectoral-trend/bltau-sec-nnew.txt, varnames(1) clear
g line=_n
save bltau-sec-nnewBT, replace
drop line
collapse (mean) _all,
g draw=_n
reshape long bs1_, i(draw) j(sec)
drop draw
outsheet b* using "../simulation/input/bltau_smdBT.txt", replace nonames delimiter(",")



local vars "schengen botheuro botheea eukor bothrta"

foreach x of local vars {

import delim ../simulation/input/BDraws/sectoral-trend/b`x'-sec-nnew.txt, varnames(1) clear
g line=_n
collapse (mean) _all
g draw=_n
*rename bs_* b*
reshape long bs, i(draw) j(sec) string
replace sec=substr(sec,3,.)
destring sec, replace
drop draw
outsheet b* using "../simulation/input/b`x'_smdBT.txt", replace nonames delimiter(",")
}






**************** prepare aggregate coefficients ******************


u drawsMANU, clear
rename _bs_* BS*_1
g line = _n
save temp, replace

u adrawsMANU, clear
rename _bs_* ABS*
forvalues l=1(1)6 {
	local j=`l'+1
	g ABS`j'_1 =ABS`l'
	drop ABS`l'
	}

g line = _n
merge 1:1 line using temp
drop _merge

egen dropme=rowmiss(BS* ABS*)
drop if dropme!=0
drop dropme

forvalues k=2(1)6 {
replace BS`k'=ABS`k' if BS1>-1
}

drop ABS*
replace BS1=-1 if BS1>-1
drop line
g line=_n
save ddrawsMANU, replace




u adrawSERV, clear
local j=2
foreach var of varlist _bs_* {
	rename `var' BS`j'_2
local j=`j'+1	
}

	
g BS1_2=-1.6189
	
egen dropme=rowmiss(BS*)
drop if dropme!=0
drop dropme
g line=_n
save ddrawsSERV, replace


u ddrawsMANU, clear

		merge 1:1 line using ddrawsSERV
		drop _merge
		keep if line<1001
		
		
		rename BS1* bltau*
rename BS2* bbothEEA*
rename BS3* bschengen*
rename BS4* beukor*
rename BS5* bbothrta*
rename BS6* bbotheuro*
rename BS7* beupta*

save ddrawsMANUSERV, replace



// this do-file compiles the 
*- 1) tariffs for the estimation datasets
*----- alltariffs_wiodsample.dta
*----- alltariffs_wiodsample_agg.dta


*- 2) tariffs for the simulation
*----- appl_tau_2014_weighted.dta
*----- mfntariffs_weight_fill_wiod.dta
*----- fin_cftau_2014_weighted_wiod.dta


* using tariff data downloaded from https://wits.worldbank.org/WITS/WITS/Restricted/Login.aspx
* and trade data downloaded from downloaded from http://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=37



cd "Repository"
mkdir temp
cd temp


************************************************************************************************
**************** TARIFFS FOR ESTIMATION ********************************************************
************************************************************************************************


*********************************************************************************
*********************************************************************************
*** Fill Preferential Tariffs
*********************************************************************************
*********************************************************************************

clear
set more off
set processors 32


global pref "../rawdata/tariffs/yearly/"
global mfn "../rawdata/tariffs/bilateralMFN/"


*-----------------------------------------------------------------------------------
*-----------------------------------------------------------------------------------
* 1) Combine MFN and Preferential  Tariffs
*-----------------------------------------------------------------------------------
*-----------------------------------------------------------------------------------

forvalues t = 1995/2014{
use "${mfn}bil_mfn`t'", clear

drop wto* eu*


merge 1:1 iso1 iso2 hs92 year using "${pref}bil_pref`t'"
rename _m merge
forvalues i = 1/2 {
rename iso`i' iso3
merge m:1 iso3 using "../rawdata/WTO_membership"
drop if _m == 2
rename wto wto`i'
rename iso3 iso`i'
drop _m
label var wto`i' "Accession-WTO `i'"
}

tab merge 
* no merge == 2!
drop if merge == 2
drop merge

gen tariff = pref
replace tariff = mfn1 if pref == . & mfn1 != .
label var tariff "Tariff (Pref or MFN)"

order iso* year hs tariff 
compress
save "tariff`t'", replace
}


*-----------------------------------------------------------------------------------
*-----------------------------------------------------------------------------------
* 2) Add EU accession year
*-----------------------------------------------------------------------------------
*-----------------------------------------------------------------------------------

forvalues t = 1995/2014{
use "tariff`t'", clear

*** generate eu variable
drop eu1 eu2
gen eu1 = .
gen eu2 = .
forvalues x = 1/2{
local iso BEL FRA ITA LUX NLD DEU
foreach i of local iso {
replace eu`x' = 1957 if iso`x' == "`i'"
}
local iso DNK IRL GBR
foreach i of local iso {
replace eu`x' = 1973 if iso`x' == "`i'"
}
local iso GRC
foreach i of local iso {
replace eu`x' = 1981 if iso`x' == "`i'"
}
local iso PRT ESP
foreach i of local iso {
replace eu`x' = 1986 if iso`x' == "`i'"
}
local iso AUT FIN SWE
foreach i of local iso {
replace eu`x' = 1995 if iso`x' == "`i'"
}
local iso CYP CZE EST HUN LVA LTU MLT POL SVK SVN
foreach i of local iso {
replace eu`x' = 2004 if iso`x' == "`i'"
}
local iso BGR ROU
foreach i of local iso {
replace eu`x' = 2007 if iso`x' == "`i'"
}
local iso HRV
foreach i of local iso {
replace eu`x' = 2012 if iso`x' == "`i'"
}
}
label var eu1 "Accession-EU 1"
label var eu2 "Accession-EU 2"

replace pref = 0 if eu1 <= year & eu2 <= year

replace iso1="ROU" if iso1=="ROM"
replace iso2="ROU" if iso2=="ROM"

save "prel1_`t'", replace
}




**** generate WIOD sample for 2000-2014
forvalues t = 2000/2014{
use ../rawdata/countrysample_wiod.dta, clear
rename iso iso1
drop if iso1=="ROW"
merge 1:m iso1 using prel1_`t'.dta
drop if _m==2
g wiod1=0
replace wiod1=1 if _m==3
drop _m eu1 eu2 wto1 wto2
sort iso2
save wiod_tariff`t', replace

use ../rawdata/countrysample_wiod.dta, clear
rename iso iso2
drop if iso2=="ROW"
merge 1:m iso2 using wiod_tariff`t'.dta
drop if _m==2 
g wiod2=0
replace wiod2=1 if _m==3
drop _m
keep if wiod1==1 & wiod2==1
sort iso1 iso2 year
save wiod_tariff`t', replace
}



* translate H0 to ISIC3
forvalues t=2000/2014{
insheet using "../rawdata/JobID-6_Concordance_H0_to_I3.csv", clear
drop hs198892productdescription isicrevision3productdescription
rename hs198892productcode hs92
rename isicrevision3productcode isicr3

merge 1:m hs92 using wiod_tariff`t'
tab hs92 if _m==2
tab iso2 if _m==2
keep if _m==3
drop _m
save help, replace

* translate ISIC3 to ISIC3.1
insheet using "../rawdata/ISIC_Rev_3-ISIC_Rev_3_1_correspondence.txt", clear
drop partial3 partial31 activity
rename rev3 isicr3
rename rev31 isicr31
drop if isicr3=="n/a"
destring isicr3, replace
joinby isicr3 using help, unm(b)
tab _m
tab isicr3 if _m==2
tab iso2 if _m==2
keep if _m==3
drop _m
save help, replace

* translate ISIC3.1 to ISIC4
insheet using "../rawdata/ISIC31_ISIC4.txt", clear
drop partialisic31 partialisic4 detail
rename isic31code isicr31
rename isic4code isicr4
joinby isicr31 using help, unm(b)
tab _m
tab isicr31 if _m==2
tab iso2 if _m==2
keep if _m==3
drop _m
save help, replace

tostring isicr4, replace
replace isicr4 = "0" + isicr4 if length(isicr4)==3
g isic=substr(isicr4,1,2)

g sec=0

*Agriculture
replace sec=1 if isic=="01" 
replace sec=2 if isic=="02" 
replace sec=3 if isic=="03" 
*Mining and Quarrying
replace sec=4 if isic=="05" |  isic=="06" | isic=="07" | isic=="08" | isic=="09"
*Manufacturing
replace sec=5 if isic=="10" |  isic=="11" | isic=="12"
replace sec=6 if isic=="13" |  isic=="14" | isic=="15"
replace sec=7 if isic=="16" 
replace sec=8 if isic=="17" 
replace sec=9 if isic=="18" 
replace sec=10 if isic=="19" 
replace sec=11 if isic=="20" 
replace sec=12 if isic=="21" 
replace sec=13 if isic=="22" 
replace sec=14 if isic=="23" 
replace sec=15 if isic=="24" 
replace sec=16 if isic=="25" 
replace sec=17 if isic=="26" 
replace sec=18 if isic=="27" 
replace sec=19 if isic=="28" 
replace sec=20 if isic=="29" 
replace sec=21 if isic=="30" 
replace sec=22 if isic=="31" |  isic=="32"
replace sec=23 if isic=="33" 
*Services
replace sec=24 if isic=="35" 
replace sec=25 if isic=="36" 
replace sec=26 if isic=="37" |  isic=="38" |  isic=="39"
replace sec=27 if isic=="41" |  isic=="42" |  isic=="43"
replace sec=28 if isic=="45" 
replace sec=29 if isic=="46" 
replace sec=30 if isic=="47"
replace sec=31 if isic=="49" 
replace sec=32 if isic=="50" 
replace sec=33 if isic=="51" 
replace sec=34 if isic=="52"
replace sec=35 if isic=="53" 
replace sec=36 if isic=="55" |  isic=="56"
replace sec=37 if isic=="58" 
replace sec=38 if isic=="59" | isic=="60"
replace sec=39 if isic=="61" 
replace sec=40 if isic=="62" | isic=="63"
replace sec=41 if isic=="64" 
replace sec=42 if isic=="65" 
replace sec=43 if isic=="66" 
replace sec=44 if isic=="68" 
replace sec=45 if isic=="69" | isic=="70"
replace sec=46 if isic=="71" 
replace sec=47 if isic=="72" 
replace sec=48 if isic=="73" 
replace sec=49 if isic=="74" | isic=="75"
replace sec=50 if isic=="77" | isic=="78" | isic=="79" | isic=="80" | isic=="81" | isic=="82"
replace sec=51 if isic=="84" 
replace sec=52 if isic=="85" 
replace sec=53 if isic=="86" | isic=="87" | isic=="88"
replace sec=54 if isic=="90" | isic=="91" | isic=="92" | isic=="93" | isic=="94" | isic=="95" | isic=="96"
replace sec=55 if isic=="97" | isic=="98" 
replace sec=56 if isic=="99" 

*** no tau in services
tab isicr4 if sec==0

drop isicr31 isicr4 isicr3 hs92 wiod1 wiod2

save ISIC3_`t'_wiodsample, replace
rename sec sec_old
merge m:1 sec_old using ../rawdata/sector_aggregation
drop if _m==2
drop _m

***collapse tariffs over WIOD sectors
collapse (mean) tariff, by(iso1 iso2 year sec_new)
*internal tariffs are zero
replace tariff=0 if iso1==iso2
*services sector tariffs are zero
replace tariff=0 if sec_new>22
save fin_`t'_wiodsample, replace

}

use fin_2000_wiodsample, clear
forvalues t=2001/2014{
append using fin_`t'_wiodsample
}
saveold ../rawdata/alltariffs_wiodsample, replace


*** aggregate tariffs for goods
use alltariffs_wiodsample, clear
keep if sec_new<23
collapse (mean) tariff, by(iso1 iso2 year)
sort iso1 iso2 year
saveold ../rawdata/alltariffs_wiodsample_agg, replace



erase help.dta




























************************************************************************************************
**************** TARIFFS FOR SIMULATION ********************************************************
************************************************************************************************



clear

u ../rawdata/hs_complete_8914, clear // downloaded from http://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=37
rename t year
keep if year>=1995 & year<=2014
sort i
saveold baci_1995-2014, replace

***merge country codes
insheet using ../rawdata/country_code_baci92.csv, clear delim(,)
drop iso2 name_french
rename iso3 iso_o
rename name_english country_o
drop if iso_o=="" | country_o==""
drop if iso_o=="YEM" | iso_o=="RVN" | iso_o=="SCG"
duplicates drop
sort i
merge 1:m i using baci_1995-2014
drop _m
*sort j
saveold baci_1995-2014, replace

insheet using ../rawdata/country_code_baci92.csv, clear delim(,)
drop iso2 name_french
rename iso3 iso_d
rename name_english country_d
drop if iso_d=="" | country_d==""
drop if iso_d=="YEM" | iso_d=="RVN" | iso_d=="SCG"
duplicates drop
rename i j
sort j
merge 1:m j using baci_1995-2014
drop _m

replace iso_o="TMP" if iso_o=="TLS"
replace iso_d="TMP" if iso_d=="TLS"
drop if iso_o=="WLD" | iso_d=="WLD"

drop if iso_o=="" | iso_d==""
replace iso_o="MKD" if iso_o=="?" & country_o=="The former Yugoslav Rep. of Macedonia"
replace iso_d="MKD" if iso_d=="?" & country_d=="The former Yugoslav Rep. of Macedonia"
destring hs6, replace


drop country*
*sort iso_o iso_d hs6 year
save baci_1995-2014, replace

egen imp_rev0_t=total(imports), by(iso_o iso_d year hs6) missing
collapse (first) imp_rev0_t, by(iso_o iso_d year hs6) fast

rename imp_rev0_t imp
label var imp "Imports (fob) in USD from BACI"
save baci_imp_9514_all_to_all_H0, replace

*split to years
forvalues x=1995(1)2014{
use baci_imp_9514_all_to_all_H0, clear
keep if year==`x'
replace iso_o="ROU" if iso_o=="ROM"
replace iso_d="ROU" if iso_d=="ROM"
saveold baci_imp_H0_`x', replace
}



forvalues x=1995(1)2014{
use "prel1_`x'", clear
rename hs92 hs6
rename iso2 iso_d
rename iso1 iso_o
merge 1:m iso_o iso_d hs6 year using baci_imp_H0_`x'
drop _m
compress
save "WITS_tariffs_imp_H0_`x'", replace
}





**** generate WIOD sample for 1995-2014
forvalues t = 1995/1997 {
use ../rawdata/countrysample_wiod.dta, clear
rename iso iso_o
*drop if iso_o=="ROW"
merge 1:m iso_o using WITS_tariffs_imp_H0_`t'
*drop if _m==2
replace iso_o="ROW" if _merge==1
g wiod1=0
replace wiod1=1 if _m==3
drop _m eu1 eu2 wto1 wto2
sort iso_d
compress
save wiod_tariff_weighted`t', replace
}

**

forvalues t = 1995/2014{
use ../rawdata/countrysample_wiod.dta, clear
rename iso iso_d
*drop if iso2=="ROW"
merge 1:m iso_d using wiod_tariff_weighted`t'.dta
*drop if _m==2 
replace iso_d="ROW" if _merge==1
g wiod2=0
replace wiod2=1 if _m==3
keep if wiod1==1 & wiod2==1
drop _merge
rename iso_d iso_d_h
drop if iso_o==iso_d_h & iso_o=="RoW"
count if iso_o~="RoW" & iso_d~="RoW" 
rename iso_o iso1
rename iso_d iso2
replace tariff=0 if pref==0
sort iso2 iso1 hs6 year
order iso* year hs6
save WITS_tariffs_H0_`t'_weighted, replace
erase wiod_tariff_weighted`t'.dta
}



* translate H0 to ISIC3
forvalues t=1995/2014{
insheet using "../rawdata/JobID-6_Concordance_H0_to_I3.csv", clear
drop hs198892productdescription isicrevision3productdescription
rename hs198892productcode hs6
rename isicrevision3productcode isicr3

merge 1:m hs6 using WITS_tariffs_H0_`t'_weighted
tab hs6 if _m==2
tab iso2 if _m==2
keep if _m==3
drop _m
save help, replace

* translate ISIC3 to ISIC3.1
insheet using "../rawdata/ISIC_Rev_3-ISIC_Rev_3_1_correspondence.txt", clear
drop partial3 partial31 activity
rename rev3 isicr3
rename rev31 isicr31
drop if isicr3=="n/a"
destring isicr3, replace
joinby isicr3 using help, unm(b)
tab _m
tab isicr3 if _m==2
tab iso2 if _m==2
keep if _m==3
drop _m
save help, replace

* translate ISIC3.1 to ISIC4
insheet using "../rawdata/ISIC31_ISIC4.txt", clear
drop partialisic31 partialisic4 detail
rename isic31code isicr31
rename isic4code isicr4
joinby isicr31 using help, unm(b)
tab _m
tab isicr31 if _m==2
tab iso2 if _m==2
keep if _m==3
drop _m
save help, replace

tostring isicr4, replace
replace isicr4 = "0" + isicr4 if length(isicr4)==3
g isic=substr(isicr4,1,2)

g sec=0

*Agriculture
replace sec=1 if isic=="01" 
replace sec=2 if isic=="02" 
replace sec=3 if isic=="03" 
*Mining and Quarrying
replace sec=4 if isic=="05" |  isic=="06" | isic=="07" | isic=="08" | isic=="09"
*Manufacturing
replace sec=5 if isic=="10" |  isic=="11" | isic=="12"
replace sec=6 if isic=="13" |  isic=="14" | isic=="15"
replace sec=7 if isic=="16" 
replace sec=8 if isic=="17" 
replace sec=9 if isic=="18" 
replace sec=10 if isic=="19" 
replace sec=11 if isic=="20" 
replace sec=12 if isic=="21" 
replace sec=13 if isic=="22" 
replace sec=14 if isic=="23" 
replace sec=15 if isic=="24" 
replace sec=16 if isic=="25" 
replace sec=17 if isic=="26" 
replace sec=18 if isic=="27" 
replace sec=19 if isic=="28" 
replace sec=20 if isic=="29" 
replace sec=21 if isic=="30" 
replace sec=22 if isic=="31" |  isic=="32"
replace sec=23 if isic=="33" 
*Services
replace sec=24 if isic=="35" 
replace sec=25 if isic=="36" 
replace sec=26 if isic=="37" |  isic=="38" |  isic=="39"
replace sec=27 if isic=="41" |  isic=="42" |  isic=="43"
replace sec=28 if isic=="45" 
replace sec=29 if isic=="46" 
replace sec=30 if isic=="47"
replace sec=31 if isic=="49" 
replace sec=32 if isic=="50" 
replace sec=33 if isic=="51" 
replace sec=34 if isic=="52"
replace sec=35 if isic=="53" 
replace sec=36 if isic=="55" |  isic=="56"
replace sec=37 if isic=="58" 
replace sec=38 if isic=="59" | isic=="60"
replace sec=39 if isic=="61" 
replace sec=40 if isic=="62" | isic=="63"
replace sec=41 if isic=="64" 
replace sec=42 if isic=="65" 
replace sec=43 if isic=="66" 
replace sec=44 if isic=="68" 
replace sec=45 if isic=="69" | isic=="70"
replace sec=46 if isic=="71" 
replace sec=47 if isic=="72" 
replace sec=48 if isic=="73" 
replace sec=49 if isic=="74" | isic=="75"
replace sec=50 if isic=="77" | isic=="78" | isic=="79" | isic=="80" | isic=="81" | isic=="82"
replace sec=51 if isic=="84" 
replace sec=52 if isic=="85" 
replace sec=53 if isic=="86" | isic=="87" | isic=="88"
replace sec=54 if isic=="90" | isic=="91" | isic=="92" | isic=="93" | isic=="94" | isic=="95" | isic=="96"
replace sec=55 if isic=="97" | isic=="98" 
replace sec=56 if isic=="99" 

*** no tau in services
tab isicr4 if sec==0

*drop isicr31 isicr4 isicr3 hs6 wiod1 wiod2

save ISIC3_`t'_weighted_wiodsample, replace
}





**************************************
*actual tariffs for the baseline *****
**************************************



*merge new aggregated sector classification
use ISIC3_2014_weighted_wiodsample, clear
rename sec sec_old
merge m:1 sec_old using ../rawdata/sector_aggregation
drop if _m==2
drop _m

*merge m:1 isic_2dig using sectors2



*internal tariffs are zero
replace tariff=0 if iso1==iso2
*services sector tariffs are zero
replace tariff=0 if sec_old>23

* set intra-EU tariffs to zero
g eu_o=0
g eu_d=0
local EUN2013 "AUT DEU BEL NLD LUX FRA ITA GBR IRL SWE DNK FIN ESP PRT GRC SVN SVK POL CZE HUN EST LVA LTU MLT CYP BGR ROM HRV"

foreach i of local EUN2013{
replace eu_o=1 if iso1=="`i'" 
replace eu_d=1 if iso2=="`i'" 
}

g eu=0
replace eu=1 if eu_o==1 & eu_d==1
replace tariff=0 if eu==1

drop eu*

g wtau=tariff*imp
egen twtau=total(wtau), by(iso1 iso2 year sec_new) missing
egen timp=total(imp), by(iso1 iso2 year sec_new) missing

g tau_weighted=twtau/timp

* simple average irrespective of whether a product is traded or not
egen tau_simple=mean(tariff), by(iso1 iso2 year sec_new) 

replace tau_weighted=0 if iso1==iso2 & iso1~="RoW"
sort iso2 iso1 sec_new year
replace timp=0 if iso1==iso2 & iso1~="RoW"

***collapse tariffs over WIOD sectors
collapse (first) tau_* (mean) tariff timp, by(iso1 iso2 sec_new year) fast

rename sec_new sector
order iso* year sec
rename timp imp
duplicates drop
sort iso2 iso1 year sector

count if tau_w==.
count if tau_w==. & tau_s==.
sum tau* tariff

replace tau_weighted=tau_simple if tau_w==.
save ../rawdata/appl_tau_2014_weighted, replace







*******************************
********* mfn tariffs *********
*******************************








*merge new aggregated sector classification
use ISIC3_2014_weighted_wiodsample, clear
rename sec sec_old
merge m:1 sec_old using ../rawdata/sector_aggregation
drop if _m==2
drop _m


*mfn
replace tariff=mfn1

*internal tariffs are zero
replace tariff=0 if iso1==iso2
*services sector tariffs are zero
replace tariff=0 if sec_old>23

g wtau=tariff*imp
egen twtau=total(wtau), by(iso1 iso2 year sec_new) missing
egen timp=total(imp), by(iso1 iso2 year sec_new) missing

g tau_weighted=twtau/timp

* simple average irrespective of whether a product is traded or not
egen tau_simple=mean(tariff), by(iso1 iso2 year sec_new) 

replace tau_weighted=0 if iso1==iso2 & iso1~="RoW"
sort iso2 iso1 sec_new year
replace timp=0 if iso1==iso2 & iso1~="RoW"

***collapse tariffs over WIOD sectors
collapse (first) tau_* (mean) tariff timp, by(iso1 iso2 sec_new year)

rename sec_new sector
order iso* year sec
rename timp imp
duplicates drop
sort iso2 iso1 year sector

count if tau_w==.
count if tau_w==. & tau_s==.
sum tau* tariff



replace tau_w=tau_w/100
replace tau_s=tau_s/100
replace tariff=tariff/100
drop imp

***fill the gaps of weighted tariffs! assuming same tariffs and same import volumes

sort iso2 iso1 sec year
replace tau_weighted=tau_simple if tau_w==.
sum if sec<23

keep iso2 iso1 sec year tau_w
rename tau_w tau_mfn_w
save ../rawdata/fin_cftau_2014_weighted_wiod,replace






*** erase data not needed

erase help.dta

forvalues t = 1995/2014 {
capture erase baci_imp_H0_`t'.dta
capture erase fin_tau_`t'_weighted.dta
capture erase ISIC3_`t'_weighted_wiodsample.dta
capture erase tariff`t'.dta
capture erase WITS_tariffs_H0_`t'_weighted.dta
capture erase WITS_tariffs_imp_H0_`t'.dta
}


forvalues t = 2000/2014{
capture erase fin_`t'_wiodsample.dta
capture erase ISIC3_`t'_wiodsample.dta
capture erase wiod_tariff`t'.dta
}



erase baci_1995-2014.dta
erase baci_imp_9514_all_to_all_H0.dta












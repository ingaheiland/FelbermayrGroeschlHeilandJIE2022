
******************** Tariff and FTA information for Schengen scenarios ************************

* no tariff changes for Schengen
* Schengen, EURO, EU, RTA info


cd "Repository/data"



u ../estimation/eucol_nnew_dataset_sec_IH, clear

keep if year==2014
drop imp  lta  sec

rename iso2 n
rename iso1 i

duplicates drop n i , force

save components, replace



******************************** pre - brexit scenarios

u WIOD_exp2014, clear
rename iso_supply n
rename iso_demand i

merge m:1 i n sector using tau_ni_2014_2digwiod
tab sector if _merge==1
drop _merge
replace tariff=0 if sec>"C22"
g tau_nik=(1+tariff)
label var tau "tariff (1+t/100) on imports of i from n in sector k"
drop tar

merge m:1 i n using components
*the _merge==1 cases all involve row
drop _merge

local vars "intst bothEEA botheuro eukor rta_c"
foreach x of local vars{
replace `x'=0 if i=="ROW"
replace `x'=0 if n=="ROW"
}
rename i iso_d
rename n iso_o

sort iso_d iso_o sector
tab iso_o if intst>0

g counter_Schengen = 0
g counter_Smarket = 0
g counter_Euro = 0
g counter_eukor =0

* undo the EU's trade agreements with outsiders
g counter_oRTA = rta_c
levelsof iso_d if bothEEA==1 & iso_d!="NOR" & iso_d!="CHE", local(EU14) clean
global EU14 "`EU14'"
di "$EU14"

foreach x in $EU14 {
replace counter_oRTA=0 if iso_d=="`x'" | iso_o=="`x'"
}

g diff_Schengen=counter_Schengen-intst
g diff_Smarket=counter_Smarket-bothEEA
g diff_Euro=counter_Euro-botheuro
g diff_oRTA=counter_oRTA-rta_c
g diff_eukor = counter_eukor-eukor

save counterfactuals, replace

u counterfactuals, clear
keep iso* sect diff*
reshape long diff_, i(iso_o sector iso_d) j(scenario) string
tab scenario
rename diff_ diff
reshape wide diff, i(scenario iso_o sector) j(iso_d) string /* note that this dataset should be sorted alike the trade share matrix (and tariffs), but b.c. dummies are "symmetric" the wrong sorting here doesn't matter*/
sort scenario sector iso_o

tab scenario

local scenarios "Schengen Smarket Euro oRTA eukor" 
foreach x of local scenarios{
outsheet diff* using "../simulation/input/dummy_diff_`x'.txt" if scenario=="`x'", replace nonames delimiter(",")
}




******************************** POST - Brexit scenarios

u WIOD_exp2014, clear
rename iso_supply n
rename iso_demand i

merge m:1 i n sector using tau_ni_2014_2digwiod
tab sector if _merge==1
drop _merge
replace tariff=0 if sec>"C22"
g tau_nik=(1+tariff)
label var tau "tariff (1+t/100) on imports of i from n in sector k"
drop tar

merge m:1 i n using components
*the _merge==1 cases all involve row
drop _merge

local vars "intst bothEEA botheuro eukor rta_c"
foreach x of local vars{
replace `x'=0 if i=="ROW"
replace `x'=0 if n=="ROW"
}

rename i iso_d
rename n iso_o

sort iso_d iso_o sector
tab iso_o if intst>0

g counter_Schengen = 0
g counter_Smarket = bothEEA
g counter_Euro = 0
g counter_eukor = eukor
g counter_oRTA = rta_c


* drop GBR
replace bothEEA=0 if iso_o=="GBR" | iso_d=="GBR"
levelsof iso_d if bothEEA==1, local(EU14pb) clean
global EU14pb "`EU14pb'"
di "$EU14pb"

foreach x in $EU14pb {
replace counter_oRTA=0 if (iso_d=="`x'" | iso_o=="`x'") & iso_d!="NOR" & iso_d!="CHE" & iso_o!="NOR" & iso_o!="CHE"
replace counter_Smarket=0 if iso_d=="`x'" | iso_o=="`x'"
replace counter_eukor=0 if iso_d=="`x'" | iso_o=="`x'"
}


g diff_Schengen=counter_Schengen-intst
g diff_Smarket=counter_Smarket-bothEEA
g diff_Euro=counter_Euro-botheuro
g diff_oRTA=counter_oRTA-rta_c
g diff_eukor=counter_eukor-eukor


keep iso* sect diff*
reshape long diff_, i(iso_o sector iso_d) j(scenario) string
tab scenario
rename diff_ diff
reshape wide diff, i(scenario iso_o sector) j(iso_d) string
sort scenario sector iso_o

tab scenario

local scenarios "Schengen Smarket Euro oRTA eukor" 
foreach x of local scenarios{
outsheet diff* using "../simulation/input/pb_dummy_diff_`x'.txt" if scenario=="`x'", replace nonames delimiter(",")
}





***************************************************************************
**************    prepare aggregate coefficients    ***********************
******                as trade(production) weighted avg.      *************
*****			         1000 cleaned sectoral draws			***********
***************************************************************************

* global production by sector

u output2014, clear
collapse (sum) out, by(sec)
sort sec
drop sec
g sec = _n
save output2014_sec, replace



* global export by sector


local vars "ltau schengen botheuro botheea eukor bothrta"

foreach x of local vars{
    
	local x = "ltau"
	import delim ../simulation/input/Bdraws/sectoral/b`x'-sec-nnew.txt, varnames(1) clear
g draw=_n
reshape long bs, i(draw) j(sec) string
replace sec = substr(sec,3,.)
destring sec, replace
merge m:1 sec using output2014_sec
drop _merge
g sec_agg = 1 if sec<23
replace sec_agg=2 if sec>22
egen tout=total(output), by(sec_agg draw)
egen ctout=total(output*bs), by(sec_agg draw)
replace bs = ctout/tout
drop ctou tout sec_agg output
reshape wide bs, i(draw) j(sec)

outsheet b* using "../simulation/input/b`x'_sec-wagg.txt", replace nonames delimiter(",")
}

***  same for the sectoral means

local vars "ltau schengen botheuro botheea eukor bothrta"

foreach x of local vars{
import delim ../simulation/input/Bdraws/sectoral/b`x'-sec-nnew.txt, varnames(1) clear
g draw=_n
reshape long bs, i(draw) j(sec) string
replace sec = substr(sec,3,.)
destring sec, replace
merge m:1 sec using output2014_sec
drop _merge
g sec_agg = 1 if sec<23
replace sec_agg=2 if sec>22
egen tout=total(output), by(sec_agg draw)
egen ctout=total(output*bs), by(sec_agg draw)
replace bs = ctout/tout
collapse (mean) bs, by(sec)
outsheet b* using "../simulation/input/b`x'_smd-wagg.txt", replace nonames delimiter(",")
}



*************************************
**tariffs: intra EU == mfn **********
*************************************


u taumfn_ni_2014_2digwiod, clear
rename tariff mfn_ink
replace mfn=1+mfn
rename sec k

merge 1:1 i n k using tau2014
drop _merge
replace mfn=1 if mfn==. /*service sectors*/
g test=mfn/tau
g taup_ink = tau_ink
g sector = k
merge m:1 n i  using components, keepusing(bothEEA)
*the missings all involve ROW
drop _merge sec test

replace  taup=mfn if bothEEA==1
save tariffs_eu_ink, replace

drop both tau_in mfn
reshape wide taup_in, i(n k) j(i) string
sort k n
outsheet taup_in* using "../simulation/input/taup_eu2014.txt", replace nonames delimiter(",")



*************************************
**tariffs: otherRTA == mfn **********
*************************************


u taumfn_ni_2014_2digwiod, clear
rename tariff mfn_ink
replace mfn=1+mfn
rename sec k

merge 1:1 i n k using tau2014
drop _merge
replace mfn=1 if mfn==. /*service sectors*/

g taup_ink = tau_ink

g sector = k

merge m:1 n i using components, keepusing(eukor bothEEA rta_c)
*the missings all involve ROW
drop _merge sec 

levelsof n if bothEEA==1 & n!="NOR" & n!="CHE", local(EU14) clean
global EU14 "`EU14'"
di "$EU14"

foreach x in $EU14 {
replace  taup=mfn if ( (n == "`x'" & rta_c ==1) | (n == "`x'" & eukor ==1) | (i == "`x'" & eukor ==1) | (i == "`x'" & rta_c ==1) ) & taup<mfn
}

save tariffs_orta_ink, replace

drop both tau_in mfn rta eukor
reshape wide taup_in, i(n k) j(i) string
sort k n
outsheet taup_in* using "../simulation/input/taup_oRTA2014.txt", replace nonames delimiter(",")





*********************************************************
**tariffs: otherRTA == mfn  +  intra EU == mfn **********
*********************************************************


u tariffs_orta_ink, clear

replace taup=mfn if bothEEA==1

save tariffs_ortaeu_ink, replace

drop both tau_in mfn eukor rta_c
reshape wide taup_in, i(n k) j(i) string
sort k n
outsheet taup_in* using "../simulation/input/taup_allEU2014.txt", replace nonames delimiter(",")








****************************************************************************
****************** Post-Brexit benchmark 



u WIOD_exp2014, clear
rename iso_supply n
rename iso_demand i

merge m:1 i n sector using tau_ni_2014_2digwiod
tab sector if _merge==1
drop _merge
replace tariff=0 if sec>"C22"
g tau_nik=(1+tariff)
label var tau "tariff (1+t/100) on imports of i from n in sector k"
drop tar

merge m:1 i n using components
*the _merge==1 cases all involve row
drop _merge

local vars "intst bothEEA botheuro  rta_c eukor"
foreach x of local vars{
replace `x'=0 if i=="ROW"
replace `x'=0 if n=="ROW"
}

rename i iso_d
rename n iso_o

sort iso_d iso_o sector
tab iso_o if intst>0

g counter_Schengen = intst
g counter_Smarket = bothEEA
replace counter_Smarket = 0 if iso_o=="GBR" | iso_d=="GBR"
g counter_Euro = botheuro
g counter_oRTA = rta_c
replace counter_oRTA=0 if iso_d=="GBR" | iso_o=="GBR"
g counter_eukor = eukor
replace counter_eukor=0 if iso_d=="GBR" | iso_o=="GBR"


g diff_Schengen=counter_Schengen-intst
g diff_Smarket=counter_Smarket-bothEEA
g diff_Euro=counter_Euro-botheuro
g diff_oRTA=counter_oRTA-rta_c
g diff_eukor=counter_eukor-eukor


keep iso* sect diff*
reshape long diff_, i(iso_o sector iso_d) j(scenario) string
tab scenario
rename diff_ diff
reshape wide diff, i(scenario iso_o sector) j(iso_d) string
sort scenario sector iso_o

tab scenario

local scenarios "Schengen Smarket Euro oRTA eukor" 
foreach x of local scenarios{
outsheet diff* using "../simulation/input/brexit_dummy_diff_`x'.txt" if scenario=="`x'", replace nonames delimiter(",")
}


*** tariffs: give UK mfn 


u taumfn_ni_2014_2digwiod, clear
rename tariff mfn_ink
replace mfn=1+mfn
rename sec k

merge 1:1 i n k using tau2014
drop _merge
replace mfn=1 if mfn==. /*service sectors*/

g taup_ink = tau_ink

g sector = k

merge m:1 n i using components, keepusing(bothEEA rta_c eukor)
*the missings all involve ROW
drop _merge sec 

rename n iso_d
rename i iso_o

levelsof iso_d if bothEEA==1, local(EU14) clean
global EU14 "`EU14'"
di "$EU14"

foreach x in $EU14 {
replace  taup=mfn if iso_o=="GBR" & iso_d=="`x'"
replace  taup=mfn if iso_d=="GBR" & iso_o=="`x'"
}

replace taup=mfn if iso_o=="GBR" & rta_c==1 /*set also tariffs of EU's RTA partners to mfn*/
replace taup=mfn if iso_d=="GBR" & rta_c==1

replace taup=mfn if iso_o=="GBR" & eukor==1 /*set also tariffs of EU's RTA partners to mfn*/
replace taup=mfn if iso_d=="GBR" & eukor==1



drop both tau_in mfn rta eukor
reshape wide taup_in, i(iso_d k) j(iso_o) string
sort k iso_d
outsheet taup_in* using "../simulation/input/taup_brexit2014.txt", replace nonames delimiter(",")










********************************************************************************
********************** prepare transfer data *********************************
******************************************************************************

import excel using ../rawdata/transfers.xls, clear cellrange(B3) firstrow sheet("2014")
keep in 103

drop Total earmarke other nonEU EU28 AK AL C

destring BE-UK, replace force

foreach var of varlist BE-UK {
rename `var' s`var'
}
drop if EUR==""

reshape long s, i(EUR) j(iso) string

replace iso ="AUT" if iso=="AT"
replace iso ="BEL" if iso=="BE"
replace iso ="BGR" if iso=="BG"
replace iso ="CYP" if iso=="CY"
replace iso ="CZE" if iso=="CZ"
replace iso ="DEU" if iso=="DE"
replace iso ="DNK" if iso=="DK"
replace iso ="GRC" if iso=="EL"
replace iso ="EST" if iso=="EE"
replace iso ="ESP" if iso=="ES"
replace iso ="FIN" if iso=="FI"
replace iso ="FRA" if iso=="FR"
replace iso ="HRV" if iso=="HR"
replace iso ="HUN" if iso=="HU"
replace iso ="IRL" if iso=="IE"
replace iso ="ITA" if iso=="IT"
replace iso ="LTU" if iso=="LT"
replace iso ="LUX" if iso=="LU"
replace iso ="LVA" if iso=="LV"
replace iso ="MLT" if iso=="MT"
replace iso ="NLD" if iso=="NL"
replace iso ="POL" if iso=="PL"
replace iso ="PRT" if iso=="PT"
replace iso ="ROU" if iso=="RO"
replace iso ="SWE" if iso=="SE"
replace iso ="SVN" if iso=="SI"
replace iso ="SVK" if iso=="SK"
replace iso ="GBR" if iso=="UK"

rename s transfer
sort iso

g year = 2014
merge m:1 year using ../rawdata/xrateEUR
keep if _merge==3
drop _merge price year

replace transfer = transfer*exc
label var transfer "net EU transfer payments in mio USD"
drop exch
save EU_transfers, replace

u EU_transfers, clear
rename iso n
merge 1:1 n using "SI_n_2014.dta" 
drop _merge

merge 1:1 n using "S_n_2014.dta" 
drop _merge

merge 1:1 n using R_n_2014
drop _merge

rename n iso
g SI_n_prime = SI_n+ transfer + R /*substract transfers*/
replace SI_n_prime = SI_n if SI_n_prime==.

g S_n_prime = S_n + transfer+R  /*substract transfers*/
replace S_n_prime = S_n if S_n_prime==.


save SI_n_prime_EU, replace

sort iso
outsheet S_n_prime using "../simulation/input/S_n_prime.txt", replace nonames delimiter(",")







*************************************************
********** make fake gamma's ********************
*** model with same-sector intermediates ********


**********************************************
*****the GAMMAs*******************************
**********************************************

u io_data_2014, clear
count if gma==. & output_nj!=0 & beta_nj!=1
collapse (mean) gma_nkj beta_nj, by(n k j)
replace gma = 1 if j==k & beta!=1
replace gma=0 if j!=k 
replace gma=0 if beta==1
drop beta
sort n j k
label var gma "cost share of intermediate from k in sector j's production in n"
reshape wide gma, i(n k) j(j) string
sort n k
outsheet gma* using "../simulation/input/G_fake.txt", delimiter(",") replace nonames




***********************************************************************
************** Input-output Matrix (J*N x N x J) **********************
***********************************************************************

* sorting: exporters in columns
* sorting rows: stack demanding sectors -> within: stack sectors -> within: sort importers


u io_data_2014, clear
replace beta_nj=1 if beta_nj==.
replace gma_nkj=0 if k!=j
replace gma_nkj=1 if k==j

*gen the new io coefficients
gen aa_inkj=(1-beta_nj)*gma_nkj*pi_ink/tau_ink
count if aa_inkj==.
keep i n k j aa_
label var aa_inkj "fob io coefficient i = iso_s, n = iso_d, j = sec_d, k = sec_s"

sort i k
reshape wide aa_inkj, i(i n k) j(j) string
reshape wide aa_inkj*, i(i k) j(n) string

drop k i
outsheet using "../simulation/input/ioc2014_fake.txt", replace nonames comma




*** this do-file reads in the world input-output table and computes the constant model parameters and baseline values for the simulation
** for the 50-sector version of our model

clear
set more off
set matsize 10000
set type double


cd "Repository/data"


*** calculation of 
*alpha_nj, beta_nj, gamma_nkj, pi_nkj, X_nk, Y_nk, I_n, VA_nk

*** moments matched in the data:
*bilateral sectoral fob trade flows  /*actually not perfectly matched bc negative gfcf is pushed to inventory*/
*tariffs
*fob intermediate usage per sector (J*K)

*** moments not matched:
*VA
*final goods vs intermediates shares in total trade
*gross production values

*** conditions observed:
*VA = output - intermediate use_cif (by sector)
*final goods expenditure = total expenditure - intermediate goods expenditure (by sector)
*income - net exports - net inventory consumption = final goods expenditure (by country)

*** critical assumption: all other trade cost are iceberg. Only then do basic prices (producer prices) incorporate the trade cost. 

u ../rawdata/WIOT2014_October16_ROW, clear // downloaded from https://www.rug.nl/ggdc/valuechain/wiod/
drop Industry*
rename Country iso_supply
rename RNr sec_supply
rename Year year
rename TOT BKG

save WIOT_STATA, replace

******original output for shares of merged sectors

preserve
keep iso_s sec_s BKG
rename iso iso 
rename sec sec
save output2014_orgsec, replace 
restore


******IO TAB*****
drop BKG
drop if sec_supply>56

recast long v*
reshape long v, i(year iso_s sec_s) j(dem) string

gen iso_demand=substr(dem,1,3)
gen sec_demand=substr(dem,-1,1) if length(dem)==4
replace sec_demand=substr(dem,-2,2) if length(dem)==5
destring sec_dem, replace
drop dem


* fix negative value addeds by merging sectors io structure
/*
BGR: average c10 w c14
MLT: average c03 w c10
LUX: average c31 w c32
*/

*average "columns"

g smu = 1 if iso_dem=="BGR" & (sec_dem==10 | sec_dem==14)
replace smu = 2 if iso_dem=="MLT" & (sec_dem==4 | sec_dem==10)
replace smu = 3 if iso_dem=="LUX" & (sec_dem==31 | sec_dem==32)

egen csum=total(v) if smu!=., by(iso* sec_s smu)

g sec = sec_d if smu!=.
g iso = iso_d if smu!=.

merge m:1 iso sec using output2014_orgsec
drop if _merge==2
drop _merge
rename BKG output
egen osum=total(output), by(iso* sec_s smu) missing
g oshare = output/osum

replace v = oshare*csum if smu!=.
drop sec iso output osum oshare smu csum


*average "rows"

g smu = 1 if iso_s=="BGR" & (sec_s==10 | sec_s==14)
replace smu = 2 if iso_s=="MLT" & (sec_s==4 | sec_s==10)
replace smu = 3 if iso_s=="LUX" & (sec_s==31 | sec_s==32)

egen rsum=total(v) if smu!=., by(iso* sec_d smu)

g sec = sec_s if smu!=.
g iso = iso_s if smu!=.

merge m:1 iso sec using output2014_orgsec
drop if _merge==2
drop _merge
rename BKG output
egen osum=total(output), by(iso* sec_d smu) missing
g oshare = output/osum

replace v = oshare*rsum if smu!=.
drop sec iso output osum oshare smu rsum


* collapse sectors with many zeros
replace sec_s=19 if sec_s==23
replace sec_s=42 if sec_s==43
replace sec_s=46 if sec_s==48 | sec_s==49
replace sec_s=54 if sec_s==55 | sec_s==56
sort sec_s
egen sec_sid=group(sec_s)


replace sec_d=19 if sec_d==23
replace sec_d=42 if sec_d==43
replace sec_d=46 if sec_d==48 | sec_d==49
replace sec_d=54 if sec_d==55 | sec_d==56
sort sec_d
egen sec_did=group(sec_d)

collapse (sum) v, by(iso_s iso_d sec_sid sec_did)

tostring sec_sid, g(sec_supply)
tostring sec_did, g(sec_demand)

drop sec_sid sec_did

replace sec_supply = "C0"+sec_supp if length(sec_s)==1
replace sec_demand = "C0"+sec_demand if length(sec_demand)==1
replace sec_supply = "C"+sec_supp if length(sec_s)==2
replace sec_demand = "C"+sec_demand if length(sec_demand)==2

rename v io

replace iso_d="ROU" if iso_d=="ROM"
replace iso_s="ROU" if iso_s=="ROM"

****** original IO TAB*****
preserve
drop if sec_demand>"C50"
drop if sec_supply>"C50"
replace io=0 if io==.
save io_orig_2014,replace
restore


**delete the inventories for now *note that this is no longer sector C61, but C55
drop if sec_demand=="C55"

count if io<0
tab sec_dem  if io<0
*some countries have negative capital formation -> I treat them alike inventories
*they are dropped here to assure that exports are always positive and will be substracted from inventory below
replace io=0 if io<0

preserve
drop if sec_demand>"C50"
*drop year
save io_tab_2014, replace
restore

*gen trade dataset
collapse (sum) io, by(iso_supply iso_demand sec_supply)
rename io exp
rename sec_su sector
label var exp "exports (final + intermediate) goods in mio current usd"
save WIOD_exp2014, replace




************ Output
u WIOT_STATA, clear

keep iso_s sec_s BKG
rename BKG output
gen year=2014
drop if sec_s>56

replace sec_s=19 if sec_s==23
replace sec_s=42 if sec_s==43
replace sec_s=46 if sec_s==48 | sec_s==49
replace sec_s=54 if sec_s==55 | sec_s==56
sort sec_s
egen sec_sid=group(sec_s)

collapse (sum) outp, by(iso sec_sid year)

tostring sec_sid, g(sec_supply)
drop sec_sid
replace sec_supply = "C0"+sec_supp if length(sec_s)==1
replace sec_supply = "C"+sec_supp if length(sec_s)==2
*replace iso_suppl="R0W" if iso_suppl=="TWN"
*collapse (sum) output, by(iso* sec* year)
replace iso_s="ROU" if iso_s=="ROM"
rename iso iso
rename sec sec
label var output "output by sector country in mio current usd"

save output2014, replace


// generate a compatible dataset with bilateral tariffs, including ROW
u ../rawdata/appl_tau_2014_weighted, clear

drop year tau_s tariff imp 
rename iso1 i /*iso demand*/
rename iso2 n /*iso_supply*/

drop if sec>23
tostring sec, replace
replace sec = "C0"+sec if length(sec)==1
replace sec = "C"+sec if length(sec)==2

renam tau_w tariff
replace tariff=tariff/100 // was in percent

save tau_norow, replace

u tau_norow, clear
rename i iso_demand
rename n iso_supply

global bricim "BRA IND RUS CHN IDN MEX"

g keepme=.
foreach x in $bricim {
replace keepme=1 if iso_d=="`x'"
}

keep if keepme==1

merge 1:1 iso* sec using WIOD_exp2014
keep  if _merge==3

drop keepme _merge
replace exp=0 if iso_s==iso_d /*dont use cty in row if iso_s is cty -> avoid downward bias for bricim countries */
egen timp=total(tariff*exp), by(sec iso_s)
egen ex = total(exp), by(sec iso_s)
g wtariff=timp/ex
replace wtariff=0 if ex==0
keep if iso_d=="BRA"
replace iso_d="ROW"
replace tariff = wtariff
drop ex timp wtar

save drowa, replace

*g ROW against ROW: using weighted average of bricim vs bricim tariffs
u drowa, clear

g keepme=.
foreach x in $bricim {
replace keepme=1 if iso_s=="`x'"
}
keep if keepme==1
egen timp=total(tariff*exp), by(sector iso_demand)
egen ex = total(exp), by(sec iso_d)
g wtariff=timp/ex
replace wtariff=0 if ex==0
keep if iso_s=="BRA"
replace iso_s="ROW"
replace tariff = wtariff
drop ex timp wtar keepme

append using drowa
save drow, replace


u tau_norow, clear
rename i iso_demand
rename n iso_supply

global bricim "BRA IND RUS CHN IDN MEX"

g keepme=.
foreach x in $bricim {
replace keepme=1 if iso_s=="`x'"
}

keep if keepme==1
merge 1:1 iso* sec using WIOD_exp2014
keep  if _merge==3

drop keepme _merge
replace exp=0 if iso_s==iso_d /**see comment above*/
egen timp=total(tariff*exp), by(sec iso_d)
egen ex = total(exp), by(sec iso_d)
g wtariff=timp/ex
replace wtariff=0 if ex==0
keep if iso_s=="BRA"
replace iso_s="ROW"
replace tariff = wtariff
drop ex timp wtar
save srow, replace

u srow,clear

append using drow

rename iso_d i
rename iso_s n

append using tau_norow

drop exp
save tau_ni_2014_2digwiod, replace

erase tau_norow.dta
erase srow.dta
erase drow.dta
erase drowa.dta


******combine IO and trade data
u WIOD_exp2014, clear
rename sector sec_supply
merge 1:m iso_su iso_dem sec_supply using io_tab_2014
drop _merge
rename exp Z_ink
label var Z_ink "fob exports from i to n in sector k in mio current usd"
rename iso_supply i
rename iso_demand n
rename sec_supply k
save io_trad_1, replace


u WIOD_exp2014, clear
rename iso_supply n
rename iso_demand i
rename sector k

merge 1:m i n k using io_trad_1
drop _merge
rename exp Z_nik
label var Z_nik "fob exports from n to i in sector k in mio current usd"


*****merge tariff data
rename k sector
merge m:1 i n sector using tau_ni_2014_2digwiod
tab sector if _merge==1
*only the service sectors have not been merged
*replace tau=0 if _merge==1
drop _merge
rename tariff tariff_ink

rename i n2
rename n i 
rename n2 n

merge m:1 i n sector using tau_ni_2014_2digwiod
*only the service sectors have not been merged
*replace tau=0 if _merge==1
drop _merge
rename tariff tariff_nik

rename i n2
rename n i 
rename n2 n
rename sector k

label var tariff_nik "tariff on imports of n from i in sector k"
label var tariff_ink "tariff on imports of i from n in sector k"

gen tau_ink=1+tariff_in
gen tau_nik=1+tariff_ni

replace tau_ink=1 if k>"C22"
replace tau_nik=1 if k>"C22"

*****cif expenditure and trade shares

gen X_ink=tau_ink*Z_ink
gen X_nik=tau_nik*Z_nik

rename sec_demand j

egen X_nk=total(X_ink), by(n j k)
egen X_ik=total(X_nik), by(i j k)
label var X_nk "total expenditure in n on goods from sector k"
label var X_ik "total expenditure in i on goods from sector k"

gen pi_ink=X_ink/X_nk
gen pi_nik=X_nik/X_ik

replace pi_ink=0 if pi_ink==.
replace pi_nik=0 if pi_nik==.


egen test=total(pi_ink), by(n k j)
summ test
drop test

rename io io_inkj
save io_exp_2014, replace


*merge output
u io_exp_2014, clear
rename n iso
rename j sec
merge m:1 iso sec using output2014
drop _merge
rename iso n
rename sec j
rename output output_nj
rename i iso
rename k sec
merge m:1 iso sec using output2014
drop _merge
rename iso i
rename sec k
rename output output_ik

*gen inventories
egen Ze_ik=total(Z_ink), by(i k j)
gen inventory_ik=output_ik-Ze_ik

preserve
collapse (mean) inventory, by(i k)
rename inventory inventory
rename i iso
rename k sec
save inventory_2014, replace
restore


*gen VA = output - intermediate inputs (cif)



*gen F_nk
egen F_nk=total(pi_ink/tau_ink), by(n k j)
label var F_nk "cif to fob correction for expenditure in n on goods from sector k"

*(1-beta)/beta*gma
egen tio_nkj_fob=total(io_inkj), by(n k j)
gen ombg_nkj=tio_nkj_fob/F_nk/output_nj
drop tio_nkj_fob


*tariff expenses for intermediate goods by demanding sectors
egen te_nj=total((tau_ink-1)/tau_ink*pi_ink*ombg_nkj*output_nj), by(n j)


**this is not the same as this
*egen te_nj2=total((tau_ink-1)*a_inkj*output_nj), by(n j)
**because the average pi's generate intermediate flows that are no longer consistent with the bilateralized io-coefficients

egen tio_nj=total(io_inkj), by(n j)

gen test=output_nj-tio_nj
summ test 
drop test

gen tio_cif_nj=tio_nj+te_nj
*output = va + cif intermediate usage
gen va_nj= output_nj-tio_cif_nj
summ va_nj

*the new beta
gen beta_nj=va_nj/output_nj
summ beta

*set mising betas to one
replace beta_nj=1 if output_nj==0

*the gamma 's
gen gma_nkj=ombg_nkj/(1-beta_nj)
sum gma
*set missing gmas to zero
replace gma_nkj=0 if beta_nj==1
egen test=total(gma_nkj), by(i n j)
summ test, detail
drop test

drop ombg tio* te_nj

*sectoral and overall trade balance
egen S_n_imp=total(Z_ink), by(n j)
egen S_n_exp=total(Z_nik), by(n j) 
gen S_n=S_n_exp-S_n_imp

egen S_i_imp=total(Z_nik), by(i j)
egen S_i_exp=total(Z_ink), by(i j) 
gen S_i=S_i_exp-S_i_imp
save io_trad_2014, replace


*gen dataset with va adjusted and beta
u io_trad_2014, clear
collapse (mean) beta_nj va_nj, by(n j)
rename beta beta 
rename va va
rename n iso
rename j sec
save va_2014, replace



*gen dataset with expenditure on intermediate goods by source sector
u io_trad_2014, clear
gen im_inkj=(1-beta_nj)*gma_nkj*output_nj
egen im_nk=total(im_inkj), by(n i k)
collapse (mean) im_nk, by(n k)
rename n iso
rename k sec
rename im im
label var im "cif expenditure for intermediates by source sector and demanding country in mio current usd"
save im2014, replace


***the alpha's
u io_trad_2014, clear
rename i iso
rename k sec
merge m:1 iso sec using im2014
drop _merge
rename im im_ik
rename iso i
rename sec k

******************** income
**tariff revenue
g r_nik=Z_nik*tariff_nik
egen R_i=total(r_nik), by(i j)
summ R_i
drop r_nik
**VA
rename i iso
rename k sec
merge m:1 iso sec using va_2014
drop _merge
rename va va_ik
rename iso i
rename sec k

egen VA_i=total(va_ik), by(n j i)



*inventory: net inventory consumption acts like an income transfer, similar to net exports

egen Inv_i=total(inventory_ik), by(n i j) 


**income
gen I_i=VA_i+R_i-S_i-Inv_i
********************************

*final goods expenditure
gen Xf_ik=X_ik-im_ik
egen Xf_i=total(Xf_ik), by(i n j)
gen test=Xf_i-I_i
summ test

gen alpha_ik=(X_ik-im)/I_i
summ alpha

egen test2=total(alph), by(i j n)
summ test2
drop test*

*u io_data_2014, clear

format Z* io_inkj-alpha %20.19f

save io_data_2014, replace




***************************************************************************************************************************************************************
********************************************DATA TRANSFER TO MATLAB********************************************************************************************
***************************************************************************************************************************************************************





************************
**trade shares**********
************************
u io_data_2014, clear
collapse (mean) pi_ink, by(i n k)
count if pi==.
label var pi_ink "country i's share in n's imports of goods from sector k"
save Din2014, replace
reshape wide pi_in, i(n k) j(i) string
sort k n
outsheet pi_in* using "../simulation/input/Din2014.txt", replace nonames delimiter(",")



************************
**tariffs**********
************************
u io_data_2014, clear
collapse (mean) tau_ink, by(i n k)
count if tau==.
label var tau_ink "n's tariffs on goods from sector k from country i"
save tau2014, replace
reshape wide tau_in, i(n k) j(i) string
sort k n
outsheet tau_in* using "../simulation/input/tariffs2014.txt", replace nonames delimiter(",")






**************************************************
**trade deficit INcluding inventory balance*******
**************************************************


u io_data_2014, clear
gen SI_i=S_i+Inv_i
collapse (mean) SI_i, by(i)
rename i n
count if SI==.
label var SI "Net exports of i + net inventory consumption (in million current US dollars)"
rename SI SI_n
save SI_n_2014, replace
sort n
outsheet SI_n using "../simulation/input/SI_n.txt", delimiter(",") replace nonames




**************************************************
**trade deficit EXcluding inventory balance*******
**************************************************

u io_data_2014, clear
collapse (mean) S_i, by(i)
rename i n
rename S_i S_n
count if S==.
label var S "Net exports of n (in million current US dollars)"
save S_n_2014, replace
sort n
outsheet S_n using "../simulation/input/S_n.txt", delimiter(",") replace nonames



****************************
*** Expenditure in 2000 ****
****************************

u io_data_2014, clear
collapse (mean) X_nk, by(n k)
count if X_nk==.
label var X_nk "expenditure in n on goods from sector k (in current mio usd)"
save X0_2014, replace
sort k n
reshape wide X_nk, i(k) j(n) string
sort k
outsheet X_n* using "../simulation/input/X0.txt", delimiter(",") replace nonames



********************************
****** Value added *************
********************************

u io_data_2014, clear
collapse (mean) VA_i, by(i)
count if VA_i==.
rename i n
rename VA VA_n
sort n
label var VA "total VA in country n (in mio current usd)"
save VA_n_2014, replace
sort n
drop n
outsheet VA using "../simulation/input/VA_n.txt", delimiter(",") replace nonames



********************************
****** Tariff revenue *************
********************************


u io_data_2014, clear
collapse (mean) R_i, by(i)
rename i n
rename R R_n
sort n
label var R "total tariff revenue in country n (in mio current usd)"
save R_n_2014, replace


********************************
****** GDP *********************
********************************

/*
u io_data_2014, clear
collapse (mean) output_nj, by(n j)
collapse (sum) output_nj, by(n)
outsheet output using "..\MATLAB\Y_n.txt", delimiter(",") replace nonames
*/


****************************************************
****** Global output by sector *********************
****************************************************

u io_data_2014, clear
collapse (sum) output_nj, by(j) fast
rename j sec
save output2014_sec,  replace




**********************************************
****** F: trade tariff weighted shares *******
**********************************************

u io_data_2014, clear
g f_nk = pi_ink/tau_ink
collapse (sum) f_nk, by(n j k)
collapse (mean) f_nk, by(n k)
label var f_nk "cif to fob correction for expenditure in n on goods from sector k"
save F_nk_2014, replace
reshape wide f_nk, i(k) j(n) string
sort k
outsheet f_nk* using "../simulation/input/F.txt", delimiter(",") replace nonames





**********************************************
****** tariff expenses for final goods********
**********************************************

/*
u io_data_2014, clear
g f_ik = pi_nik/tau_nik
collapse (sum) f_ik (mean) Xf_ik, by(i j k)
g tf_i=Xf_ik*(1-f_ik)
*there is a number of summ over pi's that is marginally greater (!) one. must be rounding errors
replace tf_i=0 if tf_i<0
collapse (mean) tf_i, by(i k)
collapse (sum) tf_i, by(i)
label var tf_i "tariff expenses for final goods consumption by demanding country"
rename i iso
save tf_i_2014, replace
*/



**********************************************
*****the GAMMAs*******************************
**********************************************

u io_data_2014, clear
count if gma==. & output_nj!=0 & beta_nj!=1
collapse (mean) gma_nkj, by(n k j)
replace gma=0 if gma==.
sort n j k
label var gma "cost share of intermediate from k in sector j's production in n"
save gma_nkj_2014, replace
reshape wide gma, i(n k) j(j) string
sort n k
outsheet gma* using "../simulation/input/G.txt", delimiter(",") replace nonames



**********************************************
****** Value added coefficient beta **********
**********************************************

u io_data_2014, clear
count if beta_nj==. & output_nj!=0
collapse (mean) beta_nj, by(n j)
replace beta=0 if beta==.
rename beta beta
rename j sector
rename n iso
label var beta "value added share in sectoral production (based on residual approach)"
save beta_2014, replace
sort sec iso
reshape wide beta, i(sec) j(iso) string
sort sec
outsheet beta* using "../simulation/input/B.txt", delimiter(",") replace nonames



******************************************
****** Expenditure shares alpha **********
******************************************


u io_data_2014, clear
collapse (mean) alpha_ik, by(i k)
count if alpha==.
rename alpha alpha
rename i iso
rename k sector
save alpha_2014, replace
sort sec iso
reshape wide alpha, i(sec) j(iso) string
sort sec
outsheet alpha* using "../simulation/input/alphas.txt", delimiter(",") replace nonames


*************************************************************************
********** Sectoral IM demand for net inventory**************************
**************demand of k-inputs for "inventory production" of all j's***
***********************same shape as alpha*******************************


u io_data_2014, clear
collapse (mean) beta_nj gma_nkj, by(n k j) 
rename n iso 
rename j sec
merge m:1 iso sec using inventory_2014
drop _merge
rename inventory inventory
rename sec j
rename iso n
gen VP_nkj=(1-beta_nj)*gma_nkj*inventory
replace VP_nkj=0 if VP_nkj==.
collapse (sum) VP_nkj, by(n k)
rename VP VP
rename n iso
rename k sector
save VP_2014, replace
sort sec iso
reshape wide VP, i(sec) j(iso) string
sort sec
outsheet VP* using "../simulation/input/VP.txt", delimiter(",") replace nonames


************************************************************************
***************sectoral inventories*************************************

u io_data_2014, clear
collapse (mean) inventory_ik, by(i k) 
rename k sec
rename i iso
sort sec iso
reshape wide inv, i(sec) j(iso) string
sort sec
outsheet inv* using "../simulation/input/inventory2014.txt", delimiter(",") replace nonames





***********************************************************************
************** Input-output Matrix (J*N x N x J) **********************
***********************************************************************

* sorting: exporters in columns
* sorting rows: stack demanding sectors -> within: stack sectors -> within: sort importers

u io_data_2014, clear
replace beta_nj=1 if beta_nj==.
replace gma_nkj=0 if gma_nkj==.
*gen the new io coefficients
gen aa_inkj=(1-beta_nj)*gma_nkj*pi_ink/tau_ink
count if aa_inkj==.
keep i n k j aa_
label var aa_inkj "fob io coefficient i = iso_s, n = iso_d, j = sec_d, k = sec_s"
save generic_iotab2014, replace



* alternative sorting
u generic_iotab2014, clear
sort k i
reshape wide aa_inkj, i(k i j) j(n) string
reshape wide aa_inkj*, i(k i) j(j) string
sort k i
drop k i
outsheet using ../simulation/input/ioc2014_2.txt, replace nonames comma







***************************************************************************************
************ send original io coeffients, corresponding va coeffients, and final goods*
* exports to MATLAB to compute original VA flows for the descriptives *****************


u io_orig_2014, clear
rename iso_d iso
rename sec_d sec

merge m:1 iso sec using output2014
drop _mer
rename iso iso_demand
rename sec sec_demand
egen IO = total(io), by(iso_d sec_d)
g beta = 1-IO/output
replace beta = 0 if beta==.
replace io = io/output
replace io = 0 if io==.

preserve
keep sec_d iso_d beta
collapse (first) beta, by(iso sec)
sort iso_d sec
export delim beta using ../simulation/input/beta4va.txt, replace novarname 
restore

keep iso* sec* io
reshape wide io, i(iso_s sec_s iso_d) j(sec_d) string
reshape wide io*, i(iso_s sec_s) j(iso_d) string
sort iso sec
export delim io* using ../simulation/input/io4va.txt, replace novarnames


u io_orig_2014, clear
collapse (sum) io, by(iso_s sec_s iso_d)
rename sec sector
merge 1:1 iso* sec using WIOD_exp2014
drop _m
g fexp=exp-io

keep iso* sec fexp
reshape wide fexp, i(iso_s sec) j(iso_d) string

export delim fexp* using ../simulation/input/fexp4va.txt, replace novarnames 









*******************************************************************************
******************** prep dataset with MFN tariffs ****************************
******************** for the counterfactual ***********************************


u ../rawdata/fin_cftau_2014_weighted_wiod, clear

keep if year==2014
drop year
rename iso1 i /*iso demand*/
rename iso2 n /*iso_supply*/

drop if sec>23

tostring sec, replace
replace sec = "C0"+sec if length(sec)==1
replace sec = "C"+sec if length(sec)==2

renam tau_mfn_w tariff


save taumfn_norow, replace

u taumfn_norow, clear
rename i iso_demand
rename n iso_supply

global bricim "BRA IND RUS CHN IDN MEX"

g keepme=.
foreach x in $bricim {
replace keepme=1 if iso_d=="`x'"
}

keep if keepme==1

merge 1:1 iso* sec using WIOD_exp2014
keep  if _merge==3

drop keepme _merge
replace exp=0 if iso_s==iso_d /*dont use Mex in row if iso_s is mex -> avoid downward bias for bricim countries */
egen timp=total(tariff*exp), by(sec iso_s)
egen ex = total(exp), by(sec iso_s)
g wtariff=timp/ex
replace wtariff=0 if ex==0
keep if iso_d=="BRA"
replace iso_d="ROW"
replace tariff = wtariff
drop ex timp wtar

save drowa_mfn, replace

*g ROW against ROW: using weighted average of bricim vs bricim tariffs
u drowa_mfn, clear

g keepme=.
foreach x in $bricim {
replace keepme=1 if iso_s=="`x'"
}
keep if keepme==1
egen timp=total(tariff*exp), by(sector iso_demand)
egen ex = total(exp), by(sec iso_d)
g wtariff=timp/ex
replace wtariff=0 if ex==0
keep if iso_s=="BRA"
replace iso_s="ROW"
replace tariff = wtariff
drop ex timp wtar keepme

append using drowa_mfn
save drow_mfn, replace


u taumfn_norow, clear
rename i iso_demand
rename n iso_supply

global bricim "BRA IND RUS CHN IDN MEX"

g keepme=.
foreach x in $bricim {
replace keepme=1 if iso_s=="`x'"
}

keep if keepme==1
merge 1:1 iso* sec using WIOD_exp2014
keep  if _merge==3

drop keepme _merge
replace exp=0 if iso_s==iso_d /**see comment above*/
egen timp=total(tariff*exp), by(sec iso_d)
egen ex = total(exp), by(sec iso_d)
g wtariff=timp/ex
replace wtariff=0 if ex==0
keep if iso_s=="BRA"
replace iso_s="ROW"
replace tariff = wtariff
drop ex timp wtar
save srow_mfn, replace

u srow_mfn,clear

append using drow_mfn
rename iso_d i
rename iso_s n
append using taumfn_norow
drop exp
save taumfn_ni_2014_2digwiod, replace



erase taumfn_norow.dta
erase drow_mfn.dta
erase srow_mfn.dta
erase drowa_mfn.dta

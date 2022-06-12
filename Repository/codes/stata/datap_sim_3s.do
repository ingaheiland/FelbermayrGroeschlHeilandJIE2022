*** this do-file reads in the world input-output table and computes the constant model parameters and baseline values for the simulation
*** for the 3-sector version of our model
*** as well as the input for the counterfactual scenarios




clear
set more off
set matsize 10000
*set processors 32
set type double

cd "Repository/data_3s"


******IO TAB*****

u ../data/WIOT_STATA, clear

drop BKG
drop if sec_supply>56

recast long v*
reshape long v, i(year iso_s sec_s) j(dem) string


gen iso_demand=substr(dem,1,3)
gen sec_demand=substr(dem,-1,1) if length(dem)==4
replace sec_demand=substr(dem,-2,2) if length(dem)==5
destring sec_dem, replace
drop dem



* collapse sectors: manufacturing 1-22 (include 23 in manuf to be consistent with the 50-sector analysis)

*tradable services: 28-49
* non-tradable services: rest

replace sec_s=1 if sec_s<24
replace sec_s=2 if sec_s>27 & sec_s<50
replace sec_s=3 if sec_s>23 & sec_s<28
replace sec_s=3 if sec_s>49 & sec_s<57

 
replace sec_d=1 if sec_d<24
replace sec_d=2 if sec_d>27 & sec_d<50
replace sec_d=3 if sec_d>23 & sec_d<28
replace sec_d=3 if sec_d>49 & sec_d<57


collapse (sum) v, by(iso_s iso_d sec_s sec_d)

tostring sec_s, replace
tostring sec_d, replace

replace sec_supply = "C0"+sec_supp if length(sec_s)==1
replace sec_demand = "C0"+sec_demand if length(sec_demand)==1

replace sec_supply = "C"+sec_supp if length(sec_s)==2
replace sec_demand = "C"+sec_demand if length(sec_demand)==2


rename v io

*replace iso_suppl="ROW" if iso_suppl=="TWN"
*replace iso_d="ROW" if iso_d=="TWN"
*collapse (sum) io, by(iso* sec* year)
replace iso_d="ROU" if iso_d=="ROM"
replace iso_s="ROU" if iso_s=="ROM"

****** original IO TAB*****
preserve
drop if sec_demand>"C03"
drop if sec_supply>"C03"
replace io=0 if io==.
save io_orig_2014,replace
restore

save io_interm, replace

u io_interm, clear
drop if sec_demand=="C61"

count if io<0
tab sec_dem  if io<0
*some countries have negative capital formation -> I treat them alike inventories
*they are dropped here to assure that exports are always positive and will be substracted from inventory below
replace io=0 if io<0

preserve
drop if sec_demand>"C03"
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
u ../data/WIOT_STATA, clear

keep iso_s sec_s BKG
rename BKG output
gen year=2014
drop if sec_s>56

replace sec_s=1 if sec_s<24
replace sec_s=2 if sec_s>27 & sec_s<50
replace sec_s=3 if sec_s>23 & sec_s<28
replace sec_s=3 if sec_s>49 & sec_s<57


collapse (sum) outp, by(iso sec_s year)

tostring sec_s, replace
replace sec_supply = "C0"+sec_supp if length(sec_s)==1
replace iso_s="ROU" if iso_s=="ROM"
rename iso iso
rename sec sec
label var output "output by sector country in mio current usd"

save output2014, replace



u ../data/tau_ni_2014_2digwiod, clear
rename i iso_supply
rename n iso_demand
merge 1:1 iso_demand iso_supply sector using ../data/WIOD_exp2014
keep if _merge==3
drop _merge
g timp = tar*exp
collapse (sum) timp exp, by(iso_s iso_d)
g tariff = timp/exp
drop timp exp
g sector = "C01"
rename iso_s i
rename iso_d n
save tau_ni_2014_2digwiod, replace


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

replace tau_ink=1 if k=="C02" | k=="C03"
replace tau_nik=1 if k=="C02" | k=="C03"

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
outsheet pi_in* using "../simulation/input_3s/Din2014.txt", replace nonames delimiter(",")



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
outsheet tau_in* using "../simulation/input_3s/tariffs2014.txt", replace nonames delimiter(",")






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
outsheet SI_n using "../simulation/input_3s/SI_n.txt", delimiter(",") replace nonames




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
outsheet S_n using "../simulation/input_3s/S_n.txt", delimiter(",") replace nonames



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
outsheet X_n* using "../simulation/input_3s/X0.txt", delimiter(",") replace nonames



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
outsheet VA using "../simulation/input_3s/VA_n.txt", delimiter(",") replace nonames





********************************
****** GDP *********************
********************************

/*
u io_data_2014, clear
collapse (mean) output_nj, by(n j)
collapse (sum) output_nj, by(n)
outsheet output using "Y_n.txt", delimiter(",") replace nonames
*/




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
outsheet f_nk* using "../simulation/input_3s/F.txt", delimiter(",") replace nonames





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
outsheet gma* using "../simulation/input_3s/G.txt", delimiter(",") replace nonames



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
outsheet beta* using "../simulation/input_3s/B.txt", delimiter(",") replace nonames



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
outsheet alpha* using "../simulation/input_3s/alphas.txt", delimiter(",") replace nonames


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
outsheet VP* using "../simulation/input_3s/VP.txt", delimiter(",") replace nonames


************************************************************************
***************sectoral inventories*************************************

u io_data_2014, clear
collapse (mean) inventory_ik, by(i k) 
rename k sec
rename i iso
sort sec iso
reshape wide inv, i(sec) j(iso) string
sort sec
outsheet inv* using "../simulation/input_3s/inventory2014.txt", delimiter(",") replace nonames


/*
************************************************************************
***************sectoral negative inventories****************************

u io_data_2014, clear
replace inventory_ik=0 if inventory_ik>0
collapse (mean) inventory_ik, by(i k) 
rename k sec
rename i iso
sort sec iso
reshape wide inv, i(sec) j(iso) string
sort sec
outsheet inv* using "inventory2014_neg.txt", delimiter(",") replace nonames



************************************************************************
***************sectoral positive inventories****************************

u io_data_2014, clear
replace inventory_ik=0 if inventory_ik<0
collapse (mean) inventory_ik, by(i k) 
rename k sec
rename i iso
sort sec iso
reshape wide inv, i(sec) j(iso) string
sort sec
outsheet inv* using "inventory2014_pos.txt", delimiter(",") replace nonames

*/







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


/*
*reshape wide aa_inkj, i(n j k) j(i) string
*sort j k n
*outsheet aa_* using "A.csv", delimiter(",") replace nonames
sort i k
reshape wide aa_inkj, i(i n k) j(j) string
reshape wide aa_inkj*, i(i k) j(n) string

drop k i
outsheet using ioc2014.txt, replace nonames comma
*/




u generic_iotab2014, clear
sort k i
reshape wide aa_inkj, i(k i j) j(n) string
reshape wide aa_inkj*, i(k i) j(j) string
sort k i
drop k i
outsheet using ../simulation/input_3s/ioc2014_2.txt, replace nonames comma




****************************************************************************************
********************** dataprep scenarios ********************************************
****************************************************************************************



******************** Tariff and FTA information for EU scenarios ************************



******************************** pre - brexit scenarios


u ../data/counterfactuals, clear
keep iso* sect diff*
g sec=1 if sect<"C23"
replace sec=2 if sect>"C26" & sect<"C46"
replace sec=3 if sect>"C22" & sect<"C27"
replace sec=3 if sect>"C45"

drop sect
rename sec sector
collapse (first) diff*, by(iso* sec)
reshape long diff_, i(iso_o sector iso_d) j(scenario) string
tab scenario
rename diff_ diff
reshape wide diff, i(scenario iso_o sector) j(iso_d) string /* note that this dataset should be sorted alike the trade share matrix (and tariffs), but b.c. dummies are "symmetric" the wrong sorting here doesn't matter*/
sort scenario sector iso_o

tab scenario

local scenarios "Schengen Smarket Euro oRTA eukor" 
foreach x of local scenarios{
outsheet diff* using "../simulation/input_3s/dummy_diff_`x'.txt" if scenario=="`x'", replace nonames delimiter(",")
}






***************************************************************************
**************    prepare aggregate coefficients    ***********************
******                as trade(production) weighted avg.      *************
*****			         1000 cleaned sectoral draws			***********
***************************************************************************

* global production by sector


local vars "ltau schengen botheuro botheea eukor bothrta"

foreach x of local vars{
    
	*local x "ltau"
import delim ../simulation/input/BDraws/sectoral/b`x'-sec-nnew.txt, varnames(1) clear
g draw=_n
reshape long bs, i(draw) j(sec) string
replace sec = substr(sec,3,.)
destring sec, replace
merge m:1 sec using ../data/output2014_sec
drop _merge
g sec_agg = 1 if sec<23
replace sec_a=2 if sec>26 & sec<46
replace sec_a=3 if sec>22 & sec<27
replace sec_a=3 if sec>45
egen tout=total(output), by(sec_agg draw)
egen ctout=total(output*bs), by(sec_agg draw)
replace bs = ctout/tout
collapse (mean) bs, by(sec_agg draw)
replace bs=0 if sec_agg==3 & "`x'"!="ltau"
reshape wide bs, i(draw) j(sec)
outsheet b* using "../simulation/input_3s/b`x'_sec-wagg.txt", replace nonames delimiter(",")
}

***  same for the sectoral means

local vars "ltau schengen botheuro botheea eukor bothrta"

foreach x of local vars{
import delim ../simulation/input/Bdraws/sectoral/b`x'-sec-nnew.txt, varnames(1) clear
g draw=_n
reshape long bs, i(draw) j(sec) string
replace sec = substr(sec,3,.)
destring sec, replace
merge m:1 sec using ../data/output2014_sec
drop _merge
g sec_agg = 1 if sec<23
replace sec_a=2 if sec>26 & sec<46
replace sec_a=3 if sec>22 & sec<27
replace sec_a=3 if sec>45
egen tout=total(output), by(sec_agg draw)
egen ctout=total(output*bs), by(sec_agg draw)
replace bs = ctout/tout
collapse (mean) bs, by(sec_agg)
replace bs=0 if sec_agg==3 & "`x'"!="ltau"
outsheet b* using "../simulation/input_3s/b`x'_smd-wagg.txt", replace nonames delimiter(",")
}


************ make similar dataset with aggregate estimates

* 1000 draws
u ../estimation/ddrawsMANUSERV, clear
keep line *_2
rename bltau bltau_3
foreach x in schengen botheuro bothEEA eukor bothrta {
g b`x'_3 = 0
drop b`x'_2
}
merge 1:1 line using ../estimation/ddrawsMANUSERV
drop _merge
rename line rep
aorder
rename bbothEEA* bbotheea*
local vars "ltau schengen botheuro botheea eukor bothrta"

foreach x of local vars{
outsheet b`x'* using "../simulation/input_3s/b`x'_sec-aest.txt", replace nonames delimiter(",")
}

collapse (mean) b*
g dum=1
reshape long bltau_ bschengen_ bbotheuro_ bbotheea_ beukor_ bbothrta_, i(dum) j(sec_agg)

foreach x of local vars{
outsheet b`x'* using "../simulation/input_3s/b`x'_smd-aest.txt", replace nonames delimiter(",")
}














*************************************
**tariffs: intra EU == mfn **********
*************************************


u ../data/taumfn_ni_2014_2digwiod, clear
rename i iso_supply
rename n iso_demand
merge 1:1 iso_demand iso_supply sector using ../data/WIOD_exp2014
keep if _merge==3
drop _merge
g timp = tar*exp
collapse (sum) timp exp, by(iso_s iso_d)
g tariff = timp/exp
drop timp exp
g sector = "C01"
rename iso_s i
rename iso_d n
save taumfn_ni_2014_2digwiod, replace


rename tariff mfn_ink
replace mfn=1+mfn
rename sec k

merge 1:1 i n k using tau2014
drop _merge
replace mfn=1 if mfn==. /*service sectors*/

g test=mfn/tau
count if test<.99 /* for a small number of obs mfn tariffs are smaller - maybe due to filling of different missings*/

g taup_ink = tau_ink

g sector = k

merge m:1 n i  using ../data/components, keepusing(bothEEA)
*the missings all involve ROW
drop _merge sec test

replace  taup=mfn if bothEEA==1
save tariffs_eu_ink, replace

drop both tau_in mfn
reshape wide taup_in, i(n k) j(i) string
sort k n
outsheet taup_in* using "../simulation/input_3s/taup_eu2014.txt", replace nonames delimiter(",")



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

g test=mfn/tau
count if test<.99 /* for a small number of obs mfn tariffs are smaller - maybe due to filling of different missings*/

g taup_ink = tau_ink

g sector = k

merge m:1 n i using ../data/components, keepusing(eukor bothEEA rta_c)
*the missings all involve ROW
drop _merge sec test


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
outsheet taup_in* using "../simulation/input_3s/taup_oRTA2014.txt", replace nonames delimiter(",")





*********************************************************
**tariffs: otherRTA == mfn  +  intra EU == mfn **********
*********************************************************


u tariffs_orta_ink, clear

replace taup=mfn if bothEEA==1

save tariffs_orta_ink, replace

drop both tau_in mfn eukor rta_c
reshape wide taup_in, i(n k) j(i) string
sort k n
outsheet taup_in* using "../simulation/input_3s/taup_allEU2014.txt", replace nonames delimiter(",")





















*************************************************
********** make fake gamma's ********************
*** model with same-sector intermediates ********



**********************************************
*****the GAMMAs*******************************
**********************************************

u io_data_2014, clear
count if gma==. & output_nj!=0 & beta_nj!=1
collapse (mean) gma_nkj, by(n k j)
replace gma = 1 if j==k
replace gma=0 if j!=k
sort n j k
label var gma "cost share of intermediate from k in sector j's production in n"
reshape wide gma, i(n k) j(j) string
sort n k
outsheet gma* using "../simulation/input_3s/G_fake.txt", delimiter(",") replace nonames




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

*reshape wide aa_inkj, i(n j k) j(i) string
*sort j k n
*outsheet aa_* using "A.csv", delimiter(",") replace nonames
sort i k
reshape wide aa_inkj, i(i n k) j(j) string
reshape wide aa_inkj*, i(i k) j(n) string

drop k i
outsheet using ../simulation/input_3s/ioc2014_fake.txt, replace nonames comma




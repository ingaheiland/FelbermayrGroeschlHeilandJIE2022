clear all
set matsize 10000


global main "Repository"
global source "$main/rawdata"
global tables "$main/tables"
global figures "$main/images"
global results "$main/simulation/Results"

cd "$results"



********* Table 5: real income changes

global scen "allEUtransfer allEU Euro Cunion Smarket oRTA  Schengen"

foreach x in $scen {
global sc "`x'"
u ${sc}1000_secB/what, clear
replace what = (what-1)*100
foreach x in 1 25 5 95 975 995 {
    g p`x'=.
}
levelsof iso, local(isolist)
foreach i in `isolist' {
_pctile what if iso=="`i'", p(.5 2.5 5 95 97.5 99.5)
replace p1=r(r1) if iso=="`i'"
replace p25=r(r2) if iso=="`i'"
replace p5=r(r3) if iso=="`i'"
replace p95=r(r4) if iso=="`i'"
replace p975=r(r5) if iso=="`i'"
replace p995=r(r6) if iso=="`i'"
}
collapse (mean) p* what, by(iso)
rename p* whatp*
rename what whatm
reshape long what, i(iso) j(unit) string
save ${sc}_smdB/what_ci, replace

u ${sc}1000_secB/wr, clear
replace wr = (wr-1)*100
foreach x in 1 25 5 95 975 995 {
    g p`x'=.
}
levelsof iso, local(isolist)
foreach i in `isolist' {
_pctile wr if iso=="`i'", p(.5 2.5 5 95 97.5 99.5)
replace p1=r(r1) if iso=="`i'"
replace p25=r(r2) if iso=="`i'"
replace p5=r(r3) if iso=="`i'"
replace p95=r(r4) if iso=="`i'"
replace p975=r(r5) if iso=="`i'"
replace p995=r(r6) if iso=="`i'"
}
collapse (mean) p* wr, by(iso)
rename p* wrp*
rename wr wrmean
reshape long wr, i(iso) j(unit) string
save ${sc}_smdB/wr_ci, replace
}


u allEU_smdB/cf_what, clear
keep if iso=="" /*make an empty dataset*/
g scen=""
foreach x in $scen {
global sc "`x'"
append using ${sc}_smdB/cf_what
append using ${sc}_smdB/what_ci
replace scen = "${sc}" if scen==""
}
drop ccode
replace what = (what-1)*100 if unit==""
replace unit="b" if unit==""
sort scen iso unit

save what_all, replace

* table for paper: real income changes

u what_all, clear
reshape wide what, i(scen iso) j(unit) string

*use mean of draws rather than sim at the mean
rename whatb whatatm
rename whatm whatb
tostring whatb, g(swhatb) format(%9.2f) force

foreach x in what {
replace s`x'b = s`x'b+"$^{*}$" if `x'b<0 & `x'p95<0 & `x'p975>0
replace s`x'b = s`x'b+"$^{**}$" if `x'b<0 & `x'p975<0 & `x'p995>0
replace s`x'b = s`x'b+"$^{***}$" if `x'b<0 & `x'p995<0
replace s`x'b = s`x'b+"$^{*}$" if `x'b>0 & `x'p5>0 & `x'p25<0
replace s`x'b = s`x'b+"$^{**}$" if `x'b>0 & `x'p25>0 & `x'p1<0
replace s`x'b = s`x'b+"$^{***}$" if `x'b>0 & `x'p1>0
}

keep iso scen s*b

reshape wide sw, i(iso) j(scen) string

global oldEU "AUT BEL DEU DNK ESP FIN FRA GBR GRC IRL ITA LUX NLD PRT SWE"
global newEU "BGR CYP CZE EST HRV HUN LTU LVA MLT POL ROU SVK SVN"

g isostar=iso
foreach x in $oldEU {
replace isostar = iso + "$^o$" if iso=="`x'" 
}
foreach x in $newEU {
replace isostar = iso + "$^n$" if iso=="`x'" 
}

listtex isostar *union *Smarket *Euro *Schengen *oRTA *allEU *transfer using "$tables/real_inc.tex", replace end(\\)

****************************************************************************************************************************
*************sectoral trade cost changes ***********************************************************************************
****************************************************************************************************************************

********** average schengen borders per group: S-S, S-n, n-n


u "../input/components", clear

*collapse sectors

collapse (first) intst, by(n i)

global schengen "AUS BEL CZE DNK EST FIN FRA DEU GRC HUN ISL ITA LVA LTU LUX MLT NLD NOR POL PRT SVN SVK ESP SWE CHE"

foreach x in $schengen {
replace i="S" if i=="`x'"
replace n="S" if n=="`x'"
}

replace i="N" if i!="S"
replace n="N" if n!="S"

collapse (mean) intst, by(i n)

*count average no of borders
*S-S: 2.7
*N-S: 1
*N-N: .2

global vars "ltau schengen botheuro botheea bothrta eukor"

foreach x in $vars{
import delim ../input/BDraws/sectoral/b`x'-sec-nnew.txt, varnames(1) clear
g rep=_n
if "`x'"!="ltau" {
merge 1:1 rep using coeffs_draws
drop _merge
}
save coeffs_draws, replace
}

rename bs1_* bltau*
rename bs2_* bbotheea*
rename bs3_* bschengen*
rename bs4_* beukor*
rename bs5_* bbothrta*
rename bs6_* bbotheuro*

save coeffs_draws, replace





/*
global labellist ""
forvalues i=1(1)38 {
global val = secshort[`i']
global labellist  $labellist  `i' "$val"
}
*/





*** btconfidence bounds

u coeffs_draws, clear
foreach x in botheuro botheea schengen {
forvalues i=1(1)50 {
gen ntb`x'`i'=(exp(b`x'`i'/bltau`i')-1)*100
}
}

forvalues i=1(1)50 {
gen ntbbothrta`i'=(exp(bbothrta`i'/bltau`i'+beukor`i'/bltau`i')-1)*100
}

forvalues i=1(1)50 {
gen ntball`i'=(exp(bschengen`i'/bltau`i'+bbotheuro`i'/bltau`i'+bbotheea`i'/bltau`i'+bbothrta`i'/bltau`i'+beukor`i'/bltau`i')-1)*100
}


forvalues i=1(1)50 {
gen ntbschengenss`i'=(exp(bschengen`i'/bltau`i'*2.7)-1)*100
gen ntbschengennn`i'=(exp(bschengen`i'/bltau`i'*.2)-1)*100
}

drop b*

reshape long ntbbotheuro ntbschengen ntbschengennn ntbschengenss ntbbotheea ntbbothrta ntball, i(rep) j(sec)

foreach var of varlist ntb* {
replace `var'=-`var'
g p95`var' = `var'
g p5`var' = `var'
}

collapse (p5) p5ntb*  (p95) p95ntb* (mean) ntb*, by(sec)

reshape long p5ntb p95ntb ntb, i(sec) j(scen) string



*cut bounds for display

replace p5 = -100 if p5<-100
replace p95 = 100 if p95>100


// Figure 3


twoway (rbar p5ntb p95ntb sec, color(gs12) barwidth(1))  (scatter ntb sec, color(black) msymbol(circle) msize()) if scen=="botheea",  legend(label(1 "90% CI") label(2 "Mean" )) ytitle("Change in NTBs (in %)") xla($labellist, angle(90)) scale(1.4)  xtitle("") graphregion(color(white))   xlabel(1(1)50) xtitle("Sector") yline(0)  xsize(10) ysize(3.5) ylabel(-100(50)100)
graph export "$figures/ntb_Smarket_bt.pdf", replace
graph export "$figures/fig3_1.pdf", replace
graph export "$figures/fig3_1.eps", replace



twoway (rbar p5ntb p95ntb sec, color(gs12) barwidth(1))  (scatter ntb sec, color(black) msymbol(circle) msize()) if scen=="schengen",  legend(label(1 "90% CI") label(2 "Mean" )) ytitle("Change in NTBs (in %)") xla($labellist, angle(90)) scale(1.4)  xtitle("") graphregion(color(white))   xlabel(1(1)50) xtitle("Sector") yline(0) xsize(10) ysize(3.5) ylabel(-100(50)100) legend(off)
graph export "$figures/ntb_Schengen_bt.pdf", replace
graph export "$figures/fig3_3.pdf", replace
graph export "$figures/fig3_3.eps", replace

twoway (rbar p5ntb p95ntb sec, color(gs12) barwidth(1))  (scatter ntb sec, color(black) msymbol(circle) msize()) if scen=="botheuro",  legend(label(1 "90% CI") label(2 "Mean" )) ytitle("Change in NTBs (in %)") xla($labellist, angle(90)) scale(1.4)  xtitle("") graphregion(color(white)) xlabel(1(1)50) xtitle("Sector") yline(0) ylabel(-100(50)100) xsize(10) ysize(3.5) legend(off)
graph export "$figures/ntb_Euro_bt.pdf", replace
graph export "$figures/fig3_2.pdf", replace
graph export "$figures/fig3_2.eps", replace

twoway (rbar p5ntb p95ntb sec, color(gs12) barwidth(1))  (scatter ntb sec, color(black) msymbol(circle) msize()) if scen=="bothrta",  legend(label(1 "90% CI") label(2 "Mean" )) ytitle("Change in NTBs (in %)") xla($labellist, angle(90)) scale(1.4)  xtitle("") graphregion(color(white)) xlabel(1(1)50) xtitle("Sector")  yline(0)  xsize(10) ysize(3.5) ylabel(-100(50)100) legend(off)
graph export "$figures/ntb_oRTA_bt.pdf", replace
graph export "$figures/fig3_4.pdf", replace
graph export "$figures/fig3_4.eps", replace

twoway (rbar p5ntb p95ntb sec, color(gs12) barwidth(1))  (scatter ntb sec, color(black) msymbol(circle) msize()) if scen=="all",  legend(label(1 "90% CI") label(2 "Mean" )) ytitle("Change in NTBs (in %)") xla($labellist, angle(90)) scale(1.4)  xtitle("") graphregion(color(white)) xlabel(1(1)50) xtitle("Sector")   yline(0) xsize(10) ysize(3.5) ylabel(-100(50)100)
graph export "$figures/ntb_allEU_bt.pdf", replace
graph export "$figures/fig3_5.pdf", replace
graph export "$figures/fig3_5.eps", replace

keep if scen=="all"
replace scen="Ball"

save scen4fig, replace



************** the same for estimates with trends


global vars "ltau schengen botheuro botheea bothrta eukor"

foreach x in $vars{
import delim ../input/BDraws/sectoral-trend/b`x'-sec-nnew.txt, varnames(1) clear
g rep=_n
if "`x'"!="ltau" {
merge 1:1 rep using coeffs_draws_trends
drop _merge
}
save coeffs_draws_trends, replace
}
rename bs1_* bltau*
rename bs2_* bbotheea*
rename bs3_* bschengen*
rename bs4_* beukor*
rename bs5_* bbothrta*
rename bs6_* bbotheuro*

save coeffs_draws_trends, replace



*** btconfidence bounds

u coeffs_draws_trends, clear
foreach x in botheuro botheea schengen {
forvalues i=1(1)50 {
gen ntb`x'`i'=(exp(b`x'`i'/bltau`i')-1)*100
}
}

forvalues i=1(1)50 {
gen ntbbothrta`i'=(exp(bbothrta`i'/bltau`i'+beukor`i'/bltau`i')-1)*100
}

forvalues i=1(1)50 {
gen ntball`i'=(exp(bschengen`i'/bltau`i'+bbotheuro`i'/bltau`i'+bbotheea`i'/bltau`i'+bbothrta`i'/bltau`i'+beukor`i'/bltau`i')-1)*100
}


forvalues i=1(1)50 {
gen ntbschengenss`i'=(exp(bschengen`i'/bltau`i'*2.7)-1)*100
gen ntbschengennn`i'=(exp(bschengen`i'/bltau`i'*.2)-1)*100
}

drop b*

reshape long ntbbotheuro ntbschengen ntbschengennn ntbschengenss ntbbotheea ntbbothrta ntball, i(rep) j(sec)

foreach var of varlist ntb* {
replace `var'=-`var'
g p95`var' = `var'
g p5`var' = `var'
}

collapse (p5) p5ntb*  (p95) p95ntb* (mean) ntb*, by(sec)

reshape long p5ntb p95ntb ntb, i(sec) j(scen) string



*cut bounds for display

replace p5 = -100 if p5<-100
replace p95 = 100 if p95>100
replace ntb=. if ntb>100 | ntb<-100

keep if scen=="all"
rename p5 tp5
rename p95 tp95
rename ntb tntb
append using scen4fig

* Figure 8
twoway (rbar p5ntb p95ntb sec  if scen=="Ball", color(gs10) barwidth(1))  (scatter ntb sec  if scen=="Ball", color(black) msymbol(circle) msize()) (rbar tp5 tp9 sec  if scen=="all", color(gs14%50) barwidth(1))  (scatter tntb sec  if scen=="all", color(orange) msymbol(circle) msize()),  legend(label(1 "Baseline: 90% CI") label(2 "Baseline: Mean" ) label(3 "Trends: 90% CI") label(4 "Trends: Mean")) ytitle("Change in NTBs (in %)") xla($labellist, angle(90)) scale(1.4)  xtitle("") graphregion(color(white))   xlabel(1(1)50) xtitle("Sector") yline(0)  xsize(10) ysize(3.5)
graph export "$figures/ntb_allEU_bt_both.pdf", replace
graph export "$figures/fig8.pdf", replace
graph export "$figures/fig8.eps", replace




********************************************************************************
************************* trade flows ******************************************

global sclist "Schengen allEU Euro Cunion Smarket oRTA allEUtransfer"

foreach sc in $sclist {
global scen ="`sc'"
cd "$results/${scen}1000_secB"
u rva, clear
merge 1:1 iso* rep sec using rvap
keep  if _merge==3
drop _merge
merge 1:1 iso* rep sec using exf
keep  if _merge==3
drop _merge
merge 1:1 iso* rep sec using exfp
keep  if _merge==3
drop _merge

rename rva va
rename rvap vap

global oldEU "AUT BEL DEU DNK ESP FIN FRA GBR GRC IRL ITA LUX NLD PRT SWE"
global newEU "BGR CYP CZE EST HRV HUN LTU LVA MLT POL ROU SVK SVN"

g reg_o=""
g reg_d=""
foreach x in $oldEU  {
replace reg_o="oldEU" if iso_o=="`x'"
replace reg_d="oldEU" if iso_d=="`x'"
}
foreach x in $newEU {
replace reg_o="newEU" if iso_o=="`x'"
replace reg_d="newEU" if iso_d=="`x'"
}
replace reg_o="nonEU" if reg_o==""
replace reg_d="nonEU" if reg_d==""

g EX=exf if iso_o!=iso_d
g EXp=exfp if iso_o!=iso_d
g VAX=va if iso_o!=iso_d
g VAXp=vap if iso_o!=iso_d

rename sec sec_id
merge m:1 sec_id using "$source/sectorlist"




*** for normalization in terms of US VA
egen vapUS = total(vap) if iso_o=="USA", by(rep)
egen vaUS = total(va) if iso_o=="USA", by(rep)

g hvaUS=vapUS/vaUS

collapse (sum) EX EXp VA* exf* va* (firstnm) hvaUS, by(reg* sec_agg rep) fast /*collapse over sectors and country groups*/

save EXVA_reps, replace


u "$results//${scen}_smdB/cf_rva", clear
merge 1:1 iso* sec using "$results//${scen}_smdB/cf_rvap"
keep if _merge==3
drop _merge
merge 1:1 iso* sec using "$results//${scen}_smdB/cf_exf"
keep if _merge==3
drop _merge
merge 1:1 iso* sec using "$results//${scen}_smdB/cf_exfp"
keep if _merge==3
drop _merge

rename rva vaa
rename rvap vaap

global oldEU "AUT BEL DEU DNK ESP FIN FRA GBR GRC IRL ITA LUX NLD PRT SWE"
global newEU "BGR CYP CZE EST HRV HUN LTU LVA MLT POL ROU SVK SVN"
g reg_o=""
g reg_d=""
foreach x in $oldEU {
replace reg_o="oldEU" if iso_o=="`x'"
replace reg_d="oldEU" if iso_d=="`x'"
}
foreach x in $newEU {
replace reg_o="newEU" if iso_o=="`x'"
replace reg_d="newEU" if iso_d=="`x'"
}
replace reg_o="nonEU" if reg_o==""
replace reg_d="nonEU" if reg_d==""

g EX=EXf if iso_o!=iso_d
g EXp=EXfp if iso_o!=iso_d
g VAX=vaa if iso_o!=iso_d
g VAXp=vaap if iso_o!=iso_d

rename sec sec_id
merge m:1 sec_id using "$source/sectorlist"


*** for normalization in terms of US VA
egen vapUS = total(vaap) if iso_o=="USA"
egen vaUS = total(vaa) if iso_o=="USA"

g hvaUS=vapUS/vaUS

collapse (sum) EX EXp VA* EXf* vaa* (firstnm) hvaUS, by(reg* sec_agg) fast /*collapse over sectors and country groups*/

save "$results//${scen}_smdB/EXVA", replace

}


* trade flows at the country group level


global sclist "Schengen allEU Euro Cunion Smarket oRTA allEUtransfer"

foreach sc in $sclist {
global scen ="`sc'"
cd "$results/${scen}1000_secB"

u EXVA_reps, clear

egen tva=total(va), by(reg_o rep)
egen tvap=total(vap), by(reg_o rep)
egen texf=total(exf), by(reg_o rep)
egen texfp=total(exfp), by(reg_o rep)

collapse (sum) EX* VA* (first) tv* tex* (firstnm) hvaUS, by(rep reg*) fast /*collapse over sectors and country groups*/

egen hVAUS = mean(hvaUS), by(rep)
drop hvaUS

g Yhat = (texfp/texf-hVAUS)*100
g VAhat = (tvap/tva-hVAUS)*100
drop tv* tex*

g Xhat=(EXp/EX-hVAUS)*100
g VAXhat=(VAXp/VAX-hVAUS)*100

*vax ratio changes
replace VAXhat = VAXhat-Xhat
replace VAhat = VAhat-Yhat


drop VAX VAXp EXp EX hVAUS

reshape wide *hat, i(rep reg_o) j(reg_d) string
drop Yhatold VAhatold Yhatnon VAhatnon

global valist ""
foreach var of varlist *EU {
foreach x in 1 25 5 95 975 995 {
    g `var'p`x'=.
}
levelsof reg_o, local(isolist)
global isolist "`isolist'"
foreach i in `isolist' {
_pctile `var' if reg=="`i'", p(.5 2.5 5 95 97.5 99.5)
replace `var'p1=r(r1) if reg=="`i'"
replace `var'p25=r(r2) if reg=="`i'"
replace `var'p5=r(r3) if reg=="`i'"
replace `var'p95=r(r4) if reg=="`i'"
replace `var'p975=r(r5) if reg=="`i'"
replace `var'p995=r(r6) if reg=="`i'"
}
global valist "$valist" "`var'"
rename `var' `var'b
}

collapse (mean) *hat*, by(reg)

g scenario="${scen}"

reshape long $valist, i(reg scen) j(unit) string

save ci_agg, replace



u "$results//${scen}_smdB/EXVA", clear

egen tva=total(vaa), by(reg_o )
egen tvap=total(vaap), by(reg_o )
egen texf=total(EXf), by(reg_o )
egen texfp=total(EXfp), by(reg_o)


collapse (sum) EX EXp VA* (first) tv* tex* (firstnm) hvaUS, by(reg*) fast /*collapse over sectors and country groups*/

egen hVAUS = mean(hvaUS)
drop hvaUS


g Yhat = (texfp/texf-hVAUS)*100
g VAhat = (tvap/tva-hVAUS)*100
drop tv* tex*

g Xhat=(EXp/EX-hVAUS)*100
g VAXhat=(VAXp/VAX-hVAUS)*100

*vax ratio changes
replace VAXhat = VAXhat-Xhat
replace VAhat = VAhat-Yhat

drop VAX VAXp EXp EX hVAUS

reshape wide *hat, i(reg_o) j(reg_d) string
drop Yhatold VAhatold Yhatnon VAhatnon


g scenario="${scen}"
g unit = "m"
append using ci_agg

save "$results//${scen}_smdB/tab_agg", replace
}



u "$results/Schengen_smdB/tab_agg", replace
keep if unit=="" /*make empty dataset*/

foreach sc in $sclist {
global scen ="`sc'"
append using "$results//${scen}_smdB/tab_agg"
}

reshape wide VAX* Xhat* Yhat VAhat, i(scen reg ) j(unit) string

tostring VAXhatnewEUb VAXhatnonEUb VAXhatoldEUb, g(sVAXhatnewEUb sVAXhatnonEUb sVAXhatoldEUb) format(%9.2f) force
tostring XhatnewEUb XhatnonEUb XhatoldEUb, g(sXhatnewEUb sXhatnonEUb sXhatoldEUb) format(%9.2f) force
tostring YhatnewEUb VAhatnewEUb, g(sYhatnewEUb sVAhatnewEUb) format(%9.2f) force

foreach x in VAhatnewEU YhatnewEU VAXhatnewEU XhatnewEU VAXhatoldEU XhatoldEU VAXhatnonEU XhatnonEU {
replace s`x'b = s`x'b+"$^{*}$" if `x'b<0 & `x'p95<0 & `x'p975>0
replace s`x'b = s`x'b+"$^{**}$" if `x'b<0 & `x'p975<0 & `x'p995>0
replace s`x'b = s`x'b+"$^{***}$" if `x'b<0 & `x'p995<0
replace s`x'b = s`x'b+"$^{*}$" if `x'b>0 & `x'p5>0 & `x'p25<0
replace s`x'b = s`x'b+"$^{**}$" if `x'b>0 & `x'p25>0 & `x'p1<0
replace s`x'b = s`x'b+"$^{***}$" if `x'b>0 & `x'p1>0
}


keep reg scen s*b

g sort = 1  if reg=="oldEU"
replace sort = 2 if reg=="newEU"
replace sort = 3 if reg=="nonEU"
sort scen sort

replace reg ="\quad old EU" if reg=="oldEU"
replace reg ="\quad new EU" if reg=="newEU"
replace reg ="\quad non-EU" if reg=="nonEU"

** Table 2

foreach x in $sclist {
listtex reg sYhat sVAhat sXhatold sVAXhatold sXhatnew sVAXhatnew sXhatnon sVAXhatnon using "$tables/aggexp`x'.tex" if scen=="`x'", replace end(\\)
}




********************************************
* trade flows at the sectoral and country-group level


foreach sc in $sclist {
global scen ="`sc'"
cd "$results//${scen}1000_secB"

u EXVA_reps, clear
drop va* exf*

egen hVAUS=mean(hvaUS), by(rep)
drop hva
drop if reg_o=="nonEU"


reshape wide EX EXp VAX VAXp, i(reg_o sec_agg rep hVAUS) j(reg_d) string

foreach x in EX EXp VAX VAXp {
g `x'EU = `x'newEU+`x'oldEU
g `x'WLD = `x'EU+`x'nonEU
}

drop  *oldEU *newEU 

g XhatWLD=(EXpWLD/EXWLD-hVAUS)*100
g XhatEU=(EXpEU/EXEU-hVAUS)*100
g XhatnonEU=(EXpnonEU/EXnonEU-hVAUS)*100

g VAXhatWLD=(VAXpWLD/VAXWLD-hVAUS)*100-XhatWLD
g VAXhatEU=(VAXpEU/VAXEU-hVAUS)*100-XhatEU
g VAXhatnonEU=(VAXpnonEU/VAXnonEU-hVAUS)*100-XhatnonEU

keep reg_o sec_agg Xhat* VAXha*

global valist ""
foreach var of varlist *EU *WLD {
foreach x in 1 25 5 95 975 995 {
    g `var'p`x'=.
}
levelsof reg_o, local(isolist)
global isolist "`isolist'"
foreach i in `isolist' {
levelsof sec_agg, local(seclist)
global seclist "`seclist'"
foreach j in `seclist' {
_pctile `var' if reg=="`i'" & sec=="`j'", p(.5 2.5 5 95 97.5 99.5)
replace `var'p1=r(r1) if reg=="`i'" & sec=="`j'"
replace `var'p25=r(r2) if reg=="`i'" & sec=="`j'"
replace `var'p5=r(r3) if reg=="`i'" & sec=="`j'"
replace `var'p95=r(r4) if reg=="`i'" & sec=="`j'"
replace `var'p975=r(r5) if reg=="`i'" & sec=="`j'"
replace `var'p995=r(r6) if reg=="`i'" & sec=="`j'"
}
}
global valist "$valist" "`var'"
rename `var' `var'b
}

collapse (mean) *hat*, by(reg sec)

g scenario="${scen}"

reshape long $valist, i(reg sec scen) j(unit) string

save ci_agg_sec, replace



u "$results//${scen}_smdB/EXVA", clear

drop vaa* EXf*
egen hVAUS = mean(hvaUS)
drop hva
reshape wide EX EXp VAX VAXp, i(reg_o sec_agg hVAUS) j(reg_d) string

foreach x in EX EXp VAX VAXp {
g `x'EU = `x'newEU+`x'oldEU
g `x'WLD = `x'EU+`x'nonEU
}

drop *oldEU *newEU 

g XhatWLD=(EXpWLD/EXWLD-hVAUS)*100
g XhatEU=(EXpEU/EXEU-hVAUS)*100
g XhatnonEU=(EXpnonEU/EXnonEU-hVAUS)*100

g VAXhatWLD=(VAXpWLD/VAXWLD-hVAUS)*100-XhatWLD
g VAXhatEU=(VAXpEU/VAXEU-hVAUS)*100-XhatEU
g VAXhatnonEU=(VAXpnonEU/VAXnonEU-hVAUS)*100-XhatnonEU


keep reg_o sec_agg Xhat* VAXha*

drop if reg_o=="nonEU"

g scenario="${scen}"
g unit = "m"
append using ci_agg_sec

replace sec_agg="Agric." if sec_agg=="Agriculture & Natural Resources"
replace sec_agg="Manuf." if sec_agg=="Manufacturing"
replace sec_agg="Serv." if sec_agg=="Services"
save "$results//${scen}_smdB/tab_agg_sec", replace
}


u "$results/Schengen_smdB/tab_agg_sec", replace
keep if unit=="" /*make empty dataset*/

foreach sc in $sclist {
global scen ="`sc'"
append using "$results//${scen}_smdB/tab_agg_sec"
}

reshape wide VAX* Xhat*, i(scen reg sec) j(unit) string

tostring VAXhatEUb VAXhatnonEUb VAXhatWLDb, g(sVAXhatEUb sVAXhatnonEUb sVAXhatWLDb) format(%9.2f) force
tostring XhatEUb XhatnonEUb XhatWLDb, g(sXhatEUb sXhatnonEUb sXhatWLDb) format(%9.2f) force

foreach x in VAXhatEU XhatEU VAXhatWLD XhatWLD VAXhatnonEU XhatnonEU {
replace s`x'b = s`x'b+"$^{*}$" if `x'b<0 & `x'p95<0 & `x'p975>0
replace s`x'b = s`x'b+"$^{**}$" if `x'b<0 & `x'p975<0 & `x'p995>0
replace s`x'b = s`x'b+"$^{***}$" if `x'b<0 & `x'p995<0
replace s`x'b = s`x'b+"$^{*}$" if `x'b>0 & `x'p5>0 & `x'p25<0
replace s`x'b = s`x'b+"$^{**}$" if `x'b>0 & `x'p25>0 & `x'p1<0
replace s`x'b = s`x'b+"$^{***}$" if `x'b>0 & `x'p1>0
}

keep reg scen s*b sec

g sort = 1  if reg=="oldEU"
replace sort = 2 if reg=="newEU"
sort scen sort sec


g reg="\quad old EU" if reg_o=="oldEU" & sec_agg=="Agric."
replace reg ="\quad new EU" if reg_o=="newEU" & sec_agg=="Agric."


*Table 3
*Table B7

foreach x in $sclist {
listtex reg sec sXhatEU sVAXhatEU sXhatnon sVAXhatnon sXhatWLD sVAXhatWLD using "$tables/aggsecexp`x'.tex" if scen=="`x'", replace end(\\)
}










/***** sectoral output growth and shares ********/

cd "$main"


* baseline values
u data/va_2014, clear
rename sec sector
merge m:1 sector using "$source/sectorlist"
drop _merge

g reg_o=""
foreach x in $oldEU {
replace reg_o="oldEU" if iso=="`x'"
}
foreach x in $newEU {
replace reg_o="newEU" if iso=="`x'"
}
replace reg_o="nonEU" if reg_o==""

collapse (sum) va, by(reg sec_ag)
egen tva=total(va), by(reg)
g sva = va/tva
replace va = va/1000 /*now in bn USD*/
replace sva = sva*100 /*now in %*/
format sva %9.1f
format va %9.0f
drop tva
replace sec_agg="Agric." if sec_agg=="Agriculture & Natural Resources"
replace sec_agg="Manuf." if sec_agg=="Manufacturing"
replace sec_agg="Serv." if sec_agg=="Services"
drop if reg_o=="nonEU"

save "$results/va_base", replace


cd "$results"

foreach sc in $sclist {
global scen ="`sc'"

u ${scen}1000_secB/EXVA_reps, clear

collapse (sum) va vap (firstnm) hvaUS, by(reg_o sec rep)
egen hVAUS = mean(hva), by(rep)
drop hvaUS
egen tva = total(va), by(reg rep)
egen tvap = total(vap), by(reg rep)
g dsva=(vap/tvap-va/tva)*100
g hva=(vap/va-hVAUS)*100

global valist ""
foreach var of varlist dsva hva {
foreach x in 1 25 5 95 975 995 {
    g `var'p`x'=.
}
levelsof reg_o, local(isolist)
global isolist "`isolist'"
foreach i in `isolist' {
levelsof sec_agg, local(seclist)
global seclist "`seclist'"
foreach j in `seclist' {
_pctile `var' if reg=="`i'" & sec=="`j'", p(.5 2.5 5 95 97.5 99.5)
replace `var'p1=r(r1) if reg=="`i'" & sec=="`j'"
replace `var'p25=r(r2) if reg=="`i'" & sec=="`j'"
replace `var'p5=r(r3) if reg=="`i'" & sec=="`j'"
replace `var'p95=r(r4) if reg=="`i'" & sec=="`j'"
replace `var'p975=r(r5) if reg=="`i'" & sec=="`j'"
replace `var'p995=r(r6) if reg=="`i'" & sec=="`j'"
}
}
global valist "$valist" "`var'"
rename `var' `var'b
}

collapse (mean) dsva* hva*, by(reg sec)

g scenario="${scen}"

reshape long $valist, i(reg sec scen) j(unit) string

save "${scen}_smdB/vahat_ci", replace


u ${scen}_smdB/EXVA, clear

collapse (sum) vaa vaap (firstnm) hvaUS, by(reg_o sec) fast
egen hVAUS = mean(hvaUS)
egen tva = total(vaa), by(reg)
egen tvap = total(vaap), by(reg)
g dsva=(vaap/tvap-vaa/tva)*100
g hva=(vaap/vaa-hVAUS)*100
keep reg sec ds hva
g scenario = "$scen"
g unit = "m"

append using  ${scen}_smdB/vahat_ci
sort reg_o sec unit

replace sec_agg="Agric." if sec_agg=="Agriculture & Natural Resources"
replace sec_agg="Manuf." if sec_agg=="Manufacturing"
replace sec_agg="Serv." if sec_agg=="Services"
save ${scen}_smdB/vahat, replace
}

* table for paper

u "$results/Schengen_smdB/vahat", replace
keep if unit=="" /*make empty dataset*/

foreach sc in $sclist {
global scen ="`sc'"
append using "$results//${scen}_smdB/vahat"
}
drop if reg_o=="nonEU"

reshape wide dsva hva, i(reg sec scen) j(unit) string

tostring dsvab, g(sdsvab) format(%9.2f) force
tostring hvab, g(shvab) format(%9.2f) force

foreach x in dsva hva {
replace s`x'b = s`x'b+"$^{*}$" if `x'b<0 & `x'p95<0 & `x'p975>0
replace s`x'b = s`x'b+"$^{**}$" if `x'b<0 & `x'p975<0 & `x'p995>0
replace s`x'b = s`x'b+"$^{***}$" if `x'b<0 & `x'p995<0
replace s`x'b = s`x'b+"$^{*}$" if `x'b>0 & `x'p5>0 & `x'p25<0
replace s`x'b = s`x'b+"$^{**}$" if `x'b>0 & `x'p25>0 & `x'p1<0
replace s`x'b = s`x'b+"$^{***}$" if `x'b>0 & `x'p1>0
}

rename sdsvab ds
rename shvab hv
keep reg scen sec hv ds
reshape wide hv ds, i(reg sec) j(scen) string

merge 1:1 reg sec using "$results/va_base"
drop _merge

g reg = reg_o if sec_agg=="Agric."
replace reg ="old EU" if reg=="oldEU"
replace reg ="new EU" if reg=="newEU"

** Table 4
listtex reg sec va hvCunion hvSmarket  hvEuro hvSchengen hvoRTA hvallEU hvallEUtran using "$tables/outputp_s.tex", replace end(\\)
listtex reg sec sva dsCunion dsSmarket dsEuro dsSchengen dsoRTA dsallEU dsallEUtran using "$tables/output_s.tex", replace end(\\)




* trade and value added at the country-sector level


*foreach sc in $sclist {
foreach sc in allEU {

global scen ="`sc'"
cd "$results//${scen}1000_secB"

u exf, clear
merge 1:1 iso* rep sec using exfp
keep  if _merge==3
drop _merge

collapse (sum) exf*, by(iso_o sec rep) /*collapse over sectors and country groups*/

global oldEU "AUT BEL DEU DNK ESP FIN FRA GBR GRC IRL ITA LUX NLD PRT SWE"
global newEU "BGR CYP CZE EST HRV HUN LTU LVA MLT POL ROU SVK SVN"

g reg_o=""
foreach x in $oldEU  {
replace reg_o="oldEU" if iso_o=="`x'"
}
foreach x in $newEU {
replace reg_o="newEU" if iso_o=="`x'"
}
replace reg_o="nonEU" if reg_o==""

rename sec sec_id
merge m:1 sec_id using "$source/sectorlist"

save EXiso_reps, replace


u "$results//${scen}_smdB/cf_EXf", clear
merge 1:1 iso* sec using "$results//${scen}_smdB/cf_EXfp"
keep if _merge==3
drop _merge

collapse (sum) EX*, by(iso_o sec)


global oldEU "AUT BEL DEU DNK ESP FIN FRA GBR GRC IRL ITA LUX NLD PRT SWE"
global newEU "BGR CYP CZE EST HRV HUN LTU LVA MLT POL ROU SVK SVN"
g reg_o=""
foreach x in $oldEU {
replace reg_o="oldEU" if iso_o=="`x'"
}
foreach x in $newEU {
replace reg_o="newEU" if iso_o=="`x'"
}
replace reg_o="nonEU" if reg_o==""


rename sec sec_id
merge m:1 sec_id using "$source/sectorlist"

drop _merge
save "$results//${scen}_smdB/EXiso", replace
}


foreach sc in $sclist {
global scen ="`sc'"
cd "$results//${scen}1000_secB"


u rva, clear
rename rva va
merge 1:1 iso* rep sec using rvap
keep  if _merge==3
drop _merge
rename rvap vap

collapse (sum) va*, by(iso_o sec rep) /*collapse over sectors and country groups*/

global oldEU "AUT BEL DEU DNK ESP FIN FRA GBR GRC IRL ITA LUX NLD PRT SWE"
global newEU "BGR CYP CZE EST HRV HUN LTU LVA MLT POL ROU SVK SVN"

g reg_o=""
foreach x in $oldEU  {
replace reg_o="oldEU" if iso_o=="`x'"
}
foreach x in $newEU {
replace reg_o="newEU" if iso_o=="`x'"
}
replace reg_o="nonEU" if reg_o==""

rename sec sec_id
merge m:1 sec_id using "$source/sectorlist"
drop _merge
save VAiso_reps, replace


u "$results//${scen}_smdB/cf_rva", 
rename rva vaa  
merge 1:1 iso* sec using "$results//${scen}_smdB/cf_rvap"
keep if _merge==3
drop _merge
rename rvap vaap
collapse (sum) va*, by(iso_o sec)


global oldEU "AUT BEL DEU DNK ESP FIN FRA GBR GRC IRL ITA LUX NLD PRT SWE"
global newEU "BGR CYP CZE EST HRV HUN LTU LVA MLT POL ROU SVK SVN"
g reg_o=""
foreach x in $oldEU {
replace reg_o="oldEU" if iso_o=="`x'"
}
foreach x in $newEU {
replace reg_o="newEU" if iso_o=="`x'"
}
replace reg_o="nonEU" if reg_o==""


rename sec sec_id
merge m:1 sec_id using "$source/sectorlist"

drop _merge
save "$results//${scen}_smdB/VAiso", replace

}

************** terms-of-trade changes

cd "$results"

foreach sc in allEU  { 
global scen ="`sc'"

u ${scen}1000_secB/EXiso_reps, clear
drop _merge
merge 1:1 iso sec_id rep using ${scen}1000_secB/VAiso_reps
	collapse (sum) exf exfp va vap, by(iso rep) fast
	g hvaUS = vap/va if iso=="USA"
	egen hVAUS=mean(hvaUS), by(rep)
	g hva=(vap/va-hVAUS)*100
	g hex=(exfp/exf-hVAUS)*100
	
	collapse (mean) hexm=hex hvam=hva (p5) hexp5=hex hvap5=hva (p95) hexp95=hex hvap95=hva, by(iso)

	g scen = "$scen"
	
reshape long hex hva, i(iso scen) j(unit) string

save ${scen}_smdB/vaghat_sd, replace


u ${scen}_smdB/EXiso, clear

merge 1:1 iso sec_id using ${scen}_smdB/VAiso

collapse (sum) vaa vaap EXf EXfp, by(iso)
	g hvaUS = vaap/vaa if iso=="USA"
	egen hVAUS=mean(hvaUS)
g hva=(vaap/vaa-hVAUS)*100
g hex = (EXfp/EXf-hVAUS)*100
keep iso hva hex
g scen = "$scen"
g unit = "b"

append using  ${scen}_smdB/vaghat_sd

sort iso unit

save ${scen}_smdB/vaghat, replace

}


u allEU_smdB/vaghat, clear
reshape wide hva hex, i(scen iso) j(unit) string

gsort scen  hvam
bysort scen: g n=_n


global labeliso ""
forvalues i=1(1)44 {
global val = iso[`i']
global labeliso  $labeliso  `i' "$val"
}

** Figure 4

twoway (rbar hvap95 hvap5 n, color(gs12) barwidth(1))  (scatter hvam n, color(black) msymbol(circle) msize(.7)) if scen=="allEU" ,  legend(off) ytitle("Wage level adjustment (in %)") xla($labeliso, angle(90)) scale(1.3)  xtitle("") graphregion(color(white)) xtitle("") yline(0)  xsize(10) ysize(3.5)
graph export "$figures/wage_allEU.pdf", replace
graph export "$figures/fig4.pdf", replace
graph export "$figures/fig4.eps", replace



****** sectoral value added - share changes -- aka employment

cd "$results"

foreach sc in $sclist {

global scen ="`sc'"

u ${scen}1000_secB/VAiso_reps, clear
replace reg_o="EU" if reg_o=="newEU" | reg_o=="oldEU"
collapse (sum) va vap, by(reg_o sec_id rep)
	egen tva=total(va), by(reg_o rep)
	egen tvap=total(vap), by(reg_o rep)
	g dsva=(vap/va-tvap/tva)*100

	
	global valist ""
foreach var of varlist dsva  {
foreach x in 1 25 5 95 975 995 {
    g `var'p`x'=.
}
levelsof reg_o, local(isolist)
global isolist "`isolist'"
foreach i in `isolist' {
forvalues j=1(1)50 {
_pctile `var' if reg=="`i'" & sec==`j', p(.5 2.5 5 95 97.5 99.5)
replace `var'p1=r(r1) if reg=="`i'" & sec==`j'
replace `var'p25=r(r2) if reg=="`i'" & sec==`j'
replace `var'p5=r(r3) if reg=="`i'" & sec==`j'
replace `var'p95=r(r4) if reg=="`i'" & sec==`j'
replace `var'p975=r(r5) if reg=="`i'" & sec==`j'
replace `var'p995=r(r6) if reg=="`i'" & sec==`j'
}
}
global valist "$valist" "`var'"
rename `var' `var'b
}

collapse (mean) dsva*, by(reg sec)

g scen="${scen}"

reshape long $valist, i(reg sec scen) j(unit) string

save ${scen}_smdB/vasds_ci, replace


u ${scen}_smdB/VAiso, clear
replace reg_o="EU" if reg_o=="newEU" | reg_o=="oldEU"
collapse (sum) vaa vaap, by(reg_o sec_id)
	egen tva=total(vaa), by(reg)
	egen tvap=total(vaap), by(reg)
g dsva=(vaap/vaa-tvap/tva)*100
keep reg sec_id ds
g scen = "$scen"
g unit = "m"

append using  ${scen}_smdB/vasds_ci
sort reg sec unit
save ${scen}_smdB/vasds, replace
}



******************************************************************************
**** sectoral value added changes 
*** note! currently the numeraire is average VA
*** everything will be scaled by US VA


foreach sc in $sclist {

global scen ="`sc'"

u ${scen}1000_secB/VAiso_reps, clear
replace reg_o="EU" if reg_o=="newEU" | reg_o=="oldEU"
egen vapUS=total(vap) if iso_o=="USA", by(rep)
egen vaUS=total(va) if iso_o=="USA", by(rep)
g hvaUS = vapUS/vaUS
collapse (sum) va vap (firstnm) hvaUS, by(reg_o sec_id rep)
egen hVAUS = mean(hvaUS), by(rep)
	g hva=(vap/va-hVAUS)*100
	
	
global valist ""
foreach var of varlist hva  {
foreach x in 1 25 5 95 975 995 {
    g `var'p`x'=.
}
levelsof reg_o, local(isolist)
global isolist "`isolist'"
foreach i in `isolist' {
forvalues j=1(1)50 {
_pctile `var' if reg=="`i'" & sec==`j', p(.5 2.5 5 95 97.5 99.5)
replace `var'p1=r(r1) if reg=="`i'" & sec==`j'
replace `var'p25=r(r2) if reg=="`i'" & sec==`j'
replace `var'p5=r(r3) if reg=="`i'" & sec==`j'
replace `var'p95=r(r4) if reg=="`i'" & sec==`j'
replace `var'p975=r(r5) if reg=="`i'" & sec==`j'
replace `var'p995=r(r6) if reg=="`i'" & sec==`j'
}
}
global valist "$valist" "`var'"
rename `var' `var'b
}

collapse (mean) hva*, by(reg sec)

g scen="${scen}"

reshape long $valist, i(reg sec scen) j(unit) string

save ${scen}_smdB/hvasec_ci, replace

u ${scen}_smdB/VAiso, clear
replace reg_o="EU" if reg_o=="newEU" | reg_o=="oldEU"
egen vapUS=total(vaap) if iso_o=="USA"
egen vaUS=total(vaa) if iso_o=="USA"
g hvaUS = vapUS/vaUS
collapse (sum) vaa vaap (firstnm) hvaUS, by(reg_o sec_id)
egen hVAUS = mean(hvaUS)
g hva=(vaap/vaa-hVAUS)*100
keep reg sec_id hva
g scen = "$scen"
g unit = "m"

append using  ${scen}_smdB/hvasec_ci
sort reg sec unit
save ${scen}_smdB/hvasec, replace
}



* table for paper

u "$results/Schengen_smdB/hvasec", replace
keep if unit=="" /*make empty dataset*/

foreach sc in $sclist {
global scen ="`sc'"
append using "$results//${scen}_smdB/hvasec"
}
drop if reg_o=="nonEU"

reshape wide hva, i(reg sec scen) j(unit) string

tostring hvab, g(shvab) format(%9.2f) force

foreach x in hva {
replace s`x'b = s`x'b+"$^{*}$" if `x'b<0 & `x'p95<0 & `x'p975>0
replace s`x'b = s`x'b+"$^{**}$" if `x'b<0 & `x'p975<0 & `x'p995>0
replace s`x'b = s`x'b+"$^{***}$" if `x'b<0 & `x'p995<0
replace s`x'b = s`x'b+"$^{*}$" if `x'b>0 & `x'p5>0 & `x'p25<0
replace s`x'b = s`x'b+"$^{**}$" if `x'b>0 & `x'p25>0 & `x'p1<0
replace s`x'b = s`x'b+"$^{***}$" if `x'b>0 & `x'p1>0
}

rename shvab hv
keep reg scen sec hv
reshape wide hv, i(reg sec) j(scen) string

merge m:1 sec_id using "$source/sectorlist"
drop _merge

merge m:1 sec_id using "$source/sectornames", keepusing(sectorname_short Isic)
drop _merge


replace Isic = subinstr(Isic,"_","\&",.)

replace sectorname = subinstr(sectorname,"&","\&",.)

* Table B8
* Table B9

listtex sectorna Isic sec_id hvCunion hvSmarket  hvEuro hvSchengen hvoRTA hvallEU hvallEUtran if sec_id<23 using "$tables/hvasec1.tex", replace end(\\)
listtex sectorna Isic sec_id hvCunion hvSmarket  hvEuro hvSchengen hvoRTA hvallEU hvallEUtran if sec_id>22 using "$tables/hvasec2.tex", replace end(\\)






*************************************************************************************************************************
********************* appendix table - additional aggregate flows at the regional level
*************************************************************************************************************************
* domestic sales and world exports


foreach sc in $sclist {
global scen ="`sc'"
cd "$results//${scen}1000_secB"

u EXVA_reps, clear

collapse (sum) ex* va* VAX* EX* (firstnm) hvaUS, by(reg_o rep)
egen hVAUS = mean(hvaUS), by(rep)

g D = exf-EX
g Dp = exfp-EXp
g DV = va - VAX
g DVp = vap - VAXp

g XhatWLD=(EXp/EX-hVAUS)*100
g XhatD=(Dp/D-hVAUS)*100

g VAXhatWLD=(VAXp/VAX-hVAUS)*100-XhatWLD
g VAXhatD=(DVp/DV-hVAUS)*100-XhatD



global valist ""
foreach var of varlist Xhat* VAXhat* {
foreach x in 1 25 5 95 975 995 {
    g `var'p`x'=.
}
levelsof reg_o, local(isolist)
global isolist "`isolist'"
foreach i in `isolist' {
_pctile `var' if reg=="`i'", p(.5 2.5 5 95 97.5 99.5)
replace `var'p1=r(r1) if reg=="`i'"
replace `var'p25=r(r2) if reg=="`i'"
replace `var'p5=r(r3) if reg=="`i'"
replace `var'p95=r(r4) if reg=="`i'"
replace `var'p975=r(r5) if reg=="`i'"
replace `var'p995=r(r6) if reg=="`i'"
}
global valist "$valist" "`var'"
rename `var' `var'b
}

collapse (mean) *hat*, by(reg)

g scenario="${scen}"

reshape long $valist, i(reg scen) j(unit) string

save ci_agg_apx, replace


u "$results//${scen}_smdB/EXVA", clear

collapse (sum) EX EXp vaa* VAX* EXf* (firstnm) hvaUS, by(reg_o)
egen hVAUS = mean(hvaUS)
drop hvaUS

g D = EXf-EX
g Dp = EXfp-EXp
g DV = vaa - VAX
g DVp = vaap - VAXp

g XhatWLD=(EXp/EX-hVAUS)*100
g XhatD=(Dp/D-hVAUS)*100
g VAXhatWLD=(VAXp/VAX-hVAUS)*100-XhatWLD
g VAXhatD=(DVp/DV-hVAUS)*100-XhatD

keep reg Xhat* VAXha*

g scenario="${scen}"
g unit = "m"
append using ci_agg_apx

save "$results//${scen}_smdB/tab_agg_apx", replace
}


cd "$results"

u Schengen_smdB/tab_agg_apx, replace
keep if unit=="" /*make empty dataset*/

foreach sc in $sclist {
global scen ="`sc'"
append using ${scen}_smdB/tab_agg_apx
}

reshape wide VAX* Xhat*, i(scen reg) j(unit) string

tostring VAXhatDb VAXhatWLDb, g(sVAXhatDb sVAXhatWLDb) format(%9.2f) force
tostring XhatDb XhatWLDb, g(sXhatDb sXhatWLDb) format(%9.2f) force

foreach x in VAXhatD XhatD VAXhatWLD XhatWLD  {
replace s`x'b = s`x'b+"$^{*}$" if `x'b<0 & `x'p95<0 & `x'p975>0
replace s`x'b = s`x'b+"$^{**}$" if `x'b<0 & `x'p975<0 & `x'p995>0
replace s`x'b = s`x'b+"$^{***}$" if `x'b<0 & `x'p995<0
replace s`x'b = s`x'b+"$^{*}$" if `x'b>0 & `x'p5>0 & `x'p25<0
replace s`x'b = s`x'b+"$^{**}$" if `x'b>0 & `x'p25>0 & `x'p1<0
replace s`x'b = s`x'b+"$^{***}$" if `x'b>0 & `x'p1>0
}


keep reg scen s*b
g sort = 1  if reg=="oldEU"
replace sort = 2 if reg=="newEU"
replace sort = 3 if reg=="nonEU"
sort scen sort 

replace reg ="\quad new EU" if reg_o=="newEU"
replace reg ="\quad old EU" if reg_o=="oldEU"
replace reg ="\quad non-EU" if reg_o=="nonEU"

* Table B6

foreach x in $sclist {
listtex reg sXhatD sVAXhatD sXhatWLD sVAXhatWLD using "$tables/agg_apx`x'.tex" if scen=="`x'", replace end(\\)
}













********  ROBUSTNESS. deficits, brexit, 



********* Tab: real income changes

global scenlist2_w "allEU1000_sec_noSB allEU1000_sec_oldB allEU1000_secBT"

foreach sc in $scenlist2_w {
global scen ="`sc'"
cd "$results//$scen"
u what, clear
replace what = (what-1)*100
*drop if rep==808 /* duplicate*/
foreach x in 1 25 5 95 975 995 {
    g p`x'=.
}
levelsof iso, local(isolist)
foreach i in `isolist' {
_pctile what if iso=="`i'", p(.5 2.5 5 95 97.5 99.5)
replace p1=r(r1) if iso=="`i'"
replace p25=r(r2) if iso=="`i'"
replace p5=r(r3) if iso=="`i'"
replace p95=r(r4) if iso=="`i'"
replace p975=r(r5) if iso=="`i'"
replace p995=r(r6) if iso=="`i'"
}
collapse (mean) p* what, by(iso)
rename p* whatp*
rename what whatb
reshape long what, i(iso) j(unit) string
save what_ci, replace
}


cd $results

u Brexit1000_secB/what, clear
replace what = (what-1)*100
*drop if rep==808 /* duplicate*/
foreach x in 1 25 5 95 975 995 {
    g p`x'=.
}
levelsof iso, local(isolist)
foreach i in `isolist' {
_pctile what if iso=="`i'", p(.5 2.5 5 95 97.5 99.5)
replace p1=r(r1) if iso=="`i'"
replace p25=r(r2) if iso=="`i'"
replace p5=r(r3) if iso=="`i'"
replace p95=r(r4) if iso=="`i'"
replace p975=r(r5) if iso=="`i'"
replace p995=r(r6) if iso=="`i'"
}
collapse (mean) p* what, by(iso)
rename p* whatp*
rename what whatb
reshape long what, i(iso) j(unit) string
rename what what_br
save Brexit_smdB/what_ci, replace




u AllEU1000_secB/what_pb, clear
merge 1:1 rep iso using AllEU1000_secB/what
replace what = (what-1)*100
replace what_pb = (what_pb-1)*100
g dwhat = what_pb-what
foreach var of varlist what_pb what dwhat {
    foreach x in 1 25 5 95 975 995 {
    g `var'p`x'=.
	}
levelsof iso, local(isolist)
foreach i in `isolist' {
_pctile `var' if iso=="`i'", p(.5 2.5 5 95 97.5 99.5)
replace `var'p1=r(r1) if iso=="`i'"
replace `var'p25=r(r2) if iso=="`i'"
replace `var'p5=r(r3) if iso=="`i'"
replace `var'p95=r(r4) if iso=="`i'"
replace `var'p975=r(r5) if iso=="`i'"
replace `var'p995=r(r6) if iso=="`i'"
}	
rename `var' `var'b
}
collapse (mean) what* dwhat*, by(iso)
reshape long what dwhat what_pb, i(iso) j(unit) string
rename what what_base
save AllEU_smdB/dwhat_ci, replace

u AllEU_smdB/dwhat_ci, clear
merge 1:1 iso unit using Brexit_smdB/what_ci
drop _mer
merge 1:1 iso unit using allEU1000_sec_oldB/what_ci
drop _mer
rename what what_oldB
merge 1:1 iso unit using allEU1000_sec_noSB/what_ci
drop _mer
rename what what_SB
merge 1:1 iso unit using allEU1000_secBT/what_ci
drop _mer
rename what what_BT


reshape wide what_BT what_SB what_oldB what_base what_br what_pb dwhat, i(iso) j(unit) string

* share of eu collapse losses accounted for by brexit

g s=what_brb/what_baseb

global oldEU "AUT BEL DEU DNK ESP FIN FRA GBR GRC IRL ITA LUX NLD PRT SWE"
global newEU "BGR CYP CZE EST HRV HUN LTU LVA MLT POL ROU SVK SVN"

g isostar=iso
foreach x in $oldEU {
replace isostar = iso + "$^o$" if iso=="`x'" 
}
foreach x in $newEU {
replace isostar = iso + "$^n$" if iso=="`x'" 
}


sum s if iso=="GBR"
sum s if iso=="IRL"
sum s if substr(isostar,-2,2)=="o$" & iso!="GBR"  & iso!="IRL"
sum s if substr(isostar,-2,2)=="n$"
drop s

foreach x in what_SB what_BT what_oldB what_base what_br what_pb dwhat {
tostring `x'b, g(s`x'b) format(%9.2f) force
}

foreach x in what_SB what_BT what_oldB what_base what_br what_pb dwhat {
replace s`x'b = s`x'b+"$^{*}$" if `x'b<0 & `x'p95<0 & `x'p975>0
replace s`x'b = s`x'b+"$^{**}$" if `x'b<0 & `x'p975<0 & `x'p995>0
replace s`x'b = s`x'b+"$^{***}$" if `x'b<0 & `x'p995<0
replace s`x'b = s`x'b+"$^{*}$" if `x'b>0 & `x'p5>0 & `x'p25<0
replace s`x'b = s`x'b+"$^{**}$" if `x'b>0 & `x'p25>0 & `x'p1<0
replace s`x'b = s`x'b+"$^{***}$" if `x'b>0 & `x'p1>0
}

keep iso*  s*b

** Table A8

listtex isostar *brb *pbb *baseb *whatb *BTb *oldBb *SBb using "$tables/real_inc_rob_new.tex", replace end(\\)













************* sectoral va changes

cd "$results//allEU1000_gfake_secB"

u rva, clear
rename rva va
merge 1:1 iso* rep sec using rvap
keep  if _merge==3
drop _merge
rename rvap vap
merge 1:1 iso* rep sec using exf
keep  if _merge==3
drop _merge
merge 1:1 iso* rep sec using exfp
keep  if _merge==3
drop _merge


global oldEU "AUT BEL DEU DNK ESP FIN FRA GBR GRC IRL ITA LUX NLD PRT SWE"
global newEU "BGR CYP CZE EST HRV HUN LTU LVA MLT POL ROU SVK SVN"

g reg_o=""
g reg_d=""
foreach x in $oldEU  {
replace reg_o="oldEU" if iso_o=="`x'"
replace reg_d="oldEU" if iso_d=="`x'"
}
foreach x in $newEU {
replace reg_o="newEU" if iso_o=="`x'"
replace reg_d="newEU" if iso_d=="`x'"
}
replace reg_o="nonEU" if reg_o==""
replace reg_d="nonEU" if reg_d==""

g EX=exf if iso_o!=iso_d
g EXp=exfp if iso_o!=iso_d
g VAX=va if iso_o!=iso_d
g VAXp=vap if iso_o!=iso_d

rename sec sec_id
merge m:1 sec_id using "$source/sectorlist"


*** for normalization in terms of US VA
egen vapUS = total(vap) if iso_o=="USA", by(rep)
egen vaUS = total(va) if iso_o=="USA", by(rep)

g hvaUS=vapUS/vaUS

collapse (sum) EX EXp VA* exf* va* (firstnm) hvaUS, by(reg* sec_agg rep) fast /*collapse over sectors and country groups*/

save EXVA_reps, replace



u "$results/allEU1000_secB/EXVA_reps", clear

collapse (sum) va vap (firstnm) hvaUS, by(reg_o sec rep)
egen hVAUS = mean(hva), by(rep)
drop hvaUS
egen tva = total(va), by(reg rep)
egen tvap = total(vap), by(reg rep)
g Bdsva=(vap/tvap-va/tva)*100
g Bhva=(vap/va-hVAUS)*100
keep reg_o sec rep Bds Bhva
save temp, replace

u EXVA_reps, clear

collapse (sum) va vap (firstnm) hvaUS, by(reg_o sec rep)
egen hVAUS = mean(hva), by(rep)
drop hvaUS
egen tva = total(va), by(reg rep)
egen tvap = total(vap), by(reg rep)
g dsva=(vap/tvap-va/tva)*100
g hva=(vap/va-hVAUS)*100

merge 1:1 sec rep reg_o using temp
g Ddsva=Bdsva-dsva
g Dhva=Bhva-hva



foreach var of varlist D* dsva hva {
    foreach x in 1 25 5 95 975 995 {
    g `var'p`x'=.
	}
levelsof reg, local(isolist)
levelsof sec, local(seclist)
foreach i in `isolist' {
    foreach j in `seclist' {
_pctile `var' if reg=="`i'" & sec=="`j'", p(.5 2.5 5 95 97.5 99.5)
replace `var'p1=r(r1) if reg=="`i'"  & sec=="`j'",
replace `var'p25=r(r2) if reg=="`i'"  & sec=="`j'",
replace `var'p5=r(r3) if reg=="`i'" & sec=="`j'",
replace `var'p95=r(r4) if reg=="`i'" & sec=="`j'",
replace `var'p975=r(r5) if reg=="`i'" & sec=="`j'",
replace `var'p995=r(r6) if reg=="`i'" & sec=="`j'",
}	
}
rename `var' `var'b
}

collapse (mean) dsva* hva* D*, by(reg sec)

g scen = "$scen"
g unit = "sd"

save vahat_ci, replace


replace sec_agg="Agric." if sec_agg=="Agriculture & Natural Resources"
replace sec_agg="Manuf." if sec_agg=="Manufacturing"
replace sec_agg="Serv." if sec_agg=="Services"


foreach x in dsva hva {
tostring `x'b, g(s`x'b) format(%9.2f) force
}

foreach x in dsva hva {
replace s`x'b = s`x'b+"$^{*}$" if `x'b<0 & `x'p95<0 & `x'p975>0
replace s`x'b = s`x'b+"$^{**}$" if `x'b<0 & `x'p975<0 & `x'p995>0
replace s`x'b = s`x'b+"$^{***}$" if `x'b<0 & `x'p995<0
replace s`x'b = s`x'b+"$^{*}$" if `x'b>0 & `x'p5>0 & `x'p25<0
replace s`x'b = s`x'b+"$^{**}$" if `x'b>0 & `x'p25>0 & `x'p1<0
replace s`x'b = s`x'b+"$^{***}$" if `x'b>0 & `x'p1>0
}

drop dsv* hva*

replace shvab = "\textbf{"+shvab+"}" if (Dhvap5>0 & Dhvap95>0) | (Dhvap5<0 & Dhvap95<0)
replace sdsvab = "\textbf{"+sdsvab+"}" if (Ddsvap5>0 & Ddsvap95>0) | (Ddsvap5<0 & Ddsvap95<0)


sort reg_o sec 
* table for paper
drop if reg_o=="nonEU"

g reg = reg_o if sec_agg=="Agric."
replace reg ="old EU" if reg=="oldEU"
replace reg ="new EU" if reg=="newEU"

* Table 6

listtex reg sec shvab sdsvab using "$tables/outputp_s_gfake.tex", replace end(\\)








cd "$results"

*** figure: comparison of transfer / no transfer scenario

import delim using AllEUtransfer_smdB/SD.txt, clear
rename v1 sd
rename v2 sd_prime
g ccode = _n
merge 1:1 cc using "$source/countrylist"
drop _merge
merge 1:1 iso using AllEUtransfer_smdB/cf_what
drop _merge
rename what what_tr
merge 1:1 iso using AllEU_smdB/cf_what
drop _merge
drop if sd==sd_prime
g Diff_what = (what_tr-what)*100
g tdiff = (sd-sd_pr)*100

*Figure 6

sort tdiff
graph hbar tdiff Diff, over(iso, sort(tdiff)) bar(1, color(black)) bar(2, color(gs12)) scale(.6) legend(label(1 "transfer cut (in % of income)") label(2 "additional income losses/gains in scenario with transfer cuts")) graphregion(color(white)) ysize(5) xsize(10) ylabel(-8(1)0)
graph export "$figures/transfer_effects.pdf", replace
graph export "$figures/fig6.pdf", replace
graph export "$figures/fig6.eps", replace




u "$source/fin_cftau_2014_weighted_wiod", clear

global EU "AUT BEL DEU DNK ESP FIN FRA GBR GRC IRL ITA LUX NLD PRT SWE BGR CYP CZE EST HRV HUN LTU LVA MLT POL ROU SVK SVN"

g eu_o=.
g eu_d=.
foreach x in $EU {
replace eu_o =1 if iso1=="`x'"
replace eu_d =1 if iso2=="`x'"
}

keep if eu_o==1 & eu_d==1
replace tau=tau*100

drop if sec>22
collapse (mean) tau, by(sector)
rename sector sec_id
merge 1:1 sec_id using "$source/sectornames", keepusing(sectorname_)
drop if _merg==2
drop _me

* Figure 2

graph bar tau, over(sectorna, sort(sec_id)  label(angle(90)) ) graphregion(color(white)) ytitle("Tariff change (in %pts.)") 
graph export "$figures/tariffsEU.pdf", replace
graph export "$figures/fig2.pdf", replace
graph export "$figures/fig2.eps", replace










**************** DESCRIPTIVES OF THE BASELINE: AGGREGATE (VA) TRADE 


cd "$main/data"


import delim ../simulation/Results/va_baseline.txt, clear
g ccode = ceil(_n/50)
g sector = _n-(ccode-1)*50

merge m:1 ccode using ../rawdata/countrylist
drop _merge ccode
rename iso iso_supply
reshape long v, i(iso_s sec) j(ccode)
merge m:1 ccode using ../rawdata/countrylist
drop _merge ccode
rename iso iso_demand
rename v va
save va_sectoral, replace

tostring sec, replace
replace sector = "C0"+sector if length(sector)==1
replace sector = "C"+sector if length(sector)==2

merge 1:1 iso* sec using WIOD_exp2014
drop _merge

rename exp EXf
rename va vaa
rename iso_s iso_o


global oldEU "AUT BEL DEU DNK ESP FIN FRA GBR GRC IRL ITA LUX NLD PRT SWE"
global newEU "BGR CYP CZE EST HRV HUN LTU LVA MLT POL ROU SVK SVN"

g reg_o=""
g reg_d=""
foreach x in $oldEU  {
replace reg_o="oldEU" if iso_o=="`x'"
replace reg_d="oldEU" if iso_d=="`x'"
}
foreach x in $newEU {
replace reg_o="newEU" if iso_o=="`x'"
replace reg_d="newEU" if iso_d=="`x'"
}
replace reg_o="nonEU" if reg_o==""
replace reg_d="nonEU" if reg_d==""


merge m:1 sector using ../rawdata/sectorlist
drop _merge

save va_baseline, replace

u va_baseline, clear
egen Y = total(EXf), by(reg_o)
g d = EXf if iso_o==iso_d
egen D = total(d), by(reg_o)

egen V= total(vaa), by(reg_o)
g da = vaa if iso_o==iso_d
egen DA = total(da), by(reg_o)

g EX = EXf if iso_o!=iso_d
g VA = vaa if iso_o!=iso_d

collapse (first) Y D V DA (sum) EX VA, by(reg_o reg_d)

reshape wide EX VA, i(reg_o Y D V DA) j(reg_d) string

foreach var of varlist Y D V DA EX* VA* {
	replace `var'=`var'/1000
}

g sort=1 if reg_o=="oldEU"
replace sort=2 if reg_o=="newEU"
replace sort=3 if reg_o=="nonEU"

replace reg_o= "old EU" if reg_o=="oldEU"
replace reg_o= "new EU" if reg_o=="newEU"
replace reg_o= "non-EU" if reg_o=="nonEU"
sort sort

format Y -VAold %9.0f


* Table A6

listtex reg Y D EXold EXnew EXnon using "../tables/aggexp_base_s.tex", replace end(\\) 
listtex reg V DA VAold VAnew VAnon using "../tables/aggVAexp_base_s.tex", replace end(\\) 

u va_baseline, clear
g EX = EXf if iso_o!=iso_d
g VA = vaa if iso_o!=iso_d
replace reg_d = "EU" if reg_d=="oldEU" | reg_d=="newEU"

collapse (sum) EX VA, by(reg_o reg_d sec_agg)


reshape wide EX VA, i(reg_o sec_ag) j(reg_d) string

g sort=1 if reg_o=="oldEU"
replace sort=2 if reg_o=="newEU"
replace sort=3 if reg_o=="nonEU"
replace reg_o= "old EU" if reg_o=="oldEU"
replace reg_o= "new EU" if reg_o=="newEU"
replace reg_o= "non-EU" if reg_o=="nonEU"

replace sec_agg="Agric." if sec_agg=="Agriculture & Natural Resources"
replace sec_agg="Manuf." if sec_agg=="Manufacturing"
replace sec_agg="Serv." if sec_agg=="Services"

sort sort sec_agg

g VAXEU=VAEU/EXEU
g VAXnonEU=VAnonEU/EXnonEU

foreach var of varlist EX* {
	replace `var'=`var'/1000
}
foreach var of varlist VAX* {
	replace `var'=`var'*100
}

format EX* %9.0f 
format VAX* %9.1f


* Table A7

listtex reg sec_a EXEU EXnonEU VAXEU VAXnonEU using "../tables/secbilexp_base_s.tex", replace end(\\) 



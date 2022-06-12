***** this do-files generates figures and table that combine the results from the 50-sector and 3-sector model

clear all
set matsize 10000


global main "Repository"
global source "$main/rawdata"
global tables "$main/tables"
global figures "$main/images"
global results "$main/simulation/Results"

cd "$results"



******************************************************
*** welfare changes in the four models
******************************************************

u allEU1000_secB/what, clear
rename what what_B
merge 1:1 iso rep using allEU1000_3secB/what
rename what what_3B
drop _merge
merge 1:1 iso rep using allEU1000_gfake_secB/what
rename what what_gfake
drop _merge
merge 1:1 iso rep using allEU1000_gfake_3secB/what
rename what what_gfake3
drop _merge
merge 1:1 iso rep using allEU1000_3EsecB/what
rename what what_3EB
drop _merge

g dw5=what_B-what_gfake

g d=what_3B-what_gfake3
g d2=what_gfake-what_gfake3
g d3=what_B-what_gfake3
g d5=what_B-what_gfake
g d6=what_3EB-what_3B

foreach var of varlist d* what* {
g p5`var'=`var'
g p95`var'=`var'
}

collapse (mean) wh* d* (p5) p5* (p95) p95*, by(iso)


foreach var of varlist what*  {
    replace `var'=(`var'-1)*100
}

foreach var of varlist d* p5* p95* {
    replace `var'=`var'*100
}

sort what_B
g n=_n



global labeliso ""
forvalues i=1(1)44 {
global val = iso[`i']
global labeliso  $labeliso  `i' "$val"
}


*Figure 9

twoway (scatter what_gfake3 n, color(black) msymbol(square_hollow)) (scatter what_gfake3 n, msymbol(square_hollow) color(black)), legend(order(1) label(1 "real consumption change"))  scale(.6) graphregion(color(white)) xlabel($labeliso, angle(90))   ytitle("percent") saving(gfake3, replace) xtitle("") title("(a) 3 sectors, simple IO") ylabel(-20(5)10)
twoway (rbar p5d p95d n , color(gs8)) (scatter d n, color(black%50) msize(.8)) (scatter what_3B n, msymbol(square_hollow) color(black)), scale(.6) graphregion(color(white)) xlabel($labeliso, angle(90)) legend(cols(3) order(3 2 1) label(1 "90% c.i. of difference") label(2 "Difference complex vs simple IO")  label(3  "real consumption change"))  ytitle("percent/percentage points") saving(sc3, replace) xtitle("") title("(b) 3 sectors, complex IO") ylabel(-20(5)10)
twoway (rbar p5d2 p95d2 n , color(gs8)) (scatter d2 n, color(black%50) msize(.8))  (scatter what_gfake n, msymbol(square_hollow) color(black)), scale(.6) graphregion(color(white)) xlabel($labeliso, angle(90)) legend(cols(3) order(3 2 1) label(1 "90% c.i. of difference") label(2 "Difference 50 vs 3 sectors") label(3 "real consumption change"))  ytitle("percent/percentage points") saving(tf, replace) xtitle("") title("(c) 50 sectors, simple IO") ylabel(-20(5)10)
twoway (rbar p5d5 p95d5 n , color(gs8)) (scatter d5 n, color(black%50) msize(.8))  (scatter what_B n, msymbol(square_hollow) color(black)), scale(.6) graphregion(color(white)) xlabel($labeliso, angle(90)) legend(cols(3) order(3 2 1) label(1 "90% c.i. of difference") label(2 "Difference complex vs simple IO")  label(3 "real consumption change"))  ytitle("percent/percentage points") saving(sc, replace) xtitle("") title("(d) 50 sectors, complex IO") ylabel(-20(5)10)

*Figure 9

graph combine gfake3.gph sc3.gph tf.gph   sc.gph, cols(2) ycommon scale(.8) graphregion(color(white))
graph export "$figures/fourways.pdf", replace
graph export "$figures/fig9.pdf", replace
graph export "$figures/fig9.eps", replace



*** table for online appendix with different model predictions ***

cd "$results"
u allEU1000_secB/what, clear
rename what what_B
merge 1:1 iso rep using allEU1000_3secB/what
rename what what_3B
drop _merge
merge 1:1 iso rep using allEU1000_gfake_secB/what
rename what what_gfake
drop _merge
merge 1:1 iso rep using allEU1000_gfake_3secB/what
rename what what_gfake3
drop _merge
merge 1:1 iso rep using allEU1000_3EsecB/what
rename what what_3EB
drop _merge

g d=what_3B-what_gfake3
g d2=what_gfake-what_gfake3
g d3=what_B-what_gfake3
g d5=what_B-what_gfake
g d6=what_3EB-what_3B

foreach var of varlist what*  {
    replace `var'=(`var'-1)*100
}

foreach var of varlist d* {
    replace `var'=`var'*100
}

foreach var of varlist d* what* {
foreach x in 1 25 5 95 975 995 {
    g p`x'`var'=.
}
levelsof iso, local(isolist)
foreach i in `isolist' {
_pctile `var' if iso=="`i'", p(.5 2.5 5 95 97.5 99.5)
replace p1`var'=r(r1) if iso=="`i'"
replace p25`var'=r(r2) if iso=="`i'"
replace p5`var'=r(r3) if iso=="`i'"
replace p95`var'=r(r4) if iso=="`i'"
replace p975`var'=r(r5) if iso=="`i'"
replace p995`var'=r(r6) if iso=="`i'"
}
}

collapse (mean) p* what* d*, by(iso)
sort what_B
g n=_n


global oldEU "AUT BEL DEU DNK ESP FIN FRA GBR GRC IRL ITA LUX NLD PRT SWE"
global newEU "BGR CYP CZE EST HRV HUN LTU LVA MLT POL ROU SVK SVN"

g isostar=iso
foreach x in $oldEU {
replace isostar = iso + "$^o$" if iso=="`x'" 
}
foreach x in $newEU {
replace isostar = iso + "$^n$" if iso=="`x'" 
}

foreach x in what_gfake3 what_3B what_gfake what_B d d2 d5 d3 what_3EB d6 {
tostring `x', g(s`x'b) format(%9.2f) force
}

foreach x in what_gfake3 what_3B what_gfake what_B d d2 d5 d3 what_3EB d6 {
rename p5`x' `x'p5 
rename  p95`x' `x'p95
rename  p25`x' `x'p25
rename p975`x' `x'p975 
rename p995`x' `x'p995 
rename p1`x' `x'p1 
}


foreach x in what_gfake3 what_3B what_gfake what_B d d2 d5 d3 what_3EB d6 {
replace s`x'b = s`x'b+"$^{*}$" if `x'<0 & `x'p95<0 & `x'p975>0
replace s`x'b = s`x'b+"$^{**}$" if `x'<0 & `x'p975<0 & `x'p995>0
replace s`x'b = s`x'b+"$^{***}$" if `x'<0 & `x'p995<0
replace s`x'b = s`x'b+"$^{*}$" if `x'>0 & `x'p5>0 & `x'p25<0
replace s`x'b = s`x'b+"$^{**}$" if `x'>0 & `x'p25>0 & `x'p1<0
replace s`x'b = s`x'b+"$^{***}$" if `x'>0 & `x'p1>0
}

keep isostar s*

*Table B10

listtex isostar swhat_gfake3 swhat_3B swhat_gfakeb swhat_B sdb sd2 sd5 sd3 swhat_3EB sd6 using "$tables/rob_models.tex", replace end(\\)





import delim using "../input_3s/bltau_sec-wagg.txt",  clear
forvalues i=1(1)3 {
replace v`i'=v`i'+1
} 
g coe="ltau"
g rep=_n
save "../input_3s/wagg_table", replace
local vars "schengen botheuro botheea eukor bothrta"
foreach x of local vars {
import delim using "../input_3s/b`x'_sec-wagg.txt", clear
g coe= "`x'"
g rep=_n
append using "../input_3s/wagg_table"
save "../input_3s/wagg_table", replace
}

reshape long v, i(coe rep) j(sec_agg)
reshape wide v, i(rep sec_agg) j(coe) string
rename v* *

foreach var of varlist both* euk lta schen {
foreach x in 1 25 5 95 975 995 {
    g p`x'`var'=.
}

forvalue i=1(1)3 {
_pctile `var' if sec==`i', p(.5 2.5 5 95 97.5 99.5)
replace p1`var'=r(r1) if sec==`i'
replace p25`var'=r(r2) if sec==`i'
replace p5`var'=r(r3) if sec==`i'
replace p95`var'=r(r4) if sec==`i'
replace p975`var'=r(r5) if sec==`i'
replace p995`var'=r(r6) if sec==`i'
}
}


collapse (mean) p*  both* euk lta schen , by(sec)


foreach x in  botheea botheuro eukor ltau schengen bothrta  {
rename p5`x' `x'p5 
rename  p95`x' `x'p95
rename  p25`x' `x'p25
rename p975`x' `x'p975 
rename p995`x' `x'p995 
rename p1`x' `x'p1 
}


foreach x in botheea botheuro eukor ltau schengen bothrta {
tostring `x', g(s`x'b) format(%9.2f) force
}

foreach x in botheea botheuro eukor schengen bothrta {
replace s`x'b = s`x'b+"$^{*}$" if `x'<0 & `x'p95<0 & `x'p975>0
replace s`x'b = s`x'b+"$^{**}$" if `x'<0 & `x'p975<0 & `x'p995>0
replace s`x'b = s`x'b+"$^{***}$" if `x'<0 & `x'p995<0
replace s`x'b = s`x'b+"$^{*}$" if `x'>0 & `x'p5>0 & `x'p25<0
replace s`x'b = s`x'b+"$^{**}$" if `x'>0 & `x'p25>0 & `x'p1<0
replace s`x'b = s`x'b+"$^{***}$" if `x'>0 & `x'p1>0
}

foreach x in botheea botheuro eukor schengen bothrta {
replace s`x'b = "--" if sec_agg==3
}

foreach x in ltau {
replace s`x'b = s`x'b+"$^{*}$" if `x'<0 & `x'p95<0 & `x'p975>0 & sec_agg==1
replace s`x'b = s`x'b+"$^{**}$" if `x'<0 & `x'p975<0 & `x'p995>0 & sec_agg==1
replace s`x'b = s`x'b+"$^{***}$" if `x'<0 & `x'p995<0 & sec_agg==1
replace s`x'b = s`x'b+"$^{*}$" if `x'>0 & `x'p5>0 & `x'p25<0 & sec_agg==1
replace s`x'b = s`x'b+"$^{**}$" if `x'>0 & `x'p25>0 & `x'p1<0 & sec_agg==1
replace s`x'b = s`x'b+"$^{***}$" if `x'>0 & `x'p1>0 & sec_agg==1
}

replace sltaub=sltaub+"$^\dagger$" if sec_agg>1

keep sec_agg s*

g sec="Manufacturing" if sec_a==1
replace sec="Services (tradable)" if sec_a==2
replace sec="Services (non-tradable)" if sec_a==3
 
 
 
* Table A5

listtex sec sbotheea sbotheur sschen sbothrt seuko   sltau using "$tables/wagg.tex", replace end(\\)




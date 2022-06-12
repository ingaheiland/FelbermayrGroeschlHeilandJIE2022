
* this do-file generates additional descriptive statistics


global main "Replication/Repository"
global source "$main/rawdata"
global tables "$main/tables"
global figures "$main/images"
global results "$main/simulation/Results"

cd "$results"


***Prepare Data for Figures

use "what_all.dta", clear
keep if unit=="m"
drop unit
rename what What_g
save What_g, replace

use "what_all.dta", clear
keep if unit=="p5"
drop unit
rename what w_r_l
merge 1:1 iso scen using What_g
drop _m
save What_g, replace

use what_all.dta, clear
keep if unit=="p95"
drop unit
rename what w_r_u
merge 1:1 iso scen using What_g
drop _m
save What_g, replace



***EEA-membership
generate EEA=0
label var EEA "member of European Economic Area"
local EEA "BEL DNK DEU FIN FRA GRC IRL ITA LUX NLD PRT ESP GBR AUT SWE NOR LIE ISL" 
foreach h of local EEA{
replace EEA=1 if iso=="`h'"
}

* EU-Osterweiterung 2004/2007/2013
local EU2004 "CZE CYP EST HUN LTU LVA MLT POL SVK SVN BGR ROU HRV" 
foreach l of local EU2004 {
replace EEA=1 if iso=="`l'"
} 

drop if EEA==0
save income_confidence_new, replace

*** Figure Appendix: Percentage Change in Real consumption relative to Status Quo, Various Scenarios
use income_confidence_new, clear
keep if scen=="Smarket"
sort What_g
g id=_n
labmask id, values(iso)
twoway (line What_g id) || (line w_r_l id, lpattern(shortdash) lc(gs10)) || (line w_r_u id, lpattern(shortdash) lc(gs10) graphregion(color(white)) legend(off) ylabel(0(5)-25, labsize(small)) xlabel(1(1)28, valuelabel angle(90) labsize(vsmall) notick) title("Single Market", size(medium)) xtitle(" ") saving(line0, replace))

use income_confidence_new, clear
keep if scen=="Cunion"
sort What_g
g id=_n
labmask id, values(iso)
twoway (line What_g id) || (line w_r_l id, lpattern(shortdash) lc(gs10)) || (line w_r_u id, lpattern(shortdash) lc(gs10) graphregion(color(white)) legend(off) ylabel(.2(.2)-1, labsize(small)) xlabel(1(1)28, valuelabel angle(90) labsize(vsmall) notick) title("Customs Union", size(medium)) xtitle(" ") saving(line1, replace))


use income_confidence_new, clear
keep if scen=="Euro"
sort What_g
g id=_n
labmask id, values(iso)
twoway (line What_g id) || (line w_r_l id, lpattern(shortdash) lc(gs10)) || (line w_r_u id, lpattern(shortdash) lc(gs10) graphregion(color(white)) legend(off) ylabel(0(1)-5, labsize(small)) xlabel(1(1)28, valuelabel angle(90) labsize(vsmall) notick) title("Euro", size(medium)) xtitle(" ") saving(line2, replace))


use income_confidence_new, clear
keep if scen=="Schengen"
sort What_g
g id=_n
labmask id, values(iso)
twoway (line What_g id) || (line w_r_l id, lpattern(shortdash) lc(gs10)) || (line w_r_u id, lpattern(shortdash) lc(gs10) graphregion(color(white)) legend(off) ylabel(3(1)-5, labsize(small)) xlabel(1(1)28, valuelabel angle(90) labsize(vsmall) notick) title("Schengen", size(medium)) xtitle(" ") saving(line3, replace))

use income_confidence_new, clear
keep if scen=="oRTA"
sort What_g
g id=_n
labmask id, values(iso)
twoway (line What_g id) || (line w_r_l id, lpattern(shortdash) lc(gs10)) || (line w_r_u id, lpattern(shortdash) lc(gs10) graphregion(color(white)) legend(off) ylabel(1.5(.5)-2, labsize(small)) xlabel(1(1)28, valuelabel angle(90) labsize(vsmall) notick) title("Other RTAs", size(medium)) xtitle(" ") saving(line4, replace))


use income_confidence_new, clear
keep if scen=="allEU"
sort What_g
g id=_n
labmask id, values(iso)
twoway (line What_g id) || (line w_r_l id, lpattern(shortdash) lc(gs10)) || (line w_r_u id, lpattern(shortdash) lc(gs10) graphregion(color(white)) legend(off) ylabel(0(5)-30, labsize(small)) xlabel(1(1)28, valuelabel angle(90) labsize(vsmall) notick) title("Complete EU", size(medium)) xtitle(" ") saving(line5, replace))

use income_confidence_new, clear
keep if scen=="allEUtransfer"
sort What_g
g id=_n
labmask id, values(iso)
twoway (line What_g id) || (line w_r_l id, lpattern(shortdash) lc(gs10)) || (line w_r_u id, lpattern(shortdash) lc(gs10) graphregion(color(white)) legend(off) ylabel(0(5)-30, labsize(small)) xlabel(1(1)28, valuelabel angle(90) labsize(vsmall) notick) title("Complete EU incl. Fiscal Transfers", size(medium)) xtitle(" ") saving(line6, replace))


* Figure B1

graph combine line0.gph line1.gph line2.gph line3.gph line4.gph line5.gph line6.gph, graphregion(color(white))
graph export "$figures/scenarios_rank.pdf", replace
graph export "$figures/figB1.pdf", replace
graph export "$figures/figB1.eps", replace




***Figure: Change in Real Consumption in % for Various Scenarios, Baseline Year 2014
use "what_all.dta", clear
keep if unit=="m"
drop unit
rename what What_g
greshape wide What_g, i(iso) j(scen)

global scen "allEU allEUtransfer Euro Cunion Smarket oRTA Schengen"

foreach x in $scen {
rename What_g`x' What_g_`x'
}

g total=What_g_Smarket+What_g_Schengen+What_g_Euro+What_g_oRTA+What_g_Cunion
sort total

** Figure 5

graph hbar What_g_Smarket What_g_Schengen What_g_Euro What_g_oRTA What_g_Cunion, over(iso, sort(total) label(labsize(vsmall))) stack ylabel(, labsize(vsmall)) bar(1, color(dknavy)) bar(2, color(eltblue)) bar(3, color(green)) bar(4, color(blue)) bar(5, color(emidblue)) graphregion(color(white)) legend(ring(0) pos(7) col(1) lab(1 "Single Market") lab(2 "Schengen") lab(3 "Euro") lab(4 "Other RTAs") lab(5 "MFN Tariffs") region(lcolor(white)) size(vsmall)) outergap(-50) xalt 
graph export "$figures/income_change_sec.pdf", replace
graph export "$figures/fig5.pdf", replace
graph export "$figures/fig5.eps", replace




***Figure: Correlation Losses and Country Characteristics: Size, the Level of Real consumption, Remoteness, and Openness
use "what_all.dta", clear
keep if unit=="m"
drop unit
rename what What_g
greshape wide What_g, i(iso) j(scen)

global scen "allEU allEUtransfer Euro Cunion Smarket oRTA Schengen"

foreach x in $scen {
rename What_g`x' What_g_`x'
}

***EEA-membership
generate EEA=0

local EEA "BEL DNK DEU FIN FRA GRC IRL ITA LUX NLD PRT ESP GBR AUT SWE NOR LIE ISL" 
foreach h of local EEA{
replace EEA=1 if iso=="`h'"
}

* EU-Osterweiterung 2004/2007/2013
local EU2004 "CZE CYP EST HUN LTU LVA MLT POL SVK SVN BGR ROU HRV" 
foreach l of local EU2004 {
replace EEA=1 if iso=="`l'"
} 

drop if EEA==0

save "scatter_new_what.dta", replace


use "$source/scatter4graphs", clear // 

drop What*
merge 1:1 iso using scatter_new_what



*** Correlation and Figure: Population (log), 1995
reg What_g_allEUtransfer lpop, rob
reg What_g_allEUtransfer lpop [aweight=pop],rob
twoway (scatter What_g_allEUtransfer lpop [aweight = pop], mcolor(black%25) msymbol(0h) ) (lfit What_g_allEUtransfer lpop, lpattern(dash)) (lfit What_g_allEUtransfer lpop [aweight = pop]),   graphregion(color(white)) xtitle("Population (log), 1995") ytitle("Change in real consumption (%)") legend(off) scale(.7) xlabel(13(1)19) ylabel(0(5)-35) saving(welfare3, replace)
*graph export pop_welfare.pdf, replace


*** Correlation and Figure: Per capita income (log), 1995
reg What_g_allEUtransfer linpc,rob
reg What_g_allEUtransfer linpc [aweight=pop],rob
twoway (scatter What_g_allEUtransfer linpc [aweight = pop], mcolor(black%25) msymbol(0h) ) (lfit What_g_allEUtransfer linpc, lpattern(dash)) (lfit What_g_allEUtransfer linpc [aweight = pop]),  graphregion(color(white)) xtitle("Per capita income (log), 1995") ytitle("Change in real consumption (%)") legend(off) scale(.7) xlabel(7(.5)11) ylabel(0(5)-35) saving(welfare1, replace)
*graph export linpc_welfare.pdf, replace

*** Correlation and Figure: Average distance from other EEA members (log)
reg What_g_allEUtransfer ldist, rob
reg What_g_allEUtransfer ldist [aweight=pop], rob
twoway (scatter What_g_allEUtransfer ldist [aweight = pop], mcolor(black%25) msymbol(0h) ) (lfit What_g_allEUtransfer ldist, lpattern(dash)) (lfit What_g_allEUtransfer ldist [aweight = pop]),  graphregion(color(white)) xtitle("Average distance from other EEA members(log)") ytitle("Change in real consumption (%)") legend(off) scale(.7) xlabel(6.7(.1)7.7) ylabel(0(5)-35) saving(welfare4, replace)
*graph export avdist_welfare.pdf, replace

*** Correlation and Figure: Trade openness (%), 1995
reg What_g_allEUtransfer open, rob
reg What_g_allEUtransfer open [aweight=pop], rob
twoway (scatter What_g_allEUtransfer open [aweight = pop], mcolor(black%25) msymbol(0h) ) (lfit What_g_allEUtransfer open, lpattern(dash)) (lfit What_g_allEUtransfer open [aweight = pop]),  graphregion(color(white)) xtitle("Trade openness (%), 1995") ytitle("Change in real consumption (%)") legend(off) scale(.7) xlabel(25(25)105) ylabel(0(5)-35) saving(welfare2, replace)
*graph export open_welfare.pdf, replace


*Figure 7

graph combine "welfare3.gph" "welfare1.gph" "welfare4.gph" "welfare2.gph", ycommon rows(2) graphregion(fcolor(white) lcolor(white) ifcolor(white) ilcolor(white))
graph export "$figures/combined.pdf", replace
graph export "$figures/fig7.pdf", replace
graph export "$figures/fig7.eps", replace


reg What_g_allEUtransfer lpop linpc ldist open, rob
reg What_g_allEUtransfer lpop linpc ldist open [aweight=pop], rob

reg What_g_allEUtransfer lpop linpc lav_dist open, rob
reg What_g_allEUtransfer lpop linpc lav_dist open [aweight=pop], rob



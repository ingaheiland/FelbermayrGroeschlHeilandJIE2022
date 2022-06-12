
** this do-file generates the regression tables 
* 1: coef_final_adj.tex
* A2: theta_org_adj_constraint.tex
* A3: coef_trend_final_adj.tex
* A4: theta_trend_org_adj_constraint.tex




cd "Repository/estimation"


u coef_final_adj,clear

*u coef_final_org,clear
reshape long bltau_ bbothEEA_ bschengen_ bbotheuro_ bbothrta_ beukor_ beupta_, i(draw) j(sec)

foreach var of varlist b* {
g p1`var'=.
g p25`var'=.
g p5`var'=.
g p95`var'=.
g p975`var'=.
g p99`var'=.
}

forvalues i=1(1)50 {
foreach var of varlist b* {
_pctile `var' if sec==`i', p(.5 2.5 5 95 97.5 99.5)
replace p1`var'=r(r1) if sec==`i'
replace p25`var'=r(r2) if sec==`i'
replace p5`var'=r(r3) if sec==`i'
replace p95`var'=r(r4) if sec==`i'
replace p975`var'=r(r5) if sec==`i'
replace p99`var'=r(r6) if sec==`i'
}
}


collapse (mean) *b* (sd) sdbltau=bltau sdbothEEA=bbothEEA sdbothrta=bbothrta sdbotheuro=bbotheuro sdschengen=bschengen sdeupta=beupta, by(sec)

foreach var of varlist bbothEEA_ bschengen_ bbotheuro_ bbothrta_ beukor_ beupta_ {
tostring `var', g(s`var') force format(%9.2f)
g ss`var'=s`var'
replace ss`var'=s`var'+"$^*$" if p5`var'>=0
replace ss`var'=s`var'+"$^{**}$" if p25`var'>=0
replace ss`var'=s`var'+"$^{***}$" if p1`var'>=0
replace ss`var'=s`var'+"$^*$" if p95`var'<=0
replace ss`var'=s`var'+"$^{**}$" if p975`var'<=0
replace ss`var'=s`var'+"$^{***}$" if p99`var'<=0
}

rename sec sec_id
merge 1:1 sec_id using ../rawdata/sectorlist
drop _merge

replace sectorname_short = subinstr(sectorname_s,"&","\&",.)

tostring p1bltau, g(c1bltau) force format(%9.2f)
tostring p99bltau, g(c99bltau) force format(%9.2f)


listtex sectorname_short sec_id ssbbothEEA ssbbotheuro ssbschengen ssbeukor ssbeupta ssbbothrta using ../tables/coef_final_adj.tex, replace end(\\)


compress
keep  *bltau* sectorname_short sec_id

replace bltau_ = bltau+1
replace p5bltau = p5bltau+1
replace p95bltau = p95bltau+1

tostring bltau, g(bltau) force format(%9.2f)

tostring p5bltau, g(c5bltau) force format(%9.2f)
tostring p95bltau, g(c95bltau) force format(%9.2f)

keep if sec_id<24

replace secto="Services (23-50)" if sec_id==23
tostring sec_id, replace
replace sec_id="23-50" if sec_id=="23"

g ci = "["+c5bltau+"; "+c95bltau+"]"

replace ci="-" if sec_id=="23-50"


save theta_adj_tab, replace



* count thetas at the constraint

u coef_final_org, clear

keep blta* draw

forvalues i=1(1)22 {
	replace bltau_`i'=. if bltau_`i'<-1
}

collapse (count) bl*

g dum=1
reshape long bltau_, i(dum) j(sec_id)

g constraint=blta/1000
drop bl*

ge sec=sec_id
tostring sec_id, replace

merge 1:1 sec_id using theta_adj_tab
drop _merge

tostring constr, format(%9.2f) replace force
replace constr="-" if sec_id=="23-50"

sort sec
listtex sectorname_short sec_id bltau ci constraint using ../tables/theta_org_adj_constraint.tex, replace end(\\)

save theta_adj_constr_tab, replace






******************* coefficients estimated with trends



u coef_trend_final_adj,clear


reshape long bltau_ bbothEEA_ bschengen_ bbotheuro_ bbothrta_ beukor_ beupta_, i(draw) j(sec)

foreach var of varlist b* {
g p1`var'=.
g p25`var'=.
g p5`var'=.
g p95`var'=.
g p975`var'=.
g p99`var'=.
}

forvalues i=1(1)50 {
foreach var of varlist b* {
_pctile `var' if sec==`i', p(.5 2.5 5 95 97.5 99.5)
replace p1`var'=r(r1) if sec==`i'
replace p25`var'=r(r2) if sec==`i'
replace p5`var'=r(r3) if sec==`i'
replace p95`var'=r(r4) if sec==`i'
replace p975`var'=r(r5) if sec==`i'
replace p99`var'=r(r6) if sec==`i'
}
}


collapse (mean) *b* (sd) sdbltau=bltau sdbothEEA=bbothEEA sdbothrta=bbothrta sdbotheuro=bbotheuro sdschengen=bschengen sdeupta=beupta, by(sec)

foreach var of varlist bbothEEA_ bschengen_ bbotheuro_ bbothrta_ beukor_ beupta_ {
tostring `var', g(s`var') force format(%9.2f)
g ss`var'=s`var'
replace ss`var'=s`var'+"$^*$" if p5`var'>=0
replace ss`var'=s`var'+"$^{**}$" if p25`var'>=0
replace ss`var'=s`var'+"$^{***}$" if p1`var'>=0
replace ss`var'=s`var'+"$^*$" if p95`var'<=0
replace ss`var'=s`var'+"$^{**}$" if p975`var'<=0
replace ss`var'=s`var'+"$^{***}$" if p99`var'<=0
}

rename sec sec_id
merge 1:1 sec_id using ../rawdata/sectorlist
drop _merge

replace sectorname_short = subinstr(sectorname_s,"&","\&",.)

tostring p1bltau, g(c1bltau) force format(%9.2f)
tostring p99bltau, g(c99bltau) force format(%9.2f)


listtex sectorname_short sec_id ssbbothEEA ssbbotheuro ssbschengen ssbeukor ssbeupta ssbbothrta using ../tables/coef_trend_final_adj.tex, replace end(\\)


compress
keep  *bltau* sectorname_short sec_id

replace bltau_ = bltau+1
replace p5bltau = p5bltau+1
replace p95bltau = p95bltau+1

tostring bltau, g(bltau) force format(%9.2f)

tostring p5bltau, g(c5bltau) force format(%9.2f)
tostring p95bltau, g(c95bltau) force format(%9.2f)

keep if sec_id<24

replace secto="Services (23-50)" if sec_id==23
tostring sec_id, replace
replace sec_id="23-50" if sec_id=="23"

g ci = "["+c5bltau+"; "+c95bltau+"]"

replace ci="-" if sec_id=="23-50"



save theta_trend_adj_tab, replace


* count thetas at the constraint

u coef_trend_final_org, clear

keep blta* draw

forvalues i=1(1)22 {
	replace bltau_`i'=. if bltau_`i'<-1
}

collapse (count) bl*

g dum=1
reshape long bltau_, i(dum) j(sec_id)

g constraint=blta/1000
drop bl*

g sec=sec_id
tostring sec_id, replace



merge 1:1 sec_id using theta_trend_adj_tab
drop _merge

tostring constr, format(%9.2f) replace force
replace constr="-" if sec_id=="23-50"

sort sec
listtex sectorname_short sec_id bltau ci constraint using ../tables/theta_trend_org_adj_constraint.tex, replace end(\\)

save theta_trend_adj_constr_tab, replace



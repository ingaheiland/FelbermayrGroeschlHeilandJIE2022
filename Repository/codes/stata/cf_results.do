clear
set more off
set matsize 10000



global main "Repository"
global source "$main/rawdata"

global scenlist "allEU_smdB Euro_smdB Cunion_smdB Smarket_smdB oRTA_smdB Schengen_smdB allEUtransfer_smdB allEU_gfake_smdB"


foreach sc in $scenlist {
global scen ="`sc'"
cd "$main/simulation/Results/$scen"

foreach x in what wrhat  {
import delim cf_`x'.txt, clear
g ccode = _n
merge 1:1 ccode using "$source/countrylist"
drop _merge
rename v `x'
save cf_`x', replace
}

foreach x in rva rvap {
import delim cf_`x'.txt, clear
g sec=ceil(_n/44)
g ccode = _n-(sec-1)*44
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename iso iso_o
reshape long v, i(iso_o sec) j(ccode)
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename iso iso_d
rename v `x'
save cf_`x', replace
}


foreach x in EXf EXfp {
import delim cf_`x'.txt, clear
g sec=ceil(_n/44)
g ccode = _n-(sec-1)*44
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename iso iso_d
reshape long v, i(iso_d sec) j(ccode)
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename iso iso_o
rename v `x'
save cf_`x', replace
}


}

global scenlist2 "CUnion1000_secB  Smarket1000_secB Schengen1000_secB oRTA1000_secB AllEU1000_secB AllEUtransfer1000_secB Euro1000_secB AllEU1000_gfake_secB"



foreach sc in $scenlist2 {
global scen ="`sc'"
cd "$main/simulation/Results/$scen"
foreach x in rva rvap {
import delim `x'.txt, clear
g rep=ceil(_n/(44*50))
g n= _n-(rep-1)*50*44
g sec=ceil(n/44)
g ccode = n-(sec-1)*44
drop n
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename iso iso_o
reshape long v, i(rep sec iso_o) j(ccode)
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename iso iso_d
rename v `x'
save `x', replace
}

foreach x in exf exfp  {
import delim `x'.txt, clear
g rep=ceil(_n/(44*50))
g n= _n-(rep-1)*50*44
g sec=ceil(n/44)
g ccode = n-(sec-1)*44
drop n
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename iso iso_d
reshape long v, i(rep sec iso_d) j(ccode)
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename iso iso_o
rename v `x'
save `x', replace
}

foreach x in what wr  {
import delim using `x'.txt, clear
destring v, replace force
g rep=ceil(_n/44)
g ccode = _n-(rep-1)*44
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename v `x'
save `x', replace
}
}




*** post-BREXIT scenarios

global scenlist_pb "AllEU_smdB"

foreach sc in $scenlist_pb {
global scen ="`sc'"
cd "$main/simulation/Results/$scen"

foreach x in what_pb {
import delim cf_`x'.txt, clear
g ccode = _n
merge 1:1 ccode using "$source/countrylist"
drop _merge
rename v `x'
save cf_`x', replace
}
}


global scenlist2_pb "AllEU1000_secB"

foreach sc in $scenlist2_pb {
global scen ="`sc'"
cd "$main/simulation/Results/$scen"

foreach x in what_pb  {
import delim using `x'.txt, clear
g rep=ceil(_n/44)
g ccode = _n-(rep-1)*44
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename v `x'
save `x', replace
}
}


***** BREXIT

global scenlist_pb "Brexit_smdB"

foreach sc in $scenlist_pb {
global scen ="`sc'"
cd "$main/simulation/Results/$scen"


foreach x in what {
import delim `x'.txt, clear
g ccode = _n
merge 1:1 ccode using "$source/countrylist"
drop _merge
rename v `x'
save cf_`x', replace
}
}


global scenlist2_pb "Brexit1000_secB"

foreach sc in $scenlist2_pb {
global scen ="`sc'"
cd "$main/simulation/Results/$scen"
foreach x in what  {
import delim using `x'.txt, clear
g rep=ceil(_n/44)
g ccode = _n-(rep-1)*44
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename v `x'
save `x', replace
}
}

*more robustness where only welfare effects are considered:

*old deficit, zero deficit

global scenlist_w "AllEU_smd_noSB AllEU_smd_oldB AllEU_smdBT"


foreach sc in $scenlist_w {
global scen ="`sc'"
cd "$main/simulation/Results/$scen"

foreach x in what {
import delim cf_`x'.txt, clear
g ccode = _n
merge 1:1 ccode using "$source/countrylist"
drop _merge
rename v `x'
save cf_`x', replace
}
}


global scenlist2_w "allEU1000_sec_noSB allEU1000_sec_oldB allEU1000_secBT"

foreach sc in $scenlist2_w {
global scen ="`sc'"
cd "$main/simulation/Results/$scen"

foreach x in what  {
import delim using `x'.txt, clear
g rep=ceil(_n/44)
g ccode = _n-(rep-1)*44
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename v `x'
save `x', replace
}
}









********************** 3-sector models ***************************

global scenlist3s "allEU_3smd allEU_gfake_3smd"

foreach sc in $scenlist3s {
global scen ="`sc'"
cd "$main/simulation/Results/$scen"


foreach x in what wrhat {
import delim cf_`x'.txt, clear
g ccode = _n
merge 1:1 ccode using "$source/countrylist"
drop _merge
rename v `x'
save cf_`x', replace
}

/*
foreach x in EXf EXfp {
import delim cf_`x'.txt, clear
g sec=ceil(_n/44)
g ccode = _n-(sec-1)*44
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename iso iso_d
reshape long v, i(iso_d sec) j(ccode)
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename iso iso_o
rename v `x'
save cf_`x', replace
}



foreach x in vaa vaap {
import delim cf_`x'.txt, clear
g sec=ceil(_n/44)
g ccode = _n-(sec-1)*44
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename iso iso_o
reshape long v, i(iso_o sec) j(ccode)
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename iso iso_d
rename v `x'
save cf_`x', replace
}
*/

}


global scenlist3sa "allEU1000_3secB allEU1000_gfake_3secB"

foreach sc in $scenlist3sa {
global scen ="`sc'"
cd "$main/simulation/Results/$scen"

/*
foreach x in va vap {
import delim `x'.txt, clear
g rep=ceil(_n/(44*2))
g n= _n-(rep-1)*3*44
g sec=ceil(n/44)
g ccode = n-(sec-1)*44
drop n
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename iso iso_o
reshape long v, i(rep sec iso_o) j(ccode)
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename iso iso_d
rename v `x'
save `x', replace
}


foreach x in exf exfp  {
import delim `x'.txt, clear
g rep=ceil(_n/(44*2))
g n= _n-(rep-1)*3*44
g sec=ceil(n/44)
g ccode = n-(sec-1)*44
drop n
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename iso iso_d
reshape long v, i(rep sec iso_d) j(ccode)
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename iso iso_o
rename v `x'
save `x', replace
}
*/

foreach x in what wr  {
import delim using `x'.txt, clear
destring v, replace force
g rep=ceil(_n/44)
g ccode = _n-(rep-1)*44
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename v `x'
save `x', replace
}
}



global scenlist3sb "allEU1000_3EsecB"

foreach sc in $scenlist3sb {
global scen ="`sc'"
cd "$main/simulation/Results/$scen"


foreach x in what wr  {
import delim using `x'.txt, clear
destring v, replace force
g rep=ceil(_n/44)
g ccode = _n-(rep-1)*44
merge m:1 ccode using "$source/countrylist"
drop _merge ccode
rename v `x'
save `x', replace
}
}


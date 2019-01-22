
use $data\finalhousechar, replace
	
	drop if  (sold_in_yr_bl==1 | constnewandprior==1)

	global xbhousing " age sa_sqft1000 sa_lotsize sa_lotsize2 distMSA"	

drop if nhousespuddist10000>1
/* want to get wells that only have 1 well drilled (at present or in the future) */
bysort apn:	egen maxwellsin10000=max(nhousespuddist10000)
keep if maxwellsin10000==1 
	

	bysort apn:	egen minmindistspud=min(mindistspud) 
	
	bysort apn:	egen minmindistdate=max(mindistdate) 
	
	format minmindistdate %td
	

bysort apn: gen st=_N
keep if st>1

 

keep if statecode==42 & Marcellus==1
egen countystateyear=group(countycod statecode year)
tab countystateyear, gen(_IyeaXcou_)
gen q=quarter(saledate)
tab quarter, gen(dq_)

****** SPECIFY TIME VARIABLE *******
global pre_var nhousespuddist10000


* 3 X Silverman's Rule of Thumb
qui: sum minmindistspud
global h = 3*1.06*r(sd)*(_N)^(-1/5)

gen temp = int(${h})
local band = temp[1]
drop temp
************************************
** Price Gradient wrt to distance (Adapted from Lala Ma's code) **
************************************


di "${h}"
cd $data

	sort minmindistspud
gen points = .
forvalues i =1/1000 {
	qui:replace points = `i'*10 in `i'
}



*qui 
reg logprice2012 $xbhousing _IyeaXcou_*  dq_* ,robust cluster(sa_census_tract)
predict price_residg, residuals


qui lpoly price_residg minmindistspud if ${pre_var} == 0  & groundwater==1, generate(yhat_before2) at(points)  bwidth(${h}) degree(1) kernel(gaussian) msymbol(oh) msize(small) mcolor(gs10) ciopts(lwidth(medium)) noscatter nograph
*692 ob
local n_beforeg = r(N)
scalar define nbefg=r(N)
qui lpoly price_residg minmindistspud if ${pre_var} == 1  & groundwater==1, generate(yhat_after2) at(points)  bwidth(${h}) degree(1) kernel(gaussian) msymbol(oh) msize(small) mcolor(gs10) ciopts(lwidth(medium)) noscatter nograph
*731 obs
local n_afterg = r(N)
scalar define naftg=r(N)


reg logprice2012 $xbhousing _IyeaXcou_*  dq_* ,robust cluster(sa_census_tract)
predict price_residp, residuals


qui lpoly price_residp minmindistspud if ${pre_var} == 0  & groundwater==0, generate(yhat_before2p) at(points)  bwidth(${h}) degree(1) kernel(gaussian) msymbol(oh) msize(small) mcolor(gs10) ciopts(lwidth(medium)) noscatter nograph
*6332 obs
local n_beforep = r(N)
scalar define nbefp=r(N)
qui lpoly price_residp minmindistspud if ${pre_var} == 1  & groundwater==0, generate(yhat_after2p) at(points)  bwidth(${h}) degree(1) kernel(gaussian) msymbol(oh) msize(small) mcolor(gs10) ciopts(lwidth(medium)) noscatter nograph
*6547 obs
local n_afterp = r(N)
scalar define naftg=r(N)

keep points yhat_after2 yhat_before2 yhat_after2p yhat_before2p
drop if points==.


*save  original_LR, replace


twoway (line yhat_before2 points, lcolor(black) lpattern(solid))  (line yhat_after2 points, lcolor(black) lpattern(dash)), /*
*/ xtitle("Distance from Well (in meters)", size(small)) note("Bandwidth= `band' meters, N1 = `n_beforeg', N2 = `n_afterg'") ytitle("Log Price Residuals", size(small)) /*
*/ legend(order(1 "Before Drilled" 2 "After Drilled") size(small)) scheme(s1mono) xlabel(, labsize(small)) ylabel(, labsize(small)) yscale(range(-1 .5))  ylabel(-1(.5).5) 
graph save $graph/LRGW, replace
graph export $graph/LRGW.eps, replace

twoway (line yhat_before2p points, lcolor(black) lpattern(solid))  (line yhat_after2p points, lcolor(black) lpattern(dash)), /*
*/ xtitle("Distance from Well (in meters)", size(small)) note("Bandwidth= `band' meters, N1 = `n_beforep', N2 = `n_afterp'") ytitle("Log Price Residuals", size(small)) /*
*/ legend(order(1 "Before Drilled" 2 "After Drilled") size(small)) scheme(s1mono) xlabel(, labsize(small)) ylabel(, labsize(small)) yscale(range(-1 .5))  ylabel(-1(.5).5) 
graph save $graph/LRPWSA, replace
graph export $graph/LRPWSA.eps, replace


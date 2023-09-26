cap program drop fedistr

program fedistr 
version 15.0

*-------------------------------------------------------------------------------
* syntax and options
*-------------------------------------------------------------------------------

#delimit ;

syntax varlist(numeric min=1 max=9)	[if] [in] [aweight fweight] , [

Controls(varlist numeric fv) fe(varlist) by(varlist min=1 max=1)
HISTogram DENSity bins(passthru) width(passthru) bwidth(passthru) COMpare
Percentiles(numlist) nosd Nobs
STANDardize log LOGVar(varlist numeric) MEANCenter COMMONSample
COLorscheme(string) PLOTScheme(string) CINTensity(numlist max=1) 
scale(string) XYSize(string)
NOLEGend LEGPos(numlist) LEGSize(numlist) NONote
opts(string asis) grcopts(string asis)

] ;

#delimit cr


*-------------------------------------------------------------------------------
* check if options are correct & output errors
*-------------------------------------------------------------------------------

// variable absorbtion
if wordcount("`varlist'")>9 {
	di as error "The command only supports up to nine variables in {it:varlist}."
}
if "`fe'`controls'"=="" & "`compare'"!="" {
	di as error "The {it:compare} option requires either {it:fe()} or {it:controls()}."
	exit 498
}

// by 
if "`by'"!="" {
	if wordcount("`varlist'")>1 {
		di as error "The by() option only works when a single variable is specified in {it:varlist}."
		exit 498
	}
	qui duplicates report `by'
	local byvalcount = r(unique_value)
	if `byvalcount'>9 {
		di as error "The command only supports up to nine distinct levels of the by-variable."
		exit 498
	}
}

// variable transformations
if "`logvar'"!="" {
	foreach x of varlist `logvar' {
		if !strpos("`varlist'","`x'") {
			di as error "Variable {it:`x'} (specified in {it:logvar()}) does not appear as part of {it:fedistr(varlist)}."
			exit 498
		}
	}
}


*-------------------------------------------------------------------------------
* prep dependencies
*-------------------------------------------------------------------------------

// install
foreach package in reghdfe ftools {
	capture which `package'
	if _rc==111 ssc install `package', replace
}
capture which colorpalette
if _rc==111 {
	ssc install colrspace, replace
	ssc install palettes, replace
}
capture findfile blindschemes.sthlp
if _rc!=0 {
	capture set scheme plotplain
	if _rc==111 ssc install blindschemes, replace
}


*-------------------------------------------------------------------------------
* prep dataset
*-------------------------------------------------------------------------------

// preserve original data
preserve

// suppress output
quietly {

// weight local
if ("`weight'"!="") local w [`weight'`exp']
if ("`weight'"!="") local weightname = subinstr("`exp'","=","",.)

// make variables numeric if required
if "`fe'`by'"!="" {
	foreach x of varlist `fe' `by' {
		capture confirm numeric variable `x'
		if _rc {
			rename `x' `x'alt
			egen `x' = group(`x'alt)
		}
	}
}

// drop observations
if "`by'"!="" | "`controls'"!="" | "`fe'"!="" {
	foreach v of varlist `by' `controls' `fe' {
		local covdrop `covdrop' `v'
	}
	egen fedistr_vardrop1 = rowmiss(`covdrop')
	drop if fedistr_vardrop1>0
}
if "`if'"!="" keep `if'
if "`in'"!="" keep `in'
if "`weight'"!="" keep if !mi(`weightname')
if "`commonsample'"!="" {
	egen fedistr_vardrop2 = rowmiss(`varlist')
	drop if fedistr_vardrop2>0
}

// drop variables
if "`fe'"!="" local vlist `vlist' `fe'
if "`controls'"!="" local vlist `vlist' `controls'
if "`weight'"!="" local vlist `vlist' `weightname'
if "`by'"!="" local vlist `vlist' `by'
keep `varlist' `vlist'


*-------------------------------------------------------------------------------
* prep variables
*-------------------------------------------------------------------------------

// by variable
local isthereby = 0
if "`by'"!="" {
	duplicates report `by'
	local byvalcount = r(unique_value)
	if `byvalcount'!=1 {
		
		* adapt locals
		local isthereby = 1
		local byparen by(`by')
		
		* separate versions of main variables for by-levels
		levelsof `by', local(bynum)
		foreach lvl in `bynum' {
			gen `varlist'`lvl' = `varlist' if `by'==`lvl'
			local byvarlist `byvarlist' `varlist'`lvl'
		}
	}
}

// main variables
if `isthereby'==0 {
	local varcount = wordcount("`varlist'")
	tokenize `varlist'
}
else {
	local varcount = `byvalcount'
	tokenize `byvarlist'
	local count = 1
	foreach n of numlist `bynum' {
		local bycode`count' = `n'
		local ++count
	}
}


*-------------------------------------------------------------------------------
* residualization & scaling
*-------------------------------------------------------------------------------

// variable transformations 
if "`log'"!="" | ("`by'"!="" & "`logvar'"!="") {
	foreach var of numlist 1/`varcount' {
		count if ``var''<=0 
		if r(N)>0 dis as error "``var'' has values <=0, which cannot be log-transformed. The values are set to missing."
		replace ``var'' = ln(``var'')
	}
}
if "`logvar'"!="" & !("`by'"!="" & "`logvar'"!="") {
	foreach x of varlist `logvar' {
		count if `x'<=0 
		if r(N)>0 dis as error "`x' has values <=0, which cannot be log-transformed. The values are set to missing."
		replace `x' = ln(`x')
	}
}
if "`standardize'"!="" {
	if `isthereby'==1 {
		sum `varlist' `w'
		local sd_bycase = r(sd)
	}
	foreach var of numlist 1/`varcount' {
		sum ``var'' `w'
		if `isthereby'==0 replace ``var'' = (``var''-r(mean))/r(sd)
		if `isthereby'==1 replace ``var'' = (``var''-r(mean))/`sd_bycase'
	}
}

if "`fe'`controls'"!="" {

// residualize
	if "`fe'"=="" local hdfeabsorb noabsorb
	if "`fe'"!="" local hdfeabsorb absorb(`fe')
	foreach var of numlist 1/`varcount' {
		reghdfe ``var'' `controls' `w', `hdfeabsorb' res(``var''_res)
		gen ``var''_smple = e(sample)
	}

// rescale original or residualized variable
	if "`meancenter'"=="" {
		if "`compare'"!="" {
			foreach var of numlist 1/`varcount' {
				sum ``var'' `w', meanonly
				replace ``var'' = ``var''-r(mean)
			}
		}
	}
	else {
		foreach var of numlist 1/`varcount' {
			sum ``var'' `w', meanonly
			replace ``var''_res = ``var''_res+r(mean)
		}
	}
}

// duplicate _res variable if nothing is absorbed 
else {
	foreach var of numlist 1/`varcount' {
		clonevar ``var''_res = ``var''
		gen ``var''_smple = 1 if !mi(``var'')
	}
}


*-------------------------------------------------------------------------------
* gather parameters
*-------------------------------------------------------------------------------

// store standard deviations
foreach var of numlist 1/`varcount' {
	sum ``var'' `w'
	local ``var''_sd = r(sd)
	sum ``var''_res if ``var''_smple==1 `w'
	local ``var''_ressd = r(sd)
	local ``var''_sdred = round((1-(```var''_ressd'/```var''_sd'))*100)
}

// store percentiles
if "`percentiles'"!="" {
	
	// count number
	local pcount = wordcount("`percentiles'")
	
	// options for line style in plot
	local pctline lp(dash) lc(gs4) lw(thin)
	
	foreach var of numlist 1/`varcount' {
		
		// store percentiles 
		_pctile ``var''_res `w', p(`percentiles')
		foreach p of numlist 1/`pcount' {
			local ``var''_p`p' = r(r`p')
			local ``var''_pct ```var''_pct' ```var''_p`p''
			local ``var''_pctlist ```var''_pctlist' ``var''_p`p'	// for later rounding
		}
		
		// prep plot line options
		local ``var''_pct xline(```var''_pct', `pctline')
	}
}

// number of observations
if "`nobs'"!="" {
	foreach var of numlist 1/`varcount' {
		count if !mi(``var'')
		local nobs`var' = r(N)
	}
}

// round parameters
*local space " "
foreach var of numlist 1/`varcount' {
	foreach par in ``var''_sd ``var''_ressd ```var''_pctlist' {
		local smallround`par' = 0
		if ``par''>=10 | ``par''<=-10 {
			local `par'round "1"
		}
		else if ``par''>=1 | ``par''<=-1 {
			local `par'round ".1"
		}
		else {
			local roundcount = 0
			local `par'string = subinstr("``par''","-","",.)
			local `par'round ".0"
			foreach rr of numlist 2/6 {
				if substr("``par'string'",`rr',1)!="0" {
					local `par'round "``par'round'1"
					continue, break
				}
				else {
					local `par'round "``par'round'0"
					local roundcount = `roundcount'+1
				}
			}
		}
		cap if strpos("``par'string'","e") & ``par''>0 local smallround`par' = 1
		cap if strpos("``par'string'","e") & ``par''<0 local smallround`par' = -1
		local `par' = round(``par'',``par'round')
		if `smallround`par''==0 {
			local `par' "`space'=`space'``par''"
		}
		else if `smallround`par''==1 {
			local `par' "`space'<`space'.00001"
		}
		else {
			local `par' "`space'{&cong}`space'0"
		}
		if strpos("``par''","000000") & "``par''"!="`space'<`space'.00001" {
			foreach zz of numlist 1/9 {
				if substr("``par''",`zz',1)=="." local dotpos = `zz'
			}
			if "`dotpos'"!="" {
				if "``par'round'"=="1" local `par' = substr("``par''",1,`dotpos'-1)
				if "``par'round'"==".1" local `par' = substr("``par''",1,`dotpos'+1)
				if "``par'round'"==".01" local `par' = substr("``par''",1,`dotpos'+2)
				if "``par'round'"==".001" local `par' = substr("``par''",1,`dotpos'+3)
				if "``par'round'"==".0001" local `par' = substr("``par''",1,`dotpos'+4)
				if "``par'round'"==".00001" local `par' = substr("``par''",1,`dotpos'+5)
			}
		}
	}
}

// prep figure text for percentiles
if "`percentiles'"!="" {
	foreach word of local percentiles {
		local lastpct `word'
	}
	foreach var of numlist 1/`varcount' {
		local count = 1
		foreach p in `percentiles' {
			if "`p'"!="`lastpct'" local comma`var'`count' ,
			local ``var''pcttxt ```var''pcttxt' pct{sub:`p'}```var''_p`count''`comma`var'`count''
			local ++count
		}
	}
}


*-------------------------------------------------------------------------------
* figure size & combination options
*-------------------------------------------------------------------------------

// scale
if "`scale'"=="" {
	if `varcount'==1 {
		local scale 1.2
	}
	else {
		if `varcount'==2 local scale 2.15
		if inlist(`varcount',3,4) local scale 1.3
		if inlist(`varcount',5,6) local scale 1.6
		if inlist(`varcount',7,8,9) local scale 1.05
	}
}
else {
	if strpos("`scale'","*") {
		local scale2 = subinstr("`scale'","*","",.)
		local scale = `scale2'
	}
}

// relative size
if "`xysize'"=="" {
	if `varcount'==1 {
		local xysize = 1.9
	}
	else {
		if `varcount'==2 local xysize = 3
		if inlist(`varcount',3,4) local xysize = 1.8
		if inlist(`varcount',5,6) local xysize = 2
		if inlist(`varcount',7,8,9) local xysize = 1.3
	}
}
else {
	if strpos("`xysize'","*") {
		local xysize2 = subinstr("`xysize'","*","",.)
		local xysize = `xysize2'
	}
}
if `xysize'<=1 {
	local ysize = 100
	local xsize = 100*`xysize'
}
if `xysize'>1 {
	local xsize = 100
	local ysize = 100*(1/`xysize')
}

// graph combine options
if `varcount'>1 {
	if inlist(`varcount',2,3,4) local grcombopts col(2)
	if inlist(`varcount',5,6,7,8,9) local grcombopts col(3)
}

// put together
if `varcount'==1 {
	local singlesize scale(*`scale') xsize(`xsize') ysize(`ysize')
}
else {
	local grcombsize iscale(*`scale') xsize(`xsize') ysize(`ysize')
}


*-------------------------------------------------------------------------------
* color palette
*-------------------------------------------------------------------------------

if "`colorscheme'"=="" {
	local cpal `" "210 0 0" "49 113 166" "'
}
else {
	local cpal `colorscheme'
}
if "`cintensity'"=="" & "`colorscheme'"=="" local cpalo_high int(1.2)
if "`cintensity'"=="" & "`colorscheme'"=="" local cpalo_low int(.7)
if "`cintensity'"=="" & "`colorscheme'"!="" local cpalo_high int(1)
if "`cintensity'"=="" & "`colorscheme'"!="" local cpalo_low int(1)
if "`cintensity'"!="" local cpalo_low int(`cintensity')
if "`cintensity'"!="" local cpalo_high int(`cintensity')

colorpalette `cpal', `cpalo_high' nograph local(,prefix(c) nonames)
colorpalette `cpal', `cpalo_low' nograph local(,prefix(c) suffix(low) nonames)
colorpalette `cpal', int(2.5) nograph local(,prefix(c) suffix(dark) nonames)
foreach i of numlist 10 20 50 80 {
	colorpalette `cpal', `cpalo' op(`i') nograph local(,prefix(c) suffix(o`i') nonames) 
}


*-------------------------------------------------------------------------------
* options for look of figure
*-------------------------------------------------------------------------------

// legend
if "`compare'"!="" {
	if "`legpos'"=="" local legpos 1
	if "`legsize'"!="" {
		if strpos("`legsize'","*") {
			local legsize2 = subinstr("`legsize'","*","",.)
			local legsize = `legsize'
		}
	}
	else {
		if `varcount'==1 local legsize 1
		if `varcount'>1 local legsize .7
	}
	local legend legend(order(1 "Original" 2 "Residual") ring(0) keygap(*.5) symxsize(*.5) bcolor(none) lcolor(none) fcolor(none) bmargin(zero) pos(`legpos') size(*`legsize') nobox)
	if "`nolegend'"!="" local legend legend(off)
}

// note 
if "`nonote'"=="" {
	if "`fe'"!="" {
		foreach x in `fe' {
			local felist `felist' `x'
		}
		local fenote `""Absorbed fixed effects: `felist'""'
	}
	if "`controls'"!="" {
		foreach x in `controls' {
			local controlslist `controlslist' `x'
		}
		local controlsnote `""Absorbed covariates: `controlslist'""'
	}
	if `varcount'==1 {
		local note `"`fenote' `controlsnote'"'
	}
	else {
		local grcnote `"`fenote' `controlsnote'"'
	}
	if `varcount'==1 local notesize size(*.75)
	if inlist(`varcount',2) local notesize size(*1.1)
	if inlist(`varcount',3,4) local notesize size(*.75)
	if inlist(`varcount',5,6) local notesize size(*.8)
	if inlist(`varcount',7,8,9) local notesize size(*.45)
}

// histgram / density
if "`density'"=="" & "`histogram'"=="" local histogram histogram
if "`compare'"=="" {	// single plot
	if "`histogram'"!="" local look lalign(inside) lc(black) lw(thin) fc(`c1') 
	if "`density'"!="" local look lalign(outside) recast(area) lc(black) lw(thin) fc(`c1low')
}
if "`compare'"!="" {	// with original below
	if "`histogram'"!="" {
		local look lalign(inside) lc(`c1o10') lw(thin) fc(`c1o80')
		local oglook lalign(inside) lc(`c2o10') lw(thin) fc(`c2o80')
	}
	if "`density'"!="" {
		local look lalign(outside) recast(area) lc(`c1dark') lw(thin) fc(`c1o50')
		local oglook lalign(outside) recast(area) lc(`c2dark') lw(thin) fc(`c2o50')
	}
}

// title with variable name and parameters
foreach var of numlist 1/`varcount' {
	if "`log'"!="" | strpos("`logvar'","``var''") local ``var''log (log)
	if "`fe'`controls'"!="" local sdreduc {it:(-```var''_sdred'%)}
	if `isthereby'==0 local title1 {bf:``var''} ```var''log'
	if `isthereby'==1 local title1 {bf:`varlist' (at `by'==`bycode`var'')} ```var''log'
	local title2 SD{sub:original}```var''_sd', SD{sub:residual}```var''_ressd' `sdreduc'
	local title3 ```var''_sdred'% of variance absorbed
	if "`nobs'"!="" local title4 `" "N = `nobs`var''"  "'
	if "`percentiles'"!="" local title5 `" "```var''pcttxt'"  "'
	local ``var''_title title("`title1'" "`title2'" `title5' `title4')
}

// overall twoway options
if "`plotscheme'"=="" {
local axisopts labs(*1.05) tlc(gs2) tlw(thin) nogrid
local commonopts scheme(plotplain) ///
	ysc(off) xsc(lc(gs2) lw(thin)) ylab(,`axisopts') xlab(#6,`axisopts') ///
	title(,size(*1)) xtitle("") note(`note', `notesize') `legend' `singlesize'
}
else {
	local commonopts scheme(`plotscheme')
}
	
// graph combine options
local grcopts `grcombopts' imargin(small) note(`grcnote', `notesize')


*-------------------------------------------------------------------------------
* plot
*-------------------------------------------------------------------------------

// decide which plot to draw
if "`histogram'"!="" {
	local plottype hist
	local plottypeopts `bins' `width'
}
if "`density'"!="" {
	local plottype kdensity
	local plottypeopts `bwidth'
}
}	// closes quietly bracket

// draw plot 
foreach var of numlist 1/`varcount' {
	
	// options & compiler for graph combine
	if `varcount'>1 {
		local sav name(h`var', replace) nodraw
		local compiler `compiler' h`var'
	}
	
	// original distribution
	if "`compare'"!="" local plot`var' `plot`var'' (`plottype' ``var'', `oglook' `plottypeopts')
	
	// residualized distribution
	local plot`var' `plot`var'' (`plottype' ``var''_res, `look' `plottypeopts')
	
	// draw with graph twoway
	tw `plot`var'', `commonopts' ```var''_pct' ```var''_title' `sav' `opts'

}

// combine multiple plots 
if `varcount'>1 {
	graph combine `compiler', `grcombsize' `grcopts'
}


restore
end



















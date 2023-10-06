*! version 1.2   Leo Ahrens   leo@ahrensmail.de

program define fedistr 
version 15.0

*-------------------------------------------------------------------------------
* syntax and options
*-------------------------------------------------------------------------------

#delimit ;

syntax varlist(numeric min=1 max=9)	[if] [in] [aweight fweight] , [

Controls(varlist numeric fv) fe(varlist) by(varlist min=1 max=1)
HISTogram DENSity bins(passthru) width(passthru) bwidth(passthru) COMpare
Percentiles(numlist) nosd Nobs
STANDardize log LOGVar(varlist numeric) MEANcenter WINSor(numlist) WINSORRes(numlist) KEEPvars REPLace
COMMONSample DROPSingle DROPZerovar
COLorscheme(string) PLOTScheme(string) CINTensity(numlist max=1) 
scale(string) XYSize(string)
NOPlot NOLEGend LEGPos(numlist) LEGSize(numlist) NONote 
opts(string asis) grcopts(string asis) slow

] ;

#delimit cr



*-------------------------------------------------------------------------------
* prep slow options
*-------------------------------------------------------------------------------
	
if "`slow'"!="" {
	local lvlsof levelsof
	local sumd sum 
	local sumd2 , d
	local ggen egen
	local duplrep duplicates report
}
else {
	local lvlsof glevelsof
	local sumd gstats sum
	local ggen gegen
	local duplrep gdistinct
}

*-------------------------------------------------------------------------------
* check if options are correct & output errors
*-------------------------------------------------------------------------------

// suppress output
quietly {

// variable absorbtion
if wordcount("`varlist'")>9 {
	di as error "The command only supports up to nine variables in {it:varlist}."
}
if "`fe'`controls'`by'"=="" & "`compare'"!="" {
	di as error "The {it:compare} option requires either {it:fe()}, {it:controls()}, or {it:by()}."
	exit 498
}

// by 
if "`by'"!="" {
	if "`if'"!="" local ifhelp &
	if "`if'"=="" local ifhelp if
	`duplrep' `by' `if' `ifhelp' !mi(`by')
	if "`slow'"=="" local byvalcount = r(ndistinct)
	if "`slow'"!="" local byvalcount = r(unique_value)
	if `byvalcount'>7 {
		di as error "The command only supports up to seven distinct levels of the by-variable."
		exit 498
	}
	if "`p'"!="" {
		di as error "The by() option cannot be combined with p()."
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
if "`slow'"=="" {
	foreach package in reghdfe ftools {
		capture which `package'
		if _rc==111 ssc install `package', replace
	}
	capture which gtools
	if _rc==111 {
		ssc install gtools, replace 
		gtools, upgrade
	}
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

// identify observations for later merge
if "`keepvars'"!=""	gen fedistr_n = _n
	
// preserve original data
preserve

// weight local
if ("`weight'"!="") local w [`weight'`exp']
if ("`weight'"!="") local weightname = subinstr("`exp'","=","",.)

// make variables numeric if required
if "`fe'`by'"!="" {
	local numcount = 1
	foreach x of varlist `fe' `by' {
		capture confirm numeric variable `x'
		if _rc {
			rename `x' fedistr_old`numcount'
			`ggen' `x' = group(fedistr_old`numcount')
			local ++numcount
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
if "`keepvars'"!="" local vlist `vlist' fedistr_n
keep `varlist' `vlist'

// possibly drop obs again while retaining vars
if "`fe'`by'"!="" {
	local numcount = 1
	foreach var of varlist `varlist' {
		`ggen' fedistr_wsd`numcount' = sd(`var'), by(`fe' `by')
		if "`dropzerovar'"!="" drop if fedistr_wsd`numcount'==0
		`ggen' fedistr_c`numcount' = count(`var'), by(`fe' `by')
		if "`dropsingle'"!="" drop if fedistr_c`numcount'==1
		local ++numcount
	}
}


*-------------------------------------------------------------------------------
* prep variables
*-------------------------------------------------------------------------------

// by variable
local isthereby = 0
if "`by'"!="" {
	`duplrep' `by'
	if "`slow'"=="" local byvalcount = r(ndistinct)
	if "`slow'"!="" local byvalcount = r(unique_value)
	if `byvalcount'!=1 local isthereby = 1
}
if `isthereby'==0 {
	gen fedistr_by = 1
	local by fedistr_by
}
`lvlsof' `by', local(bylvls)
foreach lvl of numlist `bylvls' {
	local lastby `lvl'
}

// main variables
local varcount = wordcount("`varlist'")
tokenize `varlist'

*-------------------------------------------------------------------------------
* check for singletons and no within-variance
*-------------------------------------------------------------------------------

if "`fe'"!="" | `isthereby'==1 {
	if "`fe'"!="" & `isthereby'==0 local within within `fe'
	if "`fe'"!="" & `isthereby'==1 local within within `fe' `by'
	if "`fe'"=="" & `isthereby'==1 local within within `by'
	foreach var of numlist 1/`varcount' {
		count if fedistr_wsd`var'==0 & !mi(``var'')
		if r(N)>0 di as error r(N) " observations of ``var'' do not vary `within'. Consider dropping them from FE regressions."
		count if fedistr_c`var'==1 & !mi(``var'')
		if r(N)>0 di as error "There are " r(N) " singleton observations of ``var'' `within'. Consider dropping them from FE regressions."
	}
}

*-------------------------------------------------------------------------------
* residualization & scaling
*-------------------------------------------------------------------------------

// variable transformations 
if "`winsor'"!="" {
	foreach var of numlist 1/`varcount' {
		foreach pctl of numlist `winsor' {
			if "`slow'"=="" gquantiles ``var'' `w', _pctile percentiles(`pctl')
			if "`slow'"!="" _pctile ``var'' `w', p(`pctl')
			if `pctl'>=50 replace ``var'' = r(r1) if ``var''>r(r1) & !mi(``var'')
			if `pctl'<50 replace ``var'' = r(r1) if ``var''<r(r1) & !mi(``var'')
		}
	}
}
if "`log'"!="" {
	foreach var of numlist 1/`varcount' {
		count if ``var''<=0 
		if r(N)>0 dis as error "``var'' has values <=0, which cannot be log-transformed. The values are set to missing."
		replace ``var'' = ln(``var'')
	}
}
if "`logvar'"!="" & "`log'"=="" {
	foreach x of varlist `logvar' {
		count if `x'<=0 
		if r(N)>0 dis as error "`x' has values <=0, which cannot be log-transformed. The values are set to missing."
		replace `x' = ln(`x')
	}
}
if "`standardize'"!="" {
	foreach var of numlist 1/`varcount' {
		if "`slow'"!="" {
			sum ``var'' `w'
			replace ``var'' = (``var''-r(mean))/r(sd)
		}
		if "`slow'"=="" gstats transform (standardize) ``var'' `w', replace
	}
}

// residualize
if "`fe'`controls'"!="" | `isthereby'==1 {
	if "`slow'"!="" {
		foreach fevar in `fe' {
			local fevarlist `fevarlist' i.`fevar'
		}
		foreach var of numlist 1/`varcount' {
			reg ``var'' `fevarlist' `controls' `w'
			predict fedistr_res`var', resid
			gen fedistr_smple`var' = e(sample)
		}
	}
	else if "`controls'"=="" & `isthereby'==0 {
		global GTOOLS_BETA = "I KNOW WHAT I AM DOING"
		if `isthereby'==1 local gtoolsby by(`by')
		foreach var of numlist 1/`varcount' {
			gstats residualize ``var'' `w', absorb(`fe') `gtoolsby' gen(fedistr_res`var')
			gen fedistr_smple`var' = !mi(fedistr_res`var')
		}		
	}
	else {
		if "`fe'"=="" local hdfeabsorb noabsorb
		if "`fe'"!="" | `isthereby'==1 local hdfeabsorb absorb(`fe' `by')
		foreach var of numlist 1/`varcount' {
			reghdfe ``var'' `controls' `w', `hdfeabsorb' res(fedistr_res`var')
			gen fedistr_smple`var' = e(sample)
		}
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
			replace fedistr_res`var' = fedistr_res`var'+r(mean)
		}
	}
}

// winsor residualized
if "`winsorres'"!="" {
	foreach var of numlist 1/`varcount' {
		foreach pctl of numlist `winsorres' {
			if "`slow'"=="" gquantiles fedistr_res`var' `w', _pctile percentiles(`pctl')
			if "`slow'"!="" _pctile fedistr_res`var' `w', p(`pctl')
			if `pctl'>=50 replace fedistr_res`var' = r(r1) if fedistr_res`var'>r(r1) & !mi(fedistr_res`var')
			if `pctl'<50 replace fedistr_res`var' = r(r1) if fedistr_res`var'<r(r1) & !mi(fedistr_res`var')
		}
	}
}

// duplicate _res variable if nothing is absorbed 
if !("`fe'`controls'"!="" | `isthereby'==1) {
	foreach var of numlist 1/`varcount' {
		clonevar fedistr_res`var' = ``var''
		gen fedistr_smple`var' = 1 if !mi(``var'')
	}
}


*-------------------------------------------------------------------------------
* gather parameters
*-------------------------------------------------------------------------------

if "`noplot'"=="" {
	
// store standard deviations
foreach var of numlist 1/`varcount' {
	sum ``var'' `w'
	local sd_var`var' = r(sd)
	local count = 1
	foreach lvl of numlist `bylvls' {
		sum fedistr_res`var' if fedistr_smple`var'==1 & `by'==`lvl' `w'
		local ressd_var`var'`count' = r(sd)
		local sdred_var`var'`count' = round((1-(`ressd_var`var'`count''/`sd_var`var''))*100)
		local sdlist_var`var' `sdlist_var`var'' ressd_var`var'`count'
		local ++count
	}
}

// store percentiles
if "`percentiles'"!="" {
	
	// count number
	local pcount = wordcount("`percentiles'")
	
	// options for line style in plot
	local pctline lp(dash) lc(gs4) lw(thin)
	
	foreach var of numlist 1/`varcount' {
		
		// store percentiles 
		if "`slow'"=="" gquantiles fedistr_res`var' `w', _pctile percentiles(`percentiles')
		if "`slow'"!="" _pctile fedistr_res`var' `w', p(`percentiles')
		foreach p of numlist 1/`pcount' {
			local p`p'_var`var' = r(r`p')
			local pct_var`var' `pct_var`var'' `p`p'_var`var''
			local pctlist_var`var' `pctlist_var`var'' p`p'_var`var'	// for later rounding
		}
		
		// prep plot line options
		local pct_var`var' xline(`pct_var`var'', `pctline')
	}
}

// number of observations
if "`nobs'"!="" {
	foreach var of numlist 1/`varcount' {
		local count = 1
		foreach lvl of numlist `bylvls' {
			count if !mi(``var'') & `by'==`lvl'
			local nobs`var'`count' = r(N)
			local ++count
		}
	}
}

// round parameters
*local space " "
foreach var of numlist 1/`varcount' {
	foreach par in sd_var`var' `sdlist_var`var'' `pctlist_var`var'' {
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
			local pcttxt_var`var' `pcttxt_var`var'' pct{sub:`p'}`p`count'_var`var''`comma`var'`count''
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
	local ysize = (100)/5
	local xsize = (100*`xysize')/5
}
if `xysize'>1 {
	local xsize = (100)/5
	local ysize = (100*(1/`xysize'))/5
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
	if `isthereby'==1 local cpal `" "210 0 0" "49 113 166" "15 137 1" "255 127 14" "169 58 228" "41 217 231" "250 238 22"  "222 115 50" "'
	if `isthereby'==0 local cpal `" "210 0 0" "49 113 166" "'
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
if "`compare'"!="" | `isthereby'==1 {
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
	if `isthereby'==0 {
		local legcont `"1 "Original" 2 "Residual""'
	}
	else {
		if "`compare'"!="" local legcont `"1 "Original""'
		local count = 1
		if "`compare'"!="" local ++count
		foreach lvl of numlist `bylvls' {
			local legcont `"`legcont' `count' "`by'=`lvl'""'
			local ++count
		}
	}
	local legend legend(order(`legcont') ring(0) keygap(*.5) symxsize(*.5) bcolor(none) lcolor(none) fcolor(none) bmargin(zero) pos(`legpos') size(*`legsize') nobox)
	if "`nolegend'"!="" local legend legend(off)
}
*	if `isthereby'==1 local title1 {bf:`varlist' (at `by'==`bycode`var'')} `log_var`var''
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
if "`compare'"=="" & `isthereby'==0 {
	if "`histogram'"!="" local look1 lalign(inside) lc(black) lw(thin) fc(`c1') 
	if "`density'"!="" local look1 lalign(outside) recast(area) lc(black) lw(thin) fc(`c1low')
}
else {
	local count = 1
	foreach byc of numlist `bylvls' {
		if "`histogram'"!="" local look`count' lalign(inside) lc(`c`count'o10') lw(thin) fc(`c`count'o80')
		if "`density'"!="" local look`count' lalign(outside) recast(area) lc(`c`count'dark') lw(thin) fc(`c`count'o50')
		local ++count
	}
	if "`compare'"!="" {
		if "`histogram'"!="" {
			if `isthereby'==0 local oglook lalign(inside) lc(`c2o10') lw(thin) fc(`c2o80')
			if `isthereby'==1 local oglook lalign(inside) lc(black%10) lw(thin) fc(black%70)
		}
		if "`density'"!="" {
			if `isthereby'==0 local oglook lalign(outside) recast(area) lc(`c2dark') lw(thin) fc(`c2o50')
			if `isthereby'==1 local oglook lalign(outside) recast(area) lc(black%95) lw(thin) fc(black%50)
		}
	}
}

// title with variable name and parameters
foreach var of numlist 1/`varcount' {
	* prep
	if "`log'"!="" | strpos("`logvar'","``var''") local log_var`var' (log)
	if "`fe'`controls'"!="" & `isthereby'==0 local sdreduc {it:(-`sdred_var`var'1'%)}
	* top title with variable name
	local title1_`var' {bf:``var''} `log_var`var''
	* standard deviations
	if `isthereby'==0 {
		if "`fe'`controls'"!="" local restitle , SD{sub:residual}`ressd_var`var'1' `sdreduc'
		local title2_`var' SD{sub:original}`sd_var`var''`restitle'
	}
	if `isthereby'==1 {
		local title2pre_`var' `""SD{sub:original}`sd_var`var''""'
		local count = 1
		foreach lvl of numlist `bylvls' {
			local comma 
			if `lvl'!=`lastby' local comma ,
			local title2_`var' `title2_`var'' SD{sub:`by'=`lvl'}`ressd_var`var'`count''`comma'
			local ++count
		}
	}
	* no of observations
	if "`nobs'"!="" {
		local count = 1
		foreach lvl of numlist `bylvls' {
			local comma 
			if `lvl'!=`lastby' local comma ,
			if `isthereby'==1 local sub {sub:`by'=`lvl'}
			local nobslist_`var' `nobslist_`var'' N`sub'=`nobs`var'`count''`comma'
			local ++count
		}
		local title4_`var' `" "`nobslist_`var''"  "'
	}
	* percentiles
	if "`percentiles'"!="" local title5_`var' `" "`pcttxt_var`var''"  "'
	* combine
	local title_var`var' title("`title1_`var''" `title2pre_`var'' "`title2_`var''" `title5_`var'' `title4_`var'')
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
} 	// closes plotfig brackket
}	// closes quietly bracket
if "`noplot'"=="" {

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
	local count = 1
	foreach lvl of numlist `bylvls' {
		local plot`var' `plot`var'' (`plottype' fedistr_res`var' if `by'==`lvl', `look`count'' `plottypeopts')
		local ++count
	}
		
	// draw with graph twoway
	tw `plot`var'', `commonopts' `pct_var`var'' `title_var`var'' `sav' `opts'

}

// combine multiple plots 
if `varcount'>1 graph combine `compiler', `grcombsize' `grcopts'

}

*-------------------------------------------------------------------------------
* merge back stored variables
*-------------------------------------------------------------------------------

quietly {

	if "`keepvars'"!="" {
		foreach var of numlist 1/`varcount' {
			if "`replace'"=="" {
				rename fedistr_res`var' ``var''_res
				local keeplist `keeplist' ``var''_res
			}
			if "`replace'"!="" {
				cap drop ``var''
				rename fedistr_res`var' ``var''
				local keeplist `keeplist' ``var''
			}
		}
		keep fedistr_n `keeplist'
		tempfile fedistr_merge
		save `fedistr_merge', replace
	}

	restore

	if "`keepvars'"!="" {
		if "`replace'"!="" drop `keeplist'
		merge 1:1 fedistr_n using `fedistr_merge', nogen update replace
		drop fedistr_n
	}
}


end




















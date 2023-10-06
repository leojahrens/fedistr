{smcl}
{title:help fedistr}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:fedistr}}Visualize distributions after absorbing fixed effects and covariates{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:fedistr}
varlist {ifin}
{weight} {cmd:,} [options] {p_end}

{marker opt_summary}{...}
{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}

{syntab:Absorbed variables}
{synopt :{opt fe(varlist)}}Fixed effects variables and categorical controls.{p_end}
{synopt :{opt c:ontrols(varlist)}}Continuous control variables (used as linear controls).{p_end}

{syntab:Stratify into subsamples}
{synopt :{opt by(var)}}Plots separate distributions for each level of the specified variable.{p_end}

{syntab:Plotted distribution(s)}
{synopt :{opt com:pare}}Additionally plots the original distributions before absorbing FEs/covariates.{p_end}
{synopt :{opt hist:ogram}}Plots histograms (the default).{p_end}
{synopt :{opt dens:ity}}Plots kernel densities.{p_end}
{synopt :{opt bins(num)}}Specifies the number of bins of histograms.{p_end}
{synopt :{opt width(num)}}Specifies the width of histogram bars.{p_end}
{synopt :{opt bwidth(num)}}Specifies the bandwidth of density functions.{p_end}

{syntab:Additonal parameters}
{synopt :{opt p:ercentiles(numlist)}}Includes information on specified percentiles {it:p} of the distributions (0<{it:p}<100).{p_end}
{synopt :{opt n:obs}}Includes the number of non-missing observations.{p_end}

{syntab:Variable transformations}
{synopt :{opt stand:ardize}}Bring all variables on a scale with a mean of zero and a standard deviation of one (applied before residualization).{p_end}
{synopt :{opt log}}Use the natural logarithms of all variables (applied before residualization).{p_end}
{synopt :{opt logv:var(varlist)}}Only use logged values of specified variables.{p_end}
{synopt :{opt wins:or(numlist)}}Winsorize the variables at the specified percentiles {it:p} before residualization (0<{it:p}<100).{p_end}
{synopt :{opt winsorr:es(numlist)}}Winsorize at the specified percentiles after residualization.{p_end}
{synopt :{opt mean:center}}Add back the means to all variables after residualization.{p_end}
{synopt :}{it:Note: the transformations are implemented in the following order: winsorize original variable, log, standardize, residualize, add mean back, winsor residualized}.{p_end}

{syntab:Sample manipulation}
{synopt :{opt commons:ample}}Ensures a common sample across all plots by dropping all observations with missings on any of the specified variables.{p_end}
{synopt :{opt drops:ingle}}Drop singleton observations, i.e. observations with n=1 within the specified fixed effects and/or {it:by}.{p_end}
{synopt :{opt dropz:erovar}}Drop observations that do not vary within the specified fixed effects and/or {it:by}.{p_end}

{syntab:Keep variables in dataset}
{synopt :{opt keep:vars}}Retain the residualized and possibly transformed variables in the dataset (saved with suffix {it:_res}).{p_end}
{synopt :{opt repl:ace}}Overwrite the existing version of variables instead of generating a new variable.{p_end}

{syntab:Scheme and colors}
{synopt :{opt plots:cheme(str)}}Defines an alternative graph scheme, such as {opt plotscheme(white_tableau)}. See the collection in https://github.com/asjadnaqvi/stata-schemepack.{p_end}
{synopt :{opt col:orscheme(str)}}Defines a custom color palette (e.g. {opt colorscheme(tableau)}). To define your own colors, use a list of hex colors (e.g. {opt colorscheme(#E6007E #009FE3 #312783)}).{p_end}
{synopt :{opt cint:ensity(num)}}Changes the color intensity. Higher values make colors darker.{p_end}

{syntab:Plot and element size}
{synopt :{opt scale(num)}}Resizes all elements in the plot. For example, {opt scale(1.1)} increases text, marker, and lines sizes by 10%.{p_end}
{synopt :{opt xys:ize(num)}}Specifies the relative width to height. For example, {opt xysize(1.2)} draws a plot with 20% more width than height.{p_end}
{synopt :{opt legs:ize(num)}}Resizes the legend.{p_end}

{syntab:Other}
{synopt :{opt noleg:end}}Omits the legend when the original distributions are included via {opt compare}.{p_end}
{synopt :{opt legp:os(num)}}Position of the legend (enter values 1-12, akin to the positions of clock hands).{p_end}
{synopt :{opt opts(str)}}Passes on options to the twoway plot(s). For example, {opt opts(xtitle("Example"))} changes the x axis title. See {opt help twoway}.{p_end}
{synopt :{opt grcopts(str)}}Passes on options to graph combine, which is used internally to combine twoway plots when more than one distribution is plotted. See {opt help graph combine}.{p_end}

{hline}

{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. sysuse lifeexp, clear}{p_end}

{pstd}Simple case - one variable, one fixed effect{p_end}
{phang2}{cmd:. fedistr gnppc, fe(region)}{p_end}

{pstd}Compare to the original distribution{p_end}
{phang2}{cmd:. fedistr gnppc, fe(region) compare}{p_end}

{pstd}Multiple variables{p_end}
{phang2}{cmd:. fedistr gnppc lexp, fe(region) compare}{p_end}

{pstd}Include continuous control variiables{p_end}
{phang2}{cmd:. fedistr gnppc lexp, fe(region) controls(trunk) compare}{p_end}

{pstd}Use densities instead of histograms{p_end}
{phang2}{cmd:. fedistr gnppc lexp, fe(region) density compare}{p_end}

{pstd}Include percentiles{p_end}
{phang2}{cmd:. fedistr gnppc lexp, fe(region) p(1 99) compare}{p_end}

{pstd}Separate plots for subsamples{p_end}
{phang2}{cmd:. fedistr gnppc, fe(region) by(foreign) compare}{p_end}



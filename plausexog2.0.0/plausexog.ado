*! plausexog: Estimating bounds with a plausibly exogenous exclusion restriction  
*! Version 2.0.0 junio 7, 2016 @ 10:12:45
*! Author: Damian Clarke (application of code and ideas of Conley et al., 2012)
*! Much of the heart of this code comes from the Conley et al implementation
*! Contact: damian.clarke@economics.ox.ac.uk

/*
version highlights:
0.1.0: UCI and LTZ working.  No graphing
0.2.1: Graphing added
1.0.0: Completed beta testing.  Output complete. SSC v 1.0.0
1.0.1: Minor change to graph label options.
1.1.0: Weighting incorporated
1.2.0: Bug fix: very long names passed through syntax
2.2.0: Now Allowing for all arbitrary distributions with simulation algorithm
*/

cap program drop plausexog
program plausexog, eclass
version 8.0
#delimit ;

syntax anything(name=0 id="variable list")
	[if] [in]
	[fweight pweight aweight iweight]
	[,
	grid(real 2)
	gmin(numlist)
	gmax(numlist)
	level(real 0.95)
	omega(string)
 	mu(string)
	GRAph(varlist)
	GRAPHOMega(namelist min=2 max=22)
	graphmu(namelist min=2 max=22)
	graphdelta(numlist)
	graphopts(string)
    VCE(string)
    DISTribution(string)
    seed(numlist min=1 max=1)
    iterations(integer 1000) 
	]
	;
#delimit cr


********************************************************************************
*** (1) Unpack arguments, check for valid syntax, general error capture
********************************************************************************
local 0: subinstr local 0 "=" " = ", count(local equal)
local 0: subinstr local 0 "(" " ( ", count(local lparen)
local 0: subinstr local 0 ")" " ) ", count(local rparen)

tokenize `0'

local method `1'
macro shift
dis "`0'"
	
if "`method'"!="uci"&"`method'"!="ltz"&"`method'"!="upwci" {
	dis as error "Method of estimation must be specified."
	dis "Re-specify using uci, ltz or upcwi (see help file for more detail)"
	exit 200
}

	
if `equal'!=1|`lparen'!=1|`rparen'!=1 {
	dis as error "Specification of varlist is incorrect."
	dis as error "Ensure that synatx is: method yvar [exog] (endog=iv), [opts]"	
	exit 200
}
	
local yvar `1'
macro shift

local varlist1
while regexm(`"`1'"', "\(")==0 {
	local varlist1 `varlist1' `1'
	macro shift
}

local varlist2
while regexm(`"`1'"', "=")==0 {
	local var=subinstr(`"`1'"', "(", "", 1)
	local varlist2 `varlist2' `var'
	macro shift
}	

local varlist_iv
while regexm(`"`1'"', "\)")==0 {
	local var=subinstr(`"`1'"', "=", "", 1)
	local varlist_iv `varlist_iv' `var'
	macro shift
}

foreach list in varlist1 varlist2 varlist_iv {
	fvexpand ``list''
	local `list' `r(varlist)'
}

local allout `varlist1' `varlist2' constant
local allexog `varlist1' `varlist_iv' constant
	
local count2     : word count `varlist2'
local count_iv	 : word count `varlist_iv' 
local count_all  : word count `allout'
local count_exog : word count `allexog'
local countmin   : word count `gmin'
local countmax   : word count `gmax'
	
if `count2'>`count_iv' {
	dis as error "Specify at least as many instruments as endogenous variables"
	exit 200	
}
	
if "`method'"=="uci" {
	if `countmin'!=`count_iv'|`countmax'!=`count_iv' {
		dis as error "You must define as many gamma values as instrumental variables"
		dis "If instruments are believed to be valid, specify gamma=0 for gmin and gmax"
		exit 200	
	}

	foreach item in min max {
		local count=1
		foreach num of numlist `g`item'' {
			local g`count'`item'=`num'
			local ++count
		}
	}
}
if "`method'"=="ltz" & length("`distribution'")==0 {
	if length("`omega'")==0|length("`mu'")==0 {
		dis as error "For ltz, omega and mu matrices must be defined"
		exit 200
	}
	else {
		mat def omega_in=`omega'
		mat def mu_in=`mu'
	}
}

if length("`distribution'")!=0 {
    if "`method'"!="ltz" {
        dis as error "The distribution option can only be specified with ltz"
        error 200
    }
    local distribution: subinstr local distribution "," " , ", all
    local distcnt : list sizeof distribution

    local jj=1
    foreach j of numlist 1(1)`distcnt' {
        local dist`jj': word `j' of `distribution'
        local dist`jj': subinstr local dist`jj' "," ""        
        if "`dist`jj''"!="" local ++jj
    }
    
    local derr1  "If specifying a distribution with"
    local derr2  "parameter must be specified"
    local derr2s "parameters must be specified"
    local accept "normal, uniform, chi2, poisson, t, gamma, special"
    if "`dist1'"=="normal" {
        if `jj'!=4 {
            dis as error "`derr1' normal, 2 `derr2s' (mean and standard deviation)."
            exit 200
        }
        local gammaCall rnormal(`dist2', `dist3')
    }
    else if "`dist1'"=="uniform" {
        if `jj'!=4 {
            dis as error "`derr1' uniform, 2 `derr2s' (minimum and maximum)."
            exit 200
        }
        local gammaCall `dist2'+(`dist3'-`dist2')*runiform()
    }
    else if "`dist1'"=="chi2" {
        if `jj'!=3 {
            dis as error "`derr1' chi2, 1 `derr2' (degrees of freedom)."
            exit 200
        }
        if `dist2'<1 {
            dis as error "At least 1 degree of freedom must be specified for chi2"
            exit 200
        }
        local gammaCall rchi2(`dist2')        
    }
    else if "`dist1'"=="poisson" {
        if `jj'!=3 {
            dis as error "`derr1' poisson, 1 `derr2' (distribution mean)."
            exit 200
        }
        if `dist2'<1 {
            dis as error "At least a mean of 1 must be specified for poisson"
            exit 200
        }
        local gammaCall rpoisson(`dist2')                
    }
    else if "`dist1'"=="t" {
        if `jj'!=3 {
            dis as error "`derr1' t, 1 `derr2' (degrees of freedom)."
            exit 200
        }
        if `dist2'<1 {
            dis as error "At least 1 degree of freedom must be specified for t"
            exit 200
        }
        local gammaCall rt(`dist2')
    }
    else if "`dist1'"=="gamma" {
        if `jj'!=4 {
            dis as error "`derr1' gamma, 2 `derr2s' (shape and scale)."
            exit 200
        }
        if `dist2'<=0|`dist3'<=0 {
            dis as error "The shape and scale parameter for gamma must be > 0"
            exit 200
        }
        local gammaCall rgamma(`dist2',`dist3')
    }
    else if "`dist1'"=="special" {
        local sperr "To define your own distribution you must specify"
        if `jj'!=3 {
            dis as error "`sperr' one valid variable with the empirical distribution"
            exit 200
        }
        cap sum `dist2'
        if _rc!=0 {
            dis as error "`sperr' a valid variable with the empirical distribution"
            exit 200
        }
    }
    else {
        dis as error "The distribution option can only specify: `accept'"
        error 200
    }
}

dis "Estimating Conely et al.'s `method' method"
dis "Exogenous variables: `varlist1'"
dis "Endogenous variables: `varlist2'"
dis "Instruments: `varlist_iv'"


********************************************************************************
*** (2) Estimate model under assumption of gamma=0
********************************************************************************
local cEx   : word count `varlist1'
local cEn   : word count `varlist2'
local cIV   : word count `varlist_iv'

*dis "Trying to run intial reg with `cEx' exog vars `cEn' endog and `cIV' insts"	
qui ivregress 2sls `yvar' `varlist1' (`varlist2'=`varlist_iv') `if' `in'  /*
 */ [`weight' `exp'], vce(`vce')
qui estimates store __iv
	
********************************************************************************
*** (3) Union of Confidence Intervals approach (uci)
***     Here we are creating a grid and testing for each possible gamma combo:
***     ie - {g1min,g2min,g3min}, {g1max,g2min,g3min}, ..., {g1max,g2max,g3max}
***     This conserves much of the original (Conley et al.) code, which does 
***	  this in quite a nice way.
********************************************************************************
if "`method'"=="uci" {

	local points=1
	if length("`graph'")!=0 {
		local points=4
		matrix __graphmat = J(`points',4,.)
	}		
		
	foreach gnum of numlist 1(1)`points' {
		local cRatio =	`gnum'/`points'
		foreach item in min max {
		local count=1
		foreach num of numlist `g`item'' {
			local g`count'`item'l=`num'*`cRatio'
			local ++count
			}
		}
		**************************************************************************
		*** (3a) Iterate over iter, which is each possible combination of gammas
		**************************************************************************
		local iter=1
		while `iter' <= (`grid'^`count_iv') {
			local R=`iter'-1
			local w=`count_iv'

			**Create weighting factor to grid gamma.  If grid==2, gamma={max,min}
			while `w'>0 {
				local a`w'     = floor(`R'/(`grid'^(`w'-1)))
				local R        = `R'-(`grid'^(`w'-1))*`a`w''
				local gamma`w' = `g`w'minl' + ((`g`w'maxl'-`g`w'minl')/(`grid'-1))*`a`w''
				local --w
			}
				
			tempvar Y_G
			qui gen `Y_G'=`yvar'

				
			local count=1
			foreach Z of local varlist_iv {
				qui replace `Y_G'=`Y_G'-`Z'*`gamma`count''
				local ++count
			}

			***********************************************************************
			*** (3b) Estimate model based on assumed gammas, memoize conf intervals
			***********************************************************************
			qui ivregress 2sls `Y_G' `varlist1' (`varlist2'=`varlist_iv') `if' /*
			*/ `in' [`weight' `exp'], vce(`vce')
				
			***********************************************************************
			*** (3c) Check if variable is not dropped (ie dummies) and results
			***********************************************************************
			mat b2SLS   = e(b)
			mat cov2SLS = e(V)

			local vars_final
			local counter=0
			foreach item in `e(exogr)' `e(instd)' _cons {
				if _b[`item']!=0|_se[`item']!=0 {
					local vars_final `vars_final' `item'
					local ++counter
				}
			}
			ereturn scalar numvars = `counter'
			
			mat b2SLSf   = J(1,`counter',.)
			mat se2SLSf =	J(`counter',`counter',.)
			tokenize `vars_final'

			foreach num of numlist 1(1)`counter' {
				mat b2SLSf[1,`num']=_b[``num'']
				mat se2SLSf[`num',`num']=_se[``num'']
			}
			mat CI    = -invnormal((1 - `level')/2)
			mat ltemp = vec(b2SLSf) - CI*vec(vecdiag(se2SLSf))
			mat utemp = vec(b2SLSf) + CI*vec(vecdiag(se2SLSf))

			***********************************************************************
			*** (3d) Check if CI from model is lowest/highest in union (so far)
			***********************************************************************
			foreach regressor of numlist 1(1)`counter' {
				if `iter'==1 {
					local l`regressor'=.
					local u`regressor'=.	
				}
				local l`regressor' = min(`l`regressor'',ltemp[`regressor',1])
				local u`regressor' = max(`u`regressor'',utemp[`regressor',1])
			}
			local ++iter
		}

		if `gnum'==`points' {
			dis in yellow _newline
			dis "Conley et al (2012)'s UCI results" _col(55) "Number of obs =      " e(N)
			dis in yellow in smcl "{hline 78}"
			dis "Variable" _col(13) "Lower Bound" _col(29) "Upper Bound"
			dis in yellow in smcl "{hline 78}"

			tokenize `vars_final'

			foreach regressor of numlist 1(1)`counter' {
				dis in green "``regressor''" _col(13) `l`regressor'' _col(29) `u`regressor''
				foreach var of local varlist2 {
					if `"`var'"'==`"``regressor''"' {
						ereturn scalar lb_`var'=`l`regressor''
						ereturn scalar ub_`var'=`u`regressor''
					}
				}
			}
			dis in yellow in smcl "{hline 78}"
			if length("`graph'")!=0 {		
				matrix __graphmat[`gnum',1]=`cRatio'*`g1max'
				matrix __graphmat[`gnum',2]=.
				matrix __graphmat[`gnum',3]=e(lb_`graph')
				matrix __graphmat[`gnum',4]=e(ub_`graph')
			}
		}
		else if `gnum'<`points' {
			matrix __graphmat[`gnum',1]=`cRatio'*`g1max'
			tokenize `vars_final'	
			foreach regressor of numlist 1(1)`counter' {
				if `"`graph'"'==`"``regressor''"' {
					matrix __graphmat[`gnum',3]=`l`regressor''
					matrix __graphmat[`gnum',4]=`u`regressor''
				}
			}
		}
	}
}


********************************************************************************
*** (5) Local to Zero approach (ltz)
********************************************************************************    
if "`method'"=="ltz" {
	tempvar const
	qui gen `const'=1

	*****************************************************************************
	*** (5a) Remove any colinear elements to ensure that matrices are invertible
	*** For the case of the Z vector, this requires running the first stage regs	
	*****************************************************************************		
	unab testvars: `varlist1'		
	local usevars1
	foreach var of local testvars {
		cap dis _b[`var']
		if _rc!=0 continue
		if _b[`var']!=0 local usevars1 `usevars1' `var'
	}

	unab testvars: `varlist2'
	local usevars2
	foreach var of local testvars {
		cap dis _b[`var']
		if _rc!=0 continue
		if _b[`var']!=0 local usevars2 `usevars2' `var'
	}

	*****************************************************************************
	*** (5b) Form moment matrices: Z'X and Z'Z
	*****************************************************************************
	mat ZX = J(1,`count_all',.)
	mat ZZ = J(1,`count_exog',.)

	unab allvars: `varlist_iv' `usevars1' `const'
	tokenize `allvars'
	mat vecaccum a = `1' `usevars2' `usevars1' `if' `in'
	mat ZX = a
	while length("`2'")!= 0 {		
		mat vecaccum a = `2' `usevars2' `usevars1' `if' `in'
		mat ZX         = ZX\a
		macro shift		
	}

	tokenize `allvars'
	mat vecaccum a = `1' `varlist_iv' `usevars1' `if' `in'
	mat ZZ = a

	while length("`2'")!= 0 {	
		mat vecaccum a = `2' `varlist_iv' `usevars1' `if' `in'
		mat ZZ         = ZZ\a
		macro shift
	}

	scalar s1=rowsof(ZZ)
	scalar s2=rowsof(ZX)
	scalar s3=rowsof(omega_in)

	local em "conformability errors" 

	if s1!=s3 {
		dis as err "Z'Z matrix is `=s1'*`=s1', Omega defined by user is `=s3'*`=s3'"
		dis as err "Ensure that Omega is of the same dimension as Z'Z to avoid `em'"
	}
	if s2!=s3 {
		dis as err "Z'X matrix is `=s2'*`=s2', Omega defined by user is `=s3'*`=s3'"
		dis as err "Ensure that Omega is of the same dimension as Z'X to avoid `em'"
	}

    *****************************************************************************
	*** (5c) Form estimates if non-normal distribution
	*****************************************************************************
    if length("`distribution'")!=0 {

        if length("`seed'")!= 0 {
            set seed `seed'
        }
        if "`dist1'"=="special" {
            mkmat `dist2', matrix(specialgamma) nomissing
            matrix mnsp = rowsof(specialgamma)
            local nsp   = mnsp[1,1]
            matrix gammaDonors = J(`iterations',1,.)
            foreach gd of numlist 1(1)`iterations' {
                matrix gammaDonors[`gd',1]=specialgamma[ceil(runiform()*`nsp'),1]
            }
        }
        qui estimates restore __iv
        qui estat vce
        matrix varcovar = r(V)
        mata: st_matrix("betas", select(st_matrix("e(b)"), st_matrix("e(b)") :!=0))

        matrix mnvars = colsof(betas)
        local nvars   = mnvars[1,1]
        matrix gamma  = J(`nvars',1,0)
        matrix A      = inv(ZX'*inv(ZZ)*ZX)*ZX'
        foreach num of numlist 1(1)`nvars' {
             matrix betasSim = J(`iterations',`nvars',.)
        }

        qui count
        if r(N)>=`iterations' {
            dis "Simulating.  This may take a moment..."
            local simvars
            foreach num of numlist 1(1)`nvars' {
                tempvar sim`num'
                local simvars `simvars' `sim`num''
            }
            drawnorm `simvars', cov(varcovar) means(betas)

            local iter = 1
            while `iter' <= `iterations' {
                if "`dist1'"=="special" local gammaCall = gammaDonors[`iter',1]
                matrix gamma[1,1]=`gammaCall'
                matrix F = A*gamma
                foreach num of numlist 1(1)`nvars' {
                    qui sum `sim`num'' in `iter'
                    local input = r(mean)
                    matrix betasSim[`iter',`num'] = `input'+F[`num',1] 
                }
                local ++iter
            }
        }
        else {
            dis "Simulating.  This may take a moment..."
            local iter = 1
            while `iter' <= `iterations' {
                local simvars
                foreach num of numlist 1(1)`nvars' {
                    tempvar sim`num'
                    local simvars `simvars' `sim`num''
                }
                drawnorm `simvars', cov(varcovar) means(betas)
                if "`dist1'"=="special" local gammaCall = gammaDonors[`iter',1]
                matrix gamma[1,1]=`gammaCall'
                matrix F = A*gamma
                foreach num of numlist 1(1)`nvars' {
                    qui sum `sim`num'' in 1
                    local input = r(mean)
                    matrix betasSim[`iter',`num'] = `input'+F[`num',1] 
                }
                drop `simvars'
                local ++iter
            }
        }

        foreach num of numlist 1(1)`nvars' {
            mata : st_matrix("betasSim", sort(st_matrix("betasSim"), `num'))
            
            local lowerbound = betasSim[round(`iterations'*0.025),`num']
            local upperbound = betasSim[round(`iterations'*0.975),`num']
            dis "Bound for variable `num' is [`lowerbound',`upperbound']
        }
        exit

    }
    *****************************************************************************
	*** (5ci) Form augmented var-covar and coefficient matrices for graphing
	*****************************************************************************
	if length("`graph'")!=0 {
		local Nomegas  : word count `graphomega'	
		local Nmus     : word count `graphmu'

		if `Nomegas'==0|`Nmus'==0 {
			dis as error "When specifing graph and ltz, the graphomega and graphmu options are required"
			exit 272
		}
		if `Nomegas'!=`Nmus' {
			dis as error "graphmu and graphomega must take the same number of elements"
			exit 272
		}
		matrix __graphmatLTZ = J(`Nomegas',4,.)

		local countdelta   : word count `graphdelta'
		if `countdelta'!=0 {
			local j=1	
			foreach num of numlist `graphdelta' {
					matrix __graphmatLTZ[`j',1]=`num'
				local ++j
			}
		}

		tokenize `graphmu'
		local j=1
		foreach item in `graphomega' {
			qui dis "Estimating for `item'"
			mat def omegaC=`item'
			mat def muC=``j''

			scalar s4=rowsof(omegaC)
			scalar s5=colsof(omegaC)
			scalar s6=rowsof(muC)
			scalar s7=colsof(muC)

			if s1!=s4|s1!=s5 {
				dis as err "Z'Z matrix is `=s1'*`=s1', graph omega (`j') defined by user is `=s4'*`=s5'"
				dis as err "Ensure that Omega is of the same dimension as Z'Z to avoid `em'"
			}
			if s2!=s4|s2!=s5 {
				dis as err "Z'X matrix is `=s2'*`=s2', graph omega (`j') defined by user is `=s4'*`=s5'"
				dis as err "Ensure that Omega is of the same dimension as Z'X to avoid `em'"
			}
			if s7!=1 {
				dis as error "graph mu (`j') does not have 1 column"
				dis as err "Ensure that each graph mu is 1*k vector to avoid `em'"
			}
				
			qui estimates restore __iv

			mat Vc = e(V) +  inv(ZX'*inv(ZZ)*ZX)*ZX' * omegaC * ZX*inv(ZX'*inv(ZZ)*ZX)
			mat bc = e(b) - (inv(ZX'*inv(ZZ)*ZX)*ZX' * muC)'
				
			ereturn post bc Vc
			matrix CI = -invnormal((1-`level')/2)

			if `countdelta'==0 {
				scalar delta=omegaC[1,1]
				matrix __graphmatLTZ[`j',1]=delta
			}		
			matrix __graphmatLTZ[`j',2]=_b[`graph']
			matrix __graphmatLTZ[`j',3]=_b[`graph']-_se[`graph']*CI
			matrix __graphmatLTZ[`j',4]=_b[`graph']+_se[`graph']*CI

			local ++j
		}
	}
    *****************************************************************************
	*** (5cii) Form augmented var-covar and coefficient matrices (see appendix)
	*****************************************************************************
	qui estimates restore __iv
	mat V1 = e(V) +  inv(ZX'*inv(ZZ)*ZX)*ZX' * omega_in * ZX*inv(ZX'*inv(ZZ)*ZX)
	mat b1 = e(b) - (inv(ZX'*inv(ZZ)*ZX)*ZX' * mu_in)'

	*****************************************************************************
	*** (5d) Determine lower and upper bounds
	*****************************************************************************
	mat CI  = -invnormal((1-`level')/2)
	mat lb  = b1 - vecdiag(cholesky(diag(vecdiag(V1))))*CI
	mat ub  = b1 + vecdiag(cholesky(diag(vecdiag(V1))))*CI

	dis _newline	
	dis "Conley et al. (2012)'s LTZ results" _col(55) "Number of obs =    " e(N)

	ereturn post b1 V1
	ereturn display
}

********************************************************************************
*** (6) Visualise as per Conely et al., (2012) Figures 1-2
********************************************************************************
if length("`graph'")!=0 {
	if "`method'"=="uci" {
		svmat __graphmat	

		if length("`graphopts'")==0 {
			local graphopts ytitle("Estimated Beta for `graph'") /*
			*/ xtitle("{&delta}") title("Union of Confidence Interval Approach") /*
			*/ note("Methodology described in Conley et al. (2012)")
		}
	
		twoway line __graphmat3 __graphmat1, lpattern(dash) lcolor(black) || ///
		line       __graphmat4 __graphmat1, lpattern(dash) lcolor(black)     ///
		legend(order(1 "Upper Bound (UCI)" 2 "Lower Bound (UCI)"))           ///
		scheme(s1color) `graphopts'
		
	}

	if "`method'"=="ltz" {
		svmat __graphmatLTZ

		if length("`graphopts'")==0 {
			local graphopts ytitle("{&beta}") xtitle("{&delta}") /*
			*/ title("Local to Zero Approach") /*
			*/ note("Methodology described in Conley et al. (2012)")
		}
			
		twoway line __graphmatLTZ2 __graphmatLTZ1,                        || ///
		line  __graphmatLTZ3 __graphmatLTZ1, lpattern(dash) lcolor(black) || ///
		line  __graphmatLTZ4 __graphmatLTZ1, lpattern(dash) lcolor(black)    ///
		legend(order(1 "Point Estimate (LTZ)" 2 "CI (LTZ)"))                 ///
		scheme(s1color) `graphopts'
	  
	}

	cap drop __graphmat*

}


********************************************************************************
*** (7) Clean up
********************************************************************************
estimates drop __iv

end

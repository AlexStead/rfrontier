program rfrontier_p
	version 13.1
	syntax newvarname [in] [if], [XB Residuals JLMS BC BS FN]
	local nopts : word count `residuals' `xb' `jlms' `bc' `bs' `fn'
	local neffopts: word count `jlms' `bc' `fn'
	marksample touse, novarlist
	if `nopts' > 1{
		display "{err}only one statistic may be specified"
		exit 498
	}
	local s $cost_s
	local N = _N
	local tempname_list b xb_ sv su
	if "$vdistribution"=="student" local tempname_list `tempname_list' df
	tempname `tempname_list'
	matrix `b'									=	e(b)
	matrix score double `xb_'							=	`b'
	scalar `sv'									=	exp([lnsigma_v]_cons)
	scalar `su'									=	exp([lnsigma_u]_cons)
	if "$vdistribution" == "student" scalar `df'					=	exp([lndf]_cons)
	if `nopts' == 0{
		gen `typlist' `varlist'							=	`xb_'							if `touse'
	}
	if "`residuals'" != ""{
		gen `typlist' `varlist'							=	$ML_y1-`xb_'						if `touse'
	}
	if "`xb'" != ""{
		gen `typlist' `varlist' 						=	`xb_'							if `touse'
	}
	if `neffopts' > 0{
		quietly {
			tempname numerator denominator numeratorvar denominatorvar indices
			gen double `numerator'						=	0
			gen double `denominator'					=	0
			mata hdraws(`N',1000,3,10,1,"`indices'")
			qui ds
			forval i = 1/1000{
				tempvar f1_`i'
				gen double `f1_`i'' = `:word `=`indices'[1,`i']' of `r(varlist)''
			}
			forval q = 1/1000 {
				if "$udistribution"=="hnormal" local u_`q'		=	`su'*invnormal(1/2*(1+`f1_`q''))
				if "$udistribution"=="exponential" local u_`q'		=	-`su'*ln(`f1_`q'')
				if "$udistribution"=="rayleigh" local u_`q'		=	`su'*sqrt(-2*ln(1-`f1_`q''))
				if "`jlms'"!=="" local fu `u_`q''
				if "`bc'"!=="" local fu exp(-`u_`q'')
				if "$vdistribution"=="student" replace `numerator'	=	`numerator'+`fu'/`sv'*tden(`df',($ML_y1-`xb_'+`s'*`u_`q'')/`sv')	
				if "$vdistribution"=="student" replace `denominator'	=	`denominator'+1/`sv'*tden(`df',($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="cauchy" replace `numerator'	=	`numerator'+`fu'/`sv'*tden(1,($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="cauchy" replace `denominator'	=	`denominator'+1/`sv'*tden(1,($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="logistic" replace `numerator'	=	`numerator'+`fu'/`sv'*(exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))/(1+exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))^2
				if "$vdistribution"=="logistic" replace `denominator'	=	`denominator'+1/`sv'*(exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))/(1+exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))^2
			}
			if "`jlms'"!=="" gen `typlist' `varlist'			=	exp(-`numerator'/`denominator')
			if "`bc'"!=="" gen `typlist' `varlist'				=	`numerator'/`denominator'	
		}
	}
	if "`bs'" != ""{
		quietly {
			tempname numerator denominator numeratorvar denominatorvar gu indices
				gen double `numerator'					=	0
			gen double `denominator'					=	0
			gen double `numeratorvar'					=	0
			gen double `denominatorvar'					=	0
			gen double `gu'							=	0
			mata hdraws(`N',1000,3,10,1,"`indices'")
			qui ds
			forval i = 1/1000{
				tempvar f1_`i'
				gen double `f1_`i'' = `:word `=`indices'[1,`i']' of `r(varlist)''
			}	
			forval q = 1/1000 {
				if "$udistribution"=="hnormal" local u_`q'		=	`su'*invnormal(1/2*(1+`f1_`q''))
				if "$udistribution"=="exponential" local u_`q'		=	-`su'*ln(`f1_`q'')
				if "$udistribution"=="rayleigh" local u_`q'		=	`su'*sqrt(-2*ln(1-`f1_`q''))
				local fu exp(-`u_`q'')
				if "$vdistribution"=="student" replace `numerator'	=	`numerator'+`fu'/`sv'*tden(`df',($ML_y1-`xb_'+`s'*`u_`q'')/`sv')	
				if "$vdistribution"=="student" replace `denominator'	=	`denominator'+1/`sv'*tden(`df',($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="cauchy" replace `numerator'	=	`numerator'+`fu'/`sv'*tden(1,($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="cauchy" replace `denominator'	=	`denominator'+1/`sv'*tden(1,($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="logistic" replace `numerator'	=	`numerator'+`fu'/`sv'*(exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))/(1+exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))^2
				if "$vdistribution"=="logistic" replace `denominator'	=	`denominator'+1/`sv'*(exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))/(1+exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))^2
			}
			forval q = 1/1000 {
				replace `gu'						=	(exp(-`u_`q'')-`numerator'/`denominator')^2
				if "$vdistribution"=="student" replace `numeratorvar'	=	`numeratorvar'+`gu'/`sv'*tden(`df',($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="student" replace `denominatorvar'	=	`denominatorvar'+1/`sv'*tden(`df',($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="cauchy" replace `numeratorvar'	=	`numeratorvar'+`gu'/`sv'*tden(1,($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="cauchy" replace `denominatorvar'	=	`denominatorvar'+1/`sv'*tden(1,($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="logistic" replace `numeratorvar'	=	`numeratorvar'+`gu'/`sv'*(exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))/(1+exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))^2
				if "$vdistribution"=="logistic" replace `denominatorvar'=	`denominatorvar'+1/`sv'*(exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))/(1+exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))^2
			}
			gen `typlist' `varlist'						=	`numeratorvar'/`denominatorvar'
		}
	}
end

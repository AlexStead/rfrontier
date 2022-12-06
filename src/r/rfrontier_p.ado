program rfrontier_p
	version 13.1
	syntax newvarname [in] [if], [XB Residuals JLMS BC BS FN INFluence BCINFluence JLMSINFluence NSIMulations(integer 1000) BASE(real 3) BURN(real 10)]
	local nopts : word count `residuals' `xb' `jlms' `bc' `bs' `fn' `influence' `bcinfluence' `jlmsinfluence'
	local neffopts: word count `jlms' `bc' `fn' 
	local ninfeffopts: word count `jlmsinfluence' `bcinfluence'
	local ninfopts: word count `influence' `jlmsinfluence' `bcinfluence'
	marksample touse, novarlist
	global nsimulations				=	`nsimulations'
	global base						=	`base'
	global burn						=	`burn'
	if `nopts' > 1{
		display "{err}only one statistic may be specified"
		exit 498
	}
	local s $cost_s
	local N = _N
	local npar = colsof(e(b))
	if `neffopts' > 0{
		local sv e(sigma_v)
		local su e(sigma_u)
		local df e(df)
	}
	if `neffopts' == 0{
		if "$vdistribution" == "student" {
			local su e(sigma_u)
			local df exp(xb(#3))
		}
		local sv exp(xb(#2))
		if "$vdistribution" == "student" {
			local su exp(xb(#4))
			local df exp(xb(#3))
		}
		else local su exp(xb(#3))
	}
	local sigma sqrt((`sv')^2+(`su')^2)
	local lambda (`su')/(`sv')
	local gamma (`su')^2/(`sigma')^2
	local epsilon $ML_y1-xb(#1)
	local tempname_list b xb_ 
	if "$vdistribution"=="student" local tempname_list `tempname_list'
	tempname `tempname_list' gradient gradientmat
	matrix `b'									=	e(b)
	matrix score double `xb_'							=	`b'
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
			mata hdraws(`N',`nsimulations',`base',`burn',1,"`indices'")
			qui ds
			forval i = 1/`nsimulations'{
				tempvar f1_`i'
				gen double `f1_`i'' = `:word `=`indices'[1,`i']' of `r(varlist)''
			}
			forval q = 1/`nsimulations' {
				if "$udistribution"=="hnormal"		local u_`q'		=	`su'*invnormal(1/2*(1+`f1_`q''))
				if "$udistribution"=="exponential"	local u_`q'		=	-`su'*ln(`f1_`q'')
				if "$udistribution"=="rayleigh"		local u_`q'		=	`su'*sqrt(-2*ln(1-`f1_`q''))
				if "`jlms'"!="" local fu `u_`q''
				if "`bc'"!="" local fu exp(-`u_`q'')
				if "$vdistribution"=="normal"		replace `numerator'	=	`numerator'+`fu'/`sv'*normalden(($ML_y1-`xb_'+`s'*`u_`q'')/`sv')	
				if "$vdistribution"=="normal"		replace `denominator'	=	`denominator'+1/`sv'*normalden(($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="student"		replace `numerator'	=	`numerator'+`fu'/`sv'*tden(`df',($ML_y1-`xb_'+`s'*`u_`q'')/`sv')	
				if "$vdistribution"=="student"		replace `denominator'	=	`denominator'+1/`sv'*tden(`df',($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="cauchy"		replace `numerator'	=	`numerator'+`fu'/`sv'*tden(1,($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="cauchy"		replace `denominator'	=	`denominator'+1/`sv'*tden(1,($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="logistic"		replace `numerator'	=	`numerator'+`fu'/`sv'*(exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))/(1+exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))^2
				if "$vdistribution"=="logistic"		replace `denominator'	=	`denominator'+1/`sv'*(exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))/(1+exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))^2
			}
			if "`jlms'"!="" gen `typlist' `varlist'			=	exp(-`numerator'/`denominator')
			if "`bc'"!="" gen `typlist' `varlist'				=	`numerator'/`denominator'	
		}
	}
	if `ninfeffopts' > 0{
		quietly {
			tempname numerator denominator numeratorvar denominatorvar numeratorgradient denominatorgradient indices
			mata hdraws(`N',`nsimulations',`base',`burn',1,"`indices'")
			qui ds
			forval i = 1/`nsimulations'{
				tempvar f1_`i'
				gen double `f1_`i'' = `:word `=`indices'[1,`i']' of `r(varlist)''
			}
			forval q = 1/`nsimulations' {
				if "$udistribution"=="hnormal"		local u_`q'			`su'*invnormal(1/2*(1+`f1_`q''))
				if "$udistribution"=="exponential"	local u_`q'			-`su'*ln(`f1_`q'')
				if "$udistribution"=="rayleigh"		local u_`q'			`su'*sqrt(-2*ln(1-`f1_`q''))
				if "`jlms'"!="" local fu `u_`q''
				if "`bc'"!="" local fu exp(-`u_`q'')
				if "`jlmsinfluence'"!="" local fu `u_`q''
				if "`bcinfluence'"!="" local fu exp(-`u_`q'')
				if `ninfopts' > 0{
					if "$vdistribution"=="normal"		predictnl double `numerator'`q'		=	1/$nsimulations*`fu'/`sv'*normalden((`epsilon'+`s'*`u_`q'')/`sv'), g(`numeratorgradient'`q') force
					if "$vdistribution"=="normal"		predictnl double `denominator'`q'	=	1/$nsimulations*1/`sv'*normalden((`epsilon'+`s'*`u_`q'')/`sv'), g(`denominatorgradient'`q') force
					if "$vdistribution"=="student"		predictnl double `numerator'`q'		=	1/$nsimulations*`fu'/`sv'*tden(`df',(`epsilon'+`s'*`u_`q'')/`sv'), g(`numeratorgradient'`q') force
					if "$vdistribution"=="student"		predictnl double `denominator'`q'	=	1/$nsimulations*1/`sv'*tden(`df',(`epsilon'+`s'*`u_`q'')/`sv'), g(`denominatorgradient'`q') force
					if "$vdistribution"=="cauchy"		predictnl double `numerator'`q'		=	1/$nsimulations*`fu'/`sv'*tden(1,(`epsilon'+`s'*`u_`q'')/`sv'), g(`numeratorgradient'`q') force
					if "$vdistribution"=="cauchy"		predictnl double `denominator'`q'	=	1/$nsimulations*1/`sv'*tden(1,(`epsilon'+`s'*`u_`q'')/`sv'), g(`denominatorgradient'`q') force
					if "$vdistribution"=="logistic"		predictnl double `numerator'`q'		=	1/$nsimulations*`fu'/`sv'*(exp((`epsilon'+`s'*`u_`q'')/`sv'))/(1+exp((`epsilon'+`s'*`u_`q'')/`sv'))^2, g(`numeratorgradient'`q') force
					if "$vdistribution"=="logistic"		predictnl double `denominator'`q'	=	1/$nsimulations*1/`sv'*(exp((`epsilon'+`s'*`u_`q'')/`sv'))/(1+exp((`epsilon'+`s'*`u_`q'')/`sv'))^2, g(`denominatorgradient'`q') force
				}
				else{
					if "$vdistribution"=="normal"		predictnl double `numerator'`q'		=	1/$nsimulations*`fu'/`sv'*normalden((`epsilon'+`s'*`u_`q'')/`sv')
					if "$vdistribution"=="normal"		predictnl double `denominator'`q'	=	1/$nsimulations*1/`sv'*normalden((`epsilon'+`s'*`u_`q'')/`sv')
					if "$vdistribution"=="student"		predictnl double `numerator'`q'		=	1/$nsimulations*`fu'/`sv'*tden(`df',(`epsilon'+`s'*`u_`q'')/`sv')
					if "$vdistribution"=="student"		predictnl double `denominator'`q'	=	1/$nsimulations*1/`sv'*tden(`df',(`epsilon'+`s'*`u_`q'')/`sv')
					if "$vdistribution"=="cauchy"		predictnl double `numerator'`q'		=	1/$nsimulations*`fu'/`sv'*tden(1,(`epsilon'+`s'*`u_`q'')/`sv')
					if "$vdistribution"=="cauchy"		predictnl double `denominator'`q'	=	1/$nsimulations*1/`sv'*tden(1,(`epsilon'+`s'*`u_`q'')/`sv')
					if "$vdistribution"=="logistic"		predictnl double `numerator'`q'		=	1/$nsimulations*`fu'/`sv'*(exp((`epsilon'+`s'*`u_`q'')/`sv'))/(1+exp((`epsilon'+`s'*`u_`q'')/`sv'))^2
					if "$vdistribution"=="logistic"		predictnl double `denominator'`q'	=	1/$nsimulations*1/`sv'*(exp((`epsilon'+`s'*`u_`q'')/`sv'))/(1+exp((`epsilon'+`s'*`u_`q'')/`sv'))^2
				}
			}
			egen double `numerator' = rowtotal(`numerator'*)
			egen double `denominator' = rowtotal(`denominator'*)			
			if `ninfopts' == 0 {
				if "`jlms'"!="" gen `typlist' `varlist'			=	exp(-`numerator'/`denominator')
				if "`bc'"!="" gen `typlist' `varlist'			=	`numerator'/`denominator'
			}
			if `ninfopts' > 0{
				forval i = 1/`npar'{
					egen double `numeratorgradient'`i' = rowtotal(`numeratorgradient'*`i')
					egen double `denominatorgradient'`i' = rowtotal(`denominatorgradient'*`i')
					if "`jlms'"!="" gen double `gradient'`i' = -exp(-`numerator'/`denominator')*(`numeratorgradient'`i'*`denominator'-`denominatorgradient'`i'*`numerator')/(`denominator')^2
					if "`bc'"!="" gen double `gradient'`i' = (`numeratorgradient'`i'*`denominator'-`denominatorgradient'`i'*`numerator')/(`denominator')^2
					if "`jlmsinfluence'"!="" gen double `gradient'`i' = -exp(-`numerator'/`denominator')*(`numeratorgradient'`i'*`denominator'-`denominatorgradient'`i'*`numerator')/(`denominator')^2
					if "`bcinfluence'"!="" gen double `gradient'`i' = (`numeratorgradient'`i'*`denominator'-`denominatorgradient'`i'*`numerator')/(`denominator')^2
			}
			mkmat `gradient'*, matrix(`gradientmat')
			}	
		}
	}
	if "`bs'" != ""{
		quietly {
			tempname numerator denominator numeratorvar denominatorvar indices
				gen double `numerator'					=	0
			gen double `denominator'					=	0
			gen double `numeratorvar'					=	0
			gen double `denominatorvar'					=	0
			mata hdraws(`N',`nsimulations',`base',`burn',1,"`indices'")
			qui ds
			forval i = 1/`nsimulations'{
				tempvar f1_`i'
				gen double `f1_`i'' = `:word `=`indices'[1,`i']' of `r(varlist)''
			}	
			forval q = 1/`nsimulations' {
				if "$udistribution"=="hnormal"		local u_`q'			`su'*invnormal(1/2*(1+`f1_`q''))
				if "$udistribution"=="exponential"	local u_`q'			-`su'*ln(`f1_`q'')
				if "$udistribution"=="rayleigh"		local u_`q'			`su'*sqrt(-2*ln(1-`f1_`q''))
				local fu exp(-`u_`q'')
				if "$vdistribution"=="normal"		predictnl double `numerator'`q'		=	1/$nsimulations*`fu'/`sv'*normalden(($ML_y1-`xb_'+`s'*`u_`q'')/`sv')	
				if "$vdistribution"=="normal"		predictnl double `denominator'`q'	=	1/$nsimulations*1/`sv'*normalden(($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="student"		predictnl double `numerator'`q'		=	1/$nsimulations*`fu'/`sv'*tden(`df',($ML_y1-`xb_'+`s'*`u_`q'')/`sv')	
				if "$vdistribution"=="student"		predictnl double `denominator'`q'	=	1/$nsimulations*1/`sv'*tden(`df',($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="cauchy"		predictnl double `numerator'`q'		=	1/$nsimulations*`fu'/`sv'*tden(1,($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="cauchy"		predictnl double `denominator'`q'	=	1/$nsimulations*1/`sv'*tden(1,($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="logistic"		predictnl double `numerator'`q'		=	1/$nsimulations*`fu'/`sv'*(exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))/(1+exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))^2
				if "$vdistribution"=="logistic"		predictnl double `denominator'`q'	=	1/$nsimulations*1/`sv'*(exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))/(1+exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))^2
			}
			egen double `numerator' = rowtotal(`numerator'*)
			egen double `denominator' = rowtotal(`denominator'*)
			forval q = 1/`nsimulations' {
									local `gu'						(exp(-`u_`q'')-`numerator'/`denominator')^2
				if "$vdistribution"=="normal"		predictnl double `numeratorvar'`q'		=	`gu'/`sv'*normalden(`df',($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="normal"		predictnl double `denominatorvar'`q'		=	1/`sv'*normalden(`df',($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="student"		predictnl double `numeratorvar'`q'		=	`gu'/`sv'*tden(`df',($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="student"		predictnl double `denominatorvar'`q'		=	1/`sv'*tden(`df',($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="cauchy"		predictnl double `numeratorvar'`q'		=	`gu'/`sv'*tden(1,($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="cauchy"		predictnl double `denominatorvar'`q'		=	1/`sv'*tden(1,($ML_y1-`xb_'+`s'*`u_`q'')/`sv')
				if "$vdistribution"=="logistic"		predictnl double `numeratorvar'`q'		=	`gu'/`sv'*(exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))/(1+exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))^2
				if "$vdistribution"=="logistic"		predictnl double `denominatorvar'`q'		=	1/`sv'*(exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))/(1+exp(($ML_y1-`xb_'+`s'*`u_`q'')/`sv'))^2
			}
			egen double `numeratorvar' = rowtotal(`numeratorvar'*)
			egen double `denominatorvar' = rowtotal(`denominatorvar'*)
			gen `typlist' `varlist'						=	`numeratorvar'/`denominatorvar'
		}
	}
	if `ninfopts' > 0 {
		quietly{
			local sv exp(xb(#2))
			if "$vdistribution" == "student" {
				local su exp(xb(#4))
				local df exp(xb(#3))
			}
			else local su exp(xb(#3))
			local sigma sqrt((`sv')^2+(`su')^2)
			local lambda (`su')/(`sv')
			local gamma (`su')^2/(`sigma')^2
			local epsilon $ML_y1-xb(#1)
			tempname indices scoredenom scoredenomsum scorenum scorenumsum scorepar score invinformation influencemat influencepar influenceeff influenceeffob
			mata hdraws(`=_N',`nsimulations',`base',`burn',1,"`indices'")
			qui ds
			forval q=1/`nsimulations'{
				tempvar f1_`q'
				gen double `f1_`q'' = `:word `=`indices'[1,`q']' of `r(varlist)''
			}
			forval q=1/`nsimulations'{
				if "$udistribution"=="hnormal" local u_`q'	`su'*invnormal(1/2*(1+`f1_`q''))
				if "$udistribution"=="exponential" local u_`q' -`su'*ln(`f1_`q'')
				if "$udistribution"=="rayleigh" local u_`q'	`su'*sqrt(-2*ln(1-`f1_`q''))
				if "$vdistribution"=="normal" predictnl double `scoredenom'`q' = 1/$nsimulations*1/`sv'*normalden((`epsilon'+`s'*`u_`q'')/`sv'), g(`scorenum'`q') force
				if "$vdistribution"=="student" predictnl double `scoredenom'`q' = 1/$nsimulations*1/`sv'*tden(`df',(`epsilon'+`s'*`u_`q'')/`sv'), g(`scorenum'`q') force
				if "$vdistribution"=="cauchy" predictnl double `scoredenom'`q' = 1/$nsimulations*1/`sv'*tden(1,(`epsilon'+`s'*`u_`q'')/`sv'), g(`scorenum'`q') force
				if "$vdistribution"=="logistic" predictnl double `scoredenom'`q' = 1/$nsimulations*1/`sv'*logisticden((`epsilon'+`s'*`u_`q'')/`sv'), g(`scorenum'`q') force
			}
			egen double `scoredenomsum' = rowtotal(`scoredenom'*)
			forval i = 1/`npar'{
				egen double `scorenumsum'`i' = rowtotal(`scorenum'*`i')
				gen double `scorepar'`i' = -`scorenumsum'`i'/`scoredenomsum'
			}
			mkmat `scorepar'*, matrix(`score')
			matrix `invinformation' = e(V)
			matrix `influencemat' = -(`invinformation'*`score'')'
			if "`influence'" != "" {
				svmat double `influencemat', name(`influencepar')
				if "$vdistribution" == "student" {
					replace `influencepar'`=`npar'-2' = `influencepar'`=`npar'-2'*exp(_b[lnsigma_v:_cons])
					replace `influencepar'`=`npar'' = `influencepar'`=`npar'='*exp(_b[lnsigma_u:_cons])
					replace `influencepar'`=`npar'-1' = `influencepar'`=`npar'-1'*exp(_b[lndf:_cons])
				}
				else {
					replace `influencepar'`=`npar'-1' = `influencepar'`=`npar'-1'*exp(_b[lnsigma_v:_cons])
					replace `influencepar'`=`npar'' = `influencepar'`=`npar'='*exp(_b[lnsigma_u:_cons])				
				}
				forval i = 1/`npar'{
					gen `typlist' `varlist'`i'						=	`influencepar'`i'
				}
			}
			else {
				matrix `influenceeff' = `gradientmat'*`influencemat''
				svmat double `influenceeff', name(`influenceeffob')
				forval i = 1/`N'{
					gen `typlist' `varlist'`i'						=	`influenceeffob'`i'
				}
			}
		}	
	}
end

mata matrix function hdraws(N,D,P,B,A,I) {
	if (A==1) D 		= D/2 ;
	hdrawsvec 			= J(N*D,1,.)
	hdrawsvecburn		= ghalton(N*D+B,P,0)
	for (i=1; i<=N*D; i++) {
		hdrawsvec[i,1]	= hdrawsvecburn[i+B,1]
	}
	_hdraws				= colshape(hdrawsvec,D)
	if (A==1) hdraws	= _hdraws,(J(N,D,1) - _hdraws)
	else hdraws			= _hdraws
	names				= st_tempname(cols(hdraws))
	for (i=1; i<=cols(hdraws); i++) {
		st_addvar("double", names[i])
		st_store(.,names[i], hdraws[,i])
	}
	st_matrix(I,st_varindex(names))
}

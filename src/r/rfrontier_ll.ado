program rfrontier_ll
		version 13.1
		args todo b lnf g
		local tempvar_list mu lnsigma_v lnsigma_u fj
		if "$vdistribution"=="student" local tempvar_list `tempvar_list' lndf
		tempname `tempvar_list' indices
		mleval `mu'										= 	`b',		eq(1)
		mleval `lnsigma_v' 									= 	`b',		eq(2)		scalar
		mleval `lnsigma_u'									=	`b',		eq(3)		scalar
		if "$vdistribution"=="student" mleval `lndf'						= 	`b',		eq(4)		scalar
		quietly {
		if $ld == 0{
				mata hdraws(`=_N',$nsimulations,$base,10,1,"`indices'")
				}
		if $ld == 1{
				mata hdraws(`=_N',$nsimulations,$base,10,"`indices'")
				}
		}
		qui ds
		forval i = 1/$nsimulations{
				tempvar f1_`i'
				gen double `f1_`i'' = `:word `=`indices'[1,`i']' of `r(varlist)''
				}
		quietly {
				local s									=	$cost_s
				gen double `fj'								=	0
				forval q								=	1/$nsimulations{
						if "$udistribution"=="hnormal" local u_`q'		=	exp(`lnsigma_u')*invnormal(1/2*(1+`f1_`q''))
						if "$udistribution"=="exponential" local u_`q'		=	-exp(`lnsigma_u')*ln(`f1_`q'')
						if "$udistribution"=="rayleigh" local u_`q'		=	exp(`lnsigma_u')*sqrt(-2*ln(1-`f1_`q''))
						if "$vdistribution"=="student" replace `fj'		=	`fj'+1/exp(`lnsigma_v')*tden(exp(`lndf'),($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))
						if "$vdistribution"=="logistic" replace `fj'		=	`fj'+1/exp(`lnsigma_v')*exp(-($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))/(1+exp(-($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v')))^2
						if "$vdistribution"=="cauchy" replace `fj'		=	`fj'+1/exp(`lnsigma_v')*tden(1,($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))
				}
				mlsum `lnf' 								=	ln(`fj'/$nsimulations)
		}
		if(`todo'==0 | `lnf'>=.) exit
		tempname d1 d2 d3 d4 sum_g sum_dgde sum_dgdsv sum_dgdsu sum_dgda
		quietly	{
				gen double `sum_g'							=	0
				gen double `sum_dgde'							=	0
				gen double `sum_dgdsv'							=	0
				gen double `sum_dgdsu'							=	0
				gen double `sum_dgda'							=	0
				forval q								=	1/$nsimulations{
						if "$udistribution"=="hnormal" local u_`q'		=	exp(`lnsigma_u')*invnormal(1/2*(1+`f1_`q''))
						if "$udistribution"=="exponential" local u_`q'		=	-exp(`lnsigma_u')*ln(`f1_`q'')
						if "$udistribution"=="rayleigh" local u_`q'		=	exp(`lnsigma_u')*sqrt(-2*ln(1-`f1_`q''))
						if "$udistribution"=="hnormal" local dudsu		=	invnormal(1/2*(1+`f1_`q''))
						if "$udistribution"=="exponential" local dudsu		=	-ln(`f1_`q'')
						if "$udistribution"=="rayleigh" local dudsu		=	sqrt(-2*ln(1-`f1_`q''))
						if "$vdistribution"=="student" replace `sum_g'		=	`sum_g'+(1+1/exp(`lndf')*(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))^2)^(-(exp(`lndf')+1)/2)
						if "$vdistribution"=="student" replace `sum_dgde'	=	`sum_dgde'-(exp(`lndf')+1)/exp(`lndf')*1/exp(`lnsigma_v')*($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v')*(1+1/exp(`lndf')*(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))^2)^(-(exp(`lndf')+3)/2)
						if "$vdistribution"=="student" replace `sum_dgdsv'	=	`sum_dgdsv'+(exp(`lndf')+1)/exp(`lndf')*1/exp(`lnsigma_v')*(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))^2*(1+1/exp(`lndf')*(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))^2)^(-(exp(`lndf')+3)/2)
						if "$vdistribution"=="student" replace `sum_dgdsu'	=	`sum_dgdsu'-`s'*(exp(`lndf')+1)/exp(`lndf')*1/exp(`lnsigma_v')*($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v')*(1+1/exp(`lndf')*(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))^2)^(-(exp(`lndf')+3)/2)*`dudsu'
						if "$vdistribution"=="student" replace `sum_dgda'	=	`sum_dgda'-1/2*(1+1/exp(`lndf')*(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))^2)^(-(exp(`lndf')+1)/2)*(ln(1+1/exp(`lndf')*(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))^2)-(exp(`lndf')+1)*(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))^2/(exp(`lndf')^2*(1+1/exp(`lndf')*(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))^2)))
						if "$vdistribution"=="logistic" replace `sum_g'		=	`sum_g'+exp(-($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))/(1+exp(-($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v')))^2
						if "$vdistribution"=="logistic" replace `sum_dgde'	=	`sum_dgde'+1/exp(`lnsigma_v')*(1-exp(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v')))*exp(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))/(1+exp(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v')))^3
						if "$vdistribution"=="logistic" replace `sum_dgdsv'	=	`sum_dgdsv'-($ML_y1-`mu'+`s'*`u_`q'')/(exp(`lnsigma_v'))^2*(1-exp(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v')))*exp(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))/(1+exp(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v')))^3
						if "$vdistribution"=="logistic" replace `sum_dgdsu'	=	`sum_dgdsu'+`s'*1/exp(`lnsigma_v')*(1-exp(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v')))*exp(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))/(1+exp(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v')))^3*`dudsu'						
						if "$vdistribution"=="cauchy" replace `sum_g'		=	`sum_g'+(1+(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))^2)^(-1)
						if "$vdistribution"=="cauchy" replace `sum_dgde'	=	`sum_dgde'-2/exp(`lnsigma_v')*($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v')*(1+(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))^2)^(-2)
						if "$vdistribution"=="cauchy" replace `sum_dgdsv'	=	`sum_dgdsv'+2/exp(`lnsigma_v')*(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))^2*(1+(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))^2)^(-2)
						if "$vdistribution"=="cauchy" replace `sum_dgdsu'	=	`sum_dgdsu'-`s'*2/exp(`lnsigma_v')*($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v')*(1+(($ML_y1-`mu'+`s'*`u_`q'')/exp(`lnsigma_v'))^2)^(-2)*`dudsu'
						}
				mlvecsum `lnf' `d1'							=	-`sum_dgde'/`sum_g',											eq(1)
				mlvecsum `lnf' `d2'							=	exp(`lnsigma_v')*`sum_dgdsv'/`sum_g'-1,									eq(2)
				mlvecsum `lnf' `d3'							=	exp(`lnsigma_u')*`sum_dgdsu'/`sum_g',									eq(3)
				if "$vdistribution"=="student" mlvecsum `lnf' `d4'			=	exp(`lndf')*(1/2*(digamma((exp(`lndf')+1)/2)-digamma(exp(`lndf')/2)-1/exp(`lndf'))+`sum_dgda'/`sum_g'),	eq(4)
				if "$vdistribution"=="student" matrix `g'				=	(`d1',`d2',`d3',`d4')
				if inlist("$vdistribution","logistic","cauchy") matrix `g'		=	(`d1',`d2',`d3')
		}
end

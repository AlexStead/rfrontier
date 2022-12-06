//Robust stochastic frontier models
//This programme enables the estimation of stochastic frontier models in which the noise
//term follows alternative, thick-tailed distributions.

program rfrontier, eclass
	version 13.1
	if !replay() {
		syntax varlist(numeric fv ts) [if] [in] [, Level(cilevel) noCONstant COST PRODuction Vdistribution(string) Udistribution(string) NSIMulations(integer 250) BASE(real 3) BURN(real 10) DF(real 0) LImdep SVFRONTier(numlist min=1) SVVsigma(numlist min=1 max=1) SVUsigma(numlist min=1 max=1) SVDF(numlist min=1 max=1) *]
		tempname _sigma_v _sigma_u _lndf
		local mlmodel_list /lnsigma_v
		local diparm_list diparm(lnsigma_v, exp label(sigma_v) prob)
		local vpars = 1
		local upars = 1
		if inlist(`"`vdistribution'"',"student","studen","stude","stud","stu","st","s","") {
			local vdistribution student
			local title_vdist Student's t
			local mlmodel_list `mlmodel_list' /lndf
			local diparm_list `diparm_list' diparm(lndf, exp label(df) prob)
			local vpars = 2
		}
		local diparm_list `diparm_list' diparm(lnsigma_u, exp label(sigma_u) prob)
		local mlmodel_list  `mlmodel_list' /lnsigma_u
		if inlist(`"`vdistribution'"',"normal","norma","norm","nor","no","n") {
		local vdistribution normal
		local title_vdist normal
		}
		if inlist(`"`vdistribution'"',"laplace","laplac","lapla","lapl","la") {
		local vdistribution laplace
		local title_vdist Laplace
		}
		if inlist(`"`vdistribution'"',"logistic","logisti","logist","logis","logi","log","lo") {
		local vdistribution logistic
		local title_vdist logistic
		}
		if inlist(`"`vdistribution'"',"cauchy","cauch","cauc","cau","ca","c") {
			local vdistribution cauchy
			local title_vdist Cauchy
		}
		if !inlist(`"`vdistribution'"',"student","normal","laplace","logistic","cauchy") {
			display as err "vdistribution(`vdistribution') not allowed"
			exit 198
		}
		if inlist(`"`udistribution'"',"hnormal","hnorma","hnorm","hnor","hno","hn","h","") {
			local udistribution hnormal
			local title_udist half normal
		}
		if inlist(`"`udistribution'"',"exponential","exponentia","exponenti","exponent","exponen","expone") {
			local udistribution exponential
			local title_udist exponential
		}
		if inlist(`"`udistribution'"',"expon","expo","exp","ex","e") {
			local udistribution exponential
			local title_udist exponential
		}
		if inlist(`"`udistribution'"',"rayleigh","rayleig","raylei","rayle","rayl","ray","ra","r") {
			local udistribution rayleigh
			local title_udist Rayleigh
		}
		if inlist(`"`udistribution'"',"gamma","gamm","gam","ga","g") {
			local udistribution gamma
			local title_udist gamma
			local mlmodel_list  `mlmodel_list' /lnk
			local diparm_list diparm(lnk, exp label(k) prob)
			local upars = 2
		}
		if !inlist(`"`udistribution'"',"","hnormal","exponential","rayleigh","gamma") {
			display as err "udistribution(`udistribution') not allowed"
			exit 198
		}
		global vdistribution `vdistribution'
		global udistribution `udistribution'
		local s							=	-1
		if "`cost'"=="" local s					=	1
		local ld						=	1
		if "`limdep'"=="" local ld				=	0
		global cost_s `s'
		global nsimulations `nsimulations'
		local nsimulations `:di  %9.0fc $nsimulations'
		local nsimshalf = $nsimulations/2
		mlopts mlopts, `options'
		marksample touse
		gettoken yvar xvars : varlist
		quietly {
			if "`udistribution'" == "rayleigh" capture frontier `yvar' `xvars' if `touse', `cost'
			else if "`udistribution'" == "gamma" capture frontier `yvar' `xvars' if `touse', `cost'
			else capture frontier `yvar' `xvars' if `touse', `cost' d(`udistribution')
		}
		scalar _b_cons=_b[_cons]
		scalar `_sigma_v' 					=	ln(e(sigma_v))
		if "`vdistribution'"=="logistic" scalar `_sigma_v'	=	ln(e(sigma_v))+1/2*ln(3)-ln(_pi)
		scalar `_sigma_u'					=	ln(e(sigma_u))
		tempname b0 svfront
		mat `b0' = e(b)
		local colnames			`: colnames `b0''
		if "`vdistribution'"=="student" local colnames `colnames' _cons
		if "`udistribution'"=="gamma" local colnames `colnames' _cons
		mat `b0'[1,`=colsof(`b0')'] 				=	`b0'[1,`=colsof(`b0')']/2
		mat `b0'[1,`=colsof(`b0')-1'] 				=	`b0'[1,`=colsof(`b0')-1']/2
		if `vpars'==2 mat `b0' = (`b0',0)
		if `upars'==2 mat `b0' = (`b0',0)
		if "`svusigma'"!="" mat `b0'[1,`=colsof(`b0')-`upars'+1']		=	`svusigma'
		if "`svvsigma'"!="" mat `b0'[1,`=colsof(`b0')-`upars'-`vpars'+1'] 	=	`svvsigma'
		if "`svfrontier'"!="" {
			numlist "`svfrontier'"
			local svfrontier `r(numlist)'
			foreach x in `svfrontier' {
				matrix `svfront' = nullmat(`svfront'),`x'
			}
		}
		forval i = 1/`=colsof(`b0')-`vpars'-`upars''{
			if "`svfrontier'"!="" mat `b0'[1,`i'] = `svfront'[1,`i']
			local coleq `coleq' mu
		}
		local coleq `coleq' lnsigma_v
		if "`vdistribution'"=="student" local coleq `coleq' lndf
		local coleq `coleq' lnsigma_u
		if "`udistribution'"=="gamma" local coleq `coleq' lnk
		capture if "`vdistribution'"=="student" tregress `yvar' `xvars' if `touse'
		if "`vdistribution'"=="student" scalar `_lndf'		=	_b[lndf:_cons]
		else scalar `_lndf'					=	0
		if `df' > 0 {
			local lndf = ln(`df')
			constraint 1  [lndf]_cons = `lndf'
			local constraints constraints(1)
		}
		ml model d1 rfrontier_ll (mu:`yvar'=`xvars', `constant') `mlmodel_list' if `touse', crittype("log simulated likelihood") `mlopts' `constraints'  ///
		title(Stoch. frontier `title_vdist'/`title_udist' model) `diparm_list'
		local N							=	_N
		global base						=	`base'
		global burn						=	`burn'
		local base 							`:di  %9.0fc $base'
		local burn 							`:di  %9.0fc $burn'
		global ld						=	`ld'
		tempname b V ll N df_m converged sigma_v sigma_u df
		//if "`vdistribution'"=="student" local colnames			`: colnames `b0''
		if "`vdistribution'"=="student" {
			mat `b0'[1,`=colsof(`b0')']				=	`b0'[1,`=colsof(`b0')-`upars'']
			mat `b0'[1,`=colsof(`b0')-`upars''] 	=	`_lndf'
		}
		//if "`udistribution'"=="gamma" mat `b0'		=	(`b0',1)
		matrix coleq `b0' = `coleq'
		matrix colnames `b0' = `colnames'
		ml init `b0'
		ml maximize, nooutput
		}
	matrix `b'							=	e(b)
	matrix `V'							=	e(V)
	scalar `ll'							=	e(ll)
	scalar `N'							=	e(N)
	scalar `converged'						=	e(converged)
	scalar `sigma_v'						=	exp([lnsigma_v]_cons)
	scalar `sigma_u'						=	exp([lnsigma_u]_cons)
	if "`vdistribution'"=="student"  scalar `df' 			=	exp([lndf]_cons)
	//ereturn post `b' `V', esample(`touse')
	ereturn scalar ll						=	`ll'
	ereturn scalar N						=	`N'
	ereturn scalar converged					=	`converged'
	ereturn scalar sigma_v						=	`sigma_v'
	ereturn scalar sigma_u						=	`sigma_u'
	if "`vdistribution'"=="student" ereturn scalar df		=	`df'
	ereturn local predict 	"rfrontier_p"
	ereturn local cmd		"rfrontier"
	global ML_y1 `yvar'
	Replay, level(`level')
	di as text "Integral simulated using " as result "`nsimulations'" ///
	as text " Halton draws (including antithetic draws)"
	di as text "Base " ///
	as result "`base'" as text ", first " as result "`burn'" as text " draws discarded"
end

program Replay
	syntax [, Level(cilevel)]
	ml display, level(`level')
end

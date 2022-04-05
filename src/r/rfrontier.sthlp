{smcl}
{vieweralsosee "[R] help" "help help"}{...}
{viewerjumpto "Syntax" "rfrontier##syntax"}{...}
{viewerjumpto "Description" "rfrontier##description"}{...}
{viewerjumpto "Options" "rfrontier##options"}{...}
{viewerjumpto "Remarks" "rfrontier##remarks"}{...}
{viewerjumpto "Examples" "rfrontier##examples"}{...}
{title:Title}

{phang}
{bf:rfrontier} {hline 2} Robust stochastic frontier models with non-Gaussian noise distributions


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:rfrontier}
{depvar}
[{indepvars}]
{ifin}
[{cmd:,} {it:options}]

{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opt nocons:tant}}suppress constant term{p_end}
{synopt :{cmdab:v:distribution(}{opt s:tudent)}}	Student's t distribution for the noise term, the default{p_end}
{synopt :{cmdab:v:distribution(}{opt c:auchy)}}		Cauchy distribution for the noise term{p_end}
{synopt :{cmdab:v:distribution(}{opt l:ogistic)}}	Logistic distribution for the noise term{p_end}
{synopt :{cmdab:u:distribution(}{opt h:normal)}}	Half normal distribution for the inefficiency term, the default{p_end}
{synopt :{cmdab:u:distribution(}{opt e:xponential)}}	Exponential distribution for the inefficiency term{p_end}
{synopt :{cmdab:u:distribution(}{opt r:ayleigh)}}	Rayleigh distribution for the inefficiency term{p_end}

{syntab:Starting values}
{synopt :{opt svfront:ier()}}specify a 1 X k vector of initial values for the coefficients of the frontier{p_end}
{synopt :{opt svv:sigma()}}specify an initial value for the natural logarithm of the scaling parameter of the noise distribution{p_end}
{synopt :{opt svu:sigma()}}specify an initial value for the natural logarithm of the scaling parameter of the inefficiency distribution{p_end}
{synopt :{opt svdf()}}specify an initial value for the natural logarithm of the student's t distribution; only with v(student){p_end}


{syntab:Other options}
{synopt:{opt cost}}fit cost frontier model; default is production frontier model{p_end}
{synopt :{opt nsim:ulations(#)}}# of Halton draws; default is 250{p_end}
{synopt :{opt base(#)}}base used for Halton draws; default is 3{p_end}

{syntab:Reporting}
{synopt :{opt l:evel(#)}}set confidence level; default is
{cmd:level(95)}{p_end}


{title:Description}

{pstd}
{cmd:rfrontier} Estimates a stochastic production or cost frontier model in which the noise
term follows either a Student's t distribution with scaling parameter, a Cauchy distribution,
or a logistic distribution, and the inefficiency term follows either a half normal
distribution, an exponential distribution, or a Rayleigh distribution. All models are
estimated by maximum simulated likelihood.

{title:Options}

{dlgtab:Frontier}

{phang}
{opt noconstant}; see
{helpb estimation options##noconstant:[R] estimation options}.


{dlgtab:Other options}

{phang}
{opt cost} specifies that {opt cnsf} fits a cost frontier model.

{dlgtab:Reporting}

{phang}
{opt level(#)}; see
{helpb estimation options##level():[R] estimation options}.


{title:Saved results}

{pstd}
{cmd:rfrontier} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}


{title:Authors}

{pstd}Alexander Stead{p_end}
{pstd}Institute for Transport Studies, University of Leeds{p_end}
{pstd}Leeds, UK{p_end}
{pstd}a.d.stead@leeds.ac.uk{p_end}

{title:Also see}

{psee}
{space 2}Help:  {help rfrontier_postestimation}
{p_end}
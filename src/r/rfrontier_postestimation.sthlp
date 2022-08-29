{smcl}
{cmd:help rfrontier postestimation} {...}
{right:also see:  {help rfrontier}}
{hline}

{title:Title}

{p2colset 5 32 38 2}{...}
{p2col :{hi:rfrontier postestimation} {hline 2}}Postestimation tools for
rfrontier{p_end}
{p2colreset}{...}


{title:Description}

{pstd}
The following postestimation commands are available for {opt rfrontier}:

{synoptset 11}{...}
{p2coldent :command}description{p_end}
{synoptline}
{synopt :{helpb rfrontier postestimation##predict:predict}}predictions, residuals,
efficiency scores{p_end}
{synoptline}
{p2colreset}{...}


{marker predict}{...}
{title:Syntax for predict}

{p 8 16 2}{cmd:predict} {dtype} {newvar} {ifin} [{cmd:,} {it:statistic}]

{synoptset 15 tabbed}{...}
{synopthdr :statistic}
{synoptline}
{syntab :Main}
{synopt :{opt xb}}linear prediction; the default{p_end}
{synopt :{opt residuals}}residuals{p_end}
{synopt :{opt jlms}}estimates of (technical or cost) inefficiency via exp[-{it:E}(u|e)]
(Jondrow et al., 1982){p_end}
{synopt :{opt bc}}estimates of (technical or cost) inefficiency {it:E}[exp(-u|e)]
(Battese and Coelli, 1988){p_end}
{synopt :{opt inf:luence}}estimated influence functions for the parameter vector{p_end}
{synopt :{opt jlmsinf:luence}}estimated influence functions for Jondrow et al. (1982)
efficiency predictions{p_end}
{synopt :{opt bcinf:luence}}estimated influence functions for Battese and Coelli (1988)
efficiency predictions{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
These statistics are only available for the estimation sample.

INCLUDE help menu_predict


{title:Options for predict}

{dlgtab:Main}

{phang}
{opt xb}, the default, calculates the linear prediction.

{phang}
{opt residuals}, calculates the residuals.

{phang}
{opt jlms} produces estimates of (technical or cost) efficiency via exp[-E(u|e)]. 

{phang}
{opt bc} produces estimates of (technical or cost) efficiency via E[exp(-u|e)].

{phang}
{opt inf:luence} produces K variables, with each corresponding to a model parameter,
giving the element of the estimated influence function corresponding to that parameter evaluated for each observation.
Note that we calculate influence on ancillary parameters rather than the transformations that map them to the real line.

{phang}
{opt jlmsinf:luence} produces N variables, with each corresponding to an observation in the data,
giving the estimated influences of that observation on each (technical or cost) efficiency score via exp[-{it:E}(u|e)].

{phang}
{opt bcinf:luence} produces N variables, with each corresponding to an observation in the data,
giving the estimated influences of that observation on each (technical or cost) efficiency score via {it:E}[exp(-u|e)].


{title:Authors}

{pstd}Alexander Stead{p_end}
{pstd}Institute for Transport Studies, University of Leeds{p_end}
{pstd}Leeds, UK{p_end}
{pstd}a.d.stead@leeds.ac.uk{p_end}


{title:Also see}

{psee}
{space 2}Help:  {help rfrontier}
{p_end}

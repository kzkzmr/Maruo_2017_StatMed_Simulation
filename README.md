# Maruo_2017_StatMed_Simulation
Simulation program of Maruo et al. (2017). Stat. Med., 36, 15, 2420-2434

Maruo_2017_sim_parmset_PN.R: Parameter calculation program of simulations settings for the power normal distribution. "Simulation_bcmixed_PN_parm.csv" is created.

Maruo_2017_sim_parmset_WG.R: Parameter calculation program of simulations settings for the Weibull and gamma distributions. "Simulation_bcmixed_WG_parm.csv" is created.

Maruo_2017_sim_PN.R: Simulation program of simulations settings for the power normal distribution. "Simulation_bcmixed_PN_parm.csv" is used.

Maruo_2017_sim_WG.R: Simulation program of simulations settings for the Weibull and gamma distributions. "Simulation_bcmixed_PN_parm.csv" is used.


Simulation setting variables

cov = 0: no covariate, 1: with covariate (baseline); cov = 1 for all simulations.
sgn = -1: sign of true treatment effect is minus, sgn = 1: plus; sgn = 1 for all simulations. 
lambda (only for PN) = -0.5, 0, or 0.5
wg (only for WG) = 0: Weibull distribution, 1: gamma distribution
H = 0: null hypothesis for all visit, 1: null hypothes only for last visit, 2: alternative hypothesis

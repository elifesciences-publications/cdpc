# cdpc
Bayesian method to estimate embryonic age for mouse embryos, as in Hayden et al, submitted.

## Construction of the model

We built a model to improve embryonic age estimation based on observed values of embryo weight, embryo litter, and age in dpc (for the litter). 

First, we consider a basic model (a deterministic one) where the embryonic age of embryos would be equal to their age in dpc. We assume that there is a log-linear relationship between weight and dpc specific to each strain (FVB or DUHi). The slope and offset of this relationship are estimated from the data for each strain separately and termed respectively “b” and “a” in the model. The body weight in litter “i” can then be modeled as in [1]

weight.in.logi = a + b * age.in.dpci 					[1]


In practice, age.in.dpc is not a precise measure of embryonic age. Different litters at the same age in dpc may be more or less developmentally advanced (inter-litter developmental effect) [2]. This induces a variability on weight.in.log in model [1] which is described in the model [2] by the addition of a stochastic term, eps.litter.dev, following a Gaussian distribution centered on 0 with a standard deviation, sd.litter.dev, that will be estimated from data.

weight.in.logi = a + b * age.in.dpci  + eps.litter.devi				 [2]

However, different litters at the same age in dpc, and at the same embryonic age, may have different mean weights because each pregnancy will provide a specific environment to the embryos.  The addition of another stochastic term in the model, eps.litter.preg, is needed to account for this pregnancy effect [3]. We assume that this effect follows a Gaussian distribution centered on 0 with a known standard deviation (sd.litter.preg). From our practice we consider two values for sd.litter.preg : 0.05 and 0.1. The first one (0.05) corresponds to  a 95% fluctuation interval of eps.litter.preg of [-0.1 ; 0.1] which gives in the ratio of weight [-10% ; 10%], so to a realistic maximum effect on weight of 20 mg for a 200 mg embryo.  The second value (0.1) corresponds to a 95% fluctuation interval of eps.litter.preg of [-0.2 ; 0.2] which corresponds in the ratio in weight of [-18% ; 22%] so to an excessive maximum effect on weight of 40 mg for a 200 mg embryo.

weight.in.logi = a + b * age.in.dpci  + eps.litter.devi +eps.litter.pregi 		[3]

The previous complexifications of the model are drawn at the level of the litter, presuming that all embryos of a given litter were at the same embryonic age. However, we know that 1) all embryos are not at the same embryonic age and 2) within a litter, weight is a very good indicator of relative embryonic age. To take advantage of this, we modeled that within a litter, embryonic age follows a Gaussian distribution, centered on the mean stage of the litter. So to describe the weight of the embryo j in the litter i, we add to the model [4] a last term, eps.embryo.dev, following a Gaussian distribution centered on 0 with a standard deviation of sd.embryo.dev estimated from the data and characterizing intra-litter variability on body weight.

Weight.in.logij = a + b * age.in.dpci + eps.litter.devi + eps.litter.pregi  + eps.embryo.devij 	[4]


## Estimation of embryonic age

Once the model fitted, the development age was estimated from model 4 just by removing the pregnancy effect eps.litter.pregi  and inverting the relation :
age.devij = (weight.in.logij - eps.litter.pregi – a) / b		 [5]
Note that this estimation of the embryonic age (age.devij) from the body weight in log of each embryo (weight.in.logij) requires the knowledge of parameters a and b and of random effects due to pregnancy for each litter (eps.litter.pregi). Those were previously estimated from data as explained below.


## Estimation of the parameters

Parameters a, b, sd.litter.dev, sd.embryo.dev and random effects of the model were estimated from data in two  scenarios  for two fixed sd.litter.preg values (0.05 in the realistic scenario and 0.10 in the excessive scenario) as it is not possible to dissociate only from data the two components of inter-litter variability : variability in weight due to embryonic age and to pregnancy. Vague uniform priors were assigned to the other parameters (a, b, sd.litter.dev, sd.embryo.dev) allowing variation of each within a realistic range (see Supplementary Methods Table 1)
Monte Carlo Markov-Chain (MCMC) techniques were used to estimate the joint posterior distribution of parameters from prior distributions and data. Computations were performed using the JAGS software via the R package rjags (Plummer et al., 2016) (a runnable R script is provided in Supplementary material with data corresponding to strain FVB). Three independent MCMC chains were run in parallel. For each chain, 110,000 samples were produced. The first 10,000 were considered as burn-in phase and discarded. To avoid autocorrelation, the remaining 100,000 samples were thinned by selecting one out of 20 samples, thus keeping 5000 samples per chain. We checked the convergence again by displaying MCMC chain traces and autocorrelation plots and by computing the Gelman and Rubin’s statistics as modified by Brooks and Gelman (Brooks and Gelman, 1998). For each parameter, its point estimate was defined as the median of its marginal posterior distribution, and the 95% credible interval was defined from the 2.5 and 97.5 percentiles of this distribution. The calculation of the development age for each embryo was integrated in the model was thus estimated in the same way from its posterior distribution estimated by MCMC.


## usage

To run the example, you must 
1/ install R (http://cran.r-project.org/)
2/ install JAGS (http://mcmc-jags.sourceforge.net/)
3/ install the R package rjags within R 
4/ put the file fit_model.R and data_FVB.txt in a directory that you will 
	define as the working directory in R
5/ run the script written in file fit_model.R 
(it may take a few minutes but if you only want a rough estimation
you can divide n.iter and thin by 10 to speed up computations)

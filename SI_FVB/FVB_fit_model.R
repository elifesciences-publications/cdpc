

/*
   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */


########################################################################
# Fit of the  model 
#
# Beware to install JAGS on your computer and in a second time
# the R package rjags, before using the script.
# 
# Author : Marie Laure Delignette-Muller and Marie Semon
# Email : marielaure.delignettemuller@vetagro-sup.fr
########################################################################

# Importation and visualisation of data
#######################################
rawdata <- read.table("data_FVB.txt", header = TRUE)
n <- nrow(rawdata)
ref <- 15 # value near the mean of age.in.dpc used to center data
# just to limit correlation between intercept and slope of the linear fitted relation
rawdata$age.in.dpc.c <- rawdata$age.in.dpc - ref 
plot(weight.in.log ~ age.in.dpc.c, data = rawdata,
     xlab = "litter age in dpc",
     ylab = "embryo weight in log")

# Creation of argument data for jags.model
##########################################
embryo.nb <- nrow(rawdata) # number of embryos
litter.levels <-  levels(rawdata$litter) # names of the litters
litter.nb <- length(litter.levels) # number of litters
litter.index <- match(rawdata$litter, litter.levels) # index of the litter for each embryo 

data4jags <- list(age.in.dpc.c = rawdata$age.in.dpc.c,
                  weight.in.log = rawdata$weight.in.log, 
                  embryo.nb = embryo.nb,
                  litter.nb = litter.nb, 
                  litter.index = litter.index, 
                  sd.litter.preg = 0.05, # 0.05 or 0.1 depending of the scenario
                  ref = ref)

# Bayesian estimation using JAGS
###################################
require(rjags)

# Description of the model in JAGS language
model <-
  "model
{
  # loop on litters
  for(i in 1:litter.nb)
  {
  eps.litter.preg[i] ~ dnorm(0, tau.litter.preg)
  eps.litter.dev[i] ~ dnorm(0, tau.litter.dev)
  }
  
  # loop on embryos
  for(j in 1:embryo.nb)
  {
  weight.in.log[j] ~ dnorm(a + b * age.in.dpc.c[j] + 
  eps.litter.preg[litter.index[j]] + 
  eps.litter.dev[litter.index[j]],
  tau.embryo.dev)
  
  age.dev[j] <- (weight.in.log[j] - eps.litter.preg[litter.index[j]] - a) / b + ref
  }
  
  # definition of priors
  sd.litter.dev ~ dunif(0, 1)
  sd.embryo.dev ~ dunif(0, 1)
  
  tau.litter.preg <- 1/sd.litter.preg^2  
  tau.litter.dev <- 1/sd.litter.dev^2  
  tau.embryo.dev <- 1/sd.embryo.dev^2  
  
  a ~ dunif(0, 10)
  b ~ dunif(0, 1)
  
}"

  # Inference using MCMC algorithm
  ###################################
  m <- jags.model(file = textConnection(model), data = data4jags, n.chains = 3)
  update(m, n.iter = 10000) # burnin
  
  # This step can take a few minutes 
  # you might divide by 10 n.iter and thin to obtain rough results
  mcmc <- coda.samples(m, c("a", "b", "sd.litter.dev", "sd.embryo.dev"), 
                       n.iter = 100000, thin = 20) 
  # Check of the convergence 
  ##########################
  gelman.diag(mcmc) 
  plot(mcmc, trace = TRUE, density = FALSE)
  
  # Parameter estimations
  #######################
  summary(mcmc)
  plot(mcmc, trace= FALSE, density = TRUE)
  
  # Estimation of developmental ages
  #######################
  # This step can take a few minutes 
  # you might divide by 10 n.iter and thin to obtain rough results
  mcmc.dev <- coda.samples(m, c("age.dev"), n.iter = 100000, thin = 20) 
  age.dev.estimation <- summary(mcmc.dev)$quantiles
  
  # Storage of estimations as median and 2.5 and 97.5 quantiles for each embryo
  age.dev.median <- age.dev.estimation[,3]
  age.dev.lower <- age.dev.estimation[,1]
  age.dev.upper <- age.dev.estimation[,5]
  
  # plot of estimation of age.dev as a function of weight.in.log
  pdf("FVB_reasonnable.pdf")
  plot(rawdata$weight.in.log, age.dev.median, 
       xlab = "embryo weight in log", 
       ylab = "estimated developmental age with 95% credible interval")
  segments(x0 = rawdata$weight.in.log,
           y0 = age.dev.lower, 
           x1 = rawdata$weight.in.log,
           y1 = age.dev.upper,
           col = "grey")
  points(rawdata$weight.in.log, age.dev.median, pch = 16) # to plot again points above 95% credible intervals
  dev.off()

  
  pdf("FVB_reasonnable_unlog.pdf")
  plot(exp(rawdata$weight.in.log), age.dev.median, 
       xlab = "embryo weight", 
       ylab = "estimated developmental age with 95% credible interval")
 
  points(rawdata$weight.in.log, age.dev.median, pch = 16) # to plot again points above 95% credible intervals
  dev.off()
  

# The code below can be used to reproduce the results of
# An information-theoretic approach for selecting arms in clinical trials
# by P. Mozgunov and T. Jaki (2019)
# Section 4.4 with the minimum efficacy bounds and the numerical intergration

# The "wdesign-numerical-code.R" is required to be run prior to the code below.

# Trial 1, N=423, Select-The-Best Allocation Rule (warning: long running time)
# Under the Null Hypothesis
design<-wdesign.ph2.numerical(true=c(0.3,0.3,0.3,0.3),target=0.999,n=423,prior=rep(0.99,4),
                     test="Fisher",kappa=0.50,cut.off.typeI=0.0587,beta=c(5,2,2,2),nsims=2500,
                     minimum.efficacy=0.40) # choice of the minimum efficacy
# Type I Error
design$TypeI.error
# ENS and corresponding SE
design$ENS
design$SE.ENS
# Experimentation Proportion and corresponding SE
design$Experimentation
design$Experimentation.SE

### Under the Alternative
design<-wdesign.ph2.numerical(true=c(0.3,0.3,0.3,0.5),target=0.999,n=423,prior=rep(0.99,4),
                              test="Fisher",kappa=0.50,cut.off.typeI=0.0587,beta=c(5,2,2,2),nsims=2500,
                              minimum.efficacy=0.40) # choice of the minimum efficacy

# Power
design$Power
# ENS and corresponding SE
design$ENS
design$SE.ENS
# Experimentation Proportion and corresponding SE
design$Experimentation
design$Experimentation.SE

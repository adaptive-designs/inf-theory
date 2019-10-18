# The code below can be used to reproduce the results of
# An information-theoretic approach for selecting arms in clinical trials
# by P. Mozgunov and T. Jaki (2019)
# Section 5 with co-primary efficacy endpoints

# The "wdesign-co-primary-code.R" is required to be run prior to the code below.


# WE-II Select-the-Best
# Kappa 0.54 - Robust Optimal Value under the ENS Objective Function


# Under the Null Hypothesis
design<-wdesign.co.primary.ph2(true1=c(0.10,0.10,0.10),true2=c(0.45,0.45,0.45),target1=0.999,target2=0.999,n=165,
                               prior1=rep(0.99,3),prior2=rep(0.99,3),beta1=c(5,2,2),control=1,cut.off.typeI=0.095,alternative="greater",
                               correlation=0.75,assignment="best",
                               beta2=c(5,2,2),kappa=0.54,nsims=10000,hypothesis=T,test="Fisher")
# Type I Error
design$TypeI.error
# ENS and corresponding SE
design$ENS
design$SE.ENS
design$ENS2
design$SE.ENS2
# Experimentation Proportion and corresponding SE
design$Experimentation
design$Experimentation.SE


# Under the Alternative Hypothesis
design<-wdesign.co.primary.ph2(true1=c(0.10,0.10,0.25),true2=c(0.45,0.45,0.60),target1=0.999,target2=0.999,n=165,
                               prior1=rep(0.99,3),prior2=rep(0.99,3),beta1=c(5,2,2),control=1,cut.off.typeI=0.095,alternative="greater",
                               correlation=0.75,assignment="best",
                               beta2=c(5,2,2),kappa=0.54,nsims=10000,hypothesis=T,test="Fisher")


# Power
design$Power
# ENS and corresponding SE
design$ENS
design$SE.ENS
# Experimentation Proportion and corresponding SE
design$Experimentation
design$Experimentation.SE




# Kappa 0.69 - Robust Optimal Value under the Power Objective Function

# Under the Null Hypothesis
design<-wdesign.co.primary.ph2(true1=c(0.10,0.10,0.10),true2=c(0.45,0.45,0.45),target1=0.999,target2=0.999,n=165,
                               prior1=rep(0.99,3),prior2=rep(0.99,3),beta1=c(5,2,2),control=1,cut.off.typeI=0.0980,alternative="greater",
                               correlation=0.75,assignment="best",
                               beta2=c(5,2,2),kappa=0.69,nsims=10000,hypothesis=T,test="Fisher")
# Type I Error
design$TypeI.error
# ENS and corresponding SE
design$ENS
design$SE.ENS
design$ENS2
design$SE.ENS2
# Experimentation Proportion and corresponding SE
design$Experimentation
design$Experimentation.SE


# Under the Alternative Hypothesis
design<-wdesign.co.primary.ph2(true1=c(0.10,0.10,0.25),true2=c(0.45,0.45,0.60),target1=0.999,target2=0.999,n=165,
                               prior1=rep(0.99,3),prior2=rep(0.99,3),beta1=c(5,2,2),control=1,cut.off.typeI=0.0980,alternative="greater",
                               correlation=0.75,assignment="best",
                               beta2=c(5,2,2),kappa=0.69,nsims=10000,hypothesis=T,test="Fisher")

# Power
design$Power
# ENS and corresponding SE
design$ENS
design$SE.ENS
# Experimentation Proportion and corresponding SE
design$Experimentation
design$Experimentation.SE

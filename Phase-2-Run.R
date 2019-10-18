# The code below can be used to reproduce the results of
# An information-theoretic approach for selecting arms in clinical trials
# by P. Mozgunov and T. Jaki (2019)
# Section 4 using the leading term of the asymptotic criterion

# The "wdesign-main-code.R" is required to be run prior to the code below.

######### Trial 1, N=423
######### WE-II Design; Select-the-Best Allocation Rule;

######### Kappa = 0.56 - Robust Optimal Value for ENS Objective Function
### Under the Null Hypothesis
set.seed(100)
design<-wdesign.ph2(true=c(0.3,0.3,0.3,0.3),target=0.999,n=423,prior=rep(0.99,4),
                     beta=c(5,2,2,2),kappa=0.56,nsims=10000,cut.off.typeI=0.0507,test="Fisher")
# Type I Error
design$TypeI.error
# ENS and corresponding SE
design$ENS
design$SE.ENS
# Experimentation Proportion and corresponding SE
design$Experimentation
design$Experimentation.SE


### Under the Alternative
set.seed(100)
design<-wdesign.ph2(true=c(0.3,0.3,0.3,0.5),target=0.999,n=423,prior=rep(0.99,4),
                     beta=c(5,2,2,2),kappa=0.56,nsims=10000,cut.off.typeI=0.0507,test="Fisher")
# Power
design$Power
# ENS and corresponding SE
design$ENS
design$SE.ENS
# Experimentation Proportion and corresponding SE
design$Experimentation
design$Experimentation.SE


######### Kappa = 0.65 - Robust Optimal Value for Power Objective Function
### Under the Null Hypothesis
set.seed(100)
design<-wdesign.ph2(true=c(0.3,0.3,0.3,0.3),target=0.999,n=423,prior=rep(0.99,4),
                    beta=c(5,2,2,2),kappa=0.65,nsims=10000,cut.off.typeI=0.059,test="Fisher")
# Type I Error
design$TypeI.error
# ENS and corresponding SE
design$ENS
design$SE.ENS
# Experimentation Proportion and corresponding SE
design$Experimentation
design$Experimentation.SE

### Under the Alternative
set.seed(100)
design<-wdesign.ph2(true=c(0.3,0.3,0.3,0.5),target=0.999,n=423,prior=rep(0.99,4),
                    beta=c(5,2,2,2),kappa=0.65,nsims=10000,cut.off.typeI=0.059,test="Fisher")

# Power
design$Power
# ENS and corresponding SE
design$ENS
design$SE.ENS
# Experimentation Proportion and corresponding SE
design$Experimentation
design$Experimentation.SE


######### WE-I Design; Randomization Allocation Rule;
### Under the Null Hypothesis
set.seed(100)
design<-wdesign.ph2(true=c(0.3,0.3,0.3,0.3),target=0.999,n=423,prior=rep(0.99,4),assignment="randomization",
                    beta=c(5,2,2,2),kappa=0.50,nsims=10000,cut.off.typeI=0.0824,test="Fisher")
# Type I Error
design$TypeI.error
# ENS and corresponding SE
design$ENS
design$SE.ENS
# Experimentation Proportion and corresponding SE
design$Experimentation
design$Experimentation.SE

### Under the Alternative
set.seed(100)
design<-wdesign.ph2(true=c(0.3,0.3,0.3,0.5),target=0.999,n=423,prior=rep(0.99,4),assignment="randomization",
                    beta=c(5,2,2,2),kappa=0.50,nsims=10000,cut.off.typeI=0.0824,test="Fisher")

# Power
design$Power
# ENS and corresponding SE
design$ENS
design$SE.ENS
# Experimentation Proportion and corresponding SE
design$Experimentation
design$Experimentation.SE




######### Trial 2, N=80
######### WE-II Design; Select-the-Best Allocation Rule;

######### Kappa = 0.51 - Robust Optimal Value for ENS Objective Function

### Under the Null Hypothesis
set.seed(100)
design<-wdesign.ph2(true=c(0.3,0.3,0.3,0.3),target=0.999,n=80,prior=rep(0.99,4),
                     beta=c(5,2,2,2),kappa=0.51,nsims=10000,cut.off.typeI=0.1373,test="Fisher")
# Type I Error
design$TypeI.error
# ENS and corresponding SE
design$ENS
design$SE.ENS
# Experimentation Proportion and corresponding SE
design$Experimentation
design$Experimentation.SE

### Under the Alternative
set.seed(100)
design<-wdesign.ph2(true=c(0.3,0.4,0.5,0.6),target=0.999,n=80,prior=rep(0.99,4),
                    beta=c(5,2,2,2),kappa=0.51,nsims=10000,cut.off.typeI=0.1373,test="Fisher")

# Power
design$Power
# ENS and corresponding SE
design$ENS
design$SE.ENS
# Experimentation Proportion and corresponding SE
design$Experimentation
design$Experimentation.SE

######### Kappa = 0.73 - Robust Optimal Value for Power Objective Function
### Under the Null Hypothesis
set.seed(100)
design<-wdesign.ph2(true=c(0.3,0.3,0.3,0.3),target=0.999,n=80,prior=rep(0.99,4),
                     beta=c(5,2,2,2),kappa=0.73,nsims=10000,cut.off.typeI=0.1265,test="Fisher")

# Type I Error
design$TypeI.error
# ENS and corresponding SE
design$ENS
design$SE.ENS
# Experimentation Proportion and corresponding SE
design$Experimentation
design$Experimentation.SE


### Under the Alternative
set.seed(100)
design<-wdesign.ph2(true=c(0.3,0.4,0.5,0.6),target=0.999,n=80,prior=rep(0.99,4),
                    beta=c(5,2,2,2),kappa=0.73,nsims=10000,cut.off.typeI=0.1265,test="Fisher")
# Power
design$Power
# ENS and corresponding SE
design$ENS
design$SE.ENS
# Experimentation Proportion and corresponding SE
design$Experimentation
design$Experimentation.SE


######### WE-I Design; Randomization Allocation Rule;
### Under the Null Hypothesis
set.seed(100)
design<-wdesign.ph2(true=c(0.3,0.3,0.3,0.3),target=0.999,n=80,prior=rep(0.99,4),assignment="randomization",
                    beta=c(5,2,2,2),kappa=0.50,nsims=10000,cut.off.typeI=0.0874,test="Fisher")

# Type I Error
design$TypeI.error
# ENS and corresponding SE
design$ENS
design$SE.ENS
# Experimentation Proportion and corresponding SE
design$Experimentation
design$Experimentation.SE


### Under the Alternative
set.seed(100)
design<-wdesign.ph2(true=c(0.3,0.4,0.5,0.6),target=0.999,n=80,prior=rep(0.99,4),assignment="randomization",
                    beta=c(5,2,2,2),kappa=0.50,nsims=10000,cut.off.typeI=0.0874,test="Fisher")

# Power
design$Power
# ENS and corresponding SE
design$ENS
design$SE.ENS
# Experimentation Proportion and corresponding SE
design$Experimentation
design$Experimentation.SE

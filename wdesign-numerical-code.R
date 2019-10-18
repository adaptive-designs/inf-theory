loss.uni<-function(x,target.tox,n=1,kappa=0.5,minimum.efficacy=0){
  
  target<-target.tox
  minimum<-minimum.efficacy
  
  integrand <- function(p) {-1/beta(x+1,n-x+1)*p^x*(1-p)^(n-x)*(log(1/beta(x+1,n-x+1))+x*log(p)+(n-x)*log(1-p))}
  integrand.2 <- function(p) {(-1)*as.numeric((p>minimum))*p^(target*n^kappa)*(1-p)^((1-target)*n^kappa)*p^x*(1-p)^(n-x)*(log(1/beta(x+1,n-x+1))+x*log(p)+(n-x)*log(1-p))}
  integrand.norm <- function(p) {as.numeric((p>minimum))*p^(target*n^kappa)*(1-p)^((1-target)*n^kappa)*p^x*(1-p)^(n-x)}
  
  standard<-integrate(integrand, lower = 0.01, upper =0.99)$value
  weighted<-integrate(integrand.2, lower = 0.01, upper =0.99)$value/integrate(integrand.norm, lower = 0.01, upper =0.99)$value
  y<-weighted-standard
  y<-10^(weighted-standard)
  return(y)}
library("pwr")
median.beta<-function(x,n){y<-x/n
return(y)}

wdesign.ph2.numerical<-function(true,target,n,cohort=1,assignment="best",prior=NULL,beta=1,kappa=0.5,nsims,minimum.efficacy=0,
                      hypothesis=T,alternative="greater",cut.off.typeI=0.05,control=1,test="Fisher"){
  typeI<-cut.off.typeI
  fish.mat<-mat.or.vec(2,2)
  N<-round(n/cohort)                                        
  M<-length(true)
  M2<-M-1
  p.values<-n.outcomes<-x.outcomes<-mat.or.vec(M,1)
  all.p.values<-mat.or.vec(nsims,M)
  cutoff<-typeI/(M2)
  z.norm<-qnorm(1-cutoff, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  rec.all<-result.exp<-mat.or.vec(nsims,M)
  power.final<-result.suc<-bias.final<-mat.or.vec(nsims,1)
  suc.cum<-suc<-exp<-exp.cum<-mat.or.vec(N+1,M)
  experiment<-array(0,dim=c(N+1,M,nsims))
  for (z in 1:nsims){losses<-mat.or.vec(N+1,M)
  exp[1,]<- exp.cum[1,]<-beta
  suc[1,]<-suc.cum[1,]<-prior*beta
  for (q in 1:M){
    losses[1,q]<-loss.uni(x=suc.cum[1,q],target.tox=target,n=exp.cum[1,q],kappa=kappa,minimum.efficacy=minimum.efficacy)
  }
  if(assignment=="randomization"){
    nextdose<-sample(1:M, cohort, prob = (1/(losses[1,]))/sum(1/(losses[1,])),replace = TRUE)}
  else{nextdose<-which.min(losses[1,])}
  suc[2,]<-exp[2,]<-0
  exp[2,nextdose]<-cohort
  exp.cum[2,]<-exp.cum[1,]
  exp.cum[2,nextdose]<-exp.cum[1,nextdose]+exp[2,nextdose]
  suc[2,nextdose]<-sum(rbinom(cohort,1,true[nextdose]))
  suc.cum[2,]<-suc.cum[1,] 
  suc.cum[2,nextdose]<-suc.cum[1,nextdose]+suc[2,nextdose]
  j<-2
  while (j<N+1){
    
    for ( q in 1:M){
      losses[j,q]<-loss.uni(x=suc.cum[j,q],target.tox=target,n=exp.cum[j,q],kappa=kappa,minimum.efficacy=minimum.efficacy)
      
    }
    if(assignment=="randomization"){
      nextdose<-sample(1:M, cohort, prob = (1/(losses[j,]))/sum(1/(losses[j,])),replace = TRUE)}
    else{nextdose<-which.min(losses[j,])}
    exp[j+1,]<-suc[j+1,]<-0
    exp[j+1,nextdose]<-cohort
    exp.cum[j+1,]<-exp.cum[j,]
    exp.cum[j+1,nextdose]<-exp.cum[j,nextdose]+exp[j+1,nextdose]
    suc[j+1,nextdose]<-sum(rbinom(cohort,1,true[nextdose]))
    suc.cum[j+1,]<-suc.cum[j,] 
    suc.cum[j+1,nextdose]<-suc.cum[j,nextdose]+suc[j+1,nextdose]
    j<-j+1}
  result.exp[z,]<-(exp.cum[N+1,]-beta)
  result.suc[z]<-sum(suc.cum[N+1,]-suc.cum[1,])
  j<-N+1
  for(q in 1:M){
    losses[j,q]<-loss.uni(x=suc.cum[j,q]-suc.cum[1,q],target.tox=target,n=exp.cum[j,q]-exp.cum[1,q],kappa=0.5,minimum.efficacy=minimum.efficacy)
  }
  nextdose<-which.min(losses[j,])
  rec.all[z,nextdose]<-1
  bias.final[z]<-(median.beta(suc.cum[j,control]-suc.cum[1,control],exp.cum[j,control]-exp.cum[1,control])-median.beta(suc.cum[j,which.max(true)]-suc.cum[1,which.max(true)],exp.cum[j,which.max(true)]-exp.cum[1,which.max(true)]))-(true[control]-true[which.max(true)])
  if (hypothesis==T){for (q in 1:M){x.outcomes[q]<-as.integer(suc.cum[j,q]-suc.cum[1,q])
  n.outcomes[q]<-as.integer(exp.cum[j,q]-exp.cum[1,q])}        
    if(all(true==true[control])){
      if(test=="Fisher"){for (q in 1:M){if(q!=control){fish.mat[1,1]<-x.outcomes[q]
      fish.mat[1,2]<-n.outcomes[q]-x.outcomes[q]
      fish.mat[2,1]<-x.outcomes[control]
      fish.mat[2,2]<-n.outcomes[control]-x.outcomes[control]
      my<-fisher.test(fish.mat, alternative = alternative)
      p.values[q]<-all.p.values[z,q]<-my$p.value
      
      }else{p.values[q]<-0}}}
      
      if(test=="Normal"){z.norm<-qnorm(1-(cutoff), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
      for (q in 1:M){if(q!=control){p.best<-x.outcomes[q]/n.outcomes[q]
      p0<-x.outcomes[control]/n.outcomes[control]
      pf<-(x.outcomes[control]+x.outcomes[q])/(n.outcomes[q]+n.outcomes[control])
      denom<-sqrt(pf*(1-pf)*(1/n.outcomes[control]+1/n.outcomes[q]))
      p.values[q]<-(p.best-p0)/denom}else{p.values[q]<-0}}}}else{better<-as.vector(which(true>true[control]))
      if(test=="Fisher"){for (q in 1:length(better)){fish.mat[1,1]<-x.outcomes[better[q]]
      fish.mat[1,2]<-n.outcomes[better[q]]- x.outcomes[better[q]]
      fish.mat[2,1]<-x.outcomes[control]
      fish.mat[2,2]<-n.outcomes[control]-x.outcomes[control]
      my<-fisher.test(fish.mat, alternative = alternative)
      p.values[q]<-my$p.value}}
      if(test=="Normal"){for (q in 1:length(better)){p.best<-x.outcomes[better[q]]/n.outcomes[better[q]]
      p0<-x.outcomes[control]/n.outcomes[control]
      pf<-(x.outcomes[control]+x.outcomes[better[q]])/(n.outcomes[better[q]]+n.outcomes[control])
      denom<-sqrt(pf*(1-pf)*(1/n.outcomes[control]+1/n.outcomes[better[q]]))
      p.values[q]<-(p.best-p0)/denom}}}
    if(test=="Fisher"){if(any(p.values[p.values!=0]<cutoff)){power.final[z]<-1}else{power.final[z]<-0}}
    if(test=="Normal"){if(any(p.values[p.values!=0]>z.norm)){power.final[z]<-1}else{power.final[z]<-0}}}}
  y<-colSums(rec.all)/nsims
  var.pavel<-mat.or.vec(M,1)
  for (u in 1:M){var.pavel[u]<-var(result.exp[,u]/n)}
  if(all(true==true[control])){output<-list(True.Probabilities=true,number.of.simulation=nsims,Sample.Size=mean(rowSums(result.exp)),
                                            Experimentation=colSums(result.exp)/(n*nsims),
                                            Experimentation.SE=sqrt(var.pavel),Selections=y,ENS=mean(result.suc),SE.ENS=sqrt(var(result.suc)),
                                            TypeI.error=mean(power.final),Bias=mean(bias.final),All.P.Values=all.p.values)}else{
                                              output<-list(True.Probabilities=true,number.of.simulation=nsims,Sample.Size=mean(rowSums(result.exp)),
                                                           Experimentation=colSums(result.exp)/(n*nsims),
                                                           Experimentation.SE=sqrt(var.pavel),Selections=y,ENS=mean(result.suc),SE.ENS=sqrt(var(result.suc)),
                                                           Power=mean(power.final),Bias=mean(bias.final))}
  return(output)}



#####         + ADD FIX RANDOMIZATION RESULT #####






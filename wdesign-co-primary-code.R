library("pwr")
library("mvtnorm")

median.beta<-function(x,n){y<-x/n
return(y)}
loss.uni<-function(p1,p2,target1,target2,n=1,kappa=0.5){
  alpha1<-p1*p2
  alpha2<-p1*(1-p2)
  alpha3<-(1-p1)*p2
  alpha4<-1-alpha1-alpha2-alpha3
  theta1<-target1*target2
  theta2<-target1*(1-target2)
  theta3<-(1-target1)*target2
  theta4<-1-theta1-theta2-theta3
  y1<-theta1^2/alpha1+theta2^2/alpha2+theta3^2/alpha3+theta4^2/alpha4-1
  y<-(y1)*n^(2*kappa-1)
  return(y)}

wdesign.co.primary.ph2<-function(true1,true2,target1,target2,correlation=0,n,cohort=1,assignment="best",prior1=NULL,prior2=NULL,beta1=1,beta2=1,kappa=0.5,nsims,
                      hypothesis=T,alternative="greater",cut.off.typeI=0.05,control=1,test="Fisher"){
  typeI<-cut.off.typeI
  fish.mat1<-mat.or.vec(2,2)
  fish.mat2<-mat.or.vec(2,2)
  N<-round(n/cohort)                                        
  M<-length(true1)
  M2<-M-1
  cutoff<-typeI/(M2)/2
  z.norm<-qnorm(1-cutoff, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  rec.all<-result.exp1<-result.exp2<-mat.or.vec(nsims,M)
  power.final<-result.suc1<-bias.final1<-result.suc2<-bias.final2<-mat.or.vec(nsims,1)
  suc.cum1<-suc1<-exp1<-exp.cum1<-suc.cum2<-suc2<-exp2<-exp.cum2<-mat.or.vec(N+1,M)
  experiment1<-array(0,dim=c(N+1,M,nsims))
  experiment2<-array(0,dim=c(N+1,M,nsims))
  
  for (z in 1:nsims){
    p.values<-mat.or.vec(2,M)
    losses<-mat.or.vec(N+1,M)
  exp1[1,]<- exp.cum1[1,]<-beta1
  exp2[1,]<- exp.cum2[1,]<-beta2
  suc1[1,]<-suc.cum1[1,]<-prior1*beta1
  suc2[1,]<-suc.cum2[1,]<-prior2*beta2
  p.est1<-median.beta(suc.cum1[1,],exp.cum1[1,])
  p.est2<-median.beta(suc.cum2[1,],exp.cum2[1,])
  suc.cum.new<-mat.or.vec(4,M)
  losses[1,]<-loss.uni(p1=p.est1,p2=p.est2,target1=target1,target2=target2,n=1,kappa)
  if(assignment=="randomization"){
    nextdose<-sample(1:M, cohort, prob = (1/(losses[1,]))/sum(1/(losses[1,])),replace = TRUE)}
  else{nextdose<-which.min(losses[1,])}
  suc1[2,]<-exp1[2,]<-0
  suc2[2,]<-exp2[2,]<-0
  exp1[2,nextdose]<-cohort
  exp2[2,nextdose]<-cohort
  exp.cum1[2,]<-exp.cum1[1,]
  exp.cum2[2,]<-exp.cum2[1,]
  exp.cum1[2,nextdose]<-exp.cum1[1,nextdose]+exp1[2,nextdose]
  exp.cum2[2,nextdose]<-exp.cum2[1,nextdose]+exp2[2,nextdose]
  
  y<-pnorm(rmvnorm(cohort,mean=c(0,0),sigma=rbind(c(1,correlation),c(correlation,1))))
  
  if(y[1]<true1[nextdose] & y[2]<true2[nextdose]){
    response<-1
  }else{
    if(y[1]<true1[nextdose] & y[2]>true2[nextdose]){
      response<-2
    }else{
      if(y[1]>true1[nextdose] & y[2]<true2[nextdose]){
        response<-3
      }else{
        response<-4
      }
    }
  }
  
  suc.cum.new[response,nextdose]<-suc.cum.new[response,nextdose]+1
  
  suc1[2,nextdose]<-sum(y[1]<true1[nextdose])
  suc2[2,nextdose]<-sum(y[2]<true2[nextdose])
  
  suc.cum1[2,]<-suc.cum1[1,] 
  suc.cum2[2,]<-suc.cum2[1,] 
  
  
  suc.cum1[2,nextdose]<-suc.cum1[1,nextdose]+suc1[2,nextdose]
  suc.cum2[2,nextdose]<-suc.cum2[1,nextdose]+suc2[2,nextdose]
  
  j<-2
  while (j<N+1){
    p.est1<-median.beta(suc.cum1[j,],exp.cum1[j,])
    p.est2<-median.beta(suc.cum2[j,],exp.cum2[j,])
    # cat("probability 1=",p.est1,"\n")
    # cat("probability 2=",p.est2,"\n")
    
    losses[j,]<-loss.uni(p1=p.est1,p2=p.est2,target1=target1,target2=target2,n=exp.cum1[j,],kappa)
    # cat(losses[j,],"\n")
    if(assignment=="randomization"){
      nextdose<-sample(1:M, cohort, prob = (1/(losses[j,]))/sum(1/(losses[j,])),replace = TRUE)}
    else{nextdose<-which.min(losses[j,])}
    exp1[j+1,]<-suc1[j+1,]<-0
    exp2[j+1,]<-suc2[j+1,]<-0
    
    exp1[j+1,nextdose]<-cohort
    exp2[j+1,nextdose]<-cohort
    
    exp.cum1[j+1,]<-exp.cum1[j,]
    exp.cum2[j+1,]<-exp.cum2[j,]
    
    exp.cum1[j+1,nextdose]<-exp.cum1[j,nextdose]+exp1[j+1,nextdose]
    exp.cum2[j+1,nextdose]<-exp.cum2[j,nextdose]+exp2[j+1,nextdose]
    
    y<-pnorm(rmvnorm(cohort,mean=c(0,0),sigma=rbind(c(1,correlation),c(correlation,1))))
    
    # if(y[1]<true1[nextdose] & y[2]<true2[nextdose]){
    #   response<-1
    # }else{
    #   if(y[1]<true1[nextdose] & y[2]>true2[nextdose]){
    #     response<-2
    #   }else{
    #     if(y[1]>true1[nextdose] & y[2]<true2[nextdose]){
    #       response<-3
    #     }else{
    #       response<-4
    #     }
    #   }
    # }
    
    # suc.cum.new[response,nextdose]<-suc.cum.new[response,nextdose]+1
    
    
    suc1[j+1,nextdose]<-sum(y[1]<true1[nextdose])
    suc2[j+1,nextdose]<-sum(y[2]<true2[nextdose])
    
    suc.cum1[j+1,]<-suc.cum1[j,] 
    suc.cum2[j+1,]<-suc.cum2[j,] 
    
    suc.cum1[j+1,nextdose]<-suc.cum1[j,nextdose]+suc1[j+1,nextdose]
    suc.cum2[j+1,nextdose]<-suc.cum2[j,nextdose]+suc2[j+1,nextdose]
    
    j<-j+1}
  result.exp1[z,]<-(exp.cum1[N+1,]-beta1)
  result.exp2[z,]<-(exp.cum2[N+1,]-beta2)
  
  result.suc1[z]<-sum(suc.cum1[N+1,]-suc.cum1[1,])
  result.suc2[z]<-sum(suc.cum2[N+1,]-suc.cum2[1,])
  
  j<-N+1
  p.est1<-median.beta(suc.cum1[j,]-suc.cum1[1,],exp.cum1[j,]-exp.cum1[1,])
  p.est2<-median.beta(suc.cum2[j,]-suc.cum2[1,],exp.cum2[j,]-exp.cum2[1,])
  # cat("probability 1=",p.est1,"\n")
  # cat("probability 2=",p.est2,"\n")
  
  losses[j,]<-loss.uni(p1=p.est1,p2=p.est2,target1=target1,target2=target2,n=1,kappa=0.5)
  # cat(losses[j,],"\n")
  
  nextdose<-which.min(losses[j,])
  rec.all[z,nextdose]<-1
  bias.final1[z]<-(median.beta(suc.cum1[j,control]-suc.cum1[1,control],exp.cum1[j,control]-exp.cum1[1,control])-median.beta(suc.cum1[j,which.max(true1)]-suc.cum1[1,which.max(true1)],exp.cum1[j,which.max(true1)]-exp.cum1[1,which.max(true1)]))-(true1[control]-true1[which.max(true1)])
  bias.final2[z]<-(median.beta(suc.cum2[j,control]-suc.cum2[1,control],exp.cum2[j,control]-exp.cum2[1,control])-median.beta(suc.cum2[j,which.max(true2)]-suc.cum2[1,which.max(true2)],exp.cum2[j,which.max(true2)]-exp.cum2[1,which.max(true2)]))-(true2[control]-true2[which.max(true2)])
  
  if (hypothesis==T){
    
    if(all(true1==true1[control])){
      if(test=="Fisher"){
        for (q in 1:M){if(q!=control){
      fish.mat1[1,1]<-as.integer(sum(suc.cum1[N+1,q]-suc.cum1[1,q]))
      fish.mat1[1,2]<-as.integer((exp.cum1[N+1,q]-beta1[q])-sum(suc.cum1[N+1,q]-suc.cum1[1,q]))
      fish.mat1[2,1]<-as.integer(sum(suc.cum1[N+1,control]-suc.cum1[1,control]))
      fish.mat1[2,2]<-as.integer((exp.cum1[N+1,control]-beta1[control])-sum(suc.cum1[N+1,control]-suc.cum1[1,control]))
      my<-fisher.test(fish.mat1, alternative = alternative)
      p.values[1,q]<-my$p.value
      fish.mat2[1,1]<-as.integer(sum(suc.cum2[N+1,q]-suc.cum2[1,q]))
      fish.mat2[1,2]<-as.integer((exp.cum2[N+1,q]-beta2[q])-sum(suc.cum2[N+1,q]-suc.cum2[1,q]))
      fish.mat2[2,1]<-as.integer(sum(suc.cum2[N+1,control]-suc.cum2[1,control]))
      fish.mat2[2,2]<-as.integer((exp.cum2[N+1,control]-beta2[control])-sum(suc.cum2[N+1,control]-suc.cum2[1,control]))
      my<-fisher.test(fish.mat2, alternative = alternative)
      p.values[2,q]<-my$p.value
        }else{
        p.values[1:2,q]<-0
        }}}
      
      }else{better<-as.vector(which(true1>true1[control]))
      if(test=="Fisher"){for (q in 1:length(better)){
        
        fish.mat1[1,1]<-as.integer(sum(suc.cum1[N+1,better[q]]-suc.cum1[1,better[q]]))
        fish.mat1[1,2]<-as.integer((exp.cum1[N+1,better[q]]-beta1[better[q]])-sum(suc.cum1[N+1,better[q]]-suc.cum1[1,better[q]]))
        fish.mat1[2,1]<-as.integer(sum(suc.cum1[N+1,control]-suc.cum1[1,control]))
        fish.mat1[2,2]<-as.integer((exp.cum1[N+1,control]-beta1[control])-sum(suc.cum1[N+1,control]-suc.cum1[1,control]))
        my<-fisher.test(fish.mat1, alternative = alternative)
        p.values[1,q]<-my$p.value
        fish.mat2[1,1]<-as.integer(sum(suc.cum2[N+1,better[q]]-suc.cum2[1,better[q]]))
        fish.mat2[1,2]<-as.integer((exp.cum2[N+1,better[q]]-beta2[better[q]])-sum(suc.cum2[N+1,better[q]]-suc.cum2[1,better[q]]))
        fish.mat2[2,1]<-as.integer(sum(suc.cum2[N+1,control]-suc.cum2[1,control]))
        fish.mat2[2,2]<-as.integer((exp.cum2[N+1,control]-beta2[control])-sum(suc.cum2[N+1,control]-suc.cum2[1,control]))
        my<-fisher.test(fish.mat2, alternative = alternative)
        p.values[2,q]<-my$p.value
}}
      
      }
    
    if(test=="Fisher"){
      if(any(p.values[p.values!=0]<cutoff)){
        power.final[z]<-1
      }else{
          power.final[z]<-0
      }
      }
    }
  else{
    power.final[z]<-0
  }
  
  }
  y<-colSums(rec.all)/nsims
  var.pavel<-mat.or.vec(M,1)
  for (u in 1:M){var.pavel[u]<-var(result.exp1[,u]/n)}
  if(all(true1==true1[control])){output<-list(True.Probabilities=true1,number.of.simulation=nsims,Sample.Size=mean(rowSums(result.exp1)),
                                              Experimentation=colSums(result.exp1)/(n*nsims),
                                              Experimentation.SE=sqrt(var.pavel),Selections=y,ENS=mean(result.suc1),SE.ENS=sqrt(var(result.suc1)),ENS2=mean(result.suc2),SE.ENS2=sqrt(var(result.suc2)),
                                              TypeI.error=mean(power.final),Bias=mean(bias.final1))}else{
                                                output<-list(True.Probabilities=true1,number.of.simulation=nsims,Sample.Size=mean(rowSums(result.exp1)),
                                                             Experimentation=colSums(result.exp1)/(n*nsims),
                                                             Experimentation.SE=sqrt(var.pavel),Selections=y,ENS=mean(result.suc1),SE.ENS=sqrt(var(result.suc1)),ENS2=mean(result.suc2),SE.ENS2=sqrt(var(result.suc2)),
                                                             Power=mean(power.final),Bias=mean(bias.final1))}
  return(output)}









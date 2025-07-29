distributive = function(NN,kk=2,KK=4,rate0=0.5,rate1=0.7,rate2=rate0,rate12=rate1,verbose=F,returndesign=F)
 {
  design = data.frame(t(replicate(NN,{vv=rep(0,KK);vv[sample(1:KK,kk,replace=F)]=1;vv}))) # Sampling kk among KK for each patient
  design["pp"] = rate0                                       # Setting default probability of success
  design[design[,"X1"]==1,"pp"] = rate1                      # Setting probability of success for 1 (overwrites previous)
  design[design[,"X2"]==1,"pp"] = rate2                      # Setting probability of success for 2 (overwrites for 1&2)
  design[design[,"X1"]==1 & design[,"X2"]==1,"pp"] = rate12  # Setting probability of success for 1&2 (overwrites again for 1&2)
  design["yy"] = runif(NN)<design[,"pp"]                     # Drawing success
  fitcoefs = summary(glm(paste(sep="","yy ~ ",paste(collapse="+",sep="","X",1:KK)),family=binomial,data=design))$coef # Fitting logistic model
  if(verbose) { print(design) ; print(fitcoefs) }            # Print some output
  if("X1"%in%rownames(fitcoefs)) { pp = fitcoefs["X1",4] } else { pp=NA } # extract p-value of interest, trial considered NS if not estimable
  if(returndesign) {design} else {(!is.na(pp))&(pp<(0.05/KK))} # Output, either trial data (for quality control) or discovery yes/no, with Bonferroni adjustment
  }
distributive(30,verbose=T,returndesign=T) # Giving it a shot
adsasi(distributive,0.9,kk=3,KK=7)        # Doing the empirical determination

######################################################################
######################################################################
# Adaptive simulation wrapper based on user-supplied function
######################################################################
adsasi = function(simfun,tar_power=0.9,...,nsims=5000, verbose=F, impNN=Inf, capNN=2000, trim_initiation = TRUE, savegraphs = FALSE) 
 {
  # simfun for the user-supplied simulation function that takes as first argument a sample size and returns TRUE or FALSE
  # tar_power for desired power
  # ... for additional arguments for simfun (for however the user wrote the latter)
  #    it is better to name these additional arguments to avoid confusion
  # nsims for number of simulations wanted. Because the algorithm just tries stuff around the right sample size, it yields approximately a 
  #    Monte Carlo variance for the power of the sample size it outputs. So if you want less than .5% deviation from 90% power, you need
  #    nsims = .9*.1/(.005^2) = 3600 (or more) simulated trials
  #    The algorithm does the calculations in 10% batches to parallelize things a bit so it will overshoot nsims. 
  # verbose to write some outputs in the console if you fear something is amiss
  # impNN is the impossibility threshold, used to force-exit the algorithm if the answer is more than impNN patients after a few iterations
  # capNN is a cap on the size of the trials the algorithm will simulate
  # trim_initiation is to "forget" the first simulations ; useful if the probit relationship is suspected to be misspecified 
  #    (it will generally be somewhat misspecified if the problem is complex enough to warrant use of this algorithm)
  #    also useful if one gets an early unlucky draw (e.g. a trial with 10k patients that returns FALSE)
  # savegraphs to save the diagnostic graphs instead of displaying them

  # The general way this algorithm works is we are going to compute different sample sizes only a few times and then fit a probit regression
  # this will enable picking a new range of values, and we keep doing this until we have enough simulations

  # Initializing some variables
  tarNN = round(exp(seq(log(10),log(300),length.out=48)))    # tarNN will be a vector of target sample sizes that will be simulated iteratively
  latest_estimate = sqrt(30)                                 # this is the square root of the sample size
  latest_beta=.01                                            # this is the slope nuisance parameter (how fast power drops if sample size is not adequate)
  trials = matrix(c(1,sqrt(1000000),0,1),nrow=2)             # this is where we store the TRUE/FALSE results of our simulations for each tested sample size, 
                                                             #    we are initiating assuming 1 patient is not enough and 1M patients is enough
  colnames(trials) = c("srsampsize","nullreject")            # srsampsize for "Square Root of SAMPle SIZE", nullreject for whether a discovery is made (even 
                                                             #    if it is not strictly a rejection of a well-defined null)
  se = NA                                                    # standard error of the square root of the target sample size using the Hessian matrix
  batch = 0                                                  # just an iteration of simulation batches to keep count
  
  intercept_target = qnorm(tar_power)                        # given the equation we use (see paper), we have an intercept term where we have tar_power = pnorm(intercept_target)

  # Because we use a slightly modified probit regression we have to write the likelihood by hand
  # xx is a scalar vector with two components srbeta and srsampsize. For any pair of values
  #    in xx, the partial likelihood for "successful" and "failed" trials is computed with 
  #    the normal CDF. Then we sum the partial log-likelihoods to get the total log-likelihood. 
  loglik = function(xx) {sum(
                              pnorm(intercept_target+xx["srbeta"]^2*(trials[!!trials[,"nullreject"],"srsampsize"]-xx["X_target"]),lower=T,log=T)
                             ,pnorm(intercept_target+xx["srbeta"]^2*(trials[ !trials[,"nullreject"],"srsampsize"]-xx["X_target"]),lower=F,log=T)
                             )}
  # Here is the main loop
  while(  nrow(trials)<nsims                                 # exit if enough simulations made
        & !(nrow(trials)>200 & latest_estimate>sqrt(impNN))  # exit if after 200+ simulated trials the answers seems to be >impNN
        & !(latest_estimate>314&sum(trials[,1]==sqrt(capNN))>10))  # exit in any case if answer seems >100k patients (likely a user error 
                                                                   #    leading to lots of wasted compute)
   {
    cat("-")                                                 # a homebrew progress bar
    simulations = sapply(tarNN,simfun,...)           # running the simulations, getting a vector of logicals. Note that this calls 
                                                             #    the function from the global environment so it needs to be defined there
    cat("*")                                                 # a homebrew progress bar again
    
    if(nrow(trials)<200&nrow(trials)>2) cat(round(median(tarNN)))     # output the initial tested sample sizes, usually one can see if something has gone wrong early (or if the design is very bad)
    if(nrow(trials)>200&trim_initiation) { trials=trials[-(1:100),] ; trim_initiation = FALSE } # kick out first iterations, deactivate the logical switch
    
    trials = rbind(trials,cbind(srsampsize=sqrt(tarNN),simulations))  # collating the latest simulations with the ones before
    if(savegraphs) { batch=batch+1 ; pngname = paste(sep="","adsasi",gsub(":","-",as.character(Sys.time()),fixed=T),"iteration",batch,".png") ; png(pngname,width=400*2,height=400) } # opening a device to save graphs
    par(mfrow=c(1,2),cex=1.35,cex.lab=1,cex.axis=1)           # light aesthetic setup
    plot(
          trials[,1]^2                                        # sample sizes as ordinates, the abscissae are by default 1:nrow(trials) and that's what we want (to plot them in order)
         ,xlim=c(0,nsims*1.1)
         ,ylim=median(trials[,1])^2 * c(.9,1.11)             # zoomed, not showing everything
         ,main=c("Trace of sample size exploration",paste0("Current estimate ",round(latest_estimate^2)," = (",round(latest_estimate,1),"+-",round(se,2),")²")) # writing the latest estimate and its error
         ,xlab="Simulation#",ylab="Sample size",type="p"
         ,col=paste(sep="",c(rgb(1,.5,0),rgb(0,.6,0.8)),"99")[1+trials[,2]],pch=3,cex=.5+.1*(1-trials[,2]) # successes in blue, failures in orange ; so higher ordinates are bluer, lower are oranger
         )

    
      if(verbose) { print(trials[nrow(trials)-2:0,]) ; cat("\n") ; print(c("beta"=latest_beta, "X_target"=latest_estimate)) }    # verbose descriptions of how the run is going
    initial_value = loglik(c("srbeta"=.01, "X_target"=mean(trials[,"srsampsize"])))    # log-likelihood for the initial values passed to the minimizer 
    if(initial_value!=Inf&initial_value!=-Inf&!is.na(initial_value))    # if the initial value is defined, the minimizer can get to work
     {
      try(mle <- optim(c("srbeta"=.01, "X_target"=mean(trials[,"srsampsize"])), fn=loglik, hessian = TRUE, control=list(fnscale=-1)), silent = TRUE) # using the default optimizer in R, a bit slow but the main compute use if for within-simulation inference anyway
                                                                                                                                                     # control=list(fnscale=-1)) to maximize instead of minimizing
                                                                                                                                                     # here we use try() in case it fails, gets stuck etc. in which case we keep previous values
      latest_beta = unname(mle$par["srbeta"]^2)
      latest_estimate = unname(mle$par["X_target"])
      se = latest_beta*.5
      try(se <- sqrt(diag(solve(-mle$hessian)))["X_target"], silent=TRUE)     # using the numerical Hessian approximation to get a standard error, fails from time to time, most often when not invertible
      tarNN = pmax(3,latest_estimate+c(-1,1)*min(latest_beta*.5,se,na.rm=T))  # now to compute the nest sample sizes to be tried ; first compute a window of +- 1 SE except if SE too large, minimum trial size 9 (3 is for sqrt(9))
      tarNN = round(seq(tarNN[1],tarNN[2],length.out=max(50,round(nrow(trials)/10)))^2)    # now turn those minimal and maximal values into a vector, spanning it evenly (not super elegant to overwrite the same variable)
                                                                                           # the number of simulations for the new batch is 10% of what's already been done, but at least 50
      } else {                                                                # the else loop may rarely be needed to unstick the algorithm, if the optimizer cannot be initiated
              tarNN = round(seq(quantile(trials[,1],.01),quantile(trials[,1],.99),length.out=max(50,round(nrow(trials)/10)))^2)             # we kludge something just to get to the next iteration, mostly around what we have already explored
              latest_estimate = mean(mean(trials[ !trials[,"nullreject"],"srsampsize"]),mean(trials[!!trials[,"nullreject"],"srsampsize"]))
              latest_beta=0.01
              cat("LLfail")                                                   # returning an error (super rare since using the implementation with the squares of sample size and slope)
              }       
    save(tarNN,trials,latest_beta,latest_estimate,file="inner.rda")           # for checking in case of strange behavior ; loading inner.rda puts those in the global environment              
    tarNN = pmin(capNN,tarNN)                                                 # applying the cap ; but if all simulations are made at the same sample size the optimizer cannot compute a slope and fails
    if(all(c(trials[,1]^2,tarNN)==tarNN[1])) tarNN = tarNN*rep(c(.5,1),c(round(length(tarNN)),length(tarNN)-round(length(tarNN))))    # so we divide half of values by 2 if that is the case
      if(verbose) { cat("\n") ; print(tarNN); cat("\n") ; print(se) ; cat("\n") }    # another verbose descriptions of how the run is going
    points(nrow(trials)+1:length(tarNN),tarNN,col="#55555588",pch=1,cex=.8)   # drawing current batch on the graph, result unknown so in grey
    plot(                                                                     # second graph (right panel)
          trials[,1]^2+rnorm(nrow(trials),0,.25)                              # sample size
         ,trials[,2] + (trials[,2]-.5)*2*-rexp(nrow(trials),20)               # success or failure, with some jittering to appreciate density
         ,xlim=latest_estimate^2 * c(.75,1.33)                                # slightly zoomed, not showing everything
         ,ylim=0:1,main=c("Probit regression",paste0("New estimate ",round(latest_estimate^2),", slope ",round(latest_beta,3))),xlab="Sample size",ylab="Expected power"
         ,type="p",col=paste(sep="",c(rgb(1,.5,0),rgb(0,.6,0.8)),"44")[1+trials[,2]],pch=3,cex=.5)
    lines(latest_estimate^2*c(.75,1,1),c(tar_power,tar_power,0),lty=2)
    curve(pnorm(qnorm(tar_power)+latest_beta*(sqrt(x)-latest_estimate)),add=TRUE)
    if(savegraphs) dev.off()                                                  # closing device and saving image if graphs are to be saved
    }                                                                         # end of loop
  
  cat("sample size : ", round(latest_estimate^2), " with ",nrow(trials), " simulated trials", ifelse(nrow(trials)<nsims,"(halted for inefficiency)",""), "\n")
  if( latest_estimate^2 > 4*capNN ) latest_estimate=Inf                       # there is little accuracy when trying to extrapolate that far
  round(latest_estimate^2)
  }                                                                           # end of function

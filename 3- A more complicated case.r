require(MASS)
require(lme4)
require(lmerTest)

# We are basing this example on the TRAIL-1 trial of pirfenidone in idiopathic pulmonary fibrosis (Solomon et al. Lancet Respir Med 2022)
# This is a simple two-arm randomized trial whose composite primary endpoint was not significant
# However, a secondary endpoint of slope of FVC change was significant

# We are trying to power a trial of a drug that would have the same effect as pirfenidone on FVC, as primary endpoint
# FVC or Forced Vital Capacity is the maximum volume of air that a patient can breathe out in one go (after a deep inspiration)
# It is a quantitative descriptor, where higher is better
# FVC can be expressed in liters but also in percent predicted, which is adjusted for age, sex height, weight, ethnicity. 
# This is also a quantitative descriptor where higher is better. 

# We first program the visit structure
visit_times = c(3,6,9,12)                      # These are the visits after baseline, in months
names(visit_times) = paste0("v",visit_times)   # Using some names for quick indexing

# We then program what happens with %predicted Forced Vital Capacity (%pFVC)
# In TRAIL-1, we read that %pFVC at baseline is around 70+-14.5 in both groups (Table 1, "Percent predicted FVC")
# We use a pooled estimate since the trial is randomized. 
# We first create a toy example with 10 patients in the global environment for debugging, then we will implement the same within a function
NN=10
dd = data.frame(patient = paste0("P",1:NN),baseline=rnorm(NN,70,14.5),treated=runif(NN)<0.5)
# patient is the patient id, baseline is the baselien %pFVC, treated is the allocation variable

# We then compute individual deviations from the population mean at each visit. 
# For this, we use a multivariate normal draw with a certain correlation matrix. We have to assume the correlation matrix. 
cormar = rbind(
                c(1,.8,.7,.6)
               ,c(.8,1,.8,.7)
               ,c(.7,.8,1,.8)
               ,c(.6,.7,.8,1)
               )
# We use this structure because for a linear model to make sense, correlation has to be >.5 (otherwise we should just use end of follow-up values, in the same way the variance of a (follow-up minus baseline) t-test is lower than (follow-up only) only if correlation is above .5). 
# Close visits have higher correlation than further visits. 
# These correlation coefficients can also easily be produced by someone who has the right dataset, even a non-statistician. 

# Now that we have specified the correlation matrix we can make the draw. 
dd[names(visit_times)] = mvrnorm(NN,visit_times*0,cormar)

# These deviations can then just be plugged into the population mean and standard deviation at each visit. 
# We try to read population means per visit on figure 2D. 
# Unfortunately, we only have difference from baseline, not absolute values. 
# They seem to increase linearly (this is what is computed below). 
deltas = matrix(c(visit_times/max(visit_times)*-1.02,visit_times/max(visit_times)*-3.21),nrow=2,byrow=T)
rownames(deltas) = c("treated","nontreated")
# We also have SEs here, again for the difference with baseline. They also seem to increase somewhat linearly. 
standevs = matrix(visit_times/max(visit_times)*.515*sqrt(61.5),nrow=1)

# Since we can draw the baseline for our virtual patients, we will simply add their degradation at each visit. 
# This degradation is drawn in an appropriate distribution using the correlation matrix above. 
dd[,names(visit_times)] = round( 
                                dd[,names(visit_times)]*standevs[rep(1,NN),] # We multiply individual standardized deviation by population deviation  (this is a NNx4 matrix)
                                + deltas[2-dd[,"treated"],]                  # Add the treatment arm average degradation for the corresponding visit  (this is a NNx4 matrix)
                                + dd[,rep("baseline",length(visit_times))]   # And add the personal baseline of the patient                           (this is a NNx4 matrix)
                                ,1)
# Note that the correlation matrix is for the individual degradation compared to baseline, not for the absolute values. 
print(dd[1:10,])

# Of course having population means and deviations per visit would be easier. 

# We can also plot and compare to Supplementary Figure S2 (a spaghetti plot)
plot(c(0,max(visit_times)),c(0,120),col=NA,main="Spaghetti plot check",xlab="Time",ylab="%pFVC")
for(ii in 1:nrow(dd)) { lines(c(0,visit_times),dd[ii,c("baseline",names(visit_times))],col=ii) }

# We then put the dataset in long form for the mixed model
ee = data.frame( patient=rep(dd[,"patient"],rep(1+length(visit_times),NN))
                ,baseline=rep(dd[,"baseline"],rep(1+length(visit_times),NN))
                ,time=rep(c(0,visit_times),NN)
                ,treated = rep(dd[,"treated"],rep(1+length(visit_times),NN))
                ,value = c(t(dd[,c("baseline",names(visit_times))]))
                )
print(ee[1:10,])

# We use the following model, with random patient intercepts and slopes
# We also add a time by treatment term (because the treatment acts over time on lung remodelling, not in one shot)
fit = lmer(as.formula("value ~ (time|patient) + time + treated:time"),data=ee,REML=F,control = lmerControl(optimizer ="Nelder_Mead"))
print(fit)
print(anova(fit))

# Now that we have tested things we can remove the toy example and program the function
rm(NN,dd,ee,fit)

# Everything same as above, using the "cormar" global environment variable
simulate_fvc_slope = function(NN)
 {
  dd = data.frame(patient = paste0("P",1:NN),baseline=rnorm(NN,70,14.5),treated=runif(NN)<0.5,mvrnorm(NN,visit_times*0,cormar))
  dd[,names(visit_times)] = round(dd[,names(visit_times)]*standevs[rep(1,NN),] + deltas[2-dd[,"treated"],] + dd[,rep("baseline",length(visit_times))],1)

  ee = data.frame( patient=rep(dd[,"patient"],rep(1+length(visit_times),NN))
                  ,baseline=rep(dd[,"baseline"],rep(1+length(visit_times),NN))
                  ,time=rep(c(0,visit_times),NN)
                  ,treated = rep(dd[,"treated"],rep(1+length(visit_times),NN))
                  ,value = c(t(dd[,c("baseline",names(visit_times))]))
                  )
  ee[,"changefrombaseline"] = ee[,"value"]-ee[,"baseline"]
  fit = lmer(as.formula("value ~ (time|patient) + time + treated:time"),data=ee,REML=F,control = lmerControl(optimizer ="Nelder_Mead"))
  pp = anova(fit)["time:treated","Pr(>F)"]
  output = (!is.na(pp))&pp<0.05
  output
  }

# Try it once
simulate_fvc_slope(50)

# Run adsasi, some fitting problems will happen given the high number of simulations
proposal = adsasi(simulate_fvc_slope)
# This is obviously slower than for the t-test
# Answer is usually aorund 145

# In the supplementary material, the authors give the slopes from their own modelling
# We do not have access to their modelling code & raw data, but we can check a t-test sample size derived from those values against our simulations
# Approximate slope standard deviation from supplementary material (Table S3, the given values are SE not SD)
(0.51*sqrt(63)+0.52*sqrt(60))/2  # = 4.037951
# Approximate standardized effect size, and sample size
power.t.test(NULL,(3.21-1.02)/4.04,sig=.05,pow=.9)$n*2 # = 145

# Very similar, but without their code, it is easier to trust this value after doing our own simulation


# The Adsasi algorithm (ADaptive SAmple SIze) is based on the algorithm in Haviari & Mentré BMC Med Res Methodol 2024, adapted for any clinical trial design & test
# The gist of the approach if that we try different sample sizes empirically and bidn them together with a probit regression (why probit ? read the paper ;)
# The Adsasi algorithm needs a function that takes as argument a number of subjects and outputs TRUE or FALSE per the success/failure of the clinical trial. The user has to write the function. 
# Most of the time, TRUE means rejection of some type of parametric H0 (null hypothesis), but depending on the complexity the user might want to program something more complicated (e.g. hitting multiple endpoints)
# Be careful with more complex stuff though, because power has to be a bijection of sample size for the algorithm to function properly

# For this first example, we are taking a case where the analytical formula for sample size is well known, the t-test
# The implementation uses a function in the global environment for easy fiddling
simulate_one_trial = function(NN,effsize=1) # NN for total sample size, effsize for effect size
 { 
  y0 = rnorm(round(NN/2))           # control arm
  y1 = rnorm(NN-length(y0))         # intervention arm
  y1 = y1 + effsize                 # applying the effect size to the intervention arm
  
  if(NN>=4)                         # t-test not computable for N<4. Because the Adsasi algorithm will try many things, it's usually good practice to hard-code edge cases. 
   {
    pp = t.test(y0,y1)$p.value
    } else {
            pp = 1
            } 
  pp<.05                            # nothing fancy here
  }
cat("One trial with 500 patients : ",simulate_one_trial(500),"\n")      # checking we get an output (should be TRUE unless you are cosmically unlucky)
cat("One trial with 3 patients : ",simulate_one_trial(3),"\n")          # checking we get an output (should be FALSE even if you are cosmically unlucky)
cat( "1000 trials of 60 patients : "
    ,round(100*mean(replicate(1000,simulate_one_trial(60))))            # parallelizing at one size
    ,"% TRUE"
    )
# By the way let's just check the actual (analytical) sample size real quick
power.t.test(NULL,delta=1,sig=.05,power=.9)$n*2 # times 2 because the output is per arm


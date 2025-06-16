# Checking that the base function works
readline(prompt="Press [enter] to start the tutorial")
cat("Let us first check that our individual simulator function works with command : ","\n",sep="")
readline(prompt="Press [enter] to continue")
cat("simulate_one_trial(500)","\n",sep="")
cat(simulate_one_trial(500),"\n",sep="")
cat("simulate_one_trial(500)","\n",sep="")
cat(simulate_one_trial(500),"\n",sep="")
cat("simulate_one_trial(500)","\n",sep="")
cat(simulate_one_trial(500),"\n",sep="")

readline(prompt="Press [enter] to continue")
cat("Now let's do a batch of 1000 trials and compute the observed power with","\n",sep="")
cat("mean(replicate(1000,simulate_one_trial(500)))","\n",sep="")
readline(prompt="Press [enter] to continue")
cat(mean(replicate(1000,simulate_one_trial(500))),"\n",sep="")

readline(prompt="Press [enter] to continue")
cat("Interesting, but we seem to have overshot, and we would like to avoid trial and error","\n",sep="")
readline(prompt="Press [enter] to continue")

# Running adsasi
cat("So let's try with the adsasi algorithm (check the plotting device), with command ; ","\n",sep="")
cat("adsasi(simulate_one_trial,0.90)  # 0.90 for 90% power","\n",sep="")
readline(prompt="Press [enter] to continue")
proposal = adsasi(simulate_one_trial,0.90)

readline(prompt="Press [enter] to continue")
cat("The analytical sample size for the tutorial function is 44 patients. How did adsasi do ?","\n",sep="")
readline(prompt="Press [enter] to continue")

# Running adsasi again
cat("Let's do it again in case the plots were too fast or the window dimensions too small the first time","\n",sep="")
readline(prompt="Resize the plotting window and press [enter] to continue")
proposal = adsasi(simulate_one_trial,0.90)

# Running adsasi one last time
cat("Want to change one of the arguments of the function you wrote ? Put it in there : ","\n",sep="")
cat("adsasi(simulate_one_trial,0.90,effsize=.5)  # changing the effsize parameter","\n",sep="")
readline(prompt="Press [enter] to continue")
proposal = adsasi(simulate_one_trial,0.90,effsize=.5)
readline(prompt="Press [enter] to continue")

# Scalar output
cat("The adsasi function returns a scalar (sometimes Inf), so simulations can be parallelized.","\n",sep="")
cat("proposal = adsasi(simulate_one_trial,0.90)","\n",sep="")
readline(prompt="Press [enter] to continue")
proposal = adsasi(simulate_one_trial,0.90)
cat("print(proposal)","\n",sep="")
print(proposal)
readline(prompt="Press [enter] to continue")

# End
cat("Estimations above the provided maxNN argument (default 2000) are not reliable, so adapt it to your problem.","\n",sep="")
cat("There are several other failsafes that you may want to deactivate for large (and likely slow) sample size problems.","\n",sep="")
cat("Do not hesitate to display the adsasi function itself, it is heavily annotated.","\n",sep="")
cat("Happy sample size fishing !","\n",sep="")
readline(prompt="Press [enter] to finish the tutorial")

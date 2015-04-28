## This repository contains all the necessary information and code (written in R)
## to calculate the posterior predictive mass distribution for an individual planet 
## using the result of Wolfgang, Rogers, & Ford (2015).  There are many ways in which 
## the planet radii can be specified, as listed in Step 5 below.  

## This is written in R because the posterior samples output by the MCMC algorithm we've used 
## to evaluate these hierarchical Bayesian models are R objects.  Given most astronomers'
## unfamiliarity with the R language, I have written a detailed cookbook for using this
## code below under "Usage", including all of the R scripting commands you should need.  You 
## can download R at http://www.r-project.org/; the installation is straight-forward.

## Contents of repository:
## 1) README.txt
##       A copy of this document
## 2) calc_mass_postpred.R
##       R code which computes the mass posterior predictive distribution with the physical
##       density constraints (0 < M_pl < M_pureiron(R_pl))
## 3) posterior_samples.savr
##       Some of the posterior samples presented in Figure 3 of Wolfgang, Rogers, & Ford
##       (2015), selected to be under GitHub's 100 MB size limit.

## Usage:
## 1) Download the Wolfgang15 file from GitHub and unzip.
##
## 2) Open R and type in the path to that directory:
##       > filepath = "/full/path/to/files/Wolfgang15/"
##
## 3) Read the posterior samples into R:
##       > load(paste(filepath,"posterior_samples.savr",sep=""))
##
## 4) Read in the code that computes the mass posterior predictive distribution:
##       > source(paste(filepath,"calc_mass_postpred.R",sep=""))
##    You may have to download some R packages to get these functions to compile;
##    there is an R package installer in the menu bar that makes this easy.
##
## 5) Generate samples from your planet radius distribution:
##       > numsamples = 100000 #or however many you want
##
##    a) If you only have a single radius measurement with no error bars, use:
##       > planrad = (the quoted radius measurement, in units of Earth radii)
##       > planraderr = 0
##
##    b) If you have a radius measurement with an error bar, use:
##       > planrad = (The quoted planet radius, in units of Earth radii)
##       > planraderr = (The one-sigma error of that radius measurement in Earth radii)
##       The radius distribution will be assumed to be Gaussian.
##
##    c) If you have samples from the planet radius distribution, read them into a 
##       text file (i.e. "radii.txt"), label the column of planet radius values as "Rad", 
##       then in R:
##       > data = read.table("/full/path/to/text/file/radii.txt", header=T)
##       > planrad = data$Rad
##       > planraderr = 0
##
##    d) If you have samples of the planet's transit depth from light curve modeling,
##       and an estimate of the stellar radius, read the depths into a text file (i.e. 
##       "depths.txt"), label the column of depth values as "Depth", then in R:
##       > data2 = read.table("/full/path/to/text/file/depths.txt", header=T)
##       > ratio = data2$Depth  #must be pure fractions
##       > stelrad = (stellar radius measurement, in Solar radii)
##       > planrad = trdepth2rad(stelrad,ratio,numsamp=numsamples,ratioisdepth=TRUE)
##       > planraderr = 0
##       Note that if instead you get Rplanet/Rstar from your light curve modeling,
##       you set ratioisdepth=FALSE.  Everything else is the same.
##
##    e) If you have samples of the stellar radius distribution from spectroscopy,
##       and an estimate of the transit depth, read the stellar radii into a text file 
##       (i.e. "stelrad.txt"), label the column of stellar radius values as "Srad", then in R:
##       > ratio = (transit depth measurement, as a pure fraction)
##       > data3 = read.table("/full/path/to/text/file/stelrad.txt", header=T)
##       > stelrad = data3$Srad
##       > planrad = trdepth2rad(stelrad,ratio,numsamp=numsamples,ratioisdepth=TRUE)
##       > planraderr = 0
##       Note that if you have Rplanet/Rstar instead of transit depth, you set ratioisdepth=FALSE.
##       Also, if you have an estimate of the error on the transit depth/radius ratio, you can
##       set the keyword "ratioerr" to that value, which incorporates that (Gaussian) uncertainty, too:
##       > planrad = trdepth2rad(stelrad,ratio,ratioerr=0.00001,numsamp=numsamples,ratioisdepth=TRUE)
##
##    f) If you have samples of the planet's transit depth from light curve modeling
##       *and* samples from stellar radius distribution, then read these samples into
##       a text file (jointly or separately), label the columns "Depth" and "Srad", and in R: 
##       > data2 = read.table("/full/path/to/text/file/depths.txt", header=T)
##       > ratio = data2$Depth  #must be pure fractions
##       > data3 = read.table("/full/path/to/text/file/stelrad.txt", header=T)
##       > stelrad = data3$Srad
##       > planrad = trdepth2rad(stelrad,ratio,numsamp=numsamples,ratioisdepth=TRUE)
##       > planraderr = 0
##       Note that if instead you get Rplanet/Rstar from your light curve modeling,
##       you set ratioisdepth=FALSE.  Everything else is the same.
##       
## 6) Calculate the posterior predictive mass distribution:
##       > postpred = massguess_individpl(postsamp_eqn2_baseline,numsamples,planrad,raderr=planraderr)
##       > masses = postpred$masses
##
## 7) A plot showing the joint mass-radius distribution for this individual planet should
##    pop up in a graphics window upon running the massguess_individpl() function, but if
##    you'd like to save it as a postscript file, type:
##       > postscript(file=paste(filepath,"mass_postpred_Wolfgang2015.ps",sep=""),width=10,height=8)
##       > plot_individMR(postpred$radii,postpred$masses)
##       > dev.off()
##    If you don't care about the proportion of posterior predictive mass draws that correspond
##    to possibly rocky compositions, use this instead:
##       > plot_individMR(postpred$radii,postpred$masses,rockcol="black")

## Note that in Step 6 we're using the posterior samples from the default M-R relation (Eqn 2)
## with the baseline dataset (RV only, < 4 R_Earth) to calculate the posterior predictive
## distribution.  However, "posterior_samples.savr" contains a few of the other posteriors 
## presented in Wolfgang, Rogers, & Ford (2015).  Any of these samples could be used in step 6.
## The options, with self-explanatory names, are:
## 	postsamp_eqn2_baseline
## 	postsamp_eqn2_lt8
## 	postsamp_eqn2_TTVs
##

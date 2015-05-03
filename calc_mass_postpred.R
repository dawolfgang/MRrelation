## R code that calculates the posterior predictive mass distribution for an individual planet 
## using the result of Wolfgang, Rogers, & Ford (2015).  This distribution marginalizes over
## the hyperparameter posteriors displayed in Figures 2 and 3 to produce a distribution for
## the planet\'s possible masses that incorporates these uncertainties (as opposed to 
## simply using Equation 5).  Likewise, it marginalizes over the distribution of possible 
## radii that you input.  See README.txt for usage and how to execute this code.

## Copyright (C) 2015  Angie Wolfgang

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License along
## with this program; if not, write to the Free Software Foundation, Inc.,
## 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.  You can
## contact Angie Wolfgang at awolfgan@ucsc.edu as well.

rnorm.bound <- function(Ndata, mu, stddev, lower.bound = -Inf, upper.bound = 0) {
  ## Drawing from a normal distribution that is truncated by upper and lower limits.
  x = rep(0.,Ndata)
  if (length(lower.bound) == 1) {lower.bound = rep(lower.bound,Ndata)}
  if (length(upper.bound) == 1) {upper.bound = rep(upper.bound,Ndata)}
  if (length(mu) == 1) {mu = rep(mu,Ndata)}
  if (length(stddev) == 1) {stddev = rep(stddev,Ndata)}
  for (i in 1:Ndata) {
    repeat {
      x[i] = rnorm(1,mu[i],stddev[i])
      if( (x[i]>lower.bound[i]) && (x[i]<upper.bound[i]) ) {break}
    }
  }
  x
}


genmass1 <- function(norm,powlaw,rad) {
  ## This function generates masses based on the M-R relation of Eqn 1 and the input radii.
  ## Note that the returned mass array will be the same length as the input radius array, and that if
  ## the population-wide parameters are also arrays, they need to be the same length as the radius array.
  ## The normalizing constant should NOT be exponentiated yet.

  ## Using the M-R relation to calculate masses:
  logmass = powlaw * log(rad) + norm
  mass = exp(logmass)

  mass
}


genmass2 <- function(norm,powlaw,sigmam,rad) {
  ## This function generates masses based on the M-R relation of Eqn 2 and the input radii.
  ## Note that the returned mass array will be the same length as the input radius array, and that if
  ## the population-wide parameters are also arrays, they need to be the same length as the radius array.
  ## The normalizing constant should NOT be exponentiated yet.

  ## Using the M-R relation to calculate masses:

  ## The "mean" M-R relation:
  logmass = powlaw * log(rad) + norm

  ## Drawing mass probabilistically, using the parameterized scatter around this mean relation
  ## and the physical mass constraints (mass > 0 and less than that of a 100% iron planet at the  
  ## provided radius; equation is inverted from rmf quadratic fit in erratum of Fortney et al. 2007 
  ## using rmf = 0).
  MaxM = 10^((-0.4938 + sqrt(0.4938^2-4*0.0975*(0.7932-rad)))/(2*0.0975))
  mass = rnorm.bound(length(rad),exp(logmass),sigmam,lower.bound=0,upper.bound=MaxM)

  mass
}


genmass3 <- function(norm,powlaw,varm0,varslope,rad) {
  ## This function generates masses based on the M-R relation of Eqn 3 and the input radii.
  ## Note that the returned mass array will be the same length as the input radius array, and that if
  ## the population-wide parameters are also arrays, they need to be the same length as the radius array.
  ## The normalizing constant should NOT be exponentiated yet.

  ## Using the M-R relation to calculate masses:

  ## The "mean" M-R relation:
  logmass = powlaw * log(rad) + norm

  ## Calculating the non-constant scatter:
  sigmam = sqrt(varm0 + varslope*(rad-1))

  ## Drawing mass probabilistically, using the parameterized scatter around this mean relation
  ## and the physical mass constraints (mass > 0 and less than that of a 100% iron planet at the  
  ## provided radius; equation is inverted from rmf quadratic fit in erratum of Fortney et al. 2007 
  ## using rmf = 0).
  MaxM = 10^((-0.4938 + sqrt(0.4938^2-4*0.0975*(0.7932-rad)))/(2*0.0975))
  mass = rnorm.bound(length(rad),exp(logmass),sigmam,lower.bound=0,upper.bound=MaxM)

  mass
}


isplanetrocky <- function(rad, mass, likeEarth=FALSE) {
  ## Returns boolean value: true if given masses, radii in Earth units could be consistent with a rocky planet
  ## and false if not.  From erratum of Fortney et al. 2007.
  if (likeEarth) {
    radrock = (0.0592*.67 + 0.0975)*log(mass,base=10)^2 + (0.2337*.67 + 0.4938)*log(mass,base=10) + 0.3102*.67 + 0.7932
  } else {
    ## From the ice-rock mass fraction relation (in which 0% ice means 100% rock)
    #radrock = (0.1603)*log(mass,base=10)^2 + (0.7387)*log(mass,base=10)+ 1.1193
    ## From the rock-iron mass fraction relation (in which the fraction is % rock; we are looking at 100% rock)
    radrock = (0.0592*1 + 0.0975)*log(mass,base=10)^2 + (0.2337*1 + 0.4938)*log(mass,base=10) + 0.3102*1 + 0.7932
  }
  isrocky = rad < radrock

  isrocky
}


trdepth2rad <- function(stelrad, ratio, ratioerr=0, numsamp=10000, ratioisdepth=TRUE) {
  ## This function computes samples from the planet radius distribution based on a numerical
  ## stellar radius distribution (assumed to be in units of Sun radii) and a transit depth or   
  ## radius ratio (which is specified by the ratioisdepth keyword; both inputs should be pure  
  ## fractions), which could be either a single number or posterior samples itself.
  ## Note that if ratio contains posterior samples, ratioerr should remain set to the default  
  ## value of zero.  numsamp specifies the number of samples to generate.

  ## Drawing from the stellar posterior
  drawindeces = ceiling(runif(numsamp,0,length(stelrad)))
  sradsamp = stelrad[drawindeces]

  ## Drawing from the ratio posterior and converting the planet radii
  if (ratioerr == 0) {
    drawindeces2 = ceiling(runif(numsamp,0,length(ratio)))
    depthdraw = ratio[drawindeces2]
  } else {
    depthdraw = rnorm(numsamp,ratio,ratioerr)
  }
  if (ratioisdepth) {
    prad = sqrt(depthdraw)*sradsamp*110.3175       #factor is R_Sun/R_Earth, so prad is in Earth units.
  } else {
    prad = depthdraw*sradsamp*110.3175
  }

  prad
}


plot_individMR <- function(rad, mass, plotisrocky=TRUE, rockcol="red", lhist=40){
    ## This function plots the allowed masses and radii for an individual planet, with marginalized histograms
    ## on the sides, and overplots M,R combinations consistent with rocky in the color rockcol.

    ## Set up layout and graphical parameters
    layMat <- matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
    layout(layMat, widths=c(5/7, 2/7), heights=c(2/7, 5/7))
    ospc <- 0.5 # outer space
    pextd <- 4 # par extension down
    pextl <- 6 # par extension to the left
    bspc <- 0 # space between scatter plot and bar plots
    par. <- par(mar=c(pextd, pextl, bspc, bspc), oma=rep(ospc, 4)) # plot parameters

    ## Plotting the joint posterior predictive distribution
    plot(rad, mass, xlab=expression(paste("Radius (R"["Earth"],")")), ylab=expression(paste("Mass (M"["Earth"],")")), cex=.5, pch=16, cex.axis=1.6, cex.lab=1.6)
    if (plotisrocky) {
      isrocky = isplanetrocky(rad,mass)
      whererock = which(isrocky)
      points(rad[whererock],mass[whererock],cex=.5,pch=16,col=rockcol)
    }

    ## Calculate histograms (note these use probability=TRUE)
    xhist <- hist(rad, plot=FALSE, breaks=seq(from=min(rad), to=max(rad), length.out=lhist))
    yhist <- hist(mass, plot=FALSE, breaks=seq(from=min(mass), to=max(mass), length.out=lhist))
    if (plotisrocky) {
      xhistrock <- hist(rad[whererock], plot=FALSE, breaks=xhist$breaks)
      yhistrock <- hist(mass[whererock], plot=FALSE, breaks=yhist$breaks)
    }

    ## Plot marginal histograms: top first, then right
    par(mar=c(0, pextl, 0, 0))
	if (xhist$breaks[1] == xhist$breaks[2]) {
	  plot(rep(xhist$breaks[1],2),c(0,1),type="l",yaxt="n",xaxt="n",ylab="",bty="n",lwd=4)
	  line(rep(xhist$breaks[1],2),c(0,length(whererock)/length(rad)))
	} else {
      barplot(xhist$counts, axes=FALSE, space=0, col="black", border = NA)
      if (plotisrocky) {barplot(xhistrock$counts, axes=FALSE, space=0, col=rockcol, add=TRUE, border = NA)}
	}
    par(mar=c(pextd, 0, 0, 0))
    barplot(yhist$counts, axes=FALSE,space=0, horiz=TRUE, col="black", border = NA)
    if (plotisrocky) {barplot(yhistrock$counts, axes=FALSE, space=0, horiz=TRUE, col=rockcol, add=TRUE, border = NA)}

    ## Restore full-plot parameters
    par(par.)
}


massguess_individpl <- function(sim, numdraws, rad, raderr=0, likeEarth=FALSE) {
  ## This function generates a mass posterior predictive distribution based on the radius (and optional radius errors)
  ## and calculates the fraction of simulations that fall below a 100% silicate line.

  ## Figuring out which M-R relationship to call
  modelfile = sim$model.file
  wheremodelnum = regexpr(".txt",modelfile)
  modelnum = strtoi(substr(modelfile,wheremodelnum[1]-1,wheremodelnum[1]-1))
  if (!is.finite(modelnum)) {modelnum = strtoi(substr(modelfile,wheremodelnum[1]-5,wheremodelnum[1]-5))}

  ## Generating the posterior predictive distribution:
  ## First, drawing samples for the population-wide parameters that all the M-R relationships have in common.
  ## Note that the way the mass-generating functions were defined, the normalizing constant cannot be exponentiated yet.
  normconst.post = sim$BUGSoutput$sims.matrix[,"normconst"]
  powindex.post = sim$BUGSoutput$sims.matrix[,"powindex"]
  wherepostdraw = ceiling(runif(numdraws,0,length(normconst.post)))
  normdraws = normconst.post[wherepostdraw]
  powdraws = powindex.post[wherepostdraw]

  ## Second, drawing radii from the radius distribution given.
  if (raderr > 0) {
    ## If a radius error is given, then we assume a gaussian distribution for the radius.
    raddraws = rnorm.bound(numdraws, rad, raderr, lower.bound = 0.17, upper.bound = Inf)
  } else {
    ## Otherwise, we assume that rad contains a vector with posterior samples from the radius distribution
    drawindeces = ceiling(runif(numdraws,0,length(rad)))
    raddraws = rad[drawindeces]
    while (sum(raddraws < 0.17) > 0) {
      print("Warning: extremely small radii (< 0.2 R_Earth)!!  M-R relation not valid!")
      whereradtoosmall = which(raddraws < 0.17)
      raddraws[whereradtoosmall] = rad[ceiling(runif(length(whereradtoosmall),0,length(rad)))]
    }
  }

  ## Third, drawing samples for the population-wide parameters specific to the given M-R simulation,
  ## and then generating masses from the appropriate M-R relationship.
  if (modelnum == 1) {
    masses = genmass1(normdraws,powdraws,raddraws)

  } else if (modelnum == 3) {
    sigmaMR.post = sqrt(sim$BUGSoutput$sims.matrix[,"varMR"])
    sigmaMRdraws = sigmaMR.post[wherepostdraw]
    masses = genmass2(normdraws,powdraws,sigmaMRdraws,raddraws)

  } else if (modelnum == 5) {
    varMR0.post = sim$BUGSoutput$sims.matrix[,"varMR0"]
    varslope.post = sim$BUGSoutput$sims.matrix[,"varslope"]
    varMR0draws = varMR0.post[wherepostdraw]
    varslopedraws = varslope.post[wherepostdraw]
    masses = genmass3(normdraws,powdraws,varMR0draws,varslopedraws,raddraws)

  } else {
    print("Model for mass-radius relationship not recognized.  Please define a generative model for this M-R relationship.")
    masses = NA
    isrocky = NA
  }

  ## Calculating which samples yield possibly rocky planets, and plotting a joint posterior predictive
  ## distribution with marginalized mass, radius distributions.
  if (is.numeric(masses)) {
    plot_individMR(raddraws,masses)
    isrocky = isplanetrocky(raddraws,masses,likeEarth=likeEarth)
    print(paste("Fraction of (mass,radius) samples that are more dense than 100% silicate rock:",sum(isrocky)/length(masses)))
  }
 
  ## Returning the joint posterior predictive distribution. 
  list(radii=raddraws,masses=masses,isrocky=isrocky)
}

\name{marginalized.risk}
\alias{marginalized.risk}
\alias{marginalizedRisk}
\alias{marginalized.risk.cat}
\alias{marginalized.risk.cont}
\title{
Compute Maringalized Risk
}
\description{
Computes risk of disease as a function of marker s by marginalizedizing over a covariate vector Z.
}
\usage{

marginalized.risk(fit.risk, marker.name, data, categorical.s, 
    weights = rep(1,nrow(data)), t = NULL, ss = NULL, verbose = FALSE)

marginalized.risk.cont(fit.risk, marker.name, data, 
    weights = rep(1, nrow(data)), t=NULL, ss = NULL, verbose = FALSE)

marginalized.risk.cat(fit.risk, marker.name, data, 
    weights = rep(1, nrow(data)), t =NULL, verbose = FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{fit.risk}{
A regression object where the outcome is risk of disease, e.g. y~Z+marker. Need to support predict(fit.risk)
}
  \item{marker.name}{
string
}
  \item{data}{
A data frame containing the phase 2 data
}
  \item{ss}{
A vector of marker values
}
  \item{weights}{
Inverse prob sampling weight, optional
}
  \item{t}{
If fit.risk is Cox regression, t is the time at which distribution function will be assessed
}
  \item{categorical.s}{
TRUE if the marker is categorical, FALSE otherwise
}
  \item{verbose}{
Boolean
}
}
\details{
See the vignette file for more details.
}
\value{
If ss is not NULL, a vector of probabilities are returned. If ss is NULL, a matrix of two columns are returned, where the first column is the marker value and the second column is the probabilties.
}
\examples{


#### suppose wt.loss is the marker of interest

if(requireNamespace("survival")) {

library(survival)

dat=subset(lung, !is.na(wt.loss) & !is.na(ph.ecog))

f1=Surv(time, status) ~ wt.loss + ph.ecog + age + sex

fit.risk = coxph(f1, data=dat)
     
ss=quantile(dat$wt.loss, seq(.05,.95,by=0.01))
t0=1000
prob = marginalized.risk(fit.risk, "wt.loss", dat, categorical.s=FALSE, t = t0, ss=ss)

plot(ss, prob, type="l", xlab="Weight loss", ylab=paste0("Probability of survival at day ", t0))

}


\dontrun{

#### Efron bootstrap to get confidence band

# store the current rng state 
save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }      

B=10 # bootstrap replicates, 1000 is good
numCores=1 # multiple cores can speed things up
library(doParallel)
out=mclapply(1:B, mc.cores = numCores, FUN=function(seed) {   
    set.seed(seed)    
    # a simple resampling scheme here. needs to be adapted to the sampling scheme
    dat.tmp=dat[sample(row(dat), replace=TRUE),]        
    fit.risk = coxph(f1, data=dat)
    marginalized.risk(fit.risk, "wt.loss", dat.tmp, categorical.s=FALSE, t = t0, ss=ss)
})
res=do.call(cbind, out)

# restore rng state 
assign(".Random.seed", save.seed, .GlobalEnv)    

# quantile bootstrap CI
ci.band=t(apply(res, 1, function(x) quantile(x, c(.025,.975))))

plot(ss, prob, type="l", xlab="Weight loss", ylab=paste0("Probability of survival at day ", t0), 
    ylim=range(ci.band))
lines(ss, ci.band[,1], lty=2)
lines(ss, ci.band[,2], lty=2)

}

}

\name{marginalized.risk.threshold}
\alias{marginalized.risk.threshold}
\title{
Compute Maringalized Risk as a Function of S>=s
}
\description{
Computes risk of disease conditional on S>=s by marginalizedizing over a covariate vector Z.
}
\usage{

marginalized.risk.threshold(formula, marker.name, data, weights=rep(1, nrow(data)), 
    t, ss=NULL, verbose=FALSE) 


}

\arguments{

  \item{formula}{
A formula for coxph 
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
t is the time at which survival will be assessed
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
f1=Surv(time, status) ~ ph.ecog + age + sex
ss=quantile(dat$wt.loss, seq(.05,.95,by=0.01))
t0=1000
prob = marginalized.risk.threshold(f1, "wt.loss", dat, t = t0, ss=ss)

plot(ss, prob, type="l", xlab="Weight loss (S>=s)", 
    ylab=paste0("Probability of survival at day ", t0))

}


}

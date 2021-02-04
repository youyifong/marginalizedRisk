test.marginalized.risk <- function() {

library("RUnit")
library("marginalizedRisk")
tolerance=1e-5
RNGkind("Mersenne-Twister", "Inversion")


library(survival)
dat=subset(lung, !is.na(wt.loss) & !is.na(ph.ecog))

f1=Surv(time, status) ~ wt.loss + ph.ecog + age + sex
fit.risk = coxph(f1, data=dat)     
ss=quantile(dat$wt.loss, seq(.05,.95,by=0.1))
t0=1000
prob = marginalized.risk(fit.risk, "wt.loss", dat, categorical.s=FALSE, t = t0, ss=ss)
checkEqualsNumeric (prob, c(0.9298895,0.9232542,0.9232542,0.9204857,0.9162091,0.9117831,0.9072068,0.9008702,0.8873893,0.8682598), tolerance=tolerance)


f1=Surv(time, status) ~ ph.ecog + age + sex
ss=quantile(dat$wt.loss, seq(.05,.95,by=0.05))
t0=1000
prob = marginalized.risk.threshold(f1, "wt.loss", dat, t = t0, ss=ss)
checkEqualsNumeric (prob, c(0.9069922,0.9066352,0.9092511,0.9092511,0.9092511,0.9057564,0.9113864,0.9163123,0.9249508,0.9143111,0.9192401,0.9220423,0.8951074,0.8879142,0.9174229,0.9220667,0.9327394,0.9271225,0.6151954), tolerance=tolerance)



}

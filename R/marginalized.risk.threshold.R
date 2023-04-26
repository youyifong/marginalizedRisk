# only pass ph2 data to this function
marginalized.risk.threshold=function(formula, marker.name, data, weights=rep(1, nrow(data)), t, ss=NULL, verbose=FALSE) {
    
    time.var=terms(formula)[[2]][[2]]
    outcome.var=terms(formula)[[2]][[3]]
    
    # make weights part of the data frame under a unique name
    data$wt999=weights
    
    wt999=NULL # this is to trick R CMD check. otherwise, it will complain about coxph call: Undefined global functions or variables coxph
    
    ss.is.null=is.null(ss) 
    if (ss.is.null) ss=quantile(data[[marker.name]], seq(0,.9,by=0.05))
        
    risks=sapply (ss, function(s) {    
        myprint(s)
        tmp=data[data[[marker.name]]>=s, ]
        fit.risk.1=try(survival::coxph(formula, tmp, weights=wt999, model=TRUE))
        if (inherits(fit.risk.1, "try-error")) return (NA)
        
        tmp=data # g-computation formula sums over all ph2 data instead of [data[[marker.name]]>=s,]
        tmp[[time.var]]=t
        # need to wrap predict in try because when there are no cases, even model=T won't prevent errors: 
        #           Error in predict.coxph(fit.risk.1, type = "expected") : 
        #               Data is not the same size as it was in the original fit
        pred=try(predict(fit.risk.1, newdata=tmp, type="expected"), silent=T)
        myprint(pred)
        #
        if (class(pred) != "try-error" & all(!is.na(pred))) {
            weighted.mean(1 - exp(-pred), tmp$wt999)
        } else {
            # the error we have seen is due to no cases, 0 is a reasonable output
            # but a different error may occur
            # either way we return a simple estimate of Pr(Y=1|S>=s)
            weighted.mean(tmp[[outcome.var]], tmp$wt999)
        }        
    })
    
    if (ss.is.null) cbind(marker=ss, prob=risks) else risks

}

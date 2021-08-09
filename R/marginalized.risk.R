# only pass ph2 data to these functions
marginalized.risk=function(fit.risk, marker.name, data, categorical.s, weights=rep(1, nrow(data)), t=NULL, ss=NULL, verbose=FALSE) {
    if(categorical.s) {
        marginalized.risk.cat  (fit.risk, marker.name, data, weights=weights, t=t, verbose=verbose) 
    } else {
        marginalized.risk.cont (fit.risk, marker.name, data, weights=weights, t=t, ss=ss, verbose=verbose) 
    }
}

# categorical markers
marginalized.risk.cat=function(fit.risk, marker.name, data, weights=rep(1, nrow(data)), t=NULL, verbose=FALSE) {  
    
    if("coxph" %in% class(fit.risk)) {
        time.var=as.character(fit.risk$terms[[2]][[2]])
        y.var=as.character(fit.risk$terms[[2]][[3]])
    }
    
    # ss gives the levels that S takes
    ss=unique(data[[marker.name]]); ss=sort(ss[!is.na(ss)])
    
    if (!"coxph" %in% class(fit.risk)) {
        # logistic regression
        dat.tmp.mrc=data
        risks=sapply(ss, function(s) {
            dat.tmp.mrc[[marker.name]]=s    
            risks = predict(fit.risk, newdata=dat.tmp.mrc, type="response") # glm
            sum(weights * risks) / sum(weights)    
        })
        names(risks)=levels(ss)
        risks        
            
    } else {
        # coxph or svycoxph
        if (is.null(t)) {
            # return risk versus time
            tt=sort(unique(data[[time.var]][data[[y.var]]==1]))        
            risks=sapply(tt, function (t) {
                dat.tmp.mrc=data
                dat.tmp.mrc[[time.var]]=t
                risks=sapply(ss, function(s) {        
                    dat.tmp.mrc[[marker.name]]=s    
                    risks = 1 - exp(-predict(fit.risk, newdata=dat.tmp.mrc, type="expected"))# coxph survival prob
                    sum(weights * risks) / sum(weights)
                })
            })
            risks=t(risks)
            colnames(risks)=as.character(ss)        
            list(time=tt, risk=risks)
            
        } else {
            if (verbose) print("return risk at time t")
            dat.tmp.mrc=data
            time.var=fit.risk$terms[[2]][[2]]
            dat.tmp.mrc[[time.var]]=t        
            risks=sapply(ss, function(s) {
                dat.tmp.mrc[[marker.name]]=s    
                risks = 1 - exp(-predict(fit.risk, newdata=dat.tmp.mrc, type="expected")) # coxph survival prob
                sum(weights * risks) / sum(weights)    
            })
            names(risks)=levels(ss)
            risks        
        }
    }
}

# continuous markers
marginalized.risk.cont=function(fit.risk, marker.name, data, weights=rep(1, nrow(data)), t=NULL, ss=NULL, verbose=FALSE) {
    ss.is.null=is.null(ss) 
    if (ss.is.null) ss=quantile(data[[marker.name]], seq(.05,.95,by=0.01))
        
    dat.tmp.mri=data
    if (!is.null(t)) {
        time.var=fit.risk$terms[[2]][[2]]
        dat.tmp.mri[[time.var]]=t
    }
    
    risks=sapply(ss, function(s) {
        dat.tmp.mri[[marker.name]]=s    
        risks = if(is.null(t)) {
            # glm
            predict(fit.risk, newdata=dat.tmp.mri, type="response")
        } else {
            # coxph survival prob
            1 - exp(-predict(fit.risk, newdata=dat.tmp.mri, type="expected"))
        }
        #if(any(is.na(risks))) stop("NA's found in fit.risk")
        
        sum(weights * risks) / sum(weights)    
    })
    
    if (ss.is.null) cbind(marker=ss, prob=risks) else risks
}

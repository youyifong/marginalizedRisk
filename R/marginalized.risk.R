1# only pass ph2 data to these functions
marginalized.risk=function(fit.risk, marker.name, data, categorical.s, weights=rep(1, nrow(data)), t=NULL, ss=NULL, verbose=FALSE, t.end=NULL) {
    if(categorical.s) {
        marginalized.risk.cat  (fit.risk, marker.name, data, weights=weights, t=t, verbose=verbose, t.end=t.end) 
    } else {
        marginalized.risk.cont (fit.risk, marker.name, data, weights=weights, t=t, ss=ss, verbose=verbose) 
    }
}

# categorical markers
marginalized.risk.cat=function(fit.risk, marker.name, data, weights=rep(1, nrow(data)), t=NULL, verbose=FALSE, t.end=NULL) {  
    
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
            
            # system.time({
            # tt=sort(unique(data[[time.var]][data[[y.var]]==1]))        
            # if (!is.null(t.end)) tt=unique(c(tt, t.end))
            # risks=sapply(tt, function (t) {
            #     dat.tmp.mrc=data
            #     dat.tmp.mrc[[time.var]]=t
            #     risks=sapply(ss, function(s) {        
            #         dat.tmp.mrc[[marker.name]]=s    
            #         
            #         # from survey 4.1-1 to 4.2-1, predict.svycoxph changed from
            #         # type=c("lp", "risk", "expected", "terms","curve"),
            #         # to 
            #         # type=c("lp", "risk", "terms","curve"),
            #         
            #         # one fix is to change predict to survival:::predict.coxph, which gives the same results
            #         risks = 1 - exp(-survival:::predict.coxph(fit.risk, newdata=dat.tmp.mrc, type="expected"))# coxph survival prob
            #         
            #         sum(weights * risks) / sum(weights)
            #     })
            # })
            # risks=t(risks)
            # colnames(risks)=as.character(ss)        
            # resa=list(time=tt, risk=risks)
            # })
            
            # system.time({
            # alternatively, the curve type can be used to speed up 10x
            # resb differs from resa in the support of time. in resb, time includes both censoring and events, but not t.end
            risks=sapply(ss, simplify="array", function(s) {        
              dat.tmp.mrc=data
              dat.tmp.mrc[[marker.name]]=s    
              
              risks = predict(fit.risk, newdata=dat.tmp.mrc, type="curve")
              
              tt = sapply(risks, function (x) x$time)[,1]
              risks = mysapply(risks, function (x) x$surv)
              
              cbind(tt, colSums(weights * risks) / sum(weights))
            })
            resb=list(time=risks[,"tt",1], risk=risks[,2,])
            # })
            
            # par(mfrow=c(1,2))
            # mymatplot(resa$time, resa$risk)
            # mymatplot(resb$time, 1-resb$risk, add=T)
            
            resb
            
        } else {
            if (verbose) print("return risk at time t")
            dat.tmp.mrc=data
            time.var=fit.risk$terms[[2]][[2]]
            dat.tmp.mrc[[time.var]]=t        
            risks=sapply(ss, function(s) {
                dat.tmp.mrc[[marker.name]]=s    
                
                risks = 1 - exp(-survival:::predict.coxph(fit.risk, newdata=dat.tmp.mrc, type="expected"))# coxph survival prob

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
          # survival prob
          # 1 - exp(-predict(fit.risk, newdata=dat.tmp.mri, type="expected"))
          # in survey 4.4, the above is not allowed, instead, we need
          # 1. Generate the survival curves for the new data
          sf <- survfit(fit.risk, newdata = dat.tmp.mri)
          # 2. Extract the survival probability for each person at their specific time
          predicted_surv <- diag(sf$surv[match(dat.tmp.mrc[[time.var]], sf[[time.var]]), ])
          1 - predicted_surv
          
        }
        #if(any(is.na(risks))) stop("NA's found in fit.risk")
        
        sum(weights * risks) / sum(weights)    
    })
    
    if (ss.is.null) cbind(marker=ss, prob=risks) else risks
}



marginalized.risk.cont.2=function(fit.risk, marker.name, data, weights=rep(1, nrow(data)), t, ss, marker.name.2, s.2, verbose=FALSE) {
        
    is.coxph=FALSE
    if (length(fit.risk$terms[[2]])==3) {
        # presume to be coxph
        is.coxph=TRUE
        time.var=fit.risk$terms[[2]][[2]]
        data[[time.var]]=t
    } 
    
    data[[marker.name.2]]=s.2
    
    risks=sapply(ss, function(s) {
        data[[marker.name]]=s    
        risks = if(is.coxph) {
          # 1 - exp(-predict(fit.risk, newdata=data, type="expected"))
          # in survey 4.4, the above is not allowed, instead, we need
          # 1. Generate the survival curves for the new data
          sf <- survfit(fit.risk, newdata = data)
          # 2. Extract the survival probability for each person at their specific time
          predicted_surv <- diag(sf$surv[match(dat.tmp.mrc[[time.var]], sf[[time.var]]), ])
          1 - predicted_surv
          
        } else {
            # glm
            predict(fit.risk, newdata=data, type="response")
        }
        #if(any(is.na(risks))) stop("NA's found in fit.risk")
        
        sum(weights * risks) / sum(weights)    
    })
    
    risks
}

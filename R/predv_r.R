predv_r <- function(x,prop_min,prop_max,zb=TRUE,n_iter=100000,...){
  object <- reitsma(x)
  m1 <- summary(reitsma(x))
  
  if(zb){
    if(missing(prop_min)&missing(prop_max)){
      
      PredRange <-function(object, ...) UseMethod("PredRange")
      
      PredRange.default <- function(object, mu,Sigma,alphasens = 1, alphafpr = 1,
                                    n.iter = n_iter, FUN, ...){
        samples <- rmvnorm(n.iter, mu, Sigma)
        sens <- inv.trafo(alphasens,samples[,1])
        fpr <- inv.trafo(alphafpr,samples[,2])
        out <- lapply(FUN, function(x){x(sens, fpr)})
        class(out) <- "PredRange"
        out
      }
      
      PredRange.reitsma <- function(object, prev, n.iter = n_iter, FUN = NULL, ...){
        fit <- object
        if(length(coef(fit)) > 2){
          stop("PredRange is not be used for meta-regression!")}
        if(is.null(FUN)){FUN <- list(npv= function(sens,fpr){((1 - fpr)*(1-prev))/(((1 - fpr)*(1-prev)) + (prev*(1-sens)))},
                                     ppv = function(sens, fpr){(sens*prev)/((sens*prev)+(fpr*(1-prev)))})}
        PredRange.default(mu = coef(fit)[1:2], Sigma = vcov(fit), 
                          alphasens = fit$alphasens, 
                          alphafpr = fit$alphafpr,
                          n.iter = n.iter, FUN = FUN)
      }
      
      x$prev <- (x$TP + x$FN)/(x$TP + x$FN + x$FP + x$TN)
      min_prev <- min(x$prev)
      max_prev <- max(x$prev)
      
      vec_prev <- c(round(seq(min_prev,max_prev,by=0.01),2))
      
      results_pred <- vector("list",length(vec_prev))
      for(i in seq_along(vec_prev)){
        results_pred[[i]] <- PredRange.reitsma(object,prev=vec_prev[[i]])
      }
      
      results_npv <- matrix(nrow=length(vec_prev),ncol=12)
      colnames(results_npv) <- c("prevalence","mean","sd","p2.5","p5","p10","p25","p50","p75","p90","p95","p97.5")
      results_npv[,1] <- vec_prev
      for(i in seq_along(results_pred)){
        results_npv[[i,2]] <- round(mean(results_pred[[i]]$npv),3)
        results_npv[[i,3]] <- round(sd(results_pred[[i]]$npv),3)
        results_npv[[i,4]] <- round(quantile(results_pred[[i]]$npv,p=0.025),3)
        results_npv[[i,5]] <- round(quantile(results_pred[[i]]$npv,p=0.05),3)
        results_npv[[i,6]] <- round(quantile(results_pred[[i]]$npv,p=0.10),3)
        results_npv[[i,7]] <- round(quantile(results_pred[[i]]$npv,p=0.25),3)
        results_npv[[i,8]] <- round(quantile(results_pred[[i]]$npv,p=0.5),3)
        results_npv[[i,9]] <- round(quantile(results_pred[[i]]$npv,p=0.75),3)
        results_npv[[i,10]] <- round(quantile(results_pred[[i]]$npv,p=0.90),3)
        results_npv[[i,11]] <- round(quantile(results_pred[[i]]$npv,p=0.95),3)
        results_npv[[i,12]] <- round(quantile(results_pred[[i]]$npv,p=0.975),3)
      }
      
      results_ppv <- matrix(nrow=length(vec_prev),ncol=12)
      colnames(results_ppv) <- c("prevalence","mean","sd","p2.5","p5","p10","p25","p50","p75","p90","p95","p97.5")
      results_ppv[,1] <- vec_prev
      for(i in seq_along(results_pred)){
        results_ppv[[i,2]] <- round(mean(results_pred[[i]]$ppv),3)
        results_ppv[[i,3]] <- round(sd(results_pred[[i]]$ppv),3)
        results_ppv[[i,4]] <- round(quantile(results_pred[[i]]$ppv,p=0.025),3)
        results_ppv[[i,5]] <- round(quantile(results_pred[[i]]$ppv,p=0.05),3)
        results_ppv[[i,6]] <- round(quantile(results_pred[[i]]$ppv,p=0.10),3)
        results_ppv[[i,7]] <- round(quantile(results_pred[[i]]$ppv,p=0.25),3)
        results_ppv[[i,8]] <- round(quantile(results_pred[[i]]$ppv,p=0.5),3)
        results_ppv[[i,9]] <- round(quantile(results_pred[[i]]$ppv,p=0.75),3)
        results_ppv[[i,10]] <- round(quantile(results_pred[[i]]$ppv,p=0.90),3)
        results_ppv[[i,11]] <- round(quantile(results_pred[[i]]$ppv,p=0.95),3)
        results_ppv[[i,12]] <- round(quantile(results_pred[[i]]$ppv,p=0.975),3)
      }
      
      results_npv <- as.data.frame(results_npv)
      results_ppv <- as.data.frame(results_ppv)
      
      warning("Please note that included primary studies may not be the best ones to ascertain the prevalence of a disease/condition, particularly if some of such studies have a case-control design.")
      
    }else{
      
      PredRange <-function(object, ...) UseMethod("PredRange")
      
      PredRange.default <- function(object, mu,Sigma,alphasens = 1, alphafpr = 1,
                                    n.iter = n_iter, FUN, ...){
        samples <- rmvnorm(n.iter, mu, Sigma)
        sens <- inv.trafo(alphasens,samples[,1])
        fpr <- inv.trafo(alphafpr,samples[,2])
        out <- lapply(FUN, function(x){x(sens, fpr)})
        class(out) <- "PredRange"
        out
      }
      
      PredRange.reitsma <- function(object, prev, n.iter = n_iter, FUN = NULL, ...){
        fit <- object
        if(length(coef(fit)) > 2){
          stop("PredRange is not be used for meta-regression!")}
        if(is.null(FUN)){FUN <- list(npv= function(sens,fpr){((1 - fpr)*(1-prev))/(((1 - fpr)*(1-prev)) + (prev*(1-sens)))},
                                     ppv = function(sens, fpr){(sens*prev)/((sens*prev)+(fpr*(1-prev)))})}
        PredRange.default(mu = coef(fit)[1:2], Sigma = vcov(fit), 
                          alphasens = fit$alphasens, 
                          alphafpr = fit$alphafpr,
                          n.iter = n.iter, FUN = FUN)
      }
      
      vec_prev <- c(seq(prop_min,prop_max,by=0.01))
      
      results_pred <- vector("list",length(vec_prev))
      for(i in seq_along(vec_prev)){
        results_pred[[i]] <- PredRange.reitsma(object,prev=vec_prev[[i]])
      }
      
      results_npv <- matrix(nrow=length(vec_prev),ncol=12)
      colnames(results_npv) <- c("prevalence","mean","sd","p2.5","p5","p10","p25","p50","p75","p90","p95","p97.5")
      results_npv[,1] <- vec_prev
      for(i in seq_along(results_pred)){
        results_npv[[i,2]] <- round(mean(results_pred[[i]]$npv),3)
        results_npv[[i,3]] <- round(sd(results_pred[[i]]$npv),3)
        results_npv[[i,4]] <- round(quantile(results_pred[[i]]$npv,p=0.025),3)
        results_npv[[i,5]] <- round(quantile(results_pred[[i]]$npv,p=0.05),3)
        results_npv[[i,6]] <- round(quantile(results_pred[[i]]$npv,p=0.10),3)
        results_npv[[i,7]] <- round(quantile(results_pred[[i]]$npv,p=0.25),3)
        results_npv[[i,8]] <- round(quantile(results_pred[[i]]$npv,p=0.5),3)
        results_npv[[i,9]] <- round(quantile(results_pred[[i]]$npv,p=0.75),3)
        results_npv[[i,10]] <- round(quantile(results_pred[[i]]$npv,p=0.90),3)
        results_npv[[i,11]] <- round(quantile(results_pred[[i]]$npv,p=0.95),3)
        results_npv[[i,12]] <- round(quantile(results_pred[[i]]$npv,p=0.975),3)
      }
      
      results_ppv <- matrix(nrow=length(vec_prev),ncol=12)
      colnames(results_ppv) <- c("prevalence","mean","sd","p2.5","p5","p10","p25","p50","p75","p90","p95","p97.5")
      results_ppv[,1] <- vec_prev
      for(i in seq_along(results_pred)){
        results_ppv[[i,2]] <- round(mean(results_pred[[i]]$ppv),3)
        results_ppv[[i,3]] <- round(sd(results_pred[[i]]$ppv),3)
        results_ppv[[i,4]] <- round(quantile(results_pred[[i]]$ppv,p=0.025),3)
        results_ppv[[i,5]] <- round(quantile(results_pred[[i]]$ppv,p=0.05),3)
        results_ppv[[i,6]] <- round(quantile(results_pred[[i]]$ppv,p=0.10),3)
        results_ppv[[i,7]] <- round(quantile(results_pred[[i]]$ppv,p=0.25),3)
        results_ppv[[i,8]] <- round(quantile(results_pred[[i]]$ppv,p=0.5),3)
        results_ppv[[i,9]] <- round(quantile(results_pred[[i]]$ppv,p=0.75),3)
        results_ppv[[i,10]] <- round(quantile(results_pred[[i]]$ppv,p=0.90),3)
        results_ppv[[i,11]] <- round(quantile(results_pred[[i]]$ppv,p=0.95),3)
        results_ppv[[i,12]] <- round(quantile(results_pred[[i]]$ppv,p=0.975),3)
      }
      
      results_npv <- as.data.frame(results_npv)
      results_ppv <- as.data.frame(results_ppv)
    }
  }else{
    if(missing(prop_min)&missing(prop_max)){
      
      x$prev <- (x$TP + x$FN)/(x$TP + x$FN + x$FP + x$TN)
      min_prev <- min(x$prev)
      max_prev <- max(x$prev)
     
      df1 <- data.frame("sens_lb"=1/(1+exp(-(m1$coefficients[1,5]))),"sens_ub"=1/(1+exp(-(m1$coefficients[1,6]))),"spec_lb"=1-(1/(1+exp(-(m1$coefficients[2,6])))),"spec_ub"=1-(1/(1+exp(-(m1$coefficients[2,5])))))
      sens <- beta.parms.from.quantiles(q=c(df1$sens_lb,df1$sens_ub)) 
      spec <- beta.parms.from.quantiles(q=c(df1$spec_lb,df1$spec_ub))
      
      df1 <- data.frame(df1,"sens_a"=sens$a,"sens_b"=sens$b,"spec_a"=spec$a,"spec_b"=spec$b)
      
      
      df3 <- data.frame(rsens=rbeta(n_iter,df1$sens_a,df1$sens_b),rspec=rbeta(n_iter,df1$spec_a,df1$spec_b),rprev=round(runif(n_iter,min_prev,max_prev),digits=2))
      df3$ppv <- ((df3$rsens*df3$rprev)/((df3$rsens*df3$rprev)+((1-df3$rspec)*(1-df3$rprev))))
      df3$npv <- ((df3$rspec*(1-df3$rprev))/((df3$rspec*(1-df3$rprev))+(df3$rprev*(1-df3$rsens))))
      
      quantile_more <- function(x){
        quantile(x,probs=c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975))
      }
      
      
      df4a <- aggregate(df3$npv, list(df3$rprev), FUN=mean)
      df4b <- aggregate(df3$npv, list(df3$rprev), FUN=sd)
      df4c <- aggregate(df3$npv, list(df3$rprev), FUN=quantile_more)
      
      df5a <- aggregate(df3$ppv, list(df3$rprev), FUN=mean)
      df5b <- aggregate(df3$ppv, list(df3$rprev), FUN=sd)
      df5c <- aggregate(df3$ppv, list(df3$rprev), FUN=quantile_more)
      
      results_npv <- data.frame("prevalence"=df4a[,1],"mean"=round(df4a[,2],3),"sd"=round(df4b[,2],3),"p2.5"=round(df4c$x[,"2.5%"],3),"p5"=round(df4c$x[,"5%"],3),"p10"=round(df4c$x[,"10%"],3),"p25"=round(df4c$x[,"25%"],3),"p50"=round(df4c$x[,"50%"],3),"p75"=round(df4c$x[,"75%"],3),"p90"=round(df4c$x[,"90%"],3),"p95"=round(df4c$x[,"95%"],3),"p97.5"=round(df4c$x[,"97.5%"],3))
      results_ppv <- data.frame("prevalence"=df5a[,1],"mean"=round(df5a[,2],3),"sd"=round(df5b[,2],3),"p2.5"=round(df5c$x[,"2.5%"],3),"p5"=round(df5c$x[,"5%"],3),"p10"=round(df5c$x[,"10%"],3),"p25"=round(df5c$x[,"25%"],3),"p50"=round(df5c$x[,"50%"],3),"p75"=round(df5c$x[,"75%"],3),"p90"=round(df5c$x[,"90%"],3),"p95"=round(df5c$x[,"95%"],3),"p97.5"=round(df5c$x[,"97.5%"],3))
      
      warning("Please note that included primary studies may not be the best ones to ascertain the prevalence of a disease/condition, particularly if some of such studies have a case-control design.")
      
    }else{
      
      df1 <- data.frame("sens_lb"=1/(1+exp(-(m1$coefficients[1,5]))),"sens_ub"=1/(1+exp(-(m1$coefficients[1,6]))),"spec_lb"=1-(1/(1+exp(-(m1$coefficients[2,6])))),"spec_ub"=1-(1/(1+exp(-(m1$coefficients[2,5])))))
      sens <- beta.parms.from.quantiles(q=c(df1$sens_lb,df1$sens_ub))
      spec <- beta.parms.from.quantiles(q=c(df1$spec_lb,df1$spec_ub))
      
      df1 <- data.frame(df1,"sens_a"=sens$a,"sens_b"=sens$b,"spec_a"=spec$a,"spec_b"=spec$b)
      
      
      df3 <- data.frame(rsens=rbeta(n_iter,df1$sens_a,df1$sens_b),rspec=rbeta(n_iter,df1$spec_a,df1$spec_b),rprev=round(runif(n_iter,prop_min,prop_max),digits=2))
      df3$ppv <- ((df3$rsens*df3$rprev)/((df3$rsens*df3$rprev)+((1-df3$rspec)*(1-df3$rprev))))
      df3$npv <- ((df3$rspec*(1-df3$rprev))/((df3$rspec*(1-df3$rprev))+(df3$rprev*(1-df3$rsens))))
      
      quantile_more <- function(x){
        quantile(x,probs=c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975))
      }
      
      
      df4a <- aggregate(df3$npv, list(df3$rprev), FUN=mean)
      df4b <- aggregate(df3$npv, list(df3$rprev), FUN=sd)
      df4c <- aggregate(df3$npv, list(df3$rprev), FUN=quantile_more)
      
      df5a <- aggregate(df3$ppv, list(df3$rprev), FUN=mean)
      df5b <- aggregate(df3$ppv, list(df3$rprev), FUN=sd)
      df5c <- aggregate(df3$ppv, list(df3$rprev), FUN=quantile_more)
      
      results_npv <- data.frame("prevalence"=df4a[,1],"mean"=round(df4a[,2],3),"sd"=round(df4b[,2],3),"p2.5"=round(df4c$x[,"2.5%"],3),"p5"=round(df4c$x[,"5%"],3),"p10"=round(df4c$x[,"10%"],3),"p25"=round(df4c$x[,"25%"],3),"p50"=round(df4c$x[,"50%"],3),"p75"=round(df4c$x[,"75%"],3),"p90"=round(df4c$x[,"90%"],3),"p95"=round(df4c$x[,"95%"],3),"p97.5"=round(df4c$x[,"97.5%"],3))
      results_ppv <- data.frame("prevalence"=df5a[,1],"mean"=round(df5a[,2],3),"sd"=round(df5b[,2],3),"p2.5"=round(df5c$x[,"2.5%"],3),"p5"=round(df5c$x[,"5%"],3),"p10"=round(df5c$x[,"10%"],3),"p25"=round(df5c$x[,"25%"],3),"p50"=round(df5c$x[,"50%"],3),"p75"=round(df5c$x[,"75%"],3),"p90"=round(df5c$x[,"90%"],3),"p95"=round(df5c$x[,"95%"],3),"p97.5"=round(df5c$x[,"97.5%"],3))
    }
  }
  output <- list(NPV = results_npv, PPV = results_ppv, call = match.call())
  class(output) <- "predv_r"
  output  
}

print.predv_r <- function(x,ylim_npv=c(0,1),ylim_ppv=c(0,1),...){
  cat("Estimates of predictive values")
  cat("\n")
  cat("Minimum prevalence:")
  print(min(x$NPV$prevalence))
  cat("\n")
  cat("Maximum prevalence:")
  print(max(x$NPV$prevalence))
  cat("\n")
  cat("NPV")
  cat("\n")
  print(x$NPV[,1:3])
  cat("\n")
  cat("PPV")
  cat("\n")
  print(x$PPV[,1:3])
  
  plot(x=x$NPV$prevalence,y=x$NPV$mean,type="l",ylab="Negative predictive value",ylim=ylim_npv,main="Predicted negative predictive values",xlab="Prevalence")
  lines(x$NPV$prevalence,x$NPV$p97.5, type = "l", lty=2)
  lines(x$NPV$prevalence,x$NPV$p2.5, type = "l", lty=2)
  points(x$NPV$prevalence,x$NPV$mean,pch=16)
  
  plot(x=x$PPV$prevalence,y=x$PPV$mean,type="l",ylab="Positive predictive value",ylim=ylim_ppv,main="Predicted positive predictive values",xlab="Prevalence")
  lines(x$PPV$prevalence,x$PPV$p97.5, type = "l", lty=2)
  lines(x$PPV$prevalence,x$PPV$p2.5, type = "l", lty=2)
  points(x$PPV$prevalence,x$PPV$mean,pch=16)
  
  plot(x=x$NPV$prevalence,y=x$NPV$mean,type="l",ylab="Predictive value",col=4,ylim=c(0,1),main="Projected predictive values",xlab="Prevalence")
  lines(x$NPV$prevalence,x$NPV$p97.5, type = "l", lty=2,col=4)
  lines(x$NPV$prevalence,x$NPV$p2.5, type = "l", lty=2,col=4)
  points(x$NPV$prevalence,x$NPV$mean,pch=16,col=4)
  lines(x=x$PPV$prevalence,y=x$PPV$mean,type="l",col=2)
  lines(x$PPV$prevalence,x$PPV$p97.5, type = "l", lty=2,col=2)
  lines(x$PPV$prevalence,x$PPV$p2.5, type = "l", lty=2,col=2)
  points(x$PPV$prevalence,x$PPV$mean,pch=16,col=2)
  legend("bottomright",legend = c("NPV", "PPV"),col = c(4,2), pch = 16, bty="n")
}

summary.predv_r <- function(object,ylim_npv=c(0,1),ylim_ppv=c(0,1),...){
  x <- object
  cat("Estimates of predictive values")
  cat("\n")
  cat("Minimum prevalence:")
  print(min(x$NPV$prevalence))
  cat("\n")
  cat("Maximum prevalence:")
  print(max(x$NPV$prevalence))
  cat("\n")
  cat("NPV")
  cat("\n")
  print(x$NPV)
  cat("\n")
  cat("PPV")
  cat("\n")
  print(x$PPV)
  
  
  plot(x=x$NPV$prevalence,y=x$NPV$mean,type="l",ylab="Negative predictive value",ylim=ylim_npv,main="Predicted negative predictive values",xlab="Prevalence")
  lines(x$NPV$prevalence,x$NPV$p97.5, type = "l", lty=2)
  lines(x$NPV$prevalence,x$NPV$p2.5, type = "l", lty=2)
  points(x$NPV$prevalence,x$NPV$mean,pch=16)
  
  plot(x=x$PPV$prevalence,y=x$PPV$mean,type="l",ylab="Positive predictive value",ylim=ylim_ppv,main="Predicted positive predictive values",xlab="Prevalence")
  lines(x$PPV$prevalence,x$PPV$p97.5, type = "l", lty=2)
  lines(x$PPV$prevalence,x$PPV$p2.5, type = "l", lty=2)
  points(x$PPV$prevalence,x$PPV$mean,pch=16)
  
  plot(x=x$NPV$prevalence,y=x$NPV$mean,type="l",ylab="Predictive value",col=4,ylim=c(0,1),main="Projected predictive values",xlab="Prevalence")
  lines(x$NPV$prevalence,x$NPV$p97.5, type = "l", lty=2,col=4)
  lines(x$NPV$prevalence,x$NPV$p2.5, type = "l", lty=2,col=4)
  points(x$NPV$prevalence,x$NPV$mean,pch=16,col=4)
  lines(x=x$PPV$prevalence,y=x$PPV$mean,type="l",col=2)
  lines(x$PPV$prevalence,x$PPV$p97.5, type = "l", lty=2,col=2)
  lines(x$PPV$prevalence,x$PPV$p2.5, type = "l", lty=2,col=2)
  points(x$PPV$prevalence,x$PPV$mean,pch=16,col=2)
  legend("bottomright",legend = c("NPV", "PPV"),col = c(4,2), pch = 16, bty="n")
}
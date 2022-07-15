predv_d <- function(x,prop_m,prop_sd,zb=TRUE,n_iter=100000,...){
  object <- reitsma(x)
  m1 <- summary(reitsma(x))
  
  if(zb){
    if(missing(prop_m)&missing(prop_sd)){
      
      PredDist <-function(object, ...) UseMethod("PredDist")
      
      PredDist.default <- function(object, mu,Sigma,alphasens = 1, alphafpr = 1,
                                   n.iter = n_iter, FUN, ...){
        samples <- rmvnorm(n.iter, mu, Sigma)
        sens <- inv.trafo(alphasens,samples[,1])
        fpr <- inv.trafo(alphafpr,samples[,2])
        out <- lapply(FUN, function(x){x(sens, fpr)})
        class(out) <- "PredDist"
        out
      }
      
      PredDist.reitsma <- function(object, prev, n.iter = n_iter, FUN = NULL, ...){
        fit <- object
        if(length(coef(fit)) > 2){
          stop("PredDist is not be used for meta-regression!")}
        if(is.null(FUN)){FUN <- list(npv= function(sens,fpr){((1 - fpr)*(1-prev))/(((1 - fpr)*(1-prev)) + (prev*(1-sens)))},
                                     ppv = function(sens, fpr){(sens*prev)/((sens*prev)+(fpr*(1-prev)))})}
        PredDist.default(mu = coef(fit)[1:2], Sigma = vcov(fit), 
                         alphasens = fit$alphasens, 
                         alphafpr = fit$alphafpr,
                         n.iter = n.iter, FUN = FUN)
      }
      
      m01 <- rma(measure="PLN",xi=x$TP+x$FN,ni=x$TP+x$FN+x$FP+x$TN)
      m01_i2 <- m01$I2
      m01_ub <- exp(m01$ci.ub)
      m01_lb <- exp(m01$ci.lb)
      
      prev_ab <- beta.parms.from.quantiles(q=c(m01_lb,m01_ub))
      dist_prev <- rbeta(n=n_iter,shape1 = prev_ab$a,shape2 = prev_ab$b)
      
      results_pred <- PredDist.reitsma(object,prev=dist_prev)
      
      dens_n <- density(results_pred$npv)
      dens_p <- density(results_pred$ppv)
      
      
      mean_npv <- round(mean(results_pred$npv),3)
      sd_npv <- round(sd(results_pred$npv),3)
      quantiles_npv <- round(quantile(results_pred$npv,c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975)),3)
      quantiles_npv_t <- as.data.frame(t(as.data.frame(quantiles_npv)))
      
      mean_ppv <- round(mean(results_pred$ppv),3)
      sd_ppv <- round(sd(results_pred$ppv),3)
      quantiles_ppv <- round(quantile(results_pred$ppv,c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975)),3)
      quantiles_ppv_t <- as.data.frame(t(as.data.frame(quantiles_ppv)))
      
      quantiles_t <- rbind(quantiles_npv_t,quantiles_ppv_t)
      
      df_pres <- data.frame("Mean"=c(mean_npv,mean_ppv),"SD"=c(sd_npv,sd_ppv),"p2.5"=quantiles_t[,1],"p5"=quantiles_t[,2],"p10"=quantiles_t[,3],"p25"=quantiles_t[,4],"p50"=quantiles_t[,5],"p75"=quantiles_t[,6],"p90"=quantiles_t[,7],"p95"=quantiles_t[,8],"p97.5"=quantiles_t[,9],row.names = c("NPV","PPV"))
      
      if(m01_i2>40){warning("There was substantial heterogeneity (I2 of at least 40%) in the meta-analysis of prevalences. The use of a prevalence probabilty distribution based on primary studies' data may not be fully adequate.")}
      
      warning("Please note that included primary studies may not be the best ones to ascertain the prevalence of a disease/condition, particularly if some of such studies have a case-control design.")
      
      
    }else{

      PredDist <-function(object, ...) UseMethod("PredDist")
      
      PredDist.default <- function(object, mu,Sigma,alphasens = 1, alphafpr = 1,
                                   n.iter = n_iter, FUN, ...){
        samples <- rmvnorm(n.iter, mu, Sigma)
        sens <- inv.trafo(alphasens,samples[,1])
        fpr <- inv.trafo(alphafpr,samples[,2])
        out <- lapply(FUN, function(x){x(sens, fpr)})
        class(out) <- "PredDist"
        out
      }
      
      PredDist.reitsma <- function(object, prev, n.iter = n_iter, FUN = NULL, ...){
        fit <- object
        if(length(coef(fit)) > 2){
          stop("PredDist is not be used for meta-regression!")}
        if(is.null(FUN)){FUN <- list(npv= function(sens,fpr){((1 - fpr)*(1-prev))/(((1 - fpr)*(1-prev)) + (prev*(1-sens)))},
                                     ppv = function(sens, fpr){(sens*prev)/((sens*prev)+(fpr*(1-prev)))})}
        PredDist.default(mu = coef(fit)[1:2], Sigma = vcov(fit), 
                         alphasens = fit$alphasens, 
                         alphafpr = fit$alphafpr,
                         n.iter = n.iter, FUN = FUN)
      }
      
      
      prev_a <- (((prop_m)^2)*(1-(prop_m))/((prop_sd)^2)-(prop_m))
      prev_b <- ((1-(prop_m))*(((1-(prop_m))*(prop_m))/((prop_sd)^2)-1))
      dist_prev <- rbeta(n=n_iter,shape1 = prev_a,shape2 = prev_b)
      
      results_pred <- PredDist.reitsma(object,prev=dist_prev)
      
      dens_n <- density(results_pred$npv)
      dens_p <- density(results_pred$ppv)
      
      
      mean_npv <- round(mean(results_pred$npv),3)
      sd_npv <- round(sd(results_pred$npv),3)
      quantiles_npv <- round(quantile(results_pred$npv,c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975)),3)
      quantiles_npv_t <- as.data.frame(t(as.data.frame(quantiles_npv)))
      
      mean_ppv <- round(mean(results_pred$ppv),3)
      sd_ppv <- round(sd(results_pred$ppv),3)
      quantiles_ppv <- round(quantile(results_pred$ppv,c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975)),3)
      quantiles_ppv_t <- as.data.frame(t(as.data.frame(quantiles_ppv)))
      
      quantiles_t <- rbind(quantiles_npv_t,quantiles_ppv_t)
      
      df_pres <- data.frame("Mean"=c(mean_npv,mean_ppv),"SD"=c(sd_npv,sd_ppv),"p2.5"=quantiles_t[,1],"p5"=quantiles_t[,2],"p10"=quantiles_t[,3],"p25"=quantiles_t[,4],"p50"=quantiles_t[,5],"p75"=quantiles_t[,6],"p90"=quantiles_t[,7],"p95"=quantiles_t[,8],"p97.5"=quantiles_t[,9],row.names = c("NPV","PPV"))
      
  }
    }else{
    if(missing(prop_m)&missing(prop_sd)){
      
      df1 <- data.frame("sens_lb"=1/(1+exp(-(m1$coefficients[1,5]))),"sens_ub"=1/(1+exp(-(m1$coefficients[1,6]))),"spec_lb"=1-(1/(1+exp(-(m1$coefficients[2,6])))),"spec_ub"=1-(1/(1+exp(-(m1$coefficients[2,5])))))
      sens <- beta.parms.from.quantiles(q=c(df1$sens_lb,df1$sens_ub))
      spec <- beta.parms.from.quantiles(q=c(df1$spec_lb,df1$spec_ub))
      
      m01 <- rma(measure="PLN",xi=x$TP+x$FN,ni=x$TP+x$FN+x$FP+x$TN)
      m01_i2 <- m01$I2
      m01_ub <- exp(m01$ci.ub)
      m01_lb <- exp(m01$ci.lb)
      
      prev_ab <- beta.parms.from.quantiles(q=c(m01_lb,m01_ub))
      
      df1 <- data.frame(df1,"sens_a"=sens$a,"sens_b"=sens$b,"spec_a"=spec$a,"spec_b"=spec$b,"prev_a"=prev_ab$a,"prev_b"=prev_ab$b)
      
      
      results_pred <- data.frame(rsens=rbeta(n_iter,df1$sens_a,df1$sens_b),rspec=rbeta(n_iter,df1$spec_a,df1$spec_b),rprev=rbeta(n_iter,df1$prev_a,df1$prev_b))
      results_pred$ppv <- ((results_pred$rsens*results_pred$rprev)/((results_pred$rsens*results_pred$rprev)+((1-results_pred$rspec)*(1-results_pred$rprev))))
      results_pred$npv <- ((results_pred$rspec*(1-results_pred$rprev))/((results_pred$rspec*(1-results_pred$rprev))+(results_pred$rprev*(1-results_pred$rsens))))
      
      dens_n <- density(results_pred$npv)
      dens_p <- density(results_pred$ppv)
      
      
      mean_npv <- round(mean(results_pred$npv),3)
      sd_npv <- round(sd(results_pred$npv),3)
      quantiles_npv <- round(quantile(results_pred$npv,c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975)),3)
      quantiles_npv_t <- as.data.frame(t(as.data.frame(quantiles_npv)))
      
      mean_ppv <- round(mean(results_pred$ppv),3)
      sd_ppv <- round(sd(results_pred$ppv),3)
      quantiles_ppv <- round(quantile(results_pred$ppv,c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975)),3)
      quantiles_ppv_t <- as.data.frame(t(as.data.frame(quantiles_ppv)))
      
      quantiles_t <- rbind(quantiles_npv_t,quantiles_ppv_t)
      
      df_pres <- data.frame("Mean"=c(mean_npv,mean_ppv),"SD"=c(sd_npv,sd_ppv),"p2.5"=quantiles_t[,1],"p5"=quantiles_t[,2],"p10"=quantiles_t[,3],"p25"=quantiles_t[,4],"p50"=quantiles_t[,5],"p75"=quantiles_t[,6],"p90"=quantiles_t[,7],"p95"=quantiles_t[,8],"p97.5"=quantiles_t[,9],row.names = c("NPV","PPV"))
      
      if(m01_i2>40){warning("There was substantial heterogeneity (I2 of at least 40%) in the meta-analysis of prevalences. The use of a prevalence probabilty distribution based on primary studies' data may not be fully adequate.")}
      
      warning("Please note that included primary studies may not be the best ones to ascertain the prevalence of a disease/condition, particularly if some of such studies have a case-control design.")
      
    }else{
      
      df1 <- data.frame("sens_lb"=1/(1+exp(-(m1$coefficients[1,5]))),"sens_ub"=1/(1+exp(-(m1$coefficients[1,6]))),"spec_lb"=1-(1/(1+exp(-(m1$coefficients[2,6])))),"spec_ub"=1-(1/(1+exp(-(m1$coefficients[2,5])))))
      sens <- beta.parms.from.quantiles(q=c(df1$sens_lb,df1$sens_ub))
      spec <- beta.parms.from.quantiles(q=c(df1$spec_lb,df1$spec_ub))
      
      prev_a <- (((prop_m)^2)*(1-(prop_m))/((prop_sd)^2)-(prop_m))
      prev_b <- ((1-(prop_m))*(((1-(prop_m))*(prop_m))/((prop_sd)^2)-1))
      
      df1 <- data.frame(df1,"sens_a"=sens$a,"sens_b"=sens$b,"spec_a"=spec$a,"spec_b"=spec$b,"prev_a"=prev_a,"prev_b"=prev_b)
      
      
      results_pred <- data.frame(rsens=rbeta(n_iter,df1$sens_a,df1$sens_b),rspec=rbeta(n_iter,df1$spec_a,df1$spec_b),rprev=rbeta(n_iter,df1$prev_a,df1$prev_b))
      results_pred$ppv <- ((results_pred$rsens*results_pred$rprev)/((results_pred$rsens*results_pred$rprev)+((1-results_pred$rspec)*(1-results_pred$rprev))))
      results_pred$npv <- ((results_pred$rspec*(1-results_pred$rprev))/((results_pred$rspec*(1-results_pred$rprev))+(results_pred$rprev*(1-results_pred$rsens))))
      
      
      dens_n <- density(results_pred$npv)
      dens_p <- density(results_pred$ppv)

      
      mean_npv <- round(mean(results_pred$npv),3)
      sd_npv <- round(sd(results_pred$npv),3)
      quantiles_npv <- round(quantile(results_pred$npv,c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975)),3)
      quantiles_npv_t <- as.data.frame(t(as.data.frame(quantiles_npv)))
      
      mean_ppv <- round(mean(results_pred$ppv),3)
      sd_ppv <- round(sd(results_pred$ppv),3)
      quantiles_ppv <- round(quantile(results_pred$ppv,c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975)),3)
      quantiles_ppv_t <- as.data.frame(t(as.data.frame(quantiles_ppv)))
      
      quantiles_t <- rbind(quantiles_npv_t,quantiles_ppv_t)
      
      df_pres <- data.frame("Mean"=c(mean_npv,mean_ppv),"SD"=c(sd_npv,sd_ppv),"p2.5"=quantiles_t[,1],"p5"=quantiles_t[,2],"p10"=quantiles_t[,3],"p25"=quantiles_t[,4],"p50"=quantiles_t[,5],"p75"=quantiles_t[,6],"p90"=quantiles_t[,7],"p95"=quantiles_t[,8],"p97.5"=quantiles_t[,9],row.names = c("NPV","PPV"))
    }
  }
  output <- list(results_pred=results_pred,SummaryData=df_pres,density_npv=dens_n,density_ppv=dens_p,call=match.call())
  class(output) <- "predv_d"
  output
}

print.predv_d <- function(x,xlim_npv=c(0,1),xlim_ppv=c(0,1),...){
  
  mean_npv <- x$SummaryData[1,1]
  sd_npv <- x$SummaryData[1,2]
  quantiles_npv <- x$SummaryData[1,3:11]
  
  mean_ppv <- x$SummaryData[2,1]
  sd_ppv <- x$SummaryData[2,2]
  quantiles_ppv <- x$SummaryData[2,3:11]
  
  cat("Negative predictive value")
  cat("\n")
  cat("\n")
  cat("Mean = ", mean_npv,"\n")
  cat("Standard-deviation = ", sd_npv, "\n","\n")
  cat("Quantiles","\n")
  print(quantiles_npv)
  cat("\n")
  cat("\n")
  cat("\n")
  
  cat("Positive predictive value")
  cat("\n")
  cat("\n")
  cat("Mean = ", mean_ppv,"\n")
  cat("Standard-deviation = ", sd_ppv, "\n","\n")
  cat("Quantiles","\n")
  print(quantiles_ppv)
  
  hist(x$results_pred$npv, xlab = "Negative predictive value",xlim=xlim_npv,main="Distribution of negative predictive values")
  
  plot(density(x$results_pred$npv), lwd=2, xlab = "Negative predictive value",xlim=xlim_npv,main="Distribution of negative predictive values")
  
  hist(x$results_pred$ppv, xlab = "Positive predictive value",xlim=xlim_ppv,main="Distribution of positive predictive values")
  
  plot(density(x$results_pred$ppv), lwd=2, xlab = "Positive predictive value",xlim=xlim_ppv,main="Distribution of positive predictive values")
  
  plot(density(x$results_pred$npv), col=4, lwd=2, xlab = "Predictive value",xlim=c(0,1),main="Distribution of projected predictive values",ylim=c(0,max(x$density_npv$y,x$density_ppv$y)))
  lines(density(x$results_pred$ppv),col=2, lwd=2)
  legend("topleft",legend = c("NPV", "PPV"),col = c(4,2), lwd = 2, bty="n")
  
}

summary.predv_d <- function(object,xlim_npv=c(0,1),xlim_ppv=c(0,1),...){
  x <- object
  
  print(x$SummaryData)
  
  hist(x$results_pred$npv, xlab = "Negative predictive value",xlim=xlim_npv,main="Distribution of negative predictive values")
  
  plot(density(x$results_pred$npv), lwd=2, xlab = "Negative predictive value",xlim=xlim_npv,main="Distribution of negative predictive values")
  
  hist(x$results_pred$ppv, xlab = "Positive predictive value",xlim=xlim_ppv,main="Distribution of positive predictive values")
  
  plot(density(x$results_pred$ppv), lwd=2, xlab = "Positive predictive value",xlim=xlim_ppv,main="Distribution of positive predictive values")
  
  plot(density(x$results_pred$npv), col=4, lwd=2, xlab = "Predictive value",xlim=c(0,1),main="Distribution of projected predictive values",ylim=c(0,max(x$density_npv$y,x$density_ppv$y)))
  lines(density(x$results_pred$ppv),col=2, lwd=2)
  legend("topleft",legend = c("NPV", "PPV"),col = c(4,2), lwd = 2, bty="n")
}
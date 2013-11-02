reitsma <- function(data, ...) UseMethod("reitsma")

reitsma.default <-
function(data = NULL, subset=NULL, formula = NULL,
         TP="TP", FN="FN", FP="FP", TN="TN", 
         alphasens = 1, alphafpr = 1, 
         correction = 0.5, correction.control = "all",
         method = "reml",  
         control = list(), ...)
{
  call <- match.call()
  mcall <- match.call(expand.dots = FALSE)
  
stopifnot(is.numeric(correction), 0 <= correction,  
          correction.control %in% c("all", "single", "none"),
          0 <= alphasens, alphasens <= 2, 0 <= alphafpr, alphafpr <= 2,
          method %in% c("fixed", "ml", "reml"),
          is.numeric(TP) | (is.character(TP) & length(TP) == 1),
          is.numeric(FP) | (is.character(FP) & length(FP) == 1),
          is.numeric(TN) | (is.character(TN) & length(TN) == 1),
          is.numeric(FN) | (is.character(FN) & length(FN) == 1))

if(is.null(formula)){formula <- cbind(tsens,tfpr)~1}
if(!class(formula) == "formula"){stop("formula must be of class formula")}
if(!formula[2] == (cbind(tsens,tfpr)~1)[2]){stop("The left hand side of formula must be cbind(tsens, tfpr)")}

if(!is.null(data) & is.character(c(TP,FP,TN,FN))){
X <- as.data.frame(data)
origdata <- data
TP <- getElement(X,TP)
FN <- getElement(X,FN)
FP <- getElement(X,FP)
TN <- getElement(X,TN)
}
   
if(is.null(data) & !is.character(c(TP,FP,TN,FN))){
  origdata <- data.frame(TP = TP, FN = FN, FP = FP, TN = TN)
}

freqdata <- cbind(TP,FN,FP,TN)
checkdata(freqdata)

N <- length(TP)	
	
## apply continuity correction to _all_ studies if one contains zero
if(correction.control == "all"){if(any(c(TP,FN,FP,TN) == 0))
  {TP <- TP + correction;
	 FN <- FN + correction;
	 FP <- FP + correction;
	 TN <- TN + correction}}
if(correction.control == "single"){
	correction = ((((TP == 0)|(FN == 0))|(FP == 0))| (TN == 0)) * 
    correction
  TP <- correction + TP
	FN <- correction + FN
	FP <- correction + FP
	TN <- correction + TN
	}

number.of.pos <- TP + FN
number.of.neg <- FP + TN
sens<-TP/number.of.pos
fpr <- FP/number.of.neg

senstrafo <- function(x){return(mada:::talpha(alphasens)$linkfun(x))}
fprtrafo <- function(x){return(mada:::talpha(alphafpr)$linkfun(x))}
sensinvtrafo <- function(x){return(mada:::talpha(alphasens)$linkinv(x))}
fprinvtrafo <- function(x){return(mada:::talpha(alphafpr)$linkinv(x))}


trafo.sens <- senstrafo(sens)
trafo.fpr <- fprtrafo(fpr)

var.trafo.sens <- (sens*(1-sens)/number.of.pos) * 
  (mada:::jacobiantrafo(alphasens, sens)^2)
var.trafo.fpr <- (fpr*(1-fpr)/number.of.neg) * 
  (mada:::jacobiantrafo(alphafpr, fpr)^2)
  

# preparations for call to mvmeta.fit
mvmeta_S <- cbind(var.trafo.sens, rep(0,N), var.trafo.fpr)
mvmeta_data <- data.frame(X,tsens = trafo.sens, tfpr = trafo.fpr)
mvmeta_X <- model.matrix(formula,data=mvmeta_data)  
mvmeta_y <- as.matrix(data.frame(tsens = trafo.sens, tfpr = trafo.fpr))
  
fit <- mvmeta.fit(mvmeta_X, mvmeta_y, mvmeta_S, method = method, control = control)

fit$call <- call
fit$formula <- formula
fit$control <- control
fit$method <- method
fit$alphasens <- alphasens
fit$alphafpr <- alphafpr
fit$data <- origdata
fit$freqdata <- as.data.frame(freqdata)
  
class(fit) <- "reitsma"

## modify log likelihood, especially add jacobian
attr(fit$logLik, "df") <- 2*ncol(mvmeta_X) + 3
fit$logLik  <-  fit$logLik + sum(log(mada:::jacobiantrafo(alphasens,sens)) + 
  log(mada:::jacobiantrafo(alphafpr,fpr)))
  
fit$nobs <- 2*N

return(fit)

class(output) <- "reitsma"
return(output)

}

print.reitsma <- function (x, digits = 4, ...) 
{
  methodname <- c("reml", "ml", "fixed")
  methodlabel <- c("REML", "ML", "Fixed")
  cat("Call:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  cat("Fixed-effects coefficients:", "\n", sep = "")
  table <- formatC(x$coefficients, digits = digits, format = "f")
  table <- gsub("NA"," -",table)
  print(table, quote = FALSE, right = TRUE, print.gap = 2)
  cat("\n")
  cat(x$dim$m, " studies, ", 
      x$df$fixed, " fixed and ", x$df$random, " random-effects parameters", 
      "\n", sep = "")
  table <- c(x$logLik, AIC(x), BIC(x))
  names(table) <- c("logLik", "AIC", "BIC")
  table <- formatC(table, digits = digits, format = "f")
  print(table, quote = FALSE, right = TRUE, print.gap = 2)
  cat("\n")
}

summary.reitsma <- 
function (object, level = 0.95, ...) 
{
  ci.level <- level
  if (ci.level <= 0 || ci.level >= 1) 
    stop("'ci.level' must be within 0 and 1")
  coef <- object$coefficients
  vcov <- object$vcov
  dim <- object$dim
  Psi <- object$Psi
  lab <- object$lab
  coef <- as.numeric(coef)
  coef.se <- sqrt(diag(vcov))
  zval <- coef/coef.se
  zvalci <- qnorm((1 - ci.level)/2, lower.tail = FALSE)
  pvalue <- 2 * (1 - pnorm(abs(zval)))
  ci.lb <- coef - zvalci * coef.se
  ci.ub <- coef + zvalci * coef.se
  cilab <- paste(signif(ci.level, 2) * 100, "%ci.", c("lb", 
                                                      "ub"), sep = "")
  tabfixed <- cbind(coef, coef.se, zval, pvalue, ci.lb, ci.ub)
  
    
  dimnames(tabfixed) <- list(if (dim$k > 1L) colnames(vcov) else lab$p, 
                             c("Estimate", "Std. Error", "z", "Pr(>|z|)", cilab))
  
  if(nrow(tabfixed) == 2){
    tabfixed <- rbind(tabfixed,tabfixed)
    tabfixed[3,] <- mada:::inv.trafo(object$alphasens, tabfixed[3,])
    tabfixed[4,] <- mada:::inv.trafo(object$alphafpr, tabfixed[4,])
    tabfixed[3:4,c(2,3,4)] <- NA
    rownames(tabfixed)[3:4] <- c("sensitivity","false pos. rate")
  }
    
  corFixed <- vcov/outer(coef.se, coef.se)
  ran.sd <- sqrt(diag(Psi))
  corRandom <- Psi/outer(ran.sd, ran.sd)
#  qstat <- unclass(qtest(object))  
  if(attr(object$logLik,"df")== 5){AUC <- AUC(object)}else{AUC <- NULL}
  
  if(attr(object$logLik,"df")== 5){
  ## HSROC parameters (formulae from Harbord et al)
  Theta <- 0.5*(sqrt(ran.sd[2]/ran.sd[1])*coef[1] + sqrt(ran.sd[1]/ran.sd[2])*coef[2]) ##coef[2] is the fpr, so need to change sign!
  Lambda <- sqrt(ran.sd[2]/ran.sd[1])*coef[1] - sqrt(ran.sd[1]/ran.sd[2])*coef[2] ##coef[2] is the fpr, so need to change sign!
  sigma2theta <- 0.5*(ran.sd[1]*ran.sd[2] + Psi[1,2]) ##coef[2] is the fpr, so need to change sign of Psi[1,2] as well!
  sigma2alpha <- 2*(ran.sd[1]*ran.sd[2] - Psi[1,2])  ##coef[2] is the fpr, so need to change sign of Psi[1,2] as well!
  beta <- log(ran.sd[2]/ran.sd[1])             
    coef_hsroc <- list(Theta = Theta,
                       Lambda = Lambda,
                       beta = beta,
                       sigma2theta = sigma2theta,
                       sigma2alpha = sigma2alpha)
  coef_hsroc <- lapply(coef_hsroc, function(x){attr(x, "names") <- NULL; x})
  }else{
    coef_hsroc = NULL
  }
  keep <- match(c("vcov", "Psi", "df.res", "rank", "logLik", 
                  "converged", "dim", "df", "lab", "na.action", "call", 
                  "terms", "method", "alphasens", "alphafpr"), names(object), 0L)
  out <- c(list(coefficients = tabfixed), object[keep], list(AIC = AIC(object), 
                                                             BIC = BIC(object), 
                                                             corFixed = corFixed, 
                                                             corRandom = corRandom, 
#                                                             qstat = qstat, 
                                                             ci.level = ci.level,
                                                             AUC = AUC,
                                                             coef_hsroc = coef_hsroc))
  class(out) <- "summary.reitsma"
  return(out)
}

print.summary.reitsma <- function (x, digits = 3, ...){
    
    methodname <- c("reml", "ml", "fixed")
    methodlabel <- c("REML", "ML", "Fixed")
    cat("Call:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
     cat("Bivariate diagnostic ", 
         ifelse(x$method == "fixed", "fixed", "random"), 
         "-effects meta-", ifelse(x$dim$p > 2,"regression", "analysis"), "\n", sep = "")
    if (x$method != "fixed") {
      cat("Estimation method: ", methodlabel[which(x$method == 
        methodname)], "\n", sep = "")
#      cat("Variance-covariance matrix Psi: ", "unstructured", 
#          "\n", sep = "")
    }
    cat("\n")
    cat("Fixed-effects coefficients", "\n", sep = "")
    signif <- symnum(x$coefficients[, "Pr(>|z|)"], corr = FALSE, 
                     na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 
                                               1), 
                     symbols = c("***", "**", "*", ".", " "))
    tabletot <- formatC(x$coefficients, digits = digits, format = "f")
    tabletot <- gsub("NA"," -",tabletot)
    tabletot <- cbind(tabletot, signif)
    colnames(tabletot)[7] <- ""
      print(tabletot, quote = FALSE, right = TRUE, print.gap = 1)
    cat("---\nSignif. codes: ", attr(signif, "legend"), "\n\n")
    if (!x$method == "fixed") {
      cat("Variance components: between-studies Std. Dev and correlation matrix", 
          "\n", sep = "")
      corRan <- x$corRandom
      corRan[upper.tri(x$corRan)] <- NA
      table <- cbind(`Std. Dev` = sqrt(diag(x$Psi)), corRan)
      if (x$dim$k == 1L) 
        rownames(table) <- ""
      table <- formatC(table, digits = digits, format = "f")
      table[grep("NA", table)] <- "."
      print(table, quote = FALSE, right = TRUE, na.print = "", 
            print.gap = 1)
      cat("\n")
    }
    table <- c(x$logLik, x$AIC, x$BIC)
    names(table) <- c("logLik", "AIC", "BIC")
    table <- formatC(table, digits = digits, format = "f")
    print(table, quote = FALSE, right = TRUE, print.gap = 1)
    cat("\n")

  if(!is.null(x$AUC)){
    cat(c("AUC: ",as.character(round(x$AUC$AUC,3))))
    cat("\n")  
    cat(c("Partial AUC (restricted to observed FPRs and normalized): ",as.character(round(x$AUC$pAUC,3))))
    cat("\n")  
  }
  
  if(!is.null(x$coef_hsroc)){
    cat("\n")
    cat("HSROC parameters \n")
    table <- as.numeric(x$coef_hsroc)
    names(table) <- names(x$coef_hsroc)
    table <- formatC(table, digits = digits, format = "f")
    print(table, quote = FALSE, right = TRUE, print.gap = 1)
    cat("\n")
  }
  
}
  

vcov.reitsma <- function(object, ...){object$vcov}

logLik.reitsma <- function(object, ...){object$logLik}

sroc.reitsma <- function(fit, fpr = 1:99/100, type = "ruttergatsonis", ...){
  stopifnot(type %in% c("ruttergatsonis", "naive"))
  if(type == "naive"){return(mada:::sroc2(fit, fpr=fpr, ...))}
  if(type == "ruttergatsonis"){
    sum_fit <- summary(fit)
    Lambda <- summary(fit)$coef_hsroc$Lambda    
    Beta <- summary(fit)$coef_hsroc$beta    
    return(cbind(fpr, 
#      sens = (1+exp(-Lambda*exp(-0.5*Beta)-exp(-Beta)*log(fpr/(1-fpr))))^(-1)))                 
      sens = mada:::inv.trafo(fit$alphasens, (Lambda*exp(-Beta/2) + 
          exp(-Beta)*mada:::trafo(fit$alphafpr, fpr)))))
  }
}

mcsroc.reitsma <- function(fit, fpr = 1:99/100, replications = 10000, lambda = 100, ...){
  mcsroc(fit, replications = replications, lambda = lambda, ...)
}
  
ROCellipse.reitsma <- function(x, level = 0.95, add = FALSE, pch = 1, ...){
  mada:::ROC.ellipse2(x, nobs = x$nobs/2, conf.level = level, add = add, pch = pch, ...)
}

crosshair.reitsma <- function(x, level = 0.95, length = 0.1, pch = 1, ...){
  s <- summary(x, level = level)  
  ch <- rbind(inv.trafo(x$alphasens,s$coefficients["tsens.(Intercept)",5:6]),
              inv.trafo(x$alphafpr,s$coefficients["tfpr.(Intercept)",5:6]))
  pe <- c(inv.trafo(x$alphafpr,s$coefficients["tfpr.(Intercept)",1]),
          inv.trafo(x$alphasens,s$coefficients["tsens.(Intercept)",1]))
#  if(!add){
#    plot(matrix(pe, ncol = 2), pch = pch, xlim = xlim, ylim = ylim, 
#         ylab = "Sensitivity", xlab = "False Positive Rate", ...)
#  }else{
    points(matrix(pe, ncol = 2), pch = pch, ...)
#  }
  arrows(ch[2,1], pe[2], x1 =ch[2,2],  length = length, angle = 90, code = 3, ...)
  arrows(pe[1], ch[1,1],  y1 =ch[1,2],  length = length, angle = 90, code = 3, ...)
  return(invisible(NULL))
}
  
plot.reitsma <- function(x, extrapolate = FALSE, plotsumm = TRUE, level = 0.95, 
                         ylim = c(0,1), xlim = c(0,1), pch = 1, 
                         sroclty = 1, sroclwd = 1, 
                         predict = FALSE, predlty = 3, predlwd = 1,
                         type = "ruttergatsonis",
                         ...)
{
  plot(c(2,2), ylim = ylim, xlim = xlim, 
       xlab = "False Positive Rate", ylab = "Sensitivity", ...)
  if(length(coef(x)) == 2){
  FP <- x$freqdata$FP
  negatives <- FP + x$freqdata$TN
  FPR <- FP/negatives
  
  if(extrapolate){bound = c(0,1)}
  if(!extrapolate){bound = c(min(FPR), max(FPR))}
  srocmat <- sroc(x, type = type)
  lines(srocmat[cut(srocmat[,1],bound, "withinbound") == "withinbound",], 
        lty = sroclty, lwd = sroclwd)
  }else{
    warning("Not plotting any SROC for meta-regression")
  }
  if(plotsumm){ROCellipse(x, level = level, add = TRUE, pch = pch, ...)}
  
  if(predict){
    alpha.sens <- x$alphasens
    alpha.fpr <- x$alphafpr
    mu <- x$coefficients["(Intercept)",]
    Sigma <- x$Psi + vcov(x)
    talphaellipse <- ellipse(Sigma, centre = mu, level = level)
    predellipse <- matrix(0, ncol = 2, nrow = nrow(talphaellipse))
    predellipse[,1] <- mada:::inv.trafo(alpha.fpr, talphaellipse[,2])
    predellipse[,2] <- mada:::inv.trafo(alpha.sens, talphaellipse[,1])
    lines(predellipse, lty = predlty, lwd = predlwd)
}
  
  return(invisible(NULL))
}

confint.reitsma <- function (object, parm, level = 0.95, ...) 
{
  cf <- as.numeric(coef(object))
  if(missing(parm)){parm <- colnames(vcov(object))}

  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- stats:::format.perc(a, 3)
  fac <- qnorm(a)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
                                                             pct))
  ses <- sqrt(diag(vcov(object)))
  ci[] <- cf + ses %o% fac
  ci
}


AUC.reitsma <- function(x, fpr = 1:99/100, ...){
  if(!attr(x$logLik,"df") == 5){stop("AUC can not be calculated for meta-regression")}
  estimate <- x$coefficients
  alpha.sens <- x$alphasens
  alpha.fpr <- x$alphafpr
  mu1 <- estimate[1]
  mu2 <- estimate[2]
  sigma2 <- x$Psi[2,2]
  sigma <- x$Psi[1,2]  
  rsroc <- function(x){mada:::calc.sroc(x, alpha.sens, alpha.fpr, mu1, mu2, sigma2, sigma)}
  AUC <- mada:::AUC.default(rsroc, fpr = fpr, ...)
  obsfprrange <- range(fpr(x$freqdata))
  obsfprrange[1] <- max(0.01,obsfprrange[1])
  obsfprrange[2] <- min(0.99,obsfprrange[2])
  obsfpr <- seq(from = obsfprrange[1], to = obsfprrange[2], length.out = 99)
  pAUC <- mada:::AUC.default(rsroc, fpr = obsfpr, ...)
  names(pAUC) <- c("pAUC")
  return(c(AUC,pAUC))
}



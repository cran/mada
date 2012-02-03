reitsma <- function(X, ...) UseMethod("reitsma")

reitsma.default <-
function(X = NULL, TP, FN, FP, TN, correction = 0.5, correction.control = "all",
         REML = TRUE, ...)
{
if(!is.null(X)){
X <- as.data.frame(X)
origdata <- X
TP <- X$TP
FN <- X$FN
FP <- X$FP
TN <- X$TN
}
if(is.null(X)){
  origdata <- data.frame(TP = TP, FN = FN, FP = FP, TN = TN)
}

checkdata(origdata)

N <- length(TP)	
	
## apply continuity correction to _all_ studies if one contains zero
if(correction.control == "all"){if(any(c(TP,FN,FP,TN) == 0)){TP <- TP+correction;
							FN <- FN + correction;
							FP <- FP + correction;
							TN <- TN + correction}}
if(correction.control == "single"){
	correction = ((((TP == 0)|(FN == 0))|(FP == 0))| (TN == 0))*correction
  TP <- correction + TP
	FN <- correction + FN
	FP <- correction + FP
	TN <- correction + TN
	}

number.of.pos <- TP + FN
number.of.neg <- FP + TN
sens<-TP/number.of.pos
fpr <- FP/number.of.neg

logit.sens <- logit(sens)
logit.fpr <- logit(fpr)
var.logit.sens <- 1/(sens*(1-sens)*number.of.pos)
var.logit.fpr <- 1/(fpr*(1-fpr)*number.of.neg)

# calculate individual C_i as array
C <- array(data = 0, dim = c(2,2,N))
C[1,1,] <- var.logit.sens
C[2,2,] <- var.logit.fpr

# preparations for likelihood calculations
X <- cbind(logit.sens, logit.fpr)

# negative log-likelihood
# the parameters are set up as a vector so they can be handled better by nlm



## we can reduce the likelihood to three parameters by profiling for mu!
## using the below likelihood speeds up convergence by a factor of 25% for AuditC data
## is slower for Dementia because of so many inverses! 

## function hat.mu calculate mu given Sigma

hat.mu <- function(Sigma){
	Sigmai <- array(data = 0, dim = c(2,2,N))
	
	for(i in 1:N){Sigmai[,,i] <- solve(Sigma + C[,,i])}
	
	## Q: Geht das mit apply? 
	## A: ja, zb.: Sigmai <- array(apply(C,3,solve), dim = c(2,2,N))
	## da aber C+Sigma nicht sinnvoll ist, ist das oben einfacher und vermutlich genauso schnell 
	
	# now calculate the weighted sum:  \sum_i Sigma_i^{-1} %*% X[i,]
	
	weighted.X <- c(0,0)
	
	for(i in 1:N){weighted.X <- weighted.X + (Sigmai[,,i] %*% X[i,])}
	
	## auch das obige geht vermutlich irgendwie mit apply
	
	solve(apply(Sigmai, c(1,2), sum)) %*% weighted.X
}

# later perform one step of nlm with the parametrization using Sigma and not its cholesky factorization to obtain variance estimates

negfullloglik <- function(pars)
  {
	mu <- pars[1:2]
	Sigma <- matrix( c(pars[3], pars[5], pars[5], pars[4]), nrow = 2)
  p <- numeric(N)
  
    for(i in 1:N){p[i] <- dmvnorm(X[i,], mu, Sigma+C[,,i], log= TRUE)}
	  
  Sigmai <- array(data = 0, dim = c(2,2,N))
	
	for(i in 1:N){Sigmai[,,i] <- solve(Sigma + C[,,i])}
  
  REML.adjustment <- ifelse(REML, - 0.5* log(det((apply(Sigmai, c(1,2), sum))/(2*pi))), 0)
  
  output <- -(sum(p) + REML.adjustment + sum(log(jacobiantrafo(1, sens))) + sum(log(jacobiantrafo(1, fpr))))
	
#  gradientMat <- matrix(0, ncol = 2, nrow = 2)
#  
#  for(i in 1:N)
#  {
#  gradientMat <- gradientMat + 0.5*(-Sigmai[,,i] + (Sigmai[,,i]%*%(X[i,] - mu))%*%(t(X[i,]-mu)%*%Sigmai[, ,i]))
#  }
#    
#  attr(output, "gradient") <- c(gradientMat[1,1], gradientMat[2,2], gradientMat[1,2])
  
  return(output)
  }



negloglik <- function(pars){
	L <- matrix(c(pars[1], 0, pars[3], pars[2]), ncol = 2, byrow = FALSE)
	Sigma <- t(L) %*% L

	## conditional on Sigma the mean mu depends only on data!
	## hat.mu is mu conditional on Sigma
	
	mu <- hat.mu(Sigma)
		
	output <- negfullloglik(c(mu, Sigma[1,1], Sigma[2,2], Sigma[1,2]))
  
  return(output)
  
	}

	
# compute some start values based on covariance matrix
start.Sigma <- cov(cbind(logit.sens, logit.fpr))
start.L <- chol(start.Sigma)
start.pars <- c(start.L[1,1], start.L[2,2], start.L[1,2])
# if(alpha){alpha.start <- reitsma.start(sens,fpr)$maximum}


# we do not need to compute starting values for mu, since mu is profiled

		fit <- nlm(negloglik, start.pars, hessian = FALSE)

    iterations <- fit$iterations

		
	
## use the estimates gained from cholesky run as starting values, for this compute Sigma
L <- matrix(c(fit$estimate[1], fit$estimate[3], 0, fit$estimate[2]), ncol = 2)
Sigma <- L %*% t(L)
# compute mu.hat from Sigma obtained from fit:
mu <- hat.mu(Sigma)


start.pars2 <- c(mu, Sigma[1,1], Sigma[2,2], Sigma[1,2])

output2 <- try(nlm(negfullloglik, start.pars2, hessian = TRUE))

if(class(output2) == "try-error"){cat("nlm failed, trying optim with method = SANN \n")
  output2 <- optim(start.pars2, negfullloglik, hessian = TRUE, method = "SANN")}

  if(!is.null(output2$par)){coef <- output2$par; logLik <- output2$value}
  if(is.null(output2$par)){coef <- output2$estimate; logLik <- output2$minimum}
  
  attr(logLik, "df") <- 5

  names(coef) <- c("mu1", "mu2", "Sigma11", "Sigma22", "Sigma21")

  hessian <- output2$hessian
  
  Vcov <- solve(hessian)
  colnames(Vcov) <- c("mu1", "mu2", "Sigma11", "Sigma22", "Sigma21")
  rownames(Vcov) <- c("mu1", "mu2", "Sigma11", "Sigma22", "Sigma21")

  df <- N - 5

output <- list(coefficients = coef, vcov = Vcov, df = df, nobs = 2*N, logLik = logLik,
               iterations = (iterations+1), call = match.call(), REML=REML, data = origdata)

class(output) <- "reitsma"
return(output)

}

print.reitsma <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\nCoefficients:\n")
print(x$coefficients)
}


summary.reitsma <- function(object, level = 0.95, ...)
{
confints <- cbind(object$coefficients, confint(object))
colnames(confints)[1] <- "Estimate"
confints <- rbind(inv.trafo(1,confints[1:2,]),confints)
rownames(confints)[1:2] <- c("Sens", "FPR")

res <- list(call=object$call,
confints = confints)
class(res) <- "summary.reitsma"
res
}

print.summary.reitsma <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\n")
print(x$confints)
}

vcov.reitsma <- function(object, ...){object$vcov}

logLik.reitsma <- function(object, ...){object$logLik}

sroc.reitsma <- function(fit, fpr = 1:99/100, ...){
  sroc2(fit, fpr=fpr, ...)
}

mcsroc.reitsma <- function(fit, fpr = 1:99/100, replications = 10000, lambda = 100, ...){
  mcsroc(fit, replications = replications, lambda = lambda, ...)
}
  
ROCellipse.reitsma <- function(x, level = 0.95, add = FALSE, pch = 1, ...){
  ROC.ellipse2(x, nobs = x$nobs/2, conf.level = level, add = add, pch = pch, ...)
}

crosshair.reitsma <- function(x, level = 0.95, length = 0.1, pch = 1, ...){
  ch <- inv.logit(confint(x, level = level)[c("mu1", "mu2"),])
  pe <- inv.logit(coef(x)[c("mu2", "mu1")])
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
                         sroclty = 1, sroclwd = 1, ...)
{
  FP <- x$data$FP
  negatives <- FP + x$data$TN
  FPR <- FP/negatives
  
  if(extrapolate){bound = c(0,1)}
  if(!extrapolate){bound = c(min(FPR), max(FPR))}
  plot(c(2,2), ylim = ylim, xlim = xlim, 
       xlab = "False Positive Rate", ylab = "Sensitivity", ...)
  srocmat <- sroc(x)
  lines(srocmat[cut(srocmat[,1],bound, "withinbound") == "withinbound",], lty = sroclty, lwd = sroclwd)
  if(plotsumm){ROCellipse(x, level = level, add = TRUE, pch = pch, ...)}
  return(invisible(NULL))
}




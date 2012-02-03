sroc <- function(fit, ...) UseMethod("sroc")

### calculate sroc curves
calc.sroc <- function(fpr, alpha.sens, alpha.fpr, mu1, mu2, sigma2, sigma){
  theta <- sigma/sigma2
  return(inv.trafo(alpha.sens, (mu1 - theta*mu2) + theta*trafo(alpha.fpr,fpr)))
}

sroc2 <- function(fit, fpr = 1:99/100){
  estimate <- fit$coefficients
  alpha.sens <- 1
  alpha.fpr <- 1
  if(length(estimate) == 7){alpha.sens <- estimate[6]; alpha.fpr <- estimate[7]}
  mu1 <- estimate[1]
  mu2 <- estimate[2]
  sigma2 <- estimate[4]
  sigma <- estimate[5]
  return(cbind(fpr, sens = calc.sroc(fpr, alpha.sens, alpha.fpr, mu1, mu2, sigma2, sigma)))
}

mcsroc <- function(fit, ...) UseMethod("mcsroc")

mcsroc<- function(fit, replications = 10000, lambda = 100){
 estimate <- fit$coefficients
  alpha.sens <- 1
  alpha.fpr <- 1
  if(length(estimate) == 7){alpha.sens <- estimate[6]; alpha.fpr <- estimate[7]}
  mu <- c(estimate[1], estimate[2])
  sigma1 <- estimate[3]
  sigma2 <- estimate[4]
  sigma <- estimate[5]
  Sigma <- matrix(c(sigma1, sigma, sigma, sigma2), ncol =2)
  stud.pars <- rmvnorm(replications, mu, Sigma)
  sens <- inv.trafo(alpha.sens, stud.pars[,1])
  fpr <- inv.trafo(alpha.fpr, stud.pars[,2])
  N.sens <- rpois(replications, lambda)
  N.fpr <- rpois(replications, lambda)
  TN <- rbinom(replications, N.sens, sens)
  FP <- rbinom(replications, N.fpr, fpr)
  obs.sens <- TN/N.sens
  obs.fpr <- FP/N.fpr
  return(lowess(obs.fpr, obs.sens))
}

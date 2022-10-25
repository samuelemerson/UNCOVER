###############################################
## R script for functions in UNCOVER package ##
###############################################
##
## Sam Emerson
## September 2022
##


#################################################################################
## Function for Sequential Monte Carlo Sampler to Obtain the Bayesian Evidence ##
#################################################################################


##' SMC sampler function using iterated batch importance sampling
##'
##'
##' @export
##' @name MIBIS.Z
##' @description This function uses an iterated batch importance sampling scheme with batch size one to go from one bridging distribution to another. We assume a Bayesian logistic regression model.
##'
##' The default setting is for the initial 'bridging' distribution to be the prior and the final 'bridging' distribution to be the full posterior.
##'
##' Used in UNCOVER to generate the log sub-Bayesian evidence of partitioned models, however the weighted samples can also be used for posterior inference.
##'
##'
##' @keywords sequential monte carlo
##' @param X Design matrix
##' @param y Binary response vector
##' @param sampl Matrix of samples from either the prior (default setting) or from the bridging distribution (partial posterior) formed by the observations in `in_set`
##' @param prior_pdf Probability Density Function of the prior. Must only have two arguments, `th` and `di` (a vector or matrix of regression coefficients samples and the number of dimensions of a single sample respectively).
##' @param add_set Vector of observation indices to be added to arrive at the desired bridging distribution. If the desired final distribution is the full posterior this should not be specified.
##' @param in_set Vector of observation indices already added to give the starting bridging distribution. If the starting 'bridging' distribution is the prior then this should not be specified.
##' @param w Vector of weights corresponding to the matrix of samples `sampl`. If not specified all samples are assumed to have equal weight.
##' @param logZ Log Bayesian evidence of the partial posterior formed by the observations in `in_set`. Defaults to `0` when `in_set` is not specified.
##' @return A list consisting of; weighted posterior samples (the samples and weights are given in separate lists) and the log Bayesian evidence of the full posterior.
##' @details The default setting is going from prior to full posterior. If using samples from a partial posterior instead of prior samples then `in_set` and `logZ` should be provided. If the samples from the partial posterior are weighted samples then `w` should also be provided. If the full posterior is not the desired output then the specific bridging distribution (partial posterior) required must be specified by the observations wish are to added to the starting partial posterior (i.e. `add_set` is required).
##'
##' Details of the internal mechanisms of the SMC sampler such as the Metropolis-Hastings MCMC resample move can be found in *UNCOVER paper* and Chopin (2002).
##' @references Chopin, N. (2002). A sequential particle filter method for static models. Biometrika, 89(3), 539-552.
##' @examples
##'
##' # First we generate a design matrix X and binary response vector y
##' DM <- cbind(rep(1,100),matrix(rnorm(200),100,2))
##' rv <- sample(0:1,100,replace=T)
##'
##' # We assume the prior for the regression coefficients is a standard normal
##' pr_fun <- function(th,di){return(dmvn(th,mu=rep(0,di),sigma=diag(di)))}
##'
##' # Now we can obtain 1000 samples from the posterior using only 1000 samples
##' # from the prior
##' out.1 <- MIBIS.Z(X = DM,y = rv,sampl = rmvn(1000,rep(0,3),diag(3)),prior_pdf = pr_fun)
##' pairs(out.1[[1]])
##' hist(out.1[[2]])
##' out.1[[3]]
##'
##' # If we then added 100 more observations to our data set
##' DM.2 <- rbind(DM,cbind(rep(1,100),matrix(rnorm(200),100,2)))
##' rv.2 <- c(rv,sample(0:1,100,replace=T))
##'
##' # and then wanted samples from the entire posterior instead of starting
##' # from the prior we can use the previous SMC sampler output
##' out.2 <- MIBIS.Z(X = DM.2,y = rv.2,sampl = out.1[[1]],prior_pdf = pr_fun,in_set = 1:100,w = out.1[[2]],logZ = out.1[[3]])
##' pairs(out.2[[1]])
##' hist(out.2[[2]])
##' out.2[[3]]
##'
##' # we can also observe roughly how the particles move when the 100 new
##' # observations are added in two stages
##' out.3 <- MIBIS.Z(X = DM.2,y = rv.2,sampl = out.1[[1]],prior_pdf = pr_fun,add_set = 101:150,in_set = 1:100,w = out.1[[2]],logZ = out.1[[3]])
##' out.4 <- MIBIS.Z(X = DM.2,y = rv.2,sampl = out.3[[1]],prior_pdf = pr_fun,in_set = 1:150,w = out.3[[2]],logZ = out.3[[3]])
##' pairs(rbind(out.1[[1]],out.3[[1]],out.4[[1]]),col=rep(c('red','green','blue'),each=1000))

MIBIS.Z <- function(X,y,sampl,prior_pdf,add_set=NULL,in_set=NULL,w=NULL,logZ=0){
  if(is.null(w)){
    w <- rep(1,nrow(sampl))
  }
  if(is.null(add_set)){
    add_set <- setdiff(1:length(y),in_set)
  }
  X <- as.matrix(X)
  N <- nrow(sampl)
  p <- ncol(X)
  if(anyDuplicated(sampl)==0){
    samp_t <- w
    samp_i <- 1:N
    isne <- FALSE
  } else{
    isne <- TRUE
  }
  add_set <- add_set[sample(1:length(add_set),length(add_set))]
  for(i in 1:length(add_set)){
    w_new <- w*as.vector(((1+exp(-X[add_set[i],]%*%t(sampl)))^(-y[add_set[i]]))*((1+exp(X[add_set[i],]%*%t(sampl)))^(y[add_set[i]]-1)))
    logZ <- logZ + log(sum(w_new)) - log(sum(w))
    w <- w_new
    if(isne){
      w_hat <- aggregate(w~.,data = data.frame(sampl,w),FUN=sum)$w
      wsq <- sum(w_hat^2)
    } else{
      wsq <- sum((w[samp_i]*samp_t)^2)
    }
    if(wsq==0){
      ESS <- 0
    } else{
      ESS <- (sum(w)^2)/wsq
    }
    if(ESS < N/2){
      ss <- cov.wt(sampl,wt=w,method = "ML")
      mu <- ss$center
      Sigma <- ss$cov
      samp <- sample(1:N,N,prob = w,replace = TRUE)
      sampl <- sampl[samp,]
      w <- rep(1,N)
      BC <- mvnfast::rmvn(N,mu = mu, sigma=Sigma)
      Log_1 <- log(1+exp(X[c(in_set,add_set[i]),]%*%t(BC)))*(y[c(in_set,add_set[i])]-1) + log(1+exp(-X[c(in_set,add_set[i]),]%*%t(BC)))*(-y[c(in_set,add_set[i])])
      Log_1 <- colSums(Log_1) + log(prior_pdf(th = BC,di = p)) + log(mvnfast::dmvn(sampl,mu=mu,sigma=Sigma))
      Log_2 <- log(1+exp(X[c(in_set,add_set[i]),]%*%t(sampl)))*(y[c(in_set,add_set[i])]-1) + log(1+exp(-X[c(in_set,add_set[i]),]%*%t(sampl)))*(-y[c(in_set,add_set[i])])
      Log_2 <- colSums(Log_2) + log(prior_pdf(th = sampl,di = p)) + log(mvnfast::dmvn(BC,mu=mu,sigma=Sigma))
      A <- Log_1 - Log_2 - log(runif(length(Log_1)))
      sampl[which(A>0),] <- BC[which(A>0),]
      samp[which(A>0)] <- (N+1):(N+length(which(A>0)))
      samp_t <- table(samp)
      samp_i <- match(as.integer(names(samp_t)),samp)
      if(isne){
        isne <- FALSE
      }
    }
    in_set <- c(in_set,add_set[i])
  }
  return(list(sampl,w,logZ))
}

MIBIS.Z2 <- function(X,y,sampl,prior_pdf,add_set=NULL,in_set=NULL){
  X <- as.matrix(X)
  p <- ncol(X)
  if(is.null(add_set)){
    add_set <- setdiff(1:length(y),in_set)
  }
  if(is.matrix(sampl)) {
    N <- nrow(sampl)
    w <- rep(1,nrow(sampl))
    logZ <- 0
    if(anyDuplicated(sampl)==0){
      samp_t <- w
      samp_i <- 1:N
      isne <- FALSE
    } else{
      isne <- TRUE
    }
  } else {
    isne <- FALSE
    N <- nrow(sampl$sampl)
    w <- sampl$w
    logZ <- sampl$logZ
    samp_t <- sampl$samp_t
    samp_i <- sampl$samp_i
    sampl <- sampl$sampl
  }
  add_set <- add_set[sample(1:length(add_set),length(add_set))]
  for(i in 1:length(add_set)){
    w_new <- w*as.vector(((1+exp(-X[add_set[i],]%*%t(sampl)))^(-y[add_set[i]]))*((1+exp(X[add_set[i],]%*%t(sampl)))^(y[add_set[i]]-1)))
    logZ <- logZ + log(sum(w_new)) - log(sum(w))
    w <- w_new
    if(isne){
      w_hat <- aggregate(w~.,data = data.frame(sampl,w),FUN=sum)$w
      wsq <- sum(w_hat^2)
    } else{
      wsq <- sum((w[samp_i]*samp_t)^2)
    }
    if(wsq==0){
      ESS <- 0
    } else{
      ESS <- (sum(w)^2)/wsq
    }
    if(ESS < N/2){
      ss <- cov.wt(sampl,wt=w,method = "ML")
      mu <- ss$center
      Sigma <- ss$cov
      samp <- sample(1:N,N,prob = w,replace = TRUE)
      sampl <- sampl[samp,]
      w <- rep(1,N)
      BC <- mvnfast::rmvn(N,mu = mu, sigma=Sigma)
      Log_1 <- log(1+exp(X[c(in_set,add_set[i]),]%*%t(BC)))*(y[c(in_set,add_set[i])]-1) + log(1+exp(-X[c(in_set,add_set[i]),]%*%t(BC)))*(-y[c(in_set,add_set[i])])
      Log_1 <- colSums(Log_1) + log(prior_pdf(th = BC,di = p)) + log(mvnfast::dmvn(sampl,mu=mu,sigma=Sigma))
      Log_2 <- log(1+exp(X[c(in_set,add_set[i]),]%*%t(sampl)))*(y[c(in_set,add_set[i])]-1) + log(1+exp(-X[c(in_set,add_set[i]),]%*%t(sampl)))*(-y[c(in_set,add_set[i])])
      Log_2 <- colSums(Log_2) + log(prior_pdf(th = sampl,di = p)) + log(mvnfast::dmvn(BC,mu=mu,sigma=Sigma))
      A <- Log_1 - Log_2 - log(runif(length(Log_1)))
      sampl[which(A>0),] <- BC[which(A>0),]
      samp[which(A>0)] <- (N+1):(N+length(which(A>0)))
      samp_t <- table(samp)
      samp_i <- match(as.integer(names(samp_t)),samp)
      if(isne){
        isne <- FALSE
      }
    }
    in_set <- c(in_set,add_set[i])
  }
  return(list(sampl = sampl,
              w = w,
              logZ = logZ,
              samp_t = samp_t,
              samp_i = samp_i))
}



##' Log Bayesian evidence generator
##'
##'
##' @export
##' @name lbe.gen
##' @description Estimates the log Bayesian evidence (log(Z)) of a Bayesian logistic regression model.
##'
##' @keywords sequential monte carlo
##' @param method Method to be used to estimate log(Z), can be one of `"SMC"`,`"SMC_BIC"` or `"BIC"`
##' @param thres The threshold for which the number of observations needs to exceed to consider using BIC as an estimator. Only applies if method `"SMC_BIC"` is selected. Defaults to infinity if not specified.
##' @param obs_mat Covariate matrix
##' @param res_vec Binary response vector
##' @param p_num Number of samples of the prior used for the SMC sampler. Not required if method `"BIC"` selected.
##' @param rpri Function to sample from the prior. Must only have two arguments, `p_n` and `di` (Number of prior samples to generate and the number of dimensions of a single sample respectively).
##' @param p_pdf Probability Density Function of the prior. Must only have two arguments, `th` and `di` (a vector or matrix of regression coefficients samples and the number of dimensions of a single sample respectively).
##' @return An estimation of the log Bayesian evidence
##' @details log(Z) is estimated using three possible methods:
##'
##' 1. `"SMC"`: Estimates log(Z) using an SMC sampler. See `MIBIS.Z` for more details.
##'
##' 2. `"SMC_BIC"`: Initially tries to estimate log(Z) using BIC but will revert to using an SMC sampler if the data is linearly separable. Also reverts to using an SMC sampler if the number of observations is below a certain threshold. Recommended method and is set as default.
##'
##' 3. `"BIC"`: Estimates log(Z) with the Bayesian Information Criterion (BIC). Does not require prior specification but only advisable if the number of observations is large.
##' @examples
##'
##' # First we generate a covariate matrix `obs_mat` and binary response vector `res_vec`
##' CM <- matrix(rnorm(2000),1000,2)
##' rv <- sample(0:1,1000,replace=T)
##'
##' # We can then estimate log(Z) using the `"BIC"` method
##' lbe.gen(method = "BIC",obs_mat = CM,res_vec = rv)
##'
##' # If we assume the prior for the regression coefficients is a standard
##' # normal
##' pr_samp <- function(p_n,di){return(rmvn(p_n,rep(0,di),diag(di)))}
##' pr_fun <- function(th,di){return(dmvn(th,mu=rep(0,di),sigma=diag(di)))}
##'
##' # If we suspect that the data might be linearly separable
##' lbe.gen(method = "SMC_BIC",obs_mat = CM,res_vec = rv,p_num = 500, rpri = pr_samp, p_pdf = pr_fun)
##'
##' # For a much smaller dataset
##' CM.2 <- matrix(rnorm(20),10,2)
##' rv.2 <- sample(0:1,10,replace=T)
##'
##' # it is more preferable to use method `"SMC"`
##' lbe.gen(method = "SMC",obs_mat = CM.2,res_vec = rv.2,p_num = 500, rpri = pr_samp, p_pdf = pr_fun)
##'
##' # For small datasets the estimates of log(Z) using BIC or SMC will differ
##' # significantly. This becomes relatively less problematic as the dataset
##' # becomes larger however.
##' sapply(10:11,FUN = function(u,covm,resv,ps,pf){lbe.gen(method="SMC_BIC",thres=u,obs_mat=covm,res_vec=resv,p_num=500,rpri=ps,p_pdf=pf)},covm = CM.2,resv = rv.2,ps = pr_samp,pf = pr_fun)
##' sapply(1000:1001,FUN = function(u,covm,resv,ps,pf){lbe.gen(method="SMC_BIC",thres=u,obs_mat=covm,res_vec=resv,p_num=500,rpri=ps,p_pdf=pf)},covm = CM,resv = rv,ps = pr_samp,pf = pr_fun)
##'
##'@seealso [MIBIS.Z]

lbe.gen <- function(method="SMC_BIC",thres=Inf,obs_mat,res_vec,p_num=0,rpri=NULL,p_pdf=NULL){
  opts <- c("SMC","SMC_BIC","BIC")
  if(!(method%in%opts)){
    stop("Method is not supported")
  }
  if(method!="BIC" & (is.null(rpri) | is.null(p_pdf))){
    stop("Both sampling function and probability density function of the prior are required for this method")
  }
  if(method=="BIC"){
    return(-BIC(glm(res_vec~.,family="binomial",data=data.frame(obs_mat,res_vec)))/2)
  } else{
    logZ <- T
    if(method=="SMC_BIC" & length(res_vec)>=thres){
      logZ <- tryCatch(-BIC(glm(res_vec~.,family="binomial",data=data.frame(obs_mat,res_vec)))/2,warning=function(...)T)
    }
    if(logZ==T){
      logZ <- MIBIS.Z(X = cbind(rep(1,nrow(obs_mat)),obs_mat), y = res_vec, sampl = rpri(p_n = p_num,di = (ncol(obs_mat)+1)),prior_pdf = p_pdf)[[3]]
    }
    return(logZ)
  }
}

###############################################
## R script for functions in UNCOVER package ##
###############################################
##
## Sam Emerson
## January 2023
##


#################################################################################
## Function for Sequential Monte Carlo Sampler to Obtain the Bayesian Evidence ##
#################################################################################



##' Logistic regression iterated batch importance sampling
##'
##'
##' @export
##' @name IBIS.logreg
##' @description This function uses an iterated batch importance sampling scheme
##' with batch size one to go from prior to full posterior. We
##' assume a Bayesian logistic regression model.
##'
##' @keywords sequential monte carlo
##' @param X Covariance matrix
##' @param y Binary response vector
##' @param prior Name of prior. Can be one of; `"mvn"` for multivariate normal,
##' `"mvl"` for multivariate laplace or `"mviu"` for multivariate independent
##' uniform. See details. Defaults to `"mvn"`.
##' @param N Number of prior samples. Defaults to 1000.
##' @param ess Threshold: if the effective sample size of the particle weights
##' falls below this value then a resample move step is triggered. Defaults to
##' `N/2`.
##' @param n_move Number of Metropolis-Hastings steps to apply each time a
##' resample move step is triggered. Defaults to 1.
##' @param weighted Should the outputted samples be weighted? As a default the
##' samples are not weighted.
##' @param ... Arguments for each of the priors. See details.
##' @return A list consisting of posterior samples and the
##' log Bayesian evidence of the full posterior. If `weighted==TRUE` the weights
##' of the samples is also provided.
##' @details Details of the internal mechanisms of the SMC sampler such as the
##' Metropolis-Hastings MCMC resample move can be found in *UNCOVER paper* and
##' Chopin (2002).
##'
##' When selecting a prior the arguments for said prior must be specified or
##' default values will be chosen. For a multivariate normal prior, mean `mu`
##' and covariance matrix `sigma` must be specified or `mu = rep(0,ncol(X)+1)`,
##' `sigma = diag(ncol(X)+1)` will be selected as default values. For a
##' multivariate Laplace prior, mean `mu` and covariance matrix `Sigma` must be
##' specified or `mu = rep(0,ncol(X)+1)`, `Sigma = diag(ncol(X)+1)` will be
##' selected as default values. For a multivariate independent uniform prior
##' (i.e. regression coefficients are independent under the prior), a minimum
##' value vector `a` and a maximum value vector `b` must be specified or
##' `a = rep(0,ncol(X)+1)`, `b = rep(1,ncol(X)+1)` will be selected as default
##' values.
##'
##' Note that decreasing `ess` and increasing `n_move` will lead to a more
##' accurate estimate of the Bayesian evidence, but at the cost of increased
##' computational time.
##'
##' @references Chopin, N. (2002). A sequential particle filter method for
##' static models. Biometrika, 89(3), 539-552.
##' @examples
##'
##' # First we generate a covariate matrix X and binary response vector y
##' CM <- matrix(rnorm(200),100,2)
##' rv <- sample(0:1,100,replace=T)
##'
##' # Now we can obtain 1000 samples from the posterior from a standard
##' # multivariate normal prior
##' out.1 <- IBIS.logreg(X = CM,y = rv)
##' pairs(out.1$samples)
##' out.1$log_Bayesian_evidence
##'
##' # We can specify that the samples be weighted
##' out.1.w <- IBIS.logreg(X = CM,y = rv,options = IBIS.logreg.opts(weighted = TRUE))
##' weight_fade <- colorRampPalette(c('white','black'))(100)[as.numeric(cut(out.1.w$weights,breaks = 100))]
##' pairs(out.1.w$samples,col = weight_fade)
##'
##' # We can also specify different arguments for a specific prior
##' out.2 <- IBIS.logreg(X = CM,y = rv,prior_mean = rep(-3,3),prior_var = 0.1*diag(3))
##' pairs(rbind(out.1$samples,out.2$samples),col=rep(c(1,2),each=1000))
##' out.2$log_Bayesian_evidence
##' out.3 <- IBIS.logreg(X = CM,y = rv,prior_mean = rep(3,3),prior_var = 0.1*diag(3))
##' pairs(rbind(out.1$samples,out.2$samples,out.3$samples),col=rep(c(1,2,3),each=1000))
##' out.3$log_Bayesian_evidence
##'
##' # We can also change the prior, for example a multivariate independent
##' # uniform
##' rmviu <- function(n,a,b){
##' return(mapply(FUN = function(min.vec,max.vec,pn){runif(pn,a,b)},min.vec=a,max.vec=b,MoreArgs = list(pn = n)))
##' }
##' dmviu <- function(x,a,b){
##' for(ii in 1:ncol(x)){
##'   x[,ii] <- dunif(x[,ii],a[ii],b[ii])
##' }
##' return(apply(x,1,prod))
##' }
##'
##' out.4 <- IBIS.logreg(X = CM,y = rv,options = IBIS.logreg.opts(prior.override = TRUE,rprior = rmviu,dprior = dmviu,a=rep(0,3),b=rep(1,3)))
##' pairs(rbind(out.1$samples,out.4$samples),col=rep(c(1,4),each=1000))
##' out.4$log_Bayesian_evidence
##'

IBIS.logreg <- function(X,y,options = IBIS.logreg.opts(),
                        prior_mean = rep(0,ncol(X)+1),prior_var = diag(ncol(X)+1)){
  if(options$prior.override){
    MoreArgs <- options$MoreArgs
    TotArgs <- list(options$N)
    if(length(MoreArgs)!=0){
      TotArgs <- c(TotArgs,MoreArgs)
    }
    rprior <- options$rprior
    test_rprior <- do.call(rprior,TotArgs)
    if(nrow(test_rprior)!=options$N | ncol(test_rprior)!=(ncol(X)+1)){
      stop("rprior specified does not produce an N by (ncol(X)+1) matrix")
    }
    dprior <- options$dprior
    TotArgs <- list(test_rprior)
    if(length(MoreArgs)!=0){
      TotArgs <- c(TotArgs,MoreArgs)
    }
    test_dprior <- do.call(dprior,TotArgs)
    if(length(test_dprior)!=options$N){
      stop("dprior specified does not produce a vector of densities for the
           N samples provided")
    }
  } else{
    rprior <- mvnfast::rmvn
    dprior <- mvnfast::dmvn
    MoreArgs <- list(mu = prior_mean,sigma = prior_var)
  }
  DM <- cbind(rep(1,length(y)),X)
  IBIS_out <- IBIS.Z(X = DM,y = y,rprior = rprior,N = options$N,prior_pdf = dprior,
                     ess = options$ess,n_move = options$n_move,
                     PriorArgs = MoreArgs)$output
  if(options$weighted==FALSE){
    if(length(unique(IBIS_out$weights))!=1){
      ss <- cov.wt(IBIS_out$samples,wt=IBIS_out$weights,method = "ML")
      mu <- ss$center
      Sigma <- ss$cov
      samp <- sample(as.integer(names(IBIS_out$duplication_table)),options$N,prob = IBIS_out$weights[as.integer(names(IBIS_out$duplication_table))]*IBIS_out$duplication_table,replace = TRUE)
      IBIS_out$samples <- IBIS_out$samples[samp,]
      for(j in 1:options$n_move){
        BC <- mvnfast::rmvn(options$N,mu = mu, sigma=Sigma)
        Log_1 <- log(1+exp(DM%*%t(BC)))*(y-1) + log(1+exp(-DM%*%t(BC)))*(-y)
        TotArgs <- list(BC)
        if(length(MoreArgs)!=0){
          TotArgs <- c(TotArgs,MoreArgs)
        }
        Log_1 <- colSums(Log_1) + log(do.call(dprior,TotArgs)) + log(mvnfast::dmvn(IBIS_out$samples,mu=mu,sigma=Sigma))
        Log_2 <- log(1+exp(DM%*%t(IBIS_out$samples)))*(y-1) + log(1+exp(-DM%*%t(IBIS_out$samples)))*(-y)
        TotArgs <- list(IBIS_out$samples)
        if(length(MoreArgs)!=0){
          TotArgs <- c(TotArgs,MoreArgs)
        }
        Log_2 <- colSums(Log_2) + log(do.call(dprior,TotArgs)) + log(mvnfast::dmvn(BC,mu=mu,sigma=Sigma))
        A <- Log_1 - Log_2 - log(runif(length(Log_1)))
        IBIS_out$samples[which(A>0),] <- BC[which(A>0),]
      }
    }
    res <- list(covariance_matrix = data.frame(X),
                response_vector = y,
                samples = IBIS_out$samples,
                log_Bayesian_evidence = IBIS_out$log_Bayesian_evidence)
  } else{
    res <- list(covariance_matrix = data.frame(X),
                response_vector = y,
                samples = IBIS_out$samples,
                weights = IBIS_out$weights,
                log_Bayesian_evidence = IBIS_out$log_Bayesian_evidence)
  }
  class(res) <- 'IBIS'
  res
}

print.IBIS <- function(x){
  if(length(x)==4){
    xx <- colMeans(x$samples)
    names(xx) <- c("(Intercept)",colnames(x$covariance_matrix))
    cat(nrow(x$samples),"posterior samples with mean:\n")
    cat("\n")
    print.default(xx)
  } else{
    ss <- cov.wt(x$samples,wt=x$weights,method = "ML")
    xx <- ss$center
    names(xx) <- c("(Intercept)",colnames(x$covariance_matrix))
    cat(nrow(x$samples),"posterior samples with weighted mean:\n")
    cat("\n")
    print.default(xx)
  }
  cat("\n")
  cat("Log Bayesian Evidence",x$log_Bayesian_evidence)
}


##' Log Bayesian evidence generator
##'
##'
##' @export
##' @name IBIS.logreg.opts
##' @description Estimates the log Bayesian evidence (log(Z)) of a Bayesian
##' logistic regression model.
##'
##' @keywords sequential monte carlo
##' @param obs_mat Covariate matrix
##' @param res_vec Binary response vector
##' @param obs_ind Vector of indices indicating which observations are in the
##' model. Defaults to all observations.
##' @param thres The threshold for which the number of observations needs to
##' exceed to consider using BIC as an estimator. Defaults to 30 if not
##' specified.
##' @param memo_thres_bic The threshold for when it is deemed worthwhile to
##' check the cache of function `memo.bic` for similar observation indices.
##' Defaults to never checking the cache.
##' @param memo_thres_smc The threshold for when it is deemed worthwhile to
##' check the cache of function `IBIS.Z` for similar observation indices.
##' Defaults to never checking the cache.
##' @param p_num Number of samples of the prior used for the SMC sampler.
##' Default value is 1000 samples.
##' @param rpri Function to sample from the prior. Must only have two arguments,
##' `p_n` and `di` (Number of prior samples to generate and the number of
##' dimensions of a single sample respectively).
##' @param p_pdf Probability Density Function of the prior. Must only have two
##' arguments, `th` and `di` (a vector or matrix of regression coefficients
##' samples and the number of dimensions of a single sample respectively).
##' @param efs Threshold: if the effective sample size of the particle weights
##' falls below this value then a resample move step is triggered. Defaults to
##' `p_num/2`. See `IBIS.Z` for details.
##' @param nm Number of Metropolis-Hastings steps to apply each time a
##' resample move step is triggered. Defaults to 1. See `IBIS.Z` for details.
##' @return An estimation of the log Bayesian evidence
##' @details log(Z) is estimated using three possible methods:
##'
##' 1. `"SMC"`: Estimates log(Z) using an SMC sampler. In order to achieve this
##' set `thres = Inf`. See `IBIS.Z` for more details.
##'
##' 2. `"SMC_BIC"`: Initially tries to estimate log(Z) using BIC but will revert
##' to using an SMC sampler if the data is linearly separable. Also reverts to
##' using an SMC sampler if the number of observations is below a certain
##' threshold. Recommended method and is set as default.
##'
##' 3. `"BIC"`: Always tries to estimate log(Z) using the
##' Bayesian Information Criterion (BIC). Only reverts to an SMC sampler if the
##' data is linearly separable. In order to achieve this set `thres = 0`. Only
##' advisable if the number of observations is large.
##'
##' `memo_thres_bic` is recommended to be high as the computation time spent
##' searching for a similar previous run of `memo.bic` is likely to be much
##' higher than running `memo.bic` from scratch when the number of observations
##' isn't extremely high. `memo_thres_smc` should be specified much lower as
##' the running the SMC sampler function is much more costly.
##'
##' @examples
##'
##' # First we generate a covariate matrix `obs_mat` and binary response vector
##' # `res_vec`
##' CM <- matrix(rnorm(2000),1000,2)
##' rv <- sample(0:1,1000,replace=T)
##'
##' # If we assume the prior for the regression coefficients is a standard
##' # normal
##' pr_samp <- function(p_n,di){return(mvnfast::rmvn(p_n,rep(0,di),diag(di)))}
##' pr_fun <- function(th,di){return(mvnfast::dmvn(th,mu=rep(0,di),sigma=diag(di)))}
##'
##' # We can then estimate log(Z) using the `"BIC"` method
##' lbe.gen(obs_mat = CM,res_vec = rv,thres = 0,rpri = pr_samp, p_pdf = pr_fun)
##'
##' # If we suspect that the data might be linearly separable
##' get("_cache", envir=environment(memo.bic))$reset()
##' get("_cache", envir=environment(IBIS.Z))$reset()
##' lbe.gen(obs_mat = CM,res_vec = rv,rpri = pr_samp, p_pdf = pr_fun)
##'
##' # For a much smaller dataset
##' CM.2 <- matrix(rnorm(20),10,2)
##' rv.2 <- sample(0:1,10,replace=T)
##'
##' # it is more preferable to use only the SMC sampler to estimate log(Z)
##' get("_cache", envir=environment(memo.bic))$reset()
##' get("_cache", envir=environment(IBIS.Z))$reset()
##' lbe.gen(obs_mat = CM.2,res_vec = rv.2,thres=Inf,rpri = pr_samp, p_pdf = pr_fun)
##'
##' # For small datasets the estimates of log(Z) using BIC or SMC will differ
##' # significantly. This becomes relatively less problematic as the dataset
##' # becomes larger however.
##' get("_cache", envir=environment(memo.bic))$reset()
##' get("_cache", envir=environment(IBIS.Z))$reset()
##' small.bic <- lbe.gen(thres=0,obs_mat = CM.2,res_vec = rv.2,rpri = pr_samp, p_pdf = pr_fun)
##' get("_cache", envir=environment(memo.bic))$reset()
##' get("_cache", envir=environment(IBIS.Z))$reset()
##' small.smc <- lbe.gen(thres=Inf,obs_mat = CM.2,res_vec = rv.2,rpri = pr_samp, p_pdf = pr_fun)
##'c(small.bic,small.smc)
##'
##' get("_cache", envir=environment(memo.bic))$reset()
##' get("_cache", envir=environment(IBIS.Z))$reset()
##' large.bic <- lbe.gen(thres=0,obs_mat = CM,res_vec = rv,rpri = pr_samp, p_pdf = pr_fun)
##' get("_cache", envir=environment(memo.bic))$reset()
##' get("_cache", envir=environment(IBIS.Z))$reset()
##' large.smc <- lbe.gen(thres=Inf,obs_mat = CM,res_vec = rv,rpri = pr_samp, p_pdf = pr_fun)
##' c(large.bic,large.smc)
##'
##'@seealso [IBIS.Z]

IBIS.logreg.opts <- function(N=1000,ess = N/2,n_move = 1,weighted = FALSE,
                             prior.override = FALSE,rprior = NULL,
                             dprior = NULL,...){
  if(prior.override){
    if(is.null(rprior) | is.null(dprior)){
      stop("If overriding the default prior rprior and dprior must be specified.")
    }
    MoreArgs <- list(...)
  } else{
    MoreArgs = NULL
  }
  if(N <= 1){
    stop("N must be greater than 1")
  }
  if(ess > N | ess < 0){
    stop("Effective sample size must be between 0 and N")
  }
  if(n_move < 1){
    stop("n_move must be an integer greater or equal than 1")
  }
  if(!is.logical(weighted)){
    stop("weighted must be logical")
  }
  return(list(N = 1000,ess = N/2,n_move = 1,rprior = rprior,dprior = dprior,
              prior.override = prior.override,weighted = weighted,
              MoreArgs = MoreArgs))
}


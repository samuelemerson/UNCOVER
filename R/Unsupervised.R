###############################################
## R script for functions in UNCOVER package ##
###############################################
##
## Sam Emerson
## September 2022
##


##############################################################
## Functions for Unsupervised 'Cluster then Predict' Models ##
##############################################################


##' Cluster observations through K-means and then calculate the Bayesian evidence of the resulting multiplicative Bayesian logistic regression posterior.
##'
##'
##' @export
##' @name k.means.model
##' @description This function first applies K-means clustering (with K being specified by the user) on soley the covaraite data, then a separate Bayesian logistic regression model is assumed for each cluster (with the overall model being the product of these sub-models) and the log Bayesian evidence of each of these models is calculated.
##'
##' Used as a standard method to compare UNCOVER against.
##'
##' @keywords k-means cluster-then-predict
##' @param X Covariate matrix
##' @param y Binary response vector
##' @param K Number of clusters required
##' @param SMC_method Method to be used to estimate the log Bayesian evidence of the newly formed sub-models, can be one of `"SMC"`,`"SMC_BIC"` or `"BIC"`
##' @param SMC_thres The threshold for which the number of observations needs to exceed to consider using BIC as an estimator. Only applies if method `"SMC_BIC"` is selected. Defaults to 30 if not specified.
##' @param N Number of samples of the prior used for the SMC sampler. Not required if method `"BIC"` selected.
##' @param rprior Function to sample from the prior. Must only have two arguments, `p_num` and `di` (Number of prior samples to generate and the number of dimensions of a single sample respectively).
##' @param prior_pdf Probability Density Function of the prior. Must only have two arguments, `th` and `di` (a vector or matrix of regression coefficients samples and the number of dimensions of a single sample respectively).
##' @return A list consisting of; the cluster allocation vector representing the partition of the data through K-means and the Bayesian evidence vector for the resulting sub-models.
##' @details Assumes a Bayesian logistic regression model for each cohort, with the overall model being a product of these sub-models.
##'
##' The number of clusters must be specified, and if unknown running this function multiple times with differing values of `K` and comparing the overall Bayesian evidence of the outputs is recommended.
##'
##' If `rprior` and `prior_pdf` are not specified then the default prior is a standard multivariate normal.
##'
##' For more details on the specifics of the possible values for `SMC_method`, see the help page of the function `lbe.gen`.
##'
##' @seealso [lbe.gen,UNCOVER]
##' @examples
##'
##' # First we generate a covariate matrix `X` and binary response vector `y`
##' CM <- matrix(rnorm(200),100,2)
##' rv <- sample(0:1,100,replace=T)
##'
##' # We can then run our algorithm to see what cohorts are selected for
##' # differing values of `K`
##' k.means.list <- lapply(2:5, FUN = function(u,obs,res){k.means.model(X = obs,y = res,K = u)},obs = CM,res = rv)
##'
##' # We can also discover which `K` results in the highest Bayesian evidence
##' K_opt <- (2:5)[which.max(lapply(k.means.list, FUN = function(u){sum(u[[2]])}))]
##' K_opt
##' k.means.opt <- k.means.list[[which((2:5)==K_opt)]]
##'
##' # If we don't assume the prior for the regression coefficients is a
##' # standard normal but instead a multivariate normal with mean (1,1) and the
##' # identity matrix as the covariance matrix we can specify
##' pr_samp <- function(p_n,di){return(rmvn(p_n,rep(1,di),diag(di)))}
##' pr_fun <- function(th,di){return(dmvn(th,mu=rep(1,di),sigma=diag(di)))}
##'
##' # We then can run `k.means.model` using this prior and compare to the
##' # standard result
##' k.means.opt.2 <- k.means.model(X = CM,y = rv,K = K_opt,rprior = pr_samp,prior_pdf = pr_fun)
##' c(sum(k.means.opt[[2]]),sum(k.means.opt.2[[2]]))

k.means.model <- function(X,y,K,SMC_method = "SMC_BIC",SMC_thres=30,N=1000,rprior=NULL,prior_pdf=NULL){
  if((!is.null(rprior)&is.null(prior_pdf))|(is.null(rprior)&!is.null(prior_pdf))){
    stop("Both sampling function and probability density function of the prior are required")
  }
  if(is.null(rprior)&is.null(prior_pdf)){
    rprior <- function(p_num,di){
      return(rmvn(p_num,rep(0,di),diag(di)))
    }
    prior_pdf <- function(th,di){
      return(dmvn(th,mu=rep(0,di),sigma=diag(di)))
    }
  }
  k_means_un <- kmeans(X,K)
  z <- k_means_un$cluster
  logZ <- rep(Inf,K)
  for(k in 1:K){
    logZ[k] <- lbe.gen(method = SMC_method,thres = SMC_thres,obs_mat = X[which(z==k),,drop=F],res_vec = y[which(z==k)],p_num = N,rpri = rprior,p_pdf = prior_pdf)
  }
  return(list("Cluster Allocation" = z,"Log Marginal Likelihoods" = logZ))
}

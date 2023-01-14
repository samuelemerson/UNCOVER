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
##' @name IBIS.Z
##' @description This function uses an iterated batch importance sampling scheme
##' with batch size one to go from one bridging distribution to another. We
##' assume a Bayesian logistic regression model.
##'
##' There are two initialisations for this function (each requiring different
##' inputs). You may either; start from a prior or start from a bridging
##' distribution.
##'
##' Used in UNCOVER to generate the log sub-Bayesian evidence of partitioned
##' models, however the weighted samples can also be used for posterior inference.
##'
##'
##' @keywords sequential monte carlo
##' @param X Design matrix
##' @param y Binary response vector
##' @param sampl A named list containing; `samples` - A matrix of samples from
##' the starting partial posterior, `weights` - the samples associated weights,
##' `log_Bayesian_evidence` - the log Bayesian evidence of the partial posterior
##' and `duplication_table` - a table where the names in the table refer to the
##' indices of the unique samples of the partial posterior and the elements of
##' the table indicate how many duplicates of this unique sample there are. If
##' starting from the prior this argument should be ignored.
##' @param rprior Function to sample from the prior. Must only have two
##' arguments, `p_num` and `di` (Number of prior samples to generate and the
##' number of dimensions of a single sample respectively). If starting from a
##' partial posterior this argument should be ignored.
##' @param N Number of prior samples. If starting from a partial posterior this
##' argument should be ignored.
##' @param prior_pdf Probability Density Function of the prior. Must only have
##' two arguments, `th` and `di` (a vector or matrix of regression coefficients
##' samples and the number of dimensions of a single sample respectively).
##' @param target_set Vector of observation indices for the outputted
##' distribution. If the desired final distribution is the full
##' posterior this should not be specified and will be set as `1:length(y)`.
##' @param current_set Vector of observation indices already added to give the
##' starting bridging distribution. If the starting 'bridging' distribution is
##' the prior then this should not be specified.
##' @param ess Threshold: if the effective sample size of the particle weights
##' falls below this value then a resample move step is triggered. Defaults to
##' `N/2`.
##' @param n_move Number of Metropolis-Hastings steps to apply each time a
##' resample move step is triggered. Defaults to 1.
##' @return A list consisting of output and input. Output is a list consisting
##' of; weighted posterior samples (the samples and
##' weights are given in separate lists), the log Bayesian evidence of the final
##' partial posterior (the full posterior if `target_set` is not specified) and the
##' duplication table of the outputted posterior samples. Input is the inputted
##' `target_set`.
##' @details The setting going from prior to full posterior can be achieved
##' through specification of only `X`, `y`, `rprior`, `N` and `prior_pdf`. If
##' using samples from a partial posterior instead of prior samples then
##' `rprior` and `N` can be ignored but `sampl` and `current_set` should be
##' provided. If the samples from the partial posterior are not weighted samples
##' then the `weights` element of the `sampl` list should be
##' `rep(1,nrow(sampl$samples))`. If the samples provided in `sampl` are all
##' unique then the `duplication_table` element of `sampl` should be a table of
##' `sample$weights = rep(1,nrow(sampl$samples))`. If the full posterior is not
##' the desired output then the specific bridging distribution (partial
##' posterior) required must be specified (i.e.`target_set` is required).
##'
##' Note that it is possible to specify `target_set` as a subset of
##' `current_set`. If this is selected then reverse iterated batch importance
##' sampling will be applied through removal of observations in `current_set`.
##'
##' Details of the internal mechanisms of the SMC sampler such as the
##' Metropolis-Hastings MCMC resample move can be found in *UNCOVER paper* and
##' Chopin (2002).
##'
##' Note that decreasing `ess` and increasing `n_move` will lead to a more
##' accurate estimate of the Bayesian evidence, but at the cost of increased
##' computational time.
##'
##' `IBIS.Z` is to be used as a memoised function within `lbe.gen` and as such
##' requires specification of the input as an output when accessing the cache.
##' @references Chopin, N. (2002). A sequential particle filter method for
##' static models. Biometrika, 89(3), 539-552.
##' @examples
##'
##' # First we generate a design matrix X and binary response vector y
##' DM <- cbind(rep(1,100),matrix(rnorm(200),100,2))
##' rv <- sample(0:1,100,replace=T)
##'
##' # We assume the prior for the regression coefficients is a standard normal
##' pr_samp <- function(p_num,di){return(rmvn(p_num,rep(0,di),diag(di)))}
##' pr_fun <- function(th,di){return(dmvn(th,mu=rep(0,di),sigma=diag(di)))}
##'
##' # Now we can obtain 1000 samples from the posterior using only 1000 samples
##' # from the prior
##' out.1 <- IBIS.Z(X = DM,y = rv,rprior = pr_samp,N = 1000,prior_pdf = pr_fun)
##' unique_ind <- as.numeric(names(out.1$output$duplication_table))
##' unique_samp <- out.1$output$samples[unique_ind,]
##' unique_weights <- out.1$output$weights[unique_ind]
##' weight_fade <- colorRampPalette(c('white','black'))(100)[as.numeric(cut(unique_weights,breaks = 100))]
##' pairs(unique_samp,col = weight_fade)
##' out.1$output$log_Bayesian_evidence
##'
##' # If we then added 100 more observations to our data set
##' DM.2 <- rbind(DM,cbind(rep(1,100),matrix(rnorm(200),100,2)))
##' rv.2 <- c(rv,sample(0:1,100,replace=T))
##'
##' # and then wanted samples from the entire posterior instead of starting
##' # from the prior we can use the previous SMC sampler output
##' out.2 <- IBIS.Z(X = DM.2,y = rv.2,sampl = out.1$output,prior_pdf = pr_fun,current_set = 1:100)
##' unique_ind.2 <- as.numeric(names(out.2$output$duplication_table))
##' unique_samp.2 <- out.2$output$samples[unique_ind.2,]
##' unique_weights.2 <- out.2$output$weights[unique_ind.2]
##' weight_fade.2 <- colorRampPalette(c('white','black'))(100)[as.numeric(cut(unique_weights.2,breaks = 100))]
##' pairs(unique_samp.2,col = weight_fade.2)
##' out.2$output$log_Bayesian_evidence
##'
##' # we can also observe roughly how the particles move when the 100 new
##' # observations are added in two stages
##' out.3 <- IBIS.Z(X = DM.2,y = rv.2,sampl = out.1$output,prior_pdf = pr_fun,target_set = 1:150,current_set = 1:100)
##' unique_ind.3 <- as.numeric(names(out.3$output$duplication_table))
##' unique_samp.3 <- out.3$output$samples[unique_ind.3,]
##' unique_weights.3 <- out.3$output$weights[unique_ind.3]
##' weight_fade.3 <- colorRampPalette(c('white','red'))(100)[as.numeric(cut(unique_weights.3,breaks = 100))]
##' out.4 <- IBIS.Z(X = DM.2,y = rv.2,sampl = out.3$output,prior_pdf = pr_fun,current_set = 1:150)
##' unique_ind.4 <- as.numeric(names(out.4$output$duplication_table))
##' unique_samp.4 <- out.4$output$samples[unique_ind.4,]
##' unique_weights.4 <- out.4$output$weights[unique_ind.4]
##' weight_fade.4 <- colorRampPalette(c('white','green'))(100)[as.numeric(cut(unique_weights.4,breaks = 100))]
##' pairs(rbind(unique_samp,unique_samp.3,unique_samp.4),col=c(weight_fade,weight_fade.3,weight_fade.4))
##'
##' # If we wished to remove 25 observations from our original posterior
##' out.5 <- IBIS.Z(X = DM,y = rv,sampl = out.1$output,prior_pdf = pr_fun,target_set=1:75,current_set = 1:100)
##' unique_ind.5 <- as.numeric(names(out.5$output$duplication_table))
##' unique_samp.5 <- out.5$output$samples[unique_ind.5,]
##' unique_weights.5 <- out.5$output$weights[unique_ind.5]
##' weight_fade.5 <- colorRampPalette(c('white','black'))(100)[as.numeric(cut(unique_weights.5,breaks = 100))]
##' pairs(unique_samp.5,col = weight_fade.5)
##' out.5$output$log_Bayesian_evidence
##'
##' # we can also compare this method to standard iterated batch importance
##' # sampling
##' out.6 <- IBIS.Z(X = DM,y = rv,rprior = pr_samp,N = 1000,prior_pdf = pr_fun,target_set = 1:50)
##' unique_ind.6 <- as.numeric(names(out.6$output$duplication_table))
##' unique_samp.6 <- out.6$output$samples[unique_ind.6,]
##' unique_weights.6 <- out.6$output$weights[unique_ind.6]
##' weight_fade.6 <- colorRampPalette(c('white','red'))(100)[as.numeric(cut(unique_weights.6,breaks = 100))]
##' out.7 <- IBIS.Z(X = DM,y = rv,sampl = out.1$output,prior_pdf = pr_fun,target_set = 1:50,current_set = 1:100)
##' unique_ind.7 <- as.numeric(names(out.7$output$duplication_table))
##' unique_samp.7 <- out.7$output$samples[unique_ind.7,]
##' unique_weights.7 <- out.7$output$weights[unique_ind.7]
##' weight_fade.7 <- colorRampPalette(c('white','green'))(100)[as.numeric(cut(unique_weights.7,breaks = 100))]
##' pairs(rbind(unique_samp.6,unique_samp.7),col=c(weight_fade.6,weight_fade.7))
##' out.6$output$log_Bayesian_evidence
##' out.7$output$log_Bayesian_evidence
##'
##' @seealso [lbe.gen]


IBIS.Z <- function(X,y,sampl=NULL,rprior=NULL,N=NULL,prior_pdf,
                   target_set=1:length(y),current_set=NULL,ess=NULL,n_move=1){
  if(length(target_set)==length(current_set)){
    stop("target_set and current_set cannot be the same length")
  }
  if(length(target_set)>length(current_set)){
    reverse <- FALSE
    if(length(current_set)!=0){
      if(sum(sort(target_set)[1:length(current_set)]-sort(current_set))!=0){
        stop("For standard IBIS current_set must be a subset of target_set")
      }
    }
  } else{
    reverse <- TRUE
    ess_point <- length(current_set)
    if(length(target_set)!=0){
      if(sum(sort(current_set)[1:length(target_set)]-sort(target_set))!=0){
        stop("For reverse IBIS target_set must be a subset of current_set")
      }
    }
  }
  if(is.null(sampl) & (is.null(rprior)|is.null(N))){
    stop("Either specify weighted samples and their current log Bayesian evidence or a prior function to generate samples along with the number of samples")
  }
  if(!is.null(sampl) & is.null(current_set)){
    stop("If sampl is specified you must indicate which observations have been added with current_set. If initialising with the prior please use rprior and N instead of sampl.")
  }
  X <- as.matrix(X)
  p <- ncol(X)
  if(is.null(sampl)) {
    sampl <- rprior(N,p)
    w <- rep(1,N)
    logZ <- 0
    dup_tab <- table(1:N)
  } else {
    N <- nrow(sampl$samples)
    w <- sampl$weights
    logZ <- sampl$log_Bayesian_evidence
    dup_tab <- sampl$duplication_table
    sampl <- sampl$samples
  }
  if(is.null(ess)){
    ess <- N/2
  }
  if(reverse){
    action_set <- setdiff(current_set,target_set)
  } else{
    action_set <- setdiff(target_set,current_set)
  }
  action_set <- action_set[sample(1:length(action_set),length(action_set))]
  for(i in 1:length(action_set)){
    w_new <- as.vector(((1+exp(-X[action_set[i],]%*%t(sampl)))^(-y[action_set[i]]))*((1+exp(X[action_set[i],]%*%t(sampl)))^(y[action_set[i]]-1)))
    if(reverse){
      w_new <- w/w_new
    } else{
      w_new <- w*w_new
    }
    logZ <- logZ + log(sum(w_new)) - log(sum(w))
    w <- w_new
    wsq <- sum((w[as.integer(names(dup_tab))]*dup_tab)^2)
    if(wsq==0){
      ESS <- 0
    } else{
      ESS <- (sum(w)^2)/wsq
    }
    if(ESS < ess){
      if(reverse){
        if(length(current_set)==1){
          cont_rate <- ess_point
        } else{
          cont_rate <- ess_point/(length(current_set)-1)
          ess_point <- length(current_set)-1
        }
      }
      ss <- cov.wt(sampl,wt=w,method = "ML")
      mu <- ss$center
      Sigma <- ss$cov
      samp <- sample(as.integer(names(dup_tab)),N,prob = w[as.integer(names(dup_tab))]*dup_tab,replace = TRUE)
      sampl <- sampl[samp,]
      w <- rep(1,N)
      if(reverse){
        sampl <- t(((t(sampl) - mu)*sqrt(cont_rate)) + mu)
      }
      A.all <- rep(FALSE,N)
      for(j in 1:n_move){
        BC <- mvnfast::rmvn(N,mu = mu, sigma=Sigma)
        if(reverse){
          if(length(setdiff(current_set,action_set[i]))==0){
            Log_1 <- log(prior_pdf(th = BC,di = p)) + log(mvnfast::dmvn(sampl,mu=mu,sigma=Sigma))
            Log_2 <- log(prior_pdf(th = sampl,di = p)) + log(mvnfast::dmvn(BC,mu=mu,sigma=Sigma))
          } else{
            Log_1 <- log(1+exp(X[setdiff(current_set,action_set[i]),]%*%t(BC)))*(y[setdiff(current_set,action_set[i])]-1) + log(1+exp(-X[setdiff(current_set,action_set[i]),]%*%t(BC)))*(-y[setdiff(current_set,action_set[i])])
            Log_1 <- colSums(Log_1) + log(prior_pdf(th = BC,di = p)) + log(mvnfast::dmvn(sampl,mu=mu,sigma=Sigma))
            Log_2 <- log(1+exp(X[setdiff(current_set,action_set[i]),]%*%t(sampl)))*(y[setdiff(current_set,action_set[i])]-1) + log(1+exp(-X[setdiff(current_set,action_set[i]),]%*%t(sampl)))*(-y[setdiff(current_set,action_set[i])])
            Log_2 <- colSums(Log_2) + log(prior_pdf(th = sampl,di = p)) + log(mvnfast::dmvn(BC,mu=mu,sigma=Sigma))
          }
        } else{
          Log_1 <- log(1+exp(X[c(current_set,action_set[i]),]%*%t(BC)))*(y[c(current_set,action_set[i])]-1) + log(1+exp(-X[c(current_set,action_set[i]),]%*%t(BC)))*(-y[c(current_set,action_set[i])])
          Log_1 <- colSums(Log_1) + log(prior_pdf(th = BC,di = p)) + log(mvnfast::dmvn(sampl,mu=mu,sigma=Sigma))
          Log_2 <- log(1+exp(X[c(current_set,action_set[i]),]%*%t(sampl)))*(y[c(current_set,action_set[i])]-1) + log(1+exp(-X[c(current_set,action_set[i]),]%*%t(sampl)))*(-y[c(current_set,action_set[i])])
          Log_2 <- colSums(Log_2) + log(prior_pdf(th = sampl,di = p)) + log(mvnfast::dmvn(BC,mu=mu,sigma=Sigma))
        }
        A <- Log_1 - Log_2 - log(runif(length(Log_1)))
        A.all <- A.all | (A>0)
        sampl[which(A>0),] <- BC[which(A>0),]
      }
      samp[which(A.all)] <- (N+1):(N+length(which(A.all)))
      dup_tab <- table(samp)
      names(dup_tab) <- as.character(match(as.integer(names(dup_tab)),samp))
    }
    if(reverse){
      current_set <- setdiff(current_set,action_set[i])
    } else{
      current_set <- c(current_set,action_set[i])
    }
  }
  return(list(output = list(samples = sampl,
                            weights = w,
                            log_Bayesian_evidence = logZ,
                            duplication_table = dup_tab[order(as.integer(names(dup_tab)))]),
              input = target_set))
}

##' Log Bayesian evidence estimator through BIC, built for memoisation
##'
##'
##' @export
##' @name memo.bic
##' @description Estimates the log Bayesian evidence (log(Z)) of a Bayesian
##' logistic regression model using BIC.
##'
##' @keywords memoisation BIC
##' @param X Covariate matrix
##' @param y Binary response vector
##' @param which_obs The indices of the observations in `X` that are in the
##' model. Defaults to all observations if not specified.
##' @param param_start Starting values for the model parameters used in the
##' `glm` function.
##' @return A list containing; logZ - An estimation of the log Bayesian
##' evidence, coeffs - The coefficients produced by the glm model required for
##' calculation of the BIC value and input - the input `which_obs`
##' @details This function is built with memoisation in mind and is intended to
##' be used as a supplementary function to `lbe.gen`. Producing `which_obs` as
##' an output is used to obtain input information from keys in the cache when
##' reviewing inputs which have been memoised.
##'
##' @examples
##'
##' # First we generate a covariate matrix `X` and binary response vector
##' # `y`
##' CM <- matrix(rnorm(2000000),1000000,2)
##' rv <- sample(0:1,1000000,replace=T)
##'
##' # We can then estimate log(Z) using the `"BIC"` method
##' model.1000000 <- memo.bic(X = CM,y = rv)
##' model.1000000$logZ
##' model.1000000$coeffs
##'
##' # If we are only concerned with the first 500000 observations
##' model.500000 <- memo.bic(X = CM,y = rv, which_obs = 1:500000)
##'
##' # If we are then wanting the model to include the first 500001 observations,
##' # it may be faster to use the output of model.500000
##' system.time(memo.bic(X = CM,y = rv, which_obs = 1:500001))
##' system.time(memo.bic(X = CM,y = rv, which_obs = 1:500001,param_start = model.500000$coeffs))
##'
##'@seealso [glm,BIC,lbe.gen]

memo.bic <- function(X,y,which_obs=1:length(y),param_start=NULL){
  bic_glm <- glm(y~.,family="binomial",data=data.frame(X,y)[which_obs,],start = param_start)
  lbe <- -BIC(bic_glm)/2
  return(list(logZ = lbe,coeffs = bic_glm$coefficients,input = which_obs))
}


##' Log Bayesian evidence generator
##'
##'
##' @export
##' @name lbe.gen
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
##' pr_samp <- function(p_n,di){return(rmvn(p_n,rep(0,di),diag(di)))}
##' pr_fun <- function(th,di){return(dmvn(th,mu=rep(0,di),sigma=diag(di)))}
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

lbe.gen <- function(obs_mat,res_vec,obs_ind = 1:length(res_vec),thres=30,
                    memo_thres_bic = Inf,memo_thres_smc = Inf,p_num=1000,rpri,
                    p_pdf,efs = p_num/2,nm = 1){
  cache_bic <- get("_cache", envir=environment(memo.bic))
  cache_smc <- get("_cache", envir=environment(IBIS.Z))
  cache_search_fun <- function(u,cs,oi){
    c_ind <- cs$get(u)$input
    lci <- length(c_ind)
    loi <- length(obs_ind)
    if(lci==loi){
      if(sum(c_ind - obs_ind)==0){
        return(lci-loi)
      } else{
        return(Inf)
      }
    }
    if(lci<loi){
      if(sum(c_ind - obs_ind[1:lci])==0){
        return(loi-lci)
      } else{
        return(Inf)
      }
    }
    if(lci>loi){
      if(sum(c_ind[1:loi] - obs_ind)==0){
        return(lci-loi)
      } else{
        return(Inf)
      }
    }
  }
  logZ <- T
  if(length(res_vec)>=thres){
    if(length(res_vec)<memo_thres_bic | length(cache_bic$keys())==0){
      logZ <- tryCatch(memo.bic(X = obs_mat, y = res_vec, which_obs = obs_ind)$logZ,warning=function(...)T)
    } else{
      sub_size <- sapply(cache_bic$keys(),FUN = cache_search_fun,cs = cache_bic,oi = obs_ind)
      if(min(sub_size)==0 | min(sub_size)==Inf){
        logZ <- tryCatch(memo.bic(X = obs_mat, y = res_vec, which_obs = obs_ind)$logZ,warning=function(...)T)
      } else{
        bic.start <- memo.bic(X = obs_mat, y = res_vec, which_obs = cache_bic$get(cache_bic$keys()[which.min(sub_size)])$input)$coeffs
        logZ <- tryCatch(memo.bic(X = obs_mat, y = res_vec, which_obs = obs_ind,param_start = bic.start)$logZ,warning=function(...)T)
      }
    }
  }
  if(logZ==T){
    if(length(res_vec)<memo_thres_smc | length(cache_smc$keys())==0){
      logZ <- IBIS.Z(X = cbind(rep(1,nrow(obs_mat)),obs_mat),y = res_vec,rprior = rpri,N = p_num,prior_pdf = p_pdf,target_set = obs_ind,ess = efs,n_move = nm)$output$log_Bayesian_evidence
    } else{
      sub_size <- sapply(cache_smc$keys(),FUN = cache_search_fun,cs = cache_smc,oi = obs_ind)
      if(min(sub_size)==0 | min(sub_size)==Inf){
        logZ <- IBIS.Z(X = cbind(rep(1,nrow(obs_mat)),obs_mat),y = res_vec,rprior = rpri,N = p_num,prior_pdf = p_pdf,target_set = obs_ind,ess = efs,n_move = nm)$output$log_Bayesian_evidence
      } else{
        IBIS.sub <- IBIS.Z(X = cbind(rep(1,nrow(obs_mat)),obs_mat),y = res_vec,rprior = rpri,N = p_num,prior_pdf = p_pdf,target_set = cache_smc$get(cache_smc$keys()[which.min(sub_size)])$input,ess = efs,n_move = nm)
        logZ <- IBIS.Z(X = cbind(rep(1,nrow(obs_mat)),obs_mat),y = res_vec,sampl = IBIS.sub$ouput,prior_pdf = p_pdf,target_set = obs_ind,current_set = IBIS.sub$input,ess = efs,n_move = nm)$output$log_Bayesian_evidence
      }
    }
  }
  return(logZ)
}

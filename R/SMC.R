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

##' SMC sampler function using reverse iterated batch importance sampling
##'
##'
##' @export
##' @name RIBIS.Z
##' @description This function uses a reverse iterated batch importance sampling
##' scheme with batch size one to go from one bridging distribution to another
##' by removing observations from the posterior. We assume a Bayesian logistic
##' regression model.
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
##' @param prior_pdf Probability Density Function of the prior. Must only have
##' two arguments, `th` and `di` (a vector or matrix of regression coefficients
##' samples and the number of dimensions of a single sample respectively).
##' @param drop_set Vector of observation indices to be removed to arrive at the
##' desired bridging distribution.
##' @param in_set Vector of observation indices already added to give the
##' starting bridging distribution.
##' @param ess Threshold: if the effective sample size of the particle weights
##' falls below this value then a resample move step is triggered. Defaults to
##' `N/2`.
##' @param n_move Number of Metropolis-Hastings steps to apply each time a
##' resample move step is triggered. Defaults to 1.
##' @return A list consisting of; weighted posterior samples (the samples and
##' weights are given in separate lists), the log Bayesian evidence of the final
##' partial posterior and the duplication table of the outputted posterior
##' samples.
##' @details If the samples from the partial posterior are not weighted samples
##' then the `weights` element of the `sampl` list should be
##' `rep(1,nrow(sampl$samples))`. If the samples provided in `sampl` are all
##' unique then the `duplication_table` element of `sampl` should be a table of
##' `sample$weights = rep(1,nrow(sampl$samples))`.
##'
##' Details of the internal mechanisms of the SMC sampler such as the
##' Metropolis-Hastings MCMC resample move can be found in *UNCOVER paper* and
##' Chopin (2002).
##'
##' Note that decreasing `ess` and increasing `n_move` will lead to a more
##' accurate estimate of the Bayesian evidence, but at the cost of increased
##' computational time.
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
##' unique_ind <- as.numeric(names(out.1$duplication_table))
##' unique_samp <- out.1$samples[unique_ind,]
##' unique_weights <- out.1$weights[unique_ind]
##' weight_fade <- colorRampPalette(c('white','black'))(100)[as.numeric(cut(unique_weights,breaks = 100))]
##' pairs(unique_samp,col = weight_fade)
##' out.1$log_Bayesian_evidence
##'
##' # If we then wished to remove 25 of those observations from our posterior
##' out.2 <- RIBIS.Z(X = DM,y = rv,sampl = out.1,prior_pdf = pr_fun,drop_set=100:76,in_set = 1:100)
##' unique_ind.2 <- as.numeric(names(out.2$duplication_table))
##' unique_samp.2 <- out.2$samples[unique_ind.2,]
##' unique_weights.2 <- out.2$weights[unique_ind.2]
##' weight_fade.2 <- colorRampPalette(c('white','black'))(100)[as.numeric(cut(unique_weights.2,breaks = 100))]
##' pairs(unique_samp.2,col = weight_fade.2)
##' out.2$log_Bayesian_evidence
##'
##' # we can also compare this method to standard iterated batch importance
##' # sampling
##' out.3 <- IBIS.Z(X = DM,y = rv,rprior = pr_samp,N = 1000,prior_pdf = pr_fun,add_set = 1:50)
##' unique_ind.3 <- as.numeric(names(out.3$duplication_table))
##' unique_samp.3 <- out.3$samples[unique_ind.3,]
##' unique_weights.3 <- out.3$weights[unique_ind.3]
##' weight_fade.3 <- colorRampPalette(c('white','red'))(100)[as.numeric(cut(unique_weights.3,breaks = 100))]
##' out.4 <- RIBIS.Z(X = DM,y = rv,sampl = out.1,prior_pdf = pr_fun,drop_set = 100:51,in_set = 1:100)
##' unique_ind.4 <- as.numeric(names(out.4$duplication_table))
##' unique_samp.4 <- out.4$samples[unique_ind.4,]
##' unique_weights.4 <- out.4$weights[unique_ind.4]
##' weight_fade.4 <- colorRampPalette(c('white','green'))(100)[as.numeric(cut(unique_weights.4,breaks = 100))]
##' pairs(rbind(unique_samp.3,unique_samp.4),col=c(weight_fade.3,weight_fade.4))


RIBIS.Z <- function(X,y,sampl,prior_pdf,drop_set,in_set,ess=NULL,n_move=1){
  drop_vec <- in_vec <- rep(0,length(y))
  drop_vec[drop_set] <- 1
  in_vec[in_set] <- 1
  X <- as.matrix(X)
  p <- ncol(X)
  if(sum(drop_vec*in_vec)!=length(drop_set)){
    stop("drop_set must be a subset of in_set")
  }
  N <- nrow(sampl$samples)
  if(is.null(ess)){
    ess <- N/2
  }
  w <- sampl$weights
  logZ <- sampl$log_Bayesian_evidence
  dup_tab <- sampl$duplication_table
  sampl <- sampl$samples
  drop_set <- drop_set[sample(1:length(drop_set),length(drop_set))]
  ess_point <- length(in_set)
  for(i in 1:length(drop_set)){
    w_new <- w*as.vector(((1+exp(-X[drop_set[i],]%*%t(sampl)))^(y[drop_set[i]]))*((1+exp(X[drop_set[i],]%*%t(sampl)))^(1-y[drop_set[i]])))
    logZ <- logZ + log(sum(w_new)) - log(sum(w))
    w <- w_new
    wsq <- sum((w[as.integer(names(dup_tab))]*dup_tab)^2)
    if(wsq==Inf){
      ESS <- 0
    } else{
      ESS <- (sum(w)^2)/wsq
    }
    if(ESS < ess){
      if(length(in_set)==1){
        cont_rate <- ess_point
      } else{
        cont_rate <- ess_point/(length(in_set)-1)
      }
      ess_point <- length(in_set)-1
      ss <- cov.wt(sampl,wt=w,method = "ML")
      mu <- ss$center
      Sigma <- ss$cov
      samp <- sample(as.integer(names(dup_tab)),N,prob = w[as.integer(names(dup_tab))]*dup_tab,replace = TRUE)
      sampl <- sampl[samp,]
      w <- rep(1,N)
      sampl <- t(((t(sampl) - mu)*sqrt(cont_rate)) + mu)
      A.all <- rep(FALSE,N)
      for(j in 1:n_move){
        BC <- mvnfast::rmvn(N,mu = mu, sigma=Sigma)
        if(length(setdiff(in_set,drop_set[i]))==0){
          Log_1 <- log(prior_pdf(th = BC,di = p)) + log(mvnfast::dmvn(sampl,mu=mu,sigma=Sigma))
          Log_2 <- log(prior_pdf(th = sampl,di = p)) + log(mvnfast::dmvn(BC,mu=mu,sigma=Sigma))
        } else{
          Log_1 <- log(1+exp(X[setdiff(in_set,drop_set[i]),]%*%t(BC)))*(y[setdiff(in_set,drop_set[i])]-1) + log(1+exp(-X[setdiff(in_set,drop_set[i]),]%*%t(BC)))*(-y[setdiff(in_set,drop_set[i])])
          Log_1 <- colSums(Log_1) + log(prior_pdf(th = BC,di = p)) + log(mvnfast::dmvn(sampl,mu=mu,sigma=Sigma))
          Log_2 <- log(1+exp(X[setdiff(in_set,drop_set[i]),]%*%t(sampl)))*(y[setdiff(in_set,drop_set[i])]-1) + log(1+exp(-X[setdiff(in_set,drop_set[i]),]%*%t(sampl)))*(-y[setdiff(in_set,drop_set[i])])
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
    in_set <- setdiff(in_set,drop_set[i])
  }
  return(list(samples = sampl,
              weights = w,
              log_Bayesian_evidence = logZ,
              duplication_table = dup_tab[order(as.integer(names(dup_tab)))]))
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
##' @param method Method to be used to estimate log(Z), can be one of `"SMC"`,
##' `"SMC_BIC"` or `"BIC"`
##' @param thres The threshold for which the number of observations needs to
##' exceed to consider using BIC as an estimator. Only applies if method
##' `"SMC_BIC"` is selected. Defaults to infinity if not specified.
##' @param obs_mat Covariate matrix
##' @param res_vec Binary response vector
##' @param p_num Number of samples of the prior used for the SMC sampler. Not
##' required if method `"BIC"` selected.
##' @param rpri Function to sample from the prior. Must only have two arguments,
##' `p_n` and `di` (Number of prior samples to generate and the number of
##' dimensions of a single sample respectively).
##' @param p_pdf Probability Density Function of the prior. Must only have two
##' arguments, `th` and `di` (a vector or matrix of regression coefficients
##' samples and the number of dimensions of a single sample respectively).
##' @return An estimation of the log Bayesian evidence
##' @details log(Z) is estimated using three possible methods:
##'
##' 1. `"SMC"`: Estimates log(Z) using an SMC sampler. See `IBIS.Z` for more
##' details.
##'
##' 2. `"SMC_BIC"`: Initially tries to estimate log(Z) using BIC but will revert
##' to using an SMC sampler if the data is linearly separable. Also reverts to
##' using an SMC sampler if the number of observations is below a certain
##' threshold. Recommended method and is set as default.
##'
##' 3. `"BIC"`: Estimates log(Z) with the Bayesian Information Criterion (BIC).
##' Does not require prior specification but only advisable if the number of
##' observations is large.
##' @examples
##'
##' # First we generate a covariate matrix `obs_mat` and binary response vector
##' # `res_vec`
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
##'@seealso [IBIS.Z]

lbe.gen <- function(method="SMC_BIC",thres=Inf,obs_mat,res_vec,p_num=0,rpri=NULL,
                    p_pdf=NULL,subs = NULL,efs = length(res_vec)/2,nm = 1){
  opts <- c("SMC","SMC_BIC","BIC")
  if(!(method%in%opts)){
    stop("Method is not supported")
  }
  if(method!="BIC" & (is.null(rpri) | is.null(p_pdf))){
    stop("Both sampling function and probability density function of the prior are required for this method")
  }
  cache_search_fun_bic <- function(u,cb,oi){
    c_ind <- cb$get(u)$input
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
  cache_search_fun_smc <- function(u,cs,oi){
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
        return(loi-lci)
      } else{
        return(Inf)
      }
    }
  }
  if(method=="BIC"){
    if(length(res_vec)<memo_thres | length(cache_bic$keys())==0){
      logZ <- memo.bic(X = obs_mat, y = res_vec, which_obs = obs_ind)$logZ
    } else{
      sub_size <- sapply(cache_bic$keys(),FUN = cache_search_fun_bic,cb = cache_bic,oi = obs_ind)
      if(min(sub_size)==0){
        logZ <- memo.bic(X = obs_mat, y = res_vec, which_obs = obs_ind)$logZ
      } else{
        bic.start <- memo.bic(X = obs_mat, y = res_vec, which_obs = cache_bic$get(cache_bic$keys()[which.min(sub_size)])$input)$coeffs
        logZ <- memo.bic(X = obs_mat, y = res_vec, which_obs = obs_ind,param_start = bic.start)$logZ
      }
    }
  }
  if(method=="SMC"){
    if(length(res_vec)<memo_thres | length(cache_smc$keys())==0){
      logZ <- IBIS.Z(X = cbind(rep(1,nrow(obs_mat)),obs_mat),y = res_vec,rprior = rpri,N = p_num,prior_pdf = p_pdf,add_set = obs_ind,ess = efs,n_move = nm)$log_Bayesian_evidence
    } else{
      sub_size <- sapply(cache_smc$keys(),FUN = cache_search_fun_smc,cs = cache_smc,oi = obs_ind)
      if(min(sub_size)==0){
        logZ <- IBIS.Z(X = cbind(rep(1,nrow(obs_mat)),obs_mat),y = res_vec,rprior = rpri,N = p_num,prior_pdf = p_pdf,add_set = obs_ind,ess = efs,n_move = nm)$log_Bayesian_evidence
      } else{
        if(sub_size[which.min(abs(sub_size))]>0){
          IBIS.sub <- IBIS.Z(X = cbind(rep(1,nrow(obs_mat)),obs_mat),y = res_vec,rprior = rpri,N = p_num,prior_pdf = p_pdf,target_indices = cache_smc$get(cache_smc$keys()[which.min(abs(sub_size))])$input,ess = efs,n_move = nm)$output
          logZ <- IBIS.Z(X = cbind(rep(1,nrow(obs_mat)),obs_mat),y = res_vec,sampl = IBIS.sub,prior_pdf = p_pdf,target_indices = obs_ind,ess = efs,n_move = nm)$output$log_Bayesian_evidence
        }
      }
    }
  }

  if(is.null(subs)){
    subs <- length(res_vec)
  }
  if(length(subs)<thres){
    cache_pref <- "Small"
  } else{
    ss_imp <- max(worst_set_subset[1,]) - abs(2*length(subs) - length(res_vec))
    bs_imp <- length(subs) - min(worst_set_big[1,])
    if(ss_imp < 0 & bs_imp < 0){
      cache_pref <- "None"
    } else if(ss_imp >= bs_imp){
      cache_pref <- "Subset"
    } else{
      cache_pref <- "Big"
    }
  }
  if(method=="BIC"){
    if(cache_pref=="Small"){
      return(memo.BIC.small(subs,obs_mat,res_vec)$BIC)
    } else{
      if(length(cache_subset$keys())==0){
        max_sub <- c()
      } else{
        max_sub <- c()
        for(i in cache_subset$keys()){
          mem <- cache_subset$get(u)$indices
          if(length(mem)>length(max_sub)){
            if(identical(mem,wl[1:length(mem)])){
              max_sub <- mem
            }
          }
        }
      }
    }
    if(cache_pref=="None"){
      if(length(max_sub)==0){
        return(no.memo.BIC(subs,obs_mat,res_vec)$BIC)
      } else{
        coeffs <- memo.BIC.subset(max_sub,obs_mat,res_vec)$coeffs
        return(no.memo.BIC.(subs,obs_mat,res_vec,coeffs)$BIC)
      }
    }
    if(cache_pref=="Subset"){
      if(cache_subset$size()==max_cache_size){
        drop_ind <- which.max(worst_set_subset[1,])
        cache_subset$remove(worst_set_subset[2,drop_ind])
        worst_set_subset[,drop_ind] <- c(length(res_vec) + 1,"")
      }
      if(length(max_sub)==0 | length(max_sub)==length(subs)){
        temp <- memo.BIC.subset(subs,obs_mat,res_vec)$BIC
      } else{
        coeffs <- memo.BIC.subset(max_sub,obs_mat,res_vec)$coeffs
        temp <- memo.BIC.subset(subs,obs_mat,res_vec,coeffs)$BIC
      }
      worst_set_subset[,which.max(worst_set_subset[1,])] <- c(abs(2*length(subs) - length(res_vec)),cache_subset$keys()[length(cache_subset$keys())])
      return(temp)
    }
    if(cache_pref=="Big"){
      if(cache_big$size()==max_cache_size){
        drop_ind <- which.min(worst_set_big[1,])
        cache_big$remove(worst_set_big[2,drop_ind])
        worst_set_big[,drop_ind] <- c(0,"")
      }
      if(length(max_sub)==0 | length(max_sub)==length(subs)){
        temp <- memo.BIC.subset(subs,obs_mat,res_vec)$BIC
      } else{
        coeffs <- memo.BIC.subset(max_sub,obs_mat,res_vec)$coeffs
        temp <- memo.BIC.subset(subs,obs_mat,res_vec,coeffs)$BIC
      }
      worst_set_subset[,which.min(worst_set_subset[1,])] <- c(abs(2*length(subs) - length(res_vec)),cache_subset$keys()[length(cache_subset$keys())])
      return(temp)
    }
    else{
      #keys <- get("_cache", envir=environment(memo.BIC.subset))$keys()
      if(length(cache_subset$keys())==0){
        max_sub <- c()
      } else{
        max_sub <- c()
        for(i in cache_subset$keys()){
          mem <- cache_subset$get(u)$indices
          if(length(mem)>length(max_sub)){
            if(identical(mem,wl[1:length(mem)])){
              max_sub <- mem
            }
          }
        }
      }
      if(cache_subset$size()<max_cache_size){
        if(length(max_sub)==0 | length(max_sub)==length(subs)){
          temp <- memo.BIC.subset(subs,obs_mat,res_vec)$BIC
          if(length(max_sub)==0){
            worst_set_subset[,cache_subset$size()] <- c(abs(2*length(subs) - length(res_vec)),cache_subset$keys()[length(cache_subset$keys())])
            worst_set_big[,cache_subset$size()] <- c(length(subs),cache_subset$keys()[length(cache_subset$keys())])
          }
        } else{
          coeffs <- memo.BIC.subset(max_sub,obs_mat,res_vec)$coeffs
          temp <- memo.BIC.subset(subs,obs_mat,res_vec,coeffs)$BIC
          worst_set_subset[,cache_subset$size()] <- c(abs(2*length(subs) - length(res_vec)),cache_subset$keys()[length(cache_subset$keys())])
          worst_set_big[,cache_subset$size()] <- c(length(subs),cache_subset$keys()[length(cache_subset$keys())])
        }
        if(cache_subset$size()==max_cache_size){
          cache_big <- cache_subset
          memo.BIC.big <- memoise::memoise(no.memo.BIC,cache = cache_big,omit_args = c("X","y","coeff_start"))
        }
      }
      if(max(part_out)==0){
        if(abs(2*length(subs) - length(res_vec))<subset_max){
          return(memo.BIC.subset(subs,obs_mat,res_vec)$BIC)
        } else if(length(subs) > big_min){
          return(memo.BIC.big(subs,obs_mat,res_vec)$BIC)
        } else{
          return(no.memo.BIC(subs,obs_mat,res_vec)$BIC)
        }
      } else if(max(part_out)==length(subs)){
        return(memo.BIC.subset(subs,obs_mat,res_vec)$BIC)
      } else{
        part_idx <- which.max(part_out)
        seen_before <- (get("_cache", envir=environment(memo.BIC.subset))$get(keys[part_idx]))$value[[1]]
        if(abs(2*length(subs) - length(res_vec))<subset_max){
          return(memo.BIC.subset(subs,obs_mat,res_vec,seen_before)$BIC)
        } else if(length(subs) > big_min){
          return(memo.BIC.big(subs,obs_mat,res_vec,seen_before)$BIC)
        } else{
          return(no.memo.BIC(subs,obs_mat,res_vec,seen_before)$BIC)
        }
      }
    }
    #return(-BIC(glm(res_vec~.,family="binomial",data=data.frame(obs_mat,res_vec)))/2)
  } else{
    logZ <- T
    if(method=="SMC_BIC" & length(res_vec)>=thres){
      logZ <- tryCatch(-BIC(glm(res_vec~.,family="binomial",data=data.frame(obs_mat,res_vec)))/2,warning=function(...)T)
    }
    if(logZ==T){
      logZ <- IBIS.Z(X = cbind(rep(1,nrow(obs_mat)),obs_mat), y = res_vec, sampl = rpri(p_n = p_num,di = (ncol(obs_mat)+1)),prior_pdf = p_pdf)[[3]]
    }
    return(logZ)
  }
}

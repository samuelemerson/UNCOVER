IBIS.Z <- function(X,y,sampl=NULL,rprior=NULL,N=NULL,prior_pdf,
                   target_set=1:length(y),current_set=NULL,ess=NULL,n_move=1,
                   PriorArgs,diagnostics = FALSE){
  if(length(target_set)==length(current_set)){
    stop("target_set and current_set cannot be the same length")
  }
  if(length(target_set)>length(current_set)){
    reverse <- FALSE
  } else{
    reverse <- TRUE
    if(!identical(sampl$weights,rep(1,length(sampl$weights)))){
      N <- nrow(sampl$samples)
      p <- ncol(sampl$samples)
      ss <- stats::cov.wt(sampl$samples,wt=sampl$weights,method = "ML")
      mu <- ss$center
      Sigma <- ss$cov
      samp <- sample(as.integer(names(sampl$duplication_table)),N,
                     prob = sampl$weights[as.integer(names(sampl$duplication_table))]*sampl$duplication_table,
                     replace = TRUE)
      sampl$samples <- sampl$samples[samp,]
      sampl$weights <- rep(1,N)
      A.all <- rep(FALSE,N)
      for(j in 1:n_move){
        BC <- mvnfast::rmvn(N,mu = mu, sigma=Sigma)
        Log_1 <- log(1+exp(X[current_set,]%*%t(BC)))*(y[current_set]-1) + log(1+exp(-X[current_set,]%*%t(BC)))*(-y[current_set])
        TotArgs <- list(BC)
        if(length(PriorArgs)!=0){
          TotArgs <- c(TotArgs,PriorArgs)
        }
        Log_1 <- colSums(Log_1) + log(do.call(prior_pdf,TotArgs)) + log(mvnfast::dmvn(sampl$samples,mu=mu,sigma=Sigma))
        Log_2 <- log(1+exp(X[current_set,]%*%t(sampl$samples)))*(y[current_set]-1) + log(1+exp(-X[current_set,]%*%t(sampl$samples)))*(-y[current_set])
        TotArgs <- list(sampl$samples)
        if(length(PriorArgs)!=0){
          TotArgs <- c(TotArgs,PriorArgs)
        }
        Log_2 <- colSums(Log_2) + log(do.call(prior_pdf,TotArgs)) + log(mvnfast::dmvn(BC,mu=mu,sigma=Sigma))
        A <- Log_1 - Log_2 - log(stats::runif(length(Log_1)))
        A.all <- A.all | (A>0)
        sampl$samples[which(A>0),] <- BC[which(A>0),]
      }
      samp[which(A.all)] <- (N+1):(N+length(which(A.all)))
      sampl$duplication_table <- table(samp)
      names(sampl$duplication_table) <- as.character(match(as.integer(names(sampl$duplication_table)),samp))
    }
    ess_point <- length(current_set)
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
    TotArgs <- list(N)
    if(length(PriorArgs)!=0){
      TotArgs <- c(TotArgs,PriorArgs)
    }
    sampl <- do.call(rprior,TotArgs)
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
  if(diagnostics){
    ESS.rec <- data.frame(Observations = action_set, ESS = rep(0,length(action_set)))
    LBE.rec <- ESS.rec
    colnames(LBE.rec)[2] <- "Log_Bayesian_Evidence"
    AR.rec <- data.frame()
  }
  for(i in 1:length(action_set)){
    w_new <- as.vector(((1+exp(-X[action_set[i],]%*%t(sampl)))^(-y[action_set[i]]))*((1+exp(X[action_set[i],]%*%t(sampl)))^(y[action_set[i]]-1)))
    if(reverse){
      w_new <- w/w_new
    } else{
      w_new <- w*w_new
    }
    logZ <- logZ + log(sum(w_new)) - log(sum(w))
    if(diagnostics){
      LBE.rec[i,2] <- logZ
    }
    w <- w_new
    wsq <- sum((w[as.integer(names(dup_tab))]*dup_tab)^2)
    if(wsq==0){
      ESS <- 0
    } else{
      ESS <- (sum(w)^2)/wsq
    }
    if(diagnostics){
      ESS.rec[i,2] <- ESS
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
      ss <- stats::cov.wt(sampl,wt=w,method = "ML")
      mu <- ss$center
      Sigma <- ss$cov
      samp <- sample(as.integer(names(dup_tab)),N,
                     prob = w[as.integer(names(dup_tab))]*dup_tab,
                     replace = TRUE)
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
            TotArgs <- list(BC)
            if(length(PriorArgs)!=0){
              TotArgs <- c(TotArgs,PriorArgs)
            }
            Log_1 <- log(do.call(prior_pdf,TotArgs)) + log(mvnfast::dmvn(sampl,mu=mu,sigma=Sigma))
            TotArgs <- list(sampl)
            if(length(PriorArgs)!=0){
              TotArgs <- c(TotArgs,PriorArgs)
            }
            Log_2 <- log(do.call(prior_pdf,TotArgs)) + log(mvnfast::dmvn(BC,mu=mu,sigma=Sigma))
          } else{
            Log_1 <- log(1+exp(X[setdiff(current_set,action_set[i]),]%*%t(BC)))*(y[setdiff(current_set,action_set[i])]-1) + log(1+exp(-X[setdiff(current_set,action_set[i]),]%*%t(BC)))*(-y[setdiff(current_set,action_set[i])])
            TotArgs <- list(BC)
            if(length(PriorArgs)!=0){
              TotArgs <- c(TotArgs,PriorArgs)
            }
            Log_1 <- colSums(Log_1) + log(do.call(prior_pdf,TotArgs)) + log(mvnfast::dmvn(sampl,mu=mu,sigma=Sigma))
            Log_2 <- log(1+exp(X[setdiff(current_set,action_set[i]),]%*%t(sampl)))*(y[setdiff(current_set,action_set[i])]-1) + log(1+exp(-X[setdiff(current_set,action_set[i]),]%*%t(sampl)))*(-y[setdiff(current_set,action_set[i])])
            TotArgs <- list(sampl)
            if(length(PriorArgs)!=0){
              TotArgs <- c(TotArgs,PriorArgs)
            }
            Log_2 <- colSums(Log_2) + log(do.call(prior_pdf,TotArgs)) + log(mvnfast::dmvn(BC,mu=mu,sigma=Sigma))
          }
        } else{
          Log_1 <- log(1+exp(X[c(current_set,action_set[i]),]%*%t(BC)))*(y[c(current_set,action_set[i])]-1) + log(1+exp(-X[c(current_set,action_set[i]),]%*%t(BC)))*(-y[c(current_set,action_set[i])])
          TotArgs <- list(BC)
          if(length(PriorArgs)!=0){
            TotArgs <- c(TotArgs,PriorArgs)
          }
          Log_1 <- colSums(Log_1) + log(do.call(prior_pdf,TotArgs)) + log(mvnfast::dmvn(sampl,mu=mu,sigma=Sigma))
          Log_2 <- log(1+exp(X[c(current_set,action_set[i]),]%*%t(sampl)))*(y[c(current_set,action_set[i])]-1) + log(1+exp(-X[c(current_set,action_set[i]),]%*%t(sampl)))*(-y[c(current_set,action_set[i])])
          TotArgs <- list(sampl)
          if(length(PriorArgs)!=0){
            TotArgs <- c(TotArgs,PriorArgs)
          }
          Log_2 <- colSums(Log_2) + log(do.call(prior_pdf,TotArgs)) + log(mvnfast::dmvn(BC,mu=mu,sigma=Sigma))
        }
        A <- Log_1 - Log_2 - log(stats::runif(length(Log_1)))
        A.all <- A.all | (A>0)
        sampl[which(A>0),] <- BC[which(A>0),]
        if(diagnostics){
          AR.rec <- rbind(AR.rec,c(paste0(action_set[i],".",j),sum(A>0)/N))
        }
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
  if(diagnostics){
    if(nrow(AR.rec)>0){
      colnames(AR.rec) <- c("Observations","Acceptance_Rate")
    }
    return(list(output = list(samples = sampl,
                              weights = w,
                              log_Bayesian_evidence = logZ,
                              duplication_table = dup_tab[order(as.integer(names(dup_tab)))],
                              diagnostics = list(log_Bayesian_evidence_tracker = LBE.rec,
                                                 effective_sample_size_tracker = ESS.rec,
                                                 acceptance_rate_tracker = AR.rec,
                                                 ESS_threshold = ess)),
                input = target_set))
  } else{
    return(list(output = list(samples = sampl,
                              weights = w,
                              log_Bayesian_evidence = logZ,
                              duplication_table = dup_tab[order(as.integer(names(dup_tab)))]),
                input = target_set))
  }
}

memo.bic <- function(X,y,which_obs=1:length(y),param_start=NULL){
  bic_glm <- stats::glm(y~.,family="binomial",data=data.frame(X,y)[which_obs,],start = param_start)
  lbe <- -stats::BIC(bic_glm)/2
  return(list(logZ = lbe,coeffs = bic_glm$coefficients,input = which_obs))
}

lbe.gen <- function(obs_mat,res_vec,obs_ind = 1:length(res_vec),thres=30,
                    memo_thres_bic = Inf,memo_thres_smc = Inf,p_num=1000,rpri,
                    p_pdf,efs = p_num/2,nm = 1,cache_bic,cache_smc,MA,SMC_fun,
                    BIC_fun,ribis_thres=30){
  cache_search_fun_bic <- function(u,cs,oi){
    c_ind <- cs$get(u)$value$input
    return(length(setdiff(oi,c_ind)) + length(setdiff(c_ind,oi)))
  }
  cache_search_fun_smc <- function(u,cs,oi,rt){
    c_ind <- cs$get(u)$value$input
    lci <- length(c_ind)
    loi <- length(oi)
    if(lci==loi){
      if(sum(c_ind - oi)==0){
        return(lci-loi)
      } else{
        return(Inf)
      }
    }
    if(lci<loi){
      if(all(c_ind%in%oi)){
        return(loi-lci)
      } else{
        return(Inf)
      }
    }
    if(lci>loi){
      if((lci-loi >= loi) | (loi < rt)){
        return(Inf)
      } else{
        if(all(oi%in%c_ind)){
          return(lci-loi)
        } else{
          return(Inf)
        }
      }
    }
  }
  logZ <- TRUE
  if(length(obs_ind)>=thres){
    if(length(obs_ind)<memo_thres_bic | length(cache_bic$keys())==0){
      logZ <- tryCatch(BIC_fun(X = obs_mat, y = res_vec,which_obs = obs_ind)$logZ,
                       warning=function(...)TRUE)
    } else{
      sub_size <- sapply(cache_bic$keys(),FUN = cache_search_fun_bic,cs = cache_bic,
                         oi = obs_ind)
      if(min(sub_size)==0 | min(sub_size)==Inf){
        logZ <- tryCatch(BIC_fun(X = obs_mat, y = res_vec, which_obs = obs_ind)$logZ,
                         warning=function(...)TRUE)
      } else{
        bic.start <- BIC_fun(X = obs_mat, y = res_vec, which_obs = cache_bic$get(cache_bic$keys()[which.min(sub_size)])$value$input)$coeffs
        logZ <- tryCatch(BIC_fun(X = obs_mat, y = res_vec, which_obs = obs_ind,
                                  param_start = bic.start)$logZ,
                         warning=function(...)TRUE)
      }
    }
  }
  if(logZ==TRUE){
    if(length(obs_ind)<memo_thres_smc | length(cache_smc$keys())==0){
      logZ <- SMC_fun(X = cbind(rep(1,nrow(obs_mat)),obs_mat),y = res_vec,
                     rprior = rpri,N = p_num,prior_pdf = p_pdf,
                     target_set = obs_ind,ess = efs,n_move = nm,
                     PriorArgs = MA)$output$log_Bayesian_evidence
    } else{
      sub_size <- sapply(cache_smc$keys(),FUN = cache_search_fun_smc,cs = cache_smc,
                         oi = obs_ind,rt = ribis_thres)
      if(min(sub_size)==0 | min(sub_size)==Inf){
        logZ <- SMC_fun(X = cbind(rep(1,nrow(obs_mat)),obs_mat),y = res_vec,
                       rprior = rpri,N = p_num,prior_pdf = p_pdf,
                       target_set = obs_ind,ess = efs,n_move = nm,
                       PriorArgs = MA)$output$log_Bayesian_evidence
      } else{
        IBIS.sub <- SMC_fun(X = cbind(rep(1,nrow(obs_mat)),obs_mat),y = res_vec,
                           rprior = rpri,N = p_num,prior_pdf = p_pdf,
                           target_set = cache_smc$get(cache_smc$keys()[which.min(sub_size)])$value$input,
                           ess = efs,n_move = nm,PriorArgs = MA)
        logZ <- SMC_fun(X = cbind(rep(1,nrow(obs_mat)),obs_mat),y = res_vec,
                       sampl = IBIS.sub$output,prior_pdf = p_pdf,
                       target_set = obs_ind,current_set = IBIS.sub$input,
                       ess = efs,n_move = nm,
                       PriorArgs = MA)$output$log_Bayesian_evidence
      }
    }
  }
  return(logZ)
}

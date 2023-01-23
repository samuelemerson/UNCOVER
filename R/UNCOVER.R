###############################################
## R script for functions in UNCOVER package ##
###############################################
##
## Sam Emerson
## September 2022
##


######################
## UNCOVER Function ##
######################


##' Utilising Normalisation Constant Optimisation Via Edge Removal
##'
##'
##' @export
##' @name UNCOVER
##' @description Generates cohorts for a data set through removal of edges from
##' a graphical representation of the covariates. Edges are removed (or
##' reintroduced) by considering the normalisation constant (or Bayesian
##' evidence) of a multiplicative Bayesian logistic regression model.
##'
##' The first stage of the function is concerned purely with a greedy
##' optimisation of the Bayeisan evidence through edge manipulation. The second
##' stage then addresses any other criteria (known as deforest conditions)
##' expressed by the user through reintroduction of edges.
##'
##' @keywords graph cohort cluster bayesian evidence
##' @param X Covariate matrix
##' @param y Binary response vector
##' @param mst_var A vector specifying which variables of the covariate matrix
##' will be used to form the graph. If not specified all variables will be used.
##' @param N Number of samples of the prior used for the SMC sampler. Not
##' required if method `"BIC"` selected. If required but not specified the
##' default value is set to `1000`.
##' @param stop_criterion What is the maximum number of clusters allowed before
##' we terminate the first stage and begin deforestation. If not specified the
##' algorithm continues until the Bayesian evidence cannot be improved upon,
##' however for time and memory purposes specifying a number is highly
##' recommended.
##' @param deforest_criterion Criterion in which edges are reintroduced to allow
##' the model to satisfy. Can be one of `"NoC"`,`"SoC"`,`"MaxReg"`,
##' `"Validation"`, `"Balanced"` or `"None"`.
##' @param split What fraction of the data should be used for training. Should
##' only be specified if `deforest_criterion == "Validation`. Defaults to `1`.
##' @param max_K The maximum number of clusters allowed in the final output.
##' Should only be specified if `deforest_criterion == "NoC`. Defaults to `Inf`.
##' @param min_size The minimum number of observations allowed for any cluster
##' in the final model. Should only be specified if
##' `deforest_criterion == "SoC`. Defaults to `0`.
##' @param reg Numerical natural logarithm of the tolerance parameter. Must be
##' positive. Should only be specified if `deforest_criterion == "MaxReg`.
##' Defaults to `0`.
##' @param n_min_class The minimum number of observations in any cluster that
##' has an associated response in the minority class of that cluster. Defaults
##' to `0`.
##' @param SMC_thres The threshold for which the number of observations needs to
##' exceed to consider using BIC as an estimator. Defaults to 30 if not
##' specified.
##' @param BIC_memo_thres The threshold for when it is deemed worthwhile to
##' check the cache of function `memo.bic` for similar observation indices.
##' Defaults to never checking the cache.
##' @param SMC_memo_thres The threshold for when it is deemed worthwhile to
##' check the cache of function `IBIS.Z` for similar observation indices.
##' Defaults to never checking the cache.
##' @param rprior Function to sample from the prior. Must only have two
##' arguments, `p_num` and `di` (Number of prior samples to generate and the
##' number of dimensions of a single sample respectively).
##' @param prior_pdf Probability Density Function of the prior. Must only have
##' two arguments, `th` and `di` (a vector or matrix of regression coefficients
##' samples and the number of dimensions of a single sample respectively).
##' @param ess Threshold: if the effective sample size of the particle weights
##' falls below this value then a resample move step is triggered. Defaults to
##' `N/2`. See `IBIS.Z` for details.
##' @param n_move Number of Metropolis-Hastings steps to apply each time a
##' resample move step is triggered. Defaults to 1. See `IBIS.Z` for details.
##' @param plot_progress Do you want to plot the output of the clustering each
##' time an edge is removed or reintroduced?
##' @param verbose Do you want the progress of the algorithm to be shown?
##' @return Either a list or a list of two lists. See details.
##' @details Assumes a Bayesian logistic regression model for each cohort, with
##' the overall model being a product of these sub-models.
##'
##' A minimum spanning tree graph is first constructed from a subset of the
##' covariates. Then at each iteration, each edge in the current graph is
##' checked to see if removal to split a cohort is beneficial, and then either
##' we selected the optimal edge to remove or we concluded it is not beneficial
##' to remove any more edges. At the end of each iteration we also check the set
##' of removed edges to see if it beneficial to reintroduce any previously
##' removed edges. After this process has ended we then reintroduce edges in the
##' removed set specifically to meet the criteria set by the user in the most
##' optimal manner possible through a greedy approach. For more details see the
##' help pages of `lbe.gen`,`one.stage.mst`,`two.stage.mst`,`remove.edge`,
##' `deforest.noc`,`deforest.soc`,`deforest.maxreg`,`deforest.validation`,
##' `deforest.balanced` and the paper *UNCOVER*.
##'
##' If `rprior` and `prior_pdf` are not specified then the default prior is a
##' standard multivariate normal.
##'
##' The graph can be undergo deforestation to meet 6 possible criteria:
##'
##' 1. `"NoC"`: Number of Clusters - we specify a maximum number of clusters
##' (`max_K`) we can tolerate in the final output of the algorithm. See
##' `deforest.noc` for more details.
##'
##' 2. `"SoC"`: Size of Clusters - we specify a minimum number of observations
##' (`min_size`) we can tolerate being assigned to a cluster in the final output
##' of the algorithm. See `deforest.soc` for more details.
##'
##' 3. `"MaxReg"`: Maximal Regret - we give a maximum tolerlance (`exp(reg)`)
##' that we allow the Bayesian evidence to decrease by by reintroducing an edge.
##' See `deforest.maxreg` for more details.
##'
##' 4. `"Validation"`: Validation Data - we split (using `split`) the data into
##' training and validation data, apply the first stage of the algorithm on the
##' training data and the introduce the validation data for the deforestation
##' stage. See `deforest.valaidation` for more details.
##'
##' 5. `"Balanced"`: Balanced Response Class Within Clusters - We specify a
##' minimum number of observations (`n_min_class`) in a cluster that have the
##' minority response class associated to them (the minimum response class is
##' determined for each cluster).
##'
##' 6. `"None"`: No Criteria Specified - we do not go through the second
##' deforestation stage of the algorithm.
##'
##' For more details on the specifics of the possible values for `SMC_method`,
##' see the help page of the function `lbe.gen`.
##'
##' If any deforestation criterion other than `"Validation"` is chosen, then the
##' output will be a list of the following:
##'
##' 1. Cluster Allocation - A vector indicating which cluster each observation
##' belongs to.
##'
##' 2. Log Marginal Likelihoods - A vector of the log Bayesian evidences (or log
##' marginal likelihoods) of each of the clusters. The sum of this vector will
##' be the log Bayesian evidence of the overall model.
##'
##' 3. Graph - The final graph of the data.
##'
##' 4. Number of Clusters - Total number of clusters (or cohorts) in the final
##' output.
##'
##' 5. Edges Removed - A matrix of the edges removed, expressed as the two
##' vertices in the graph that the edge connected.
##'
##' If the deforestation criterion chosen is `"Validation"`, then we produce a
##' list of two lists. The first is the list given above for only the training
##' data, and the second is the list given above for both the training data and
##' the valaidation data combined.
##'
##' @seealso [lbe.gen,one.stage.mst,two.stage.mst,remove.edge,deforest.noc,deforest.soc,deforest.maxreg,deforest.validation,deforest.balanced,IBIS.Z,memo.bic]
##' @examples
##'
##' # First we generate a covariate matrix `X` and binary response vector `y`
##' CM <- matrix(rnorm(200),100,2)
##' rv <- sample(0:1,100,replace=T)
##'
##' # We can then run our algorithm to see what cohorts are selected for each
##' # of the different deforestation criteria
##' UN.none <- UNCOVER(X = CM,y = rv, stop_criterion = 8, deforest_criterion = "None", verbose = F)
##' UN.noc <- UNCOVER(X = CM,y = rv, stop_criterion = 8, deforest_criterion = "NoC", max_K = 3, verbose = F)
##' UN.soc <- UNCOVER(X = CM,y = rv, stop_criterion = 8, deforest_criterion = "SoC", min_size = 10, verbose = F)
##' UN.maxreg <- UNCOVER(X = CM,y = rv, stop_criterion = 8, deforest_criterion = "MaxReg", reg = 1, verbose = F)
##' UN.validation <- UNCOVER(X = CM,y = rv, stop_criterion = 8, deforest_criterion = "Validation", split = 0.8, verbose = F)
##' UN.balanced <- UNCOVER(X = CM,y = rv, stop_criterion = 8, deforest_criterion = "Balanced", n_min_class = 2, verbose = F)
##' clu_al_mat <- rbind(UN.none[[1]],UN.noc[[1]],UN.soc[[1]],UN.maxreg[[1]],UN.validation[[2]][[1]],UN.balanced[[1]])
##' # We can create a matrix where each entry shows in how many of the methods
##' # did the indexed observations belong to the same cluster
##' obs_con_mat <- matrix(0,100,100)
##' for(i in 1:100){
##' for(j in 1:100){
##' obs_con_mat[i,j] <- length(which(clu_al_mat[,i]-clu_al_mat[,j]==0))/6
##' obs_con_mat[j,i] <- obs_con_mat[i,j]
##' }
##' }
##' obs_con_mat
##' # We can also view the outputted overall Bayesian evidence of the five
##' # models as well
##' c(sum(UN.none[[2]]),sum(UN.noc[[2]]),sum(UN.soc[[2]]),sum(UN.maxreg[[2]]),sum(UN.validation[[2]][[2]]),sum(UN.balanced[[2]]))
##'
##' # If we don't assume the prior for the regression coefficients is a
##' # standard normal but instead a multivariate normal with mean (1,1) and the
##' # identity matrix as the covariance matrix we can specify
##' pr_samp <- function(p_n,di){return(mvnfast::rmvn(p_n,rep(1,di),diag(di)))}
##' pr_fun <- function(th,di){return(mvnfast::dmvn(th,mu=rep(1,di),sigma=diag(di)))}
##'
##' # We then can run UNCOVER using this prior and compare to the standard result
##' UN.none.2 <- UNCOVER(X = CM,y = rv, stop_criterion = 8, deforest_criterion = "None", rprior = pr_samp,prior_pdf = pr_fun,verbose = F)
##' c(sum(UN.none[[2]]),sum(UN.none.2[[2]]))


UNCOVER <- function(X,y,mst_var=NULL,N=1000,stop_criterion=Inf,
                    deforest_criterion="None",split=1,max_K=Inf,min_size=0,
                    reg=0,n_min_class = 0,SMC_thres=30,BIC_memo_thres = Inf,
                    SMC_memo_thres = Inf,rprior=NULL,prior_pdf=NULL,ess=N/2,
                    n_move=1,plot_progress=F,verbose = T){
  if((!is.null(rprior)&is.null(prior_pdf))|(is.null(rprior)&!is.null(prior_pdf))){
    stop("Both sampling function and probability density function of the prior are required")
  }
  if(is.null(rprior)&is.null(prior_pdf)){
    rprior <- function(p_num,di){
      return(mvnfast::rmvn(p_num,rep(0,di),diag(di)))
    }
    prior_pdf <- function(th,di){
      return(mvnfast::dmvn(th,mu=rep(0,di),sigma=diag(di)))
    }
  }
  if(is.null(mst_var)){
    mst_var <- 1:ncol(X)
  }
  if(deforest_criterion=="Validation"){
    samp <- sort(sample((1:length(y)),round(length(y)*split)))
    g_val <- two.stage.mst(obs_mat = X,tr_ind = samp,mst_sub = mst_var)
    g <- g_val[[1]]
    g_all <- g_val[[2]]
    X_all <- X
    y_all <- y
    X <- X[samp,,drop=F]
    y <- y[samp]
  } else{
    g <- one.stage.mst(obs = X,rho = mst_var)
  }
  depth_g <- V(g)$name[dfs(g,get_diameter(g)[1])$order]
  edge_rank <- matrix(0,length(E(g)),2)
  for(i in 2:length(depth_g)){
    pot_edg <- V(g)$name[neighbors(g,depth_g[i])]
    edge_rank[i-1,] <- c(depth_g[i],pot_edg[which(pot_edg %in% depth_g[1:i])])
  }
  n <- length(y)
  K <- 1
  z <- rep(1,n)
  logZ <- lbe.gen(thres = SMC_thres,obs_mat = X,res_vec = y,memo_thres_bic = BIC_memo_thres,memo_thres_smc = SMC_memo_thres,p_num = N,rpri = rprior,p_pdf = prior_pdf,efs = ess,nm = n_move)
  edge_rem <- c()
  system_save <- vector(mode = "list",length=1)
  combine_save <- list()
  if(deforest_criterion=="NoC" | deforest_criterion=="SoC" | deforest_criterion=="Balanced"){
    model_selection <- list(z,logZ,g,K,edge_rem)
  }
  m <- 1
  while(T){
    if(verbose){
      message("")
      message(paste0("Iteration ",m))
    }
    for(k in 1:K){
      if(length(system_save[[k]][[1]])!=0){
        next
      }
      if(verbose){
        message("")
        message(green(paste0("Assessing the removal of edges from connected component ",k)))
      }
      if(deforest_criterion=="Validation"){
        system_save[[k]] <- list(z,logZ,g,c(),g_all)
      } else{
        system_save[[k]] <- list(z,logZ,g,c())
      }
      for(q in 1:nrow(edge_rank)){
        i <- get.edge.ids(g,edge_rank[q,])
        if(i==0){
          if(verbose){
            txtProgressBar(min=1,max=nrow(edge_rank),initial=q,style=3)
          }
          next
        }
        if(z[get.edgelist(g,names=F)[i,]][1]!=k){
          if(verbose){
            txtProgressBar(min=1,max=nrow(edge_rank),initial=q,style=3)
          }
          next
        }
        if(deforest_criterion=="Validation"){
          g_temp <- delete_edges(g_all,E(g_all)[get.edge.ids(g_all,get.edgelist(g)[i,])])
          z_temp <- components(g_temp)$membership
          if(!identical(sort(unique(z_temp[-samp])),as.numeric(1:(K+1)))){
            if(verbose){
              txtProgressBar(min=1,max=nrow(edge_rank),initial=q,style=3)
            }
            next
          }
        }
        er_temp <- remove.edge(gra = g,j = i,clu_al = z,lbe = logZ,obs = X,res = y,est_thres = SMC_thres,mtb = BIC_memo_thres,mts = SMC_memo_thres,par_no = N,rfun = rprior,pdf_fun = prior_pdf,efsamp = ess,methas = n_move)
        if(sum(er_temp[[2]])>sum(system_save[[k]][[2]])){
          if(deforest_criterion=="Validation"){
            system_save[[k]] <- c(er_temp,list(get.edgelist(g)[i,]),list(g_temp))
          } else{
            system_save[[k]] <- c(er_temp,list(get.edgelist(g)[i,]))
          }
        }
        if(verbose){
          txtProgressBar(min=1,max=nrow(edge_rank),initial=q,style=3)
        }
      }
      if(verbose){
        message("")
      }
    }
    logZ_regions <- sapply(system_save,FUN = function(u){sum(u[[2]])})
    if(all(logZ_regions<=sum(logZ))){
      break
    }
    k_change <- which.max(logZ_regions)
    if(verbose){
      message("")
      message(blue(paste0("Removing edge from connected component ",k_change)))
    }
    combine_save[[K]] <- list(which(z==k_change),logZ[k_change])
    z <- system_save[[k_change]][[1]]
    logZ <- system_save[[k_change]][[2]]
    g <- system_save[[k_change]][[3]]
    edge_rem <- rbind(edge_rem,system_save[[k_change]][[4]])
    if(deforest_criterion=="Validation"){
      g_all <- system_save[[k_change]][[5]]
    }
    if(plot_progress){
      pairs(X[,mst_var],pch=as.character(y),col=z,cex=0.5)
    }
    for(k in 1:K){
      if(k==k_change){
        system_save[[k]] <- vector(mode = "list",length=1)
      } else{
        system_save[[k]][[1]][which(system_save[[k]][[1]]==(K+1))] <- K+2
        system_save[[k]][[1]][which(z==(K+1))] <- K+1
        system_save[[k]][[2]] <- c(system_save[[k]][[2]],system_save[[k]][[2]][K+1])
        system_save[[k]][[2]][K+1] <- logZ[K+1]
        system_save[[k]][[2]][k_change] <- logZ[k_change]
        system_save[[k]][[3]] <- delete_edges(system_save[[k]][[3]],E(system_save[[k]][[3]])[get.edge.ids(system_save[[k]][[3]],edge_rem[nrow(edge_rem),])])
        if(deforest_criterion=="Validation"){
          system_save[[k]][[5]] <- delete_edges(system_save[[k]][[5]],E(system_save[[k]][[5]])[get.edge.ids(system_save[[k]][[5]],edge_rem[nrow(edge_rem),])])
        }
      }
    }
    K <- K+1
    system_save[[K]] <- vector(mode = "list",length=1)
    if(deforest_criterion=="NoC" | deforest_criterion=="SoC" | deforest_criterion=="Balanced"){
      unbal_clu <- sapply(1:K,FUN = function(u,res,clu_al,ups){all(table(factor(res[which(clu_al==u)],levels=0:1))>=ups)},res=y,clu_al=z,ups=n_min_class)
      if((K<=max_K & deforest_criterion=="NoC" & sum(logZ)>sum(model_selection[[2]])) | (all(table(z)>=min_size) & deforest_criterion=="SoC" & sum(logZ)>sum(model_selection[[2]])) | (all(unbal_clu) & deforest_criterion=="Balanced" & sum(logZ)>sum(model_selection[[2]]))){
        model_selection <- list(z,logZ,g,K,edge_rem)
      }
    }
    j <- 1
    if(verbose){
      message("")
      message(green("Assessing the reintroduction of edges which have been removed"))
    }
    while(j <= nrow(edge_rem)){
      edge_z <- sort(z[match(edge_rem[j,],V(g)$name)])
      if(!identical(combine_save[[j]][[1]],which(z==edge_z[1]| z==edge_z[2]))){
        combine_save[[j]] <- list(which(z==edge_z[1]| z==edge_z[2]),lbe.gen(thres = SMC_thres,obs_mat = X,res_vec = y,obs_ind = which(z==edge_z[1]| z==edge_z[2]),p_num = N,rpri = rprior,p_pdf = prior_pdf,memo_thres_bic = BIC_memo_thres,memo_thres_smc = SMC_memo_thres,efs = ess,nm = n_move))
      }
      if(verbose){
        txtProgressBar(min=1,max=nrow(edge_rem)+1,initial=j+1,style=3)
      }
      if(j==nrow(edge_rem)){
        message("")
        logZ_comb.1 <- sapply(combine_save,FUN = function(u){u[[2]]})
        logZ_comb.2 <- apply(edge_rem,MARGIN = 1,FUN = function(u,lZ,ca,cg){sum(lZ[-ca[match(u,V(cg)$name)]])},lZ = logZ,ca = z, cg = g)
        logZ_comb <- logZ_comb.1 + logZ_comb.2
        if(any(logZ_comb>sum(logZ))){
          cut_comb <- which.max(logZ_comb)
          edge_z <- sort(z[match(edge_rem[cut_comb,],V(g)$name)])
          if(verbose){
            message("")
            message(blue(paste0("Combining connected components ",edge_z[1]," and ",edge_z[2])))
          }
          system_save[[edge_z[1]]] <- vector(mode = "list",length=1)
          system_save[[edge_z[2]]] <- c()
          z[which(z==edge_z[2])] <- edge_z[1]
          for(k in (edge_z[2]+1):K){
            z[which(z==k)] <- k-1
          }
          for(k in 1:length(system_save)){
            if(is.null(system_save[[k]][[1]])){
              next
            }
            change_z <- which(system_save[[k]][[1]]==(K+1))
            system_save[[k]][[1]] <- z
            system_save[[k]][[1]][change_z] <- K
            system_save[[k]][[2]][edge_z[1]] <- combine_save[[cut_comb]][[2]]
            system_save[[k]][[2]] <- system_save[[k]][[2]][-edge_z[2]]
            system_save[[k]][[3]] <- add_edges(system_save[[k]][[3]],as.numeric(edge_rem[cut_comb,]))
            if(deforest_criterion=="Validation"){
              system_save[[k]][[5]] <- add_edges(system_save[[k]][[5]],as.numeric(edge_rem[cut_comb,]))
            }
          }
          logZ[edge_z[1]] <- combine_save[[cut_comb]][[2]]
          logZ <- logZ[-edge_z[2]]
          g <- add_edges(g,edge_rem[cut_comb,])
          if(deforest_criterion=="Validation"){
            g_all <- add_edges(g_all,edge_rem[cut_comb,])
          }
          edge_rem <- edge_rem[-cut_comb,,drop=F]
          combine_save[[cut_comb]] <- c()
          K <- K-1
          if(plot_progress){
            pairs(X[,mst_var],pch=as.character(y),col=z,cex=0.5)
          }
          j <- 0
          if(verbose){
            message("")
            message(green("Reassessing the reintroduction of edges which have been removed"))
          }
          if(deforest_criterion=="NoC" | deforest_criterion=="SoC" | deforest_criterion=="Balanced"){
            unbal_clu <- sapply(1:K,FUN = function(u,res,clu_al,ups){all(table(factor(res[which(clu_al==u)],levels=0:1))>=ups)},res=y,clu_al=z,ups=n_min_class)
            if((K<=max_K & deforest_criterion=="NoC" & sum(logZ)>sum(model_selection[[2]])) | (all(table(z)>=min_size) & deforest_criterion=="SoC" & sum(logZ)>sum(model_selection[[2]]))| (all(unbal_clu) & deforest_criterion=="Balanced" & sum(logZ)>sum(model_selection[[2]]))){
              model_selection <- list(z,logZ,g,K,edge_rem)
            }
          }
        }
      }
      j <- j+1
    }
    if(K==stop_criterion){
      break
    }
    m <- m + 1
  }
  if(verbose){
    message("")
    message("Deforestation Stage")
  }
  if(deforest_criterion=="NoC"){
    pnoc <- deforest.noc(obs = X,res = y,gra = g,lbe = logZ,eps = edge_rem,K_dag = max_K,clu_al = z,c_s = combine_save,est_thres = SMC_thres,mtb = BIC_memo_thres,mts = SMC_memo_thres,par_no = N,rfun = rprior,pdf_fun = prior_pdf,efsamp = ess,methas = n_move,p_p = plot_progress,rho = mst_var,vb = verbose)
    get("_cache", envir=environment(memo.bic))$reset()
    get("_cache", envir=environment(IBIS.Z))$reset()
    if(sum(model_selection[[2]])>sum(pnoc[[2]])){
      if(plot_progress){
        pairs(X[,mst_var],pch=as.character(y),col=model_selection[[1]],cex=0.5)
      }
      return(list("Cluster Allocation" = model_selection[[1]],
                  "Log Marginal Likelihoods" = model_selection[[2]],
                  "Graph" = model_selection[[3]],
                  "Number of Clusters" = model_selection[[4]],
                  "Edges Removed" = model_selection[[5]]))
    } else{
      if(plot_progress){
        pairs(X[,mst_var],pch=as.character(y),col=pnoc[[1]],cex=0.5)
      }
      return(pnoc)
    }
  }
  if(deforest_criterion=="SoC"){
    psoc <- deforest.soc(obs = X,res = y,gra = g,lbe = logZ,eps = edge_rem,n_dag = min_size,clu_al = z,c_s = combine_save,est_thres = SMC_thres,mtb = BIC_memo_thres,mts = SMC_memo_thres,par_no = N,rfun = rprior,pdf_fun = prior_pdf,efsamp = ess,methas = n_move,p_p = plot_progress,rho = mst_var,vb = verbose)
    get("_cache", envir=environment(memo.bic))$reset()
    get("_cache", envir=environment(IBIS.Z))$reset()
    if(sum(model_selection[[2]])>sum(psoc[[2]])){
      if(plot_progress){
        pairs(X[,mst_var],pch=as.character(y),col=model_selection[[1]],cex=0.5)
      }
      return(list("Cluster Allocation" = model_selection[[1]],
                  "Log Marginal Likelihoods" = model_selection[[2]],
                  "Graph" = model_selection[[3]],
                  "Number of Clusters" = model_selection[[4]],
                  "Edges Removed" = model_selection[[5]]))
    } else{
      if(plot_progress){
        pairs(X[,mst_var],pch=as.character(y),col=psoc[[1]],cex=0.5)
      }
      return(psoc)
    }
  }
  if(deforest_criterion=="Balanced"){
    pbal <- deforest.balanced(obs = X,res = y,gra = g,lbe = logZ,eps = edge_rem,ups = n_min_class,clu_al = z,c_s = combine_save,est_thres = SMC_thres,mtb = BIC_memo_thres,mts = SMC_memo_thres,par_no = N,rfun = rprior,pdf_fun = prior_pdf,efsamp = ess,methas = n_move,p_p = plot_progress,rho = mst_var,vb = verbose)
    get("_cache", envir=environment(memo.bic))$reset()
    get("_cache", envir=environment(IBIS.Z))$reset()
    if(sum(model_selection[[2]])>sum(pbal[[2]])){
      if(plot_progress){
        pairs(X[,mst_var],pch=as.character(y),col=model_selection[[1]],cex=0.5)
      }
      return(list("Cluster Allocation" = model_selection[[1]],
                  "Log Marginal Likelihoods" = model_selection[[2]],
                  "Graph" = model_selection[[3]],
                  "Number of Clusters" = model_selection[[4]],
                  "Edges Removed" = model_selection[[5]]))
    } else{
      if(plot_progress){
        pairs(X[,mst_var],pch=as.character(y),col=pbal[[1]],cex=0.5)
      }
      return(pbal)
    }
  }
  if(deforest_criterion=="None"){
    get("_cache", envir=environment(memo.bic))$reset()
    get("_cache", envir=environment(IBIS.Z))$reset()
    if(plot_progress){
      pairs(X[,mst_var],pch=as.character(y),col=model_selection[[1]],cex=0.5)
    }
    return(list("Cluster Allocation" = z,
                "Log Marginal Likelihoods" = logZ,
                "Graph" = g,
                "Number of Clusters" = K,
                "Edges Removed" = edge_rem))
  }
  if(deforest_criterion=="MaxReg"){
    pmaxreg <- deforest.maxreg(obs = X,res = y,gra = g,lbe = logZ,eps = edge_rem,tau = reg,clu_al = z,c_s = combine_save,est_thres = SMC_thres,mtb = BIC_memo_thres,mts = SMC_memo_thres,par_no = N,rfun = rprior,pdf_fun = prior_pdf,efsamp = ess,methas = n_move,p_p = plot_progress,rho = mst_var,vb = verbose)
    get("_cache", envir=environment(memo.bic))$reset()
    get("_cache", envir=environment(IBIS.Z))$reset()
    if(plot_progress){
      pairs(X[,mst_var],pch=as.character(y),col=pmaxreg[[1]],cex=0.5)
    }
    return(pmaxreg)
  }
  if(deforest_criterion=="Validation"){
    pval <- deforest.validation(obs = X,obs_all = X_all,res = y,res_all = y_all,gra = g,lbe = logZ,eps = edge_rem,gra_all = g_all,clu_al = z,c_s = combine_save,est_thres = SMC_thres,mtb = BIC_memo_thres,mts = SMC_memo_thres,par_no = N,rfun = rprior,pdf_fun = prior_pdf,efsamp = ess,methas = n_move,p_p = plot_progress,rho = mst_var,vb = verbose)
    get("_cache", envir=environment(memo.bic))$reset()
    get("_cache", envir=environment(IBIS.Z))$reset()
    if(plot_progress){
      pairs(X_all[,mst_var],pch=as.character(y_all),col=pval[[2]][[1]],cex=0.5)
    }
    return(pval)
  }
}

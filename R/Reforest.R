###############################################
## R script for functions in UNCOVER package ##
###############################################
##
## Sam Emerson
## September 2022
##


################################
## Pruning Criteria Functions ##
################################


##' Reintroducing edges to a graph to reduce the number of clusters to the maximal amount specified
##'
##'
##' @export
##' @name reforest.noc
##' @description Adds edges to the graph to reduce the number of clusters until the criteria is met and it is no longer beneficial to reintroduce any more edges. This is done in a greedy manner to optimise the log Bayesian evidence of the resulting models.
##'
##' Used in UNCOVER if the reforest condition is set to "NoC".
##'
##' @keywords reforest number
##' @param obs Covariate matrix
##' @param res Binary response vector
##' @param gra `igraph` object which contains the information of the graph of the current model
##' @param lbe A vector detailing the log Bayesian evidences of all the sub-models defined by the separated components of `gra`
##' @param eps A 2-column matrix of edges previously removed. Rows correspond to edges and edges should be expressed as the two vertices the edge connects.
##' @param K_dag The maximum number of clusters allowed in the final output
##' @param clu_al A vector detailing the cluster allocation of each observation. If not specified the function will generate this vector.
##' @param c_s A list of length `nrow(eps)`. See details more information. Does not need to be specified.
##' @param est_method Method to be used to estimate the log Bayesian evidence of the newly formed sub-models, can be one of `"SMC"`,`"SMC_BIC"` or `"BIC"`
##' @param est_thres The threshold for which the number of observations needs to exceed to consider using BIC as an estimator. Only applies if method `"SMC_BIC"` is selected. Defaults to `Inf` if not specified.
##' @param par_no Number of samples of the prior used for the SMC sampler. Not required if method `"BIC"` selected.
##' @param rfun Function to sample from the prior. Must only have two arguments, `p_n` and `di` (Number of prior samples to generate and the number of dimensions of a single sample respectively).
##' @param pdf_fun Probability Density Function of the prior. Must only have two arguments, `th` and `di` (a vector or matrix of regression coefficients samples and the number of dimensions of a single sample respectively).
##' @param p_p Do you want to plot the output of the clustering each time an edge is reintroduced?
##' @param rho Only applies if `p_p` is `TRUE`. A vector specifying which variables of the covariate matrix to plot.
##' @param vb Do you want the progress of the algorithm to be shown?
##' @return A list consisting of; the cluster allocation vector of the new model, the resulting Bayesian evidence vector for the new model, an `igraph` object containing information on the new graph, the number of clusters in the model and the edges that have been removed from the graph to achieve this model.
##' @details Requires a minimum spanning forest graph which defines components for a multiplicative Bayesian logistic regression model, and the edges removed to achieve this graph.
##'
##' For each edge previously removed the log Bayesian evidence `lbe` is estimated for the model obtained by reintroducing said edge. The optimal edge of this group is then selected and if this either showcases an improvement on the current model or the number of clusters in the current model exceeds the maximum number of clusters allowed then we reintroduce this edge to the graph and update the model. If added the remaining edges are then reconsidered and this process is repeated until either; the criterion is met and it is not beneficial to add anymore edges to the graph or there are no longer any edges to reintroduce.
##'
##' If the vertices of graph `gra` are specifically named then elements of the matrix `eps` should not be numeric and instead be the names of the vertices involved.
##'
##' If the clusters specified by the initial model have fixed labels then this should be specified by `clu_al`. `clu_al` must be a one of the possible labelling of the observations defined by the clusters of the graph. For example for a graph where there is only one connected component, if `clu_al` is specified it cannot be anything other than `rep(1,length(res)`.
##'
##' `c_s` allows the user to specify information about the edges removed. For example if `c_s` is specified it must be of the form of list with each element representing information on the reintroduction of an edge. The index of this list corresponds to the index of the edges in `eps`. Furthermore, each element of the list will itself be a list of two elements, the first being the indices of the observations combined by introducing this edge and the second being the sub log Bayesian evidence of the cluster formed through this edge reintroduction. `c_s` is intended to be used to reduce computation time, and so whilst incorrect information on the observations involved in particular lists of `c_s` will not produce incorrect results, it will not have the desired time saving effect.
##'
##' For more details on the specifics of the possible values for `est_method`, see the help page of the function `lbe.gen`.
##' @seealso [lbe.gen]
##' @examples
##'
##' # First we generate a covariate matrix `obs` and binary response vector `res`
##' CM <- matrix(rnorm(20),10,2)
##' rv <- sample(0:1,10,replace=T)
##'
##' # We then generate the minimum spanning tree graph from CM
##' gra_CM <- one.stage.mst(obs = CM)
##' plot(gra_CM,layout=CM)
##'
##' # If we assume the prior for the regression coefficients is a standard normal
##' pr_samp <- function(p_n,di){return(rmvn(p_n,rep(0,di),diag(di)))}
##' pr_fun <- function(th,di){return(dmvn(th,mu=rep(0,di),sigma=diag(di)))}
##'
##' # Then we can generate the log Bayesian evidence for this one component model
##' lZ <- lbe.gen(method = "SMC_BIC",obs_mat = CM,res_vec = rv,p_num = 500, rpri = pr_samp, p_pdf = pr_fun)
##' lZ
##'
##' # If we wished to remove the first edge in the graph we would have
##' er <- remove.edge(gra = gra_CM, j = 1, lbe = lZ, obs = CM, res = rv, est_method = "SMC_BIC", par_no = 500, rfun = pr_samp, pdf_fun = pr_fun)
##' clu_al.2 <- er[[1]]
##' lZ.2 <- er[[2]]
##' gra_CM.2 <- er[[3]]
##' plot(gra_CM.2,layout=CM)
##' lZ.2
##' sum(lZ.2)
##'
##' # If we then wanted to remove the first edge of the newly formed graph gra_CM.2 we would have
##' er.2 <- remove.edge(gra = gra_CM.2, j = 1, clu_al = clu_al.2, lbe = lZ.2, obs = CM, res = rv, est_method = "SMC_BIC", par_no = 500, rfun = pr_samp, pdf_fun = pr_fun)
##' plot(er.2[[3]],layout=CM)
##' er.2[[2]]
##' sum(er.2[[2]])

reforest.noc <- function(obs,res,gra,lbe,eps,K_dag,clu_al=NULL,c_s=NULL,est_method="SMC_BIC",
                         est_thres=Inf,par_no=0,rfun=NULL,pdf_fun=NULL,p_p=F,rho=NULL,vb = F){
  K <- igraph::count_components(gra)
  if(is.null(clu_al)){
    clu_al <- igraph::components(gra)$membership
  }
  if(is.null(c_s)){
    c_s <- vector(mode="list",length=1)
  }
  j <- 1
  if(vb){
    message("Assessing the reintroduction of edges which have been removed")
  }
  while(j <= nrow(eps)){
    edge_clu_al <- sort(clu_al[match(eps[j,],V(gra)$name)])
    if(!identical(c_s[[j]][[1]],which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]))){
      c_s[[j]] <- list(which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]),lbe.gen(method = est_method,thres = est_thres,obs_mat = obs[which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]),,drop=F],res_vec = res[which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2])],p_num = par_no,rpri = rfun,p_pdf = pdf_fun))
    }
    if(vb){
      txtProgressBar(min=1,max=nrow(eps)+1,initial=j+1,style=3)
    }
    if(j==nrow(eps)){
      if(vb){
        message("")
      }
      lbe_comb.1 <- sapply(c_s,FUN = function(u){u[[2]]})
      lbe_comb.2 <- apply(eps,MARGIN = 1,FUN = function(u,lZ,ca,cg){sum(lZ[-ca[match(u,V(cg)$name)]])},lZ = lbe,ca = clu_al,cg = gra)
      lbe_comb <- lbe_comb.1 + lbe_comb.2
      if(K>K_dag | any(lbe_comb>sum(lbe))){
        cut_comb <- which.max(lbe_comb)
        edge_clu_al <- sort(clu_al[match(eps[cut_comb,],V(gra)$name)])
        if(vb){
          message("")
          message(blue(paste0("Combining connected components ",edge_clu_al[1]," and ",edge_clu_al[2])))
        }
        clu_al[which(clu_al==edge_clu_al[2])] <- edge_clu_al[1]
        for(k in (edge_clu_al[2]+1):K){
          clu_al[which(clu_al==k)] <- k-1
        }
        lbe[edge_clu_al[1]] <- c_s[[cut_comb]][[2]]
        lbe <- lbe[-edge_clu_al[2]]
        gra <- igraph::add_edges(gra,eps[cut_comb,])
        eps <- eps[-cut_comb,,drop=F]
        c_s[[cut_comb]] <- c()
        K <- K-1
        if(p_p){
          if(is.null(rho)){
            rho <- 1:ncol(obs)
          }
          pairs(obs[,rho],pch=as.character(res),col=clu_al,cex=0.5)
        }
        j <- 0
        if(nrow(eps)!=0 & vb){
          message("")
          message(green("Reassessing the reintroduction of edges which have been removed"))
        }
      }
    }
    j <- j+1
  }
  return(list("Cluster Allocation" = clu_al,
              "Log Marginal Likelihoods" = lbe,
              "Graph" = gra,
              "Number of Clusters" = K,
              "Edges Removed" = eps))
}

##' Reintroducing edges to a graph in order to ensure the size of all clusters is above the specified threshold
##'
##'
##' @export
##' @name reforest.soc
##' @description Adds edges to the graph to increase the minimum number of observations in models clusters, until the criteria is met and it is no longer beneficial to reintroduce any more edges. This is done in a greedy manner to optimise the log Bayesian evidence of the resulting models.
##'
##' Used in UNCOVER if the reforest condition is set to "SoC".
##'
##' @keywords reforest size
##' @param obs Covariate matrix
##' @param res Binary response vector
##' @param gra `igraph` object which contains the information of the graph of the current model
##' @param lbe A vector detailing the log Bayesian evidences of all the sub-models defined by the separated components of `gra`
##' @param eps A 2-column matrix of edges previously removed. Rows correspond to edges and edges should be expressed as the two vertices the edge connects.
##' @param n_dag The minimum number of observations allowed for any cluster in the final model
##' @param clu_al A vector detailing the cluster allocation of each observation. If not specified the function will generate this vector.
##' @param c_s A list of length `nrow(eps)`. See details more information. Does not need to be specified.
##' @param est_method Method to be used to estimate the log Bayesian evidence of the newly formed sub-models, can be one of `"SMC"`,`"SMC_BIC"` or `"BIC"`
##' @param est_thres The threshold for which the number of observations needs to exceed to consider using BIC as an estimator. Only applies if method `"SMC_BIC"` is selected. Defaults to `Inf` if not specified.
##' @param par_no Number of samples of the prior used for the SMC sampler. Not required if method `"BIC"` selected.
##' @param rfun Function to sample from the prior. Must only have two arguments, `p_n` and `di` (Number of prior samples to generate and the number of dimensions of a single sample respectively).
##' @param pdf_fun Probability Density Function of the prior. Must only have two arguments, `th` and `di` (a vector or matrix of regression coefficients samples and the number of dimensions of a single sample respectively).
##' @param p_p Do you want to plot the output of the clustering each time an edge is reintroduced?
##' @param rho Only applies if `p_p` is `TRUE`. A vector specifying which variables of the covariate matrix to plot.
##' @param vb Do you want the progress of the algorithm to be shown?
##' @return A list consisting of; the cluster allocation vector of the new model, the resulting Bayesian evidence vector for the new model, an `igraph` object containing information on the new graph, the number of clusters in the model and the edges that have been removed from the graph to achieve this model.
##' @details Requires a minimum spanning forest graph which defines components for a multiplicative Bayesian logistic regression model, and the edges removed to achieve this graph.
##'
##' For each edge previously removed the log Bayesian evidence `lbe` is estimated for the model obtained by reintroducing said edge. The optimal edge of this group is then selected and if this either showcases an improvement on the current model or the minimum size of all clusters in the current model is below the minimum size of cluster allowed then we reintroduce this edge to the graph and update the model. If added the remaining edges are then reconsidered and this process is repeated until either; the criterion is met and it is not beneficial to add anymore edges to the graph or there are no longer any edges to reintroduce. We highlight here that if an edge is reintroduced and it does not benefit the overall model then at least one cluster which does not meet the criteria must be being combined.
##'
##' If the vertices of graph `g` are specifically named then elements of the matrix `eps` should not be numeric and instead be the names of the vertices involved.
##'
##' If the clusters specified by the initial model have fixed labels then this should be specified by `clu_al`. `clu_al` must be a one of the possible labelling of the observations defined by the clusters of the graph. For example for a graph where there is only one connected component, if `clu_al` is specified it cannot be anything other than `rep(1,length(res)`.
##'
##' `c_s` allows the user to specify information about the edges removed. For example if `c_s` is specified it must be of the form of list with each element representing information on the reintroduction of an edge. The index of this list corresponds to the index of the edges in `eps`. Furthermore, each element of the list will itself be a list of two elements, the first being the indices of the observations combined by introducing this edge and the second being the sub log Bayesian evidence of the cluster formed through this edge reintroduction. `c_s` is intended to be used to reduce computation time, and so whilst incorrect information on the observations involved in particular lists of `c_s` will not produce incorrect results, it will not have the desired time saving effect.
##'
##' For more details on the specifics of the possible values for `est_method`, see the help page of the function `lbe.gen`.
##' @seealso [lbe.gen]
##' @examples
##'

reforest.soc <- function(obs,res,gra,lbe,eps,n_dag,clu_al=NULL,c_s=NULL,est_method="SMC_BIC",
                         est_thres=Inf,par_no=0,rfun=NULL,pdf_fun=NULL,p_p=F,rho=NULL,vb = F){
  K <- igraph::count_components(gra)
  if(is.null(clu_al)){
    clu_al <- igraph::components(gra)$membership
  }
  if(is.null(c_s)){
    c_s <- vector(mode="list",length=1)
  }
  j <- 1
  if(vb){
    message("Assessing the reintroduction of edges which have been removed")
  }
  while(j <= nrow(eps)){
    edge_clu_al <- sort(clu_al[match(eps[j,],V(gra)$name)])
    if(!identical(c_s[[j]][[1]],which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]))){
      c_s[[j]] <- list(which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]),lbe.gen(method = est_method,thres = est_thres,obs_mat = obs[which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]),,drop=F],res_vec = res[which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2])],p_num = par_no,rpri = rfun,p_pdf = pdf_fun))
    }
    if(vb){
      txtProgressBar(min=1,max=nrow(eps)+1,initial=j+1,style=3)
    }
    if(j==nrow(eps)){
      if(vb){
        message("")
      }
      lbe_comb.1 <- sapply(c_s,FUN = function(u){u[[2]]})
      lbe_comb.2 <- apply(eps,MARGIN = 1,FUN = function(u,lZ,ca,cg){sum(lZ[-ca[match(u,V(cg)$name)]])},lZ = lbe, ca = clu_al,cg = gra)
      lbe_comb <- lbe_comb.1 + lbe_comb.2
      small_clu <- which(table(clu_al)<n_dag)
      small_bool <- apply(eps,MARGIN = 1,FUN = function(u,ca,sc,cg){ca[match(u,V(cg)$name)[1]]%in%sc | ca[match(u,V(cg)$name)[2]]%in%sc},ca=clu_al,sc=small_clu,cg=gra)
      small_bool <- small_bool | lbe_comb>sum(lbe)
      if(any(small_bool)){
        lbe_comb[which(!small_bool)] <- -Inf
        cut_comb <- which.max(lbe_comb)
        edge_clu_al <- sort(clu_al[match(eps[cut_comb,],V(gra)$name)])
        if(vb){
          message("")
          message(blue(paste0("Combining connected components ",edge_clu_al[1]," and ",edge_clu_al[2])))
        }
        clu_al[which(clu_al==edge_clu_al[2])] <- edge_clu_al[1]
        for(k in (edge_clu_al[2]+1):K){
          clu_al[which(clu_al==k)] <- k-1
        }
        lbe[edge_clu_al[1]] <- c_s[[cut_comb]][[2]]
        lbe <- lbe[-edge_clu_al[2]]
        gra <- igraph::add_edges(gra,eps[cut_comb,])
        eps <- eps[-cut_comb,,drop=F]
        c_s[[cut_comb]] <- c()
        K <- K-1
        if(p_p){
          if(is.null(rho)){
            rho <- 1:ncol(obs)
          }
          pairs(obs[,rho],pch=as.character(res),col=clu_al,cex=0.5)
        }
        j <- 0
        if(nrow(eps)!=0 & vb){
          message("")
          message(green("Reassessing the reintroduction of edges which have been removed"))
        }
      }
    }
    j <- j+1
  }
  return(list("Cluster Allocation" = clu_al,
              "Log Marginal Likelihoods" = lbe,
              "Graph" = gra,
              "Number of Clusters" = K,
              "Edges Removed" = eps))
}

##' Reintroducing edges to a graph that decrease the Bayesian evidence within a specified tolerance
##'
##'
##' @export
##' @name reforest.maxreg
##' @description Selects the previously removed edge with the highest Bayesian evidence, and adds this edge to the graph if this edges' Bayesian evidence multiplied by a tolerance parameter (the natural logarithm of which is specified by the user) is larger than the current models Bayesian evidence. This process then repeats until it is no longer beneficial (with consideration to the tolerance parameter) to add an edge to the graph.
##'
##' Used in UNCOVER if the reforest condition is set to "MaxReg".
##'
##' @keywords reforest maximal regret
##' @param obs Covariate matrix
##' @param res Binary response vector
##' @param gra `igraph` object which contains the information of the graph of the current model
##' @param lbe A vector detailing the log Bayesian evidences of all the sub-models defined by the separated components of `gra`
##' @param eps A 2-column matrix of edges previously removed. Rows correspond to edges and edges should be expressed as the two vertices the edge connects.
##' @param tau Numerical natural logarithm of the tolerance parameter. Must be positive.
##' @param clu_al A vector detailing the cluster allocation of each observation. If not specified the function will generate this vector.
##' @param c_s A list of length `nrow(eps)`. See details more information. Does not need to be specified.
##' @param est_method Method to be used to estimate the log Bayesian evidence of the newly formed sub-models, can be one of `"SMC"`,`"SMC_BIC"` or `"BIC"`
##' @param est_thres The threshold for which the number of observations needs to exceed to consider using BIC as an estimator. Only applies if method `"SMC_BIC"` is selected. Defaults to `Inf` if not specified.
##' @param par_no Number of samples of the prior used for the SMC sampler. Not required if method `"BIC"` selected.
##' @param rfun Function to sample from the prior. Must only have two arguments, `p_n` and `di` (Number of prior samples to generate and the number of dimensions of a single sample respectively).
##' @param pdf_fun Probability Density Function of the prior. Must only have two arguments, `th` and `di` (a vector or matrix of regression coefficients samples and the number of dimensions of a single sample respectively).
##' @param p_p Do you want to plot the output of the clustering each time an edge is reintroduced?
##' @param rho Only applies if `p_p` is `TRUE`. A vector specifying which variables of the covariate matrix to plot.
##' @param vb Do you want the progress of the algorithm to be shown?
##' @return A list consisting of; the cluster allocation vector of the new model, the resulting Bayesian evidence vector for the new model, an `igraph` object containing information on the new graph, the number of clusters in the model and the edges that have been removed from the graph to achieve this model.
##' @details Requires a minimum spanning forest graph which defines components for a multiplicative Bayesian logistic regression model, and the edges removed to achieve this graph.
##'
##' The tolerance parameter `exp(tau)` should be interpreted as the maximum the user is willing to regret by reintroducing an edge and decreasing the Bayesian evidence, similar to the concept of minimum improvement. The reintroduction of edges to the graph which are beneficial to the model are naturally accepted as well.
##'
##' For each edge previously removed the log Bayesian evidence `lbe` is estimated for the model obtained by reintroducing said edge. The optimal edge of this group is then selected and if this log Bayesian evidence plus `tau` is greater than `lbe` of the current model then we reintroduce this edge to the graph and update the model. If added the remaining edges are then reconsidered and this process is repeated until either; it is not beneficial (considering `tau`) to add anymore edges to the graph or there are no longer any edges to reintroduce.
##'
##' If the vertices of graph `gra` are specifically named then elements of the matrix `eps` should not be numeric and instead be the names of the vertices involved.
##'
##' If the clusters specified by the initial model have fixed labels then this should be specified by `clu_al`. `clu_al` must be a one of the possible labelling of the observations defined by the clusters of the graph. For example for a graph where there is only one connected component, if `clu_al` is specified it cannot be anything other than `rep(1,length(res)`.
##'
##' `c_s` allows the user to specify information about the edges removed. For example if `c_s` is specified it must be of the form of list with each element representing information on the reintroduction of an edge. The index of this list corresponds to the index of the edges in `eps`. Furthermore, each element of the list will itself be a list of two elements, the first being the indices of the observations combined by introducing this edge and the second being the sub log Bayesian evidence of the cluster formed through this edge reintroduction. `c_s` is intended to be used to reduce computation time, and so whilst incorrect information on the observations involved in particular lists of `c_s` will not produce incorrect results, it will not have the desired time saving effect.
##'
##' For more details on the specifics of the possible values for `est_method`, see the help page of the function `lbe.gen`.
##' @seealso [lbe.gen]
##' @examples
##'

reforest.maxreg <- function(obs,res,gra,lbe,eps,tau,clu_al=NULL,c_s=NULL,est_method="SMC_BIC",
                            est_thres=Inf,par_no=0,rfun=NULL,pdf_fun=NULL,p_p=F,rho=NULL,vb = F){
  K <- igraph::count_components(gra)
  if(is.null(clu_al)){
    clu_al <- igraph::components(gra)$membership
  }
  if(is.null(c_s)){
    c_s <- vector(mode="list",length=1)
  }
  j <- 1
  if(vb){
    message("Assessing the reintroduction of edges which have been removed")
  }
  while(j <= nrow(eps)){
    edge_clu_al <- sort(clu_al[match(eps[j,],V(gra)$name)])
    if(!identical(c_s[[j]][[1]],which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]))){
      c_s[[j]] <- list(which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]),lbe.gen(method = est_method,thres = est_thres,obs_mat = obs[which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]),,drop=F],res_vec = res[which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2])],p_num = par_no,rpri = rfun,p_pdf = pdf_fun))
    }
    if(vb){
      txtProgressBar(min=1,max=nrow(eps)+1,initial=j+1,style=3)
    }
    if(j==nrow(eps)){
      if(vb){
        message("")
      }
      lbe_comb.1 <- sapply(c_s,FUN = function(u){u[[2]]})
      lbe_comb.2 <- apply(eps,MARGIN = 1,FUN = function(u,lZ,ca,cg){sum(lZ[-ca[match(u,V(cg)$name)]])},lZ = lbe,ca = clu_al,cg = gra)
      lbe_comb <- lbe_comb.1 + lbe_comb.2
      if(any(tau + lbe_comb>sum(lbe))){
        cut_comb <- which.max(lbe_comb)
        edge_clu_al <- sort(clu_al[match(eps[cut_comb,],V(gra)$name)])
        if(vb){
          message("")
          message(blue(paste0("Combining connected components ",edge_clu_al[1]," and ",edge_clu_al[2])))
        }
        clu_al[which(clu_al==edge_clu_al[2])] <- edge_clu_al[1]
        for(k in (edge_clu_al[2]+1):K){
          clu_al[which(clu_al==k)] <- k-1
        }
        lbe[edge_clu_al[1]] <- c_s[[cut_comb]][[2]]
        lbe <- lbe[-edge_clu_al[2]]
        gra <- igraph::add_edges(gra,eps[cut_comb,])
        eps <- eps[-cut_comb,,drop=F]
        c_s[[cut_comb]] <- c()
        K <- K-1
        if(p_p){
          if(is.null(rho)){
            rho <- 1:ncol(obs)
          }
          pairs(obs[,rho],pch=as.character(res),col=clu_al,cex=0.5)
        }
        j <- 0
        if(nrow(eps)!=0 & vb){
          message("")
          message(green("Reassessing the reintroduction of edges which have been removed"))
        }
      }
    }
    j <- j+1
  }
  return(list("Cluster Allocation" = clu_al,
              "Log Marginal Likelihoods" = lbe,
              "Graph" = gra,
              "Number of Clusters" = K,
              "Edges Removed" = eps))
}

##' Adding validation data to the training data graph and model, then reintroducing edges that increase a robustness statistic involving the overall model and the training model
##'
##'
##' @export
##' @name reforest.validation
##' @description Adds validation to the training data graph and updates the log Bayesian evidence, then selects the previously removed edge with the highest Bayesian evidence ratio, and adds this edge to the graph if this edges' addition gives a larger Bayesian evidence ratio than the current Bayesian evidence ratio. This process then repeats until it is no longer beneficial to add an edge to the graph.
##'
##' Used in UNCOVER if the reforest condition is set to "Validation".
##'
##' @keywords reforest validation
##' @param obs Covariate matrix of the training data
##' @param obs_all Covariate matrix of training and validation data
##' @param res Binary response vector of the training data
##' @param res_all Binary response vector of training and validation data
##' @param gra `igraph` object which contains the information of the graph of the current model for the training data
##' @param lbe A vector detailing the log Bayesian evidences of all the sub-models defined by the separated components of `gra`
##' @param eps A 2-column matrix of edges previously removed from the graph `gra`. Rows correspond to edges and edges should be expressed as the two vertices the edge connects.
##' @param gra_all `igraph` object which contains the information of the graph of the current model for the all the data. If not specified `tr_ind` and `rho` must be specified. See details.
##' @param tr_ind Only applies if `gra_all` is `NULL`. Vector of indices of the training observations.
##' @param rho Only applies if `gra_all` is `NULL` or `p_p` is `TRUE`. A vector specifying which variables of the covariate matrix were used to construct the graph `gra`.
##' @param clu_al A vector detailing the cluster allocation of each observation. If not specified the function will generate this vector.
##' @param c_s A list of length `nrow(eps)`. See details more information. Does not need to be specified.
##' @param est_method Method to be used to estimate the log Bayesian evidence of the newly formed sub-models, can be one of `"SMC"`,`"SMC_BIC"` or `"BIC"`
##' @param est_thres The threshold for which the number of observations needs to exceed to consider using BIC as an estimator. Only applies if method `"SMC_BIC"` is selected. Defaults to `Inf` if not specified.
##' @param par_no Number of samples of the prior used for the SMC sampler. Not required if method `"BIC"` selected.
##' @param rfun Function to sample from the prior. Must only have two arguments, `p_n` and `di` (Number of prior samples to generate and the number of dimensions of a single sample respectively).
##' @param pdf_fun Probability Density Function of the prior. Must only have two arguments, `th` and `di` (a vector or matrix of regression coefficients samples and the number of dimensions of a single sample respectively).
##' @param p_p Do you want to plot the output of the clustering each time an edge is reintroduced?
##' @param vb Do you want the progress of the algorithm to be shown?
##' @return A list of two lists, one containing information of the training data model and one containing information on the overall model. Each list consists of a further list containing; the cluster allocation vector of the new model, the resulting Bayesian evidence vector for the new model, an `igraph` object containing information on the new graph, the number of clusters in the model and the edges that have been removed from the graph to achieve this model.
##' @details Requires a minimum spanning forest graph which defines components for a multiplicative Bayesian logistic regression model, and the edges removed to achieve this graph. This will correspond to the training model.
##'
##' `gra_all` represents the graph obtained by adding the validation data to the initial minimum spanning tree graph of the training data and then removing the edges removed from the training data graph. See details in `two.stage.mst` for more information. If `gra_all` is not specified but `tr_ind` and `rho` are, then `tr_ind` must be the row indices of `obs_all` that are not present in `obs` and `rho` must be the same variables that constucted the training graph.
##'
##' The names of the vertices in `gra` (and `gra_all` if provided) must correspond to the row index of the observation in `obs_all`.
##'
##' `reforest.validation` first assigns the validation data to a cluster and subsequently updates each sub models log Bayesian evidence, summing to make an updated overall log Bayesian evidence. The difference between the overall log Bayesian evidence for the updated model and the training model is the new robustness statistic (`RobS`), and the edge (if any) whose reintroduction increases `RobS` the most is added back to both graphs (`gra` and `gra_all`) and subsequently both models including `RobS` are updated. If an edge is added the remaining edges are then reconsidered (with respect to `RobS`) and this process is repeated until either; it is not beneficial to add anymore edges to the graph or there are no longer any edges to reintroduce.
##'
##' If the clusters specified by the initial model have fixed labels then this should be specified by `clu_al`. `clu_al` must be a one of the possible labelling of the observations defined by the clusters of the graph. For example for a graph where there is only one connected component, if `clu_al` is specified it cannot be anything other than `rep(1,length(res)`.
##'
##' `c_s` allows the user to specify information about the edges removed. For example if `c_s` is specified it must be of the form of list with each element representing information on the reintroduction of an edge. The index of this list corresponds to the index of the edges in `eps`. Furthermore, each element of the list will itself be a list of two elements, the first being the indices of the observations combined by introducing this edge and the second being the sub log Bayesian evidence of the cluster formed through this edge reintroduction. `c_s` is intended to be used to reduce computation time, and so whilst incorrect information on the observations involved in particular lists of `c_s` will not produce incorrect results, it will not have the desired time saving effect.
##'
##' For more details on the specifics of the possible values for `est_method`, see the help page of the function `lbe.gen`.
##' @seealso [lbe.gen],[two.stage.mst]
##' @examples
##'


reforest.validation <- function(obs,obs_all,res,res_all,gra,lbe,eps,gra_all = NULL,tr_ind=NULL,
                                rho=NULL,clu_al = NULL,c_s=NULL,est_method="SMC_BIC",
                                est_thres=Inf,par_no=0,rfun=NULL,pdf_fun=NULL,p_p=F,vb = F){
  K <- igraph::count_components(gra)
  if(is.null(clu_al)){
    clu_al <- igraph::components(gra)$membership
  }
  if(is.null(c_s)){
    c_s <- vector(mode="list",length=1)
  }
  if(is.null(gra_all)){
    gra_val <- two.stage.mst(obs_mat = obs_all,tr_ind = tr_ind,mst_sub = rho)
    gra_all <- delete_edges(gra_val[[2]],E(gra_val[[2]])[get.edge.ids(gra_val[[2]],eps)])
  }
  clu_al_all <- components(gra_all)$membership
  for(k in 1:K){
    clu_al_all[which(clu_al_all==clu_al_all[as.numeric(V(gra)$name[which(clu_al==k)[1]])])] <- K+k
  }
  clu_al_all <- clu_al_all - K
  lbe_all <- lbe
  for(k in 1:K){
    lbe_all[k] <- lbe.gen(method = est_method,thres = est_thres,obs_mat = obs_all[which(clu_al_all==k),,drop=F],res_vec = res_all[which(clu_al_all==k)],p_num = par_no,rpri = rfun,p_pdf = pdf_fun)
  }
  RobS <- sum(lbe_all-lbe)
  j <- 1
  c_s_all <- vector(mode="list",length=nrow(eps))
  if(vb){
    message("")
    message("Assessing the reintroduction of edges which have been removed")
  }
  while(j <= nrow(eps)){
    edge_clu_al <- sort(clu_al_all[as.numeric(eps[j,])])
    if(!identical(c_s[[j]][[1]],which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]))){
      c_s[[j]] <- list(which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]),lbe.gen(method = est_method,thres = est_thres,obs_mat = obs[which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]),,drop=F],res_vec = res[which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2])],p_num = par_no,rpri = rfun,p_pdf = pdf_fun))
    }
    if(!identical(c_s_all[[j]][[1]],which(clu_al_all==edge_clu_al[1]| clu_al_all==edge_clu_al[2]))){
      c_s_all[[j]] <- list(which(clu_al_all==edge_clu_al[1]| clu_al_all==edge_clu_al[2]),lbe.gen(method = est_method,thres = est_thres,obs_mat = obs_all[which(clu_al_all==edge_clu_al[1]| clu_al_all==edge_clu_al[2]),,drop=F],res_vec = res_all[which(clu_al_all==edge_clu_al[1]| clu_al_all==edge_clu_al[2])],p_num = par_no,rpri = rfun,p_pdf = pdf_fun))
    }
    if(vb){
      txtProgressBar(min=1,max=nrow(eps)+1,initial=j+1,style=3)
    }
    if(j==nrow(eps)){
      lbe_comb.1 <- sapply(c_s,FUN = function(u){u[[2]]})
      lbe_comb.2 <- apply(eps,MARGIN = 1,FUN = function(u,lZ,ca,cg){sum(lZ[-ca[match(u,V(cg)$name)]])},lZ = lbe,ca = clu_al,cg = gra)
      lbe_comb <- lbe_comb.1 + lbe_comb.2
      lbe_comb_all.1 <- sapply(c_s_all,FUN = function(u){u[[2]]})
      lbe_comb_all.2 <- apply(eps,MARGIN = 1,FUN = function(u,lZa,caa){sum(lZa[-caa[as.numeric(u)]])},lZa = lbe_all,caa = clu_al_all)
      lbe_comb_all <- lbe_comb_all.1 + lbe_comb_all.2
      RobS_comb <- lbe_comb_all - lbe_comb
      if(any(RobS_comb>RobS)){
        cut_comb <- which.max(RobS_comb)
        edge_clu_al <- sort(clu_al_all[as.numeric(eps[cut_comb,])])
        if(vb){
          message("")
          message(blue(paste0("Combining connected components ",edge_clu_al[1]," and ",edge_clu_al[2])))
        }
        clu_al[which(clu_al==edge_clu_al[2])] <- edge_clu_al[1]
        for(k in (edge_clu_al[2]+1):K){
          clu_al[which(clu_al==k)] <- k-1
        }
        clu_al_all[which(clu_al_all==edge_clu_al[2])] <- edge_clu_al[1]
        for(k in (edge_clu_al[2]+1):K){
          clu_al_all[which(clu_al_all==k)] <- k-1
        }
        lbe[edge_clu_al[1]] <- c_s[[cut_comb]][[2]]
        lbe <- lbe[-edge_clu_al[2]]
        lbe_all[edge_clu_al[1]] <- c_s_all[[cut_comb]][[2]]
        lbe_all <- lbe_all[-edge_clu_al[2]]
        RobS <- sum(lbe_all-lbe)
        gra <- add_edges(gra,eps[cut_comb,])
        gra_all <- add_edges(gra_all,eps[cut_comb,])
        eps <- eps[-cut_comb,,drop=F]
        c_s[[cut_comb]] <- c()
        c_s_all[[cut_comb]] <- c()
        K <- K-1
        if(p_p){
          if(is.null(rho)){
            rho <- 1:ncol(obs)
          }
          pairs(obs_all[,rho],pch=as.character(res_all),col=clu_al_all,cex=0.5)
        }
        j <- 0
        if(nrow(eps)!=0 & vb){
          message("")
          message(green("Reassessing the reintroduction of edges which have been removed"))
        }
      }
    }
    j <- j+1
  }
  return(list("Training Data" = list("Cluster Allocation" = clu_al,
                                     "Log Marginal Likelihoods" = lbe,
                                     "Graph" = gra,
                                     "Number of Clusters" = K,
                                     "Edges Removed" = eps),
              "All Data" = list("Cluster Allocation" = clu_al_all,
                                "Log Marginal Likelihoods" = lbe_all,
                                "Graph" = gra_all,
                                "Number of Clusters" = K,
                                "Edges Removed" = eps)))
}


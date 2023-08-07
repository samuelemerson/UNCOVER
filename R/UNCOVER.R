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
##' a graphical representation of the co-variates. Edges are removed (or
##' reintroduced) by considering the normalisation constant (or Bayesian
##' evidence) of a multiplicative Bayesian logistic regression model.
##'
##' The first stage of the function is concerned purely with a greedy
##' optimisation of the Bayesian evidence through edge manipulation. The second
##' stage then addresses any other criteria (known as deforestation conditions)
##' expressed by the user through reintroduction of edges.
##'
##' @keywords graph cohort cluster Bayesian evidence
##' @param X Co-variate matrix
##' @param y Binary response vector
##' @param mst_var A vector specifying which variables of the co-variate matrix
##' will be used to form the graph. If not specified all variables will be used.
##' @param options Additional arguments that can be specified for `UNCOVER`.
##' See [UNCOVER.opts()] for details. Can be ignored.
##' @param stop_criterion What is the maximum number of clusters allowed before
##' we terminate the first stage and begin deforestation. Defaults to 5.
##' @param deforest_criterion Constraint type which the final model must satisfy.
##' Can be one of `"NoC"`, `"SoC"`, `"MaxReg"`, `"Validation"`, `"Diverse"` or
##' `"None"`. See details. Defaults to `"None"`.
##' @param prior_mean Mean for the multivariate normal prior used in the SMC
##' sampler. See details. Defaults to the origin.
##' @param prior_var Variance matrix for the multivariate normal prior used in
##' the SMC sampler. See details. Defaults to the identity matrix.
##' @param verbose Do you want the progress of the algorithm to be shown?
##' Defaults to `TRUE`.
##' @return An object of class `"UNCOVER"`, which is a list consisting of:
##'
##' \describe{
##' \item{`Covariate_Matrix`}{The co-variate matrix provided.}
##' \item{`Response_Vector`}{The binary response vector provided.}
##' \item{`Minimum_Spanning_Tree_Variables`}{A vector of indices for the
##' co-variates used to construct the minimum spanning tree.}
##' \item{`Control`}{A list of the additional arguments specified by `options`.}
##' \item{`Deforestation_Criterion`}{The deforestation criterion specified.}
##' \item{`Prior_Mean`}{The mean of multivariate normal prior. Meaningless if
##' prior is overridden in `options`.}
##' \item{`Prior_Variance`}{The variance of multivariate normal prior.
##' Meaningless if prior is overridden in `options`.}
##' \item{`Model`}{List containing; the cluster allocation of the training data,
##' the log Bayesian evidences of the sub-models, the final graph of the
##' clustered data, the number of clusters, the edges which were removed from
##' the graph and a diagnostics data frame (the contents of which vary depending
##' on the deforestation criterion).}
##' }
##'
##' If `deforest_criterion=="Validation"` then `Model` is instead a list of two
##' lists; one containing the model information for the training data
##' (`Training_Data`) and the other containing model information for all of the
##' data (`All_Data`). Diagnostic information is only included in the `All_Data`
##' list.
##'
##' @details Assumes a Bayesian logistic regression model for each cohort, with
##' the overall model being a product of these sub-models.
##'
##' A minimum spanning tree graph is first constructed from a subset of the
##' co-variates. Then at each iteration, each edge in the current graph is
##' checked to see if removal to split a cohort is beneficial, and then either
##' we selected the optimal edge to remove or we conclude it is not beneficial
##' to remove any more edges. At the end of each iteration we also check the set
##' of removed edges to see if it is beneficial to reintroduce any previously
##' removed edges. After this process has ended we then reintroduce edges in the
##' removed set specifically to meet the criteria set by the user in the most
##' optimal manner possible through a greedy approach. For more details see the
##' Emerson and Aslett (2023).
##'
##' The graph can be undergo deforestation to meet 6 possible criteria:
##'
##' \enumerate{
##'  \item `"NoC"`: Number of Clusters - we specify a maximum number of clusters
##' (`options$max_K`) we can tolerate in the final output of the algorithm.
##' \item `"SoC"`: Size of Clusters - we specify a minimum number of
##' observations (`options$min_size`) we can tolerate being assigned to a
##' cluster in the final output of the algorithm.
##' \item `"MaxReg"`: Maximal Regret - we give a maximum tolerance
##' (`exp(options$reg)`) that we allow the Bayesian evidence to decrease by
##' reintroducing an edge.
##' \item `"Validation"`: Validation Data - we split (using
##' `options$train_frac`) the data into training and validation data, apply the
##' first stage of the algorithm on the training data and the introduce the
##' validation data for the deforestation stage. Edges are reintroduced if they
##' lead to improved prediction of the validation data using the training data
##' model (i.e. we aim to maximise a robustness statistic).
##' \item `"Diverse"`: Diverse Response Class Within Clusters - We specify a
##' minimum number of observations (`options$n_min_class`) in a cluster that
##' have the minority response class associated to them (the minimum response
##' class is determined for each cluster).
##' \item `"None"`: No Criteria Specified - we do not go through the second
##' deforestation stage of the algorithm.
##' }
##'
##' All deforestation criteria other than `"None"` require additional arguments
##' to be specified in `options`. See examples and
##' [UNCOVER.opts()] for more information. It is never
##' recommended to use anything other than
##' `UNCOVER.opts` to provide the `options` argument.
##'
##' The prior used for the UNCOVER procedure will take the form of a
##' multivariate normal, where the parameters can be specified directly by the
##' user. It is however possible to override this default prior distributional
##' form by specifying `prior.override=TRUE` and providing the relevant prior
##' functions in `UNCOVER.opts`.
##'
##' The diagnostic data frames will track various outputs of the UNCOVER
##' procedure depending on the deforestation criterion. All data frames will
##' contain an action (removal or addition of an edge to the graph) and the
##' total log Bayesian evidence of the model gained through deployment of that
##' action (for `"Validation"` this will be two columns, one for the training
##' data and one for all of the data). `"NoC"` will also track the number of
##' clusters, `"SoC"` will track the minimum cluster size and the number of
##' criterion breaking clusters, `"Validation"` will track the robustness
##' statistic and "`Diverse"` will track the minimum minority class across all
##' clusters alongside the number of criterion breaking clusters.
##'
##' @seealso [UNCOVER.opts()], [print.UNCOVER()], [predict.UNCOVER()], [plot.UNCOVER()]
##' @references \itemize{
##' \item Emerson, S.R. and Aslett, L.J.M. (2023). Joint cohort and prediction
##' modelling through graphical structure analysis (to be released)
##' }
##' @examples
##'
##' \donttest{
##' # First we generate a co-variate matrix and binary response vector
##' CM <- matrix(rnorm(200),100,2)
##' rv <- sample(0:1,100,replace=TRUE)
##'
##' # We can then run our algorithm to see what cohorts are selected for each
##' # of the different deforestation criteria
##' UN.none <- UNCOVER(X = CM,y = rv, deforest_criterion = "None",
##'                    verbose = FALSE)
##' UN.noc <- UNCOVER(X = CM,y = rv, deforest_criterion = "NoC",
##'                   options = UNCOVER.opts(max_K = 3), verbose = FALSE)
##' UN.soc <- UNCOVER(X = CM,y = rv, deforest_criterion = "SoC",
##'                   options = UNCOVER.opts(min_size = 10), verbose = FALSE)
##' UN.maxreg <- UNCOVER(X = CM,y = rv, deforest_criterion = "MaxReg",
##'                      options = UNCOVER.opts(reg = 1), verbose = FALSE)
##' UN.validation <- UNCOVER(X = CM,y = rv, deforest_criterion = "Validation",
##'                          options = UNCOVER.opts(train_frac = 0.8),
##'                          verbose = FALSE)
##' UN.diverse <- UNCOVER(X = CM,y = rv, deforest_criterion = "Diverse",
##'                        options = UNCOVER.opts(n_min_class = 2),
##'                        verbose = FALSE)
##' clu_al_mat <- rbind(UN.none$Model$Cluster_Allocation,
##'                     UN.noc$Model$Cluster_Allocation,
##'                     UN.soc$Model$Cluster_Allocation,
##'                     UN.maxreg$Model$Cluster_Allocation,
##'                     UN.validation$Model$All_Data$Cluster_Allocation,
##'                     UN.diverse$Model$Cluster_Allocation)
##' # We can create a matrix where each entry shows in how many of the methods
##' # did the indexed observations belong to the same cluster
##' obs_con_mat <- matrix(0,100,100)
##' for(i in 1:100){
##' for(j in 1:100){
##' obs_con_mat[i,j] <- length(which(clu_al_mat[,i]-clu_al_mat[,j]==0))/6
##' obs_con_mat[j,i] <- obs_con_mat[i,j]
##' }
##' }
##' head(obs_con_mat)
##'
##' # We can also view the outputted overall Bayesian evidence of the five
##' # models as well
##' c(sum(UN.none$Model$Log_Marginal_Likelihoods),
##'   sum(UN.noc$Model$Log_Marginal_Likelihoods),
##'   sum(UN.soc$Model$Log_Marginal_Likelihoods),
##'   sum(UN.maxreg$Model$Log_Marginal_Likelihoods),
##'   sum(UN.validation$Model$All_Data$Log_Marginal_Likelihoods),
##'   sum(UN.diverse$Model$Log_Marginal_Likelihoods))
##'
##' # If we don't assume the prior for the regression coefficients is a
##' # standard multivariate normal but instead a multivariate normal with
##' # different parameters
##' UN.none.2 <- UNCOVER(X = CM,y = rv, deforest_criterion = "None",
##'                      prior_mean = rep(1,3), prior_var = 0.5*diag(3),
##'                      verbose = FALSE)
##' c(sum(UN.none$Model$Log_Marginal_Likelihoods),
##'   sum(UN.none.2$Model$Log_Marginal_Likelihoods))
##' # We can also specify a completely different prior, for example a
##' # multivariate independent uniform
##' rmviu <- function(n,a,b){
##' return(mapply(FUN = function(min.vec,max.vec,pn){
##'                       stats::runif(pn,a,b)},min.vec=a,max.vec=b,
##'                                      MoreArgs = list(pn = n)))
##' }
##' dmviu <- function(x,a,b){
##' for(ii in 1:ncol(x)){
##'   x[,ii] <- dunif(x[,ii],a[ii],b[ii])
##' }
##' return(apply(x,1,prod))
##' }
##' UN.none.3 <- UNCOVER(X = CM,y = rv,deforest_criterion = "None",
##'                      options = UNCOVER.opts(prior.override = TRUE,
##'                                             rprior = rmviu,
##'                                             dprior = dmviu,a=rep(0,3),
##'                                             b=rep(1,3)),verbose = FALSE)
##' c(sum(UN.none$Model$Log_Marginal_Likelihoods),
##'   sum(UN.none.2$Model$Log_Marginal_Likelihoods),
##'   sum(UN.none.3$Model$Log_Marginal_Likelihoods))
##'
##' # We may only wish to construct our minimum spanning tree based on the first
##' # variable
##' UN.none.4 <- UNCOVER(X = CM,y = rv,mst_var = 1,deforest_criterion = "None",
##'                      verbose = FALSE)
##' c(sum(UN.none$Model$Log_Marginal_Likelihoods),
##'   sum(UN.none.4$Model$Log_Marginal_Likelihoods))
##'
##' # Increasing the stop criterion may uncover more clustering structure within
##' # the data, but comes with a time cost
##' system.time(UNCOVER(X = CM,y = rv,stop_criterion = 4,verbose = FALSE))
##' system.time(UNCOVER(X = CM,y = rv,stop_criterion = 6,verbose = FALSE))
##' }
##'



UNCOVER <- function(X,y,mst_var=NULL,options = UNCOVER.opts(),stop_criterion=5,
                    deforest_criterion="None",prior_mean=rep(0,ncol(X)+1),
                    prior_var=diag(ncol(X)+1),verbose = TRUE){
  memo.bic <- memoise::memoise(memo.bic,
                               cache = options$BIC_cache,
                               omit_args = c("param_start"))
  IBIS.Z <- memoise::memoise(IBIS.Z,cache = options$SMC_cache,
                              omit_args = c("sampl","current_set"))
  X_names <- colnames(X)
  X <- model.matrix(~.,data = as.data.frame(X))[,-1]
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
    prior_mean = "overridden"
    prior_var = "overridden"
  } else{
    rprior <- mvnfast::rmvn
    dprior <- mvnfast::dmvn
    MoreArgs <- list(mu = prior_mean,sigma = prior_var)
  }
  if(is.null(mst_var)){
    mst_var <- 1:ncol(X)
  }
  if(deforest_criterion=="Validation"){
    if(options$train_frac==1){
      stop("train_frac must be specified as less than 1 if Validation criterion selected")
    }
    samp <- sort(sample((1:length(y)),round(length(y)*options$train_frac)))
    g_val <- two.stage.mst(obs_mat = X,tr_ind = samp,mst_sub = mst_var)
    g <- g_val[[1]]
    g_all <- g_val[[2]]
    beid <- g_val[[3]]
    X_all <- X
    y_all <- y
    X <- X[samp,,drop=FALSE]
    y <- y[samp]
  } else{
    osm <- one.stage.mst(obs = X,rho = mst_var)
    g <- osm[[1]]
    beid <- osm[[2]]
  }
  depth_g <- igraph::V(g)$name[igraph::dfs(g,igraph::get_diameter(g)[1])$order]
  edge_rank <- matrix(0,length(igraph::E(g)),2)
  for(i in 2:length(depth_g)){
    pot_edg <- igraph::V(g)$name[igraph::neighbors(g,depth_g[i])]
    edge_rank[i-1,] <- c(depth_g[i],pot_edg[which(pot_edg %in% depth_g[1:i])])
  }
  n <- length(y)
  K <- 1
  z <- rep(1,n)
  logZ <- lbe.gen(thres = options$SMC_thres,obs_mat = X,res_vec = y,
                  memo_thres_bic = options$BIC_memo_thres,
                  memo_thres_smc = options$SMC_memo_thres,
                  p_num = options$N,rpri = rprior,p_pdf = dprior,MA = MoreArgs,
                  efs = options$ess,nm = options$n_move,
                  cache_bic = options$BIC_cache,cache_smc = options$SMC_cache,
                  SMC_fun = IBIS.Z,BIC_fun = memo.bic,
                  ribis_thres = options$RIBIS_thres)
  if(options$diagnostics){
    n_min_size <- sum(table(z)<options$min_size)
    num_min_cla <- sum(table(y)<options$n_min_class)
    Track <- data.frame(Action = "Init.",Log_Bayesian_Evidence = logZ,
                        Number_Of_Clusters = K,Minimum_Cluster_Size = n,
                        Number_Of_Small_Clusters = n_min_size,
                        Minimum_Minority_Class = min(table(y)),
                        Number_Of_Undiverse_Clusters = num_min_cla)
  }
  edge_rem <- c()
  system_save <- vector(mode = "list",length=1)
  combine_save <- list()
  if(deforest_criterion=="NoC" | deforest_criterion=="SoC" | deforest_criterion=="Diverse"){
    model_selection <- list(z,logZ,g,K,edge_rem)
  }
  m <- 1
  while(TRUE){
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
        message(crayon::green(paste0("Assessing the removal of edges from connected component ",k)))
      }
      if(deforest_criterion=="Validation"){
        system_save[[k]] <- list(z,logZ,g,c(),g_all)
      } else{
        system_save[[k]] <- list(z,logZ,g,c())
      }
      for(q in 1:nrow(edge_rank)){
        i <- igraph::get.edge.ids(g,edge_rank[q,])
        if(i==0){
          if(verbose){
            utils::txtProgressBar(min=1,max=nrow(edge_rank),initial=q,style=3)
          }
          next
        }
        if(z[igraph::get.edgelist(g,names=FALSE)[i,]][1]!=k | E(g)[igraph::get.edge.ids(g,igraph::get.edgelist(g)[i,])]$weight<beid){
          if(verbose){
            utils::txtProgressBar(min=1,max=nrow(edge_rank),initial=q,style=3)
          }
          next
        }
        if(deforest_criterion=="Validation"){
          g_temp <- igraph::delete_edges(g_all,igraph::E(g_all)[igraph::get.edge.ids(g_all,igraph::get.edgelist(g)[i,])])
          z_temp <- igraph::components(g_temp)$membership
          if(!identical(sort(unique(z_temp[-samp])),as.numeric(1:(K+1)))){
            if(verbose){
              utils::txtProgressBar(min=1,max=nrow(edge_rank),initial=q,style=3)
            }
            next
          }
        }
        er_temp <- remove.edge(gra = g,j = i,clu_al = z,lbe = logZ,obs = X,
                               res = y,est_thres = options$SMC_thres,
                               mtb = options$BIC_memo_thres,
                               mts = options$SMC_memo_thres,
                               par_no = options$N,rfun = rprior,
                               pdf_fun = dprior,efsamp = options$ess,
                               methas = options$n_move,cb = options$BIC_cache,
                               cs = options$SMC_cache,PA = MoreArgs,
                               SMC_f = IBIS.Z,BIC_f = memo.bic,
                               rt = options$RIBIS_thres)
        if(sum(er_temp[[2]])>sum(system_save[[k]][[2]])){
          if(deforest_criterion=="Validation"){
            system_save[[k]] <- c(er_temp,list(igraph::get.edgelist(g)[i,]),
                                  list(g_temp))
          } else{
            system_save[[k]] <- c(er_temp,list(igraph::get.edgelist(g)[i,]))
          }
        }
        if(verbose){
          utils::txtProgressBar(min=1,max=nrow(edge_rank),initial=q,style=3)
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
      message(crayon::blue(paste0("Removing edge from connected component ",k_change)))
    }
    combine_save[[K]] <- list(which(z==k_change),logZ[k_change])
    z <- system_save[[k_change]][[1]]
    logZ <- system_save[[k_change]][[2]]
    g <- system_save[[k_change]][[3]]
    edge_rem <- rbind(edge_rem,system_save[[k_change]][[4]])
    if(deforest_criterion=="Validation"){
      g_all <- system_save[[k_change]][[5]]
    }
    if(options$diagnostics){
      tab_z <- table(z)
      sap_div <- sapply(1:(K+1),FUN = function(u,res,clu_al){min(table(factor(res[which(clu_al==u)],levels=0:1)))},res=y,clu_al=z)
      Track <- rbind(Track,
                     c(paste0("Rem.",system_save[[k_change]][[4]][1],"-",system_save[[k_change]][[4]][2]),
                       sum(logZ),K+1,min(tab_z),sum(tab_z<options$min_size),
                       min(sap_div),sum(sap_div<options$n_min_class)))
    }
    for(k in 1:K){
      if(k==k_change){
        system_save[[k]] <- vector(mode = "list",length=1)
      } else{
        system_save[[k]][[1]][which(system_save[[k]][[1]]==(K+1))] <- K+2
        system_save[[k]][[1]][which(z==(K+1))] <- K+1
        system_save[[k]][[2]] <- c(system_save[[k]][[2]],
                                   system_save[[k]][[2]][K+1])
        system_save[[k]][[2]][K+1] <- logZ[K+1]
        system_save[[k]][[2]][k_change] <- logZ[k_change]
        system_save[[k]][[3]] <- igraph::delete_edges(system_save[[k]][[3]],igraph::E(system_save[[k]][[3]])[igraph::get.edge.ids(system_save[[k]][[3]],edge_rem[nrow(edge_rem),])])
        if(deforest_criterion=="Validation"){
          system_save[[k]][[5]] <- igraph::delete_edges(system_save[[k]][[5]],igraph::E(system_save[[k]][[5]])[igraph::get.edge.ids(system_save[[k]][[5]],edge_rem[nrow(edge_rem),])])
        }
      }
    }
    K <- K+1
    system_save[[K]] <- vector(mode = "list",length=1)
    if(deforest_criterion=="NoC" & K<=options$max_K){
      if(sum(logZ)>sum(model_selection[[2]])){
        model_selection <- list(z,logZ,g,K,edge_rem)
      }
    }
    if(deforest_criterion=="SoC" & all(table(z)>=options$min_size)){
      if(sum(logZ)>sum(model_selection[[2]])){
        model_selection <- list(z,logZ,g,K,edge_rem)
      }
    }
    if(deforest_criterion=="Diverse"){
      undiv_clu <- sapply(1:K,FUN = function(u,res,clu_al,ups){all(table(factor(res[which(clu_al==u)],levels=0:1))>=ups)},res=y,clu_al=z,ups=options$n_min_class)
      if(all(undiv_clu) & sum(logZ)>sum(model_selection[[2]])){
        model_selection <- list(z,logZ,g,K,edge_rem)
      }
    }
    j <- 1
    if(verbose){
      message("")
      message(crayon::green("Assessing the reintroduction of edges which have been removed"))
    }
    while(j <= nrow(edge_rem)){
      edge_z <- sort(z[match(edge_rem[j,],igraph::V(g)$name)])
      if(!identical(combine_save[[j]][[1]],which(z==edge_z[1]| z==edge_z[2]))){
        combine_save[[j]] <- list(which(z==edge_z[1]| z==edge_z[2]),
                                  lbe.gen(thres = options$SMC_thres,obs_mat = X,
                                          res_vec = y,
                                          obs_ind = which(z==edge_z[1]| z==edge_z[2]),
                                          p_num = options$N,rpri = rprior,
                                          p_pdf = dprior,
                                          memo_thres_bic = options$BIC_memo_thres,
                                          memo_thres_smc = options$SMC_memo_thres,
                                          efs = options$ess,nm = options$n_move,
                                          cache_bic = options$BIC_cache,
                                          cache_smc = options$SMC_cache,
                                          MA = MoreArgs,SMC_fun = IBIS.Z,
                                          BIC_fun = memo.bic,
                                          ribis_thres = options$RIBIS_thres))
      }
      if(verbose){
        utils::txtProgressBar(min=1,max=nrow(edge_rem)+1,initial=j+1,style=3)
      }
      if(j==nrow(edge_rem)){
        logZ_comb.1 <- sapply(combine_save,FUN = function(u){u[[2]]})
        logZ_comb.2 <- apply(edge_rem,MARGIN = 1,FUN = function(u,lZ,ca,cg){sum(lZ[-ca[match(u,igraph::V(cg)$name)]])},lZ = logZ,ca = z, cg = g)
        logZ_comb <- logZ_comb.1 + logZ_comb.2
        if(any(logZ_comb>sum(logZ))){
          cut_comb <- which.max(logZ_comb)
          edge_z <- sort(z[match(edge_rem[cut_comb,],igraph::V(g)$name)])
          if(verbose){
            message("")
            message(crayon::blue(paste0("Combining connected components ",edge_z[1]," and ",edge_z[2])))
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
            system_save[[k]][[3]] <- igraph::add_edges(system_save[[k]][[3]],
                                               edge_rem[cut_comb,],weight = beid)
            if(deforest_criterion=="Validation"){
              system_save[[k]][[5]] <- igraph::add_edges(system_save[[k]][[5]],
                                                 edge_rem[cut_comb,])
            }
          }
          logZ[edge_z[1]] <- combine_save[[cut_comb]][[2]]
          logZ <- logZ[-edge_z[2]]
          g <- igraph::add_edges(g,edge_rem[cut_comb,],weight = beid)
          if(deforest_criterion=="Validation"){
            g_all <- igraph::add_edges(g_all,edge_rem[cut_comb,])
          }
          if(options$diagnostics){
            tab_z <- table(z)
            sap_div <- sapply(1:(K-1),FUN = function(u,res,clu_al){min(table(factor(res[which(clu_al==u)],levels=0:1)))},res=y,clu_al=z)
            Track <- rbind(Track,c(paste0("Add.",edge_rem[cut_comb,1],"-",edge_rem[cut_comb,2]),
                                   sum(logZ),K-1,min(tab_z),sum(tab_z<options$min_size),
                                   min(sap_div),sum(sap_div<options$n_min_class)))
          }
          edge_rem <- edge_rem[-cut_comb,,drop=FALSE]
          combine_save[[cut_comb]] <- c()
          K <- K-1
          j <- 0
          if(verbose){
            message("")
            message(crayon::green("Reassessing the reintroduction of edges which have been removed"))
          }
          if(deforest_criterion=="NoC" & K<=options$max_K){
            if(sum(logZ)>sum(model_selection[[2]])){
              model_selection <- list(z,logZ,g,K,edge_rem)
            }
          }
          if(deforest_criterion=="SoC" & all(table(z)>=options$min_size)){
            if(sum(logZ)>sum(model_selection[[2]])){
              model_selection <- list(z,logZ,g,K,edge_rem)
            }
          }
          if(deforest_criterion=="Diverse"){
            undiv_clu <- sapply(1:K,FUN = function(u,res,clu_al,ups){all(table(factor(res[which(clu_al==u)],levels=0:1))>=ups)},res=y,clu_al=z,ups=options$n_min_class)
            if(all(undiv_clu) & sum(logZ)>sum(model_selection[[2]])){
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
    if(K==1){
      if(options$diagnostics){
        pnoc <- list(Cluster_Allocation = z,
                    Log_Marginal_Likelihoods = logZ,
                    Graph = g,
                    Number_of_Clusters = K,
                    Edges_Removed = edge_rem,
                    Diagnostics = Track[,1:3])
      } else{
        pnoc <- list(Cluster_Allocation = z,
                    Log_Marginal_Likelihoods = logZ,
                    Graph = g,
                    Number_of_Clusters = K,
                    Edges_Removed = edge_rem)
      }
    } else{
      if(options$diagnostics){
        Tr_in <- Track[,1:3]
      } else{
        Tr_in <- NULL
      }
      pnoc <- deforest.noc(obs = X,res = y,gra = g,lbe = logZ,eps = edge_rem,
                           K_dag = options$max_K,clu_al = z,c_s = combine_save,
                           est_thres = options$SMC_thres,
                           mtb = options$BIC_memo_thres,
                           mts = options$SMC_memo_thres,par_no = options$N,
                           rfun = rprior,pdf_fun = dprior,efsamp = options$ess,
                           methas = options$n_move,vb = verbose,
                           cb = options$BIC_cache,cs = options$SMC_cache,
                           PA = MoreArgs,diagnostics = options$diagnostics,
                           Tr = Tr_in,SMC_f = IBIS.Z,BIC_f = memo.bic,
                           rt = options$RIBIS_thres)
    }
    memoise::forget(memo.bic)
    memoise::forget(IBIS.Z)
    if(sum(model_selection[[2]])>sum(pnoc[[2]])){
      if(options$diagnostics){
        Track <- pnoc$Diagnostics
        Track <- rbind(Track,c("Rev.",sum(model_selection[[2]]),model_selection[[4]]))
        res <- list(Cluster_Allocation = model_selection[[1]],
                    Log_Marginal_Likelihoods = model_selection[[2]],
                    Graph = model_selection[[3]],
                    Number_of_Clusters = model_selection[[4]],
                    Edges_Removed = model_selection[[5]],
                    Diagnostics = Track)
      } else{
        res <- list(Cluster_Allocation = model_selection[[1]],
                    Log_Marginal_Likelihoods = model_selection[[2]],
                    Graph = model_selection[[3]],
                    Number_of_Clusters = model_selection[[4]],
                    Edges_Removed = model_selection[[5]])
      }
    } else{
      res <- pnoc
    }
  }
  if(deforest_criterion=="SoC"){
    if(K==1){
      if(options$diagnostics){
        psoc <- list(Cluster_Allocation = z,
                    Log_Marginal_Likelihoods = logZ,
                    Graph = g,
                    Number_of_Clusters = K,
                    Edges_Removed = edge_rem,
                    Diagnostics = Track[,c(1:2,4:5)])
      } else{
        psoc <- list(Cluster_Allocation = z,
                    Log_Marginal_Likelihoods = logZ,
                    Graph = g,
                    Number_of_Clusters = K,
                    Edges_Removed = edge_rem)
      }
    } else{
      if(options$diagnostics){
        Tr_in <- Track[,c(1:2,4:5)]
      } else{
        Tr_in <- NULL
      }
      psoc <- deforest.soc(obs = X,res = y,gra = g,lbe = logZ,eps = edge_rem,
                           n_dag = options$min_size,clu_al = z,c_s = combine_save,
                           est_thres = options$SMC_thres,
                           mtb = options$BIC_memo_thres,
                           mts = options$SMC_memo_thres,par_no = options$N,
                           rfun = rprior,pdf_fun = dprior,efsamp = options$ess,
                           methas = options$n_move,vb = verbose,
                           cb = options$BIC_cache,cs = options$SMC_cache,
                           PA = MoreArgs,diagnostics = options$diagnostics,
                           Tr = Tr_in,SMC_f = IBIS.Z,
                           BIC_f = memo.bic, rt = options$RIBIS_thres)
    }
    memoise::forget(memo.bic)
    memoise::forget(IBIS.Z)
    if(sum(model_selection[[2]])>sum(psoc[[2]])){
      if(options$diagnostics){
        Track <- psoc$Diagnostics
        tab_z <- table(model_selection[[1]])
        Track <- rbind(Track,c("Rev.",sum(model_selection[[2]]),
                               min(tab_z),
                               sum(tab_z<options$min_size)))
        res <- list(Cluster_Allocation = model_selection[[1]],
                    Log_Marginal_Likelihoods = model_selection[[2]],
                    Graph = model_selection[[3]],
                    Number_of_Clusters = model_selection[[4]],
                    Edges_Removed = model_selection[[5]],
                    Diagnostics = Track)
      } else{
        res <- list(Cluster_Allocation = model_selection[[1]],
                    Log_Marginal_Likelihoods = model_selection[[2]],
                    Graph = model_selection[[3]],
                    Number_of_Clusters = model_selection[[4]],
                    Edges_Removed = model_selection[[5]])
      }
    } else{
      res <- psoc
    }
  }
  if(deforest_criterion=="Diverse"){
    if(K==1){
      if(options$diagnostics){
        pdiv <- list(Cluster_Allocation = z,
                    Log_Marginal_Likelihoods = logZ,
                    Graph = g,
                    Number_of_Clusters = K,
                    Edges_Removed = edge_rem,
                    Diagnostics = Track[,c(1:2,6:7)])
      } else{
        pdiv <- list(Cluster_Allocation = z,
                    Log_Marginal_Likelihoods = logZ,
                    Graph = g,
                    Number_of_Clusters = K,
                    Edges_Removed = edge_rem)
      }
    } else{
      if(options$diagnostics){
        Tr_in <- Track[,c(1:2,6:7)]
      } else{
        Tr_in <- NULL
      }
      pdiv <- deforest.diverse(obs = X,res = y,gra = g,lbe = logZ,eps = edge_rem,
                                ups = options$n_min_class,clu_al = z,
                                c_s = combine_save,est_thres = options$SMC_thres,
                                mtb = options$BIC_memo_thres,
                                mts = options$SMC_memo_thres,par_no = options$N,
                                rfun = rprior,pdf_fun = dprior,
                                efsamp = options$ess,methas = options$n_move,
                                vb = verbose,cb = options$BIC_cache,
                                cs = options$SMC_cache,PA = MoreArgs,
                                diagnostics = options$diagnostics,
                                Tr = Tr_in,SMC_f = IBIS.Z,
                                BIC_f = memo.bic, rt = options$RIBIS_thres)
    }
    memoise::forget(memo.bic)
    memoise::forget(IBIS.Z)
    if(sum(model_selection[[2]])>sum(pdiv[[2]])){
      if(options$diagnostics){
        Track <- pdiv$Diagnostics
        sap_div <- sapply(1:(model_selection[[4]]),FUN = function(u,res,clu_al){min(table(factor(res[which(clu_al==u)],levels=0:1)))},res=y,clu_al=model_selection[[1]])
        Track <- rbind(Track,c("Rev.",sum(model_selection[[2]]),min(sap_div),sum(sap_div<options$n_min_class)))
        res <- list(Cluster_Allocation = model_selection[[1]],
                    Log_Marginal_Likelihoods = model_selection[[2]],
                    Graph = model_selection[[3]],
                    Number_of_Clusters = model_selection[[4]],
                    Edges_Removed = model_selection[[5]],
                    Diagnostics = Track)
      } else{
        res <- list(Cluster_Allocation = model_selection[[1]],
                    Log_Marginal_Likelihoods = model_selection[[2]],
                    Graph = model_selection[[3]],
                    Number_of_Clusters = model_selection[[4]],
                    Edges_Removed = model_selection[[5]])
      }
    } else{
      res <- pdiv
    }
  }
  if(deforest_criterion=="None"){
    memoise::forget(memo.bic)
    memoise::forget(IBIS.Z)
    if(options$diagnostics){
      res <- list(Cluster_Allocation = z,
                  Log_Marginal_Likelihoods = logZ,
                  Graph = g,
                  Number_of_Clusters = K,
                  Edges_Removed = edge_rem,
                  Diagnostics = Track[,1:2])
    } else{
      res <- list(Cluster_Allocation = z,
                  Log_Marginal_Likelihoods = logZ,
                  Graph = g,
                  Number_of_Clusters = K,
                  Edges_Removed = edge_rem)
    }
  }
  if(deforest_criterion=="MaxReg"){
    if(K==1){
      if(options$diagnostics){
        pmaxreg <- list(Cluster_Allocation = z,
                    Log_Marginal_Likelihoods = logZ,
                    Graph = g,
                    Number_of_Clusters = K,
                    Edges_Removed = edge_rem,
                    Diagnostics = Track[,1:2])
      } else{
        pmaxreg <- list(Cluster_Allocation = z,
                    Log_Marginal_Likelihoods = logZ,
                    Graph = g,
                    Number_of_Clusters = K,
                    Edges_Removed = edge_rem)
      }
    } else{
      if(options$diagnostics){
        Tr_in <- Track[,1:2]
      } else{
        Tr_in <- NULL
      }
      pmaxreg <- deforest.maxreg(obs = X,res = y,gra = g,lbe = logZ,
                                 eps = edge_rem, tau = options$reg,clu_al = z,
                                 c_s = combine_save,
                                 est_thres = options$SMC_thres,
                                 mtb = options$BIC_memo_thres,
                                 mts = options$SMC_memo_thres,par_no = options$N,
                                 rfun = rprior,pdf_fun = dprior,
                                 efsamp = options$ess,methas = options$n_move,
                                 vb = verbose,cb = options$BIC_cache,
                                 cs = options$SMC_cache,PA = MoreArgs,
                                 diagnostics = options$diagnostics,
                                 Tr = Tr_in,SMC_f = IBIS.Z,BIC_f = memo.bic,
                                 rt = options$RIBIS_thres)
    }
    memoise::forget(memo.bic)
    memoise::forget(IBIS.Z)
    res <- pmaxreg
  }
  if(deforest_criterion=="Validation"){
    if(K==1){
      if(options$diagnostics){
        Tr <- data.frame(Track[,1:2],Log_Bayesian_Evidence_All = rep(NA,nrow(Track[,1:2])),
                           Robustness_Statistic = rep(NA,nrow(Track[,1:2])))
      }
      z_all <- igraph::components(g_all)$membership
      for(k in 1:K){
        z_all[which(z_all==z_all[as.numeric(igraph::V(g)$name[which(z==k)[1]])])] <- K+k
      }
      z_all <- z_all - K
      logZ_all <- logZ
      for(k in 1:K){
        logZ_all[k] <- lbe.gen(thres = options$SMC_thres,obs_mat = X_all,
                               res_vec = y_all,obs_ind = which(z_all==k),
                               memo_thres_bic = options$BIC_memo_thres,
                               memo_thres_smc = options$SMC_memo_thres,
                               p_num = options$N,rpri = rprior,p_pdf = dprior,
                               efs = options$ess,nm = options$n_move,
                               cache_bic = options$BIC_cache,
                               cache_smc = options$SMC_cache,MA = MoreArgs,
                               SMC_fun = IBIS.Z,BIC_fun = memo.bic,
                               ribis_thres = options$RIBIS_thres)
      }
      RobS <- sum(logZ_all-logZ)
      if(options$diagnostics){
        Tr[nrow(Tr),3:4] <- c(sum(logZ_all),RobS)
      }
      TD = list(Cluster_Allocation = z,
                Log_Marginal_Likelihoods = logZ,
                Graph = g,
                Number_of_Clusters = K,
                Edges_Removed = edge_rem)
      if(options$diagnostics){
        AD = list(Cluster_Allocation = z_all,
                  Log_Marginal_Likelihoods = logZ_all,
                  Graph = g_all,
                  Number_of_Clusters = K,
                  Edges_Removed = edge_rem,
                  Diagnostics = Tr)
      } else{
        AD = list(Cluster_Allocation = z_all,
                  Log_Marginal_Likelihoods = logZ_all,
                  Graph = g_all,
                  Number_of_Clusters = K,
                  Edges_Removed = edge_rem)
      }
      pval <- list(Training_Data = TD,
                   All_Data = AD)
    } else{
      if(options$diagnostics){
        Tr_in <- Track[,1:2]
      } else{
        Tr_in <- NULL
      }
      pval <- deforest.validation(obs = X,obs_all = X_all,res = y,res_all = y_all,
                                  gra = g,lbe = logZ,eps = edge_rem,gra_all = g_all,
                                  clu_al = z,c_s = combine_save,
                                  est_thres = options$SMC_thres,
                                  mtb = options$BIC_memo_thres,
                                  mts = options$SMC_memo_thres,par_no = options$N,
                                  rfun = rprior,pdf_fun = dprior,
                                  efsamp = options$ess,methas = options$n_move,
                                  rho = mst_var,vb = verbose,
                                  cb = options$BIC_cache,cs = options$SMC_cache,
                                  PA = MoreArgs,diagnostics = options$diagnostics,
                                  Tr = Tr_in,SMC_f = IBIS.Z,
                                  BIC_f = memo.bic,rt = options$RIBIS_thres)
    }
    memoise::forget(memo.bic)
    memoise::forget(IBIS.Z)
    res <- pval
    X <- X_all
    y <- y_all
  }
  X <- data.frame(X)
  if(!is.null(X_names)){
    colnames(X) <- X_names
  }
  res <- append(list(Covariate_Matrix = X,Response_Vector = y,
                     Minimum_Spanning_Tree_Variables = mst_var,
                     Control = options,
                     Deforestation_Criterion = deforest_criterion,
                     Prior_Mean = prior_mean,Prior_Variance = prior_var),
                list(Model = res))
  class(res) <- 'UNCOVER'
  res
}

##' Print UNCOVER
##'
##'
##' @export
##' @keywords print UNCOVER
##' @name print.UNCOVER
##' @description Prints summary information from an UNCOVER object.
##'
##' @param x Object of class `"UNCOVER"`
##' @param ... Further arguments passed to or from other methods
##' @return No return value, called for side effects
##' @details When running the function [UNCOVER()] the printed
##' information will contain information regarding; the number of clusters, the
##' cluster sizes, the sub-model log Bayesian evidences and the total model log
##' Bayesian evidence.
##' @seealso [UNCOVER()]
##'

print.UNCOVER <- function(x,...){
  if(x$Deforestation_Criterion=="Validation"){
    x$Model <- x$Model$All_Data
  }
  cat(x$Model$Number_of_Clusters,"clusters were uncovered \n")
  cat("\n")
  cat("Cluster sizes:\n")
  cat("\n")
  xz <- as.vector(table(x$Model$Cluster_Allocation))
  names(xz) <- paste0("Cluster ",1:x$Model$Number_of_Clusters)
  print.default(xz,print.gap = 2)
  cat("\n")
  cat("Sub-model log Bayesian evidences:\n")
  cat("\n")
  xx <- x$Model$Log_Marginal_Likelihoods
  names(xx) <- paste0("Cluster ",1:x$Model$Number_of_Clusters)
  print.default(xx)
  cat("\n")
  cat("Total model log Bayesian evidence:",sum(x$Model$Log_Marginal_Likelihoods))
}

##' Prediction method for UNCOVER
##'
##'
##' @export
##' @name predict.UNCOVER
##' @description Predicts the response of new observations and their cluster
##' assignment using an object of class `"UNCOVER"`.
##'
##' @keywords predict UNCOVER
##' @param object Object of class `"UNCOVER"`
##' @param newX Data frame containing new observations to predict. If not
##' specified the fitted values will be returned instead.
##' @param type Either `"prob"` for a probabilistic response prediction or
##' `"response"` for a hard output of the predicted response
##' @param ... Additional arguments affecting the predictions produced
##' @return Either a data frame of response probabilities with cluster
##' assignment for each observation or a data frame of predicted responses with
##' cluster assignment for each observation.
##' @details Note that this is a Bayesian prediction method and so samples of
##' the posterior, defined by `"UNCOVER"` object provided, will be obtained
##' through SMC methods for prediction. See [IBIS.logreg()] for
##' more details.
##' @seealso [UNCOVER()], [IBIS.logreg()]
##' @examples
##'
##' \donttest{
##' # First we generate a co-variate matrix and binary response vector
##' CM <- data.frame(X1 = rnorm(100),X2 = rnorm(100))
##' rv <- sample(0:1,100,replace=TRUE)
##'
##' # We can then run UNCOVER with no deforestation criteria
##' UN.none <- UNCOVER(X = CM,y = rv, deforest_criterion = "None", verbose = FALSE)
##'
##' # The fitted values of UN.none can then be obtained
##' predict(UN.none)
##' predict(UN.none,type = "response")
##'
##' # We can also predict the response for new data
##' CM.2 <- data.frame(X1 = rnorm(10),X2 = rnorm(10))
##' cbind(CM.2,predict(UN.none,newX = CM.2))
##' }
##'

predict.UNCOVER <- function(object,newX=NULL,type="prob",...){
  if(is.null(newX)){
    newX = object$Covariate_Matrix
  }
  newX <- as.matrix(newX)
  if(type!="prob" & type!="response"){
    stop("type not supported")
  }
  object.infer <- UNCOVER.infer(object)
  newX.assign <- UNCOVER.assign(object,newX)
  newX <- as.matrix(newX,ncol = ncol(object$Covariate_Matrix))
  DM <- cbind(rep(1,nrow(newX)),newX)
  p1 <- rep(0,nrow(newX))
  for(i in 1:nrow(newX)){
    p1[i] <- rowMeans(1/(1+exp(-DM[i,,drop=FALSE]%*%t(object.infer[[newX.assign[i]]]))))
  }
  if(type=="prob"){
    res <- data.frame(1-p1,p1,newX.assign)
    colnames(res) <- c(0:1,"Cluster")
  } else{
    res <- rep(0,length(p1))
    res[which(p1>=0.5)] <- 1
    res <- data.frame(newX.assign,res)
    colnames(res) <- c("Cluster","Predicted Response")
  }
  return(res)
}

##' Plot various outputs of UNCOVER
##'
##'
##' @export
##' @name plot.UNCOVER
##' @description Allows visualisation of many aspects of UNCOVER, including
##' co-variate, posterior and diagnostic plots.
##'
##' @keywords plot UNCOVER
##' @param x Object of class `"UNCOVER"`
##' @param type Can be one of; `"covariates"` for cluster assignment
##' visualisation for the co-variates, `"fitted"` for co-variate visualisation
##' with respect to their fitted values, `"samples"` for posterior visualisation
##' or `"diagnostics"` for diagnostic plots. See details. Defaults to
##' `"covariates"`.
##' @param plot_var Vector specifying which columns (or associated logistic
##' regression coefficients) of the co-variate matrix should be plotted. Does not
##' apply when `type=="diagnostics"`. Defaults to all columns being selected.
##' @param diagnostic_x_axis Only applies if `"type=="diagnostics"`. Either
##' `"full"` (default) for all observations indices to be plotted on the x-axis
##' or `"minimal"` for only some observations indices to be plotted on the
##' x-axis.
##' @param ... Arguments to be passed to methods
##' @return No return value, called for side effects
##' @details If `type=="covariates"`, the resulting plot will be a ggpairs plot.
##' The diagonal entries will be a plot of K density plots (K being the number
##' of clusters). The off-diagonal elements are scatter-plots of the
##' observations, given a label according to their true response and a colour
##' based on their assigned cluster. If `length(plot_var)==1` then the density
##' plot and the scatter-plot are combined. If a cluster contains less than two
##' data points the density will not be plotted.
##'
##' If `"type=="fitted"`, the resulting plot will be a ggpairs plot. The
##' diagonal entries will be two density plots, one for data predicted to have
##' response 0 by the model (red) and one for training data predicted to have
##' response 1 by the model (green). The off-diagonal elements are
##' scatter-plots of the observations, given a label according to their actual
##' response and a colour scale based on their predicted response. If
##' `length(plot_var)==1` then the density plot and the scatter-plot are
##' combined. If a predicted class (0 or 1) contains less than two data points
##' the density will not be plotted.
##'
##' If `type=="samples"`, the resulting plot will be a ggpairs plot of the
##' clusters posteriors, giving the coefficient densities in the diagonal and
##' scatter-plots of the posterior samples in the off-diagonal. The transparency
##' is increased in the upper triangle for scenarios when posteriors overlap.
##'
##' If `"type==diagnostics"`, the resulting plot depends on the deforestation
##' criterion used to create the `"UNCOVER"` object:
##' \describe{
##' \item{`"None"`}{A plot tracking the overall log Bayesian evidence every time
##' an action is executed.}
##' \item{`"NoC"`}{A plot tracking the overall log Bayesian evidence after every
##' action and a plot tracking the number of clusters after every action.}
##' \item{`"SoC"`}{Three plots; one tracking the overall log Bayesian evidence
##' after every action, one tracking the number of criterion breaking clusters
##' after every action and one tracking the minimum cluster size after every
##' action.}
##' \item{`"MaxReg"`}{A plot tracking the overall log Bayesian evidence every
##' time an action is executed. Actions are coloured and each action has an
##' associated coloured dashed line indicating the log Bayesian evidence plus
##' the logarithm of the maximum tolerance provided.}
##' \item{`"Validation"`}{A plot tracking the overall log Bayesian evidence
##' after every action (for both the training data and all of the data) and a
##' plot tracking the robustness statistic after every deforestation action.}
##' \item{`"Diverse"`}{Three plots; one tracking the overall log Bayesian
##' evidence after every action, one tracking the number of criterion breaking
##' clusters after every action and one tracking the minimum minority class
##' across clusters after every action.}
##' }
##' Actions are defined as either edge removals, edge additions or edge
##' additions in the deforestation stage. The syntax for an action will be the
##' 'type_of_action.edge'. For example the removal of an edge connecting
##' observation 1 and observation 2 will be displayed 'Rem.1-2'. If the edge was
##' being added this would be displayed 'Def.Add.1-2' if in the deforestation
##' stage and 'Add.1-2' otherwise. When the data for the `"UNCOVER"` object
##' created is large setting `diagnostic_x_axis=="minimal"` is recommended as it
##' gives a more visually appealing output.
##'
##' @seealso [UNCOVER()]
##' @examples
##'
##' \donttest{
##' require(graphics)
##' # First we generate a co-variate matrix and binary response vector
##' CM <- matrix(rnorm(200),100,2)
##' rv <- sample(0:1,100,replace=TRUE)
##'
##' # We can then run our algorithm for each of the different deforestation
##' # criteria
##' UN.none <- UNCOVER(X = CM,y = rv, deforest_criterion = "None", verbose = FALSE)
##' UN.noc <- UNCOVER(X = CM,y = rv, deforest_criterion = "NoC",
##'                   options = UNCOVER.opts(max_K = 3), verbose = FALSE)
##' UN.soc <- UNCOVER(X = CM,y = rv, deforest_criterion = "SoC",
##'                   options = UNCOVER.opts(min_size = 10), verbose = FALSE)
##' UN.maxreg <- UNCOVER(X = CM,y = rv, deforest_criterion = "MaxReg",
##'                      options = UNCOVER.opts(reg = 1), verbose = FALSE)
##' UN.validation <- UNCOVER(X = CM,y = rv, deforest_criterion = "Validation",
##'                          options = UNCOVER.opts(train_frac = 0.8),
##'                          verbose = FALSE)
##' UN.diverse <- UNCOVER(X = CM,y = rv, deforest_criterion = "Diverse",
##'                        options = UNCOVER.opts(n_min_class = 2), verbose = FALSE)
##' plot(UN.none,type = "covariates")
##' plot(UN.none,type = "fitted")
##' plot(UN.none,type = "samples")
##' plot(UN.none,type = "diagnostics",diagnostic_x_axis = "minimal")
##' plot(UN.noc,type = "diagnostics",diagnostic_x_axis = "minimal")
##' plot(UN.soc,type = "diagnostics",diagnostic_x_axis = "minimal")
##' plot(UN.maxreg,type = "diagnostics",diagnostic_x_axis = "minimal")
##' plot(UN.validation,type = "diagnostics",diagnostic_x_axis = "minimal")
##' plot(UN.diverse,type = "diagnostics",diagnostic_x_axis = "minimal")
##'
##' # If we only wanted to view the second co-variate
##' plot(UN.none,type = "covariates",plot_var=2)
##' plot(UN.none,type = "fitted",plot_var=2)
##' plot(UN.none,type = "samples",plot_var=2)
##' }
##'


plot.UNCOVER <- function(x,type = "covariates",
                         plot_var = x$Minimum_Spanning_Tree_Variables,
                         diagnostic_x_axis = "full",...){
  if(diagnostic_x_axis!="full" & diagnostic_x_axis!="minimal"){
    stop("diagnostic_x_axis must be either full or minimal")
  }
  if(any(is.na(match(plot_var,1:ncol(x$Covariate_Matrix))))){
    stop("cannot subset the co-variate matrix with plot_var provided")
  }
  if(type!="covariates" & type!="fitted" & type!="samples" & type!="diagnostics"){
    stop("type not supported")
  }
  if(type=="covariates"){
    if(x$Deforestation_Criterion=="Validation"){
      x$Model <- x$Model$All_Data
    }
    x$Covariate_Matrix <- x$Covariate_Matrix[,plot_var,drop=FALSE]
    if(ncol(x$Covariate_Matrix)==1){
      plotdf <- data.frame(x$Covariate_Matrix,y.axis = rep(0,nrow(x$Covariate_Matrix)),
                           y = as.character(x$Response_Vector),
                           Cluster = as.factor(x$Model$Cluster_Allocation))
      overall_plot <- ggplot2::ggplot(plotdf) +
        ggplot2::geom_density(show.legend = FALSE,
                              ggplot2::aes_string(colnames(plotdf)[1],
                                                  colour = "Cluster",
                                                  fill = "Cluster"),alpha = 0.25) +
        ggplot2::geom_point(alpha=0,
                            ggplot2::aes_string(colnames(plotdf)[1],
                                                "y.axis",colour = "Cluster")) +
        ggplot2::geom_text(alpha = 1,show.legend = FALSE,size = 3,
                           ggplot2::aes_string(colnames(plotdf)[1],
                                               "y.axis",
                                               label = "y",colour = "Cluster")) +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha=1),
                                                       title = "Cluster")) +
        ggplot2::labs(y = "density") +
        ggplot2::theme(legend.position = "bottom")
    } else{
      plotdf <- data.frame(x$Covariate_Matrix[,c(1,2)],
                           Cluster = as.factor(x$Model$Cluster_Allocation))
      legend_plot <- ggplot2::ggplot(plotdf,
                                     ggplot2::aes_string(colnames(plotdf)[1],
                                                         colnames(plotdf)[2],
                                                         colour = "Cluster")) +
        ggplot2::geom_point() + ggplot2::theme(legend.position = "bottom")
      label_plot <- function(data,mapping,...){
        ggplot2::ggplot(data = data,mapping = mapping) +
        ggplot2::geom_point(alpha=0,...) +
          ggplot2::geom_text(alpha = 1,size = 3,...)
      }
      overall_plot <- GGally::ggpairs(data.frame(x$Covariate_Matrix,
                                                 y = as.character(x$Response_Vector),
                                                 Cluster = as.factor(x$Model$Cluster_Allocation)),
                                      ggplot2::aes_string(colour = "Cluster",
                                                          alpha = 0.5,
                                                          label = "y",
                                                          fill = "Cluster"),
                                      columns = 1:ncol(x$Covariate_Matrix),
                                      upper = list(continuous = label_plot),
                                      lower = list(continuous = label_plot),
                                      legend = GGally::grab_legend(legend_plot)) +
        ggplot2::theme(legend.position = "bottom")
    }
  }
  if(type=="fitted"){
    X_prob <- predict.UNCOVER(object = x)[,2]
    x$Covariate_Matrix <- x$Covariate_Matrix[,plot_var,drop=FALSE]
    if(ncol(x$Covariate_Matrix)==1){
      plotdf <- data.frame(x$Covariate_Matrix,
                           y.axis = rep(0,nrow(x$Covariate_Matrix)),
                           y = as.character(x$Response_Vector),
                           Fitted.Probabilities = X_prob,
                           Hard.Assignment = X_prob>=0.5)
      overall_plot <- ggplot2::ggplot(plotdf) +
        ggplot2::geom_density(show.legend = FALSE,
                              ggplot2::aes_string(colnames(plotdf)[1],
                                                  colour = "Hard.Assignment",
                                                  fill = "Hard.Assignment"),alpha=0.25) +
        ggplot2::scale_colour_manual(values = c("red","green")) +
        ggplot2::scale_fill_manual(values = c("red","green")) +
        ggnewscale::new_scale_color() +
        ggplot2::geom_point(alpha=0,
                            ggplot2::aes_string(colnames(plotdf)[1],
                                                "y.axis",
                                                colour = "Fitted.Probabilities")) +
        ggplot2::geom_text(alpha = 1,show.legend = FALSE,size = 3,
                           ggplot2::aes_string(colnames(plotdf)[1],
                                               "y.axis",
                                               label = "y",
                                               colour = "Fitted.Probabilities")) +
        ggplot2::scale_color_gradient(low = "red",high = "green",name = "Fitted Probabilities",limits = c(0,1)) +
        ggplot2::labs(y = "density") +
        ggplot2::theme(legend.position = "bottom")
    } else{
      plotdf <- data.frame(x$Covariate_Matrix[,c(1,2)],
                           Fitted.Probabilities = X_prob)
      legend_plot <- ggplot2::ggplot(plotdf,
                                     ggplot2::aes_string(colnames(plotdf)[1],
                                         colnames(plotdf)[2],
                                         colour = "Fitted.Probabilities")) +
        ggplot2::geom_point() + ggplot2::theme(legend.position = "bottom") +
        ggplot2::scale_color_gradient(low = "red",high = "green",name = "Fitted Probabilities",limits = c(0,1))
      diag_plot <- function(data,mapping,...){
        ggplot2::ggplot(data = data,mapping = mapping) +
          ggplot2::geom_density(ggplot2::aes_string(colour = "Hard.Assignment",
                                                    fill = "Hard.Assignment")) +
          ggplot2::scale_colour_manual(values = c("red","green")) +
          ggplot2::scale_fill_manual(values = c("red","green"))
      }
      label_plot_fitted <- function(data,mapping,...){
        ggplot2::ggplot(data = data,mapping = mapping) +
          ggplot2::geom_point(alpha=0,
                              ggplot2::aes_string(colour = "Fitted.Probabilities")) +
          ggplot2::geom_text(alpha = 1,size = 3,
                             ggplot2::aes_string(colour = "Fitted.Probabilities",
                                                 label = "y")) +
          ggplot2::scale_color_gradient(low = "red",high = "green",
                                        name = "Fitted Probabilities",
                                        limits = c(0,1))
      }
      overall_plot <- GGally::ggpairs(data.frame(x$Covariate_Matrix,
                                                 y = as.character(x$Response_Vector),
                                                 Fitted.Probabilities = X_prob,
                                                 Hard.Assignment = X_prob>=0.5),
                                      ggplot2::aes_string(alpha = 0.5),
                                      columns = 1:ncol(x$Covariate_Matrix),
                                      diag = list(continuous = diag_plot),
                                      upper = list(continuous = label_plot_fitted),
                                      lower = list(continuous = label_plot_fitted),
                                      legend = GGally::grab_legend(legend_plot)) +
        ggplot2::theme(legend.position = "bottom")
    }
  }
  if(type=="samples"){
    samp.df <- data.frame(do.call(rbind,UNCOVER.infer(x))[,c(1,plot_var+1)])
    colnames(samp.df) <- paste0("beta[",c(0,plot_var),"]")
    if(x$Deforestation_Criterion=="Validation"){
      x$Model <- x$Model$All_Data
    }
    plotdf <- data.frame(x$Covariate_Matrix[,c(1,2)],
                         Cluster = as.factor(x$Model$Cluster_Allocation))
    legend_plot <- ggplot2::ggplot(plotdf,
                                   ggplot2::aes_string(colnames(plotdf)[1],
                                                       colnames(plotdf)[2],
                                                       colour = "Cluster")) +
      ggplot2::geom_point() + ggplot2::theme(legend.position = "bottom")
    overall_plot <- GGally::ggpairs(samp.df,ggplot2::aes_string(alpha = 0.5,
                                                         color = as.factor(rep(1:x$Model$Number_of_Clusters,each = x$Control$N)),
                                                         fill = as.factor(rep(1:x$Model$Number_of_Clusters,each = x$Control$N))),
                                    labeller = "label_parsed",
                                    upper = list(continuous = GGally::wrap("points",size=0.5,
                                                                           alpha = 0.25)),
                                    lower = list(continuous = GGally::wrap("points",size=0.5)),
                                    legend = GGally::grab_legend(legend_plot)) +
      ggplot2::theme(legend.position = "bottom")
  }
  if(type=="diagnostics"){
    if(x$Deforestation_Criterion=="Validation"){
      modelnames <- names(x$Model$All_Data)
    } else{
      modelnames <- names(x$Model)
    }
    if(!("Diagnostics" %in% modelnames)){
      stop("For diagnostics plot diagnostics=TRUE needs to be specified when creating UNCOVER object")
    }
    if(x$Deforestation_Criterion=="None"){
      obs_lev <- x$Model$Diagnostics$Action
      obs_lev_tics <- rep("",length(obs_lev))
      obs_lev_tics[c(1,round(length(obs_lev)*(1:4)/4))] <- obs_lev[c(1,round(length(obs_lev)*(1:4)/4))]
      x$Model$Diagnostics$Action <- factor(x$Model$Diagnostics$Action,levels = obs_lev)
      x$Model$Diagnostics$Log_Bayesian_Evidence <- as.numeric(x$Model$Diagnostics$Log_Bayesian_Evidence)
      overall_plot <- ggplot2::ggplot(x$Model$Diagnostics,
                                      ggplot2::aes_string(x = "Action",
                                                   y = "Log_Bayesian_Evidence",group=1)) +
        ggplot2::geom_point() + ggplot2::geom_line() +
        ggplot2::labs(x = "Action",y = "Log Bayesian Evidence")
      if(diagnostic_x_axis=="minimal"){
        overall_plot <- overall_plot + ggplot2::scale_x_discrete(labels = obs_lev_tics)
      }
    }
    if(x$Deforestation_Criterion=="NoC"){
      obs_lev <- x$Model$Diagnostics$Action
      obs_lev_tics <- rep("",length(obs_lev))
      obs_lev_tics[c(1,round(length(obs_lev)*(1:4)/4))] <- obs_lev[c(1,round(length(obs_lev)*(1:4)/4))]
      x$Model$Diagnostics$Action <- factor(x$Model$Diagnostics$Action,levels = obs_lev)
      x$Model$Diagnostics$Log_Bayesian_Evidence <- as.numeric(x$Model$Diagnostics$Log_Bayesian_Evidence)
      plot_1 <- ggplot2::ggplot(x$Model$Diagnostics,
                                      ggplot2::aes_string(x = "Action",
                                                   y = "Log_Bayesian_Evidence",group=1)) +
        ggplot2::geom_point() + ggplot2::geom_line()
      if(diagnostic_x_axis=="minimal"){
        plot_1 <- plot_1 + ggplot2::scale_x_discrete(labels = obs_lev_tics)
      }
      x$Model$Diagnostics$Number_Of_Clusters <- as.integer(x$Model$Diagnostics$Number_Of_Clusters)
      plot_2 <- ggplot2::ggplot(x$Model$Diagnostics,
                                ggplot2::aes_string(x = "Action",
                                             y = "Number_Of_Clusters",group=1)) +
        ggplot2::geom_point() + ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = x$Control$max_K,linetype = 2) +
        ggplot2::scale_y_continuous(breaks = unique(sort(c(round(seq(min(x$Model$Diagnostics$Number_Of_Clusters),max(x$Model$Diagnostics$Number_Of_Clusters),length.out = 5)),x$Control$max_K))))
      if(diagnostic_x_axis=="minimal"){
        plot_2 <- plot_2 + ggplot2::scale_x_discrete(labels = obs_lev_tics)
      }
      overall_plot <- GGally::ggmatrix(list(plot_1,plot_2),nrow = 2,ncol = 1,
                                       xlab = "Action",
                                       yAxisLabels = c("Log Bayesian Evidence",
                                                       "Number Of Clusters"))
    }
    if(x$Deforestation_Criterion=="SoC"){
      obs_lev <- x$Model$Diagnostics$Action
      obs_lev_tics <- rep("",length(obs_lev))
      obs_lev_tics[c(1,round(length(obs_lev)*(1:4)/4))] <- obs_lev[c(1,round(length(obs_lev)*(1:4)/4))]
      x$Model$Diagnostics$Action <- factor(x$Model$Diagnostics$Action,levels = obs_lev)
      x$Model$Diagnostics$Log_Bayesian_Evidence <- as.numeric(x$Model$Diagnostics$Log_Bayesian_Evidence)
      plot_1 <- ggplot2::ggplot(x$Model$Diagnostics,
                                ggplot2::aes_string(x = "Action",
                                             y = "Log_Bayesian_Evidence",group=1)) +
        ggplot2::geom_point() + ggplot2::geom_line() +
        ggplot2::labs(x = "Action",y = "Log Bayesian Evidence")
      if(diagnostic_x_axis=="minimal"){
        plot_1 <- plot_1 + ggplot2::scale_x_discrete(labels = obs_lev_tics)
      }
      x$Model$Diagnostics$Minimum_Cluster_Size <- as.integer(x$Model$Diagnostics$Minimum_Cluster_Size)
      plot_2 <- ggplot2::ggplot(x$Model$Diagnostics,
                                ggplot2::aes_string(x = "Action",
                                             y = "Minimum_Cluster_Size",group=1)) +
        ggplot2::geom_point() + ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = x$Control$min_size,linetype = 2) +
        ggplot2::scale_y_continuous(breaks = unique(sort(c(round(seq(min(x$Model$Diagnostics$Minimum_Cluster_Size),max(x$Model$Diagnostics$Minimum_Cluster_Size),length.out = 5)),x$Control$min_size)))) +
        ggplot2::labs(x = "Action",y = "Minimum Cluster Size")
      if(diagnostic_x_axis=="minimal"){
        plot_2 <- plot_2 + ggplot2::scale_x_discrete(labels = obs_lev_tics)
      }
      x$Model$Diagnostics$Number_Of_Small_Clusters <- as.integer(x$Model$Diagnostics$Number_Of_Small_Clusters)
      plot_3 <- ggplot2::ggplot(x$Model$Diagnostics,
                                ggplot2::aes_string(x = "Action",
                                             y = "Number_Of_Small_Clusters",group=1)) +
        ggplot2::geom_point() + ggplot2::geom_line() +
        ggplot2::labs(x = "Action",y = "Number of Criterion Breaking Clusters")
      if(diagnostic_x_axis=="minimal"){
        plot_3 <- plot_3 + ggplot2::scale_x_discrete(labels = obs_lev_tics)
      }
      overall_plot <- ggpubr::ggarrange(ggpubr::ggarrange(plot_1,plot_3,ncol=2),plot_2,nrow=2)
    }
    if(x$Deforestation_Criterion=="MaxReg"){
      obs_lev <- x$Model$Diagnostics$Action
      obs_lev_tics <- rep("",length(obs_lev))
      obs_lev_tics[c(1,round(length(obs_lev)*(1:4)/4))] <- obs_lev[c(1,round(length(obs_lev)*(1:4)/4))]
      x$Model$Diagnostics$Action <- factor(x$Model$Diagnostics$Action,levels = obs_lev)
      x$Model$Diagnostics$Log_Bayesian_Evidence <- as.numeric(x$Model$Diagnostics$Log_Bayesian_Evidence)
      overall_plot <- ggplot2::ggplot(x$Model$Diagnostics,
                                ggplot2::aes_string(x = "Action",
                                             y = "Log_Bayesian_Evidence",group=1,
                                             colour = "Action")) +
        ggplot2::geom_point(show.legend = FALSE) + ggplot2::geom_line(show.legend = FALSE) +
        ggplot2::labs(x = "Action",y = "Log Bayesian Evidence")
      for(i in 1:nrow(x$Model$Diagnostics)){
        yint <- x$Model$Diagnostics$Log_Bayesian_Evidence[i] + x$Control$reg
        ycol <- scales::hue_pal()(nrow(x$Model$Diagnostics))[i]
        overall_plot <- overall_plot + ggplot2::geom_hline(yintercept = yint, colour = ycol,linetype = 2)
      }
      if(diagnostic_x_axis=="minimal"){
        overall_plot <- overall_plot + ggplot2::scale_x_discrete(labels = obs_lev_tics)
      }
    }
    if(x$Deforestation_Criterion=="Validation"){
      obs_lev <- x$Model$All_Data$Diagnostics$Action
      obs_lev_tics <- rep("",length(obs_lev))
      obs_lev_tics[c(1,round(length(obs_lev)*(1:4)/4))] <- obs_lev[c(1,round(length(obs_lev)*(1:4)/4))]
      p1df <- x$Model$All_Data$Diagnostics
      p1df <- rbind(p1df,p1df)
      p1df[(nrow(x$Model$All_Data$Diagnostics)+1):nrow(p1df),2] <- p1df[(nrow(x$Model$All_Data$Diagnostics)+1):nrow(p1df),3]
      p1df <- p1df[,1:2]
      p1df$Data <- rep(c("Training","All"),each=nrow(x$Model$All_Data$Diagnostics))
      p1df$Action <- factor(p1df$Action,levels = obs_lev)
      p1df$Log_Bayesian_Evidence <- as.numeric(p1df$Log_Bayesian_Evidence)
      plot_1 <- ggplot2::ggplot(p1df,
                                      ggplot2::aes_string(x = "Action",
                                                   y = "Log_Bayesian_Evidence",
                                                   group="Data",color = "Data")) +
        ggplot2::geom_point() + ggplot2::geom_line(show.legend = FALSE) +
        ggplot2::theme(legend.position = "top")
      if(diagnostic_x_axis=="minimal"){
        plot_1 <- plot_1 + ggplot2::scale_x_discrete(labels = obs_lev_tics)
      }
      x$Model$All_Data$Diagnostics$Action <- factor(x$Model$All_Data$Diagnostics$Action,levels = obs_lev)
      x$Model$All_Data$Diagnostics$Robustness_Statistic <- as.numeric(x$Model$All_Data$Diagnostics$Robustness_Statistic)
      plot_2 <- ggplot2::ggplot(x$Model$All_Data$Diagnostics,
                                ggplot2::aes_string(x = "Action",
                                             y = "Robustness_Statistic",
                                             group=1)) +
        ggplot2::geom_point(show.legend = FALSE) + ggplot2::geom_line(show.legend = FALSE)
      if(diagnostic_x_axis=="minimal"){
        plot_2 <- plot_2 + ggplot2::scale_x_discrete(labels = obs_lev_tics)
      }
      overall_plot <- GGally::ggmatrix(list(plot_1,plot_2),nrow = 2,ncol = 1,
                                       xlab = "Action",
                                       yAxisLabels = c("Log Bayesian Evidence",
                                                       "Log Robustness Statistic"),
                                       legend = 1) + ggplot2::theme(legend.position = "top")
    }
    if(x$Deforestation_Criterion=="Diverse"){
      obs_lev <- x$Model$Diagnostics$Action
      obs_lev_tics <- rep("",length(obs_lev))
      obs_lev_tics[c(1,round(length(obs_lev)*(1:4)/4))] <- obs_lev[c(1,round(length(obs_lev)*(1:4)/4))]
      x$Model$Diagnostics$Action <- factor(x$Model$Diagnostics$Action,levels = obs_lev)
      x$Model$Diagnostics$Log_Bayesian_Evidence <- as.numeric(x$Model$Diagnostics$Log_Bayesian_Evidence)
      plot_1 <- ggplot2::ggplot(x$Model$Diagnostics,
                                ggplot2::aes_string(x = "Action",
                                             y = "Log_Bayesian_Evidence",group=1)) +
        ggplot2::geom_point() + ggplot2::geom_line() +
        ggplot2::labs(x = "Action",y = "Log Bayesian Evidence")
      if(diagnostic_x_axis=="minimal"){
        plot_1 <- plot_1 + ggplot2::scale_x_discrete(labels = obs_lev_tics)
      }
      x$Model$Diagnostics$Minimum_Minority_Class <- as.integer(x$Model$Diagnostics$Minimum_Minority_Class)
      plot_2 <- ggplot2::ggplot(x$Model$Diagnostics,
                                ggplot2::aes_string(x = "Action",
                                             y = "Minimum_Minority_Class",group=1)) +
        ggplot2::geom_point() + ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = x$Control$n_min_class,linetype = 2) +
        ggplot2::scale_y_continuous(breaks = unique(sort(c(round(seq(min(x$Model$Diagnostics$Minimum_Minority_Class),max(x$Model$Diagnostics$Minimum_Minority_Class),length.out = 5)),x$Control$n_min_class)))) +
        ggplot2::labs(x = "Action",y = "Minimum Minority Class")
      if(diagnostic_x_axis=="minimal"){
        plot_2 <- plot_2 + ggplot2::scale_x_discrete(labels = obs_lev_tics)
      }
      x$Model$Diagnostics$Number_Of_Undiverse_Clusters <- as.integer(x$Model$Diagnostics$Number_Of_Undiverse_Clusters)
      plot_3 <- ggplot2::ggplot(x$Model$Diagnostics,
                                ggplot2::aes_string(x = "Action",
                                             y = "Number_Of_Undiverse_Clusters",group=1)) +
        ggplot2::geom_point() + ggplot2::geom_line() +
        ggplot2::labs(x = "Action",y = "Number of Criterion Breaking Clusters")
      if(diagnostic_x_axis=="minimal"){
        plot_3 <- plot_3 + ggplot2::scale_x_discrete(labels = obs_lev_tics)
      }
      overall_plot <- ggpubr::ggarrange(ggpubr::ggarrange(plot_1,plot_3,ncol=2),plot_2,nrow=2)
    }
  }
  suppressWarnings(print(overall_plot))
}

##' Additional argument generator for [UNCOVER()]
##'
##'
##' @export
##' @name UNCOVER.opts
##' @description This function is used to specify additional arguments to
##' `UNCOVER`.
##'
##' @keywords UNCOVER options control
##' @param N Number of particles for the SMC sampler. Defaults to 1000.
##' @param train_frac What fraction of the data should be used for training.
##' Should only be directly specified if `deforest_criterion == "Validation`.
##' Defaults to `1`.
##' @param max_K The maximum number of clusters allowed in the final output.
##' Should only be directly specified if `deforest_criterion == "NoC`. Defaults
##' to `Inf`.
##' @param min_size The minimum number of observations allowed for any cluster
##' in the final model. Should only be directly specified if
##' `deforest_criterion == "SoC`. Defaults to 0.
##' @param reg Numerical natural logarithm of the tolerance parameter. Must be
##' positive. Should only be directly specified if
##' `deforest_criterion == "MaxReg`. Defaults to 0.
##' @param n_min_class Each cluster will have an associated minority class.
##' `n_min_class` specifies a minimum number of observations that should have
##' that class for each and every cluster. Should only be directly specified if
##' `deforest_criterion == "Diverse`. Defaults to 0.
##' @param SMC_thres The threshold for which the number of observations needs to
##' exceed to consider using BIC as an estimator. Defaults to 30 if not
##' specified.
##' @param BIC_memo_thres Only used when estimating the log Bayesian evidence of
##' a cluster using BIC. When the number of observations exceeds `BIC_memo_thres`
##' the function checks for similar inputs evaluated previously. See details.
##' Defaults to never checking.
##' @param SMC_memo_thres Only used when estimating the log Bayesian evidence of
##' a cluster using SMC. When the number of observations exceeds `SMC_memo_thres`
##' the function checks for similar inputs evaluated previously. See details.
##' Defaults to never checking.
##' @param ess Effective Sample Size Threshold: If the effective sample size of
##' the particles falls below this value then a resample move step is
##' triggered. Defaults to `N/2`.
##' @param n_move Number of Metropolis-Hastings steps to apply each time a
##' resample move step is triggered. Defaults to 1.
##' @param prior.override Are you overriding the default multivariate normal
##' form of the prior? Defaults to `FALSE`.
##' @param rprior Function which produces samples from your prior if the default
##' prior form is to be overridden. If using the default prior form this does
##' not need to be specified.
##' @param dprior Function which produces your specified priors density for
##' inputted samples if the default prior form is to be overridden. If using the
##' default prior form this does not need to be specified.
##' @param diagnostics Should diagnostic data be recorded and outputted?
##' Defaults to `TRUE`.
##' @param RIBIS_thres The threshold for which the number of observations needs
##' to exceed to consider ever using RIBIS as an estimator. Defaults to 30 if
##' not specified. See details.
##' @param BIC_cache The cache for the function which estimates the log Bayesian
##' evidence using BIC. Defaults to a cache with standard size and least
##' recently used eviction policy.
##' @param SMC_cache The cache for the function which estimates the log Bayesian
##' evidence using SMC. Defaults to a cache with standard size and least
##' recently used eviction policy.
##' @param ... Additional arguments required for complete specification of the
##' two prior functions given, if the default prior form is to be overridden.
##' @return A list consisting of:
##'
##' \describe{
##' \item{`N`}{Number of particles for the SMC sampler}
##' \item{`train_frac`}{Training data fraction}
##' \item{`max_K`}{Maximum number of clusters allowed}
##' \item{`min_size`}{Minimum size of clusters allowed}
##' \item{`reg`}{Log of the maximum regret tolerance parameter}
##' \item{`n_min_class`}{Minimum size of cluster minority class allowed}
##' \item{`SMC_thres`}{Threshold for when estimation with BIC is attempted}
##' \item{`BIC_memo_thres`}{Threshold for when we review previous inputs of the
##' BIC function for similarities}
##' \item{`SMC_memo_thres`}{Threshold for when we review previous inputs of the
##' SMC function for similarities}
##' \item{`ess`}{Effective Sample Size Threshold}
##' \item{`n_move`}{Number of Metropolis-Hastings steps}
##' \item{`rprior`}{Function which produces samples from your prior. `NULL` if
##' `prior.override==FALSE`.}
##' \item{`dprior`}{Function which produces your specified priors density for
##' inputted samples. `NULL` if `prior.override==FALSE`.}
##' \item{`prior.override`}{Logical value indicating if the prior has been
##' overridden or not}
##' \item{`diagnostics`}{Logical value indicating whether diagnostic information
##' should be included in the output of `UNCOVER`}
##' \item{`RIBIS_thres`}{The threshold for allowing the use of RIBIS}
##' \item{`BIC_cache`}{Cache for the memoised function which estimates the log
##' Bayesian evidence using BIC}
##' \item{`SMC_cache`}{Cache for the memoised function which estimates the log
##' Bayesian evidence using SMC}
##' \item{`MoreArgs`}{A list of the additional arguments required for `rprior`
##' and `dprior`. `NULL` if `prior.override==FALSE`.}
##' }
##'
##' @details This function should only be used to provide additional control
##' arguments to `UNCOVER`. Arguments that are for a particular deforestation
##' criteria should not be altered from the defaults for other deforestation
##' criteria.
##'
##' BIC refers to the Bayesian Information Criterion. The use of BIC when
##' estimating the log Bayesian evidence is valid assuming the number of
##' observations is large, and if specifying `SMC_thres` this should be balanced
##' with computational expense (as the function which relies
##' on BIC values is much faster than the SMC sampler).
##'
##' In an attempt to improve computational time, the SMC sampler along with the
##' function which uses BIC values are memoised, with the cache for each of
##' these memoised functions be specified by `SMC_cache` and `BIC_cache`
##' respectively. See [memoise::memoise()] for more details. If we do
##' not get and each match from the function input to a previously evaluated
##' input, we may wish to search the cache for similar inputs which could
##' provide a reasonable starting point. Checking the cache however takes time,
##' and so we allow the user to specify at which size of cluster to they deem it
##' worthwhile to check. Which value threshold to select to optimise run time is
##' problem specific, however for `BIC_memo_thres` it is almost always
##' beneficial to never check the cache (the exception for this being when the
##' cluster sizes are extremely large, for example containing a million
##' observations). `SMC_memo_thres` can be much lower as the SMC sampler is a
##' much more expensive function to run. See Emerson and Aslett (2023) for more
##' details.
##'
##' `RIBIS_thres` can be specified to have a higher value to ensure that the
##' asymptotic properties which Reverse Iterated Batch Importance Sampling
##' (RIBIS) relies upon hold. See Emerson and Aslett (2023) for more details.
##'
##' Specifying `rprior` and `dprior` will not override the default prior form
##' unless `prior.override=TRUE`. If a multivariate normal form is required then
##' the arguments for this prior should be specified in `UNCOVER`.
##' @seealso [UNCOVER()]
##' @references \itemize{
##' \item Emerson, S.R. and Aslett, L.J.M. (2023). Joint cohort and prediction
##' modelling through graphical structure analysis (to be released)
##' }
##' @examples
##'
##' #Specifying a multivariate independent uniform prior
##'
##' rmviu <- function(n,a,b){
##' return(mapply(FUN = function(min.vec,max.vec,pn){stats::runif(pn,a,b)},
##'               min.vec=a,max.vec=b,MoreArgs = list(pn = n)))
##' }
##' dmviu <- function(x,a,b){
##' for(ii in 1:ncol(x)){
##'   x[,ii] <- dunif(x[,ii],a[ii],b[ii])
##' }
##' return(apply(x,1,prod))
##' }
##'
##' UNCOVER.opts(prior.override = TRUE,rprior = rmviu,
##'                  dprior = dmviu,a=rep(0,3),b=rep(1,3))
##'
##' \donttest{
##' # If we generate a co-variate matrix and binary response vector
##' CM <- matrix(rnorm(200),100,2)
##' rv <- sample(0:1,100,replace=TRUE)
##'
##' # We can then run our algorithm with a SMC threshold of 50 and a SMC cache
##' # checking threshold of 25 to see if this is quicker than the standard
##' # version
##' system.time(UNCOVER(X = CM,y = rv,verbose = FALSE))
##' system.time(UNCOVER(X = CM,y = rv,
##'                     options = UNCOVER.opts(SMC_thres = 50),
##'                     verbose = FALSE))
##' system.time(UNCOVER(X = CM,y = rv,
##'                     options = UNCOVER.opts(SMC_thres = 50,
##'                                            SMC_memo_thres = 25),
##'                     verbose = FALSE))
##' }
##'

UNCOVER.opts <- function(N = 1000,train_frac = 1,max_K = Inf,min_size = 0,
                        reg = 0,n_min_class = 0,SMC_thres = 30,
                        BIC_memo_thres = Inf,SMC_memo_thres = Inf,ess = N/2,
                        n_move = 1,prior.override = FALSE,rprior = NULL,
                        dprior = NULL,diagnostics = TRUE,
                        RIBIS_thres = 30,
                        BIC_cache = cachem::cache_mem(max_size = 1024 * 1024^2,
                                                      evict = "lru"),
                        SMC_cache = cachem::cache_mem(max_size = 1024 * 1024^2,
                                                      evict = "lru"),...){
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
  if(train_frac <=0 | train_frac > 1){
    stop("train_frac must be between 0 and 1 if validation data is required")
  }
  if(max_K < 1){
    stop("There must be at least 1 cluster")
  }
  if(min_size < 0){
    stop("min_size must be non-negative")
  }
  if(reg < 0){
    stop("reg must be non-negative")
  }
  if(n_min_class < 0){
    stop("n_min_class must be non-negative")
  }
  if(BIC_memo_thres < 0 ){
    stop("BIC_memo_thres must be non-negative")
  }
  if(SMC_memo_thres < 0 ){
    stop("SMC_memo_thres must be non-negative")
  }
  if(ess > N | ess < 0){
    stop("Effective sample size must be between 0 and N")
  }
  if(n_move < 1){
    stop("n_move must be an integer greater or equal than 1")
  }
  if(!inherits(BIC_cache,"cache_mem") | !inherits(SMC_cache,"cache_mem")){
    stop("Cache objects must be of class cache_mem")
  }
  if(SMC_thres < 0){
    stop("SMC_thres must be non-negative")
  }
  if(!is.logical(diagnostics)){
    stop("diagnostics must be logical")
  }
  return(list(N = N,train_frac = train_frac,max_K = max_K,min_size = min_size,
              reg = reg,n_min_class = n_min_class,SMC_thres = SMC_thres,
              BIC_memo_thres = BIC_memo_thres,
              SMC_memo_thres = SMC_memo_thres,ess = ess,n_move = n_move,
              rprior = rprior,dprior = dprior,prior.override = prior.override,
              diagnostics = diagnostics, RIBIS_thres = RIBIS_thres,
              BIC_cache = BIC_cache,SMC_cache = SMC_cache,MoreArgs = MoreArgs))
}

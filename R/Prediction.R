UNCOVER.infer <- function(x){
  if(!inherits(x,"UNCOVER")){
    stop("This function is only for outputs of UNCOVER")
  }
  if(x$Deforestation_Criterion=="Validation"){
    x$Model <- x$Model$All_Data
  }
  res.list <- vector(mode = "list",length = x$Model$Number_of_Clusters)
  for(k in 1:x$Model$Number_of_Clusters){
    if(typeof(x$Prior_Mean)=="character"){
      x$Prior_Mean = rep(0, ncol(x$Covariate_Matrix) + 1)
      x$Prior_Variance = diag(ncol(x$Covariate_Matrix) + 1)
    }
    res.list[[k]] <- IBIS.logreg(X = x$Covariate_Matrix[which(x$Model$Cluster_Allocation==k),,drop=FALSE],
                                 y = x$Response_Vector[which(x$Model$Cluster_Allocation==k)],
                                 options = do.call(IBIS.logreg.opts,
                                                   do.call(c,list(list(N=x$Control$N),
                                                                  ess = x$Control$ess,
                                                                  n_move = x$Control$n_move,
                                                                  prior.override = x$Control$prior.override,
                                                                  rprior = x$Control$rprior,
                                                                  dprior = x$Control$dprior,
                                                                  x$Control$MoreArgs))),
                                 prior_mean = x$Prior_Mean,
                                 prior_var = x$Prior_Variance)$samples
  }
  names(res.list) <- paste0("Cluster ",1:x$Model$Number_of_Clusters)
  return(res.list)
}

UNCOVER.assign <- function(x,nX){
  if(!inherits(x,"UNCOVER")){
    stop("This function is only for outputs of UNCOVER")
  }
  if(x$Deforestation_Criterion=="Validation"){
    x$Model <- x$Model$All_Data
  }
  nX <- as.matrix(nX,ncol = ncol(x$Covariate_Matrix))
  conn <- as.matrix(stats::dist(rbind(as.matrix(x$Covariate_Matrix),nX),method = "euclidean"))
  nn <- function(u){
    nns <- which(u==min(u))
    if(length(nns)==1){
      return(nns)
    } else{
      return(sample(nns,1))
    }
  }
  if(nrow(nX)==1){
    n.assign <- nn(conn[(nrow(x$Covariate_Matrix)+1):nrow(conn),1:nrow(x$Covariate_Matrix)])
  } else{
    n.assign <- apply(conn[(nrow(x$Covariate_Matrix)+1):nrow(conn),1:nrow(x$Covariate_Matrix)],1,nn)
  }
  c.assign <- x$Model$Cluster_Allocation[n.assign]
  return(c.assign)
}

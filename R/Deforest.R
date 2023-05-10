deforest.noc <- function(obs,res,gra,lbe,eps,K_dag,clu_al=NULL,c_s=NULL,
                         est_thres=30,mtb = Inf,mts = Inf,par_no=1000,rfun=NULL,
                         pdf_fun=NULL,efsamp = par_no/2,methas = 1,vb = FALSE,
                         cb,cs,PA,diagnostics = FALSE,Tr=NULL,SMC_f,BIC_f,rt = 30){
  K <- igraph::count_components(gra)
  if(is.null(clu_al)){
    clu_al <- igraph::components(gra)$membership
  }
  if(is.null(c_s)){
    c_s <- vector(mode="list",length=nrow(eps))
  }
  j <- 1
  if(vb){
    message("")
    message(crayon::green("Assessing the reintroduction of edges which have been removed"))
  }
  while(j <= nrow(eps)){
    edge_clu_al <- sort(clu_al[match(eps[j,],igraph::V(gra)$name)])
    if(!identical(c_s[[j]][[1]],which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]))){
      c_s[[j]] <- list(which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]),
                       lbe.gen(thres = est_thres,obs_mat = obs,res_vec = res,
                               obs_ind = which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]),
                               memo_thres_bic = mtb,memo_thres_smc = mts,
                               p_num = par_no,rpri = rfun,p_pdf = pdf_fun,
                               efs = efsamp,nm = methas,cache_bic = cb,
                               cache_smc = cs,MA = PA,SMC_fun = SMC_f,
                               BIC_fun = BIC_f, ribis_thres = rt))
    }
    if(vb){
      utils::txtProgressBar(min=1,max=nrow(eps)+1,initial=j+1,style=3)
    }
    if(j==nrow(eps)){
      lbe_comb.1 <- sapply(c_s,FUN = function(u){u[[2]]})
      lbe_comb.2 <- apply(eps,MARGIN = 1,FUN = function(u,lZ,ca,cg){sum(lZ[-ca[match(u,igraph::V(cg)$name)]])},lZ = lbe,ca = clu_al,cg = gra)
      lbe_comb <- lbe_comb.1 + lbe_comb.2
      if(K>K_dag | any(lbe_comb>sum(lbe))){
        cut_comb <- which.max(lbe_comb)
        edge_clu_al <- sort(clu_al[match(eps[cut_comb,],igraph::V(gra)$name)])
        if(vb){
          message("")
          message(crayon::blue(paste0("Combining connected components ",edge_clu_al[1]," and ",edge_clu_al[2])))
        }
        clu_al[which(clu_al==edge_clu_al[2])] <- edge_clu_al[1]
        for(k in (edge_clu_al[2]+1):K){
          clu_al[which(clu_al==k)] <- k-1
        }
        lbe[edge_clu_al[1]] <- c_s[[cut_comb]][[2]]
        lbe <- lbe[-edge_clu_al[2]]
        gra <- igraph::add_edges(gra,eps[cut_comb,])
        if(diagnostics){
          Tr <- rbind(Tr,c(paste0("Def.Add.",eps[cut_comb,1],"-",eps[cut_comb,2]),sum(lbe),K-1))
        }
        eps <- eps[-cut_comb,,drop=FALSE]
        c_s[[cut_comb]] <- c()
        K <- K-1
        j <- 0
        if(nrow(eps)!=0 & vb){
          message("")
          message(crayon::green("Reassessing the reintroduction of edges which have been removed"))
        }
      }
    }
    j <- j+1
  }
  if(diagnostics){
    return(list(Cluster_Allocation = clu_al,
                Log_Marginal_Likelihoods = lbe,
                Graph = gra,
                Number_of_Clusters = K,
                Edges_Removed = eps,
                Diagnostics = Tr))
  } else{
    return(list(Cluster_Allocation = clu_al,
                Log_Marginal_Likelihoods = lbe,
                Graph = gra,
                Number_of_Clusters = K,
                Edges_Removed = eps))
  }
}

deforest.soc <- function(obs,res,gra,lbe,eps,n_dag,clu_al=NULL,c_s=NULL,
                         est_thres=30,mtb = Inf,mts = Inf,par_no=1000,rfun=NULL,
                         pdf_fun=NULL,efsamp = par_no/2,methas = 1,vb = FALSE,
                         cb,cs,PA,diagnostics=FALSE,Tr=NULL,SMC_f,BIC_f,
                         rt = 30){
  K <- igraph::count_components(gra)
  if(is.null(clu_al)){
    clu_al <- igraph::components(gra)$membership
  }
  if(is.null(c_s)){
    c_s <- vector(mode="list",length=nrow(eps))
  }
  j <- 1
  if(vb){
    message("")
    message(crayon::green("Assessing the reintroduction of edges which have been removed"))
  }
  while(j <= nrow(eps)){
    edge_clu_al <- sort(clu_al[match(eps[j,],igraph::V(gra)$name)])
    if(!identical(c_s[[j]][[1]],which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]))){
      c_s[[j]] <- list(which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]),
                       lbe.gen(thres = est_thres,obs_mat = obs,res_vec = res,
                               obs_ind = which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]),
                               memo_thres_bic = mtb,memo_thres_smc = mts,
                               p_num = par_no,rpri = rfun,p_pdf = pdf_fun,
                               efs = efsamp,nm = methas,cache_bic = cb,
                               cache_smc = cs,MA = PA,SMC_fun = SMC_f,
                               BIC_fun = BIC_f, ribis_thres = rt))
    }
    if(vb){
      utils::txtProgressBar(min=1,max=nrow(eps)+1,initial=j+1,style=3)
    }
    if(j==nrow(eps)){
      lbe_comb.1 <- sapply(c_s,FUN = function(u){u[[2]]})
      lbe_comb.2 <- apply(eps,MARGIN = 1,FUN = function(u,lZ,ca,cg){sum(lZ[-ca[match(u,igraph::V(cg)$name)]])},lZ = lbe, ca = clu_al,cg = gra)
      lbe_comb <- lbe_comb.1 + lbe_comb.2
      small_clu <- which(table(clu_al)<n_dag)
      small_bool <- apply(eps,MARGIN = 1,FUN = function(u,ca,sc,cg){ca[match(u,igraph::V(cg)$name)[1]]%in%sc | ca[match(u,igraph::V(cg)$name)[2]]%in%sc},ca=clu_al,sc=small_clu,cg=gra)
      small_bool <- small_bool | lbe_comb>sum(lbe)
      if(any(small_bool)){
        lbe_comb[which(!small_bool)] <- -Inf
        cut_comb <- which.max(lbe_comb)
        edge_clu_al <- sort(clu_al[match(eps[cut_comb,],igraph::V(gra)$name)])
        if(vb){
          message("")
          message(crayon::blue(paste0("Combining connected components ",edge_clu_al[1]," and ",edge_clu_al[2])))
        }
        clu_al[which(clu_al==edge_clu_al[2])] <- edge_clu_al[1]
        for(k in (edge_clu_al[2]+1):K){
          clu_al[which(clu_al==k)] <- k-1
        }
        lbe[edge_clu_al[1]] <- c_s[[cut_comb]][[2]]
        lbe <- lbe[-edge_clu_al[2]]
        gra <- igraph::add_edges(gra,eps[cut_comb,])
        if(diagnostics){
          tab_clu_al <- table(clu_al)
          Tr <- rbind(Tr,c(paste0("Def.Add.",eps[cut_comb,1],"-",eps[cut_comb,2]),
                           sum(lbe),min(tab_clu_al),sum(tab_clu_al<n_dag)))
        }
        eps <- eps[-cut_comb,,drop=FALSE]
        c_s[[cut_comb]] <- c()
        K <- K-1
        j <- 0
        if(nrow(eps)!=0 & vb){
          message("")
          message(crayon::green("Reassessing the reintroduction of edges which have been removed"))
        }
      }
    }
    j <- j+1
  }
  if(diagnostics){
    return(list(Cluster_Allocation = clu_al,
                Log_Marginal_Likelihoods = lbe,
                Graph = gra,
                Number_of_Clusters = K,
                Edges_Removed = eps,
                Diagnostics = Tr))
  } else{
    return(list(Cluster_Allocation = clu_al,
                Log_Marginal_Likelihoods = lbe,
                Graph = gra,
                Number_of_Clusters = K,
                Edges_Removed = eps))
  }
}

deforest.maxreg <- function(obs,res,gra,lbe,eps,tau,clu_al=NULL,c_s=NULL,
                            est_thres=30,mtb = Inf,mts = Inf,par_no=1000,
                            rfun=NULL,pdf_fun=NULL,efsamp = par_no/2,methas = 1,
                            vb = FALSE,cb,cs,PA,diagnostics=FALSE,Tr=NULL,
                            SMC_f,BIC_f,rt = 30){
  K <- igraph::count_components(gra)
  if(is.null(clu_al)){
    clu_al <- igraph::components(gra)$membership
  }
  if(is.null(c_s)){
    c_s <- vector(mode="list",length=nrow(eps))
  }
  j <- 1
  if(vb){
    message("")
    message(crayon::green("Assessing the reintroduction of edges which have been removed"))
  }
  while(j <= nrow(eps)){
    edge_clu_al <- sort(clu_al[match(eps[j,],igraph::V(gra)$name)])
    if(!identical(c_s[[j]][[1]],which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]))){
      c_s[[j]] <- list(which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]),
                       lbe.gen(thres = est_thres,obs_mat = obs,res_vec = res,
                               obs_ind = which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]),
                               memo_thres_bic = mtb,memo_thres_smc = mts,
                               p_num = par_no,rpri = rfun,p_pdf = pdf_fun,
                               efs = efsamp,nm = methas,cache_bic = cb,
                               cache_smc = cs,MA = PA,SMC_fun = SMC_f,
                               BIC_fun = BIC_f,ribis_thres = rt))
    }
    if(vb){
      utils::txtProgressBar(min=1,max=nrow(eps)+1,initial=j+1,style=3)
    }
    if(j==nrow(eps)){
      lbe_comb.1 <- sapply(c_s,FUN = function(u){u[[2]]})
      lbe_comb.2 <- apply(eps,MARGIN = 1,FUN = function(u,lZ,ca,cg){sum(lZ[-ca[match(u,igraph::V(cg)$name)]])},lZ = lbe,ca = clu_al,cg = gra)
      lbe_comb <- lbe_comb.1 + lbe_comb.2
      if(any(tau + lbe_comb>sum(lbe))){
        cut_comb <- which.max(lbe_comb)
        edge_clu_al <- sort(clu_al[match(eps[cut_comb,],igraph::V(gra)$name)])
        if(vb){
          message("")
          message(crayon::blue(paste0("Combining connected components ",edge_clu_al[1]," and ",edge_clu_al[2])))
        }
        clu_al[which(clu_al==edge_clu_al[2])] <- edge_clu_al[1]
        for(k in (edge_clu_al[2]+1):K){
          clu_al[which(clu_al==k)] <- k-1
        }
        lbe[edge_clu_al[1]] <- c_s[[cut_comb]][[2]]
        lbe <- lbe[-edge_clu_al[2]]
        gra <- igraph::add_edges(gra,eps[cut_comb,])
        if(diagnostics){
          Tr <- rbind(Tr,c(paste0("Def.Add.",eps[cut_comb,1],"-",eps[cut_comb,2]),sum(lbe)))
        }
        eps <- eps[-cut_comb,,drop=FALSE]
        c_s[[cut_comb]] <- c()
        K <- K-1
        j <- 0
        if(nrow(eps)!=0 & vb){
          message("")
          message(crayon::green("Reassessing the reintroduction of edges which have been removed"))
        }
      }
    }
    j <- j+1
  }
  if(diagnostics){
    return(list(Cluster_Allocation = clu_al,
                Log_Marginal_Likelihoods = lbe,
                Graph = gra,
                Number_of_Clusters = K,
                Edges_Removed = eps,
                Diagnostics = Tr))
  } else{
    return(list(Cluster_Allocation = clu_al,
                Log_Marginal_Likelihoods = lbe,
                Graph = gra,
                Number_of_Clusters = K,
                Edges_Removed = eps))
  }
}

deforest.validation <- function(obs,obs_all,res,res_all,gra,lbe,eps,
                                gra_all=NULL,which_tr=NULL,rho=NULL,clu_al=NULL,
                                c_s=NULL,est_thres=30,mtb=Inf,mts=Inf,
                                par_no=1000,rfun,pdf_fun,efsamp=par_no/2,
                                methas=1,vb=FALSE,cb,cs,PA,diagnostics=FALSE,
                                Tr=NULL,SMC_f,BIC_f,rt=30){
  if(diagnostics){
    if(is.null(Tr)){
      stop("If diagnostics=TRUE then Tr needs to be specified")
    } else{
      Tr <- data.frame(Tr,Log_Bayesian_Evidence_All = rep(NA,nrow(Tr)),
                       Robustness_Statistic = rep(NA,nrow(Tr)))
    }
  }
  K <- igraph::count_components(gra)
  if(is.null(clu_al)){
    clu_al <- igraph::components(gra)$membership
  }
  if(is.null(c_s)){
    c_s <- vector(mode="list",length=nrow(eps))
  }
  if(is.null(gra_all)){
    gra_val <- two.stage.mst(obs_mat = obs_all,tr_ind = which_tr,mst_sub = rho)
    gra_all <- igraph::delete_edges(gra_val[[2]],igraph::E(gra_val[[2]])[igraph::get.edge.ids(gra_val[[2]],c(t(eps)))])
  }
  clu_al_all <- igraph::components(gra_all)$membership
  for(k in 1:K){
    clu_al_all[which(clu_al_all==clu_al_all[as.numeric(igraph::V(gra)$name[which(clu_al==k)[1]])])] <- K+k
  }
  clu_al_all <- clu_al_all - K
  lbe_all <- lbe
  for(k in 1:K){
    lbe_all[k] <- lbe.gen(thres = est_thres,obs_mat = obs_all,res_vec = res_all,
                          obs_ind = which(clu_al_all==k),memo_thres_bic = mtb,
                          memo_thres_smc = mts,p_num = par_no,rpri = rfun,
                          p_pdf = pdf_fun,efs = efsamp,nm = methas,
                          cache_bic = cb,cache_smc = cs,MA = PA,SMC_fun = SMC_f,
                          BIC_fun = BIC_f,ribis_thres = rt)
  }
  RobS <- sum(lbe_all-lbe)
  if(diagnostics){
    Tr[nrow(Tr),3:4] <- c(sum(lbe_all),RobS)
  }
  j <- 1
  c_s_all <- vector(mode="list",length=nrow(eps))
  if(vb){
    message("")
    message(crayon::green("Assessing the reintroduction of edges which have been removed"))
  }
  while(j <= nrow(eps)){
    edge_clu_al <- sort(clu_al_all[as.numeric(eps[j,])])
    if(!identical(c_s[[j]][[1]],which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]))){
      c_s[[j]] <- list(which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]),
                       lbe.gen(thres = est_thres,obs_mat = obs,res_vec = res,
                               obs_ind = which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]),
                               memo_thres_bic = mtb,memo_thres_smc = mts,
                               p_num = par_no,rpri = rfun,p_pdf = pdf_fun,
                               efs = efsamp,nm = methas,cache_bic = cb,
                               cache_smc = cs,MA = PA,SMC_fun = SMC_f,
                               BIC_fun = BIC_f,ribis_thres = rt))
    }
    if(!identical(c_s_all[[j]][[1]],which(clu_al_all==edge_clu_al[1]| clu_al_all==edge_clu_al[2]))){
      c_s_all[[j]] <- list(which(clu_al_all==edge_clu_al[1]| clu_al_all==edge_clu_al[2]),
                           lbe.gen(thres = est_thres,obs_mat = obs_all,
                                   res_vec = res_all,
                                   obs_ind = which(clu_al_all==edge_clu_al[1]| clu_al_all==edge_clu_al[2]),
                                   memo_thres_bic = mtb,memo_thres_smc = mts,
                                   p_num = par_no,rpri = rfun,p_pdf = pdf_fun,
                                   efs = efsamp,nm = methas,cache_bic = cb,
                                   cache_smc = cs,MA = PA,SMC_fun = SMC_f,
                                   BIC_fun = BIC_f,ribis_thres = rt))
    }
    if(vb){
      utils::txtProgressBar(min=1,max=nrow(eps)+1,initial=j+1,style=3)
    }
    if(j==nrow(eps)){
      lbe_comb.1 <- sapply(c_s,FUN = function(u){u[[2]]})
      lbe_comb.2 <- apply(eps,MARGIN = 1,FUN = function(u,lZ,ca,cg){sum(lZ[-ca[match(u,igraph::V(cg)$name)]])},lZ = lbe,ca = clu_al,cg = gra)
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
          message(crayon::blue(paste0("Combining connected components ",edge_clu_al[1]," and ",edge_clu_al[2])))
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
        gra <- igraph::add_edges(gra,eps[cut_comb,])
        gra_all <- igraph::add_edges(gra_all,eps[cut_comb,])
        if(diagnostics){
          Tr <- rbind(Tr,c(paste0("Def.Add.",eps[cut_comb,1],"-",eps[cut_comb,2]),sum(lbe),sum(lbe_all),RobS))
        }
        eps <- eps[-cut_comb,,drop=FALSE]
        c_s[[cut_comb]] <- c()
        c_s_all[[cut_comb]] <- c()
        K <- K-1
        j <- 0
        if(nrow(eps)!=0 & vb){
          message("")
          message(crayon::green("Reassessing the reintroduction of edges which have been removed"))
        }
      }
    }
    j <- j+1
  }
  TD = list(Cluster_Allocation = clu_al,
                       Log_Marginal_Likelihoods = lbe,
                       Graph = gra,
                       Number_of_Clusters = K,
                       Edges_Removed = eps)
  if(diagnostics){
    AD = list(Cluster_Allocation = clu_al_all,
                    Log_Marginal_Likelihoods = lbe_all,
                    Graph = gra_all,
                    Number_of_Clusters = K,
                    Edges_Removed = eps,
                    Diagnostics = Tr)
  } else{
    AD = list(Cluster_Allocation = clu_al_all,
                    Log_Marginal_Likelihoods = lbe_all,
                    Graph = gra_all,
                    Number_of_Clusters = K,
                    Edges_Removed = eps)
  }
  return(list(Training_Data = TD,
              All_Data = AD))
}

deforest.diverse <- function(obs,res,gra,lbe,eps,ups,clu_al=NULL,c_s=NULL,
                              est_thres=30,mtb = Inf,mts = Inf,par_no=1000,
                              rfun=NULL,pdf_fun=NULL,efsamp=par_no/2,methas=1,
                              vb = FALSE,cb,cs,PA,diagnostics=FALSE,Tr=NULL,
                              SMC_f,BIC_f,rt=30){
  K <- igraph::count_components(gra)
  if(is.null(clu_al)){
    clu_al <- igraph::components(gra)$membership
  }
  if(is.null(c_s)){
    c_s <- vector(mode="list",length=nrow(eps))
  }
  j <- 1
  if(vb){
    message("")
    message(crayon::green("Assessing the reintroduction of edges which have been removed"))
  }
  while(j <= nrow(eps)){
    edge_clu_al <- sort(clu_al[match(eps[j,],igraph::V(gra)$name)])
    if(!identical(c_s[[j]][[1]],which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]))){
      c_s[[j]] <- list(which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]),
                       lbe.gen(thres = est_thres,obs_mat = obs,res_vec = res,
                               obs_ind = which(clu_al==edge_clu_al[1]| clu_al==edge_clu_al[2]),
                               memo_thres_bic = mtb,memo_thres_smc = mts,
                               p_num = par_no,rpri = rfun,p_pdf = pdf_fun,
                               efs = efsamp,nm = methas,cache_bic = cb,
                               cache_smc = cs,MA = PA,SMC_fun = SMC_f,
                               BIC_fun = BIC_f,ribis_thres = rt))
    }
    if(vb){
      utils::txtProgressBar(min=1,max=nrow(eps)+1,initial=j+1,style=3)
    }
    if(j==nrow(eps)){
      lbe_comb.1 <- sapply(c_s,FUN = function(u){u[[2]]})
      lbe_comb.2 <- apply(eps,MARGIN = 1,FUN = function(u,lZ,ca,cg){sum(lZ[-ca[match(u,igraph::V(cg)$name)]])},lZ = lbe, ca = clu_al,cg = gra)
      lbe_comb <- lbe_comb.1 + lbe_comb.2
      undiv_clu <- c()
      for(k in 1:K){
        if(any(table(factor(res[which(clu_al==k)],levels = 0:1))<ups)){
          undiv_clu <- c(undiv_clu,k)
        }
      }
      undiv_bool <- apply(eps,MARGIN = 1,FUN = function(u,ca,uc,cg){ca[match(u,igraph::V(cg)$name)[1]]%in%uc | ca[match(u,igraph::V(cg)$name)[2]]%in%uc},ca=clu_al,uc=undiv_clu,cg=gra)
      undiv_bool <- undiv_bool | lbe_comb>sum(lbe)
      if(any(undiv_bool)){
        lbe_comb[which(!undiv_bool)] <- -Inf
        cut_comb <- which.max(lbe_comb)
        edge_clu_al <- sort(clu_al[match(eps[cut_comb,],igraph::V(gra)$name)])
        if(vb){
          message("")
          message(crayon::blue(paste0("Combining connected components ",edge_clu_al[1]," and ",edge_clu_al[2])))
        }
        clu_al[which(clu_al==edge_clu_al[2])] <- edge_clu_al[1]
        for(k in (edge_clu_al[2]+1):K){
          clu_al[which(clu_al==k)] <- k-1
        }
        lbe[edge_clu_al[1]] <- c_s[[cut_comb]][[2]]
        lbe <- lbe[-edge_clu_al[2]]
        gra <- igraph::add_edges(gra,eps[cut_comb,])
        if(diagnostics){
          sap_div <- sapply(1:(K-1),FUN = function(u,y,z){min(table(factor(y[which(z==u)],levels=0:1)))},y=res,z=clu_al)
          Tr <- rbind(Tr,c(paste0("Def.Add.",eps[cut_comb,1],"-",eps[cut_comb,2]),sum(lbe),min(sap_div),sum(sap_div<ups)))
        }
        eps <- eps[-cut_comb,,drop=FALSE]
        c_s[[cut_comb]] <- c()
        K <- K-1
        j <- 0
        if(nrow(eps)!=0 & vb){
          message("")
          message(crayon::green("Reassessing the reintroduction of edges which have been removed"))
        }
      }
    }
    j <- j+1
  }
  if(diagnostics){
    return(list(Cluster_Allocation = clu_al,
                Log_Marginal_Likelihoods = lbe,
                Graph = gra,
                Number_of_Clusters = K,
                Edges_Removed = eps,
                Diagnostics = Tr))
  } else{
    return(list(Cluster_Allocation = clu_al,
                Log_Marginal_Likelihoods = lbe,
                Graph = gra,
                Number_of_Clusters = K,
                Edges_Removed = eps))
  }
}

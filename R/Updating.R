remove.edge <- function(gra,j,clu_al=NULL,lbe,obs,res,est_thres=30,mtb=Inf,
                        mts=Inf,par_no=1000,rfun,pdf_fun,efsamp=par_no/2,
                        methas=1,cb,cs,PA,SMC_f,BIC_f){
  K <- igraph::count_components(gra)
  if(is.null(clu_al)){
    clu_al <- igraph::components(gra)$membership
  }
  k <- clu_al[igraph::get.edgelist(gra,names=FALSE)[j,]][1]
  gra_rem <- igraph::delete_edges(gra,igraph::get.edge.ids(gra,igraph::get.edgelist(gra)[j,]))
  clu_al_rem <- igraph::components(gra_rem)$membership
  change_set <- which(clu_al_rem==clu_al_rem[igraph::get.edgelist(gra,names=FALSE)[j,2]])
  clu_al[change_set] <- K+1
  for(l in c(k,K+1)){
    lbe[l] <- lbe.gen(thres = est_thres,obs_mat = obs,res_vec = res,
                      obs_ind = which(clu_al==l),memo_thres_bic = mtb,
                      memo_thres_smc = mts,p_num = par_no,rpri = rfun,
                      p_pdf = pdf_fun,efs = efsamp,nm = methas,
                      cache_bic = cb,cache_smc = cs,MA = PA,SMC_fun = SMC_f,
                      BIC_fun = BIC_f)
  }
  return(list(clu_al,lbe,gra_rem))
}

###############################################
## R script for functions in UNCOVER package ##
###############################################
##
## Sam Emerson
## September 2022
##


#######################################################################################################################
## Functions for Extracting the Log Bayesian Evidence of a Model Obtained by Removing or Adding an Edge to the Graph ##
#######################################################################################################################


##' Log Bayesian evidence estimation of a model obtained by removing an edge from the models associated graph
##'
##'
##' @export
##' @name remove.edge
##' @description This function removes an edge from the graph, creating a further partition of the data, then calculates the Bayesian evidence of the new model created by this removal.
##'
##' Used in UNCOVER to check if the removal of certain edges is beneficial to the model.
##'
##' @keywords edge removal
##' @param gra `igraph` object which contains the information of the graph of the current model
##' @param j Index of the edge which will be removed from the graph
##' @param clu_al A vector detailing the cluster allocation of each observation. If not specified the function will generate this vector.
##' @param lbe A vector detailing the log Bayesian evidences of all the sub-models defined by the separated components of `gra`
##' @param obs Covariate matrix
##' @param res Binary response vector
##' @param est_method Method to be used to estimate the log Bayesian evidence of the newly formed sub-models, can be one of `"SMC"`,`"SMC_BIC"` or `"BIC"`
##' @param est_thres The threshold for which the number of observations needs to exceed to consider using BIC as an estimator. Only applies if method `"SMC_BIC"` is selected. Defaults to infinity if not specified.
##' @param par_no Number of samples of the prior used for the SMC sampler. Not required if method `"BIC"` selected.
##' @param rfun Function to sample from the prior. Must only have two arguments, `p_n` and `di` (Number of prior samples to generate and the number of dimensions of a single sample respectively).
##' @param pdf_fun Probability Density Function of the prior. Must only have two arguments, `th` and `di` (a vector or matrix of regression coefficients samples and the number of dimensions of a single sample respectively).
##' @return A list consisting of; the cluster allocation vector of the new model, the resulting Bayesian evidence vector for the new model and an `igraph` object containing information on the new graph.
##' @details Requires a minimum spanning forest graph which defines components for a multiplicative Bayesian logistic regression model. `lbe.gen` should be used to obtain the log Bayesian evidence of this initial model before this function is used.
##'
##' If the clusters specified by the initial model have fixed labels then this should be specified by `clu_al`. `clu_al` must be a one of the possible labelling of the observations defined by the clusters of the graph. For example for a graph where there is only one connected component, if `clu_al` is specified it cannot be anything other than `rep(1,length(res)`.
##'
##' For more details on the specifics of the possible values for `est_method`, see the help page of the function `lbe.gen`.
##'
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
##' # If we assume the prior for the regression coefficients is a standard
##' # normal
##' pr_samp <- function(p_n,di){return(rmvn(p_n,rep(0,di),diag(di)))}
##' pr_fun <- function(th,di){return(dmvn(th,mu=rep(0,di),sigma=diag(di)))}
##'
##' # Then we can generate the log Bayesian evidence for this one component
##' # model
##' lZ <- lbe.gen(method = "SMC_BIC",thres = 30,obs_mat = CM,res_vec = rv,p_num = 500, rpri = pr_samp, p_pdf = pr_fun)
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
##' # If we then wanted to remove the first edge of the newly formed graph
##' # gra_CM.2 we would have
##' er.2 <- remove.edge(gra = gra_CM.2, j = 1, clu_al = clu_al.2, lbe = lZ.2, obs = CM, res = rv, est_method = "SMC_BIC", par_no = 500, rfun = pr_samp, pdf_fun = pr_fun)
##' plot(er.2[[3]],layout=CM)
##' er.2[[2]]
##' sum(er.2[[2]])

remove.edge <- function(gra,j,clu_al=NULL,lbe,obs,res,est_method="SMC_BIC",est_thres=Inf,par_no=NULL,rfun=NULL,pdf_fun=NULL){
  K <- igraph::count_components(gra)
  if(is.null(clu_al)){
    clu_al <- igraph::components(gra)$membership
  }
  k <- clu_al[igraph::get.edgelist(gra,names=F)[j,]][1]
  gra_rem <- igraph::delete_edges(gra,E(gra)[.env$j])
  clu_al_rem <- igraph::components(gra_rem)$membership
  change_set <- which(clu_al_rem==clu_al_rem[igraph::get.edgelist(gra,names=F)[j,2]])
  clu_al[change_set] <- K+1
  for(l in c(k,K+1)){
    lbe[l] <- lbe.gen(method = est_method,thres = est_thres,obs_mat = obs[which(clu_al==l),,drop=F],res_vec = res[which(clu_al==l)],p_num = par_no,rpri = rfun,p_pdf = pdf_fun)
  }
  return(list(clu_al,lbe,gra_rem))
}

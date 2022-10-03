################################################
## R script for functions in UNCOVER package ##
###############################################
##
## Sam Emerson
## September 2022
##


############################################################
## Function for Generating Types of Minimum Spanning Tree ##
############################################################


##' One stage minimum spanning tree generation
##'
##'
##' @export
##' @name one.stage.mst
##' @description Generation of a minimum spanning tree induced subgraph using the covariate matrix provided.
##'
##' Used in UNCOVER to generate the initial graph for the algorithm if the pruning method does not involve validation data.
##'
##' @keywords minimum spanning tree
##' @param obs Covariate matrix
##' @param rho Vector of indices which indicate the variables of `obs` to use to construct the minimum spanning tree
##' @return An `igraph` object of the minimum spanning tree induced subgraph
##' @examples
##'
##' # First we generate a covariate matrix `obs`
##' CM <- matrix(rnorm(300),100,3)
##'
##' # Assuming we require a minimum spanning tree induced subgraph using all
##' # variables
##' one.stage.mst(obs = CM)
##'
##' # If we only require the minimum spanning tree to be constructed from the
##' # first two variables
##' one.stage.mst(obs = CM,rho = 1:2)
##'

one.stage.mst <- function(obs,rho=NULL){
  if(is.null(rho)){
    rho <- 1:ncol(obs)
  }
  conn <- as.matrix(dist(obs[,rho],method = "euclidean"))
  g <- igraph::graph_from_adjacency_matrix(conn,weighted=T,mode="undirected")
  g <- igraph::mst(g,algorithm="prim")
  return(g)
}

##' Two stage minimum spanning tree generation
##'
##'
##' @export
##' @name two.stage.mst
##' @description Generation of a two stage minimum spanning tree induced subgraph using the covariate matrix provided and an indictor of which rows are to be used as training data.
##'
##' Used in UNCOVER to generate the initial graph for the algorithm if the pruning method involves validation data.
##'
##' @keywords minimum spanning tree
##' @param obs_mat Covariate matrix
##' @param tr_ind Vector of indices of the observations that will be used as training data
##' @param mst_sub Vector of indices which indicate the variables of `obs_mat` to use to construct the minimum spanning tree
##' @return A list of two `igraph` objects; the minimum spanning tree induced subgraph of the training data and the two stage minimum spanning tree induced subgraph obtained by adding the validation data.
##' @details A minimum spanning tree is first constructed on the training data, these edges are then fixed. Validation data is then added in the second stage and we continue with Prim's algorithm of constructing a minimum spanning tree (given the current edges) until all observations are connected.
##' @examples
##'
##' # First we generate a covariate matrix `obs_mat` and assign observations as
##' # training data or validation data
##' CM <- matrix(rnorm(300),100,3)
##' ti <- sort(sample(1:100,50))
##'
##' # Assuming we require a two stage minimum spanning tree induced subgraph
##' # using all variables
##' tsm <- two.stage.mst(obs_mat = CM,tr_ind = ti)
##'
##' # The edges from the first graph should also be present in the second graph
##' E(tsm[[1]])
##' E(tsm[[2]])
##'
##' # If we only require the two stage minimum spanning tree to be constructed
##' # from the first two variables
##' tsm.2 <- two.stage.mst(obs_mat = CM,tr_ind = ti,mst_sub = 1:2)
##' E(tsm.2[[1]])
##' E(tsm.2[[2]])
##'

two.stage.mst <- function(obs_mat,tr_ind,mst_sub=NULL){
  if(length(tr_ind)==nrow(obs_mat)){
    stop("There is no validation data to produce a two stage minimum spanning tree. Consider using one.stage.mst instead.")
  }
  if(is.null(mst_sub)){
    mst_sub <- 1:ncol(obs_mat)
  }
  conn <- as.matrix(dist(obs_mat[,mst_sub],method = "euclidean"))
  conn_temp <- conn
  conn_tr <- conn[tr_ind,tr_ind]
  g_tr <- igraph::graph_from_adjacency_matrix(conn_tr,weighted=T,mode="undirected")
  g_tr <- igraph::mst(g_tr,algorithm="prim")
  igraph::V(g_tr)$name <- as.character(tr_ind)
  e <- max(igraph::E(g_tr)$weight)
  conn_temp <- conn_temp + e
  edg_ind <- igraph::get.edgelist(g_tr)
  class(edg_ind) <- "numeric"
  conn_temp[edg_ind] <- conn_temp[edg_ind] - e
  conn_temp[edg_ind[,2:1]] <- conn_temp[edg_ind[,2:1]] - e
  diag(conn_temp) <- 0
  g <- igraph::graph_from_adjacency_matrix(conn_temp,weighted = T,mode="undirected")
  g <- igraph::mst(g,algorithm = "prim")
  igraph::E(g)$weight <- conn[igraph::get.edgelist(g,names=F)]
  return(list(g_tr,g))
}

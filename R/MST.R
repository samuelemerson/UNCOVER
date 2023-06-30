one.stage.mst <- function(obs,rho=NULL){
  if(is.null(rho)){
    rho <- 1:ncol(obs)
  }
  conn <- as.matrix(stats::dist(obs[,rho],method = "euclidean"))
  beid <- min(conn[conn>0])
  if(any(which(conn==0,arr.ind = T)[,1] != which(conn==0,arr.ind = T)[,2])){
    conn[conn==0] <- beid/2
    diag(conn) <- 0
  }
  g <- igraph::graph_from_adjacency_matrix(conn,weighted=TRUE,mode="undirected")
  g <- igraph::mst(g,algorithm="prim")
  return(list(g,beid))
}

two.stage.mst.defunct <- function(obs_mat,tr_ind,mst_sub=NULL){
  if(length(tr_ind)==nrow(obs_mat)){
    stop("There is no validation data to produce a two stage minimum spanning tree. Consider using one.stage.mst instead.")
  }
  if(is.null(mst_sub)){
    mst_sub <- 1:ncol(obs_mat)
  }
  conn <- as.matrix(stats::dist(obs_mat[,mst_sub],method = "euclidean"))
  conn_temp <- conn
  conn_tr <- conn[tr_ind,tr_ind]
  g_tr <- igraph::graph_from_adjacency_matrix(conn_tr,weighted=TRUE,mode="undirected")
  g_tr <- igraph::mst(g_tr,algorithm="prim")
  igraph::V(g_tr)$name <- as.character(tr_ind)
  e <- max(igraph::E(g_tr)$weight)
  conn_temp <- conn_temp + e
  edg_ind <- igraph::get.edgelist(g_tr)
  class(edg_ind) <- "numeric"
  conn_temp[edg_ind] <- conn_temp[edg_ind] - e
  conn_temp[edg_ind[,2:1]] <- conn_temp[edg_ind[,2:1]] - e
  diag(conn_temp) <- 0
  g <- igraph::graph_from_adjacency_matrix(conn_temp,weighted = TRUE,mode="undirected")
  g <- igraph::mst(g,algorithm = "prim")
  igraph::E(g)$weight <- conn[igraph::get.edgelist(g,names=FALSE)]
  return(list(g_tr,g))
}

two.stage.mst <- function(obs_mat,tr_ind,mst_sub=NULL){
  if(length(tr_ind)==nrow(obs_mat)){
    stop("There is no validation data to produce a two stage minimum spanning tree. Consider using one.stage.mst instead.")
  }
  if(is.null(mst_sub)){
    mst_sub <- 1:ncol(obs_mat)
  }
  conn <- as.matrix(stats::dist(obs_mat[,mst_sub],method = "euclidean"))
  beid <- min(conn[conn>0])
  if(any(which(conn==0,arr.ind = T)[,1] != which(conn==0,arr.ind = T)[,2])){
    conn[conn==0] <- beid/2
    diag(conn) <- 0
  }
  conn_temp <- conn
  conn_tr <- conn[tr_ind,tr_ind]
  g_tr <- igraph::graph_from_adjacency_matrix(conn_tr,weighted=TRUE,mode="undirected")
  g_tr <- igraph::mst(g_tr,algorithm="prim")
  igraph::V(g_tr)$name <- as.character(tr_ind)
  val_ind <- setdiff(1:nrow(obs_mat),tr_ind)
  conn_temp[val_ind,val_ind] <- Inf
  diag(conn_temp) <- Inf
  add_edges_val <- matrix(0,nrow = length(val_ind),ncol = 2)
  add_edges_val[,1] <- as.character(val_ind)
  add_edges_val[,2] <- as.character(apply(conn_temp[val_ind,],1,which.min))
  conn_temp <- matrix(0,nrow(conn),ncol(conn),dimnames = dimnames(conn))
  conn_temp[igraph::get.edgelist(g_tr)] <- conn[igraph::get.edgelist(g_tr)]
  conn_temp[igraph::get.edgelist(g_tr)[,2:1]] <- conn[igraph::get.edgelist(g_tr)]
  conn_temp[add_edges_val] <- conn[add_edges_val]
  conn_temp[add_edges_val[,2:1]] <- conn[add_edges_val]
  g <- igraph::graph_from_adjacency_matrix(conn_temp,weighted = TRUE,mode="undirected")
  return(list(g_tr,g,beid))
}

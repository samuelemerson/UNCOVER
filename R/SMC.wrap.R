###############################################
## R script for functions in UNCOVER package ##
###############################################
##
## Sam Emerson
## January 2023
##


#################################################################################
## Function for Sequential Monte Carlo Sampler to Obtain the Bayesian Evidence ##
#################################################################################



##' Logistic regression iterated batch importance sampling
##'
##'
##' @export
##' @name IBIS.logreg
##' @description This function uses an Iterated Batch Importance Sampling (IBIS)
##' scheme with batch size one to go from prior to full posterior. We
##' assume a Bayesian logistic regression model.
##'
##' @keywords sequential monte carlo IBIS
##' @param X Co-variate matrix
##' @param y Binary response vector
##' @param options Additional arguments that can be specified for `IBIS.logreg`.
##' See [IBIS.logreg.opts()] for details. Can be ignored.
##' @param prior_mean Mean for the multivariate normal prior used in the SMC
##' sampler. See details. Defaults to the origin.
##' @param prior_var Variance matrix for the multivariate normal prior used in
##' the SMC sampler. See details. Defaults to the identity matrix.
##' @return An object of class `"IBIS"`, which is a list consisting of:
##'
##' \describe{
##' \item{`covariate_matrix`}{The co-variate matrix provided.}
##' \item{`response_vector`}{The binary response vector provided.}
##' \item{`samples`}{A matrix of samples from the posterior.}
##' \item{`log_Bayesian_evidence`}{An estimate of the log Bayesian evidence (or
##' normalisation constant) of the posterior.}
##' \item{`diagnostics`}{A data frame recording the features of the SMC sampler
##' as the observations were added.}
##' }
##'
##' If `weighted==TRUE` then an additional element of the list (`weights`) is
##' added detailing the weights of the posterior samples.
##' @details Details of the internal mechanisms of the SMC sampler such as the
##' Metropolis-Hastings MCMC resample move can be found in Emerson and Aslett
##' (2023) and Chopin (2002).
##'
##' It is never recommended to use anything other than
##' `IBIS.logreg.opts` to provide the `options` argument. See
##' examples and [IBIS.logreg.opts()] for more information.
##'
##' The prior used for the IBIS procedure will take the form of a multivariate
##' normal, where the parameters can be specified directly by the user. It is
##' however possible to override this default prior distributional form by
##' specifying `prior.override=TRUE` and providing the relevant prior functions
##' in `IBIS.logreg.opts`.
##'
##' @seealso [IBIS.logreg.opts()], [print.IBIS()], [predict.IBIS()], [plot.IBIS()]
##' @references \itemize{
##' \item Emerson, S.R. and Aslett, L.J.M. (2023). Joint cohort and prediction
##' modelling through graphical structure analysis (to be released)
##' \item Chopin, N. (2002). A sequential particle filter method for
##' static models. Biometrika, 89(3), 539-552, \doi{10.1093/biomet/89.3.539}
##' }
##' @examples
##'
##' \donttest{
##' require(graphics)
##' # First we generate a co-variate matrix X and binary response vector y
##' CM <- matrix(rnorm(200),100,2)
##' rv <- sample(0:1,100,replace=TRUE)
##'
##' # Now we can obtain 1000 samples from the posterior from a standard
##' # multivariate normal prior
##' out.1 <- IBIS.logreg(X = CM,y = rv)
##' plot(out.1)
##' out.1$log_Bayesian_evidence
##'
##' # We can specify that the samples be weighted
##' out.1.w <- IBIS.logreg(X = CM,y = rv,
##'                        options = IBIS.logreg.opts(weighted = TRUE))
##' out.1.w$weights
##' plot(out.1.w)
##'
##' # We can also specify different arguments for a specific prior
##' out.2 <- IBIS.logreg(X = CM,y = rv,prior_mean = rep(-3,3),
##'                      prior_var = 0.1*diag(3))
##' samp.df <- data.frame(rbind(out.1$samples,out.2$samples))
##' colnames(samp.df) <- paste0("beta[",c(0:2),"]")
##' GGally::ggpairs(samp.df,
##'                 labeller = "label_parsed",
##'                 ggplot2::aes(color = as.factor(rep(c(1,2),each=1000))),
##'                 upper = list(continuous = GGally::wrap("density")),
##'                 lower = list(continuous = GGally::wrap("points",size=0.5)))
##' out.2$log_Bayesian_evidence
##' out.3 <- IBIS.logreg(X = CM,y = rv,prior_mean = rep(3,3),
##'                      prior_var = 0.1*diag(3))
##' samp.df <- data.frame(rbind(out.1$samples,out.2$samples,out.3$samples))
##' colnames(samp.df) <- paste0("beta[",c(0:2),"]")
##' GGally::ggpairs(samp.df,
##'                 labeller = "label_parsed",
##'                 ggplot2::aes(color = as.factor(rep(c(1,2,3),each=1000))),
##'                 upper = list(continuous = GGally::wrap("density")),
##'                 lower = list(continuous = GGally::wrap("points",size=0.5)))
##' out.3$log_Bayesian_evidence
##'
##' # We can also change the prior, for example a multivariate independent
##' # uniform
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
##' out.4 <- IBIS.logreg(X = CM,y = rv,
##'                      options = IBIS.logreg.opts(prior.override = TRUE,
##'                                                 rprior = rmviu,
##'                                                 dprior = dmviu,a=rep(0,3),
##'                                                 b=rep(1,3)))
##' samp.df <- data.frame(rbind(out.1$samples,out.4$samples))
##' colnames(samp.df) <- paste0("beta[",c(0:2),"]")
##' GGally::ggpairs(samp.df,
##'                 labeller = "label_parsed",
##'                 ggplot2::aes(color = as.factor(rep(c(1,4),each=1000))),
##'                 upper = list(continuous = GGally::wrap("points",size=0.5)),
##'                 lower = list(continuous = GGally::wrap("points",size=0.5)))
##' out.4$log_Bayesian_evidence
##' }
##'

IBIS.logreg <- function(X,y,options = IBIS.logreg.opts(),
                        prior_mean = rep(0,ncol(X)+1),
                        prior_var = diag(ncol(X)+1)){
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
  } else{
    rprior <- mvnfast::rmvn
    dprior <- mvnfast::dmvn
    MoreArgs <- list(mu = prior_mean,sigma = prior_var)
  }
  DM <- as.matrix(cbind(rep(1,length(y)),X))
  IBIS_out <- IBIS.Z(X = DM,y = y,rprior = rprior,N = options$N,
                     prior_pdf = dprior,
                     ess = options$ess,n_move = options$n_move,
                     PriorArgs = MoreArgs,diagnostics=TRUE)$output
  if(options$weighted==FALSE){
    if(length(unique(IBIS_out$weights))!=1){
      ss <- stats::cov.wt(IBIS_out$samples,wt=IBIS_out$weights,method = "ML")
      mu <- ss$center
      Sigma <- ss$cov
      samp <- sample(as.integer(names(IBIS_out$duplication_table)),options$N,prob = IBIS_out$weights[as.integer(names(IBIS_out$duplication_table))]*IBIS_out$duplication_table,replace = TRUE)
      IBIS_out$samples <- IBIS_out$samples[samp,]
      for(j in 1:options$n_move){
        BC <- mvnfast::rmvn(options$N,mu = mu, sigma=Sigma)
        Log_1 <- log(1+exp(DM%*%t(BC)))*(y-1) + log(1+exp(-DM%*%t(BC)))*(-y)
        TotArgs <- list(BC)
        if(length(MoreArgs)!=0){
          TotArgs <- c(TotArgs,MoreArgs)
        }
        Log_1 <- colSums(Log_1) + log(do.call(dprior,TotArgs)) + log(mvnfast::dmvn(IBIS_out$samples,mu=mu,sigma=Sigma))
        Log_2 <- log(1+exp(DM%*%t(IBIS_out$samples)))*(y-1) + log(1+exp(-DM%*%t(IBIS_out$samples)))*(-y)
        TotArgs <- list(IBIS_out$samples)
        if(length(MoreArgs)!=0){
          TotArgs <- c(TotArgs,MoreArgs)
        }
        Log_2 <- colSums(Log_2) + log(do.call(dprior,TotArgs)) + log(mvnfast::dmvn(BC,mu=mu,sigma=Sigma))
        A <- Log_1 - Log_2 - log(stats::runif(length(Log_1)))
        IBIS_out$samples[which(A>0),] <- BC[which(A>0),]
      }
    }
    res <- list(covariate_matrix = data.frame(X),
                response_vector = y,
                samples = IBIS_out$samples,
                log_Bayesian_evidence = IBIS_out$log_Bayesian_evidence,
                diagnostics = IBIS_out$diagnostics)
  } else{
    res <- list(covariate_matrix = data.frame(X),
                response_vector = y,
                samples = IBIS_out$samples,
                weights = IBIS_out$weights,
                log_Bayesian_evidence = IBIS_out$log_Bayesian_evidence,
                diagnostics = IBIS_out$diagnostics)
  }
  class(res) <- 'IBIS'
  res
}

##' Print IBIS
##'
##'
##' @export
##' @keywords print IBIS
##' @name print.IBIS
##' @description Prints summary information from an IBIS object.
##'
##' @param x Object of class `"IBIS"`
##' @param ... Further arguments passed to or from other methods
##' @return No return value, called for side effects
##' @details When running the function [IBIS.logreg()] the printed
##' information will contain information regarding; the number of samples, the
##' mean of those samples and the log Bayesian evidence of the posterior.
##' @seealso [IBIS.logreg()]
##'

print.IBIS <- function(x,...){
  if("weights" %in% names(x)){
    ss <- stats::cov.wt(x$samples,wt=x$weights,method = "ML")
    xx <- ss$center
    names(xx) <- c("(Intercept)",colnames(x$covariate_matrix))
    cat(nrow(x$samples),"posterior samples with weighted mean:\n")
    cat("\n")
    print.default(xx)
  } else{
    xx <- colMeans(x$samples)
    names(xx) <- c("(Intercept)",colnames(x$covariate_matrix))
    cat(nrow(x$samples),"posterior samples with mean:\n")
    cat("\n")
    print.default(xx)
  }
  cat("\n")
  cat("Log Bayesian Evidence:",x$log_Bayesian_evidence)
}

##' Prediction method for IBIS
##'
##'
##' @export
##' @name predict.IBIS
##' @description Predicts the response of new observations using an object of
##' class `"IBIS"`.
##'
##' @keywords predict IBIS
##' @param object Object of class `"IBIS"`
##' @param newX Data frame containing new observations to predict. If not
##' specified the fitted values will be returned instead.
##' @param type Either `"prob"` for a probabilistic output or `"response"` for
##' a hard output of the predicted response
##' @param ... Additional arguments affecting the predictions produced
##' @return Either a matrix of response probabilities for each observation or
##' a vector of predicted responses for each observation.
##' @details Note that this is a Bayesian prediction method as objects with
##' class `"IBIS"` will provide samples from a posterior.
##' @seealso [IBIS.logreg()]
##' @examples
##'
##' \donttest{
##' # First we generate a co-variate matrix X and binary response vector y
##' CM <- data.frame(X1 = rnorm(100),X2 = rnorm(100))
##' rv <- sample(0:1,100,replace=TRUE)
##'
##' # Now we can obtain 1000 samples from the posterior from a standard
##' # multivariate normal prior
##' out <- IBIS.logreg(X = CM,y = rv)
##'
##' # The fitted values of out can be obtained
##' predict(out)
##' predict(out,type = "response")
##'
##' # We can also predict the response for new data
##' CM.2 <- data.frame(X1 = rnorm(10),X2 = rnorm(10))
##' cbind(CM.2,predict(out,newX = CM.2))
##' }
##'


predict.IBIS <- function(object,newX = NULL,type = "prob",...){
  if(type!="prob" & type!="response"){
    stop("type not supported")
  }
  if(is.null(newX)){
    newX = object$covariate_matrix
  }
  newX <- as.matrix(newX,ncol = ncol(object$covariate_matrix))
  DM <- cbind(rep(1,nrow(newX)),newX)
  p1 <- rowMeans(1/(1+exp(-DM%*%t(object$samples))))
  if(type=="prob"){
    res <- data.frame(1-p1,p1)
    colnames(res) <- 0:1
  } else{
    res <- rep(0,length(p1))
    res[which(p1>=0.5)] <- 1
  }
  return(res)
}

##' Plot various outputs of IBIS
##'
##'
##' @export
##' @name plot.IBIS
##' @description Allows visualisation of many aspects of IBIS, including
##' co-variate, posterior and diagnostic plots.
##'
##' @keywords plot IBIS
##' @param x Object of class `"IBIS"`
##' @param type Can be one of; `"samples"` for posterior visualisation,
##' `"fitted"` for co-variate visualisation or `"diagnostics"` for diagnostic
##' plots. See details. Defaults to `"samples"`.
##' @param plot_var Vector specifying which columns (or associated logistic
##' regression coefficients) of the co-variate matrix should be plotted. Does not
##' apply when `type=="diagnostics"`. Defaults to all columns being selected.
##' @param diagnostic_x_axis Only applies if `"type=="diagnostics"`. Either
##' `"full"` (default) for all observations indices to be plotted on the x-axis
##' or `"minimal"` for only some observations indices to be plotted on the
##' x-axis.
##' @param ... Arguments to be passed to methods
##' @return No return value, called for side effects
##' @details If `type=="samples"`, the resulting plot will be a ggpairs plot
##' giving the coefficient densities in the diagonal, points plots of the
##' posterior samples in the lower triangle and contour plots in the upper
##' triangle.
##'
##' If `"type=="fitted"`, the resulting plot will be a ggpairs plot. The
##' diagonal entries will be two density plots, one for training data predicted
##' to have response 0 by the model (red) and one for training data predicted
##' to have response 1 by the model (green). The off-diagonal elements are
##' scatter-plots of the observations, given a label according to their actual
##' response and a colour scale based on their predicted response. If
##' `length(plot_var)==1` then the co-variate variable is plotted against it's
##' index and a density plot is not provided. If `length(plot_var)==1` then the
##' density plot and the scatter-plot are combined. If a predicted class (0 or
##' 1) contains less than two data points the density will not be plotted.
##'
##' If `"type==diagnostics"`, the resulting plot will be a combination of three
##' plots; one tracking the log Bayesian evidence as observations are added to
##' the posterior, one tracking the effective sample size of the particles for
##' each step of the SMC sampler and one tracking the acceptance rate of the
##' Metropolis-Hastings step when a resample-move is triggered. See
##' Emerson and Aslett (2023) and Chopin (2002) for more details. Multiple
##' Metropolis-Hastings steps can be performed when a resample-move step is
##' triggered, and so for the acceptance rate plot observations are suffixed
##' with "." and then the index of current Metropolis-Hastings step. For example
##' the x-axis label for the acceptance rate of the 2nd Metropolis-Hastings step
##' which was triggered by adding observation 1 to the posterior would be
##' labelled "1.2". When the training data for the `"IBIS"` object created is
##' large setting `diagnostic_x_axis=="minimal"` is recommended as it gives a
##' more visually appealing output.
##'
##' @seealso [IBIS.logreg()]
##' @references \itemize{
##' \item Emerson, S.R. and Aslett, L.J.M. (2023). Joint cohort and prediction
##' modelling through graphical structure analysis (to be released)
##' \item Chopin, N. (2002). A sequential particle filter method for
##' static models. Biometrika, 89(3), 539-552, \doi{10.1093/biomet/89.3.539}
##' }
##' @examples
##'
##' \donttest{
##' require(graphics)
##' # First we generate a co-variate matrix X and binary response vector y
##' CM <- matrix(rnorm(200),100,2)
##' rv <- sample(0:1,100,replace=TRUE)
##'
##' # Now we can obtain 1000 samples from the posterior from a standard
##' # multivariate normal prior and plot the results
##' out <- IBIS.logreg(X = CM,y = rv)
##' plot(out,type = "samples")
##' plot(out,type = "fitted")
##' plot(out,type = "diagnostics",diagnostic_x_axis = "minimal")
##'
##' # If we only wanted to view the second co-variate
##' plot(out,type = "samples",plot_var = 2)
##' plot(out,type = "fitted",plot_var = 2)
##' }
##'

plot.IBIS <- function(x,type = "samples",plot_var = NULL,
                      diagnostic_x_axis = "full",...){
  if(diagnostic_x_axis!="full" & diagnostic_x_axis!="minimal"){
    stop("diagnostic_x_axis must be either full or minimal")
  }
  if(is.null(plot_var)){
    plot_var <- 1:ncol(x$covariate_matrix)
  }
  if(any(is.na(match(plot_var,1:ncol(x$covariate_matrix))))){
    stop("cannot subset the co-variate matrix with plot_var provided")
  }
  if(type!="samples" & type!="fitted" & type!="diagnostics"){
    stop("type not supported")
  }
  if(type=="samples"){
    x$samples <- x$samples[,c(1,plot_var+1)]
    samp.df <- data.frame(x$samples)
    colnames(samp.df) <- paste0("beta[",c(0,plot_var),"]")
    if("weights" %in% names(x)){
      overall_plot <- GGally::ggpairs(samp.df,labeller = "label_parsed",
                                      upper = list(continuous = GGally::wrap("density",
                                                                             colour = "black")),
                                      lower = list(continuous = GGally::wrap("points",
                                                                             size=0.5,
                                                                             alpha = x$weights)))
    } else{
      overall_plot <- GGally::ggpairs(samp.df,labeller = "label_parsed",
                                      upper = list(continuous = GGally::wrap("density",
                                                                             colour = "black")),
                                      lower = list(continuous = GGally::wrap("points",
                                                                             size=0.5)))
    }
  }
  if(type=="fitted"){
    X_prob <- predict.IBIS(object = x)[,2]
    x$covariate_matrix <- x$covariate_matrix[,plot_var,drop=FALSE]
    if(ncol(x$covariate_matrix)==1){
      overall_plot <- ggplot2::ggplot(data.frame(x$covariate_matrix,
                                                 y.axis = rep(0,nrow(x$covariate_matrix)),
                                                 y = as.character(x$response_vector),
                                                 Fitted.Probabilities = X_prob,
                                                 Hard.Assignment = X_prob>=0.5)) +
        ggplot2::geom_density(show.legend = FALSE,
                              ggplot2::aes_string(colnames(x$covariate_matrix)[1],
                                                  colour = "Hard.Assignment",
                                                  fill = "Hard.Assignment"),alpha=0.25) +
        ggplot2::scale_colour_manual(values = c("red","green")) +
        ggplot2::scale_fill_manual(values = c("red","green")) +
        ggnewscale::new_scale_color() +
        ggplot2::geom_point(alpha=0,
                            ggplot2::aes_string(colnames(x$covariate_matrix)[1],
                                                "y.axis",
                                                colour = "Fitted.Probabilities")) +
        ggplot2::geom_text(alpha = 1,show.legend = FALSE,size = 3,
                           ggplot2::aes_string(colnames(x$covariate_matrix)[1],
                                               "y.axis",
                                               label = "y",
                                               colour = "Fitted.Probabilities")) +
        ggplot2::scale_color_gradient(low = "red",high = "green",name = "Fitted Probabilities",limits = c(0,1)) +
        ggplot2::labs(y = "density") +
        ggplot2::theme(legend.position = "bottom")
    } else{
      legend_plot <- ggplot2::ggplot(data.frame(x$covariate_matrix[,c(1,2)],
                                                Fitted.Probabilities = X_prob),
                                     ggplot2::aes_string(colnames(x$covariate_matrix)[1],
                                                         colnames(x$covariate_matrix)[2],
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
      overall_plot <- GGally::ggpairs(data.frame(x$covariate_matrix,
                                                 y = as.character(x$response_vector),
                                                 Fitted.Probabilities = X_prob,
                                                 Hard.Assignment = X_prob>=0.5),
                                      ggplot2::aes_string(alpha = 0.5),
                                      columns = 1:ncol(x$covariate_matrix),
                                      diag = list(continuous = diag_plot),
                                      upper = list(continuous = label_plot_fitted),
                                      lower = list(continuous = label_plot_fitted),
                                      legend = GGally::grab_legend(legend_plot)) +
        ggplot2::theme(legend.position = "bottom")
    }
  }
  if(type=="diagnostics"){
    obs_lev <- x$diagnostics$log_Bayesian_evidence_tracker$Observations
    obs_lev_tics <- rep("",length(obs_lev))
    obs_lev_tics[c(1,round(length(obs_lev)*(1:4)/4))] <- obs_lev[c(1,round(length(obs_lev)*(1:4)/4))]
    x$diagnostics$log_Bayesian_evidence_tracker$Observations <- factor(x$diagnostics$log_Bayesian_evidence_tracker$Observations,levels = obs_lev)
    plot_1 <- ggplot2::ggplot(x$diagnostics$log_Bayesian_evidence_tracker,
                              ggplot2::aes_string(x = "Observations",
                                  y = "Log_Bayesian_Evidence",group=1)) +
      ggplot2::geom_line() + ggplot2::labs(x = "Observations", y = "Log Bayesian Evidence")
    if(diagnostic_x_axis=="minimal"){
      plot_1 <- plot_1 + ggplot2::scale_x_discrete(labels = obs_lev_tics)
    }
    x$diagnostics$effective_sample_size_tracker$Observations <- factor(x$diagnostics$effective_sample_size_tracker$Observations,levels = obs_lev)
    plot_2 <- ggplot2::ggplot(x$diagnostics$effective_sample_size_tracker,
                              ggplot2::aes_string(x = "Observations",
                                  y = "ESS",group=1)) +
      ggplot2::geom_line() + ggplot2::labs(x = "Observations") +
      ggplot2::geom_hline(yintercept = x$diagnostics$ESS_threshold,linetype = 2) +
      ggplot2::scale_y_continuous(breaks = unique(sort(c(round(nrow(x$samples)*(2:5)/5),x$diagnostics$ESS_threshold))))
    if(diagnostic_x_axis=="minimal"){
      plot_2 <- plot_2 + ggplot2::scale_x_discrete(labels = obs_lev_tics)
    }
    if(nrow(x$diagnostics$acceptance_rate_tracker)==0){
      overall_plot <- ggpubr::ggarrange(plot_1,plot_2)
    } else{
      obs_lev_tics2 <- rep("",length(x$diagnostics$acceptance_rate_tracker$Observations))
      obs_lev_tics2[c(1,round(length(obs_lev_tics2)*(1:4)/4))] <- x$diagnostics$acceptance_rate_tracker$Observations[c(1,round(length(obs_lev_tics2)*(1:4)/4))]
      plot_3 <- ggplot2::ggplot(x$diagnostics$acceptance_rate_tracker,
                                ggplot2::aes_string(x = "Observations",
                                    y = "Acceptance_Rate",group=1)) +
        ggplot2::geom_line() + ggplot2::labs(y = "Acceptance Rate")
      if(diagnostic_x_axis=="minimal"){
        plot_3 <- plot_3 + ggplot2::scale_x_discrete(labels = obs_lev_tics2)
      }
      overall_plot <- ggpubr::ggarrange(ggpubr::ggarrange(plot_1,plot_2,ncol=2),plot_3,nrow=2)
    }
  }
  suppressWarnings(print(overall_plot))
}

##' Additional argument generator for [IBIS.logreg()]
##'
##'
##' @export
##' @name IBIS.logreg.opts
##' @description This function is used to specify additional arguments to
##' `IBIS.logreg`.
##'
##' @keywords IBIS options control
##' @param N Number of particles for the SMC sampler. Defaults to 1000.
##' @param ess Effective Sample Size Threshold: If the effective sample size of
##' the particles falls below this value then a resample move step is
##' triggered. Defaults to `N/2`.
##' @param n_move Number of Metropolis-Hastings steps to apply each time a
##' resample move step is triggered. Defaults to 1.
##' @param weighted Should the outputted samples be weighted? Defaults to
##' `FALSE`.
##' @param prior.override Are you overriding the default multivariate normal
##' form of the prior? Defaults to `FALSE`.
##' @param rprior Function which produces samples from your prior if the default
##' prior form is to be overridden. If using the default prior form this does
##' not need to be specified.
##' @param dprior Function which produces your specified priors density for
##' inputted samples if the default prior form is to be overridden. If using the
##' default prior form this does not need to be specified.
##' @param ... Additional arguments required for complete specification of the
##' two prior functions given, if the default prior form is to be overridden.
##' @return A list consisting of:
##'
##' \describe{
##' \item{`N`}{Number of particles for the SMC sampler}
##' \item{`ess`}{Effective Sample Size Threshold}
##' \item{`n_move`}{Number of Metropolis-Hastings steps}
##' \item{`rprior`}{Function which produces samples from your prior. `NULL` if
##' `prior.override==FALSE`.}
##' \item{`dprior`}{Function which produces your specified priors density for
##' inputted samples. `NULL` if `prior.override==FALSE`.}
##' \item{`prior.override`}{Logical value indicating if the prior has been
##' overridden or not.}
##' \item{`weighted`}{Logical value indicating if the outputted particles of
##' `IBIS.logreg` should be weighted or not.}
##' \item{`MoreArgs`}{A list of the additional arguments required for `rprior`
##' and `dprior`. `NULL` if `prior.override==FALSE`.}
##' }
##'
##' @details This function should only be used to provide additional control
##' arguments to `IBIS.logreg`.
##'
##' Specifying `rprior` and `dprior` will not override the default prior form
##' unless `prior.override=TRUE`. If a multivariate normal form is required then
##' the arguments for this prior should be specified in `IBIS.logreg`.
##' @seealso [IBIS.logreg()]
##'
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
##' IBIS.logreg.opts(prior.override = TRUE,rprior = rmviu,
##'                  dprior = dmviu,a=rep(0,3),b=rep(1,3))
##'
##'

IBIS.logreg.opts <- function(N=1000,ess = N/2,n_move = 1,weighted = FALSE,
                             prior.override = FALSE,rprior = NULL,
                             dprior = NULL,...){
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
  if(ess > N | ess < 0){
    stop("Effective sample size must be between 0 and N")
  }
  if(n_move < 1){
    stop("n_move must be an integer greater or equal than 1")
  }
  if(!is.logical(weighted)){
    stop("weighted must be logical")
  }
  return(list(N = N,ess = ess,n_move = n_move,rprior = rprior,dprior = dprior,
              prior.override = prior.override,weighted = weighted,
              MoreArgs = MoreArgs))
}


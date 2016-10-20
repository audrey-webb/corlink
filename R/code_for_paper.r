

# alldata is a table with match patterns for all rows being linked - the last pattern being the counts of each pattern.  0 means fields observed for both records, no match; 1 fields observed on both records, match; 2 fields missing
# alg takes 3 values "i" for original Feligi Sunter EM algorithm, "m" for correlation within agreement patterns for 'matches' only, "b" for general model and "a" for all models.
#' Function to impute missing agreement patterns and then to link data
#'
#' @param d Matrix of agreement patterns with final column counting the number of times that pattern was observed.  See Details
#' @param initial_m starting probabilities for per-field agreement in record pairs, both records being generated from the same individual.  Defaults to NULL
#' @param initial_u starting probabilities for per-field agreement in record pairs, with the two records being generated from differing individuals  Defaults to NULL
#' @param p_init starting probability that both records for a randomly selected record pair is associated with the same individual
#' @param fixed_col vector indicating columns that are not to be updated in initial EM algorithm.  Useful if good prior estimates of the mis-match probabilities.  See details
#' @param alg character; see Details
#' @keywords EM algorithm, probabilistic linkage, Feligi/Sunter, latent class, correlation
#' @details \code{d} is a numeric matrix with N rows corresponding to N record pairs, and L+1 columns the first L of which show the field agreement patterns observed over the record pairs, and the last column the total number of times that pattern was observed in the database.  The code 0 is used for a field that differs for two record, 1 for a field that agrees, and 2 for a missing field. \code{fixed_col} indicates the components of the \code{u} vector (per field probabilities of agreement for 2 records from differing individuals) that are not to be updated when applying the EM algorithm to estimate components of the Feligi Sunter model.  \code{alg} has four possible values.  The default \code{'m'} fits a log-linear model for the agreement counts only within the record pairs that corresponds to the same individual, \code{'b'} fits differing log-linear models for the 2 clusters, \code{'i'} corresponds to the original Feligi Sunter algorithm, with probabilities estimated via the EM algorithm, \code{'a'} fits all the previously listed models
#' @return A list, the first component is a matrix -  the posterior probabilities of being a true match is the last column, the second component are the fitted models used to generate the predicted probabilities
#' @importFrom stats qnorm rbinom rnorm
#' @importFrom utils flush.console
#' @export
#' @examples
#'
#' # Simulate data
#' m_probs <- rep(0.8,6)
#' u_probs <- rep(0.2,6)
#' means_match <- -1*qnorm(1-m_probs)
#' means_mismatch <- -1*qnorm(1-u_probs)
#' missingprobs <- rep(.2,6)
#' thedata <- do_sim(cor_match=0.2,cor_mismatch=0,nsample=10^4,pi_match=.5,
#' m_probs=rep(0.8,5),u_probs=rep(0.2,5),missingprobs=rep(0.4,5))
#' colnames(thedata) <- c(paste("V",1:5,sep="_"),"count")
#' output <- linkd(thedata)
#' output$fitted_probs

linkd <- function(d, initial_m=NULL, initial_u=NULL, p_init=0.5,fixed_col=NULL,alg="m"){
  colnames(d)[ncol(d)] <- "totalcount"
  L <- ncol(d) - 1
  if(is.null(initial_m)) initial_m <- rep(0.8,L)
  if(is.null(initial_u)) initial_u <- rep(0.2,L)
  if(is.null(fixed_col)) fixed_col <- c()
  themat <- create_01mat(L)
  out <- imputemissing(d,themat)
  out <- cbind(themat,counts=out$counts)
  out <- as.data.frame(out)
  L <- ncol(d)-1
  allterms <- c()
  for (lower in 1:(L-2)) {
    for(upper in (lower+1):L)
      allterms <- c(allterms,paste("var",lower, "*","var",upper, sep = ""))
  }
  allterms <- c(allterms,paste("var",(L-1),"*","var",L , sep = ""))
  if(alg == "i" | alg == "a"){
  results_independence <- EM_match_independence_v2(out,m_1=initial_m,u_1=initial_u,m_2=rep(0,L),u_2=rep(0,L),p_init=p_init,tol=10^-5, fixedcol=fixed_col)
  probs_independence <- reassign_probs(d[,1:(ncol(d)-1)], out, results_independence$probs)
  }
  if(alg == "b" | alg == "a"){
  results_loglinear <- EM_match_modelsearch(out,m_1=initial_m,u_1=initial_u,m_2=rep(0,L),u_2=rep(0,L),p_init=p_init,tol=10^-5, fixedcol=fixed_col,allterms=allterms)
  probs_loglinear <- reassign_probs(d[,1:(ncol(d)-1)], out,results_loglinear$probs)
  }
  if(alg == "m" | alg == "a"){
   results_loglinear_iu <- EM_match_modelsearch_iu(out,m_1=initial_m,u_1=initial_u,m_2=rep(0,L),u_2=rep(0,L),p_init=p_init,tol=10^-5, fixedcol=fixed_col,allterms=allterms)
   probs_loglinear_iu <- reassign_probs(d[,1:(ncol(d)-1)], out,results_loglinear_iu$probs)
  }
  if(alg == "i") return(list(fitted_probs=cbind(d,"fitted_prob_match"=probs_independence),fitted_models=results_independence))
  if(alg == "b") return(list(fitted_probs=cbind(d,"fitted_prob_match"=probs_loglinear),fitted_models=results_loglinear))
  if(alg == "m") return(list(fitted_probs=cbind(d,"fitted_prob_match"=probs_loglinear_iu),fitted_models=results_loglinear_iu))
  if(alg == "a") return(list(fitted_probs=cbind(d,"independence"=probs_independence,"log_linear"=probs_loglinear,"log_linear_iu"=probs_loglinear_iu),model_loglinear_iu=results_loglinear_iu,model_loglinear=results_loglinear,model_independence=results_independence))

}



#' corlink: Record linkage, Incorporating Imputation for Missing Agreement Patterns, and Modeling correlation patterns between fields 
#'
#' A matrix of agreement patterns and counts for record pairs is the input for the procedure.  An EM algorithm is used to impute plausible values for missing record pairs.  A second EM algorithm, incorporating possible correlations between per-field agreement, is used to estimate posterior probabilites that each pair is a true match - i.e. constitutes the same individual.
#' 
#' @section corlink functions:
#' linkd 
#'
#' @docType package
#' @name corlink
NULL
#> NULL
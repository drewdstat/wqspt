#' Full wrapper WQS permutation test 
#' 
#' \code{wqs_full_perm} is a full wrapper function that is a full implementation 
#' of the Weighted Quantile Sum (WQS) regression method followed by the 
#' permutation test to determine the significance of the WQS coefficient. 
#' 
#' @param formula An object of class formula. The wqs term must be included in 
#' the formula (e.g., y ~ wqs + ...).
#' @param data The \code{data.frame} to be used in the WQS regression run. This
#' can be of class \code{data.frame} or it can be a tibble from the tidyverse.
#' @param mix_name A vector with the mixture column names. 
#' @param q An integer to indicate the number of quantiles to split the mixture 
#' variables. 
#' @param b_main The number of bootstraps for the main WQS regression run. 
#' @param b_perm The number of bootstraps for the iterated permutation test 
#' WQS regression runs and the reference WQS regression run (only for linear
#' WQS regression and only when b_mean != b_perm). 
#' @param b1_pos A logical value that indicates whether beta values should be 
#' positive or negative.
#' @param rs A logical value indicating whether random subset implementation 
#' should be performed. 
#' @param niter Number of permutation test iterations. 
#' @param seed An integer to fix the seed. This will only impact the the initial 
#' WQS regression run and not the permutation test iterations. The default 
#' setting is NULL, which means no seed is used for the initial WQS regression. 
#' The seed will be saved in the "gwqs_main" object as "gwqs_main$seed".
#' @param plan_strategy Evaluation strategy for the plan function. You can choose 
#' among "sequential", "transparent", "multisession", "multicore", "multiprocess", 
#' "cluster" and "remote." See future::plan documentation for full details. 
#' @param b_constr Logical value that determines whether to apply positive or 
#' negative constraints in the optimization function for the weight optimization.
#' Note that this won't guarantee that the iterated b1 values in the 
#' weight optimization are only positive (if \code{b1_pos = TRUE}) or only 
#' negative (if \code{b1_pos = FALSE}) as seen in the \code{bres} matrix output 
#' by the \code{gwqs} models (i.e., column \code{bres$b1}), but it does 
#' substantially increase the probability that those b1 values will be 
#' constrained to be either positive or negative. This defaults to \code{FALSE}.
#' @param family A description of the error distribution and link function to be 
#' used in the model. This can be a character string naming a family function 
#' (e.g., \code{"binomial"}) or a family object (e.g., 
#' \code{binomial(link="logit")}). Currently validated families include 
#' \code{gaussian(link="identity")} for linear regression, \code{binomial()} 
#' with any accepted link function (e.g., \code{"logit"} or \code{"probit"}), 
#' \code{poisson(link = "log")}, \code{quasipoisson(link = "log")}, or 
#' \code{"negbin"} for negative binomial. The \code{"multinomial"} family is 
#' not yet supported.
#' @param stop_if_nonsig if TRUE, the function will not proceed with the 
#' permutation test if the main WQS regression run produces nonsignificant 
#' p-value.
#' @param stop_thresh numeric p-value threshold required in order to proceed 
#' with the permutation test, if \code{stop_if_nonsig = TRUE}.
#' @param ... Other parameters to put into the gwqs function call.
#'
#' @return \code{wqs_full_perm} returns an object of class \code{wqs_perm}, 
#' which contains three sublists: 
#' 
#' \item{perm_test}{List containing: 
#' \itemize{
#'   \item \code{pval}: permutation test p-value 
#'   \item (linear regression only) \code{testbeta1}: reference WQS regression 
#'   coefficient beta1 value 
#'   \item (linear regression only) \code{betas}: Vector of beta values from 
#'   each permutation test run 
#'   \item (logistic regression only) \code{testpval}: test reference p-value 
#'   \item (logistic regression only) \code{permpvals}: p-values from the null models
#' }
#' }
#' \item{gwqs_main}{Main gWQS object (same as model input). This will now 
#' include an additional object "seed" that returns the seed used for this main 
#' WQS regression.}
#' \item{gwqs_perm}{Permutation test reference gWQS object (NULL if model 
#' \code{family != "gaussian"} or if same number of bootstraps are used in 
#' permutation test WQS regression runs as in the main run).}
#' @import gWQS
#' @export wqs_full_perm
#'
#' @examples
#' library(gWQS)
#'
#' # mixture names
#' PCBs <- names(wqs_data)[1:5] 
#' # Only using 1st 5 of the original 34 exposures for this quick example
#' 
#' # quick example with only 4 bootstraps each WQS regression iteration, and 
#' # only 3 iterations
#' 
#' perm_test_res <- wqs_full_perm(formula = yLBX ~ wqs, data = wqs_data, 
#'                                 mix_name = PCBs, q = 10, b_main = 4, 
#'                                 b_perm = 4, b1_pos = TRUE, b_constr = FALSE, 
#'                                 niter = 3, seed = 16, 
#'                                 plan_strategy = "multicore", 
#'                                 stop_if_nonsig = FALSE)
#' 
#' # Note: The default values of b_main = 1000, b_perm = 200, and niter = 200 
#' # are the recommended parameter values. This example has a lower b_main, 
#' # b_perm, and niter in order to serve as a shorter test run. 
#'  
wqs_full_perm <- function(formula, data, mix_name, q = 10, b_main = 1000, 
                          b_perm = 200, b1_pos = TRUE, b_constr = FALSE, 
                          rs = FALSE, niter = 200, seed = NULL, 
                          family = "gaussian", plan_strategy = "multicore",
                          stop_if_nonsig = FALSE, stop_thresh = 0.05, ...){
  
  if (is.character(family)) {
    if (family=="multinomial"){
      stop("This simulation function doesn't yet accomodate 
           multinomial WQS regression.")
    }
  }
  
  # run main WQS regression
  gwqs_res_main <- gWQS::gwqs(formula = formula, data = data, mix_name = mix_name, 
                              q = q, b = b_main, b1_pos = b1_pos, 
                              b_constr = b_constr, rs = rs, seed = seed, 
                              validation = 0, family = family, 
                              plan_strategy = plan_strategy, ...) 
  
  gwqs_res_main$seed<-seed
  naive_p <- summary(gwqs_res_main)$coefficients["wqs", 4]

  if (stop_if_nonsig == TRUE & naive_p > stop_thresh){
    message(sprintf("The main WQS regression run did not give a significant 
                    result (p = %s)", 
                    naive_p))
    
    results <- list(gwqs_main = gwqs_res_main, 
                    family = gwqs_res_main$family$family,
                    gwqs_perm = NULL, 
                    perm_test = NULL)
  } else {
    # run permutation test (using wqs_perm function) 
    results <- wqs_pt(gwqs_res_main, niter = niter, boots = b_perm, 
                        b1_pos = b1_pos, b_constr = b_constr, rs = rs, 
                        plan_strategy = plan_strategy, seed = seed)
  }
  
  class(results) <- "wqs_pt"
  
  results
}

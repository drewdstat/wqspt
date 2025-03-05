#' Full wrapper WQS permutation test 
#' 
#' \code{wqs_full_perm} is a full wrapper function that is a full implementation 
#' of the Weighted Quantile Sum (WQS) regression method followed by the 
#' permutation test to determine the significance of the WQS coefficient. 
#' 
#' @param formula An object of class formula. The wqs term must be included in 
#' the formula (e.g., \code{y ~ wqs + ...}).
#' @param data The \code{data.frame} to be used in the WQS regression run. This
#' can be of class \code{data.frame} or it can be a tibble from the tidyverse.
#' @param mix_name A vector with the mixture column names. 
#' @param q An integer to indicate the number of quantiles to split the mixture 
#' variables. 
#' @param b_main The number of bootstraps for the main WQS regression run. 
#' @param b_perm The number of bootstraps for the iterated permutation test 
#' WQS regression runs and the reference WQS regression run (only for linear
#' WQS regression and only when \code{b_mean != b_perm}). 
#' @param b1_pos A logical value that indicates whether beta values should be 
#' positive or negative.
#' @param rs A logical value indicating whether random subset implementation 
#' should be performed. 
#' @param niter Number of permutation test iterations. 
#' @param seed An integer to fix the seed. This will only impact the the initial 
#' WQS regression run and not the permutation test iterations. The default 
#' setting is NULL, which means no seed is used for the initial WQS regression. 
#' The seed will be saved in the \code{"gwqs_main"} object as 
#' \code{"gwqs_main$seed"}.
#' @param plan_strategy Evaluation strategy for the plan function. You can choose 
#' among \code{"sequential"}, \code{"multisession"}, \code{"multicore"}, and  
#' \code{"cluster"}. This defaults to \code{"multicore"}. See the 
#' \code{future::plan} documentation for full details. 
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
#' @param nworkers (optional) If the \code{plan_strategy} is not 
#' \code{"sequential"}, this argument defines the number of parallel processes 
#' to use, which can be critical when using a high-performance computing (HPC) 
#' cluster. This should be an integer value. The default behavior for 
#' \code{gWQS::gwqs} is to use all detected cores on a machine, but for many 
#' HPC use scenarios, this would call in cores that have not been allotted by 
#' the HPC scheduler, resulting in the submitted job being halted. For example, 
#' if one has requested 14 cores on a 28-core HPC queue, one would want to set 
#' \code{nworkers = 14}. If \code{nworkers} was greater than 14 in that case, 
#' the HPC job would be terminated. This argument defaults to \code{NULL}, in 
#' which case \code{length(future::availableWorkers())} will be used to 
#' determine the number of parallel processes to use. 
#' @param ... (optional) Additional arguments to pass to the \code{gWQS::gwqs} 
#' function.
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
#' @import gWQS future
#' @export wqs_full_perm
#'
#' @examples
#' library(gWQS)
#'
#' # mixture names
#' PCBs <- names(wqs_data)[1:5] 
#' # Only using 1st 5 of the original 34 exposures for this quick example
#' 
#' # quick example with only 3 bootstraps each WQS regression iteration, and 
#' # only 2 iterations
#' 
#' perm_test_res <- wqs_full_perm(formula = yLBX ~ wqs, data = wqs_data, 
#'                                 mix_name = PCBs, q = 10, b_main = 3, 
#'                                 b_perm = 3, b1_pos = TRUE, b_constr = FALSE, 
#'                                 niter = 2, seed = 16, 
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
                          stop_if_nonsig = FALSE, stop_thresh = 0.05, 
                          nworkers = NULL, ...){
  
  if (is.character(family)) {
    if (family=="multinomial"){
      stop("This simulation function doesn't yet accomodate 
           multinomial WQS regression.")
    }
  }
  
  if(is.null(nworkers)) nworkers = availableWorkers()
  if(length(nworkers) > 1) nworkers = length(nworkers)
  
  # run main WQS regression
  gwqs_res_main <- gwqs_hpc(formula = formula, data = data, mix_name = mix_name,
                            q = q, b = b_main, b1_pos = b1_pos,
                            b_constr = b_constr, rs = rs, seed = seed,
                            validation = 0, family = family,
                            plan_strategy = plan_strategy, 
                            n_workers = nworkers, ...) 
  
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
                      plan_strategy = plan_strategy, seed = seed, 
                      nworkers = nworkers, ...)
  }
  
  class(results) <- "wqs_pt"
  
  results
}

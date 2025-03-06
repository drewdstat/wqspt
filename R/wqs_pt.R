#' WQS permutation test
#' 
#' \code{wqs_pt} takes a \code{gwqs} object as an input and runs the permutation 
#' test (Day et al. 2022) to obtain an estimate for the p-value significance for 
#' the WQS coefficient.  
#' 
#' To use \code{wqs_pt}, we first need to run an initial WQS regression run while 
#' setting \code{validation = 0}. We will use this \code{gwqs} object as the 
#' model argument for the \code{wqs_pt} function. Note that permutation test 
#' has so far only been validated for linear WQS regression (i.e., 
#' \code{family = "gaussian"}) or logistic WQS regression (i.e., 
#' \code{family = binomial(link = "logit")}), though the
#' permutation test algorithm should also work for all WQS GLMs. Therefore,
#' this function accepts \code{gwqs} objects made with the following families: 
#' \code{"gaussian"} or \code{gaussian(link = "identity")}, \code{"binomial"} or 
#' \code{binomial()} with any accepted link function (e.g., \code{"logit"} or 
#' \code{"probit"}), \code{"poisson"} or \code{poisson(link="log")}, 
#' \code{"negbin"} for negative binomial, and \code{"quasipoisson"} or
#' \code{quasipoisson(link="log")}. This function cannot currently accommodate 
#' \code{gwqs} objects made with the \code{"multinomial"} family, and it is not 
#' currently able to accommodate stratified weights or WQS interaction terms 
#' (e.g., \code{y ~ wqs * sex}).
#' 
#' The argument \code{boots} is the number of bootstraps for the WQS regression 
#' run in each permutation test iteration. Note that we may elect a bootstrap 
#' count \code{boots} lower than that specified in the model object for the 
#' sake of efficiency. If \code{boots} is not specified, then we will use the 
#' same bootstrap count in the permutation test WQS regression runs as that 
#' specified in the model argument.
#'
#' The arguments \code{b1_pos} and \code{rs} should be consistent with the 
#' inputs chosen in the model object. The seed should ideally be consistent 
#' with the seed set in the model object for consistency, though this is not 
#' required.
#'
#' @param model A \code{gwqs} object as generated from the \code{gWQS} package.  
#' @param niter Number of permutation test iterations. 
#' @param boots Number of bootstrap samples for each permutation test WQS 
#' regression iteration. If \code{boots} is not specified, then we will use the 
#' same bootstrap count for each permutation test WQS regression iteration as  
#' was specified in the main WQS regression run.
#' @param b1_pos A logical value that indicates whether beta values should be 
#' positive or negative.
#' @param b_constr Logical value that determines whether to apply positive or 
#' negative constraints in the optimization function for the weight optimization. 
#' Note that this won't guarantee that the iterated b1 values in the 
#' weight optimization are only positive (if \code{b1_pos = TRUE}) or only 
#' negative (if \code{b1_pos = FALSE}) as seen in the \code{bres} matrix output 
#' by the \code{gwqs} models (i.e., column \code{bres$b1}), but it does 
#' substantially increase the probability that those b1 values will be 
#' constrained to be either positive or negative. This defaults to \code{FALSE}.
#' @param rs A logical value indicating whether random subset implementation 
#' should be performed. 
#' @param plan_strategy Evaluation strategy for the \code{plan} function. You 
#' can choose among \code{"sequential"}, \code{"transparent"}, 
#' \code{"multisession"}, \code{"multicore"}, and \code{"cluster"}. This 
#' defaults to \code{"multicore"}. See the \code{future::plan} documentation 
#' for full details. 
#' @param seed (optional) Random seed for the permutation test WQS reference run. 
#' This should be the same random seed as used for the main WQS regression run. 
#' This seed will be saved in the \code{"gwqs_perm"} object as 
#' \code{gwqs_perm$seed}. This defaults to \code{NULL}.
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
#' @return \code{wqs_pt} returns an object of class \code{"wqs_pt"}, which 
#' contains: 
#' 
#' \item{perm_test}{List containing: (1) \code{pval}: permutation test p-value, 
#' (2) (linear WQS regression only) \code{testbeta1}: reference WQS coefficient 
#' beta1 value, 
#' (3) (linear WQS regression only) \code{betas}: Vector of beta values from 
#' each permutation test run, (4) (WQS GLM only) \code{testpval}: test reference 
#' p-value, (5) (WQS GLM only) \code{permpvals}: p-values from the null models.}
#' \item{gwqs_main}{Main gWQS object (same as model input).}
#' \item{gwqs_perm}{Permutation test reference gWQS object (NULL if model 
#' \code{family != "gaussian"} or if same number of bootstraps are used in 
#' permutation test WQS regression runs as in the main run).}
#' @import gWQS ggplot2 viridis cowplot stats methods future
#' @export wqs_pt
#' 
#' @examples
#' library(gWQS)
#' 
#' # mixture names
#' PCBs <- names(wqs_data)[1:5] 
#'  # Only using 1st 5 of the original 34 exposures for this quick example
#' 
#' # create reference wqs object with 4 bootstraps
#' wqs_main <- gwqs(yLBX ~ wqs, mix_name = PCBs, data = wqs_data, q = 10, 
#'                  validation = 0, b = 3, b1_pos = TRUE, b_constr = FALSE,
#'                  plan_strategy = "multicore", family = "gaussian", seed = 16)
#' # Note: We recommend niter = 1000 for the main WQS regression. This example
#' # has a lower number of bootstraps to serve as a shorter test run.
#' 
#' # run the permutation test
#' perm_test_res <- wqs_pt(wqs_main, niter = 2, b1_pos = TRUE)
#' 
#' # Note: The default value of niter = 200 is the recommended parameter value. 
#' # This example has a lower niter in order to serve as a shorter test run. 
#' 
#' @references
#' 
#' Day, D. B., Sathyanarayana, S., LeWinn, K. Z., Karr, C. J., Mason, W. A., & 
#' Szpiro, A. A. (2022). A permutation test-based approach to strengthening 
#' inference on the effects of environmental mixtures: comparison between single 
#' index analytic methods. Environmental Health Perspectives, 130(8).
#' 
#' Day, D. B., Collett, B. R., Barrett, E. S., Bush, N. R., Swan, S. H., Nguyen, 
#' R. H., ... & Sathyanarayana, S. (2021). Phthalate mixtures in pregnancy, 
#' autistic traits, and adverse childhood behavioral outcomes. Environment 
#' International, 147, 106330.
#' 
#' Loftus, C. T., Bush, N. R., Day, D. B., Ni, Y., Tylavsky, F. A., Karr, C. J., 
#' ... & LeWinn, K. Z. (2021). Exposure to prenatal phthalate mixtures and 
#' neurodevelopment in the Conditions Affecting Neurocognitive Development and 
#' Learning in Early childhood (CANDLE) study. Environment International, 150, 
#' 106409.
#' 
wqs_pt <- function(model, niter = 200, boots = NULL, b1_pos = TRUE, 
                     b_constr = FALSE, rs = FALSE, plan_strategy = "multicore", 
                     seed = NULL, nworkers = NULL, ...) {
  
  pbapply::pboptions(type="timer")
  
  if (is(model, "gwqs")) {
    if (model$family$family == "multinomial"){
      stop("The permutation test is not currently set up to accomodate the 
           multinomial WQS regressions.")
    }
  } else stop("'model' must be of class 'gwqs' (see gWQS package).")
  
  mm <- model$fit
  formchar <- as.character(formula(mm))
  
  if(is.null(nworkers)) nworkers = availableWorkers()
  if(length(nworkers) > 1) nworkers = length(nworkers)

  if (length(formchar) == 1) {
    tempchar <- rep(NA, 3)
    tempchar[1] <- "~"
    tempchar[2] <- gsub("\\ ~.*", "", formchar)
    tempchar[3] <- gsub(".*\\~ ", "", formchar)
    formchar <- tempchar
    rm(tempchar)
  }
  
  if (!is.null(model$stratified) | grepl("wqs:", formchar[3], fixed = TRUE))
  {
    stop("This permutation test is not yet set up to accomodate stratified 
         weights or WQS interaction terms.")
  }  
  
  cl = match.call()
  yname <- formchar[2]
  mix_name <- names(model$bres)[names(model$bres) %in% 
                                  model$final_weights$mix_name]
  
  if (!is.null(model$qi)) {
    nq <- max(sapply(model$qi, length)) - 1
  } else {
    # this is for cases when there is no quantile transformation or it's 
    # already been done in the data frame
    nq <- NULL
  }
  
  if (is.null(boots)){
    boots <- length(model$bindex)
  }
  
  if (model$family$family == "gaussian"){
    Data <- model$data[, -which(names(model$data) %in% c("wqs", "wghts"))]
    
    # reference WQS regression run 
    if(boots == length(model$bindex)){
      perm_ref_wqs <- model
      ref_beta1 <- mm$coef[2]
    } else {
      perm_ref_wqs <- gwqs_hpc(formula = formula(mm), data = Data, 
                               mix_name = mix_name, q = nq, b = boots, rs = rs, 
                               validation = 0, plan_strategy = plan_strategy, 
                               b1_pos = b1_pos, b_constr = b_constr, 
                               seed = seed, n_workers = nworkers, ...)
      
      ref_beta1 <- perm_ref_wqs$fit$coef[2]
    }
    
    
    if (length(mm$coef) > 2) {
      # This is the permutation test algorithm when there are multiple independent 
      # variables in the model
      lm_form <- formula(paste0(formchar[2], formchar[1], 
                                gsub("wqs + ", "", formchar[3], fixed = TRUE)))
      fit.partial <- lm(lm_form, data = Data)
      partial.yhat <- predict(fit.partial)
      partial.resid <- resid(fit.partial)
      reorgmat <- matrix(NA, dim(Data)[1], niter)
      reorgmat <- apply(reorgmat, 2, function(x) 
        partial.yhat + sample(partial.resid, replace = F))
    } else {
      # This is the permutation test algorithm when there is only one independent 
      # variable in the model
      if(any(grepl("tbl", class(Data)))){
        Data <- as.data.frame(Data)
      }
      reorgmat <- matrix(NA, dim(Data)[1], niter)
      reorgmat <- apply(reorgmat, 2, function(x) sample(Data[, yname]))
    }
    
    # This function is the main iterated worker function for wqs_pt
    getbetas <- function(x) {
      
      newDat <- Data
      newDat[, yname] <- x
      names(newDat) <- c(names(Data))
      
      if (length(mm$coef) > 2) {
        form1 <- formula(paste0(formchar[2], formchar[1], formchar[3]))
      } else {
        form1 <- formula(paste0(formchar[2], formchar[1], "wqs"))
      }
      
      gwqs1 <- tryCatch({
        suppressWarnings(gwqs_hpc(formula = form1, data = newDat, mix_name = mix_name, 
                              q = nq, b = boots, rs = rs, validation = 0, 
                              plan_strategy = plan_strategy, b1_pos = b1_pos, 
                              b_constr = b_constr, n_workers = nworkers, ...))
      }, error = function(e) NULL)
      
      if (is.null(gwqs1))
        lm1 <- NULL else lm1 <- gwqs1$fit
      if (is.null(lm1)) {
        retvec <- NA
      } else {
        retvec <- lm1$coef[2]
      }
      return(retvec)
    }
    
    betas <- pbapply::pbapply(reorgmat, 2, getbetas)
    
    calculate_pval <- function(x, true, posb1 = b1_pos) {
      if (posb1) {
        length(which(x > true))/length(betas)
      } else {
        length(which(x < true))/length(betas)
      }
    }
    
    pval <- calculate_pval(betas, ref_beta1, b1_pos)
    
    perm_retlist <- list(pval = pval, testbeta1 = ref_beta1, betas = betas, 
                         call = cl)
    
    model$b1_pos <- b1_pos
    perm_ref_wqs$b1_pos <- b1_pos
    perm_ref_wqs$seed <- seed
    
    if (boots == length(model$bindex)){
      ret_ref_wqs <- NULL
    } 
    else{
      ret_ref_wqs <- perm_ref_wqs
    }
  } else {
    Data <- model$data[, -which(names(model$data) %in% c("wqs", "wghts"))]
    
    initialfit <- function(m) {
      if(length(mm$coef) > 2){
        newform <- formula(paste0(m, "~", gsub("wqs + ", "", formchar[3], 
                                               fixed = TRUE)))
      } else {
        newform <- formula(paste0(m, "~1"))
      }
      
      fit.x1 <- lm(newform, data = Data)
      return(resid(fit.x1))
    }
    
    residmat <- sapply(model$mix_name, initialfit)
    Data[, model$mix_name] <- residmat
    
    lwqs1 <- tryCatch({
      suppressWarnings(gwqs_hpc(formula = formula(mm), data = Data, 
                            mix_name = model$mix_name, q = nq, b = boots, 
                            rs = rs, validation = 0, 
                            plan_strategy = plan_strategy, b1_pos = b1_pos, 
                            family = model$family, seed = seed,
                            b_constr = b_constr, n_workers = nworkers, ...))
    }, error = function(e) NULL)
    
    fit1 <- lwqs1$fit
    if(length(mm$coef) > 2){
      fit2form <- formula(paste0(yname, "~", gsub("wqs + ", "", formchar[3], 
                                             fixed = TRUE)))
    } else {
      fit2form <- formula(paste0(yname, "~1"))
    }
    fit2 <- glm(fit2form, data = Data, family = model$family$family)
    
    p.value.obs <- 1 - pchisq(abs(fit1$deviance - fit2$deviance), 1)
    
    reorgmatlist <- lapply(1:niter, function(x) residmat[sample(1:nrow(residmat), 
                                                                replace = F), ])
    
    getperms <- function(x) {
      newDat <- Data
      newDat[, model$mix_name] <- x
      formchar <- as.character(formula(mm))
      if (length(mm$coef) > 2) {
        form1 <- formula(paste0(formchar[2], formchar[1], formchar[3]))
      } else {
        form1 <- formula(paste0(formchar[2], formchar[1], "wqs"))
      }
      
      
      gwqs1 <- tryCatch({
        suppressWarnings(gwqs_hpc(formula = form1, data = newDat, 
                              mix_name = mix_name, q = model$q, b = boots, 
                              rs = rs, validation = 0, 
                              plan_strategy = plan_strategy, 
                              b1_pos = b1_pos, family = model$family$family,
                              b_constr = b_constr, n_workers = nworkers, ...)
          )}, error = function(e) NULL)
      
      if (is.null(gwqs1))
        lm1 <- NULL
      else
        lm1 <- gwqs1$fit
      if (is.null(lm1)) {
        retvec <- NA
      } else {
        devi <- lm1$deviance
        pperm <- 1 - pchisq(abs(devi - fit2$deviance), 1)
        retvec <- pperm
      }
      return(retvec)
    }
    
    permstats <- pbapply::pbsapply(reorgmatlist, getperms)
    
    p0 <- length(permstats[which(permstats <= p.value.obs)]) / niter

    perm_retlist <- list(pval = p0, 
                         testpval = p.value.obs,
                         permpvals = permstats)
    
    model$b1_pos <- b1_pos

    ret_ref_wqs <- NULL 
  }
  
  results <- list(gwqs_main = model, 
                  family = model$family$family,
                  gwqs_perm = ret_ref_wqs, 
                  perm_test = perm_retlist)
  
  class(results) <- "wqs_pt"
  
  results
}

#' @rawNamespace S3method(print, wqs_pt)
print.wqs_pt <- function(x, ...){
  
  cat("Permutation test WQS coefficient p-value: \n", 
      x$perm_test$pval,
      "\n")
  
  main_sum <- summary(x$gwqs_main)
  
  print(main_sum)
  
}

#' @rawNamespace S3method(summary, wqs_pt)
summary.wqs_pt <- function(object, ...){
  
  message("Permutation test WQS coefficient p-value: \n", 
      object$perm_test$pval,
      "\n")
  
  main_sum <- summary(object$gwqs_main)
  
  main_sum
  
}

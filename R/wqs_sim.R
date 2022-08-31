#' WQS simulated dataset generator
#' 
#' \code{wqs_sim} generates a simulated dataset of mixture components, covariates, 
#' and outcomes based on an initial set of specifications. 
#'
#' @param nmix Number of mixture components in simulated dataset.
#' @param ncovrt Number of covariates in simulated dataset.
#' @param nobs Number of observations in simulated dataset.
#' @param ntruewts Number of mixture components that have a non-zero association 
#' with the outcome (i.e., are not noise). 
#' @param ntruecovrt Number of covariates that have a non-zero association with 
#' the outcome (i.e., are not noise).
#' @param corrstruct Correlation matrix.
#' @param eps Dispersion parameter. If the type is "gaussian", this corresponds
#' to the residual standard deviation. If the type is "binomial",
#' this parameter is ignored.
#' @param truewqsbeta Simulated WQS beta_1 value. If NULL, then this value will 
#' be randomly sampled depending on the parameter rnd_wqsbeta_dir. 
#' @param truebeta0 Simulated beta_0 value. If NULL, then this value will be
#' randomly sampled from a standard normal distribution. 
#' @param truewts Simulated vector of mixture weights. If NULL, then this value 
#' will be randomly sampled from a Dirichlet distribution with a vector of alpha
#' values all equal to 1. 
#' @param truegamma Simulated gamma vector. If NULL, then this value will be
#' randomly sampled from a standard normal distribution. 
#' @param rnd_wqsbeta_dir Direction of randomly sampled truewqsbeta (if 
#' truewqsbeta = NULL). You can choose between "positive", "negative", or NULL. 
#' If "positive" or "negative", the truewqsbeta will be sampled from a standard
#' half normal distribution in either of those respective directions. If NULL, 
#' then truewqsbeta will be sampled from a standard normal distribution. 
#' @param seed Random seed.
#' @param q Number of quantiles. 
#' @param type Outcome family ("gaussian" for continuous outcomes or "binomial" 
#' for binary outcomes).
#'
#' @return \code{wqs_perm} returns a list of:
#' \item{weights}{Simulated weights.}
#' \item{coef}{Simulated beta coefficients.}
#' \item{Data}{Simulated dataset.}
#' \item{yhat}{simulated predicted y values from the data generating model.}
#' \item{wqs}{Weighted quantile sum vector (quantile-transformed mixture 
#' components multiplied by weights).}
#' \item{modmat}{Model matrix.}
#' \item{Xq}{Quantile-transformed mixture components.}
#' 
#' @import mvtnorm extraDistr
#' @export wqs_sim
#'
#' @examples
#'
#' # For these examples, we only run a GLM using the simulated dataset
#' # including the simulated WQS vector just to show that the user-specified
#' # coefficients for beta1 and beta0 are returned. An example of running
#' # the full permutation test WQS regression for the simulated dataset
#' # (for which the WQS vector would be determined by the model)
#' # with the "gaussian" family is shown as well.
#' 
#' wqsform<-formula(paste0("y~wqs+",paste(paste0("C",1:10),collapse="+")))
#' 
#' testsim_gaussian<-
#'   wqs_sim(truewqsbeta=0.2,truebeta0=-2,
#'           truewts=c(rep(0.15,5),rep(0.05,5)),type="gaussian")
#' Dat<-testsim_gaussian$Data
#' Dat$wqs<-testsim_gaussian$wqs
#' summary(glm(wqsform,data=Dat,family="gaussian"))$coef[1:2,]
#' perm_test_res <- wqs_full_perm(formula = wqsform, data = testsim_gaussian$Data, 
#'                                mix_name = paste0("T",1:10), q = 10, b_main = 5, 
#'                                b_perm = 5, b1_pos = TRUE, b1_constr = FALSE, 
#'                                niter = 4, seed = 16, plan_strategy = "multicore", 
#'                                stop_if_nonsig = FALSE)
#' # Note: The default values of b_main = 1000, b_perm = 200, and niter = 200 
#' # are the recommended parameter values. This example has a lower b_main, 
#' # b_perm, and niter in order to serve as a shorter test run. 
#' 
#' testsim_logit<-
#'   wqs_sim(truewqsbeta=0.2,truebeta0=-2,
#'           truewts=c(rep(0.15,5),rep(0.05,5)),type="binomial")
#' Dat<-testsim_logit$Data
#' Dat$wqs<-testsim_logit$wqs
#' summary(glm(wqsform,data=Dat,family="binomial"))$coef[1:2,]
#'
#'
wqs_sim <- function(nmix = 10, ncovrt = 10, nobs = 500, ntruewts = 10, 
                    ntruecovrt = 5, corrstruct = 0, eps = 1, truewqsbeta = NULL, 
                    truebeta0 = NULL, truewts = NULL, truegamma = NULL, 
                    rnd_wqsbeta_dir = "none", seed = 101, q = 10, 
                    type = "gaussian") {
  
  if (!type %in% c("gaussian", "binomial")){
    stop("This simulation function can only continuous (type = 'gaussian') or 
         binary (type = 'binomial') outcomes.")
  }
  
  if (length(corrstruct) == 1) {
    Rho <- diag(nmix + ncovrt)
    Rho[upper.tri(Rho)] <- Rho[lower.tri(Rho)] <- corrstruct
  } else {
    Rho <- corrstruct
  }
  
  weights <- rep(0, nmix)
  if (is.null(truewts)) {
    set.seed(seed)
    truewts <- extraDistr::rdirichlet(1, rep(1, ntruewts))
    weights[1:ntruewts] <- truewts
  } else {
    if (length(truewts) == nmix & sum(abs(truewts)) != 1) {
      truewts <- truewts / sum(truewts)
    }
    if (length(truewts) < nmix) {
      weights[1:length(truewts)] <- truewts
      weights[(length(truewts) + 1):nmix] <-
        (1 - sum(truewts)) / (nmix - length(truewts))
    } else {
      weights[1:length(truewts)] <- truewts
    }
  }
  
  if (round(sum(weights), 3) != 1.0) {
    warning(print(paste0("weights add up to ", sum(weights))))
  }
  
  set.seed(seed)
  Xmat <- mvtnorm::rmvnorm(nobs, mean = rep(0, nmix + ncovrt), sigma = Rho)
  if (is.null(q)) {
    Xmatquant <- Xmat
  } else {
    Xmatquant <- Xmat
    Xmatquant[, 1:nmix] <- apply(
      Xmatquant[, 1:nmix],
      2,
      FUN = function(x) {
        as.numeric(as.character(cut(
          x,
          breaks = quantile(x, probs = seq(0, 1, by = (1 / q))),
          include.lowest = T,
          labels = 0:(q - 1)
        )))
      }
    )
  }
  
  if (ncovrt < ntruecovrt) {
    ntruecovrt <- ncovrt
  }
  if (!is.null(truegamma)) {
    if (length(truegamma) == 1) {
      covrtbetas <- rep(truegamma, ncovrt)
    } else {
      covrtbetas <- truegamma
    }
  } else {
    set.seed(seed)
    covrtbetas <- c(rnorm(ntruecovrt), rep(0, length = ncovrt - ntruecovrt))
  }
  
  set.seed(seed)
  if (!is.null(truebeta0)) {
    beta0 <- truebeta0
  } else {
    beta0 <- rnorm(1)
  }
  
  if (!is.null(truewqsbeta)) {
    wqsbeta <- truewqsbeta
  } else {
    set.seed(seed)
    if (rnd_wqsbeta_dir == "positive") {
      wqsbeta <- extraDistr::rhnorm(1)
    } else if (rnd_wqsbeta_dir == "negative") {
      wqsbeta <- extraDistr::rhnorm(1) * -1
    } else {
      wqsbeta <- rnorm(1)
    }
  }
  
  wqs <- Xmatquant[, 1:nmix] %*% weights
  if (ncovrt > 0) {
    modmat <- cbind(1, wqs, Xmat[, c((nmix + 1):(nmix + ncovrt))])
    dimnames(modmat)[[2]] <-
      c("Intercept", "wqs", paste0("C", 1:ncovrt))
    betas <- c(beta0, wqsbeta, covrtbetas)
    names(betas) <- c("beta0", "beta1", paste0("gamma", 1:ncovrt))
  } else {
    modmat <- cbind(1, wqs)
    dimnames(modmat)[[2]] <-
      dimnames(modmatq)[[2]] <- c("Intercept", "wqs")
    betas <- c(beta0, wqsbeta)
    names(betas) <- c("beta0", "beta1")
  }
  
  if (type == "binomial"){
    etahat <- modmat %*% betas
    probs <- 1 / (1 + exp(-etahat))
    set.seed(seed)
    y <- rbinom(nobs, size = 1, prob = probs)
    Data <- data.frame(cbind(y, Xmat))
    yhat <- NULL
  }
  else{
    yhat <- modmat %*% betas
    set.seed(seed)
    epsilon <- rnorm(nobs, sd = eps)
    y <- yhat + epsilon
    Data <- data.frame(cbind(y, Xmat))    
    etahat <- NULL 
  }
  
  if (ncovrt > 0) {
    names(Data) <- c("y", paste0("T", 1:nmix), paste0("C", 1:ncovrt))
    colnames(Xmatquant) <-
      c(paste0("T", 1:nmix), paste0("C", 1:ncovrt))
    colnames(Xmat) <- c(paste0("T", 1:nmix), paste0("C", 1:ncovrt))
  } else {
    names(Data) <- c("y", paste0("T", 1:nmix))
    colnames(Xmatquant) <- c(paste0("T", 1:nmix))
    colnames(Xmat) <- c(paste0("T", 1:nmix))
  }
  
  wtmat <- data.frame(mix_name = paste0("T", 1:nmix), true_weight = weights)
  
  retlist <-
    list(
      weights = wtmat,
      coef = betas,
      Data = Data,
      yhat = yhat,
      etahat = etahat,
      wqs = wqs,
      modmat = modmat,
      Xq = Xmatquant
    )
  return(retlist)
}

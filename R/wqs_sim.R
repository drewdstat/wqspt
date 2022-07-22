#' WQS simulated dataset generator
#' 
#' \code{wqs_sim} generates a simulated dataset of mixture components, covariates, 
#' and outcomes based on an initial set of specifications. 
#'
#' @param nmix Number of mixture components in simulated dataset.
#' @param ncovrt Number of covariates in simulated dataset.
#' @param nobs Number of observations in simulated dataset.
#' @param ntruewts Number of mixture components that have a non-zero association 
#' with the outcome (i.e. are not noise). 
#' @param ntruecovrt Number of covariates that have a non-zero association with 
#' the outcome (i.e. are not noise).
#' @param corrstruct Correlation matrix.
#' @param eps Dispersion parameter. If the family is gaussian(), this corresponds
#' to the residual standard deviation. If the family is binomial() or poisson(),
#' this parameter is ignored. If the family is "negbin", then this will be the 
#' "size" parameter in the rnbinom() function. If the family is quasipoisson(),
#' an rnbinom() distribution with "size" equal to exp(etahat)/(eps-1) is used to
#' generate the outcome variable.
#' @param truewqsbeta Simulated WQS beta_1 value. If NULL, then this value will 
#' be randomly sampled. 
#' @param truebeta0 Simulated beta_0 value. If NULL, then this value will be
#' randomly sampled. 
#' @param truewts Simulated vector of mixture weights. If NULL, then this value 
#' will be randomly sampled. 
#' @param truegamma Simulated gamma (i.e., covariate coefficients) vector. 
#' If NULL, then these values will be randomly sampled. 
#' @param rnd_wqsbeta_dir Direction of randomly sampled truewqsbeta (if 
#' truewqsbeta = NULL). You can choose between "positive", "negative", or NULL. 
#' If "positive" or "negative", the truewqsbeta will be sampled from a half 
#' normal distribution in either of those respective directions. If NULL, then 
#' truewqsbeta will be sampled from a normal distribution. 
#' @param seed Random seed.
#' @param q Number of quantiles. 
#' @param family Family for the generative model creating the outcome vector. 
#' Options include "gaussian" or gaussian(link = "identity") for a continuous
#' outcome, "binomial" or binomial() with any accepted link function for a binary
#' outcome, and finally for count outcomes this can be "poisson" or 
#' poisson(link="log"), "negbin" for negative binomial, and "quasipoisson" or
#' quasipoisson(link="log").
#'
#' @return \code{wqs_perm} returns a list of:
#' \item{weights}{Simulated weights.}
#' \item{coef}{Simulated beta coefficients.}
#' \item{Data}{Simulated dataset.}
#' \item{etahat}{Linear predictor values from the data generating model.}
#' \item{wqs}{Weighted quantile sum vector (quantile-transformed mixture 
#' components multiplied by weights and then summed).}
#' \item{modmat}{Model matrix.}
#' \item{Xq}{Quantile-transformed mixture components.}
#' 
#' @import mvtnorm extraDistr
#' @export wqs_sim
#'
#' @examples
#' 
#' # For these examples, I'll just run a GLM using the simulated dataset
#' # including the simulated WQS vector just to show that the user-specified
#' # coefficients for beta1 and beta0 are returned. An example of running
#' # the full permutation test WQS regression for the simulated dataset
#' # (for which the WQS vector would be determined by the model)
#' # with the gaussian() family is shown as well.
#' 
#' wqsform<-formula(paste0("y~wqs+",paste(paste0("C",1:10),collapse="+")))
#' 
#' testsim_gaussian<-
#'   wqs_sim(truewqsbeta=0.2,truebeta0=-2,
#'           truewts=c(rep(0.15,5),rep(0.05,5)),family="gaussian")
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
#'           truewts=c(rep(0.15,5),rep(0.05,5)),family="binomial")
#' Dat<-testsim_logit$Data
#' Dat$wqs<-testsim_logit$wqs
#' summary(glm(wqsform,data=Dat,family="binomial"))$coef[1:2,]
#' 
#' testsim_probit<-wqs_sim(truewqsbeta=0.2,truebeta0=-2,
#'                         truewts=c(rep(0.15,5),rep(0.05,5)),
#'                         family=binomial(link="probit"))
#' wqsform<-formula(paste0("y~wqs+",paste(paste0("C",1:10),collapse="+")))
#' Dat<-testsim_probit$Data
#' Dat$wqs<-testsim_probit$wqs
#' summary(glm(wqsform,data=Dat,family=binomial(link="probit")))$coef[1:2,]
#' 
#' testsim_nb<-wqs_sim(truewqsbeta=0.2,truebeta0=0,
#'                     truewts=c(rep(0.15,5),rep(0.05,5)),
#'                     family="negbin",eps=1)
#' #var = exp(eta) + (exp(eta)^2/eps)
#' Dat<-testsim_nb$Data
#' Dat$wqs<-testsim_nb$wqs
#' summary(glm.nb(wqsform,data=Dat))$coef[1:2,]
#' 
#' testsim_qp<-wqs_sim(truewqsbeta=0.2,truebeta0=0,
#'                     truewts=c(rep(0.15,5),rep(0.05,5)),
#'                     family="quasipoisson",eps=1)
#' #var = exp(eta) + (exp(eta)^2/(exp(eta)/(eps-1)))
#' Dat<-testsim_qp$Data
#' Dat$wqs<-testsim_qp$wqs
#' summary(glm(wqsform,data=Dat,family="quasipoisson"))$coef[1:2,]
#' 
wqs_sim <- function(nmix = 10, ncovrt = 10, nobs = 500, ntruewts = 10, 
                    ntruecovrt = 5, corrstruct = 0, eps = 1, truewqsbeta = NULL, 
                    truebeta0 = NULL, truewts = NULL, truegamma = NULL, 
                    rnd_wqsbeta_dir = "none", seed = 101, q = 10, 
                    family = gaussian(link="identity")) {
  
  if (is.character(family)) {
    if (family=="multinomial"){
      stop("This simulation function doesn't yet accomodate 
           multinomial WQS regression.")
    }
    if (family %in% c("negbin")) 
      family <- list(family = family)
    else family <- get(family, mode = "function", 
                       envir = parent.frame())
  }
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized\n")
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
  Xmat <- rmvnorm(nobs, mean = rep(0, nmix + ncovrt), sigma = Rho)
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
  
  etahat <- modmat %*% betas
  if(family$family=="binomial"){
    probs <- family$linkinv(etahat)
    set.seed(seed)
    y <- rbinom(nobs, size = 1, prob = probs)
  } else if(family$family=="gaussian"){
    if(family$link!="identity") { 
      stop("The gaussian() family is only supported for 
           link='identity' for WQS regression")
    }
    set.seed(seed)
    epsilon <- rnorm(nobs, sd = eps)
    y <- etahat + epsilon
  } else if(family$family=="poisson"){
    if(family$link!="log") { 
      stop("The poisson() family is only supported for 
           link='log' for WQS regression")
    }
    y<-rpois(nobs,lambda = exp(etahat))
  } else if(family$family=="quasipoisson"){
    if(family$link!="log") { 
      stop("The quasipoisson() family is only supported for 
           link='log' for WQS regression")
    }
    rqpois <- function(n, mu, theta) {
      rnbinom(n = n, mu = mu, size = mu/(theta-1))
    }
    y <- rqpois(nobs, exp(etahat), eps)
  } else if(family$family=="negbin"){
    y <- rnbinom(n = nobs, mu = exp(etahat), size = eps)
  }
  
  Data <- data.frame(cbind(y, Xmat))
  
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
      etahat = etahat,
      wqs = wqs,
      modmat = modmat,
      Xq = Xmatquant
    )
  return(retlist)
}

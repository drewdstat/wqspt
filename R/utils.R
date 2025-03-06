#' @importFrom MASS glm.nb
#' @importFrom pscl zeroinfl
#' @importFrom car vif
#' @importFrom future.apply future_lapply
#' @importFrom reshape2 melt
#' @importFrom grid is.grob
#' @importFrom utils flush.console
#' @importFrom nnet multinom
#' @import future gWQS

#work around for errors in cowplot v1.1.3 get_legend when legend is 
# not in the 'right' position
get_legend2 <- function(plot, legend = NULL) {
  if (is.ggplot(plot)) {
    gt <- ggplotGrob(plot)
  } else {
    if (is.grob(plot)) {
      gt <- plot
    } else {
      stop("Plot object is neither a ggplot nor a grob.")
    }
  }
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }
  indices <- grep(pattern, gt$layout$name)
  not_empty <- !vapply(
    gt$grobs[indices], 
    inherits, what = "zeroGrob", 
    FUN.VALUE = logical(1)
  )
  indices <- indices[not_empty]
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}

# gWQS::gwqs sets the plan() strategy internally in the function, but doesn't 
# provide any arguments for changing the maximum number of workers or other 
# important plan() parameters in the function call. This has the effect of 
# having no cap on the number of workers or memory usage when running WQSPT on 
# a high-performance computing queue, which can cause the whole job to quit if 
# the HPC-allotted core or memory constraints are exceeded. gwqs_hpc adds 
# an argument specifying the number of workers to use with the plan_strategy to 
# gWQS::gwqs to help ensure that these HPC issues are avoided, though excessive 
# memory usage is still a concern.
gwqs_hpc <- function (formula, data, na.action, weights, mix_name, stratified, 
                      rh = 1, b = 100, b1_pos = TRUE, bint_cont_pos = NULL, 
                      bint_cat_pos = NULL, b_constr = FALSE, zero_infl = FALSE, 
                      q = 4, validation = 0.6, validation_rows = NULL, 
                      family = gaussian, 
                      signal = c("t2", "t3", "one", "abst", "expt"), 
                      rs = FALSE, n_vars = NULL, 
                      zilink = c("logit", "probit", "cloglog", "cauchit", "log"), 
                      seed = NULL, wp = NULL, wn = NULL, 
                      plan_strategy = "sequential", 
                      lambda = 0, 
                      optim.method = c("BFGS", "Nelder-Mead", "CG", "SANN"), 
                      control = list(trace = FALSE, maxit = 2000, 
                                     reltol = 1e-09), 
                      b1_constr = NULL, n_workers = NULL, ...) 
{
  wqsGaussBin <- function(initp, kw, bXm, bY, boffset, bQ, 
                          kx, Xnames, level_names, wqsvars, pwqsvars, nwqsvars, 
                          family, dwqs, zilink, zero_infl, formula, ff, bwghts, 
                          stratified, stratlev, b1_pos, b_constr, lambda, intpos, 
                          wqsint, intvars, bint_cont_pos, bint_cat_pos, intcont) {
    if (b_constr) {
      initp["wqs"] <- ifelse(b1_pos, abs(initp["wqs"]), 
                             -abs(initp["wqs"]))
      if (wqsint) {
        if (intcont) 
          initp[intvars] <- sapply(1:length(intvars), 
                                   function(i) ifelse(bint_cont_pos[i], abs(initp[intvars[i]]), 
                                                      -abs(initp[intvars[i]])))
        else initp[intvars] <- sapply(1:length(intvars), 
                                      function(i) ifelse(bint_cat_pos[i] & initp[intvars[i]] - 
                                                           initp["wqs"] < 0, abs(initp[intvars[i]]), 
                                                         ifelse(!bint_cat_pos[i] & initp[intvars[i]] - 
                                                                  initp["wqs"] > 0, -abs(initp[intvars[i]]), 
                                                                abs(initp[intvars[i]]))))
      }
    }
    w <- initp[(kx + 1):(kx + kw)]
    pen <- sum(abs(w))
    w <- w^2
    if (!is.null(stratified)) {
      w <- matrix(w, kw/stratlev, stratlev)
      w <- apply(w, MARGIN = 2, FUN = function(i) i/sum(i))
      w <- as.numeric(w)
    }
    else w <- w/sum(w)
    bXm[, "wqs"] <- as.numeric(bQ %*% w)
    if (wqsint) 
      bXm[, intpos] <- bXm[, "wqs"] * bXm[, intvars]
    b_covs <- initp[1:kx]
    term <- as.numeric(bXm %*% b_covs) + boffset
    f = sum(family$dev.resids(y = bY, mu = family$linkinv(term), 
                              wt = bwghts)) + lambda * pen
    return(f)
  }
  wqsGaussBin_dwqs <- function(initp, kw, bXm, bY, boffset, 
                               bQ, kx, Xnames, level_names, wqsvars, pwqsvars, nwqsvars, 
                               family, dwqs, zilink, zero_infl, formula, ff, bwghts, 
                               stratified, stratlev, b1_pos, b_constr, lambda, intpos, 
                               wqsint, intvars, bint_cont_pos, bint_cat_pos, intcont) {
    initp["pwqs"] <- abs(initp["pwqs"])
    initp["nwqs"] <- -abs(initp["nwqs"])
    wp <- initp[(kx + 1):(kx + kw)]
    wn <- initp[(kx + kw + 1):(kx + 2 * kw)]
    pen <- sum(abs(wp)) + sum(abs(wn))
    wp <- (wp^2)/sum(wp^2)
    wn <- (wn^2)/sum(wn^2)
    bXm[, "pwqs"] <- as.numeric(bQ %*% wp)
    bXm[, "nwqs"] <- as.numeric(bQ %*% wn)
    b_covs <- initp[1:kx]
    term <- as.numeric(bXm %*% b_covs) + boffset
    f = sum(family$dev.resids(y = bY, mu = family$linkinv(term), 
                              wt = bwghts)) + lambda * pen
    return(f)
  }
  wqsPoisson <- function(initp, kw, bXm, bY, boffset, bQ, 
                         kx, Xnames, n_levels, level_names, wqsvars, pwqsvars, 
                         nwqsvars, family, dwqs, zilink, zero_infl, formula, 
                         ff, bwghts, stratified, stratlev, b1_pos, b_constr, 
                         lambda, intpos, wqsint, intvars, bint_cont_pos, bint_cat_pos, 
                         intcont) {
    if (b_constr) {
      initp["wqs"] <- ifelse(b1_pos, abs(initp["wqs"]), 
                             -abs(initp["wqs"]))
      if (wqsint) {
        if (intcont) 
          initp[intvars] <- sapply(1:length(intvars), 
                                   function(i) ifelse(bint_cont_pos[i], abs(initp[intvars[i]]), 
                                                      -abs(initp[intvars[i]])))
        else initp[intvars] <- sapply(1:length(intvars), 
                                      function(i) ifelse(bint_cat_pos[i] & initp[intvars[i]] - 
                                                           initp["wqs"] < 0, abs(initp[intvars[i]]), 
                                                         ifelse(!bint_cat_pos[i] & initp[intvars[i]] - 
                                                                  initp["wqs"] > 0, -abs(initp[intvars[i]]), 
                                                                abs(initp[intvars[i]]))))
      }
    }
    w <- initp[(kx + 1):(kx + kw)]
    pen <- sum(abs(w))
    w <- w^2
    if (!is.null(stratified)) {
      w <- matrix(w, kw/stratlev, stratlev)
      w <- apply(w, MARGIN = 2, FUN = function(i) i/sum(i))
      w <- as.numeric(w)
    }
    else w <- w/sum(w)
    bXm[, "wqs"] <- as.numeric(bQ %*% w)
    if (wqsint) 
      bXm[, intpos] <- bXm[, "wqs"] * bXm[, intvars]
    b_covs <- initp[1:kx]
    term <- as.numeric(bXm %*% b_covs) + boffset
    f = -sum(dpois(bY, lambda = exp(term), log = TRUE) * 
               bwghts) + lambda * pen
    return(f)
  }
  wqsPoisson_dwqs <- function(initp, kw, bXm, bY, boffset, 
                              bQ, kx, Xnames, n_levels, level_names, wqsvars, pwqsvars, 
                              nwqsvars, family, dwqs, zilink, zero_infl, formula, 
                              ff, bwghts, stratified, stratlev, b1_pos, b_constr, 
                              lambda, intpos, wqsint, intvars, bint_cont_pos, bint_cat_pos, 
                              intcont) {
    initp["pwqs"] <- abs(initp["pwqs"])
    initp["nwqs"] <- -abs(initp["nwqs"])
    wp <- initp[(kx + 1):(kx + kw)]^2
    wn <- initp[(kx + kw + 1):(kx + 2 * kw)]^2
    pen <- sum(wp) + sum(wn)
    wp <- wp/sum(wp)
    wn <- wn/sum(wn)
    bXm[, "pwqs"] <- as.numeric(bQ %*% wp)
    bXm[, "nwqs"] <- as.numeric(bQ %*% wn)
    b_covs <- initp[1:kx]
    term <- as.numeric(bXm %*% b_covs) + boffset
    f = -sum(dpois(bY, lambda = exp(term), log = TRUE) * 
               bwghts) + lambda * pen
    return(f)
  }
  wqsNegBin <- function(initp, kw, bXm, bY, boffset, bQ, kx, 
                        Xnames, n_levels, level_names, wqsvars, pwqsvars, nwqsvars, 
                        family, dwqs, zilink, zero_infl, formula, ff, bwghts, 
                        stratified, stratlev, b1_pos, b_constr, lambda, intpos, 
                        wqsint, intvars, bint_cont_pos, bint_cat_pos, intcont) {
    if (b_constr) {
      initp["wqs"] <- ifelse(b1_pos, abs(initp["wqs"]), 
                             -abs(initp["wqs"]))
      if (wqsint) {
        if (intcont) 
          initp[intvars] <- sapply(1:length(intvars), 
                                   function(i) ifelse(bint_cont_pos[i], abs(initp[intvars[i]]), 
                                                      -abs(initp[intvars[i]])))
        else initp[intvars] <- sapply(1:length(intvars), 
                                      function(i) ifelse(bint_cat_pos[i] & initp[intvars[i]] - 
                                                           initp["wqs"] < 0, abs(initp[intvars[i]]), 
                                                         ifelse(!bint_cat_pos[i] & initp[intvars[i]] - 
                                                                  initp["wqs"] > 0, -abs(initp[intvars[i]]), 
                                                                abs(initp[intvars[i]]))))
      }
    }
    w <- initp[(kx + 1):(kx + kw)]
    pen <- sum(abs(w))
    w <- w^2
    if (!is.null(stratified)) {
      w <- matrix(w, kw/stratlev, stratlev)
      w <- apply(w, MARGIN = 2, FUN = function(i) i/sum(i))
      w <- as.numeric(w)
    }
    else w <- w/sum(w)
    bXm[, "wqs"] <- as.numeric(bQ %*% w)
    if (wqsint) 
      bXm[, intpos] <- bXm[, "wqs"] * bXm[, intvars]
    b_covs <- initp[1:kx]
    theta <- exp(initp[length(initp)])
    term <- as.numeric(bXm %*% b_covs) + boffset
    f = -sum((suppressWarnings(dnbinom(bY, size = theta, 
                                       mu = exp(term), log = TRUE))) * bwghts) + lambda * 
      pen
    return(f)
  }
  wqsNegBin_dwqs <- function(initp, kw, bXm, bY, boffset, 
                             bQ, kx, Xnames, n_levels, level_names, wqsvars, pwqsvars, 
                             nwqsvars, family, dwqs, zilink, zero_infl, formula, 
                             ff, bwghts, stratified, stratlev, b1_pos, b_constr, 
                             lambda, intpos, wqsint, intvars, bint_cont_pos, bint_cat_pos, 
                             intcont) {
    initp["pwqs"] <- abs(initp["pwqs"])
    initp["nwqs"] <- -abs(initp["nwqs"])
    wp <- initp[(kx + 1):(kx + kw)]^2
    wn <- initp[(kx + kw + 1):(kx + 2 * kw)]^2
    pen <- sum(wp) + sum(wn)
    wp <- wp/sum(wp)
    wn <- wn/sum(wn)
    bXm[, "pwqs"] <- as.numeric(bQ %*% wp)
    bXm[, "nwqs"] <- as.numeric(bQ %*% wn)
    b_covs <- initp[1:kx]
    theta <- exp(initp[length(initp)])
    term <- as.numeric(bXm %*% b_covs) + boffset
    f = -sum((suppressWarnings(dnbinom(bY, size = theta, 
                                       mu = exp(term), log = TRUE))) * bwghts) + lambda * 
      pen
    return(f)
  }
  if (!is.null(b1_constr)) {
    b_constr <- b1_constr
    warning("The argument 'b1_constr' is deprecated, use 'b_constr' instead.")
  }
  if(is.null(n_workers)){
    n_workers <- availableWorkers()
    if(length(n_workers) > 1) n_workers <- length(n_workers)
  }
  one <- function(x) rep(1, length(x))
  t2 <- function(x) x^2
  t3 <- function(x) x^3
  expt <- function(x) exp(abs(x)) - 1
  if (is.character(family)) {
    if (family == "negbin") 
      family <- list(family = family)
    else if (family == "multinomial") 
      stop("'family = multinomial' is not supported by the 'gwqs' function. Please use 'gwqs_multinom' function.\n")
    else family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized\n")
  }
  wqs_in_formula <- any(grepl("wqs", rownames(attr(terms(formula), 
                                                   "factors"))))
  pwqs_in_formula <- any(grepl("pwqs", rownames(attr(terms(formula), 
                                                     "factors"))))
  nwqs_in_formula <- any(grepl("nwqs", rownames(attr(terms(formula), 
                                                     "factors"))))
  if (!wqs_in_formula & (!pwqs_in_formula | !nwqs_in_formula)) 
    stop("'formula' must contain either 'wqs' term or 'pwqs' and 'nwqs' terms: e.g. y ~ wqs + ...\n y ~ pwqs + nwqs + ...\n")
  dwqs <- pwqs_in_formula & nwqs_in_formula
  if (dwqs & rs) 
    stop("The 2iWQS does not support yet the random subset method\n")
  if (dwqs) 
    objfn <- switch(family$family, gaussian = wqsGaussBin_dwqs, 
                    binomial = wqsGaussBin_dwqs, poisson = wqsPoisson_dwqs, 
                    quasipoisson = wqsPoisson_dwqs, negbin = wqsNegBin_dwqs)
  else objfn <- switch(family$family, gaussian = wqsGaussBin, 
                       binomial = wqsGaussBin, poisson = wqsPoisson, quasipoisson = wqsPoisson, 
                       negbin = wqsNegBin)
  optim.method <- match.arg(optim.method)
  signal <- match.arg(signal)
  signal <- switch(signal, one = one, abst = abs, t2 = t2, 
                   t3 = t3, expt = expt)
  if (zero_infl) {
    zilink <- make.link(match.arg(zilink))
    ff = formula
    if (length(formula[[3]]) > 1 && identical(formula[[3]][[1]], 
                                              as.name("|"))) 
      formula = as.formula(paste0(ff[[2]], " ~ ", deparse(ff[[3]][[2]])))
  }
  else zilink <- ff <- NULL
  cl <- match.call()
  mc <- match.call(expand.dots = FALSE)
  m <- match(c("data", "na.action", "formula", "mix_name", 
               "weights", "stratified"), names(mc), 0)
  dtf <- mc[c(1, m)]
  dtf[[2]] <- data
  dtf[[1]] <- as.name("selectdatavars")
  dtf <- eval(dtf, parent.frame())
  l <- list(...)
  solve_dir_issue <- ifelse(is.null(l$solve_dir_issue), FALSE, 
                            l$solve_dir_issue)
  if (is.null(q)) {
    Q = as.matrix(dtf[, mix_name])
    qi = NULL
  }
  else {
    q_f = gwqs_rank(dtf, mix_name, q)
    Q = q_f$Q
    qi = q_f$qi
  }
  m <- match(c("stratified"), names(mc), 0)
  if (!is.null(l$stratified_rh)) 
    m <- 0
  if (m) {
    strtfd_out <- stratified_f(Q, dtf, stratified, mix_name)
    Q <- strtfd_out$Q
    mix_name <- strtfd_out$mix_name
  }
  else stratified <- NULL
  if (!is.null(l$stratified_rh)) 
    stratified <- l$stratified_rh
  stratlev <- ifelse(is.null(stratified), 1, nlevels(unlist(dtf[, 
                                                                stratified])))
  if (dwqs & !is.null(stratified)) 
    stop("The 2iWQS does not support yet the presence of stratified weights\n")
  if (!is.null(seed) & !is.numeric(seed)) 
    stop("seed must be numeric or NULL\n")
  if (!is.null(seed)) 
    set.seed(seed)
  fseed <- TRUE
  N <- nrow(dtf)
  if (!is.numeric(validation) | validation < 0 | validation >= 
      1) 
    stop("'validation' must be numeric >= 0 and < 1\n")
  if (is.null(validation_rows)) 
    validation_rows <- lapply(1:rh, function(i) sample(c(F, 
                                                         T), N, replace = T, prob = c(1 - validation, validation)))
  else if (length(validation_rows) != rh | !any(unique(unlist(validation_rows)) %in% 
                                                c(F, T))) 
    stop("'validation_rows' must be a list of length 'rh' of logical vectors\n")
  if (dwqs) {
    mc_pwqs <- mc_nwqs <- mc
    mc_pwqs$data <- mc_nwqs$data <- data
    mc_pwqs$mix_name <- mc_nwqs$mix_name <- mix_name
    mc_pwqs$validation <- mc_nwqs$validation <- 0
    mc_pwqs$plan_strategy <- mc_nwqs$plan_strategy <- "sequential"
    mc_pwqs$b_constr <- mc_nwqs$b_constr <- TRUE
    mc_pwqs$b1_pos <- TRUE
    mc_nwqs$b1_pos <- FALSE
    mc_pwqs$validation_rows <- mc_nwqs$validation_rows <- NULL
    if (rs) 
      mc_pwqs$b <- mc_nwqs$b <- 100
    else mc_pwqs$b <- mc_nwqs$b <- 1
    mc_pwqs$rh <- mc_nwqs$rh <- 1
    mc_pwqs$formula <- mc_nwqs$formula <- as.formula(gsub("pwqs", 
                                                          "wqs", deparse(formula)))
    mc_pwqs$formula <- mc_nwqs$formula <- remove_terms(mc_pwqs$formula, 
                                                       "nwqs")
    mc_pwqs$solve_dir_issue <- mc_nwqs$solve_dir_issue <- "inverse"
    if (is.null(wp)) {
      pwqs0 <- tryCatch(eval(mc_pwqs), error = function(e) NULL)
      if (is.null(pwqs0)) 
        wp <- rep(1/length(mix_name), length(mix_name))
      else wp <- pwqs0$final_weights$mean_weight[match(mix_name, 
                                                       pwqs0$final_weights$mix_name)]
    }
    if (is.null(wn)) {
      nwqs0 <- tryCatch(eval(mc_nwqs), error = function(e) NULL)
      if (is.null(nwqs0)) 
        wn <- rep(1/length(mix_name), length(mix_name))
      else wn <- nwqs0$final_weights$mean_weight[match(mix_name, 
                                                       nwqs0$final_weights$mix_name)]
    }
  }
  if (is.null(n_vars)) 
    n_vars = round(sqrt(length(mix_name)))
  m <- match(c("weights"), names(mc), 0)
  if (m[1]) 
    dtf$wghts <- wghts <- unlist(dtf[, weights, drop = FALSE])
  else dtf$wghts <- wghts <- rep(1, N)
  if (family$family %in% c("gaussian", "quasipoisson")) 
    ts = "t"
  else if (family$family %in% c("binomial", "poisson", "multinomial", 
                                "negbin")) 
    ts = "z"
  if (!is.numeric(b)) 
    stop("'b' must be a number\n")
  Xnames <- parnames(dtf, formula, NULL)
  kx <- length(Xnames)
  kw <- ncol(Q)
  initp <- values.0(kw, Xnames, kx, formula, ff, wghts, dtf, 
                    stratified, stratlev, b1_pos, family, wp, wn, dwqs, 
                    zilink, zero_infl)
  mf <- model.frame(formula, dtf)
  Y <- model.response(mf, "any")
  if (family$family == "binomial" & any(class(Y) %in% c("factor", 
                                                        "character"))) {
    if (inherits(Y, "character")) 
      Y = factor(Y)
    Y <- as.numeric(Y != levels(Y)[1])
  }
  offset <- model.offset(mf)
  if (is.null(offset)) 
    offset <- rep(0, nrow(dtf))
  Xm <- model.matrix(formula, dtf)
  intpos <- grepl("wqs", colnames(Xm)) & grepl(":", colnames(Xm))
  intvarname <- colnames(Xm)[intpos]
  wqsint <- any(intpos)
  intvars <- gsub(":wqs", "", gsub("wqs:", "", colnames(Xm)[intpos]))
  if (wqsint) {
    intvar0 <- gsub(":wqs", "", gsub("wqs:", "", attr(terms(formula), 
                                                      "term.labels")[grepl("wqs:", attr(terms(formula), 
                                                                                        "term.labels")) | grepl(":wqs", attr(terms(formula), 
                                                                                                                             "term.labels"))]))
    intcont <- is.numeric(dtf[, intvar0])
    if (wqsint & ((is.null(bint_cont_pos) & intcont) | (is.null(bint_cat_pos) & 
                                                        !intcont))) 
      stop(paste0("Please specify the direction of the interaction term(s) thourgh the parameter ", 
                  ifelse(intcont, "bint_cont_pos\n", "bint_cat_pos\n")))
  }
  else bint_cont_pos <- bint_cat_pos <- intcont <- NULL
  if (dwqs & wqsint) 
    stop("The 2iWQS does not support yet the presence of an interaction between a covariate and the WQS index\n")
  plan_strategy_rh <- ifelse(rh == 1, "sequential", plan_strategy)
  plan_strategy_b <- ifelse(rh == 1, plan_strategy, "sequential")
  if(plan_strategy_rh == "sequential") plan(plan_strategy_rh) else 
    plan(plan_strategy_rh, workers = n_workers)
  gwqslist <- future_lapply(X = 1:rh, FUN = function(i) {
    if (control$trace) 
      cat("start opt\n")
    if(plan_strategy_b == "sequential") plan(plan_strategy_b) else 
      plan(plan_strategy_b, workers = n_workers)
    param <- future_lapply(X = 1:b, FUN = optim.f, objfn = objfn, 
                           Y = Y[!validation_rows[[i]]], Xm = Xm[!validation_rows[[i]], 
                           ], Q = Q[!validation_rows[[i]], ], offset = offset[!validation_rows[[i]]], 
                           wghts = wghts[!validation_rows[[i]]], initp = initp, 
                           b1_pos = b1_pos, b_constr = b_constr, n_vars = n_vars, 
                           dwqs = dwqs, family = family, rs = rs, zilink = zilink, 
                           zero_infl = zero_infl, formula = formula, ff = ff, 
                           kx = kx, kw = kw, Xnames = Xnames, stratified = stratified, 
                           b = b, optim.method = optim.method, control = control, 
                           lambda = lambda, stratlev = stratlev, intpos = intpos, 
                           wqsint = wqsint, intvars = intvars, bint_cont_pos = bint_cont_pos, 
                           bint_cat_pos = bint_cat_pos, intcont = intcont, 
                           future.seed = fseed)
    conv <- c(sapply(param, function(i) i$conv))
    counts <- c(sapply(param, function(i) i$counts))
    val <- c(sapply(param, function(i) i$val))
    mex <- lapply(param, function(i) i$mex)
    bindex <- lapply(param, function(i) i$bindex)
    slctd_vars <- lapply(param, function(i) i$slctd_vars)
    if (rs) {
      if(plan_strategy == "sequential") plan(plan_strategy) else 
        plan(plan_strategy, workers = n_workers)
      param <- future_lapply(X = 1:b, FUN = set_par_names, 
                             slctd_vars, param, q_name = colnames(Q), family = family, 
                             dwqs = dwqs, future.seed = FALSE)
    }
    build_bres <- function(wqs_var_name, interval) {
      wght_matrix <- do.call("rbind", lapply(param, function(i) i$par_opt[interval]))
      b1 <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients[wqs_var_name, 
                                                                       "Estimate"])
      if (wqsint) {
        bint <- as.data.frame(do.call("rbind", lapply(param, 
                                                      function(i) summary(i$mfit$m_f)$coefficients[intpos, 
                                                                                                   "Estimate"])))
        if (!intcont) 
          bint <- b1 + bint
        names(bint) <- intvarname
      }
      se <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients[wqs_var_name, 
                                                                       "Std. Error"])
      stat <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients[wqs_var_name, 
                                                                         paste0(ts, " value")])
      p_val <- sapply(param, function(i) summary(i$mfit$m_f)$coefficients[wqs_var_name, 
                                                                          gsub("x", ts, "Pr(>|x|)")])
      if (dwqs) {
        tol <- sapply(param, function(i) {
          tmp <- vif(i$mfit$m_f)
          if ("matrix" %in% class(tmp)) 
            tmp <- tmp[, 1]
          1/tmp[wqs_var_name]
        })
      }
      n_non_conv = sum(conv == 1)
      if (n_non_conv == 0 & control$trace) 
        cat(paste0("The optimization function always converged\n"))
      else if (n_non_conv == b & rh == 1) 
        stop("The optimization function never converged\n")
      else if (n_non_conv == b & rh > 1) 
        return(NULL)
      else if (control$trace) 
        cat(paste0("The optimization function did not converge ", 
                   n_non_conv, " time/times\n"))
      if (control$trace) 
        cat(paste0("There are ", ifelse(b1_pos, sum(b1 >= 
                                                      0, na.rm = T), sum(b1 <= 0, na.rm = T)), ifelse(b1_pos, 
                                                                                                      " positive", " negative"), " bootstrapped b1 out of ", 
                   b, "\n"))
      bres <- as.data.frame(cbind(wght_matrix, b1, se, 
                                  stat, p_val))
      names(bres) <- c(colnames(Q), "b1", "Std_Error", 
                       "stat", "p_val")
      if (wqsint) 
        bres <- cbind(bres, bint)
      if (dwqs) 
        bres$tol <- tol
      strata_names <- NULL
      return(bres)
    }
    if (dwqs) {
      bresp <- build_bres("pwqs", 1:ncol(Q))
      bresn <- build_bres("nwqs", (ncol(Q) + 1):(2 * ncol(Q)))
      bres <- list(bresp = bresp, bresn = bresn)
    }
    else {
      bresp <- bresn <- NULL
      bres <- build_bres("wqs", 1:ncol(Q))
    }
    if (is.null(bres)) 
      return(NULL)
    mean_weight_f <- function(bres, b1_pos, par) {
      if (rs) 
        bres[mix_name][is.na(bres[mix_name])] <- 0
      filter_direction <- apply(as.matrix(bres[, c("b1", 
                                                   intvarname)] > 0), 1, function(i) all(i == c(b1_pos, 
                                                                                                ifelse(intcont, bint_cont_pos, bint_cat_pos))))
      mean_weight = apply(bres[filter_direction & conv == 
                                 0, mix_name], 2, weighted.mean, signal(bres[filter_direction & 
                                                                               conv == 0, par]))
      if (all(is.nan(mean_weight))) {
        if (dwqs) {
          Sys.sleep(0.01)
          message(paste0("There are no ", ifelse(b1_pos, 
                                                 "positive", "negative"), " b1 and/or bint values in the specified direction(s) among the bootstrapped models.\nRunning the model in a single ", 
                         ifelse(!b1_pos, "positive", "negative"), 
                         " direction...\n"))
          flush.console()
          mc$formula <- as.formula(gsub(ifelse(b1_pos, 
                                               "nwqs", "pwqs"), "wqs", deparse(formula)))
          mc$formula <- remove_terms(mc$formula, ifelse(b1_pos, 
                                                        "pwqs", "nwqs"))
          mc$b1_pos <- !b1_pos
          single_wqs <- eval(mc)
          return(single_wqs)
        }
        else if (solve_dir_issue == "inverse") 
          mean_weight <- 1/apply(bres[, mix_name], 2, 
                                 weighted.mean, signal(bres[, par]))/sum(1/apply(bres[, 
                                                                                      mix_name], 2, weighted.mean, signal(bres[, 
                                                                                                                               par])))
        else if (solve_dir_issue == "average") 
          mean_weight <- rep(1/length(mix_name), length(mix_name))
        else if (!(solve_dir_issue %in% c("average", 
                                          "inverse")) & rh == 1) 
          stop("There are no b1 or bint values in the specified direction(s) among the bootstrapped models\n")
        else if (!(solve_dir_issue %in% c("average", 
                                          "inverse")) & rh > 1) 
          return(NULL)
      }
      return(mean_weight)
    }
    if (dwqs) {
      mean_weight_p <- mean_weight_f(bresp, b1_pos = TRUE, 
                                     par = "tol")
      if (inherits(mean_weight_p, "gwqs")) 
        return(mean_weight_p)
      mean_weight_n <- mean_weight_f(bresn, b1_pos = FALSE, 
                                     par = "tol")
      if (inherits(mean_weight_n, "gwqs")) 
        return(mean_weight_n)
      mean_weight <- c(mean_weight_p, mean_weight_n)
    }
    else mean_weight <- mean_weight_f(bres, b1_pos = b1_pos, 
                                      par = "stat")
    if (is.null(mean_weight)) 
      return(NULL)
    if (validation == 0) 
      validation_rows <- lapply(validation_rows, function(i) !i)
    wqs_model = model.fit(mean_weight, dtf[validation_rows[[i]], 
    ], Q[validation_rows[[i]], ], Y[validation_rows[[i]]], 
    family, dwqs, zilink, formula, ff, wghts[validation_rows[[i]]], 
    offset[validation_rows[[i]]], initp, Xnames, stratified, 
    b1_pos, zero_infl, kx, kw, intpos, wqsint, intvars)
    if (dwqs) {
      data_plot <- data.frame(mix_name, mean_weight_p, 
                              mean_weight_n, stringsAsFactors = TRUE)
      data_plot = data_plot[order(data_plot$mean_weight_p, 
                                  decreasing = TRUE), ]
      pwqs_index = as.numeric(unlist(wqs_model$pwqs))
      nwqs_index = as.numeric(unlist(wqs_model$nwqs))
      dtf$pwqs <- as.numeric(Q %*% mean_weight_p)
      dtf$nwqs <- as.numeric(Q %*% mean_weight_n)
      dtf$wqs <- wqs_index <- NULL
    }
    else {
      data_plot <- data.frame(mix_name, mean_weight, stringsAsFactors = TRUE)
      data_plot <- data_plot[order(data_plot$mean_weight, 
                                   decreasing = TRUE), ]
      wqs_index <- as.numeric(unlist(wqs_model$wqs))
      dtf$wqs <- as.numeric(Q %*% mean_weight)
      form_terms <- attr(terms.formula(formula), "term.labels")
      form_terms <- gsub("^`|`$", "", form_terms)
      wqsint <- any(grepl("wqs:", form_terms))
      intwqs <- any(grepl(":wqs", form_terms))
      intpos <- ifelse(wqsint, which(grepl("wqs:", form_terms)), 
                       ifelse(intwqs, which(grepl(":wqs", form_terms)), 
                              NA))
      if (wqsint | intwqs) {
        int_vars <- unlist(strsplit(form_terms[intpos], 
                                    ":"))
        intvar <- int_vars[int_vars != "wqs"]
      }
      dtf$pwqs <- dtf$nwqs <- pwqs_index <- nwqs_index <- NULL
    }
    if (rh == 1) {
      if (all(attr(terms(formula), "term.labels") %in% 
              "wqs") | all(attr(terms(formula), "term.labels") %in% 
                           c("pwqs", "nwqs"))) 
        y_plot <- model.response(model.frame(formula, 
                                             dtf[validation_rows[[i]], ]), "any")
      else {
        if (dwqs) {
          formula_wo_wqs <- remove_terms(formula, "pwqs")
          formula_wo_wqs <- remove_terms(formula_wo_wqs, 
                                         "nwqs")
        }
        else formula_wo_wqs <- remove_terms(formula, 
                                            "wqs")
        if (zero_infl) {
          if (length(ff[[3]]) > 1 && identical(ff[[3]][[1]], 
                                               as.name("|"))) {
            if (all(grepl("wqs", attr(terms(as.formula(paste0(ff[[2]], 
                                                              " ~ ", deparse(ff[[3]][[2]])))), "term.labels")))) 
              f1 <- as.formula(paste0(ff[[2]], " ~ ", 
                                      1))
            else f1 <- remove_terms(as.formula(paste0(ff[[2]], 
                                                      " ~ ", deparse(ff[[3]][[2]]))), "wqs")
            if (all(grepl("wqs", attr(terms(as.formula(paste0(ff[[2]], 
                                                              " ~ ", deparse(ff[[3]][[3]])))), "term.labels")))) 
              f2 <- as.formula(paste0(ff[[2]], " ~ ", 
                                      1))
            else f2 <- remove_terms(as.formula(paste0(ff[[2]], 
                                                      " ~ ", deparse(ff[[3]][[3]]))), "wqs")
            formula_wo_wqs <- as.formula(paste0(deparse(f1), 
                                                " | ", f2[[3]]))
          }
          fit <- zeroinfl(formula_wo_wqs, dtf[validation_rows[[i]], 
          ], dist = family$family, link = zilink$name)
        }
        else {
          if (family$family == "negbin") 
            fit = glm.nb(formula_wo_wqs, dtf[validation_rows[[i]], 
            ])
          else fit = glm(formula_wo_wqs, dtf[validation_rows[[i]], 
          ], family = family)
        }
        if (family$family == "binomial") 
          y_plot = fit$y
        else y_plot = mean(fit$y) + resid(fit, type = "pearson")
      }
      if (dwqs) {
        y_adj_wqs_df = as.data.frame(cbind(y_plot, wqs_model$pwqs, 
                                           wqs_model$nwqs))
        names(y_adj_wqs_df) = c(ifelse(family$family == 
                                         "binomial", "y", "y_adj"), "pwqs", "nwqs")
        y_adj_wqs_df <- melt(y_adj_wqs_df, measure.vars = c("pwqs", 
                                                            "nwqs"))
      }
      else {
        y_adj_wqs_df <- as.data.frame(cbind(y_plot, 
                                            wqs_model$wqs))
        names(y_adj_wqs_df) <- c(ifelse(family$family == 
                                          "binomial", "y", "y_adj"), "wqs")
      }
    }
    results = list(fit = wqs_model$m_f, final_weights = data_plot, 
                   conv = conv, bres = bres, wqs = wqs_index, pwqs = pwqs_index, 
                   nwqs = nwqs_index, qi = qi, bindex = bindex, slctd_vars = slctd_vars, 
                   validation_rows = validation_rows[[i]], vindex = validation_rows[[i]], 
                   y_wqs_df = if (rh == 1) y_adj_wqs_df else NULL, 
                   family = family, call = cl, formula = formula, mix_name = mix_name, 
                   stratified = stratified, q = q, zero_infl = zero_infl, 
                   zilink = zilink, dwqs = dwqs, data = dtf, objfn_values = val, 
                   optim_messages = mex, rh = rh)
    if (zero_infl) 
      results$formula <- ff
    class(results) <- "gwqs"
    return(results)
  }, future.seed = fseed)
  if (rh == 1) 
    return(gwqslist[[1]])
  gwqslist <- gwqslist[!vapply(gwqslist, is.null, logical(1))]
  if (length(gwqslist) == 0) 
    stop("The model never converged. Try to increase the number of repeated holdout 'rh' or to change the direction of the association.\n")
  if (length(gwqslist) == 1) 
    stop("The model converged for only a single iteration. Try to increase the number of repeated holdout 'rh'\n")
  if (control$trace) 
    cat(paste0("The model converged ", length(gwqslist), 
               " times out of ", rh, " repeated holdout iterations.\n"))
  coeflist <- lapply(gwqslist, function(i) {
    if (zero_infl) 
      i$fit$coefficients$count
    else i$fit$coefficients
  })
  coefmat <- do.call("rbind", coeflist)
  coefmean <- colMeans(coefmat)
  coefsd <- apply(coefmat, 2, sd)
  coefCInorm <- cbind(coefmean - 1.96 * coefsd, coefmean + 
                        1.96 * coefsd)
  coefmedian <- apply(coefmat, 2, median)
  coefCIperc <- apply(coefmat, 2, function(i) quantile(i, 
                                                       probs = c(0.025, 0.975)))
  coefest <- cbind(coefmean, coefsd, coefCInorm, coefmedian, 
                   t(coefCIperc))
  colnames(coefest) <- c("Mean Est.", "Std. Error", "norm 2.5 %", 
                         "norm 97.5 %", "Median Est.", "perc 2.5 %", "perc 97.5 %")
  if (zero_infl) {
    countcoefmat <- coefmat
    countcoefest <- coefest
    zerocoeflist <- lapply(gwqslist, function(i) i$fit$coefficients$zero)
    zerocoefmat <- do.call("rbind", zerocoeflist)
    zerocoefmean <- colMeans(zerocoefmat)
    zerocoefsd <- apply(zerocoefmat, 2, sd)
    zerocoefCInorm <- cbind(zerocoefmean - 1.96 * zerocoefsd, 
                            zerocoefmean + 1.96 * zerocoefsd)
    zerocoefmedian <- apply(zerocoefmat, 2, median)
    zerocoefCIperc <- apply(zerocoefmat, 2, function(i) quantile(i, 
                                                                 probs = c(0.025, 0.975)))
    zerocoefest <- cbind(zerocoefmean, zerocoefsd, zerocoefCInorm, 
                         zerocoefmedian, t(zerocoefCIperc))
    colnames(zerocoefest) <- c("Mean Est.", "Std. Error", 
                               "norm 2.5 %", "norm 97.5 %", "Median Est.", "perc 2.5 %", 
                               "perc 97.5 %")
    coefest <- list(countcoefest = countcoefest, zerocoefest = zerocoefest)
    coefmat <- list(countcoefmat = countcoefmat, zerocoefmat = zerocoefmat)
  }
  wl <- lapply(gwqslist, function(i) i$final_weights)
  if (dwqs) {
    for (i in 1:length(wl)) {
      if (i == 1) {
        wmatpos <- wl[[1]][, 1:2]
        wmatneg <- wl[[1]][, c(1, 3)]
      }
      else {
        wmatpos <- suppressWarnings(merge(wmatpos, wl[[i]][, 
                                                           1:2], by = "mix_name"))
        wmatneg <- suppressWarnings(merge(wmatneg, wl[[i]][, 
                                                           c(1, 3)], by = "mix_name"))
      }
    }
    nameswmat <- as.character(wmatpos[, 1])
    wmatpos <- t(wmatpos[, -1])
    wmatneg <- t(wmatneg[, -1])
    colnames(wmatpos) <- colnames(wmatneg) <- nameswmat
    rownames(wmatpos) <- rownames(wmatneg) <- NULL
    wmeanpos <- colMeans(wmatpos)
    wmeanneg <- colMeans(wmatneg)
    wCIpos <- apply(wmatpos, 2, function(i) quantile(i, 
                                                     probs = c(0.025, 0.975)))
    wCIneg <- apply(wmatneg, 2, function(i) quantile(i, 
                                                     probs = c(0.025, 0.975)))
    final_weights_pos <- as.data.frame(cbind(wmeanpos, t(wCIpos)))
    final_weights_neg <- as.data.frame(cbind(wmeanneg, t(wCIneg)))
    names(final_weights_pos) <- c("Estimate pos", "2.5% pos", 
                                  "97.5% pos")
    names(final_weights_neg) <- c("Estimate neg", "2.5% neg", 
                                  "97.5% neg")
    final_weights <- cbind(final_weights_pos, final_weights_neg)
    wmat <- list(wmatpos = wmatpos, wmatneg = wmatneg)
  }
  else {
    for (i in 1:length(wl)) {
      if (i == 1) 
        wmat <- wl[[1]]
      else wmat <- suppressWarnings(merge(wmat, wl[[i]], 
                                          by = "mix_name"))
    }
    nameswmat <- as.character(wmat[, 1])
    wmat <- t(wmat[, -1])
    colnames(wmat) <- nameswmat
    rownames(wmat) <- NULL
    wmean <- colMeans(wmat)
    wCI <- apply(wmat, 2, function(i) quantile(i, probs = c(0.025, 
                                                            0.975)))
    final_weights <- as.data.frame(cbind(wmean, t(wCI)))
    names(final_weights) <- c("Estimate", "2.5%", "97.5%")
  }
  final_weights <- cbind(rownames(final_weights), final_weights)
  names(final_weights)[1] <- "mix_name"
  if (dwqs) {
    dtf$pwqs <- pwqs <- as.numeric(Q %*% final_weights[match(mix_name, 
                                                             final_weights[, 1]), 2])
    dtf$nwqs <- nwqs <- as.numeric(Q %*% final_weights[match(mix_name, 
                                                             final_weights[, 1]), 5])
    dtf$wqs <- wqs <- NULL
  }
  else {
    dtf$wqs <- wqs <- as.numeric(Q %*% final_weights[match(mix_name, 
                                                           final_weights[, 1]), 2])
    dtf$pwqs <- pwqs <- dtf$nwqs <- nwqs <- NULL
  }
  if (all(grepl("wqs", attr(terms(formula), "term.labels")))) 
    y_plot <- model.response(model.frame(formula, dtf), 
                             "any")
  else {
    formula_wo_wqs = remove_terms(formula, "wqs")
    if (family$family != "multinomial") {
      if (zero_infl) {
        if (length(ff[[3]]) > 1 && identical(ff[[3]][[1]], 
                                             as.name("|"))) {
          if (all(grepl("wqs", attr(terms(as.formula(paste0(ff[[2]], 
                                                            " ~ ", deparse(ff[[3]][[2]])))), "term.labels")))) 
            f1 <- as.formula(paste0(ff[[2]], " ~ ", 
                                    1))
          else f1 <- remove_terms(as.formula(paste0(ff[[2]], 
                                                    " ~ ", deparse(ff[[3]][[2]]))), "wqs")
          if (all(grepl("wqs", attr(terms(as.formula(paste0(ff[[2]], 
                                                            " ~ ", deparse(ff[[3]][[3]])))), "term.labels")))) 
            f2 <- as.formula(paste0(ff[[2]], " ~ ", 
                                    1))
          else f2 <- remove_terms(as.formula(paste0(ff[[2]], 
                                                    " ~ ", deparse(ff[[3]][[3]]))), "wqs")
          formula_wo_wqs <- as.formula(paste0(deparse(f1), 
                                              " | ", f2[[3]]))
        }
        fit <- zeroinfl(formula_wo_wqs, dtf, dist = family$family, 
                        link = zilink$name)
      }
      else {
        if (family$family == "negbin") 
          fit = glm.nb(formula_wo_wqs, dtf)
        else fit = glm(formula_wo_wqs, dtf, family = family)
      }
      if (family$family == "binomial") 
        y_plot = fit$y
      else y_plot = mean(fit$y) + resid(fit, type = "pearson")
    }
  }
  final_weights <- final_weights[order(final_weights[, 2], 
                                       decreasing = TRUE), ]
  if (dwqs) {
    y_adj_wqs_df <- as.data.frame(cbind(y_plot, pwqs, nwqs))
    names(y_adj_wqs_df) <- c(ifelse(family$family == "binomial", 
                                    "y", "y_adj"), "pwqs", "nwqs")
  }
  else {
    y_adj_wqs_df <- as.data.frame(cbind(y_plot, wqs))
    names(y_adj_wqs_df) <- c(ifelse(family$family == "binomial", 
                                    "y", "y_adj"), "wqs")
  }
  fit <- list(coefficients = NULL, aic = NULL, dispersion = NULL, 
              deviance = NULL, df.residual = NULL, null.deviance = NULL, 
              df.null = NULL, theta = NULL, SE.theta = NULL)
  fit$coefficients <- coefest
  if (zero_infl) 
    fit$aic <- mean(2 * (nrow(coefest[[1]]) + nrow(coefest[[2]]) + 
                           ifelse(family$family == "negbin", 1, 0)) - 2 * sapply(gwqslist, 
                                                                                 function(i) i$fit$loglik), na.rm = TRUE)
  else {
    fit$dispersion <- mean(sapply(gwqslist, function(i) summary(i$fit)$dispersion), 
                           na.rm = TRUE)
    fit$deviance <- mean(sapply(gwqslist, function(i) i$fit$deviance), 
                         na.rm = TRUE)
    fit$df.residual <- gwqslist[[1]]$fit$df.residual
    fit$null.deviance <- mean(sapply(gwqslist, function(i) i$fit$null.deviance), 
                              na.rm = TRUE)
    fit$df.null <- gwqslist[[1]]$fit$df.null
    fit$aic <- mean(sapply(gwqslist, function(i) i$fit$aic), 
                    na.rm = TRUE)
  }
  if (family$family == "negbin") {
    fit$theta <- mean(sapply(gwqslist, function(i) i$fit$theta), 
                      na.rm = TRUE)
    fit$SE.theta <- sd(sapply(gwqslist, function(i) i$fit$theta), 
                       na.rm = TRUE)
  }
  results <- list(fit = fit, final_weights = final_weights, 
                  wqs = wqs, qi = qi, y_wqs_df = y_adj_wqs_df, family = family, 
                  call = cl, formula = formula, mix_name = mix_name, stratified = stratified, 
                  q = q, zero_infl = zero_infl, zilink = zilink, dwqs = dwqs, 
                  data = dtf, gwqslist = gwqslist, coefmat = coefmat, 
                  wmat = wmat, rh = rh)
  if (zero_infl) 
    results$formula <- ff
  class(results) <- "gwqs"
  return(results)
}

#reproducing gWQS:::optim.f below to get gwqs_hpc to run
optim.f <- function (i, objfn, Y, Xm, Q, offset, wghts, initp, b1_pos, 
                     b_constr, n_vars, dwqs, family, rs, zilink, zero_infl, formula, 
                     ff, kx, kw, Xnames, stratified, b, optim.method, control, 
                     lambda, stratlev, intpos, wqsint, intvars, bint_cont_pos, 
                     bint_cat_pos, intcont){
  if (rs) {
    bindex <- 1:nrow(Xm)
    slctd_vars <- sample(colnames(Q), n_vars, replace = FALSE)
    if (dwqs) 
      initp <- initp[c(rep(T, ncol(Xm)), rep(colnames(Q) %in% 
                                               slctd_vars, 2))]
    else initp <- initp[c(rep(T, ncol(Xm)), colnames(Q) %in% 
                            slctd_vars)]
    kw <- length(slctd_vars)
  }
  else {
    if (b == 1) 
      bindex <- 1:nrow(Xm)
    else bindex <- sample(1:nrow(Xm), nrow(Xm), replace = TRUE)
    slctd_vars <- colnames(Q)
  }
  bXm <- Xm[bindex, ]
  bQ <- Q[bindex, slctd_vars]
  bY <- Y[bindex]
  bwghts <- wghts[bindex]
  boffset <- offset[bindex]
  opt_res <- optim(par = initp, fn = objfn, method = optim.method, 
                   control = control, kw = kw, bXm = bXm, bY = bY, boffset = boffset, 
                   bQ = bQ, kx = kx, Xnames = Xnames, family = family, 
                   dwqs = dwqs, zilink = zilink, zero_infl = zero_infl, 
                   formula = formula, ff = ff, bwghts = bwghts, stratified = stratified, 
                   stratlev = stratlev, b1_pos = b1_pos, b_constr = b_constr, 
                   lambda = lambda, intpos = intpos, wqsint = wqsint, intvars = intvars, 
                   bint_cont_pos = bint_cont_pos, bint_cat_pos = bint_cat_pos, 
                   intcont = intcont)
  if (!is.null(opt_res)) {
    if (dwqs) 
      par_opt <- c(opt_res$par[(kx + 1):(kx + kw)]^2/sum(opt_res$par[(kx + 
                                                                        1):(kx + kw)]^2), opt_res$par[(kx + kw + 1):(kx + 
                                                                                                                       2 * kw)]^2/sum(opt_res$par[(kx + kw + 1):(kx + 
                                                                                                                                                                   2 * kw)]^2))
    else {
      if (!is.null(stratified)) {
        par_opt <- lapply(1:stratlev, function(i) {
          ncolQ <- ncol(bQ)/stratlev
          tmp <- (opt_res$par[(kx + 1 + (i - 1) * ncolQ):(kx + 
                                                            kw/stratlev * i)]^2)/sum(opt_res$par[(kx + 
                                                                                                    1 + (i - 1) * ncolQ):(kx + kw/stratlev * 
                                                                                                                            i)]^2)
          return(tmp)
        })
        par_opt <- do.call("c", par_opt)
      }
      else par_opt <- (opt_res$par[(kx + 1):(kx + kw)]^2)/sum(opt_res$par[(kx + 
                                                                             1):(kx + kw)]^2)
    }
    conv <- opt_res$convergence
    counts <- opt_res$counts
    val <- opt_res$val
    mex <- opt_res$message
  }
  else {
    if (dwqs) 
      par_opt <- c(initp[(kx + 1):(kx + kw)]^2/sum(initp[(kx + 
                                                            1):(kx + kw)]^2), initp[(kx + kw + 1):(kx + 
                                                                                                     2 * kw)]^2/sum(initp[(kx + kw + 1):(kx + 2 * 
                                                                                                                                           kw)]^2))
    else par_opt <- (initp[(kx + 1):(kx + kw)]^2)/sum(initp[(kx + 
                                                               1):(kx + kw)]^2)
    conv <- 1
    counts <- val <- mex <- NULL
  }
  if (any(is.infinite(par_opt))) {
    if (dwqs) {
      par_opt[is.infinite(par_opt[1:kw])] <- 1/sum(is.infinite(par_opt[1:kw]))
      par_opt[!(is.infinite(par_opt[1:kw]))] <- 0
      par_opt[(kw + 1):2 * kw][is.infinite(par_opt[(kw + 
                                                      1):2 * kw])] <- 1/sum(is.infinite(par_opt[(kw + 
                                                                                                   1):2 * kw]))
      par_opt[(kw + 1):2 * kw][!(is.infinite(par_opt[(kw + 
                                                        1):2 * kw]))] <- 0
    }
    else {
      par_opt[is.infinite(par_opt)] <- 1/sum(is.infinite(par_opt))
      par_opt[!(is.infinite(par_opt))] <- 0
    }
  }
  mfit <- model.fit(w = par_opt, dt = bXm, bQ = bQ, Y = bY, 
                    family = family, dwqs = dwqs, zilink = zilink, formula = formula, 
                    ff = ff, wghts = bwghts, offset = boffset, initp = initp, 
                    Xnames = Xnames, stratified = stratified, b1_pos = b1_pos, 
                    zero_infl = zero_infl, kx = kx, kw = kw, intpos = intpos, 
                    wqsint = wqsint, intvars = intvars)
  out <- list(par_opt = par_opt, conv = conv, counts = counts, 
              val = val, mex = mex, mfit = mfit, bindex = bindex, 
              slctd_vars = slctd_vars)
  return(out)
}

#reproducing gWQS:::model.fit below to get gwqs_hpc to run
model.fit <- function (w, dt, bQ, Y, family, dwqs, zilink, formula, ff, wghts, 
                       offset, initp, Xnames, stratified, b1_pos, zero_infl, kx, 
                       kw, intpos, wqsint, intvars){
  if (is.matrix(dt)) {
    dtf <- as.data.frame(dt[, which(colnames(dt) != "(Intercept)")])
    formula <- ff <- Y ~ .
    zero_infl <- FALSE
    intpos <- grepl("wqs", colnames(dtf)) & grepl(":", colnames(dtf))
  }
  else dtf <- dt
  if (dwqs) {
    pwqs <- dtf[, "pwqs"] <- as.numeric(bQ %*% w[1:kw])
    nwqs <- dtf[, "nwqs"] <- as.numeric(bQ %*% w[(kw + 1):(2 * 
                                                             kw)])
    wqs <- NULL
  }
  else {
    wqs <- dtf[, "wqs"] <- as.numeric(bQ %*% w)
    if (wqsint & is.matrix(dt)) 
      dtf[, intpos] <- wqs * dtf[, intvars]
    pwqs <- nwqs <- NULL
  }
  if (zero_infl) 
    m_f <- zeroinfl(ff, dtf, dist = family$family, link = zilink$name)
  else {
    if (family$family == "negbin") 
      m_f = glm.nb(formula, data = dtf, weights = wghts)
    else m_f = glm(formula, data = dtf, family = family, 
                   weights = wghts)
  }
  mf_out = list(wqs = wqs, pwqs = pwqs, nwqs = nwqs, m_f = m_f)
  return(mf_out)
}

#reproducing gWQS:::parnames below to get gwqs_hpc to run
parnames <- function (df, formula, form2){
  if (!is.null(form2)) {
    mf <- model.frame(form2, df)
    Y <- model.response(mf, "any")
    df$yz <- ifelse(Y == 0, 0, 1)
  }
  mm <- model.matrix(formula, df)
  colnames(mm)
}

#reproducing gWQS:::remove_terms below to get gwqs_hpc to run
remove_terms <- function (form, term){
  fterms <- terms(form)
  fac <- attr(fterms, "factors")
  if (term %in% rownames(fac)) {
    fac_wqs <- fac[grep(term, rownames(fac)), ]
    if (NCOL(fac_wqs) == 1) 
      idx <- which(as.logical(fac[term, ]))
    else idx <- which(apply(fac_wqs, 2, function(i) any(i == 
                                                          1)))
    new_fterms <- drop.terms(fterms, dropx = idx, keep.response = TRUE)
    return(formula(new_fterms))
  }
  else return(form)
}

#reproducing gWQS:::set_par_names below to get gwqs_hpc to run
set_par_names <- function (i, slctd_vars, param, q_name, family, dwqs){
  temp <- param[[i]]$par_opt
  if (family$family == "multinomial") {
    param[[i]]$par_opt <- matrix(NA, length(q_name), dim(temp)[2])
    param[[i]]$par_opt[which(q_name %in% slctd_vars[[i]]), 
    ] <- temp
    rownames(param[[i]]$par_opt) <- q_name
  }
  else {
    param[[i]]$par_opt <- rep(NA, ifelse(dwqs, 2, 1) * length(q_name))
    names(param[[i]]$par_opt) <- rep(q_name, ifelse(dwqs, 
                                                    2, 1))
    if (dwqs) {
      param[[i]]$par_opt[1:length(q_name)][slctd_vars[[i]]] <- temp[1:(length(temp)/2)]
      param[[i]]$par_opt[(length(q_name) + 1):(2 * length(q_name))][slctd_vars[[i]]] <- temp[(length(temp)/2 + 
                                                                                                1):length(temp)]
    }
    else param[[i]]$par_opt[slctd_vars[[i]]] <- temp
  }
  return(param[[i]])
}

#reproducing gWQS:::stratified_f below to get gwqs_hpc to run
stratified_f <- function (Q, dtf, stratified, mix_name){
  ls = levels(unlist(dtf[, stratified, drop = FALSE]))
  if (is.null(ls)) 
    stop("'stratified' must be factor\n")
  llsm = lapply(ls, function(ls) {
    mat = diag(as.numeric(dtf[, stratified] == ls))
    sub_dt = mat %*% Q[, mix_name]
    colnames(sub_dt) = paste(mix_name, ls, sep = "_")
    return(sub_dt)
  })
  Q = do.call("cbind", llsm)
  mix_name = colnames(Q)
  strtfd_out = list(Q, mix_name)
  names(strtfd_out) = c("Q", "mix_name")
  return(strtfd_out)
}

#reproducing gWQS:::values.0 below to get gwqs_hpc to run
values.0 <- function (kw, Xnames, kx, formula, ff, wghts, bdtf, stratified, 
                      stratlev, b1_pos, family, wp, wn, dwqs, zilink, zero_infl){
  w = rep(1/kw * stratlev, kw)
  w <- sqrt(w)
  if (family$family == "multinomial") {
    n_levels <- nlevels(eval(formula[[2]], envir = bdtf))
    fit = multinom(formula, bdtf, trace = F, weights = wghts)
    bj = c(sapply(1:(n_levels - 1), function(i) coef(fit)[i, 
    ]))
    names(bj) <- Xnames
    if (dwqs) 
      val.0 = c(bj, sqrt(wp), sqrt(wn))
    else {
      w = rep(w, (n_levels - 1))
      val.0 = c(bj, w)
    }
  }
  else {
    if (family$family == "negbin") 
      fit = suppressWarnings(glm.nb(formula, bdtf, weights = wghts))
    else fit = glm(formula, bdtf, family = family, weights = wghts)
    bj = coef(fit)
    if (dwqs) 
      val.0 = c(bj, sqrt(wp), sqrt(wn))
    else val.0 = c(bj, w)
    if (family$family == "negbin") {
      if (length(attr(terms(formula), "term.labels")) == 
          1) 
        val.0 = c(val.0, 1)
      else val.0 = c(val.0, log(fit$theta))
    }
  }
  if (dwqs) {
    pwqs_site <- which(grepl("pwqs", Xnames))
    nwqs_site <- which(grepl("nwqs", Xnames))
    val.0[pwqs_site] = 1e-04
    val.0[nwqs_site] = -1e-04
  }
  else {
    wqs_site <- which(grepl("wqs", Xnames))
    val.0[wqs_site] = sapply(b1_pos, function(i) ifelse(i, 
                                                        1e-04, -1e-04))
  }
  return(val.0)
}

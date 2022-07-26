#' WQS permutation test
#' 
#' \code{wqs_pt} takes a `gwqs` object as an input and runs the permutation 
#' test (Day et al. 2022) to obtain an estimate for the p-value significance for 
#' the WQS coefficient.  
#' 
#' To use `wqs_pt`, we first need to run an initial WQS regression run while 
#' setting `validation = 0`. We will use this `gwqs` object as the model argument 
#' for the `wqs_pt` function. Note that permutation test has so far only been
#' validated for linear WQS regression (i.e., `family = "gaussian"`) or logistic
#' WQS regression (i.e., `family = binomial(link = "logit")`), though the
#' permutation test algorithm should also work for all WQS GLMs. Therefore,
#' this function accepts `gwqs` objects made with the following families: 
#' "gaussian" or gaussian(link = "identity"), "binomial" or binomial() with 
#' any accepted link function (e.g., "logit" or "probit"), "poisson" or 
#' poisson(link="log"), "negbin" for negative binomial, and "quasipoisson" or
#' quasipoisson(link="log"). This function cannot currently accomodate `gwqs`
#' objects made with the "multinomial" family, and it is not currently able to 
#' accomodate stratified weights or WQS interaction terms (e.g., `y ~ wqs * sex`).
#' 
#' The argument `boots` is the number of bootstraps for the WQS regression run 
#' in each permutation test iteration. Note that we may elect a bootstrap count 
#' `boots` lower than that specified in the model object for the sake of 
#' efficiency. If `boots` is not specified, then we will use the same bootstrap 
#' count in the permutation test WQS regression runs as that specified in the 
#' model argument.
#'
#' The arguments `b1_pos` and `rs` should be consistent with the inputs chosen 
#' in the model object. The seed should ideally be consistent with the seed set 
#' in the model object for consistency, though this is not required.
#'
#' @param model A \code{gwqs} object as generated from the \code{gWQS} package.  
#' @param niter Number of permutation test iterations. 
#' @param boots Number of bootstrap samples for each permutation test WQS 
#' regression iteration. If `boots` is not specified, then we will use the same 
#' bootstrap count for each permutation test WQS regression iteration as that 
#' specified in the main WQS regression run.
#' @param b1_pos A logical value that indicates whether beta values should be 
#' positive or negative.
#' @param b1_constr Logical value that determines whether to apply positive or 
#' negative constraints in the optimization function for the weight optimization.
#' @param rs A logical value indicating whether random subset implementation 
#' should be performed. 
#' @param plan_strategy Evaluation strategy for the plan function. You can choose 
#' among "sequential", "transparent", "multisession", "multicore", 
#' "multiprocess", "cluster" and "remote." See future::plan documentation for full 
#' details. 
#' @param seed (optional) Random seed for the permutation test WQS reference run. 
#' This should be the same random seed as used for the main WQS regression run. 
#' This seed will be saved in the "gwqs_perm" object as "gwqs_perm$seed".
#' 
#' @return \code{wqs_pt} returns an object of class `wqs_pt`, which contains: 
#' 
#' \item{perm_test}{List containing: (1) `pval`: permutation test p-value, 
#' (2) (linear WQS regression only) `testbeta1`: reference WQS coefficient beta1 value, 
#' (3) (linear WQS regression only) `betas`: Vector of beta values from 
#' each permutation test run, (4) (WQS GLM only) `testpval`: test reference 
#' p-value, (5) (WQS GLM only) `permpvals`: p-values from the null models.}
#' \item{gwqs_main}{Main gWQS object (same as model input).}
#' \item{gwqs_perm}{Permutation test reference gWQS object (NULL if model 
#' `family != "gaussian"` or if same number of bootstraps are used in permutation 
#' test WQS regression runs as in the main run).}
#' @import gWQS ggplot2 viridis cowplot stats methods
#' @export wqs_pt
#' 
#' @examples
#' library(gWQS)
#' 
#' # mixture names
#' PCBs <- names(wqs_data)[1:17] #half of the original 34 for quick computation
#' 
#' # create reference wqs object with 5 bootstraps
#' wqs_main <- gwqs(yLBX ~ wqs, mix_name = PCBs, data = wqs_data, q = 10, 
#'                  validation = 0, b = 5, b1_pos = TRUE, b1_constr = FALSE,
#'                  plan_strategy = "multicore", family = "gaussian", seed = 16)
#' # Note: We recommend niter = 1000 for the main WQS regression. This example
#' # has a lower number of bootstraps to serve as a shorter test run.
#' 
#' # run permutation test
#' perm_test_res <- wqs_pt(wqs_main, niter = 4, b1_pos = TRUE)
#' 
#' # Note: The default value of niter = 200 is the recommended parameter values. 
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
                     b1_constr = FALSE, rs = FALSE, plan_strategy = "multicore", 
                     seed = NULL) {
  
  pbapply::pboptions(type="timer")
  
  if (is(model, "gwqs")) {
    if (model$family$family == "multinomial"){
      stop("The permutation test is not currently set up to accomodate the 
           multinomial WQS regressions.")
    }
  } else stop("'model' must be of class 'gwqs' (see gWQS package).")
  
  mm <- model$fit
  formchar <- as.character(formula(mm))

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
  yname <- as.character(formula(mm))[2]
  mix_name <- names(model$bres)[names(model$bres) %in% model$final_weights$mix_name]
  
  if (!is.null(model$qi)) {
    nq <- max(sapply(model$qi, length)) - 1
  } else {
    # this is for cases when there is no quantile transformation or it's already been
    # done in the data frame
    nq <- NULL
  }
  
  if (is.null(boots)){
    boots <- length(model$bindex)
  }
  
  if (model$family$family == "gaussian"){
    Data <- model$data[, -which(names(model$data) %in% c("wqs", "wghts"))]
    
    # reference WQS regression run 
    if (boots == length(model$bindex)){
      perm_ref_wqs <- model
      ref_beta1 <- mm$coef[2]
    }
    
    else{
      perm_ref_wqs <- gwqs(formula = formula(mm), data = Data, mix_name = mix_name, 
                           q = nq, b = boots, rs = rs, validation = 0, 
                           plan_strategy = plan_strategy, b1_pos = b1_pos, 
                           b1_constr = b1_constr, seed = seed)
      
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
      reorgmat <- matrix(NA, dim(Data)[1], niter)
      reorgmat <- apply(reorgmat, 2, function(x) sample(Data[, yname]))
    }
    
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
        suppressWarnings(gwqs(formula = form1, data = newDat, mix_name = mix_name, 
                              q = nq, b = boots, rs = rs, validation = 0, 
                              plan_strategy = plan_strategy, b1_pos = b1_pos, 
                              b1_constr = b1_constr))
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
    
    if (any(is.na(betas))) {
      print(paste0(length(which(is.na(betas))), " failed model attempts"))
    }
    
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
  } 
  
  else {
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
      suppressWarnings(gwqs(formula = formula(mm), data = Data, 
                            mix_name = model$mix_name, q = nq, b = boots, 
                            rs = rs, validation = 0, 
                            plan_strategy = plan_strategy, b1_pos = b1_pos, 
                            family = model$family, seed = seed,
                            b1_constr = b1_constr))
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
        suppressWarnings(gwqs(formula = form1, data = newDat, 
                              mix_name = mix_name, q = model$q, b = boots, 
                              rs = rs, validation = 0, 
                              plan_strategy = plan_strategy, 
                              b1_pos = b1_pos, family = model$family$family,
                              b1_constr = b1_constr)
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
    if (any(is.na(permstats))) {
      print(paste0(length(which(is.na(permstats))), " failed model attempts"))
    }
    
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
  
  cat("Permutation test WQS coefficient p-value: \n", 
      object$perm_test$pval,
      "\n")
  
  main_sum <- summary(object$gwqs_main)
  
  print(main_sum)
  
}

#' Plotting method for wqspt object 
#' 
#' Generates plots to help visualize and summarize WQS permutation test results. 
#'
#' @param wqsptresults An object of class `wqs_pt`.
#' @param FixedPalette If TRUE, the heatmap color key for the mixture weights has 
#' categorical cutoffs with the following categories: <0.1, 0.1 - <0.2, 0.2 - <0.3, 
#' and >= 0.3. If false, the heatmap color key is continuous and dependent on the 
#' weight values.
#' @param InclKey If TRUE, a horizontal heatmap legend is included at the bottom 
#' of the full plot.
#' @param AltMixName Defaults to NULL. If not NULL, these are alternative names 
#' for the mixture components to be displayed on the heatmap y axis.
#' @param AltOutcomeName Defaults to NULL. If not NULL, this is an alternative 
#' name for the outcome to be displayed on the heatmap x axis.
#' @param ViridisPalette  Color palette to be used for the viridisLite 
#' package-based coloring of the heatmap, with possible values from 'A' to 'E'. 
#' Defaults to 'D'.
#' @param StripTextSize Text size for the plot strip labels. Defaults to 14.
#' @param AxisTextSize.Y Text size for the y axis text. Defaults to 12.
#' @param AxisTextSize.X Text size for the x axis text. Defaults to 12.
#' @param LegendTextSize Text text size for the legend text. Defaults to 14.
#' @param PvalLabelSize The geom_text size for the permutation test p-value 
#' label. Defaults to 5.
#' @param HeatMapTextSize The geom_text size for the mixture weight heatmap 
#' labels. Defaults to 5.
#'
#' @return Returns a list with 4 objects. 
#' 
#' \item{FullPlot}{Two plots stacked vertically: (1) Forest plot of the beta WQS 
#' coefficient with the naive confidence intervals as well as the permutation 
#' test p-value (2) A heatmap of the WQS weights for each mixture component.}
#' \item{CoefPlot}{Forest plot of the beta WQS 
#' coefficient with the naive confidence intervals as well as the permutation 
#' test p-value.}
#' \item{WtPlot}{A heatmap of the WQS weights for each mixture component.}
#' \item{WtLegend}{A legend for the weights in the WtPlot heatmap.}
#' 
#' @importFrom rlang .data
#' 
#' @export
#'
wqspt_plot <- function(wqsptresults, FixedPalette = FALSE, InclKey = FALSE, 
                          AltMixName = NULL, AltOutcomeName = NULL, 
                          ViridisPalette = "D", StripTextSize = 14,
                          AxisTextSize.Y = 12, AxisTextSize.X = 12, 
                          LegendTextSize = 14, PvalLabelSize = 5, 
                          HeatMapTextSize = 5) {
  
  wqs_fam <- wqsptresults$family
  if(!is.character(wqs_fam)) wqs_fam <- wqs_fam$family

  thisfit <- wqsptresults$gwqs_main$fit
  b1pos <- wqsptresults$gwqs_main$b1_pos
  if (b1pos)
    thisdir <- "Positive"
  else
    thisdir <- "Negative"
  if (!is.null(AltOutcomeName))
    outname <- AltOutcomeName
  else
    outname <- as.character(attr(thisfit$terms, "variables")[[2]])
  
  if (wqs_fam == "gaussian"){
    pval <- summary(thisfit)$coef["wqs", "Pr(>|t|)"]
  }
  else {
    pval <- summary(thisfit)$coef["wqs", "Pr(>|z|)"]
  }
  
  WQSResults <-
    data.frame(
      Outcome = outname,
      Direction = thisdir,
      Beta = thisfit$coef['wqs'],
      LCI = suppressMessages(confint(thisfit)[2, 1]),
      UCI = suppressMessages(confint(thisfit)[2, 2]),
      pval = pval,
      PTp = wqsptresults$perm_test$pval
    )
  WQSResults$PTlabel <- paste0("PTp=", signif(WQSResults$PTp, 3))
  WQSResults$FacetLabel <- "Coefficient"
  cirange <- WQSResults$UCI - WQSResults$LCI
  widercirange <-
    c(WQSResults$LCI - (WQSResults$LCI / 10),
      WQSResults$UCI + (WQSResults$UCI / 10))
  if (widercirange[1] < 0 & widercirange[2] > 0) {
    gg1 <-
      ggplot(WQSResults, aes(x = .data$Outcome, y = .data$Beta)) + geom_point(size = 3) +
      theme_bw() +
      geom_errorbar(aes(ymin = .data$LCI, ymax = .data$UCI),
                    size = 1,
                    width = 0.75) +
      geom_hline(yintercept = 0) +
      geom_text(aes(label = .data$PTlabel, y = .data$UCI + cirange / 10), 
                size = PvalLabelSize) +
      facet_grid(FacetLabel ~ Direction) +
      theme(
        strip.text = element_text(size = StripTextSize),
        axis.text.y = element_text(size = AxisTextSize.Y),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank()
      )
  } else {
    gg1 <-
      ggplot(WQSResults, aes(x = .data$Outcome, y = .data$Beta)) + geom_point(size = 3) +
      theme_bw() +
      geom_errorbar(aes(ymin = .data$LCI, ymax = .data$UCI),
                    size = 1,
                    width = 0.75) +
      geom_text(aes(label = .data$PTlabel, y = .data$UCI + cirange / 10), 
                size = PvalLabelSize) +
      facet_grid(FacetLabel ~ Direction) +
      theme(
        strip.text = element_text(size = StripTextSize),
        axis.text.y = element_text(size = AxisTextSize.Y),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
  
  WQSwts <-
    wqsptresults$gwqs_main$final_weights[wqsptresults$gwqs_main$mix_name, ]
  WQSwts$FacetLabel <- "Weights"
  WQSwts$Outcome <- WQSResults$Outcome
  WQSwts$Direction <- WQSResults$Direction
  WQSwts$mix_name <-
    factor(as.character(WQSwts$mix_name), 
           levels = wqsptresults$gwqs_main$mix_name)
  if (!is.null(AltMixName))
    levels(WQSwts$mix_name) <- AltMixName
  WQSwts$mix_name <-
    factor(WQSwts$mix_name, levels = rev(levels(WQSwts$mix_name)))
  names(WQSwts)[1:2] <- c("Exposure", "Weight")
  if (FixedPalette) {
    mypal <- viridis::viridis_pal(option = ViridisPalette)(4)
    WQSwts$Wt <- WQSwts$Weight
    WQSwts$Weight <- factor(
      ifelse(
        WQSwts$Wt < 0.1,
        "<0.1",
        ifelse(
          WQSwts$Wt >= 0.1 & WQSwts$Wt < 0.2,
          "0.1-<0.2",
          ifelse(
            WQSwts$Wt >= 0.2 & WQSwts$Wt < 0.3,
            "0.2-<0.3",
            paste0("\u2265", "0.3")
          )
        )
      ),
      levels = c("<0.1", "0.1-<0.2", "0.2-<0.3", paste0("\u2265", "0.3"))
    )
    Virclr <- ifelse(
      WQSwts$Weight == "<0.1",
      mypal[1],
      ifelse(
        WQSwts$Weight == "0.1-0.2",
        mypal[2],
        ifelse(
          WQSwts$Weight == "0.2-0.3",
          mypal[3],
          ifelse(is.na(WQSwts$Weight) == T, "grey50", mypal[4])
        )
      )
    )
    names(Virclr) <- as.character(WQSwts$Weight)
    legplot <-
      ggplot(data.frame(Weight = factor(
        levels(WQSwts$Weight), levels = levels(WQSwts$Weight)
      )),
      aes(x = 1, y = .data$Weight)) +
      geom_tile(aes(fill = .data$Weight)) + scale_fill_manual(values = mypal) +
      theme(
        legend.position = "bottom",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14)
      )
    l1 <- cowplot::get_legend(legplot)
    
    gg2 <- ggplot(WQSwts, aes(x = .data$Outcome, y = .data$Exposure)) + 
      theme_classic() +
      geom_tile(aes(fill = .data$Weight), alpha = 0.7) + 
      geom_text(aes(label =round(.data$Wt, 2)), size = HeatMapTextSize) +
      scale_fill_manual(values = Virclr) +
      facet_grid(FacetLabel ~ Direction) +
      theme(
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = StripTextSize),
        axis.text.x = element_text(size = AxisTextSize.X),
        axis.text.y = element_text(size = AxisTextSize.Y),
        axis.title = element_blank(),
        strip.background.y = element_rect(fill = "grey85", colour = "grey20"),
        legend.position = "bottom",
        legend.title = element_text(size = LegendTextSize, face = "bold"),
        legend.text = element_text(size = LegendTextSize)
      )
  } else {
    gg2 <- ggplot(WQSwts, aes(x = .data$Outcome, y = .data$Exposure)) + 
      theme_classic() +
      geom_tile(aes(fill = .data$Weight), alpha = 0.7) + 
      geom_text(aes(label = round(.data$Weight, 2)), size = HeatMapTextSize) +
      scale_fill_viridis_c(option = ViridisPalette) +
      facet_grid(FacetLabel ~ Direction) +
      theme(
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = StripTextSize),
        axis.text.x = element_text(size = AxisTextSize.X),
        axis.text.y = element_text(size = AxisTextSize.Y),
        axis.title = element_blank(),
        strip.background.y = element_rect(fill = "grey85", colour = "grey20"),
        legend.position = "bottom",
        legend.title = element_text(size = LegendTextSize, face = "bold"),
        legend.text = element_text(size = LegendTextSize),
        legend.key.size = unit(0.4, units = 'in')
      )
    l1 <- cowplot::get_legend(gg2)
  }
  
  if (InclKey) {
    gg2 <- gg2 + theme(legend.position = "none")
    fullplot <-
      cowplot::plot_grid(
        cowplot::plot_grid(
          gg1,
          gg2,
          ncol = 1,
          align = "v",
          rel_heights = c(0.4, 0.6)
        ),
        l1,
        ncol = 1,
        rel_heights = c(1, 0.1)
      )
  } else {
    gg2 <- gg2 + theme(legend.position = "none")
    fullplot <-
      cowplot::plot_grid(
        gg1,
        gg2,
        ncol = 1,
        rel_heights = c(0.4, 0.6),
        align = "v"
      )
  }
  
  return(list(
    FullPlot = fullplot,
    CoefPlot = gg1,
    WtPlot = gg2,
    WtLegend = l1
  ))
  
}

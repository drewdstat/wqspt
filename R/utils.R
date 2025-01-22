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
    suppressWarnings(gwqs(formula = form1, data = newDat, mix_name = mix_name, 
                          q = nq, b = boots, rs = rs, validation = 0, 
                          plan_strategy = plan_strategy, b1_pos = b1_pos, 
                          b_constr = b_constr))
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
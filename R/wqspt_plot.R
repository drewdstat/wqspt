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

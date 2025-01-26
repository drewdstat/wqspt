# wqspt 1.0.2

* This function previously did not accommodate tibbles for the `data` argument 
of `wqs_full_perm`, but now `data` can be a `data.frame` object or a tibble.

* Added a `utils.R` script with a function `get_legend2()` to replace the 
malfunctioning `cowplot::get_legend()` call in the function `wqspt_plot`. In 
`cowplot` v1.1.3, `get_legend()` returns an error if the legend position is 
anywhere besides on the right of the plot. Since `wqspt_plot` extracts a legend 
positioned on the bottom of a plot to be included in the plot output if 
`InclKey = TRUE`, this was previously throwing an error.

* Added the arguments `LegendWidthIn` and `LegendHeightIn` to `wqspt_plot` to 
control the bottom legend width and height, respectively.

* Added some examples to the documentation for `wqspt_plot` to better illustrate 
its use.

* Replaced all `b1_constr` arguments with `b_constr` to match the change in the 
name of this argument in the latest version of the `gWQS` package (v3.0.5).

* Replaced URL links to referenced papers with doi.org links.
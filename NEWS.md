# wqspt 1.0.2

* This function previously did not accommodate tibbles for the `data` argument 
of `wqs_full_perm`, but now `data` can be a `data.frame` object or a tibble.

* Added a `utils.R` script with a function `get_legend2()` to replace the 
malfunctioning `cowplot::get_legend()` call in the function `wqspt_plot`. In 
`cowplot` v1.1.3, `get_legend()` returns an error if the legend position is 
anywhere besides on the right of the plot. Since `wqspt_plot` extracts a legend 
positioned on the bottom of a plot to be included in the plot output if 
`InclKey = TRUE`, this was previously throwing an error.

* Added a `gwqs_hpc` function to `utils.R` along with its necessary hidden 
functions imported from the `gWQS` package to add an option to specify a 
number of workers for parallel processes for `wqs_pt` and `wqs_full_perm`. The 
`gwqs` function in the `gWQS` package (v3.0.5) uses as many parallel processes 
as there are cores detected, which can be problematic with high-performance 
computing (HPC) environments. HPC schedulers allocate a specific number of 
cores, and if `gwqs` tries to use more than that number, it will terminate the 
HPC job. Therefore, `gwqs_hpc` was created to allow for the number of parallel 
processes to be specified, thereby avoiding this problem.

* Added `...` arguments to `wqs_pt` and `wqs_full_perm` that ensure that the 
additional arguments passed to `gwqs_hpc` are passed to every iteration of that 
function, whereas previously those additional arguments only passed to the 
main WQS regression in `wqs_full_perm`.

* Added the arguments `LegendWidthIn` and `LegendHeightIn` to `wqspt_plot` to 
control the bottom legend width and height, respectively.

* Added some examples to the documentation for `wqspt_plot` to better illustrate 
its use.

* Replaced all `b1_constr` arguments with `b_constr` to match the change in the 
name of this argument in the latest version of the `gWQS` package (v3.0.5).

* Replaced URL links to referenced papers with doi.org links.

* Added additional examples of papers using the WQSPT method.
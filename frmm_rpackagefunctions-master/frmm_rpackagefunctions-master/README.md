# README.md, 22-6-20

Files and directories forming the skeleton (and its contents) of an R package:

- files: DESCRIPTION, FRMM2.Rproj, NAMESPACE;

- directories: R, man, data.

At the moment, in the R directory, we have

- functions1_3-6-20.R with the code of the functions, plus
- y.summer.m.transf.loess-data.R, x1.spring.m.transf.loess-data.R, x2.autumn.m.transf.loess-data.R, x3.winter.m.transf.loess-data.R.

The man directory has:

- binary_matrix_truth_function.Rd, 
- binary_vectors_lth_partition_function.Rd, 
- consensus_matrix_fun.Rd,
- iterative_clustering_function.Rd (this function contains operations in 1 run of FRMM),
- prop_correct_pairs_lth_partition_and_truth_function.Rd (this function computes the adjusted Rand Index plus
the Rand index, True Positive Rate and True Negative Rate between a partition and an index partition; all measures are 
multiplied by 100).

The functions names have the underscores of the manual files above substituted by dots. E.g. iterative.clustering.function 
corresponds to iterative_clustering_function.Rd, and so on.

The data directory has 4 rda files, each with the (transformed+loess) observed variables from the real seasonal data set
in the Nature Plants article. These are:

- y.summer.m.transf.loess.rda,
- x1.spring.m.transf.loess.rda,
- x2.autumn.m.transf.loess.rda,
- x3.winter.m.transf.loess.rda.

Please see the explanations and examples in the documentation files (in man), plus the explanatory comments and other 
code in functions1_3-6-20.R.




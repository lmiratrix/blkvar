
# The blkvar package

This code bundle is a package for estimating impacts and cross-site variation in multi-site or blocked RCTs.
It also has a bunch of routines for running multisite simulation scenarios, in particular code for generating multisite data parameterized by parameters such as the intraclass correlation coefficient and degree of cross-site impact variation.
It is used in a variety of research projects connected with the Miratrix CARES Lab (https://cares.gse.harvard.edu/) and other research that Miratrix is part of.

A core method of the package is `compare_methods()`, which uses a bunch of different estimation strategies on multisite data to estimate the ATE, and then makes a table of the results.

Also see `block_estimator()` which implements the various variance estimators for blocked experiments as discussed in the Pashley & Miratrix JEBS papers.

Changes in version 0.99.106 (2020-06-29)
+ Submitted to Bioconductor

Changes in version 1.0.0 (2020-10-27)
+ Accepted by Bioconductor. Included in RELEASE 3.12

Changes in version 1.0.1 (2020-01-26)
+ Fixed bug. Tests calling gprofiler now fixed to avoir error based on databases updates
+ Corrected undocumented params

Changes in version 1.1.5
+ Added a sub module detection function : get_sub_clusters()
+ Added a coloration argument on plot_module() to color graph according to given groups
+ Added warnings on suspicious power fitted and FAQ to help their resolution (get_fit.cor(), get_fit.expr(), build_net())
+ Added handling of NA in model matrix of associate_phenotype() and warnings when they are present
+ Added use of other correlation function in associate_phenotype()
+ Added a ggplot2 like palette generator (internal)
+ Removed hard coded arguments : 
  o alpha for plot_expression_profiles()
  o text angle for plot_comparison_stats(), plot_modules_phenotype()
  o deepSplit for detect_modules()
+ Changed WGCNA::cor() + WGCNA::corPvalueStudent() functions for WGCNA::corAndPvalue()
+ Changed erroneous info on SVD/PCA in FAQ

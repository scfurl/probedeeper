# probedeeper
Making the Analysis of Microarray Data More Efficient

##Welcome to the probedeeper library

Probedeeper is simply a library that allows for encapsulation of the following key elements for expression analysis:

1. An ExpressionSet [eset] (See Bioconductor "Biobase" package)
2. A Color Object [ColObj] - this is a probedeeper object that holds:
  a. An "assign" data-frame of two columns, "Color" and "Group".  This data-frame is a repository for color assignments
  b. A "classvec" factor of samples in the eset (this should be equivalent to as.factor(colnames(exprs(eset)))
  c. The set of colors for each sample (full).
  d. The set of colors for each type of sample (match).
3. A object holding differential expression output from the Limma package.

As this package is still in development, there is little to no documentation beyond which is written here.  Good Luck!

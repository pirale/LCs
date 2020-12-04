# LCs

The folder LCs contains the R script LCs.R and all files required for the computation of localized clusters. To run LCs.R, you will need R 3.4.0 or higher and the libraries reshape2, igraph, and parallel. The results of the computations will be stored in the variables _allScores_, containing the scoring for each of the clusterings; _cutoffs_, the species-specific SNV cutoffs with which the highest scoring clusterings were computed; and _clutering_, the highest-scoring clustering.

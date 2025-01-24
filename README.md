![](https://github.com/alexyermanos/AntibodyForests/blob/main/vignettes/imgs/main/logo.png)

# Delineating inter- and intra-antibody repertoire evolution with AntibodyForests

The generated wealth of immune repertoire sequencing data requires software to investigate and quantify inter- and intra-antibody repertoire evolution to uncover how B cells evolve during immune responses. Here, we present AntibodyForests, a software to investigate and quantify inter- and intra-antibody repertoire evolution.  

This R package is currently composed of a pipeline to reconstruct lineage trees from 10x single-cell V(D)J sequencing data preprocessed with the [Platypus package](https://github.com/alexyermanos/Platypus) and compare trees within and across repertoires. Furthermore, it has modalities to integrate bulk RNA sequencing data, features of protein 3D structure, and evolutionary likelihoods generated with protein language models.  

# Installation
Both Platypus and AntibodyForests can be installed from CRAN.

```r
#Install from CRAN
install.packages("Platypus")
install.packages("AntibodyForests")
```

# Quick Start

This quick start gives a short use case of AntibodyForests. Single-cell V(D)J sequencing 10x output of five mice immunized with Ovalbumin (OVA) from [Neumeier et al. (2022)](https://doi.org/10.1073/pnas.2113766119) are used to create a VDJ dataframe with [Platypus](https://github.com/alexyermanos/Platypus). AntibodyForests is used to create lineage trees for each B cell clonotype using an MST-like algorithm.

```r
#Load the libraries
library(Platypus)
library(AntibodyForests)

# Import 10x Genomics output files into VDJ dataframe, only keep cells with one VDJ and one VJ transcript, and trim the germline sequences
VDJ_OVA <- VDJ_build(VDJ.directory = "10x_output/VDJ/",
                     remove.divergent.cells = TRUE,
                     complete.cells.only = TRUE,
                     trim.germlines = TRUE)

# Build lineage trees for all clones present in the VDJ dataframe with the default algorithm
AntibodyForests_OVA <- Af_build(VDJ = VDJ_OVA, construction.method = "phylo.network.default")

# Plot one of the lineage trees as an example
Af_plot_tree(AntibodyForests_object = AntibodyForests_OVA, sample = "S1", clonotype = "clonotype3")
```
![](https://github.com/alexyermanos/AntibodyForests/blob/main/vignettes/imgs/QuickStart/Tree_OVA_s1_clonotype3.png)

Now we cluster the trees in this AntibodyForests object based on the Jensen-Shannon divergence between the Spectral Density profiles. We visualize the results in a heatmap and observe two clusters.

```r
# Cluster the trees that contain at least 8 nodes
out <- Af_compare_within_repertoires(input = AntibodyForests_OVA
                                     min.nodes = 8,
                                     distance.method = "jensen-shannon",
                                     clustering.method = "mediods",
                                     visualization.methods = "heatmap")
# Plot the heatmap
out$plots$heatmap_clusters
```
![](https://github.com/alexyermanos/AntibodyForests/blob/main/vignettes/imgs/QuickStart/quick-start_heatmap.png)

When we analyze the difference between the clusters, we observe that trees in cluster 2 have deep branching events indicated by the negative asymmetry index and contain multiple spectral density modalities. This indicates that various events of diversification took place during the evolution of these clonotypes and that cells with a small amount of SHM were recovered.

```r
# Analyze the difference between the clusters
plots <- Af_cluster_metrics(input = AntibodyForests_OVA_default,
                   clusters = out$clustering,
                   metrics = "spectral.density",
                   min.nodes = 8,
                   significance = T)

plots$spectral.asymmetry
plots$modalities
```
![](https://github.com/alexyermanos/AntibodyForests/blob/main/vignettes/imgs/QuickStart/QuickStart_boxplot.png)

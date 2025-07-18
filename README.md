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

This quick start gives a short use case of AntibodyForests. Single-cell V(D)J sequencing 10x output of five mice immunized with Ovalbumin (OVA) from [Neumeier et al. (2022)](https://doi.org/10.1073/pnas.2113766119) are used to create a VDJ dataframe with [Platypus](https://github.com/alexyermanos/Platypus). AntibodyForests is used to create lineage trees for each B cell clonotype using an MST-like algorithm and for downstream analysis.

## Construct lineage trees

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

## Quantify evolution

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
plots <- Af_cluster_metrics(input = AntibodyForests_OVA,
                   clusters = out$clustering,
                   metrics = "spectral.density",
                   min.nodes = 8,
                   significance = T)

plots$spectral.asymmetry
plots$modalities
```
![](https://github.com/alexyermanos/AntibodyForests/blob/main/vignettes/imgs/QuickStart/QuickStart_boxplot.png)

## Protein Language Model Likelihoods

We can calculate PLM likelihoods of the sequences in the lineage trees with the [PLM-pipeline](https://github.com/dvginneken/PLM-pipeline)

```r
# Extract the heavy chain sequences from the AntibodyForests objects
df <- Af_get_sequences(AntibodyForests_object = AntibodyForests_OVA, 
                       sequence.name = "VDJ_sequence_aa_trimmed")

#Save the data frame as CSV to serve as input for the PLM-pipeline
write.csv(df, file = "path/to/PLM_input.csv", row.names = FALSE)
```

Clone the repository of the [PLM-pipeline](https://github.com/dvginneken/PLM-pipeline) and run these example in the command line.

```bash
# Get the per-position likelihoods of each amino acid
python3 path/to/PLM-pipeline/scripts/pipeline.py --model_name ESMC \
--file_path path/to/PLM_input.csv --sequence_column sequence \
--calc_list probability_matrix

# Get the per-sequence pseudolikelihoods
python3 path/to/PLM-pipeline/scripts/pipeline.py --model_name ESMC \
--file_path path/to/VDJ_OVA.csv --sequence_column VDJ_sequence_aa_trimmed \
--calc_list pseudolikelihood
```

Next, we use AntibodyForests to integrate the calculated likelihoods and create a PLM dataframe with substitution likelihoods over each edge in the trees.

```r
# Add pseudolikelihood as node feature to the AntibodyForests object
AntibodyForests_OVA_default <- Af_add_node_feature(AntibodyForests_object = AntibodyForests_OVA, 
                                                   feature.df = VDJ_OVA,
                                                   feature.names = "evo_likelihood")

#Create a PLM dataframe
PLM_dataframe <- Af_PLM_dataframe(AntibodyForests_object =  AntibodyForests_OVA, 
                                  sequence.name = "VDJ_sequence_aa_trimmed", 
                                  path_to_probabilities = "path/to/probability_matrices")
```

We can analyze the correlation between the distance to the germline and the PLM pseudolikelihood. Here we observe a minor negative correlation, meaning that sequences further away from the germline have an overall lower pseudolikelihood.

```r
# Plot pseudolikelihood against distance to the germline
Af_distance_scatterplot(AntibodyForests_object = AntibodyForests_OVA, 
                        node.features = "evo_likelihood",
                        distance = "edge.length",
                        min.nodes = 1,
                        correlation = "pearson",
                        color.by = "sample")
```
![](https://github.com/alexyermanos/AntibodyForests/blob/development/vignettes/imgs/QuickStart/quick-start_pseudolikelihood.png)

We can also analyze the rank of the subsitutions in the trees. This reveals that with this specific PLM (ESM Cambrian) does not capture clear patterns of SHM in our lineage trees.

```r
# Plot the substitution rank
Af_plot_PLM(PLM_dataframe = PLM_dataframe, 
            values = "substitution_rank", 
            group_by = "sample_id")
```
![](https://github.com/alexyermanos/AntibodyForests/blob/development/vignettes/imgs/QuickStart/quick-start_substitutionRank.png)


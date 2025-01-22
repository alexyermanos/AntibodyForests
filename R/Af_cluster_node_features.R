#' Function to create a barplot of the cluster composition of selected features from each tree in an AntibodyForests-object
#' @description Function to create a barplot of the cluster composition of selected features from each tree in an AntibodyForests-object
#' @param input AntibodyForests-object(s), output from Af_build()
#' @param features Character vector of features to include in the barplot. (these features need to be present in the nodes of the trees)
#' @param clusters Named vector with the cluster assignments of the trees, output from Af_compare_within_repertoires().
#' @param fill identify each unique feature per tree (unique, default), or assign the most observed feature to the tree (max)
#' @param colors Color palette to use for the features.
#' @param text.size Size of the text in the plot. Default is 12.
#' @param significance Logical, whether to add Chi-squared Test p-value to the plot. Default is FALSE.
#' @return A list with barplots for each provided feature.
#' @export
#' @examples
#' plot <- Af_cluster_node_features(input = AntibodyForests::small_af,
#'                                  clusters = AntibodyForests::compare_repertoire[["clustering"]],
#'                                  features = "isotype",
#'                                  fill = "max")
#' plot$isotype

Af_cluster_node_features <- function(input,
                                    features,
                                    clusters,
                                    fill,
                                    colors,
                                    text.size,
                                    significance){

  #Stop when no input is provided
  if(missing(input)){stop("Please provide a valid input object.")}
  if(missing(features)){stop("Please provide a valid features object.")}
  if(missing(clusters)){stop("Please provide a valid clusters object.")}
  #Set defaults
  if(missing(fill)){fill <- "unique"}
  if(missing(text.size)){text.size <- 12}
  if(missing(significance)){significance <- FALSE}
  #Check if input is valid
  if(!(fill %in% c("unique", "max"))){stop("Please provide a valid fill argument.")}
  if(!all(features %in% names(input[[1]][[1]]$nodes[[1]]))){stop("Features are not present in AntibodyForests-object.")}

  #Set global variable for CRAN
  count <- NULL

  #Function to get the node features for each tree
  get_node_features <- function(clonotype, features, fill){
    node_features <- c()
    for (feature in features){
      if(fill == "unique"){
        labels <- paste(sort(unique(unlist(lapply(clonotype$nodes, function(c){c[[feature]]})))), collapse = ";")
        }
      if(fill == "max"){
        if(all(is.na(unlist(lapply(clonotype$nodes, function(c){c[[feature]]}))))){labels <- NA
        }else{labels <- names(sort(table(unlist(lapply(clonotype$nodes, function(c){c[[feature]]}))), decreasing = TRUE))[1]}
      }
      node_features <- c(node_features, labels)
    }
    names(node_features) <- features
    return(node_features)
  }

  #Create a dataframe with the node features for each tree
  df <- t(as.data.frame(lapply(input, function(sample){
    lapply(sample, function(clonotype){
      get_node_features(clonotype, features, fill)
    })
  })))

  #Add column to this dataframe with the assigned clusters
  df <- as.data.frame(df)
  df$tree <- gsub("^X", "", rownames(df))
  clusters <- as.data.frame(clusters)
  clusters$tree <- rownames(clusters)
  df <- dplyr::left_join(df, clusters, by = "tree")
  #Only keep trees with an assigned cluster
  df <- df[!is.na(df$clusters),]
  #Stop if there is no match
  if(nrow(df) == 0){stop("Tree names of the clusters are not in the AntibodyForests-object.")}

  #Create barplots
  output_list <- list()
  for(feature in features){


    p <- ggplot2::ggplot(df, ggplot2::aes(x=as.factor(clusters), fill=!!rlang::sym(feature))) +
      ggplot2::geom_bar(position = "fill", stat = "count", width = 0.8) +
      ggplot2::theme_classic() +
      ggplot2::theme(text = ggplot2::element_text(size = text.size),
                     axis.ticks.x=ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_blank()) +
      ggplot2::xlab("Cluster") + ggplot2::ylab("Percentage of cells")
    if(!missing(colors)){p <- p + ggplot2::scale_fill_manual(values = colors)}
    if(significance){
      #Create a contigency table for the Chi-squared test
      df_temp <- df[,c("clusters", feature)]
      df_temp$count <- 1
      contigency_table <- as.matrix(tidyr::pivot_wider(df_temp, names_from = clusters, values_from = count, values_fn = length, values_fill = 0))
      contigency_table <- apply(contigency_table[,2:ncol(contigency_table)], 2, as.numeric)
      #Perform the Chi-squared test
      chi2 <- stats::chisq.test(contigency_table)$p.value
      #Add the p-value to the plot
      p <- p + ggplot2::ggtitle(paste("Chi-squared test p-value: ", round(chi2, 3)))
    }
    #Add to the output list
    output_list[[feature]] <- p
  }

  return(output_list)

}









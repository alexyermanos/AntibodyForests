#' Function to make a grouped boxplot of the node sizes (number of cells with the exact same sequence) from specific groups of lineage trees constructed with AntibodyForests.
#' @description Function to compare trees.
#' @param AntibodyForests_object AntibodyForests-object, output from Af_build()
#' @param node.feature Node feature in the AntibodyForests-object to compare.
#' @param min.nodes The minimum number of nodes for a tree to be included in this analysis (this included the germline)
#' @param groups Which groups to compare. These groups need to be in the node features of the AntibodyForests-object. Set to NA if all features should displayed. (default is NA)
#' If you want to compare IgM and IgG for example, groups should be c("IgM, "IgG") (not "Isotypes")
#' @param colors Optionally specific colors for the group (Will be matched to the groups/names on alphabetical order).
#' @param text.size Font size in the plot (default 20).
#' @param x.label Label for the x-axis (default is the node feature).
#' @param group.order Order of the groups on the x-axis. (default is alphabetical/numerical)
#' @param significance If TRUE, the significance of the difference (paired t-test) between the groups is plotted. (default FALSE)
#' @param parallel If TRUE, the metric calculations are parallelized across clonotypes. (default FALSE)
#' @param output.file string - specifies the path to the output file (PNG of PDF). Defaults to NULL.
#' @return A ggplot2 object with the boxplot.
#' @export
#' @examples
#' Af_distance_boxplot(AntibodyForests::small_af,
#'                     min.nodes = 5,
#'                     groups = c("IGHA", "IgG1"),
#'                     node.feature = "isotype",
#'                     unconnected = TRUE)
#'

Af_node_size_boxplot <- function(AntibodyForests_object,
                                min.nodes,
                                groups,
                                node.feature,
                                colors,
                                text.size,
                                x.label,
                                group.order,
                                significance,
                                parallel,
                                output.file){
  
  #Set defaults and check for missing input
  if(missing(AntibodyForests_object)){stop("Please provide an AntibodyForests-object as input.")}
  if(missing(node.feature)){stop("Please provide a node feature to compare.")}
  if(missing(groups)){groups = NA}
  if(missing(text.size)){text.size = 20}
  if(missing(min.nodes)){min.nodes = 0}
  if(missing(parallel)){parallel <- F}
  if(missing(x.label)){x.label = node.feature}
  if(missing(significance)){significance = F}
  if(missing(group.order)){group.order = NA}
  if(missing(output.file)){output.file <- NULL}
  #Check if group are in the metric dataframe
  #if(!(all(groups %in% colnames(metric_df)))){stop("Groups are not in the column names of the metric dataframe.")}
  
  message("This function takes a long runtime if the AntibodyForests-object is large.")
  
  #Global variable definitions for CRAN checks
  depth <- NULL
  clonotype <- NULL
  png <- NULL
  pdf <- NULL
  group <- NULL
  
  #Calculate the average node size per group
  metric_df <- Af_metrics(AntibodyForests_object,
                          parallel = parallel,
                          min.nodes = min.nodes,
                          metrics = "group.node.size",
                          node.feature = node.feature,
                          group.node.feature = groups)
  #Remove column with sample names
  df <- metric_df[,colnames(metric_df) != "sample"]
  
  #Error if zero or only one tree is in the metric_df
  if(is.null(nrow(df))){stop("Your AntibodyForests-object does not have enough trees that pass the min.nodes threshold.")}
  
  #Add clonotype as column
  df$clonotype <- rownames(df)
  
  #Check if there are clonotypes left after NA removal
  if(nrow(df) == 0){stop("No trees contain nodes from all groups.")}
  
  #Transform dataframe for visualization
  df <- tidyr::pivot_longer(df, cols=colnames(df)[1:ncol(df)-1],
                            names_to='group',
                            values_to='depth')
  df$group <- gsub(paste0(".node.size"), "", df$group)
  
  #Select all groups if groups is NA
  if(all(is.na(groups))){groups <- unique(df$group)}
  
  #Set colors if not provided
  if(missing(colors)){colors = scales::hue_pal()(length(groups))}
  
  #Set order of groups if provided
  if(!all(is.na(group.order))){
    df$group <- factor(df$group, levels = group.order)
  }
  
  #Plot the grouped boxplots with lines
  p <- ggplot2::ggplot(df, ggplot2::aes(x=group, y=depth, fill=group)) +
    ggplot2::geom_boxplot()+
    ggplot2::scale_fill_manual(values=colors) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = text.size),
                   legend.position = "none")  +
    ggplot2::xlab(x.label) + ggplot2::ylab(paste0("Average node size"))
  
  #Add significance to the plot
  if(significance){
    #Get the unique combinations of groups if there are more than 2 groups
    if(length(groups) > 2){
      #Get the unique combinations of clusters
      combinations <- combinat::combn(unique(df$group), 2)
      combinations_list <- split(combinations, col(combinations))
    }else{
      combinations_list <- list(unique(df$group))
    }
    #Add to the existing plot
    p <- p + ggsignif::geom_signif(comparisons=combinations_list, step_increase = 0.1, test = "t.test",
                                   test.args = list(paired = T))
  }
  
  if(!is.null(output.file)){
    # Check if the output.file is png or pdf
    if (grepl(pattern = ".png$", output.file)){
      png(file = output.file)
      print(p)
      grDevices::dev.off()
    }else if (grepl(pattern = ".pdf$", output.file)){
      pdf(file = output.file)
      print(p)
      grDevices::dev.off()
    }
  }
  
  return(p)
  
  
}
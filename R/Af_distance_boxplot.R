#' Function to make a grouped boxplot of distance between nodes from specific groups and the germline of lineage trees constructed with AntibodyForests.
#' @description Function to compare trees.
#' @param AntibodyForests_object AntibodyForests-object, output from Af_build()
#' @param node.feature Node feature in the AntibodyForests-object to compare.
#' @param distance - string - How to calculate the distance to the germline.
#' 'node.depth'     : Average of the sum of edges on the shortest parth between germline and nodes from this group. 
#' 'edge.length'    : Average of the sum of edge length of the shortest path between germline and nodes from this group. (Default)
#' @param min.nodes The minimum number of nodes for a tree to be included in this analysis (this included the germline)
#' @param groups Which groups to compare. These groups need to be in the node features of the AntibodyForests-object. Set to NA if all features should displayed. (default is NA)
#' If you want to compare IgM and IgG for example, groups should be c("IgM, "IgG") (not "Isotypes")
#' @param unconnected If TRUE, trees that don't have all groups will be plotted, but not included in significance analysis. (default FALSE)
#' @param colors Optionally specific colors for the group (Will be matched to the groups/names on alphabetical order).
#' @param text.size Font size in the plot (default 20).
#' @param x.label Label for the x-axis (default is the node feature).
#' @param group.order Order of the groups on the x-axis. (default is alphabetical/numerical)
#' @param significance If TRUE, the significance of the difference (paired t-test) between the groups is plotted. (default FALSE)
#' @param parallel If TRUE, the metric calculations are parallelized across clonotypes. (default FALSE)
#' @param output.file string - specifies the path to the output file (PNG of PDF). Defaults to NULL.
#' @export
#' @examples
#' Af_distance_boxplot(AntibodyForests_object, 
#'   distance = "edge.length", 
#'   min.nodes = 5, 
#'   groups = c("IgM", "IgG"), 
#'   node.feature = "isotype", 
#'   unconnected = TRUE, 
#'   colors = c("red", "blue"), 
#'   text.size = 20, 
#'   x.label = "Isotype", 
#'   group.order = c("IgM", "IgG"), 
#'   significance = TRUE, 
#'   parallel = FALSE, 
#'   output.file = "output.png")
#' 

Af_distance_boxplot <- function(AntibodyForests_object,
                                     distance,
                                     min.nodes,
                                     groups,
                                     node.feature,
                                     unconnected,
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
  if(missing(distance)){distance = "edge.length"}
  if(missing(text.size)){text.size = 20}
  if(missing(min.nodes)){min.nodes = 0}
  if(missing(parallel)){parallel <- F}
  if(missing(x.label)){x.label = node.feature}
  if(missing(significance)){significance = F}
  if(missing(unconnected)){unconnected = F}
  if(missing(group.order)){group.order = NA}
  if(missing(output.file)){output.file <- NULL}  
  #Check if group are in the metric dataframe
  #if(!(all(groups %in% colnames(metric_df)))){stop("Groups are not in the column names of the metric dataframe.")}
  
  print("Warning: This function takes a long runtime if the AntibodyForests-object is large.")
  
  #Calculate the average distance to the germline per group
  metric_df <- Af_metrics(AntibodyForests_object,
                                       parallel = parallel,
                                       min.nodes = min.nodes,
                                       metrics = paste0("group.",distance),
                                       node.feature = node.feature,
                                       group.node.feature = groups)
  #Remove column with sample names
  df_all <- metric_df[,colnames(metric_df) != "sample"]
  
  #Error if zero or only one tree is in the metric_df
  if(is.null(nrow(df_all))){stop("Your AntibodyForests-object does not have enough trees that pass the min.nodes threshold.")}
  
  #Add clonotype as column
  df_all$clonotype <- rownames(df_all)
  
  #Only keep clonotypes that have nodes of all groups
  df <- as.data.frame(na.omit(df_all))
  
  #Check if there are clonotypes left after NA removal
  if(nrow(df) == 0){stop("No trees contain nodes from all groups.")}
  
  #Transform dataframe for visualization
  df <- tidyr::pivot_longer(df, cols=colnames(df)[1:ncol(df)-1],
                             names_to='group',
                             values_to='depth')
  
  #Select all groups if groups is NA
  if(all(is.na(groups))){groups <- gsub(paste0(".",distance), "", unique(df$group))}
  
  #Set colors if not provided
  if(missing(colors)){colors = scales::hue_pal()(length(groups))}
  
  #Set order of groups if provided
  if(!all(is.na(group.order))){
    df$group <- factor(df$group, levels = paste0(group.order,".", distance))
  }
  
  #Plot the grouped boxplots with lines
  p <- ggplot2::ggplot(df, ggplot2::aes(x=group, y=depth, fill=group)) + 
    ggplot2::geom_boxplot()+ 
    ggplot2::geom_point()+ 
    ggplot2::scale_fill_manual(values=colors) +
    ggplot2::geom_line(ggplot2::aes(group=clonotype)) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = text.size),
                   legend.position = "none")  +
    # ggplot2::scale_x_discrete(breaks=paste0(groups,".", distance),
    #                  labels=groups) +
    ggplot2::ggtitle(paste0("Distance (", distance, ") to germline")) +
    ggplot2::xlab(x.label)
  
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
  
  #Add unconnected points
  if(unconnected){
    #Create dataframe with trees that don't have all groups
    df_na <- df_all[rowSums(is.na(df_all)) > 0,]
    #Transform dataframe for visualization
    df_na <- tidyr::pivot_longer(df_na, cols=colnames(df_na)[1:ncol(df_na)-1],
                              names_to='group',
                              values_to='depth')
    #Add to the plot
    p <- p + ggplot2::geom_point(data = df_na, color = "darkgrey", ggplot2::aes(x=group, y=depth))
  }
  
  if(!is.null(output.file)){
    # Check if the output.file is png or pdf
    if (grepl(pattern = ".png$", output.file)){
      png(file = output.file)
      print(p)
      dev.off()
    }else if (grepl(pattern = ".pdf$", output.file)){
      pdf(file = output.file)
      print(p)
      dev.off()
    }
  }else{print(p)}
  

}

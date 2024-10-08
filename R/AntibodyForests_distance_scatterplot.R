#' Function to scatterplot the distance to the germline to a numerical node feature of the AntibodyForests-object
#' @description Function to scatterplot the distance to the germline to a numerical node feature of the AntibodyForests-object
#' @param input AntibodyForests-object(s), output from AntibodyForests()
#' @param node.features Node features in the AntibodyForests-object to compare (needs to be numerical)
#' @param min.nodes The minimum number of nodes for a tree to be included in this analysis (this included the germline). Default is 2.
#' @param color.by Color the scatterplot by a node.feature in the AntibodyForests-object, by the sample, or no color ("none). Default is "none".
#' @param color.by.numeric Logical. If TRUE, the color.by feature is treated as a numerical feature. Default is FALSE.
#' @param correlation "pearson", "spearman", "kendall", or "none"
#' @param color.palette The color palette to use for the scatterplot. Default for numerical color.by is "viridis".
#' @param font.size The font size of the plot. Default is 12.
#' @param output.file string - specifies the path to the output file (PNG of PDF). Defaults to NULL.
#' @return 
#' @export

AntibodyForests_distance_scatterplot <- function(input,
                                                 node.features,
                                                 min.nodes,
                                                 color.by,
                                                 color.by.numeric,
                                                 correlation,
                                                 color.palette,
                                                 font.size,
                                                 output.file){
  
  #Set defaults and check for missing input
  if(missing(input)){stop("Please provide an AntibodyForests-object as input.")}
  if(missing(node.features)){stop("Please provide a node feature to compare.")}
  # If the node.features could not be found for all nodes, a message is returned and execution is stopped
  for(feature in node.features){if(!(all(sapply(names(input[[sample]][[clonotype]][["nodes"]])[!names(input[[sample]][[clonotype]][["nodes"]]) == "germline"], function(x) feature %in% names(input[[sample]][[clonotype]][["nodes"]][[x]]))))){stop("The feature specified with the 'node.features' parameter could not be found for all nodes.")}}  
  if(missing(min.nodes)){min.nodes <- 2}
  if(missing(color.by)){color.by <- "none"}
  # If the color.by feature could not be found for all nodes, a message is returned and execution is stopped
  if(color.by != "sample" & !(all(sapply(names(input[[sample]][[clonotype]][["nodes"]])[!names(input[[sample]][[sample]][[clonotype]][["nodes"]]) == "germline"], function(x) color.by %in% names(input[[sample]][[clonotype]][["nodes"]][[x]]))))){stop("The feature specified with the 'color.by' parameter could not be found for all nodes.")}
  if(missing(color.by.numeric)){color.by.numeric <- F}
  if(missing(correlation)){correlation <- "none"}
  if(!(correlation %in% c("none", "spearman", "pearson", "kendall"))){stop("Please provide a valid correlation method ('none', 'spearman', 'pearson', or 'kendall').")}
  if(missing(color.palette)){color.palette <- NULL}
  if(missing(font.size)){font.size <- 12}
  if(missing(output.file)){output.file <- NULL}  
  
  #Create input for ggplot
  #Create empty dataframe
  df <- data.frame()
  #Loop though each sample and clonotype
  for (sample in names(input)){
    for (clonotype in names(input[[sample]])){
      #Get the igraph object of the tree
      tree <- input[[sample]][[clonotype]][["igraph"]]
      
      if (igraph::vcount(tree) >= min.nodes){
        #Calculate distances to germline
        distances_df <- t(igraph::distances(tree,
                                            v = "germline",
                                            to = igraph::V(tree)[names(igraph::V(tree)) != "germline"],
                                            weights = as.numeric(igraph::edge_attr(tree)$edge.length)))
        
        #Add sample name, clonotype, and node number
        distances_df <- cbind(distances_df, sample = replicate(nrow(distances_df), sample), 
                              clonotype = replicate(nrow(distances_df), clonotype),
                              node = rownames(distances_df))
        
        #Add color.by to the node features
        if (!color.by %in% c("none", "sample")){adding.features <- c(node.features, color.by)}
        
        #Add the node features to the distances dataframe
        for (feature in adding.features){
          #Get the node feature
          feature_list <- lapply(input[[sample]][[clonotype]][["nodes"]], function(x){
            value = unique(x[[feature]])
            if(length(value) == 1){return(value)}else{return(NA)}
          })
          feature_list <- feature_list[which(names(feature_list) != "germline")]
          feature_df <- t(as.data.frame(feature_list))
          colnames(feature_df) <- feature
          feature_df <- cbind(feature_df, sample = replicate(nrow(feature_df), sample), 
                              clonotype = replicate(nrow(feature_df), clonotype),
                              node = rownames(feature_df))
          
          #Merge the feature dataframe with the distances dataframe
          distances_df <- dplyr::left_join(as.data.frame(distances_df), as.data.frame(feature_df), by = c("sample", "clonotype", "node"))
        }
        
        #Add to the final dataframe
        df <- rbind(df, distances_df)
        
        #Make the node features numerical
        if(length(node.features) > 1){df[, node.features] <- apply(df[, node.features], 2, function(x) as.numeric(as.character(x)))}
        else if (length(node.features == 1)){df[, node.features] <- as.numeric(as.character(df[,node.features]))}
        
        
        #Make the color.by feature numerical
        if(color.by.numeric){df[, color.by] <- as.numeric(as.character(df[[color.by]]))}
      }
    }
  }
  
  #Create the plots
  for(feature in node.features){
    #Create the scatterplot
    p <- ggplot2::ggplot(df, ggplot2::aes(x = as.numeric(germline), y = .data[[feature]])) +
      ggplot2::geom_point() +
      ggplot2::xlab("Distance to germline") +
      ggplot2::theme_classic() +
      ggplot2::theme(text = ggplot2::element_text(size = font.size))
    
    #Color by color.by feature
    if (color.by != "none"){
      p <- p + ggplot2::aes(color = .data[[color.by]])
    }
    
    #If color.palette is provided, use that
    if(!is.null(color.palette)){
      #If the color.by is numeric, create a gradient between the two color.palette colors
      if(color.by.numeric){p <- p + ggplot2::scale_colour_gradient(low = color.palette[1], high = color.palette[2])}
      #If the color.by is not numeric, use the color.palette as is
      else{p <- p + ggplot2::scale_colour_manual(values = color.palette)}
    }
    #If no color.palette is provided and the color.by is numeric, color by viridis palette
    else if(is.null(color.palette) & color.by.numeric){p <- p + viridis::scale_color_viridis()}
    
    if (correlation != "none"){
      cor <- stats::cor.test(as.numeric(df$germline), df[,feature], method = correlation, exact = F)
      p <- p + ggplot2::ggtitle(paste0(correlation, " R\u00b2 = ", round(cor$estimate, digits = 2)))
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
  
}
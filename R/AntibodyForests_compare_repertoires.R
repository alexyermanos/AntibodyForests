AntibodyForests_compare_repertoires <- function(AntibodyForests_list, #list of antibodyforests objects to compare
                                                metrics, #which metrics to use
                                                plot){ #freqpoly or boxplot
  
  #Set defaults and check for input validity
  if(missing(AntibodyForests_list)){stop("A list of AntibodyForests objects must be provided")}
  if(missing(metrics)){metrics <- c("betweenness", "degree")}
  if(missing(plot)){plot <- "freqpoly"}
  if (plot != "freqpoly" & plot != "boxplot"){stop("plot must be either 'freqpoly' or 'boxplot'")}
  if (!(all(metrics %in% c("nr.nodes", "nr.cells", "mean.depth", "mean.edge.length", 
                           "sackin.index", "spectral.density", "degree", "betweenness")))){
    stop("Unknown metric")
  }
  
  output_plot <- list()
  
  #Name the list of AntibodyForests objects if not named
  if (is.null(names(AntibodyForests_list))){names(AntibodyForests_list) <- paste0("group",1:length(AntibodyForests_list))}
  
  # If betweenness or degree is in the metrics, all igraph object must be combined
  if ("betweenness" %in% metrics | "degree" %in% metrics){
    igraph_list <- list()
    #Combine all the igraph objects in each group
    for (group in names(AntibodyForests_list)){
      #Store all edges of this group
      edges <- c()
      #Loop over every tree of this group
      for (sample in names(AntibodyForests_list[[group]])){
        for (clonotype in names(AntibodyForests_list[[group]][[sample]])){
          #Get the igraph object
          tree <- AntibodyForests_list[[group]][[sample]][[clonotype]]$igraph
          #Change node names to include sample and clonotype information
          tree <- igraph::set_vertex_attr(tree, "name", value=paste0(sample, "_", clonotype,"_",igraph::V(tree)$name))
          #Add the tree to the list
          edges <- rbind(edges, igraph::as_edgelist(tree, names = TRUE))
        }
      }
      #Combine all trees of this group
      combined_igraph <- igraph::graph_from_edgelist(edges, directed = TRUE)
      #Add the igraph of this group to the igraph_list
      igraph_list[[group]] <- combined_igraph
    }
    
    #Calculate and plot the betweenness
    if ("betweenness" %in% metrics){
      #Store the betweenness values per group
      betweenness_list <- list()
      for (group in names(AntibodyForests_list)){
        #Calculate the betweenness and write to dataframe
        betweenness_df <- data.frame(igraph::betweenness(igraph_list[[group]]))
        colnames(betweenness_df) <- "betweenness"
        betweenness_df$group <- group
        #Add to the betweenness_list
        betweenness_list[[group]] <- betweenness_df
      }
      #Combine the betweenness dataframes
      betweennness_df <- do.call("rbind", betweenness_list)
      
      #Plot
      if (plot == "freqpoly"){
        output_plot[["betweenness"]] <- ggplot2::ggplot(betweennness_df, ggplot2::aes(betweenness)) +
          ggplot2::geom_freqpoly(binwidth = max(as.numeric(betweenness_df$betweenness))/1000, ggplot2::aes(colour = group)) +
          ggplot2::scale_x_continuous(trans='log10') +
          ggplot2::theme_minimal()
      }
      if (plot == "boxplot"){
        output_plot[["betweenness"]] <- ggplot2::ggplot(betweennness_df, ggplot2::aes(x = group, y = betweenness, fill = group)) +
          ggplot2::geom_boxplot() +
          ggplot2::scale_y_continuous(trans='log10') +
          ggplot2::theme_minimal()
      }
      
    }
    #Calculate and plot the degree
    if ("degree" %in% metrics){
      #Store the degree values per group
      degree_list <- list()
      for (group in names(AntibodyForests_list)){
        #Calculate the degree and write to dataframe
        degree_df <- data.frame(igraph::degree(igraph_list[[group]]))
        colnames(degree_df) <- "degree"
        degree_df$group <- group
        #Add to the degree_list
        degree_list[[group]] <- degree_df
      }
      #Combine the degree dataframes
      degree_df <- do.call("rbind", degree_list)
      
      #Plot
      if (plot == "freqpoly"){
        output_plot[["degree"]] <- ggplot2::ggplot(degree_df, ggplot2::aes(degree)) +
          ggplot2::geom_freqpoly(binwidth = max(as.numeric(betweenness_df$betweenness))/1000, ggplot2::aes(colour = group)) +
          ggplot2::scale_x_continuous(trans='log10') +
          ggplot2::theme_minimal()
      }
      if (plot == "boxplot"){
        output_plot[["degree"]] <- ggplot2::ggplot(degree_df, ggplot2::aes(x = group, y = degree, fill = group)) +
          ggplot2::geom_boxplot() +
          ggplot2::scale_y_continuous(trans='log10') +
          ggplot2::theme_minimal()
      }
    }
    
    
  }
  
  # For the metrics "nr.nodes", "nr.cells", "mean.depth", "mean.edge.length", "sackin.index" and "spectral.density" run AntibodyForests_metrics
  if (any(metrics %in% c("nr.nodes", "nr.cells", "mean.depth", "mean.edge.length", "sackin.index", "spectral.density"))){
    #Store metrics per group
    metrics_list <- list()
    for (group in names(AntibodyForests_list)){
      #Calculate the metrics
      if("spectral.density" %in% metrics){min.nodes = 3}else{min.nodes = 1}
      metrics_df <- AntibodyForests_metrics(input = AntibodyForests_list[[group]], metrics = metrics, min.nodes)
      metrics_df$group <- group
      #Add to the metrics_list
      metrics_list[[group]] <- metrics_df
    }
    #Combine the metrics dataframes
    metrics_df <- do.call("rbind", metrics_list)
    
    #Set new metric names
    metrics <- colnames(metrics_df)[!(colnames(metrics_df) %in% c("sample", "group"))]
    
    #Plot
    for (metric in metrics){
      temp_metrics_df <- metrics_df[!is.na(metrics_df[,metric]),]
      if (plot == "freqpoly"){
        output_plot[[metric]] <- ggplot2::ggplot(temp_metrics_df, ggplot2::aes(.data[[metric]])) +
          ggplot2::geom_freqpoly(binwidth = max(as.numeric(temp_metrics_df[,metric]))/30, ggplot2::aes(colour = group)) +
          #ggplot2::scale_x_continuous(trans='log10') +
          ggplot2::theme_minimal() +
          ggplot2::xlab(metric)
      }
      if (plot == "boxplot"){
        output_plot[[metric]] <- ggplot2::ggplot(metrics_df, ggplot2::aes(x = group, y = .data[[metric]], fill = group)) +
          ggplot2::geom_boxplot() +
          #ggplot2::scale_y_continuous(trans='log10') +
          ggplot2::theme_minimal() +
          ggplot2::ylab(metric)
      }
    }
  }
  return(output_plot)
}

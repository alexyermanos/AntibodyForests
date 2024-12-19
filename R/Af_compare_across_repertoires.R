#' A function to compare dynamics of B cell evolution across different repertoires.
#' @description Compare tree topology metrics across different (groups) of AntibodyForests objects.
#' @param AntibodyForests_list A list of AntibodyForests objects to compare.
#' @param metrics Which metrics to use for comparison. Options are:
#'. betweenness       : The number of shortest paths that pass through each node (Default)
#'. degree            : The number of edges connected to each node (Default)
#' 'nr.nodes'         : The total number of nodes
#' 'nr.cells'         : The total number of cells in this clonotype
#' 'mean.depth'       : Mean of the number of edges connecting each node to the germline
#' 'mean.edge.length' : Mean of the edge lengths between each node and the germline
#' 'sackin.index'     : Sum of the number of nodes between each terminal node and the germline, normalized by the number of terminal nodes
#' 'spectral.density' : Metrics of the spectral density profiles (calculated with package RPANDA)
#'    - peakedness            : Tree balance
#'    - asymmetry             : Shallow or deep branching events
#'    - principal eigenvalue  : Phylogenetic diversity
#'    - modalities            : The number of different structures within the tree
#' @param plot What kind of plot to make.
#'  boxplot (default)
#'  freqpoly
#' @param text.size Font size in the plot (default 20).
#' @param colors Optionally specific colors for the groups. If not provided, the default ggplot2 colors are used.
#' @param significane If TRUE, the significance of a T test between the groups is plotted in the boxplot (default FALSE)
#' @param parallel If TRUE, the metric calculations are parallelized (default FALSE)
#' @param num.cores Number of cores to be used when parallel = TRUE. (Defaults to all available cores - 1)
#' @return Plots to compare the repertoires on the supplied metrics.
#' @export
#' @examples
#' boxplots <- Af_compare_across_repertoires(list("S1" = AntibodyForests::small_af[1], "S2" = AntibodyForests::small_af[2]),
#'             metrics = c("sackin.index", "betweenness", "degree"),
#'             plot = "boxplot")
#' boxplots$betweenness

Af_compare_across_repertoires <- function(AntibodyForests_list,
                                          metrics,
                                          plot,
                                          text.size,
                                          colors,
                                          significance,
                                          parallel,
                                          num.cores){

  #Set defaults and check for input validity
  if(missing(AntibodyForests_list)){stop("A list of AntibodyForests objects must be provided")}
  if(missing(metrics)){metrics <- c("betweenness", "degree")}
  if(missing(plot)){plot <- "boxplot"}
  if(missing(text.size)){text.size <- 20}
  if(missing(significance)){significance <- FALSE}
  if(missing(colors)){colors = scales::hue_pal()(length(AntibodyForests_list))}
  if (plot != "freqpoly" & plot != "boxplot"){stop("plot must be either 'freqpoly' or 'boxplot'")}
  if (!(all(metrics %in% c("nr.nodes", "nr.cells", "mean.depth", "mean.edge.length",
                           "sackin.index", "spectral.density", "degree", "betweenness")))){
    stop("Unknown metric")
  }
  if(missing(parallel)){parallel <- FALSE}
  # If 'parallel' is set to TRUE but 'num.cores' is not specified, the number of cores is set to all available cores - 1
  if(parallel == TRUE && missing(num.cores)){num.cores <- parallel::detectCores() -1}

  output_plot <- list()

  #Name the list of AntibodyForests objects if not named
  if (is.null(names(AntibodyForests_list))){names(AntibodyForests_list) <- paste0("group",1:length(AntibodyForests_list))}

  # If betweenness or degree is in the metrics, all igraph object must be combined
  if ("betweenness" %in% metrics | "degree" %in% metrics){

    # If 'parallel' is set to FALSE, the VDJ dataframe construction is not parallelized
    if(!parallel){
      igraph_list <- lapply(names(AntibodyForests_list), function(group){
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
        return(combined_igraph)
      })
    }
    # If 'parallel' is set to TRUE, the VDJ dataframe construction is parallelized
    if(parallel){
      # Retrieve the operating system
      operating_system <- Sys.info()[['sysname']]
      # If the operating system is Linux or Darwin, 'mclapply' is used for parallelization
      if(operating_system %in% c('Linux', 'Darwin')){
        igraph_list <- parallel::mclapply(mc.cores = num.cores, names(AntibodyForests_list), function(group){
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
          return(combined_igraph)
        })
      }
      if(operating_system == "Windows"){
        # Create cluster
        cluster <- parallel::makeCluster(num.cores)

        igraph_list <- parallel::parLapply(cluster, names(AntibodyForests_list), function(group){
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
          return(combined_igraph)
        })
      }
      }

    names(igraph_list) <- names(AntibodyForests_list)
    #Calculate and plot the betweenness
    if ("betweenness" %in% metrics){
      #Store the betweenness values per group
      #betweenness_list <- list()
      betweenness_list <- lapply(names(AntibodyForests_list), function(group){
        #Calculate the betweenness and write to dataframe
        betweenness_df <- data.frame(igraph::betweenness(igraph_list[[group]]))
        colnames(betweenness_df) <- "betweenness"
        betweenness_df$group <- group
        return(betweenness_df)
      })

      #Combine the betweenness dataframes
      betweenness_df <- do.call("rbind", betweenness_list)

      #Plot
      if (plot == "freqpoly"){
        output_plot[["betweenness"]] <- ggplot2::ggplot(betweenness_df, ggplot2::aes(betweenness)) +
          ggplot2::geom_freqpoly(binwidth = max(as.numeric(betweenness_df$betweenness))/1000, ggplot2::aes(colour = group)) +
          ggplot2::scale_colour_manual(values=colors) +
          ggplot2::scale_x_continuous(trans='log10') +
          ggplot2::theme_classic() +
          ggplot2::theme(text = ggplot2::element_text(size = text.size),
                         legend.position = "none")
      }
      if (plot == "boxplot"){
        p <- ggplot2::ggplot(betweenness_df, ggplot2::aes(x = group, y = betweenness, fill = group)) +
          ggplot2::geom_boxplot() +
          ggplot2::scale_fill_manual(values=colors) +
          ggplot2::theme_classic() +
          ggplot2::theme(text = ggplot2::element_text(size = text.size),
                         legend.position = "none")

        #Add significance to the plot
        if(significance){
          #Get the unique combinations of groups if there are more than 2 clusters
          if(length(unique(betweenness_df$group)) > 2){
            #Get the unique combinations
            combinations <- combinat::combn(unique(betweenness_df$group), 2)
            combinations_list <- split(combinations, col(combinations))
          }else{combinations_list <- list(unique(betweenness_df$group))}
          #Add to the existing plot
          p <- p + ggsignif::geom_signif(comparisons=combinations_list, step_increase = 0.1, test = "t.test")
        }

        output_plot[["betweenness"]] <- p
      }

    }
    #Calculate and plot the degree
    if ("degree" %in% metrics){
      #Store the degree values per group
      #degree_list <- list()

      degree_list <- lapply(names(AntibodyForests_list), function(group){
        #Calculate the betweenness and write to dataframe
        degree_df <- data.frame(igraph::degree(igraph_list[[group]]))
        colnames(degree_df) <- "degree"
        degree_df$group <- group
        return(degree_df)
      })

      #Combine the degree dataframes
      degree_df <- do.call("rbind", degree_list)

      #Plot
      if (plot == "freqpoly"){
        output_plot[["degree"]] <- ggplot2::ggplot(degree_df, ggplot2::aes(degree)) +
          ggplot2::geom_freqpoly(binwidth = max(as.numeric(betweenness_df$betweenness))/1000, linewidth = 1, ggplot2::aes(colour = group)) +
          ggplot2::scale_colour_manual(values=colors) +
          ggplot2::scale_x_continuous(trans='log10') +
          ggplot2::theme_classic() +
          ggplot2::theme(text = ggplot2::element_text(size = text.size))
      }
      if (plot == "boxplot"){
        p <- ggplot2::ggplot(degree_df, ggplot2::aes(x = group, y = degree, fill = group)) +
          ggplot2::geom_boxplot() +
          ggplot2::scale_fill_manual(values=colors) +
          ggplot2::theme_classic() +
          ggplot2::theme(text = ggplot2::element_text(size = text.size),
                         legend.position = "none")

        #Add significance to the plot
        if(significance){
          #Get the unique combinations of groups if there are more than 2 clusters
          if(length(unique(degree_df$group)) > 2){
            #Get the unique combinations
            combinations <- combinat::combn(unique(degree_df$group), 2)
            combinations_list <- split(combinations, col(combinations))
          }else{combinations_list <- list(unique(degree_df$group))}
          #Add to the existing plot
          p <- p + ggsignif::geom_signif(comparisons=combinations_list, step_increase = 0.1, test = "t.test", y_position = max(degree_df$degree, na.rm = T))
        }
        #Log scale the y axis for visibility
        p <- p + ggplot2::scale_y_continuous(trans='log10')
        output_plot[["degree"]] <- p
      }
    }


  }

  # For the metrics "nr.nodes", "nr.cells", "mean.depth", "mean.edge.length", "sackin.index" and "spectral.density" run AntibodyForests_metrics
  if (any(metrics %in% c("nr.nodes", "nr.cells", "mean.depth", "mean.edge.length", "sackin.index", "spectral.density"))){
    #Store metrics per group
    #metrics_list <- list()

    metrics_list <- lapply(names(AntibodyForests_list), function(group){
      #Calculate the metrics
      if("spectral.density" %in% metrics){min.nodes = 3}else{min.nodes = 1}
      metrics_df <- Af_metrics(input = AntibodyForests_list[[group]], metrics = metrics[which(!metrics %in% c("betweenness", "degree"))], min.nodes = min.nodes, parallel = parallel)
      metrics_df$group <- group
      return(metrics_df)
    })

    #Combine the metrics dataframes
    metrics_df <- do.call("rbind", metrics_list)

    #Set new metric names
    metrics <- colnames(metrics_df)[!(colnames(metrics_df) %in% c("sample", "group"))]

    #Plot
    for (metric in metrics){
      temp_metrics_df <- metrics_df[!is.na(metrics_df[,metric]),]
      if (plot == "freqpoly"){
        output_plot[[metric]] <- ggplot2::ggplot(temp_metrics_df, ggplot2::aes(.data[[metric]])) +
          ggplot2::geom_freqpoly(binwidth = max(as.numeric(temp_metrics_df[,metric]))/30, linewidth = 1, ggplot2::aes(colour = group)) +
          #ggplot2::scale_x_continuous(trans='log10') +
          ggplot2::theme_classic() +
          ggplot2::theme(text = ggplot2::element_text(size = text.size)) +
          ggplot2::xlab(metric)
      }
      if (plot == "boxplot"){
        p <- ggplot2::ggplot(metrics_df, ggplot2::aes(x = group, y = .data[[metric]], fill = group)) +
          ggplot2::geom_boxplot() +
          ggplot2::scale_fill_manual(values=colors) +
          #ggplot2::scale_y_continuous(trans='log10') +
          ggplot2::theme_classic() +
          ggplot2::theme(text = ggplot2::element_text(size = text.size),
                         legend.position = "none") +
          ggplot2::ylab(metric)

        #Add significance to the plot
        if(significance){
          #Get the unique combinations of groups if there are more than 2 clusters
          if(length(unique(metrics_df$group)) > 2){
            #Get the unique combinations
            combinations <- combinat::combn(unique(metrics_df$group), 2)
            combinations_list <- split(combinations, col(combinations))
          }else{combinations_list <- list(unique(metrics_df$group))}
          #Add to the existing plot
          p <- p + ggsignif::geom_signif(comparisons=combinations_list, step_increase = 0.1, test = "t.test")
        }

        output_plot[[metric]] <- p
      }
    }
  }
  return(output_plot)
}

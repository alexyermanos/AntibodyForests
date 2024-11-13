#' Function to calculate metrics for each tree in an AntibodyForests-object
#' @description Function to calculate metrics for each tree in an AntibodyForests-object
#' @param input AntibodyForests-object(s), output from AntibodyForests()
#' @param min.nodes The minimum number of nodes in a tree to calculate metrics (including the germline).
#' @param multiple.objects If TRUE: input should contain multiple AntibodyForests-objects (default FALSE)
#' @param metrics The metrics to be calculated (default mean.depth and nr.nodes)
#' 'nr.nodes'         : The total number of nodes
#' 'nr.cells'         : The total number of cells in this clonotype
#' 'mean.depth'       : Mean of the number of edges connecting each node to the germline
#' 'mean.edge.length' : Mean of the edge lengths between each node and the germline
#' 'group.node.depth'      : Mean of the number of edges connecting each node per group (node.features of the AntibodyForests-object) to the germline. (default FALSE)
#' 'group.edge.length'    : Mean of the sum of edge length of the shortest path between germline and nodes per group (node.features of the AntibodyForests-object) 
#' 'sackin.index'     : Sum of the number of nodes between each node and the germline
#' 'spectral.density' : Metrics of the spectral density profiles (calculated with package RPANDA)
#'    - peakedness            : Tree balance
#'    - asymmetry             : Shallow or deep branching events
#'    - principal eigenvalue  : Phylogenetic diversity
#'    - modalities            : The number of different structures within the tree
#' @param node.feature The node feature to be used for the group.depth metric.
#' @param group.node.feature The groups in the node feature to be plotted. Set to NA if all features should displayed. (default NA)
#' @param parallel If TRUE, the metric calculations are parallelized (default FALSE)
#' @param num.cores Number of cores to be used when parallel = TRUE. (Defaults to all available cores - 1)
#' @param output.format The format of the output. If set to "dataframe", a dataframe is returned. If set to "AntibodyForests", the metrics are added to the AntibodyForests-object. (default "dataframe")
#' @return Returns either a dataframe where the rows are trees and the columns are metrics or an AntibodyForests-object with the metrics added to trees
#' @export

AntibodyForests_metrics <- function(input,
                                    min.nodes,
                                    node.feature,
                                    group.node.feature,
                                    multiple.objects,
                                    metrics,
                                    parallel,
                                    num.cores,
                                    output.format){
  #Stop when no input is provided
  if(missing(input)){stop("Please provide a valid input object.")}
  #Set defaults
  if(missing(multiple.objects)){multiple.objects = F}
  if(missing(min.nodes)){min.nodes <- 1}
  if(missing(metrics)){metrics <- c("mean.depth", "nr.nodes")}
  if(missing(parallel)){parallel <- FALSE}
  if(missing(group.node.feature)){group.node.feature = NA}
  if(missing(output.format)){output.format = "dataframe"}
  #Check if the input is in the correct format.
  #If multiple.objects is TRUE, multiple AntibodyForests-objects should be in the input list, where the third item in the nested AntibodyForest-object should be of class "igraph"
  if((multiple.objects == F && class(input[[1]][[1]][["igraph"]]) != "igraph")|| (multiple.objects == T && class(input[[1]][[1]][[1]][["igraph"]]) != "igraph")){
    stop("The input is not in the correct AntibodyForests-object format.")}
  # If 'parallel' is set to TRUE but 'num.cores' is not specified, the number of cores is set to all available cores - 1
  if(parallel == TRUE && missing(num.cores)){num.cores <- parallel::detectCores() -1}
  #If metric is group.node.depth or group.edge.length, node.features should be provided
  if("group.node.depth" %in% metrics || "group.edge.length" %in% metrics){
    if(missing(node.feature)){stop("Please provide node features.")}
  }
  # If 'output.format' is not "AntibodyForests" or "dataframe", an error is thrown
  if(!(output.format %in% c("AntibodyForests", "dataframe"))){stop("Please provide a valid output format ('dataframe' or 'AntibodyForests').")}
  
  #Functions to calculate metrics
  calculate_mean_depth <- function(tree, nodes){
    #Get the shortest paths between each node and the germline
    paths <- igraph::shortest_paths(tree, from = "germline", to = nodes, output = "both")
    #Take the mean of the number of edges on each path
    depth <- mean(unlist(lapply(paths$epath, length)))
    return(depth)
  }
  
  calculate_mean_edge_length <- function(tree, nodes){
    #Get the total length of shortest paths between each node and the germline
    distance <- igraph::distances(tree, v = "germline", to = nodes, algorithm = "dijkstra",
                                    weights = as.numeric(igraph::edge_attr(tree)$edge.length))
    #Take the mean of these distances
    distance <- mean(distance)
    return(distance)
  }
    
  
  
  calculate_sackin_index <- function(tree){
    #Get the shortest paths between each node and the germline
    paths <- igraph::shortest_paths(tree, from = "germline", output = "both")
    #Sum the number of nodes in the paths
    depth <- sum(unlist(lapply(paths$vpath, length)))
    return(depth)
  }
  
  # calculate_colless_number <- function(tree){
  #   #transform igraph network into bifurcating phylo tree
  #   phylo_tree <- AntibodyForests_phylo(tree, solve_multichotomies = T)
  #   #calculate the colless number
  #   colless <- phyloTop::colless.phylo(phylo_tree, normalise = T)
  #   return(colless)
  # }
  
  calculate_spectral_density <- function(tree){
    #transform igraph network into bifurcating phylo tree
    phylo_tree <- AntibodyForests_phylo(tree, solve_multichotomies = F)
    #Calculate the spectral density of the tree
    sd <- RPANDA::spectR(phylo_tree, meth = "normal")
    return(sd)
  }
  
  #Calculate the metrics for a clonotype
  calculate_metrics <- function(clonotype, min.nodes, metrics){
    
    if (igraph::vcount(clonotype$igraph) >= min.nodes){
      #Create empty vector to store metrics
      metrics_vector <- c()
      
      if ("mean.depth" %in% metrics){
        #Calculate the mean depth for all nodes except the germline
        depth <- calculate_mean_depth(clonotype$igraph, 
                                      nodes = igraph::V(clonotype$igraph)[names(igraph::V(clonotype$igraph)) != "germline"])
        #Add to the metrics vector
        metrics_vector["mean.depth"] <- depth
      }
      
      
      if ("mean.edge.length" %in% metrics){
        mean_edge_length <- calculate_mean_edge_length(clonotype$igraph, 
                                                       nodes = igraph::V(clonotype$igraph)[names(igraph::V(clonotype$igraph)) != "germline"])
        #Add to the metrics vector
        metrics_vector["mean.edge.length"] <- mean_edge_length
      }
      
      
      if ("sackin.index" %in% metrics){
        si <- calculate_sackin_index(clonotype$igraph)
        #Add to the metrics vector
        metrics_vector["sackin.index"] <- si
      }
      
      if ("spectral.density" %in% metrics){
        if (igraph::vcount(clonotype$igraph) > 3){
          spectr <- calculate_spectral_density(clonotype$igraph)
          #Add to the metrics vector
          metrics_vector["spectral.peakedness"] <- spectr$peakedness
          metrics_vector["spectral.asymmetry"] <- spectr$asymmetry
          metrics_vector["spectral.principal.eigenvalue"] <- spectr$principal_eigenvalue
          metrics_vector["modalities"] <- spectr$eigengap
        }else {
          warning("Tree needs at least 2 nodes (additional to germline) to calculate spectral density.")
          #Add NA to the metrics vector
          metrics_vector["spectral.peakedness"] <- NA
          metrics_vector["spectral.asymmetry"] <- NA
          metrics_vector["spectral.principal.eigenvalue"] <- NA
          metrics_vector["modalities"] <-NA
        }
      }
      
      if ("group.node.depth" %in% metrics){
        #Get all the unique elements in this group
        if (multiple.objects){
          groups <- unique(unlist(lapply(input[[1]],function(a){lapply(a,function(b){lapply(b$nodes, function(c){c[[node.feature]]})})})))
        }else{groups <- unique(unlist(lapply(input,function(a){lapply(a,function(b){lapply(b$nodes, function(c){c[[node.feature]]})})})))}
        
        #If group.node.features is specified, check wether they are all present in the node features of the AntibodyForests-object.
        if(!all(is.na(group.node.feature))){
          if(all(group.node.feature %in% groups) == FALSE){stop("The groups are not in the node features of the AntibodyForests-object.")
          }else{groups <- group.node.feature}}

        for (group in groups){
          #Take the nodes have a cell of this group
          nodes <- names(which(lapply(clonotype$nodes, function(x){group %in% x[node.feature]}) == TRUE))
          if (identical(nodes, character(0))){
            #Add NA to the metrics vector
            metrics_vector[paste0(group,".node.depth")] <- NA
          }else{
            #Calcute the mean depth for the nodes that have this group
            depth <- calculate_mean_depth(clonotype$igraph, nodes = igraph::V(clonotype$igraph)[nodes])
            #Add to the metrics vector
            metrics_vector[paste0(group,".node.depth")] <- depth
          }
        }
      }
      
      if ("group.edge.length" %in% metrics){
        #Get all the unique elements in this group
        if (multiple.objects){
          groups <- unique(unlist(lapply(input[[1]],function(a){lapply(a,function(b){lapply(b$nodes, function(c){c[[node.feature]]})})})))
        }else{groups <- unique(unlist(lapply(input,function(a){lapply(a,function(b){lapply(b$nodes, function(c){c[[node.feature]]})})})))}

        #If group.node.features is specified, check wether they are all present in the node features of the AntibodyForests-object.
        if(!all(is.na(group.node.feature))){
          if(all(group.node.feature %in% groups) == FALSE){stop("The groups are not in the node features of the AntibodyForests-object.")
          }else{groups <- group.node.feature}}
          
        for (group in groups){
          #Take the nodes have a cell of this group
          nodes <- names(which(lapply(clonotype$nodes, function(x){group %in% x[node.feature]}) == TRUE))
          if (identical(nodes, character(0))){
            #Add NA to the metrics vector
            metrics_vector[paste0(group,".edge.length")] <- NA
          }else{
            #Calcute the mean edge length for the nodes that have this group
            el <- calculate_mean_edge_length(clonotype$igraph, nodes = igraph::V(clonotype$igraph)[nodes])
            #Add to the metrics vector
            metrics_vector[paste0(group,".edge.length")] <- el
          }
        }
      }

      #Total number of nodes
      if ("nr.nodes" %in% metrics){
        nodes <- igraph::vcount(clonotype$igraph)
        metrics_vector["nr.nodes"] <- nodes
      }
      
      #Total number of cells
      if ("nr.cells" %in% metrics){
        cells <- sum(unlist(lapply(clonotype$nodes, function(x){length(x$barcodes)})))
        metrics_vector["nr.cells"] <- cells
      }
      
      return(metrics_vector)
    }else{
      return(NA)
    }
  }
  


  # If 'parallel' is set to TRUE, the metric calculation is parallelized across the clonotypes
  if(parallel){
    # Retrieve the operating system
    operating_system <- Sys.info()[['sysname']]
    # If the operating system is Linux or Darwin, 'mclapply' is used for parallelization
    if(operating_system %in% c('Linux', 'Darwin')){
      #Go over each tree in the AntibodyForests object and add the metrics to the AntibodyForests-object
      if(output.format == "AntibodyForests"){
        antibodyforests_object <- parallel::mclapply(mc.cores = num.cores,input, function(sample){
          parallel::mclapply(mc.cores = num.cores,sample, function(clonotype){
            metrics <- calculate_metrics(clonotype, min.nodes, metrics)
            clonotype$metrics <- metrics
            return(clonotype)
          })
        })
        antibodyforests_object <- base::structure(antibodyforests_object, class = "AntibodyForests")
        return(antibodyforests_object)
      }
      #Go over each tree in the AntibodyForests object and create a metric dataframe
      if(output.format == "dataframe"){
        metric_df <- t(as.data.frame(parallel::mclapply(mc.cores = num.cores,input, function(sample){
          parallel::mclapply(mc.cores = num.cores,sample, function(clonotype){
            calculate_metrics(clonotype, min.nodes, metrics)
          })
        })))
        #Only keep rows that have enough nodes (nr_nodes != NA)
        metric_df <- metric_df[which(is.na(metric_df[,ncol(metric_df)]) == FALSE),]
        #Transform to dataframe
        metric_df <- as.data.frame(metric_df)
        #Add sample names to the dataframe
        metric_df$sample <- gsub(".clonotype(.*)$", "", rownames(metric_df))
        return(metric_df)
      }
    }
    # If the operating system is Windows, "parLapply" is used for parallelization
    if(operating_system == "Windows"){
      # Create cluster
      cluster <- parallel::makeCluster(num.cores)
      
      # If the operating system is Linux or Darwin, 'mclapply' is used for parallelization
      if(output.format == "AntibodyForests"){
        #Go over each tree in the AntibodyForests object and add the metrics to the AntibodyForests-object
        antibodyforests_object <- parallel::parLapply(cluster,input, function(sample){
          parallel::parLapply(cluster,sample, function(clonotype){
            metrics <- calculate_metrics(clonotype, min.nodes, metrics)
            clonotype$metrics <- metrics
            return(clonotype)
          })
        })
        antibodyforests_object <- base::structure(antibodyforests_object, class = "AntibodyForests")
        return(antibodyforests_object)
      }
      #Go over each tree in the AntibodyForests object and create a metric dataframe
      if(output.format == "dataframe"){
        metric_df <- t(as.data.frame(parallel::parLapply(cluster,input, function(sample){
          parallel::parLapply(cluster,sample, function(clonotype){
            calculate_metrics(clonotype, min.nodes, metrics)
          })
        })))
        #Only keep rows that have enough nodes (nr_nodes != NA)
        metric_df <- metric_df[which(is.na(metric_df[,ncol(metric_df)]) == FALSE),]
        #Transform to dataframe
        metric_df <- as.data.frame(metric_df)
        #Add sample names to the dataframe
        metric_df$sample <- gsub(".clonotype(.*)$", "", rownames(metric_df))
        return(metric_df)
      }
    }
  }
  
  
  # If 'parallel' is set to FALSE, the network inference is not parallelized
  if(!parallel){
    #Go over each tree in the AntibodyForests object and add the metrics to the AntibodyForests-object
    if(output.format == "AntibodyForests"){
      antibodyforests_object <- lapply(input, function(sample){
        lapply(sample, function(clonotype){
          metrics <- calculate_metrics(clonotype, min.nodes, metrics)
          clonotype$metrics <- metrics
          return(clonotype)
        })
      })
      antibodyforests_object <- base::structure(antibodyforests_object, class = "AntibodyForests")
      return(antibodyforests_object)
    }
    # Go over each tree in the AntibodyForests object and create a metric dataframe
    if(output.format == "dataframe"){
      metric_df <- t(as.data.frame(lapply(input, function(sample){
        lapply(sample, function(clonotype){
          calculate_metrics(clonotype, min.nodes, metrics)
        })
      })))
      #Only keep rows that have enough nodes (nr_nodes != NA)
      metric_df <- metric_df[which(is.na(metric_df[,ncol(metric_df)]) == FALSE),]
      #Transform to dataframe
      metric_df <- as.data.frame(metric_df)
      #Add sample names to the dataframe
      metric_df$sample <- gsub(".clonotype(.*)$", "", rownames(metric_df))
      return(metric_df)
    }
    
  }
  
  
}






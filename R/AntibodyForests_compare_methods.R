#' Function to compare trees created with different algorithms from the same clonotype.
#' @description Function to compare different trees from the same clonotype to compare various graph construction and phylogenetic reconstruction methods.
#' @param input A list of AntibodyForests-objects as output from the function AntibodyForests(). These objects should contain the same samples/clonotypes. For easy interpretation of the results, please name the objects in the list according to their tree-construction method.
#' @param min.nodes The minimum number of nodes in a tree to include in the comparison, this includes the germline. Default is 2 (this includes all trees).
#' @param include.average If TRUE, the average distance matrix and visualizations between the trees is included in the output (default FALSE)
#' @param distance.method The method to calculate the distance between trees (default euclidean)
#' 'euclidean' : Euclidean distance between the depth of each node in the tree
#' 'GBLD'      : Generalized Branch Length Distance, derived from Mahsa Farnia & Nadia Tahiri, Algorithms Mol Biol 19, 22 (2024). https://doi.org/10.1186/s13015-024-00267-1
#' @param depth If distance.methods is 'euclidean', method to calculate the germline-to-node depth (default edge.count)
#' 'edge.count'   : The number of edges between each node and the germline
#' 'edge.length'  : The sum of edge lengths between each node and the germline
#' @param clustering.method Method to cluster trees (default NULL)
#' NULL             : No clustering
#' 'mediods'        : Clustering based on the k-mediods method. The number of clusters is estimated based on the optimum average silhouette.
#' @param visualization.methods The methods to analyze similarity (default NULL)
#' NULL             : No visualization
#' 'PCA'            : Scatterplot of the first two principal components.
#' 'MDS'            : Scatterplot of the first two dimensions using multidimensional scaling.
#' "heatmap'        : Heatmap of the distance
#' @param parallel If TRUE, the depth calculations are parallelized across clonotypes (default FALSE)
#' @param num.cores Number of cores to be used when parallel = TRUE. (Defaults to all available cores - 1)
#' @return A list with all clonotypes that pass the min.nodes threshold including the distance matrix, possible clustering and visualization
#' @export

AntibodyForests_compare_methods <- function(input,
                                    min.nodes,
                                    include.average,
                                    distance.method,
                                    depth,
                                    clustering.method,
                                    visualization.methods,
                                    parallel,
                                    num.cores){
  
  #1. Set defaults and check for missing input
  if(missing(input)){stop("Please provide a valid input object.")}
  if(missing(min.nodes)){min.nodes = 2}
  if(missing(distance.method)){distance.method = "euclidean"}
  if(missing(depth)){depth = "edge.count"}
  if(missing(visualization.methods)){visualization.methods = NULL}
  if(missing(clustering.method)){clustering.method = NULL}
  if(missing(parallel)){parallel <- FALSE}
  if(missing(include.average)){include.average <- FALSE}
  # If 'parallel' is set to TRUE but 'num.cores' is not specified, the number of cores is set to all available cores - 1
  if(parallel == TRUE && missing(num.cores)){num.cores <- parallel::detectCores() -1}
  
  #2. Check if the input is correct
  if (class(input[[1]][[1]][[1]][["igraph"]]) != "igraph"){stop("The input is not in the correct format.")}
  if(!(all(equal <- sapply(input, function(x) all.equal(names(x), names(input[[1]])))))){stop("The input list does not contain the same clonotypes.")}
  if(min.nodes > max(unlist(lapply(input, function(x){lapply(x,function(y){lapply(y, function(y){lapply(y$nodes, function(z){z$size})})})})))){
    stop("min.nodes is larger than the biggest clonotype.")}
  if(!(depth %in% c("edge.count", "edge.length"))){stop("Unvalid depth argument.")}
  if(!(distance.method %in% c("euclidean", "GBLD"))){stop("Unvalid distance.method argument.")}
  if(!(is.null(clustering.method)) && clustering.method != "mediods"){stop("Unvalid clustering method.")}
  if(!(is.null(visualization.methods)) && all(!(visualization.methods %in% c("PCA", "MDS", "heatmap")))){stop("Unvalid visualization methods.")}
  
  
  #3. Define functions
  #Calculate the depths of each node in each tree
  calculate_depth_list <- function(input, min.nodes, parallel, num.cores){
    # If 'parallel' is set to TRUE, the depth calculation is parallelized across the trees
    if(parallel){
      # Retrieve the operating system
      operating_system <- Sys.info()[['sysname']]
      # If the operating system is Linux or Darwin, 'mclapply' is used for parallelization
      if(operating_system %in% c('Linux', 'Darwin')){
        #Go over each tree in clonotype and create a depth list
        depth_list <- parallel::mclapply(mc.cores = num.cores, input, function(object){
          parallel::mclapply(mc.cores = num.cores, object, function(sample){
            parallel::mclapply(mc.cores = num.cores, sample, function(clonotype){
              #Only calculate depth for trees with at least min.nodes
              if (is.null(clonotype$igraph)){return(NA)}
              else{
                if(igraph::vcount(clonotype$igraph) >= min.nodes){
                  #Calculate the depth for all nodes except the germline
                  depth_vector <- calculate_depth(clonotype$igraph,
                                                  nodes = igraph::V(clonotype$igraph)[names(igraph::V(clonotype$igraph)) != "germline"],
                                                  depth)
                  return(depth_vector)
                }else{
                  return(NA)
                }
              }
              
            })
          })
        })
      }
      # If the operating system is Windows, "parLapply" is used for parallelization
      if(operating_system == "Windows"){
        # Create cluster
        cluster <- parallel::makeCluster(num.cores)
        #Go over each tree in each of the AntibodyForests objects and create a list of depths
        depth_list <- parallel::parLapply(cluster, input, function(object){
          parallel::parLapply(cluster, object, function(sample){
            parallel::parLapply(cluster, sample, function(clonotype){
              #Only calculate depth for trees with at least min.nodes
              if (is.null(clonotype$igraph)){return(NA)}
              else{
                if(igraph::vcount(clonotype$igraph) >= min.nodes){
                  #Calculate the depth for all nodes except the germline
                  depth_vector <- calculate_depth(clonotype$igraph,
                                                  nodes = igraph::V(clonotype$igraph)[names(igraph::V(clonotype$igraph)) != "germline"],
                                                  depth)
                  return(depth_vector)
                }else{
                  return(NA)
                }
              }
            })
          })
        })
        # Stop cluster
        parallel::stopCluster(cluster)
      }
    }
    # If 'parallel' is set to FALSE, the network inference is not parallelized
    if(!parallel){
      #Go over each tree in each of the AntibodyForests objects and create a list of depths
      depth_list <- lapply(input, function(object){
        lapply(object, function(sample){
          lapply(sample, function(clonotype){
            #Only calculate depth for trees with at least min.nodes
            if (is.null(clonotype$igraph)){
              return(NA)}
            else{
              if(igraph::vcount(clonotype$igraph) >= min.nodes){
                #Calculate the depth for all nodes except the germline
                depth_vector <- calculate_depth(clonotype$igraph,
                                                nodes = igraph::V(clonotype$igraph)[names(igraph::V(clonotype$igraph)) != "germline"],
                                                depth)
                return(depth_vector)
              }else{
                return(NA)
              }
            }
          })
        })
      })
    }
    
    #Remove NA (clonotypes with less then min.nodes nodes)
    depth_list <- remove_na_entries(depth_list)
    return(depth_list)
  }
  
  #Calculate the number of edges between certain nodes and the germline of a single tree
  calculate_depth <- function(tree, nodes, depth){
    if (depth == "edge.count"){
      #Get the shortest paths between each node and the germline
      paths <- igraph::shortest_paths(tree, from = "germline", to = nodes, output = "both")
      #Set names to the list of vpath
      names(paths$epath) <- names(unlist(lapply(paths$vpath, function(x){tail(x,n=1)})))
      #Reorder according to node number
      paths$epath <- paths$epath[paste0("node",sort(as.numeric(stringr::str_sub(names(paths$epath), start = 5))))]
      #Get the number of edges along the shortes paths
      depths <- unlist(lapply(paths$epath, length))
      #Set names to the vector of edge counts
      names(depths) <- names(paths$epath)
    }
    if (depth == "edge.length"){
      #Get the total length of shortest paths between each node and the germline
      depths <- igraph::distances(tree, v = "germline", to = nodes, algorithm = "dijkstra",
                                    weights = as.numeric(igraph::edge_attr(tree)$edge.length))
      #Reorder according to node number
      depths <- depths[,paste0("node",sort(stringr::str_sub(colnames(depths), start = 5)))]
    }
    #Return the named vector of depths per node
    return(depths)
  }
  
  #Remove entries from (nested) lists that contain NA
  remove_na_entries <- function(lst) {
    if (!is.list(lst)) {
      return(lst)
    } else {
      cleaned <- lapply(lst, remove_na_entries)
      cleaned <- cleaned[sapply(cleaned, function(x) !all(is.na(x)))]
      return(cleaned)
    }
  }
  
  #Calculate the euclidean distance between node depth for each clonotype in the depth_list
  calculate_euclidean <- function(depth_list){
    distance_list <- list()
    samples <- unique(unlist(lapply(depth_list, names)))
    for (sample in samples){
      #Get the clonotype names in this sample
      clonotypes <- names(depth_list[[1]][[sample]])
      for (clonotype in clonotypes){
        #Create a dataframe where the rows are the tree construction methods and the columns are the depth per node
        nodes <- names(depth_list[[1]][[sample]][[clonotype]])
        depth_df <- matrix(ncol = length(nodes), nrow = 0)
        colnames(depth_df) <- nodes
        
        skip = F
        for (method in names(depth_list)){
          #If the number of columns is not the same, skip this clonotype
          if(ncol(depth_df) != length(depth_list[[method]][[sample]][[clonotype]])){skip = T;break}
          #Add the depth to the dataframe
          depth_df <- rbind(depth_df, depth_list[[method]][[sample]][[clonotype]])
        }
        if(skip){next}
        rownames(depth_df) <- names(depth_list)
        
        #Calculate the euclidean distance between the tree methods for this clonotype
        euclidean_matrix <- as.matrix(stats::dist(depth_df, method = "euclidean", diag = T, upper = T))
        
        #Add distance matrix to the list
        distance_list[[paste0(sample,".",clonotype)]] <- euclidean_matrix
      }
    }
    return(distance_list)
  }
  
  cluster_mediods <- function(df){
    distance_matrix <- df
    #Define the max number of clusters
    max_cluster <- dim(as.matrix(distance_matrix))[1]-1
    #Perform clustering
    mediods <- fpc::pamk(distance_matrix,krange=1:max_cluster)
    clusters <- mediods$pamobject$clustering
    #Assign the clusters to the tree names
    names(clusters) <- rownames(df)
    return(clusters)
  }
  
  #Calculate principle components
  calculate_PC <- function(df, to.scale){
    names <- rownames(as.data.frame(df))
    #Run a PCA and save the PCs in a dataframe
    pca_results <- as.data.frame(stats::prcomp(df, scale. = to.scale)$x)
    #Keep the first two PCs
    pca_results <- pca_results[,1:2]
    colnames(pca_results) <- c("Dim1", "Dim2")
    #Add tree names to the dataframe
    pca_results$tree <- names
    return(pca_results)
  }
  
  #Multidimensional scaling
  calculate_MDS <- function(df){
    names <- rownames(as.data.frame(df))
    #Compute classical metric multidimensional scaling
    results <- as.data.frame(stats::cmdscale(df))
    colnames(results) <- c("Dim1", "Dim2")
    #Add tree names to the dataframe
    results$tree <- names
    
    return(results)
  }
  
  #Plot the first two dimensions of a PCA or MDS
  plot <- function(df, color, name){
    p <- ggplot2::ggplot(as.data.frame(df), ggplot2::aes(x=Dim1,y=Dim2, color=as.factor(.data[[color]]))) +
      ggplot2::geom_point(size=5) +
      ggrepel::geom_label_repel(ggplot2::aes(label = tree))+
      ggplot2::theme_minimal() +
      ggplot2::theme(text = ggplot2::element_text(size = 20)) +
      ggplot2::ggtitle(name)
    return(p)
  }
  
  #Plot a heatmap of the distance matrix
  heatmap <- function(df, clusters){
    #If there are clusters make clustered heatmap
    if(!(is.null(clusters))){
      #Make df of the clusters
      clusters <- as.data.frame(clusters)
      colnames(clusters) <- "cluster"
      #rownames(clusters) <- labels(df)[[1]]
      clusters$cluster <- as.factor(clusters$cluster)
      
      #print(labels(df)[[1]])
      #print(clusters)
      #Plot the clustered heatmap
      p <- pheatmap::pheatmap(as.matrix(df),
                         annotation_col = clusters)
    }
    #If there are no clusters
    else{
      p <- pheatmap::pheatmap(as.matrix(df))}
    
    return(p)
  }
  
  #For GBLD, Code is derived from https://github.com/tahiri-lab/ClonalTreeClustering/blob/main/src/Python/GBLD_Metric_Final.ipynb 
  #Algorithms Mol Biol 19, 22 (2024). https://doi.org/10.1186/s13015-024-00267-1
  #Calculate the distance (sum of branch lengths) between terminal nodes
  calculate_distances <- function(tree_igraph, method){
    #Get the terminal node(s)
    terminals <- igraph::degree(tree_igraph)[which(igraph::degree(tree_igraph) == 1 & names(igraph::degree(tree_igraph)) != "germline")]
    #Create empty distance matrix
    distances = matrix(0, nrow = length(terminals), ncol = length(terminals))
    #Calculate the distance between terminal nodes and add to the distance matrix
    for (i in seq(1:length(terminals))){
      for (j in seq(1:length(terminals))){
        if (i != j){
          distances[i, j] <- igraph::distances(tree_igraph, v = names(terminals)[i], to = names(terminals)[j], algorithm = "dijkstra",
                                               weights = as.numeric(igraph::edge_attr(tree_igraph)$edge.length))[1,1]
        }
      }
    }
    #Add the node names to the distance matrix
    colnames(distances) <- paste0(method, "_", names(terminals))
    rownames(distances) <- paste0(method, "_", names(terminals))
    return(distances)
  }
  #Normalize a matrix
  normalize_matrix <- function(matrix){
    if (min(matrix) == max(matrix)){return(matrix(0, nrow(matrix), ncol(matrix)))}
    return((matrix - min(matrix))/(max(matrix) - min(matrix)))
  }
  #Normalize the weights
  normalize_weights <- function(weights){
    weights <- weights[[1]]
    non_zero_weights = weights[weights != 0]
    if(is.null(non_zero_weights)){return(weights)}
    else{
      min_weight = min(non_zero_weights)
      max_weight = max(non_zero_weights)
      if (min_weight == max_weight){
        weights[weights != 0] = 1
      }else{
        weights[weights != 0] = (weights[weights != 0] - min_weight)/(max_weight - min_weight)
      }
      return(weights)
    }
  }
  #Calculate the distance between two trees
  calculate_D <- function(dist1, dist2, max_nodes){
    #Make dist1 and dist2 the same size
    if (ncol(dist1) < max_nodes){
      dist1 <- cbind(dist1, matrix(rep(rep(0, nrow(dist1)), (max_nodes - nrow(dist1))), nrow = nrow(dist1)))
      dist1 <- rbind(dist1, matrix(rep(rep(0, ncol(dist1)), (max_nodes - nrow(dist1))), ncol = ncol(dist1)))
    }
    if (ncol(dist2) < max_nodes){
      dist2 <- cbind(dist2, matrix(rep(rep(0, nrow(dist2)), (max_nodes - nrow(dist2))), nrow = nrow(dist2)))
      dist2 <- rbind(dist2, matrix(rep(rep(0, ncol(dist2)), (max_nodes - nrow(dist2))), ncol = ncol(dist2)))
    }
    #Calculate the distance between the trees
    diff <- dist1 - dist2
    diff_squared <- diff ^ 2
    sum_diff_squared <- sum(diff_squared)
    root_sum_squared <- sqrt(sum_diff_squared)
    
    return(root_sum_squared / max_nodes)
  }
  #Calculate the penalty of uncommon nodes between two trees
  calculate_penalty <- function(tree1, tree2){
    common_nodes = sum(tree1 != 0 & tree2 != 0)
    total_nodes = sum(tree1 != 0 | tree2 != 0)
    return(1 - (common_nodes / total_nodes))
  }
  #Calculate the GBLD distance between trees constructed with different methods
  calculate_GBLD_methods <- function(AntibodyForests_object, min.nodes){
    all_weights = list()
    all_distances = list()
    all_seq_names = c()
    
    tree_count = 0
    
    for(method in names(AntibodyForests_object)){
      tree_igraph <- AntibodyForests_object[[method]][['igraph']]
      #Only keep trees with a minimum number of nodes (min.nodes)
      if (igraph::vcount(tree_igraph) >= min.nodes){
        #Keep track of the number of trees
        tree_count = tree_count + 1
        
        #Set the tree names
        tree_label = method
        
        #Get the node sizes
        weights = lapply(AntibodyForests_object[[method]][['nodes']], function(x){if (is.null(x$size)){return(1)}else{return(x$size)}})
        #names(weights) <- paste0(method, "_", names(weights))
        
        for(node in names(weights)){
          if(!(node %in% names(all_weights))){
            all_weights[[node]] = rep(0, (tree_count - 1))
          }
          all_weights[[node]] = c(all_weights[[node]], weights[[node]])
        }
        for (node in names(all_weights)){
          if (length(all_weights[[node]]) < tree_count){
            all_weights[[node]] <- c(all_weights[[node]], rep(0, (tree_count - length(all_weights[[node]]))))
          }
        }
        
        #Calculate the distance (sum of branch lengths) between terminal nodes
        distances <- calculate_distances(tree_igraph, method)
        normalized_distances <- normalize_matrix(distances)
        #Store te normalized distances and terminal nodes
        all_distances[[tree_label]] <- normalized_distances
        all_seq_names <- c(all_seq_names, colnames(normalized_distances))
      }
    }
    
    all_nodes <- sort(names(all_weights))
    original_weights_matrix <- matrix(unlist(all_weights), nrow = length(all_nodes), byrow = T)
    nodes_per_tree <- colSums(original_weights_matrix != 0)
    
    normalized_weights_matrix <- matrix(0, nrow = nrow(original_weights_matrix), ncol = ncol(original_weights_matrix))
    i = 0
    for (node in all_nodes){
      i = i + 1
      normalized_weights <- normalize_weights(all_weights[node])
      for (j in seq(1:length(all_distances))){
        normalized_weights_matrix[i,j] <- normalized_weights[j]
      }
    }
    
    num_trees = length(all_distances)
    diff_columns = list()
    normalized_diff_sums <- c()
    
    for (i in seq(from = 1, to = num_trees, by = 1)){
      if (i != num_trees){
        for (j in seq(from = i+1, to = num_trees, by = 1)){
          diff <- abs(normalized_weights_matrix[,i] - normalized_weights_matrix[,j])
          diff_columns[[paste0("Diff ", names(all_distances)[i], "-", names(all_distances)[j])]] <- diff
          
          max_nodes <- max(nodes_per_tree[i], nodes_per_tree[j])
          normalized_sum = sum(diff)/max_nodes
          normalized_diff_sums <- c(normalized_diff_sums, normalized_sum)
        }
      }
    }
    
    #Create empty distance matrix
    distance_matrix <- matrix(0, nrow = num_trees, ncol = num_trees)
    #Change dashes and underscores to dots to match the tree names in downstream AntibodyForests functions
    colnames(distance_matrix) <- gsub("[-_]", ".", names(all_distances)) 
    rownames(distance_matrix) <- gsub("[-_]", ".", names(all_distances))
    
    W_values = c()
    for (i in seq(from = 1, to = num_trees, by = 1)){
      if (i != num_trees){
        for (j in seq(from = i+1, to = num_trees, by = 1)){
          #Calculate the distance between the trees
          max_nodes <- max(nodes_per_tree[i], nodes_per_tree[j])
          D = calculate_D(all_distances[[i]], all_distances[[j]], max_nodes)
          
          #Calculate the penalty of uncommon nodes
          P = calculate_penalty(normalized_weights_matrix[,i], normalized_weights_matrix[,j])
          
          #Get the weights of the nodes (abundance)
          W = normalized_diff_sums[length(W_values)+1]
          W_values <- c(W_values, W)
          
          #Calculate the GBLD distance between the trees
          GBLD = P + W + D
          
          #Add the GBLD distance to the distance matrix
          distance_matrix[i, j] <- GBLD
          distance_matrix[j, i] <- GBLD
        }
      }
    }
    
    return(distance_matrix)
  }
  
  #Calculate the GBLD distance between trees for each clonotype per sample
  calculate_GBLD_list <- function(input, min.nodes){
    distance_list <- list()
    samples <- unique(unlist(lapply(input, names)))
    methods <- names(input)
    for (sample in samples){
      #Get the clonotype names in this sample
      clonotypes <- names(input[[1]][[sample]])
      for (clonotype in clonotypes){
        #Create antibodyforests object for this clonotype in every method
        clonotype_antibodyforests <- lapply(names(input), function(x){return(input[[x]][[sample]][[clonotype]])})
        names(clonotype_antibodyforests) <- paste0(methods)
        
        #Calculate the GBLD
        distance_matrix <- tryCatch({calculate_GBLD_methods(clonotype_antibodyforests, min.nodes = min.nodes)},
                                    error = function(e) {return(NA)})
        #Add to the distance list
        if(is.matrix(distance_matrix)){distance_list[[paste0(sample,".",clonotype)]] <- distance_matrix}
        }
    }
    return(distance_list)
    }
  


  #4. Calculate the distance between trees
  #Euclidean distance between node depths
  if(distance.method == "euclidean"){
    #Calculate the depth of each node in each tree
    depth_list <- calculate_depth_list(input, min.nodes, parallel, num.cores)
    #Calculate the euclidean distance between trees for each clonotype per sample
    output_list <- calculate_euclidean(depth_list)
  }
  #Generalized Branch Length Distance
  if(distance.method == "GBLD"){
    #Calculate the GBLD list between trees for each clonotype per sample
    output_list <- calculate_GBLD_list(input, min.nodes)
  }
  
  #5. Clustering and visualization for each clonotype
  output_list <- lapply(output_list, function(distance_matrix){
    #Create inner list
    temp_list <- list()
    #Add the distance matrix to this list
    temp_list[["distance.matrix"]] <- distance_matrix
    
    #Clustering
    if(!(is.null(clustering.method))){
      #K-mediods clustering
      if (clustering.method == "mediods"){
        #Get clusters
        clusters <- cluster_mediods(distance_matrix)
      }
      #Add to the list
      temp_list[["clusters"]] <- clusters
    }
    
    #Visualization
    if(!(is.null(visualization.methods))){
      #PCA analysis
      if("PCA" %in% visualization.methods){
        #Get PCA dimensions
        pca <- calculate_PC(distance_matrix, to.scale = F)
        #If there are clusters calculated
        if(!(is.null(clustering.method))){
          #Add clusters to the pca dataframe
          pca <- cbind(pca, "clusters" = clusters)
          #Create plot
          plot <- plot(pca, color = "clusters", name = "PCA")
        }
        #If there are no clusters
        else{
          #Create plot
          plot <- plot(pca, color = "tree", name = "PCA")
        }
        #Add to the list
        temp_list[["PCA"]] = plot
      }
      #MDS analysis
      if("MDS" %in% visualization.methods){
        #Only do MDS when there are more than 2 trees to compare
        if (ncol(as.matrix(distance_matrix)) > 2){
          #Get MDS dimensions
          mds <- calculate_MDS(distance_matrix)
          #If there are clusters calculated
          if(!(is.null(clustering.method))){
            #Add clusters to the pca dataframe
            mds <- cbind(mds, "clusters" = clusters)
            #Create plot
            plot <- plot(mds, color = "clusters", name = "MDS")
          }
          #If there are no clusters
          else{
            #Create plot
            plot <- plot(mds, color = "tree", name = "MDS")
          }
          #Add to the list
          temp_list[["MDS"]] = plot
        }else{
          temp_list[["MDS"]] = "Need at least 3 trees to compute MDS"
        }
      }
      #Heatmap
      if("heatmap" %in% visualization.methods){
        #If there are clusters calculated
        if(!(is.null(clustering.method))){
          #Create plot
          hm <- heatmap(distance_matrix, clusters = clusters)
        }
        #If there are no clusters
        else{
          #Create plot
          hm <- heatmap(distance_matrix, clusters = NULL)
        }
        #Add to the list
        temp_list[["Heatmap"]] <- hm
      }
      
    }
    
    return(temp_list)
  })
  
  #6. If include.average is TRUE, calculate the average distance matrix and visualizations
  if(include.average){
    temp_list <- list()
    #Calculate the average distance matrix over all clontoypes
    distance_list <- lapply(output_list, function(x) as.matrix(x$distance.matrix))
    mean_matrix <- Rcompadre::mat_mean(distance_list)
    temp_list[["distance.matrix"]] <- mean_matrix
    
    #Clustering
    if(!(is.null(clustering.method))){
      #K-mediods clustering
      if (clustering.method == "mediods"){
        #Get clusters
        clusters <- cluster_mediods(mean_matrix)
      }
      #Add to the list
      temp_list[["clusters"]] <- clusters
    }
    
    #Visualization
    if(!(is.null(visualization.methods))){
      #PCA analysis
      if("PCA" %in% visualization.methods){
        #Get PCA dimensions
        pca <- calculate_PC(mean_matrix, to.scale = F)
        #If there are clusters calculated
        if(!(is.null(clustering.method))){
          #Add clusters to the pca dataframe
          pca <- cbind(pca, "clusters" = clusters)
          #Create plot
          plot <- plot(pca, color = "clusters", name = "PCA")
        }
        #If there are no clusters
        else{
          #Create plot
          plot <- plot(pca, color = "tree", name = "PCA")
        }
        #Add to the list
        temp_list[["PCA"]] = plot
      }
      #MDS analysis
      if("MDS" %in% visualization.methods){
        #Only do MDS when there are more than 2 trees to compare
        if (ncol(as.matrix(mean_matrix)) > 2){
          #Get MDS dimensions
          mds <- calculate_MDS(mean_matrix)
          #If there are clusters calculated
          if(!(is.null(clustering.method))){
            #Add clusters to the pca dataframe
            mds <- cbind(mds, "clusters" = clusters)
            #Create plot
            plot <- plot(mds, color = "clusters", name = "MDS")
          }
          #If there are no clusters
          else{
            #Create plot
            plot <- plot(mds, color = "tree", name = "MDS")
          }
          #Add to the list
          temp_list[["MDS"]] = plot
        }else{
          temp_list[["MDS"]] = "Need at least 3 trees to compute MDS"
        }
      }
      #Heatmap
      if("heatmap" %in% visualization.methods){
        #If there are clusters calculated
        if(!(is.null(clustering.method))){
          #Create plot
          hm <- heatmap(as.dist(mean_matrix), clusters = clusters)
        }
        #If there are no clusters
        else{
          #Create plot
          hm <- heatmap(as.dist(mean_matrix), clusters = NULL)
        }
        #Add to the list
        temp_list[["Heatmap"]] <- hm
      }
    }
    #Add to the output list
    output_list[["average"]] <- temp_list
  }
  
  return(output_list)


}
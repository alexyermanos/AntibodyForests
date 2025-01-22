#' Calculate the GBLD distance between trees in an AntibodyForests object. Code is derived from https://github.com/tahiri-lab/ClonalTreeClustering/blob/main/src/Python/GBLD_Metric_Final.ipynb
#' Farnia, M., Tahiri, N. New generalized metric based on branch length distance to compare B cell lineage trees. Algorithms Mol Biol 19, 22 (2024). https://doi.org/10.1186/s13015-024-00267-1
#' @description Calculate the GBLD distance between trees in an AntibodyForests object. Code is derived from https://github.com/tahiri-lab/ClonalTreeClustering/blob/main/src/Python/GBLD_Metric_Final.ipynb
#' @param AntibodyForests_object AntibodyForests-object, output from AntibodyForests()
#' @param min.nodes - integer - The minimum number of nodes (including the germline) in a tree to include in the analysis. Default is 3.
#' @return A matrix with the GBLD distances between trees in the AntibodyForests object.
#' @export
#' @examples
#' GBLD_matrix <- calculate_GBLD(AntibodyForests_object = AntibodyForests::small_af)
#' GBLD_matrix[1:5, 1:5]

calculate_GBLD <- function(AntibodyForests_object,
                           min.nodes){

  if(missing(AntibodyForests_object)){stop("Please provide an AntibodyForests object")}
  if(missing(min.nodes)){min.nodes = 3}

  message("Calculating the GBLD can have a long runtime for large AntibodyForests-objects.")

  #Calculate the distance (sum of branch lengths) between terminal nodes
  calculate_distances <- function(tree_igraph, sample, clonotype){
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
    colnames(distances) <- paste0(sample, "_", clonotype, "_", names(terminals))
    rownames(distances) <- paste0(sample, "_", clonotype, "_", names(terminals))
    return(distances)
  }

  #Normalize a matrix
  normalize_matrix <- function(matrix){
    if (min(matrix) == max(matrix)){return(matrix(0, nrow(matrix), ncol(matrix)))}
    return((matrix - min(matrix))/(max(matrix) - min(matrix)))
  }

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

  all_weights = list()
  all_distances = list()
  all_seq_names = c()

  tree_count = 0

  for(sample in names(AntibodyForests_object)){
    for(clonotype in names(AntibodyForests_object[[sample]])){
      tree_igraph <- AntibodyForests_object[[sample]][[clonotype]][['igraph']]
      #Only keep trees with a minimum number of nodes (min.nodes)
      if (igraph::vcount(tree_igraph) >= min.nodes){
        #Keep track of the number of trees
        tree_count = tree_count + 1

        #Set the tree names
        tree_label = paste0(sample, "_", clonotype)

        #Get the node sizes
        weights = lapply(AntibodyForests_object[[sample]][[clonotype]][['nodes']], function(x){if (is.null(x$size)){return(1)}else{return(x$size)}})
        names(weights) <- paste0(sample, "_", clonotype, "_", names(weights))

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
        distances <- calculate_distances(tree_igraph, sample, clonotype)
        normalized_distances <- normalize_matrix(distances)
        #Store te normalized distances and terminal nodes
        all_distances[[tree_label]] <- normalized_distances
        all_seq_names <- c(all_seq_names, colnames(normalized_distances))
      }
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

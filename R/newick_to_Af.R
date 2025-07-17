#'Converts files with phylogenetic trees in newick format into an AntibodyForests object.
#'@description Converts files with phylogenetic trees in newick format into an AntibodyForests object. Make sure that the germline node is called "germline" and that every line represents a new tree in the newick file. All trees in the same file should be from the same sample.
#' @param file.list list - list of newick files to be converted to AntibodyForests object. Could be a named list where the names correspond to sample IDs.
#' @param file.dir directory - directory where the newick files are stored. If provided, the function will read all newick files in the directory.
#' @return AntibodyForests object
#' @export
#' @examples
#' \dontrun{
#' af <- newick_to_Af(file.list = list("S1" = "path/to/sample1.nwk", "S2" = "path/to/sample2.nwk"))
#' }


newick_to_Af <- function(file.list,
                         file.dir){
  # Check the input
  if(missing(file.list) && !missing(file.dir)){stop("Please provide either a directory or a list of newick files to be converted to AntibodyForests object.")}
  if(missing(file.list)){file.list <- list.files(file.dir, pattern = "\\.newick|\\.nwk$", full.names = TRUE)}
  if(length(file.list) == 0){stop("No newick files found in the provided directory. Please provide a valid directory with newick files.")}
  
  reorder_edges <- function(edges){
    
    # Reorganizes the edges in a dataframe to enable the construction of a directed graph, ensuring that the 'germline' node is placed in the first row and in the 'upper.node' column, and all its descendants are in subsequent rows.
    # Arguments:
    # - edges: dataframe that contains two columns ('upper.node' and 'lower.node') that contain names of the nodes of the network/tree, whereby each row represent an edge
    
    # Retrieve the names of all the nodes in the input dataframe
    nodes <- unique(c(edges$upper.node, edges$lower.node))
    
    # Create new dataframe 'edges_reorganized' to store reorganized edges and select the edges containing the germline node
    edges_reorganized <- edges[edges$upper.node == "germline" | edges$lower.node == "germline", ]
    
    # Make sure that the germline node is in the 'upper.column' and swap the nodes if necessary
    edges_reorganized[edges_reorganized$lower.node == "germline", ] <- edges_reorganized[edges_reorganized$lower.node == "germline", c("lower.node", "upper.node", colnames(edges_reorganized[3:ncol(edges_reorganized)]))]
    
    # Start reordering the nodes in the 'edges' dataframe from the 'germline'
    current_upper_nodes <- edges_reorganized[edges_reorganized$upper.node == "germline", "lower.node"]
    processed_nodes <- c("germline", current_upper_nodes)
    
    # Keep reordering the nodes in the 'edges' dataframe until all nodes are processed and present in the 'processed_nodes' vector
    while(all(nodes %in% processed_nodes) == FALSE){
      
      # Create empty vector to store nodes that will be connected to the nodes in the 'current_upper_nodes' vector
      processed_lower_nodes <- c()
      
      # Iterate through nodes in the 'current_upper_nodes' vector
      for(upper_node in current_upper_nodes){
        
        # Select the rows/edges from 'edges' dataframe that contain the current 'upper_node'
        selected_edges <- edges[edges$upper.node == upper_node | edges$lower.node == upper_node, ]
        
        # Retrieve all the nodes that are present in this selection of edges
        selected_nodes <- unique(c(selected_edges$upper.node, selected_edges$lower.node))
        
        # Remove the nodes that are already processed, the remaining nodes will be the lower nodes of the current 'upper_node'
        current_lower_nodes <- selected_nodes[!selected_nodes %in% processed_nodes]
        
        # Iterate through nodes in 'current_lower_nodes'
        for(lower_node in current_lower_nodes){
          
          # Append edge to the 'edges_reorganized' dataframe
          edges_reorganized <- rbind(edges_reorganized, c(upper_node, lower_node, edges[(edges$upper.node == upper_node & edges$lower.node == lower_node) | (edges$upper.node == lower_node & edges$lower.node == upper_node), colnames(edges_reorganized[3:ncol(edges_reorganized)])]))
        }
        
        # After iterating trough nodes in 'current_lower_nodes', ...
        processed_lower_nodes <- c(processed_lower_nodes, current_lower_nodes)
      }
      
      # All the nodes in 'processed_lower_nodes' will be the 'current_upper_nodes' in the next iteration
      current_upper_nodes <- processed_lower_nodes
      
      # Update 'processed_nodes' vector by appending the nodes in the 'processed_lower_nodes' vector to the 'processed_nodes' vector
      processed_nodes <- c(processed_nodes, processed_lower_nodes)
    }
    
    # Return reorganized dataframe
    return(edges_reorganized)
  }
  
  # Construct the list of AntibodyForests objects from the newick files
  output_list <- lapply(file.list, function(file) {
    # Read the raw file and remove RTF formatting
    lines <- readLines(file, warn = F)
    rtf_text <- paste(lines, collapse = " ")
    plain_text <- gsub("\\\\[a-zA-Z0-9]+[ ]?|\\{[^}]*\\}|\\}", "", rtf_text)
    plain_text <- gsub("\\\\'", "", plain_text)
    plain_text <- trimws(plain_text)
    raw_file <- unlist(strsplit(plain_text, "\\ "))
    raw_file <- gsub("\\\\", "", raw_file)
    
    sample_list <- list()
    for (line in seq(1:length(raw_file))) {
      # Read the newick file
      tree_data <- ape::read.tree(text = raw_file[line])
      
      if("germline" %in% tree_data$tip.label == FALSE){
        stop("The newick file does not contain a 'germline' node, please make sure that the germline node is called 'germline'.")
      }
      
      # Create edge dataframe
      edges <- data.frame(upper.node = tree_data$edge[, 1],
                          lower.node = tree_data$edge[, 2],
                          edge.length = tree_data$edge.length)
      
      # Rename the nodes
      n_nodes <- length(tree_data$tip.label) + tree_data$Nnode
      nodes <- seq(1:n_nodes)
      nodes <- nodes
      names(nodes) <- seq(1:n_nodes)
      names(nodes)[1:length(tree_data$tip.label)] <- tree_data$tip.label
      names(nodes) <- paste0("node",names(nodes))
      names(nodes)[names(nodes) == "nodegermline"] <- "germline"
      edges$upper.node <- names(nodes)[match(edges$upper.node, nodes)]
      edges$lower.node <- names(nodes)[match(edges$lower.node, nodes)]
      
      # Create igraph object
      edges <- reorder_edges(edges) #germline as first row for correct directed graph construction
      igraph_object <- igraph::graph_from_data_frame(edges, directed = T)
      
      # Create nodes list
      nodes_list <- list()
      for (node in names(nodes)) {
        if(node != "germline" && node %in% paste0("node",tree_data$tip.label)){nodes_list[[node]] <- list("barcodes" = node, "size" = 1)}
      }
      
      # Add to the sample
      sample_list[[paste0("clonotype", line)]] <- list("nodes" = nodes_list,
                                                       "phylo" = tree_data,
                                                       "igraph" = NULL,
                                                       "igraph.with.inner.nodes" = igraph_object,
                                                       "edges" = NULL,
                                                       "edges.with.inner.nodes" = edges)
    }
    return(sample_list)
    
  })
  
  # Convert 'reorganized_output_list' of class 'list' into object of class 'AntibodyForests'
  AntibodyForests_object <- base::structure(output_list, class = "AntibodyForests")
  
  return(AntibodyForests_object)
  
}
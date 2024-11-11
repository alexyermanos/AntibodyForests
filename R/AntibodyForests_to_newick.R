#'Saves an AntibodyForests-object into a newick file
#'@description Saves an AntibodyForests-object into a newick file. The node labels will have the format node\@size where size is the size of the node.
#' @param AntibodyForests_object AntibodyForests-object, output from AntibodyForests()
#' @param min.nodes The minimum number of nodes in a tree to calculate metrics (including the germline).
#' @param output.file string - specifies the path to the output file
#' @export
#' @examples

AntibodyForests_to_newick <- function(AntibodyForests_object,
                                      min.nodes,
                                      output.file){
  
  if(missing(AntibodyForests_object)){stop("Please provide an AntibodyForests object")}
  if(missing(output.file)){stop("Please provide an output file")}
  if(missing(min.nodes)){min.nodes <- 2}
  
  #Delete file if it exists
  if (file.exists(output.file)){file.remove(output.file)}
  
  #loop over each tree
  for (sample in names(AntibodyForests_object)){
    for (tree in names(AntibodyForests_object[[sample]])){
      if(igraph::vcount(AntibodyForests_object[[sample]][[tree]]$igraph) >= min.nodes){
        #Transform igraph into phylo object
        phylo <- AntibodyForests_phylo(AntibodyForests_object[[sample]][[tree]]$igraph, solve_multichotomies = F)
        
        #Add the node sizes to the tip and node labels
        size_list <- lapply(AntibodyForests_object[[sample]][[tree]]$nodes, function(x){if (is.null(x$size)){return(1)}else{return(x$size)}})
        phylo$tip.label <- paste0(sample, "_", tree, "_", phylo$tip.label, "@", size_list[phylo$tip.label])
        phylo$node.label <- paste0(sample, "_", tree, "_", phylo$node.label, "@", size_list[phylo$node.label])
        
        #Write to newick format
        phylo$edge.length <- as.numeric(phylo$edge.length)
        ape::write.tree(phylo, file = output.file,
                        tree.names = paste0(sample, "_", tree, ":"), append = T)
      }
    }
  }
}
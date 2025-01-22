#' Function to get the sequences from the nodes in an AntibodyForest object
#' @description Function to get the sequences from the nodes in an AntibodyForest object
#' @param AntibodyForests_object AntibodyForests-object, output from Af_build()
#' @param sequence.name character, name of the sequence column in the AntibodyForests object (example VDJ_sequence_aa_trimmed)
#' @param min.nodes integer, minimum number of nodes in the tree (not including germline)
#' @param min.edges integer, minimum number of edges in the tree (not including edges to the germline)
#' @return A dataframe with the sequences and sequence identifiers
#' @export
#' @examples
#' sequence_df <- Af_get_sequences(AntibodyForests::small_af,
#'                sequence.name = "VDJ_sequence_aa_trimmed")

Af_get_sequences <- function(AntibodyForests_object,
                                          sequence.name,
                                          min.nodes,
                                          min.edges){

  if(missing(AntibodyForests_object)){stop("Please provide an AntibodyForests object")}
  if(missing(sequence.name)){stop("Please provide the name of the sequences in the AntibodyForests object")}
  if(missing(min.nodes)){min.nodes <- 1}
  if(missing(min.edges)){min.edges <- 0}

  message("This function can have a long runtime for large AntibodyForests-objects.")


  #Function to get the sequences from a tree
  df_per_clone <- function(AntibodyForests_object, sample, clonotype){
    #Igraph object
    tree <- AntibodyForests_object[[sample]][[clonotype]][["igraph"]]

    #Check the number of edges
    #Get edgelist
    edges <- igraph::as_edgelist(tree, names = T)
    edges <- as.data.frame(edges)
    #Remove germline from the edge list
    edges <- edges[edges$V1 != "germline" & edges$V1 != "germline",]
    nr_edges = nrow(edges)

    #Check the number of nodes
    nr_nodes <- length(AntibodyForests_object[[sample]][[clonotype]][["nodes"]]) - 1

    #If there are not enough edges or nodes, return NA
    if (nr_edges >= min.edges & nr_nodes >= min.nodes){
      seqs <- c()
      seq_ids <- c()
      #Get the sequences
      for (node in names(AntibodyForests_object[[sample]][[clonotype]][["nodes"]])){
        if (node != "germline"){
          seq = AntibodyForests_object[[sample]][[clonotype]][["nodes"]][[node]][[sequence.name]]
          seqs <- c(seqs, seq)
          seq_ids <- c(seq_ids, paste0(sample, "_", clonotype, "_", node))
        }
      }
      df <- data.frame(sequence = seqs, sequence_id = seq_ids)
    }else{
      df <- data.frame(sequence = NA, sequence_id = NA)
    }
    return(df)
  }

  #Initiate the output dataframe
  output_df <- data.frame()
  for (sample in names(AntibodyForests_object)){
    for (clonotype in names(AntibodyForests_object[[sample]])){
      #Get the sequences per tree
      df <- df_per_clone(AntibodyForests_object, sample, clonotype)

      #Add to the output dataframe
      output_df <- rbind(output_df, df)
    }
  }
  #Remove trees that did not pass the requirements
  output_df <- stats::na.omit(output_df)

  return(output_df)

}

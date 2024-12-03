#' Function to create a dataframe of the Protein Language Model probabilities and ranks of the mutations along the edges of B cell lineage trees.
#' @description Function to create a dataframe of the Protein Language Model probabilities and ranks of the mutations along the edges of B cell lineage trees.
#' @param AntibodyForests_object AntibodyForests-object, output from Af_build()
#' @param sequence.name character, name of the sequence column in the AntibodyForests object (example VDJ_sequence_aa_trimmed)
#' @param path_to_probabilities character, path to the folder containing probability matrices for all sequences. Probability matrices should be in CSV format and the filename should include sampleID_clonotypeID_nodeNR, matching the AntibodyForests-object.
#' @return a dataframe with the sample, clonotype, node numbers, number of substitutions, mean substitution rank and mean substitution probability

Af_PLM_dataframe <- function(AntibodyForests_object,
                                          sequence.name,
                                          path_to_probabilities){
  
  if(missing(AntibodyForests_object)){stop("Please provide an AntibodyForests object")}
  if(missing(sequence.name)){stop("Please provide the name of the sequences in the AntibodyForests object")}
  if(missing(path_to_probabilities)){stop("Please provide the path to the folder containing the probability matrices.")}
  
  #Check if probability matrices are present
  prob_matrix_list <- list.files(path_to_probabilities)
  if (length(prob_matrix_list) == 0){stop("No probability matrices found in the specified folder.")}
  if (length(grep(".csv", prob_matrix_list)) == 0){stop("Probabilty matrices should be in CSV format.")}
  
  message("This function can have a long runtime for large AntibodyForests-objects.")
  
  #Create dataframe per clonotype
  df_per_clone <- function(sample, clonotype){
    #Igraph object
    tree <- AntibodyForests_object[[sample]][[clonotype]][["igraph"]]
    
    #Get edgelist
    edges <- igraph::as_edgelist(tree, names = T)
    edges <- as.data.frame(edges)
    
    #Remove germline from the edge list
    edges <- edges[edges$V1 != "germline" & edges$V1 != "germline",]
    colnames(edges) <- c("node1", "node2")
    
    #If there are not enough edges, return NA
    if (nrow(edges) > 0){
      #Get the sequences
      nodes <- AntibodyForests_object[[sample]][[clonotype]][["nodes"]]
      
      #Initiate output dataframe
      clonotype_df <- data.frame(sample = character(), clonotype = character(), node1 = character(), node2 = character(),
                       n_subs = numeric(), mean_substitution_rank = numeric(), mean_substitution_probability = numeric())
      
      for (row in 1:nrow(edges)){
        node1 <- edges[row, "node1"]
        node2 <- edges[row, "node2"]
        seq1 <- strsplit(nodes[[node1]][[sequence.name]], "")[[1]]
        seq2 <- strsplit(nodes[[node2]][[sequence.name]], "")[[1]]
        
        #Only for sequences of the same length
        if (length(seq1) == length(seq2)){
          #Get the probability matrix
          prob_file <- prob_matrix_list[grep(paste0(sample,"_",clonotype,"_",node1,"_"), prob_matrix_list)]
          if (length(prob_file) == 1){
            prob_matrix <- read.csv(paste0(path_to_probabilities, "/", prob_file), header = T)
            
            #Get the mutating positions
            diff_positions <- c()
            for (k in 1:length(seq1)){if (seq1[k] != seq2[k]){diff_positions <- c(diff_positions, k)}}
            
            #Initiate rank vector and probabilities vector
            substitute_ranks <- c()
            substitute_probabilities <- c()
            original_ranks <- c()
            original_probabilities <- c()
            
            #Loop over the mutational positions
            for (pos in diff_positions){
              #Get the aa probabilities for this position
              likelihood_values <- prob_matrix[pos,]
              
              #Get the rank for each aa
              ranks <- rank(-likelihood_values)
              
              #Get the rank of the mutations
              mut_rank <- ranks[seq2[pos]]
              substitute_ranks <- c(substitute_ranks, mut_rank)
              
              #Get the probability of the mutation
              probability <- likelihood_values[seq2[pos]]
              substitute_probabilities <- c(substitute_probabilities, probability)
              
              #Get the rank of the original residue
              orig_rank <- ranks[seq1[pos]]
              original_ranks <- c(original_ranks, orig_rank)
              
              #Get the probability of the original
              orig_probability <- likelihood_values[seq1[pos]]
              original_probabilities <- c(original_probabilities, orig_probability)
            }
            #If there are multiple mutations in a sequence, get the average
            mean_substitution_rank <- mean(substitute_ranks)
            mean_original_rank <- mean(original_ranks)
            mean_substitution_probability <- mean(unlist(substitute_probabilities))
            mean_original_probability <- mean(unlist(original_probabilities))
            
            #Add to the output dataframe
            edge_df <- data.frame(sample = sample, clonotype = clonotype, n_subs = length(diff_positions), node1 = node1, node2 = node2,
                                  mean_substitution_rank = mean_substitution_rank, mean_substitution_probability = mean_substitution_probability,
                                  mean_original_rank = mean_original_rank, mean_original_probability = mean_original_probability)
            clonotype_df <- rbind(clonotype_df, edge_df)
          }
          else{
            print("No probability matrix found for ", paste0(sample,"_",clonotype,"_",node1))
            clonotype_df <- rbind(clonotype_df, data.frame(sample = NA, clonotype = NA, node1 = NA, node2 = NA,
                                 n_subs = NA, mean_substitution_rank = NA, mean_substitution_probability = NA,
                                 mean_original_rank = NA, mean_original_probability = NA))
          }
        }
      }
    }
    else{
      clonotype_df <- data.frame(sample = NA, clonotype = NA, node1 = NA, node2 = NA,
                                 n_subs = NA, mean_substitution_rank = NA, mean_substitution_probability = NA,
                                 mean_original_rank = NA, mean_original_probability = NA)
    }
    return(clonotype_df)
  }
  
  output_df <- data.frame(sample = character(), clonotype = character(), node1 = character(), node2 = character(),
                          n_subs = numeric(), mean_substitution_rank = numeric(), mean_substitution_probability = numeric(),
                          mean_original_rank = numeric(), mean_original_probability = numeric())
  for (sample in names(AntibodyForests_object)){
    for (clonotype in names(AntibodyForests_object[[sample]])){
      tree_df <- df_per_clone(sample, clonotype)
      output_df <- rbind(output_df, tree_df)
    }
  }
  
  return(output_df)
  
}



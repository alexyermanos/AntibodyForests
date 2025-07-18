#' Function to compare the distributions of the Protein Language Model probabilities or ranks of the mutations along the edges of B cell lineage trees across repertoires using the Jensen-Shannon divergence.
#' @description Function to compare the distributions of the Protein Language Model probabilities or ranks of the mutations along the edges of B cell lineage trees across repertoires using the Jensen-Shannon divergence.
#' @param PLM_dataframe Dataframe resulting from Af_PLM_dataframe(). This contains the Protein Language Model probabilities and ranks of the mutations along the edges of B cell lineage trees.
#' @param values What values to compare
#' "substitution_rank" will compare the rank of the mutation along the edge of the tree (Highest probability is rank 1).
#' "substitution_probability" will copmare the probability of the mutation along the edge of the tree.
#' "original_rank" will compare the rank of the original amino acid at the site of mutation along the edge of the tree (Highest probability is rank 1).
#' "original_probability" will compare the probability of the original amino acid at the site of mutation along the edge of the tree.
#' "unmutated_rank" will compare the rank of the unmutated amino acids along the edge of the tree (Highest probability is rank 1).
#' "unmutated_probability" will compare the probabilities of the unmutated amino acids along the edge of the tree.
#' @param font.size Font size for the plot. Default is 16.
#' @param output.file string - specifies the path to the output file (PNG of PDF). Defaults to NULL.
#' @return A pheatmap of the Jensen-Shannon distance between repertoires
#' @export
#' @importFrom dplyr .data
#' @examples
#' Af_compare_PLM(PLM_dataframe = AntibodyForests::PLM_dataframe,
#'             values = "original_probability")
          
Af_compare_PLM <- function(PLM_dataframe,
                           values,
                           font.size,
                           output.file){
  
  #Check input
  if(missing(PLM_dataframe)){stop("Please provide a PLM dataframe resulting from Af_PLM_dataframe function.")}
  if(length(unique(PLM_dataframe$sample)) < 3){stop("Please provide a PLM dataframe with at least three samples.")}
  if(all(colnames(PLM_dataframe) %in% c("sample", "clonotype", "n_subs", "node1", "node2", "mean_substitution_rank", "
                                  mean_substitution_probability"))){stop("Please provide a PLM dataframe resulting from Af_PLM_dataframe function.")}
  
  #Set defaults
  if(missing(values)){values <- "substitution_rank"}
  if(missing(font.size)){font.size <- 16}
  if(missing(output.file)){output.file <- NULL}
  
  if(values == "substitution_rank"){plot_values <- "mean_substitution_rank"}
  if(values == "substitution_probability"){plot_values <- "mean_substitution_probability"}
  if(values == "original_rank"){plot_values <- "mean_original_rank"}
  if(values == "original_probability"){plot_values <- "mean_original_probability"}
  if(values == "unmutated_rank"){plot_values <- "mean_unmutating_rank"}
  if(values == "unmutated_probability"){plot_values <- "mean_unmutating_probability"}
  
  #Set global variables for CRAN check
  png <- NULL
  pdf <- NULL
  
  PLM_dataframe <- stats::na.omit(PLM_dataframe)
  
  # Set the bins
  if (length(grep("rank", values))!=0){bins <- seq(0, 22, by = 1)}
  if (length(grep("probability", values))!=0){bins <- seq(0, 1, by = 0.1)}
  
  # Function to create a probability vector from the values
  make_prob_vector <- function(values, bins) {
    hist_vals <- graphics::hist(unlist(values), breaks = bins, plot = FALSE)$counts
    prob <- hist_vals / sum(hist_vals)
    return(prob)
  }
  
  # Create a list of samples and their corresponding values
  samples = list()
  for (sample in unique(PLM_dataframe$sample)){
    sample_values <- as.vector(PLM_dataframe[PLM_dataframe$sample == sample, plot_values])
    if (length(sample_values) > 0){
      samples[[sample]] <- sample_values
    }
  }
  
  # Create a probability matrix from the samples
  prob_matrix <- do.call(rbind, lapply(samples, make_prob_vector, bins = bins))
  rownames(prob_matrix) <- names(samples)
  
  # compute the JS-distance matrix of a given probability matrix
  JSDMatrix <- philentropy::JSD(prob_matrix)
  rownames(JSDMatrix) <- rownames(prob_matrix)
  colnames(JSDMatrix) <- rownames(prob_matrix)
  

  #Create the plot
  p <- pheatmap::pheatmap(JSDMatrix, 
                          fontsize = font.size,
                          main = paste0("JS Distance between PLM ", gsub("_", " ", values), " distribution"))
  
  if(!is.null(output.file)){
    # Check if the output.file is png or pdf
    if (grepl(pattern = ".png$", output.file)){
      png(file = output.file)
      print(p)
      grDevices::dev.off()
    }else if (grepl(pattern = ".pdf$", output.file)){
      pdf(file = output.file)
      print(p)
      grDevices::dev.off()
    }
  }
  
  return(p)
}

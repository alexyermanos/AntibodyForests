#' Function to create a distribution plot of the Protein Language Model probabilities and ranks of the mutations along the edges of B cell lineage trees.
#' @description Function to create a distribution plot of the Protein Language Model probabilities and ranks of the mutations along the edges of B cell lineage trees.
#' @param PLM_dataframe Dataframe resulting from Af_PLM_dataframe(). This contains the Protein Language Model probabilities and ranks of the mutations along the edges of B cell lineage trees.
#' @param values What values to plot. Can be "rank" (default) or "probability".
#' "substitution_rank" will plot the rank of the mutation along the edge of the tree (Highest probability is rank 1).
#' "substitution_probability" will plot the probability of the mutation along the edge of the tree.
#' "original_rank" will plot the rank of the original amino acid at the site of mutation along the edge of the tree (Highest probability is rank 1).
#' "original_probability" will plot the probability of the original amino acid at the site of mutation along the edge of the tree.
#' @param group_by Plot a seperate line per sample or everything together (default).
#' "sample_id"
#' "none"
#' @param colors Color to use for the lines. When group_by = "sample_id": This should be a vector of the same length as the number of samples.
#' @param font.size Font size for the plot. Default is 16.
#' @param output.file string - specifies the path to the output file (PNG of PDF). Defaults to NULL.
#' @export
#' @importFrom dplyr .data
#' @examples
#' \dontrun{
#' Af_plot_PLM(PLM_dataframe = PLM_dataframe,
#'             values = "original_probability",
#'             group_by = "sample_id")
#' }

Af_plot_PLM <- function(PLM_dataframe,
                                     values,
                                     group_by,
                                     colors,
                                     font.size,
                                     output.file){

  #Check input
  if(missing(PLM_dataframe)){stop("Please provide a PLM dataframe resulting from Af_PLM_dataframe function.")}
  if(all(colnames(PLM_dataframe) %in% c("sample", "clonotype", "n_subs", "node1", "node2", "mean_substitution_rank", "
                                  mean_substitution_probability"))){stop("Please provide a PLM dataframe resulting from Af_PLM_dataframe function.")}

  #Set defaults
  if(missing(values)){values <- "substitution_rank"}
  if(missing(group_by)){group_by <- "none"}
  if(missing(colors)){colors <- NULL}
  if(missing(font.size)){font.size <- 16}
  if(missing(output.file)){output.file <- NULL}

  if(values == "substitution_rank"){plot_values <- "mean_substitution_rank"}
  if(values == "substitution_probability"){plot_values <- "mean_substitution_probability"}
  if(values == "original_rank"){plot_values <- "mean_original_rank"}
  if(values == "original_probability"){plot_values <- "mean_original_probability"}

  #Set global variables for CRAN check
  png <- NULL
  pdf <- NULL


  #Create the plot
  p <- ggplot2::ggplot(PLM_dataframe, ggplot2::aes(.data[[plot_values]])) +
    ggplot2::theme_minimal() +
    ggplot2::theme(text = ggplot2::element_text(size = font.size))

  if (values == "substitution_rank"){
    #Set binwidt to 1 for freqpoly
    bin_width <- 1
    #Set the x-axis
    p <- p + ggplot2::scale_x_continuous(breaks = c(1,5,10,15,20), limits = c(1,NA)) +
      ggplot2::xlab("Substitution Rank") + ggplot2::ylab("Number of edges")
  }
  if (values == "substitution_probability"){
    #Set binwidt to 1 for freqpoly
    bin_width <- 0.1
    #Set the x-axis
    p <- p + ggplot2::scale_x_continuous(limits = c(0,1)) +
      ggplot2::xlab("Substitution Likelihood") + ggplot2::ylab("Number of edges")
  }
  if (values == "original_rank"){
    #Set binwidt to 1 for freqpoly
    bin_width <- 1
    #Set the x-axis
    p <- p + ggplot2::scale_x_continuous(breaks = c(1,5,10,15,20), limits = c(1,NA)) +
      ggplot2::xlab("Original Rank") + ggplot2::ylab("Number of edges")
  }
  if (values == "original_probability"){
    #Set binwidt to 1 for freqpoly
    bin_width <- 0.1
    #Set the x-axis
    p <- p + ggplot2::scale_x_continuous(limits = c(0,1)) +
      ggplot2::xlab("Original Likelihood") + ggplot2::ylab("Number of edges")
  }

  #Plot the lines
  if (group_by == "none"){
    if(is.null(colors)){colors <- "black"}
    p <- p +  ggplot2::geom_freqpoly(binwidth = bin_width, linewidth = 1, color = colors)}
  if (group_by == "sample_id"){
    p <- p +  ggplot2::geom_freqpoly(ggplot2::aes(colour = sample), binwidth = bin_width, linewidth = 1)
    if (!is.null(colors)){p <- p + ggplot2::scale_color_manual(values = colors)}
    }

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
  }else{print(p)}
}

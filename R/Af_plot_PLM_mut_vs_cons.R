#' Function to create a boxplot of the Protein Language Model probabilities
#' @description Function to create a boxplot of the Protein Language Model probabilities and ranks of the mutating vs. conserved residues along the edges of B cell lineage trees.
#' @param PLM_dataframe Dataframe resulting from Af_PLM_dataframe(). This contains the Protein Language Model probabilities and ranks of the mutations along the edges of B cell lineage trees.
#' @param values What values to plot. Can be "rank" (default) or "probability".
#' "rank" will plot the rank of the amino acid (Highest probability is rank 1).
#' "probability" will plot the probability of the amino acid.
#' @param dots Whether to plot the individual points. Can be "none" (default), "all_edges", "sample_average"
#' @param group_by Color the dots on a group. Can be "none" (default), "sample_id", or "n_subs".
#' @param colors Color to use for the dots. When group_by = "sample_id": This should be a vector of the same length as the number of samples.
#' @param font.size Font size for the plot. Default is 16.
#' @param output.file string - specifies the path to the output file (PNG of PDF). Defaults to NULL.
#' @return A ggplot2 object of the PLM boxplot
#' @export
#' @importFrom dplyr .data
#' @examples
#' Af_plot_PLM_mut_vs_cons(PLM_dataframe = AntibodyForests::PLM_dataframe,
#'             values = "probability")


Af_plot_PLM_mut_vs_cons <- function(PLM_dataframe,
                        values,
                        dots,
                        group_by,
                        colors,
                        font.size,
                        output.file){
  
  #Check input
  if(missing(PLM_dataframe)){stop("Please provide a PLM dataframe resulting from Af_PLM_dataframe function.")}
  if(all(colnames(PLM_dataframe) %in% c("sample", "clonotype", "n_subs", "node1", "node2", "mean_original_rank", "
                                  mean_original_probability", "mean_unmutating_rank"))){
    stop("Please provide a PLM dataframe resulting from Af_PLM_dataframe function.")}
  if(!any(values %in% c("rank", "probability"))){stop("Please provide a valid value for 'values'. Can be 'rank' or 'probability'.")}
  if(missing(dots)){dots <- "none"}
  if(!any(dots %in% c("none", "sample_average", "clonotype_average", "all_edges"))){
    stop("Please provide a valid value for 'dots'. Can be 'none', 'sample_average', or 'all_edges'.")}
  
  #Set defaults
  if(missing(values)){values <- "rank"}
  if(missing(group_by)){group_by <- "none"}
  if(missing(colors)){colors <- NULL}
  if(missing(font.size)){font.size <- 16}
  if(missing(output.file)){output.file <- NULL}
  
  if(values == "rank"){
    PLM_dataframe <- tidyr::pivot_longer(PLM_dataframe, names_to = "residue", values_to = "rank",
                        cols = c("mean_original_rank", "mean_unmutating_rank"))}
  if(values == "probability"){
    PLM_dataframe <- tidyr::pivot_longer(PLM_dataframe, names_to = "residue", values_to = "probability",
                                                                   cols = c("mean_original_probability", "mean_unmutating_probability"))}
  #Change names
  PLM_dataframe <- dplyr::mutate(PLM_dataframe, 
                             residue = dplyr::case_when(
                               residue == "mean_original_rank" ~ "Mutating",
                               residue == "mean_unmutating_rank" ~ "Conserved",
                               residue == "mean_original_probability" ~ "Mutating",
                               residue == "mean_unmutating_probability" ~ "Conserved"))
  
  #Set global variables for CRAN check
  png <- NULL
  pdf <- NULL
  residue <- NULL
  
  PLM_dataframe <- stats::na.omit(PLM_dataframe)
  
  if (group_by == "n_subs"){
    PLM_dataframe$n_subs <- dplyr::case_match(PLM_dataframe$n_subs,
                                              1 ~ "1",
                                              seq(2,max(PLM_dataframe$n_subs)) ~ ">1")
  }
  if (group_by == "sample_id"){group_by <- "sample"}
  
  
  #Create the plot
  p <- ggplot2::ggplot(PLM_dataframe, ggplot2::aes(x = residue, y = .data[[values]])) +
    ggplot2::theme_minimal() +
    ggplot2::theme(text = ggplot2::element_text(size = font.size),
                   axis.title.y = ggplot2::element_blank())
  #Set the y-axis
  if (values == "rank"){p <- p + ggplot2::ylab("Average Likelihood Rank")}
  if (values == "probability"){p <- p + ggplot2::ylab("Average Residue Likelihood")}
  
  #Plot the individual points (mean per sample)
  if (dots != "none"){
    if (dots == "all_edges"){
      if (group_by == "none"){
        if(is.null(colors)){colors <- "black"}
        p <- p +  ggplot2::geom_jitter(color = colors, size = 0.1)}
      else{
        p <- p +  ggplot2::geom_jitter(ggplot2::aes(colour = .data[[group_by]]), size = 0.1)
        if (!is.null(colors)){p <- p + ggplot2::scale_color_manual(values = colors)}
      }
    }
    if (dots == "sample_average"){
      p <- p + ggplot2::stat_summary(ggplot2::aes(colour = sample), 
                                     fun = "mean", geom = "point", size = 1.5)
      if (!is.null(colors)){p <- p + ggplot2::scale_color_manual(values = colors)}
    }
  }
  
  #Plot the boxes
  p <- p +  ggplot2::geom_boxplot(color = "black", fill = NA)
  
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

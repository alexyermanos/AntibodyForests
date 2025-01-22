#' Function to calculate the RMSD between sequences over each edge in the AntibodyForest object
#' @description This function calculates the RMSD between sequences over each edge in the AntibodyForest object.
#' @param AntibodyForests_object AntibodyForests-object, output from Af_build()
#' @param VDJ The dataframe with V(D)J information such as the output of Platypus::VDJ_build() that was used to create the AntibodyForests-object. Must contain columns sample_id, clonotype_id, barcode.
#' @param pdb.dir a directory containing PDB files.
#' @param file.df a dataframe of pdb filenames (column file_name) to be used and sequence IDs (column sequence) corresponding to the the barcodes in the AntibodyForests-object
#' @param sequence.region a character vector of the sequence region to be used to calculate properties. Default is "full.sequence".
#' - full.sequence: the full sequence(s) in the PDB file
#' - sub.sequence: part of the full sequence, for example the CDR3 region in the PDB file. This sub sequence must be a column in the VDJ dataframe.
#' - binding.residues: the binding residues in the PDB file
#' @param sub.sequence.column a character vector of the column name in the VDJ dataframe containing the sub sequence to be used to calculate properties. Default is NULL.
#' @param chain a character vector of the chain to be used to calculate properties. Default is both heavy and light chain
#' Assuming chain "A" is heavy chain, chain "B" is light chain, and possible chain "C" is the antigen.
#' - HC+LC: both heavy and light chain
#' - HC: heavy chain, assuming chain A is the heavy chain.
#' - LC: light chain, assuming chain B is the light chain.
#' - AG: antigen, assuming chain C is the antigen.
#' - whole.complex: the whole complex of antibody-antigen (all available chains in the pdb file).
#' @param font.size The font size of the plot. Default is 12.
#' @param point.size The size of the points in the scatterplot. Default is 1.
#' @param color The color of the dots in the scatterplot. Default is "black".
#' @param output.file string - specifies the path to the output file (PNG of PDF). Defaults to NULL.
#' @return A list with the edge dataframe and a ggplot2 object
#' @export
#' @examples
#' \dontrun{
#' rmsd_df <- Af_edge_RMSD(AntibodyForests::small_af,
#'                        VDJ = AntibodyForests::small_vdj,
#'                        pdb.dir = "~/path/PDBS_superimposed/",
#'                        file.df = files,
#'                        sequence.region = "full.sequence",
#'                        chain = "HC+LC")}


Af_edge_RMSD <- function(AntibodyForests_object,
                         VDJ,
                         pdb.dir,
                         file.df,
                         sequence.region,
                         sub.sequence.column,
                         chain,
                         font.size,
                         point.size,
                         color,
                         output.file){

  if(missing(AntibodyForests_object)){stop("Please provide an AntibodyForests object")}
  if(missing(VDJ)){stop("Please provide the VDJ dataframe")}
  if(!any(c("sample_id", "clonotype_id", "barcode") %in% colnames(VDJ))){stop("VDJ dataframe must contain columns sample_id, clonotype_id, barcode")}
  if(missing(pdb.dir)){stop("pdb.dir is missing")}
  if(!dir.exists(pdb.dir)){stop("pdb.dir does not exist")}
  if(missing(file.df)){stop("file.df directory is missing")}
  if(!any(c("file_name", "sequence") %in% colnames(file.df))){stop("file.df dataframe must contain columns file_name, sequence")}
  if(!(all(file.df$file_name %in% list.files(pdb.dir)))){stop("file.df contains files that do not exist in pdb.dir")}
  if(missing(sequence.region)){sequence.region = "full.sequence"}
  if(!sequence.region %in% c("full.sequence", "sub.sequence", "binding.residues")){stop("sequence.region must be one of 'full.sequence', 'sub.sequence', 'binding.residues'")}
  if(missing(chain)){chain = "HC+LC"}
  if(!chain %in% c("HC+LC", "HC", "LC", "AG", "whole.complex")){stop("chain must be one of 'HC+LC', 'HC', 'LC', 'AG', 'whole.complex'")}
  if(missing(sub.sequence.column)){sub.sequence.column = NULL}
  if(!is.null(sub.sequence.column)){if(!sub.sequence.column %in% colnames(VDJ)){stop("sub.sequence.column does not exist in VDJ dataframe")}}
  if(sequence.region == "sub.sequence" & is.null(sub.sequence.column)){stop("sub.sequence.column is missing")}
  if(missing(font.size)){font.size = 12}
  if(missing(point.size)){point.size = 1}
  if(missing(output.file)){output.file = NULL}
  if(missing(color)){color = "black"}

  #Global variable definitions for CRAN checks
  resno <- NULL
  resid <- NULL
  str_count <- NULL
  n_subs <- NULL
  edge_RMSD <- NULL
  png <- NULL
  pdf <- NULL

  output_df <- data.frame(sample = character(), clonotype = character(), node1 = character(), node2 = character(),
                          n_subs = numeric(), edge_RMSD = numeric())


  #Set the chains
  if(chain == "HC+LC"){chains = c("A", "B")}
  if(chain == "HC"){chains = "A"}
  if(chain == "LC"){chains = "B"}
  if(chain == "AG"){chains = "C"}
  if(chain == "whole.complex"){chains = c("A", "B", "C")}

  #Adapted from Lucas' Structure_utils.py
  #Function to get the binding site residues in a PDB file
  #* pdb object
  #* chains_binder: The chain or a list of chains of the binders
  #* chains_ligand: The chain or a list of chains of the ligands
  #* cutoff: The cutoff in distance (A) between atoms to be part of the binding site.
  get_bs_res <- function(pdb, chains_binder, chains_ligand, cutoff = 5) {
    # Extract atoms for binders and ligands
    atoms_binder <- pdb$atom[pdb$atom$chain %in% chains_binder, ]
    atoms_ligand <- pdb$atom[pdb$atom$chain %in% chains_ligand, ]
    # Extract coordinates
    coord_binder <- atoms_binder[, c("x", "y", "z")]
    coord_ligand <- atoms_ligand[, c("x", "y", "z")]
    # Calculate the distance matrix
    dist_mat <- bio3d::dist.xyz(coord_binder, coord_ligand)
    # Apply cutoff to distance matrix
    dist_mat_bin <- dist_mat <= cutoff
    # Get binding site residues
    bs_res_binder <- unique(atoms_binder[rowSums(dist_mat_bin) >= 1, "resno"])
    bs_res_ligand <- unique(atoms_ligand[colSums(dist_mat_bin) >= 1, "resno"])
    return(c(bs_res_binder, bs_res_ligand))
  }

  #Get the resnos from a subsequences in the PDB file
  get_subseq_res <- function(pdb, sub_seq, CHAIN){
    #Get the residues in the PDB file
    pdb$atom |> dplyr::filter(chain %in% CHAIN) |> dplyr::distinct(resno, .keep_all = TRUE) |> dplyr::pull(resid, resno) -> resids
    resids_df <- as.data.frame(resids)
    resids_df$resno <- rownames(resids_df)
    resids_df$resids <- sapply(resids_df$resids, function(x){x <- stringr::str_to_title(x);return(seqinr::a(x))})
    #Align the subsequence to the PDB sequence
    pwa <- pwalign::pairwiseAlignment(sub_seq, paste(resids_df$resids, collapse = ""))
    aligned_subseq <- as.character(pwalign::alignedPattern(pwa))
    #Get the matching indices
    ind <- which(strsplit(aligned_subseq, "")[[1]] != "-")
    #Get the corresponding resnos
    resnos <- resids_df$resno[ind]

    return(as.numeric(resnos))
  }

  SHM_per_alignment <- function(seq_1, seq_2){
    alignment <- pwalign::pairwiseAlignment(seq_1, seq_2)
    seq_1 <- as.character(alignment@pattern)
    seq_2 <- as.character(alignment@subject)

    hamming_dist <- stringdist::stringdist(seq_1,seq_2, method = "hamming")
    gaps_1 <- str_count(seq_1,"-")
    gaps_2 <- str_count(seq_2,"-")


    return(hamming_dist - gaps_1 - gaps_2)

  }

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
                                 n_subs = numeric(), edge_RMSD = numeric())

      for (row in 1:nrow(edges)){
        node1 <- edges[row, "node1"]
        node2 <- edges[row, "node2"]
        barcode1 <- unique(nodes[[node1]][["barcodes"]])
        barcode2 <- unique(nodes[[node2]][["barcodes"]])
        #If multiple barcodes, take the first one
        if (length(barcode1) > 1){barcode1 <- barcode1[1]}
        if (length(barcode2) > 1){barcode2 <- barcode2[1]}



        #Read in the pdb file
        file1 <- file.df[file.df$sequence == barcode1,"file_name"]
        file2 <- file.df[file.df$sequence == barcode2,"file_name"]
        pdb1 <- bio3d::read.pdb(paste0(pdb.dir, file1))
        pdb2 <- bio3d::read.pdb(paste0(pdb.dir, file2))

        #Get the residues for RMSD calculations
        if(sequence.region == "full.sequence"){
          resnos1 <- unique(pdb1$atom[pdb1$atom$chain %in% chains, ]$resno)
          resnos2 <- unique(pdb2$atom[pdb2$atom$chain %in% chains, ]$resno)
        }
        if(sequence.region == "binding.residues"){
          resnos1 <- get_bs_res(pdb1, chains_binder = c("A", "B"), chains_ligand = "C")
          resnos2 <- get_bs_res(pdb2, chains_binder = c("A", "B"), chains_ligand = "C")
        }
        if(sequence.region == "sub.sequence"){
          #get the subsequence
          sub_seq1 <- VDJ[VDJ$barcode == barcode1,sub.sequence.column]
          sub_seq2 <- VDJ[VDJ$barcode == barcode2,sub.sequence.column]
          resnos1 <- get_subseq_res(pdb1, sub_seq1, chains)
          resnos2 <- get_subseq_res(pdb2, sub_seq2, chains)
        }

        #Get the mutating positions
        total_subs = 0
        for (i in chains){
          seq1 <- paste(pdb1$atom |> dplyr::filter(chain == i) |> dplyr::filter(resno %in% resnos1) |>
                          dplyr::distinct(resno, .keep_all = T) |> dplyr::pull(resid) |> stringr::str_to_title() |> seqinr::a(), collapse = "")
          seq2 <- paste(pdb2$atom |> dplyr::filter(chain == i) |> dplyr::filter(resno %in% resnos2) |>
                          dplyr::distinct(resno, .keep_all = T) |> dplyr::pull(resid) |> stringr::str_to_title() |> seqinr::a(), collapse = "")
          n_subs <- SHM_per_alignment(seq1, seq2)
          total_subs <- total_subs + n_subs
        }


        #Get the RMSD between the two sequences
        ca.inds1 <- bio3d::atom.select(pdb1, "calpha", chain = chains, resno = resnos1)
        ca.inds2 <- bio3d::atom.select(pdb2, "calpha", chain = chains, resno = resnos2)
        rmsd <- bio3d::rmsd(pdb1, pdb2, fit = T, a.inds = ca.inds1, b.inds = ca.inds2)

        #Add to the output dataframe
        edge_df <- data.frame(sample = sample, clonotype = clonotype, n_subs = total_subs, node1 = node1, node2 = node2, edge_RMSD = rmsd)
        clonotype_df <- rbind(clonotype_df, edge_df)
      }
    }
    else{
      clonotype_df <- data.frame(sample = NA, clonotype = NA, node1 = NA, node2 = NA,n_subs = NA, edge_RMSD = NA)
    }
    return(clonotype_df)
  }

  for (sample in names(AntibodyForests_object)){
    for (clonotype in names(AntibodyForests_object[[sample]])){
      tree_df <- df_per_clone(sample, clonotype)
      output_df <- rbind(output_df, tree_df)
    }
  }

  #Calculate the correlation between number of substitutions and edge RMSD
  cor <- stats::cor.test(as.numeric(output_df$n_subs), output_df$edge_RMSD, method = "pearson", exact = F)

  p <- ggplot2::ggplot(output_df, ggplot2::aes(x = as.numeric(n_subs), y = edge_RMSD)) +
    ggplot2::geom_point(size = point.size, color = color) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = font.size)) +
    ggplot2::geom_smooth(method = "lm", color = "black") +
    ggplot2::ylab("Edge RMSD") + ggplot2::xlab("Number of substitutions") +
    ggplot2::ggtitle(paste0("R\u00b2 = ", round(cor$estimate, digits = 2)))

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

  return(list(dataframe = output_df, plot = p))


}



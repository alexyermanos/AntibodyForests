#' A function to integrate bulk and single cell data
#' @description Integrate bulk and single-cell data by reannotating the germline genes and integrating the bulk sequences into the existing single-cell clonotypes.
#' @param sc.VDJ VDJ dataframe of the single cell data created with Platypus VDJ_build function.
#' @param bulk.tsv A tab separated file of the bulk sequences with the at least columns containing the sequence, a sample ID, a barcode, and the isotype.
#' @param bulk.tsv.sequence.column column name of the bulk tsv that contains the nucleotide sequence
#' @param bulk.tsv.sample.column column name of the bulk tsv that contains the sample_id that matches the sample_id in sc_VDJ
#' @param bulk.tsv.barcode.column column name of the bulk tsv that contains the barcode/identifier of the recovered sequence
#' @param bulk.tsv.isotype.column column name of the bulk tsv that contains the isotype of the recovered sequence
#' @param organism "human" or "mouse"
#' @param scRNA_seqs_annotations A tab separated file of the reannotated single-cell sequences using Change-O AssignGenes.py. If NULL, this function will run Change-O AssignGenes.py (Make sure to have this installed, including igblast.dir). Default is NULL.
#' @param bulkRNA_seqs_annotations A tab separated file of the reannotated bulk sequences using Change-O AssignGenes.py. If NULL, this function will run Change-O AssignGenes.py (Make sure to have this installed, including igblast.dir). Default is NULL.
#' @param igblast.dir directory where the igblast executables are located. For example: use the instruction to set up IgPhyML environment in the AntibodyForests vignette ($(conda info --base)/envs/igphyml/share/igblast)
#' @param trim.FR1 - boolean - whether to trim the FR1 region from the sequences and germline, this is recommended to account for variation in primer design during sequencing (Default is TRUE)
#' @param tie.resolvement How to resolve a bulk sequence for which multiple clonotypes match.
#'  "all" - assign the bulk sequence to all matching clonotypes (Default)
#'  "none" - do not assign the bulk sequence to any clonotype
#'  "random" - randomly assign the bulk sequence to one of the matching clonotypes
#' @param seq.identity sequence identity threshold for clonotype assignment (Default: 0.85)
#' @return The VDJ dataframe of both the bulk and single-cell data
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#' VDJ <- VDJ_integrate_bulk(sc_VDJ = VDJ,
#'   bulk_tsv = "bulk_rna.tsv",
#'   bulk_tsv_sequence_column = "sequence",
#'   bulk_tsv_sample_column = "sample_id",
#'   bulk_tsv_barcode_column = "barcode",
#'   bulk_tsv_isotype_column = "isotype",
#'   organism = "human",
#'   igblast_dir = "anaconda3/envs/igphyml/share/igblast",
#'   tie_resolvement = "random",
#'   seq_identity = 0.85)
#'}

VDJ_integrate_bulk <- function(sc.VDJ,
                               bulk.tsv,
                               bulk.tsv.sequence.column,
                               bulk.tsv.sample.column,
                               bulk.tsv.barcode.column,
                               bulk.tsv.isotype.column,
                               organism,
                               scRNA_seqs_annotations,
                               bulkRNA_seqs_annotations,
                               igblast.dir,
                               trim.FR1,
                               tie.resolvement,
                               seq.identity
                               ){

  #Check missing input
  if(missing(sc.VDJ)){stop("sc.VDJ must be specified")}
  if(missing(bulk.tsv)){stop("bulk.tsv must be specified")}
  if(missing(bulk.tsv.sequence.column)){stop("bulk.tsv.sequence.column must be specified")}
  if(missing(bulk.tsv.sample.column)){stop("bulk.tsv.sample.column must be specified")}
  if(missing(bulk.tsv.barcode.column)){stop("bulk.tsv.barcode.column must be specified")}
  if(missing(bulk.tsv.isotype.column)){stop("bulk.tsv.isotype.column must be specified")}
  if(missing(organism)){stop("organism must be specified")}
  if(missing(bulkRNA_seqs_annotations)){bulkRNA_seqs_annotations <- NULL}
  if(missing(scRNA_seqs_annotations)){scRNA_seqs_annotations <- NULL}
  if(missing(igblast.dir)){igblast.dir <- NULL}
  if(missing(tie.resolvement)){tie.resolvement <- "all"}
  if(missing(seq.identity)){seq.identity <- 0.85}
  if(missing(trim.FR1)){trim.FR1 <- TRUE}

  #Check input validity
  if(!(organism %in% c("human", "mouse"))){stop("organism must be either human or mouse")}
  if(!(tie.resolvement %in% c("all", "none", "random"))){stop("tie.resolvement must be either all, none, or random")}
  if(!(is.numeric(seq.identity) & seq.identity <= 1 & seq.identity >= 0)){stop("seq.identity must be a numeric value between 0 and 1.")}
  if(is.null(bulkRNA_seqs_annotations) & is.null(bulkRNA_seqs_annotations) & is.null(igblast.dir)){stop("igblast.dir must be specified")}

  message("This function can have a long runtime for large datasets.")

  switch(Sys.info()[['sysname']],
         Windows= {message("Windows system detected")
           operating.system <- "Windows"},
         Linux  = {message("Linux system detected")
           operating.system <- "Linux"},
         Darwin = {message("MAC system detected")
           operating.system <- "Darwin"})
  if(operating.system == 'Windows'){
    os_prefix <- 'cmd.exe '
  }else{
    os_prefix <- ''
  }

  #Set global variables for CRAN
  clonotype_id <- NULL

  #Function to build a VDJ dataframe from the IgBLAST annotations
  Transform_to_VDJ <- function(annotations_combined, sc.VDJ, bulk.tsv, bulk.tsv.sample.column, bulk.tsv.barcode.column,
                               bulk.tsv.sequence.column, bulk.tsv.isotype.column, trim.FR1){
    # Process each sample separately and build a combined VDJ dataframe
    VDJ_combined_unclonotyped_unfiltered <- do.call(rbind, lapply(unique(sc.VDJ$sample_id), function(sample){
      # Select relevant columns from single-cell and bulk VDJ dataframe for the current sample
      sc_columns <- c("sample_id", "barcode", "VDJ_contig_id", "VDJ_sequence_nt_raw", "isotype")
      sub_bulk <- bulk.tsv[, c(bulk.tsv.sample.column, bulk.tsv.barcode.column, "VDJ_contig_id", bulk.tsv.sequence.column, bulk.tsv.isotype.column)]
      colnames(sub_bulk) <- c("sample_id", "barcode", "VDJ_contig_id", "VDJ_sequence_nt_raw", "isotype")
      sample_combined <- rbind(sc.VDJ[sc.VDJ$sample_id == sample & sc.VDJ$VDJ_chain_count == 1, sc_columns], sub_bulk[sub_bulk$sample_id == sample,])
      # Process each row to annotate the combined data with IgBLAST annotations
      VDJ_combined_unclonotyped_unfiltered_sample <- lapply(1:nrow(sample_combined), function(row){
        # Extract barcode and contig name for the current sequence
        barcode <- sample_combined[row, "barcode"]
        contig_name <- sample_combined[row, "VDJ_contig_id"]
        # Determine dataset type (single-cell or bulk)
        if(grepl(pattern = "^bulk", contig_name)){dataset <- "bulk"}else{dataset <- "single-cell"}
        # Get the raw VDJ nucleotide sequence
        raw_sequence <- sample_combined[row, "VDJ_sequence_nt_raw"]
        # Retrieve IgBLAST annotations for the current sequence
        IgBLAST_annotations <- annotations_combined[annotations_combined$sequence == raw_sequence,]
        # If no annotations are found, skip this sequence
        if(nrow(IgBLAST_annotations) == 0){return(NULL)}
        # Check for completeness of IgBLAST annotations and skip incomplete entries
        if("" %in% c(IgBLAST_annotations$v_call, IgBLAST_annotations$j_call, IgBLAST_annotations$fwr1, IgBLAST_annotations$cdr1, IgBLAST_annotations$fwr2, IgBLAST_annotations$cdr2, IgBLAST_annotations$fwr3, IgBLAST_annotations$cdr3, IgBLAST_annotations$fwr4, IgBLAST_annotations$fwr1_aa, IgBLAST_annotations$cdr1_aa, IgBLAST_annotations$fwr2_aa, IgBLAST_annotations$cdr2_aa, IgBLAST_annotations$fwr3_aa, IgBLAST_annotations$cdr3_aa, IgBLAST_annotations$fwr4_aa)){return(NULL)}
        # Process V and J gene names by removing allele information to ensure consistency
        vgene <- gsub(pattern = "\\*.*$", replacement = "", strsplit(IgBLAST_annotations$v_call, split = ",")[[1]])
        jgene <- gsub(pattern = "\\*.*$", replacement = "", strsplit(IgBLAST_annotations$j_call, split = ",")[[1]])
        # If multiple gene names exist, check if one gene name is found within the others, and if a single consistent gene name is identified, return it
        if(length(vgene) > 1){if(all(sapply(vgene[2:length(vgene)], function(gene) grepl(pattern = vgene[1], gene)))){vgene <- vgene[1]}else{vgene <- unique(vgene)}}
        if(length(jgene) > 1){if(all(sapply(jgene[2:length(jgene)], function(gene) grepl(pattern = jgene[1], gene)))){jgene <- jgene[1]}else(jgene <- unique(jgene))}
        # Skip entry if V or J genes are not consistent
        if(length(vgene) != 1 | length(jgene) != 1){return(NULL)}
        # Initialize amino acid sequence variable
        aa_sequence <- NULL
        # Identify the start codon and translate the raw nucleotide sequence to amino acid sequence
        for(i in 1:(nchar(raw_sequence)-2)){
          if(substr(raw_sequence, start = i, stop = i+2) == "ATG"){
            aa_sequence <- suppressWarnings(as.character(Biostrings::translate(Biostrings::DNAString(substr(raw_sequence, start = i, stop = nchar(raw_sequence))), if.fuzzy.codon = "solve")))
            #Check if the sequence contains a stop codon
            if(!grepl(pattern = "\\*", aa_sequence)){start_codon_pos <- i-1; break}
            else{aa_sequence <- NULL}
          }
        }
        # Use IgBLAST annotation if translation is not possible
        if(is.null(aa_sequence)){aa_sequence <- IgBLAST_annotations$sequence_aa; start_codon_pos <- NULL}
        # Ensure the longest amino acid sequence is used
        if(nchar(aa_sequence) < nchar(IgBLAST_annotations$sequence_aa)){aa_sequence <- IgBLAST_annotations$sequence_aa; start_codon_pos <- NULL}


        #If FR1 is trimmed, trim this region from the germline and raw sequence
        if(trim.FR1){
          #Construct the trimmed sequences without FR1 by pasting together all annotated regions
          VDJ_sequence_nt_trimmed = paste0(IgBLAST_annotations$cdr1,IgBLAST_annotations$fwr2,IgBLAST_annotations$cdr2,
                                           IgBLAST_annotations$fwr3,IgBLAST_annotations$cdr3,IgBLAST_annotations$fwr4)
          VDJ_sequence_aa_trimmed = paste0(IgBLAST_annotations$cdr1_aa,IgBLAST_annotations$fwr2_aa,IgBLAST_annotations$cdr2_aa,
                                           IgBLAST_annotations$fwr3_aa,IgBLAST_annotations$cdr3_aa,IgBLAST_annotations$fwr4_aa)

          #Trim the germline sequence by aligning the CDR1 region to the raw germline sequence
          #NT sequence
          VDJ_germline_nt_raw = gsub(pattern = "N", replacement = "", IgBLAST_annotations$germline_alignment)
          alignment_nt <- pwalign::pairwiseAlignment(pattern = IgBLAST_annotations$cdr1, subject = VDJ_germline_nt_raw)
          cdr1_alignment <- as.character(pwalign::alignedPattern(alignment_nt))
          start_pos <- nchar(substr(cdr1_alignment, 1, regexpr("[^-]", cdr1_alignment)))
          VDJ_germline_nt_trimmed = substring(VDJ_germline_nt_raw, start_pos)
          #AA sequence
          VDJ_germline_aa_raw = gsub(pattern = "X", replacement = "", IgBLAST_annotations$germline_alignment_aa)
          alignment_aa <- pwalign::pairwiseAlignment(pattern = IgBLAST_annotations$cdr1_aa, subject = VDJ_germline_aa_raw)
          cdr1_alignment <- as.character(pwalign::alignedPattern(alignment_aa))
          start_pos <- nchar(substr(cdr1_alignment, 1, regexpr("[^-]", cdr1_alignment)))
          VDJ_germline_aa_trimmed = substring(VDJ_germline_aa_raw, start_pos)

        }else{
          #Construct the trimmed sequences with FR1 by pasting together all annotated regions
          VDJ_sequence_nt_trimmed = paste0(IgBLAST_annotations$fwr1,IgBLAST_annotations$cdr1,IgBLAST_annotations$fwr2,
                                           IgBLAST_annotations$cdr2,IgBLAST_annotations$fwr3,IgBLAST_annotations$cdr3,
                                           IgBLAST_annotations$fwr4)
          VDJ_sequence_aa_trimmed = paste0(IgBLAST_annotations$fwr1_aa,IgBLAST_annotations$cdr1_aa,IgBLAST_annotations$fwr2_aa,
                                           IgBLAST_annotations$cdr2_aa,IgBLAST_annotations$fwr3_aa,IgBLAST_annotations$cdr3_aa,
                                           IgBLAST_annotations$fwr4_aa)

          #Don't trim the germline sequences
          VDJ_germline_nt_raw = gsub(pattern = "N", replacement = "", IgBLAST_annotations$germline_alignment)
          VDJ_germline_nt_trimmed = NA
          VDJ_germline_aa_raw = gsub(pattern = "X", replacement = "", IgBLAST_annotations$germline_alignment_aa)
          VDJ_germline_aa_trimmed = NA

        }

        # Construct the row for the combined VDJ dataframe
        VDJ_combined_unclonotyped_row <- data.frame(sample_id = sample,
                                                    barcode = barcode,
                                                    dataset = dataset,
                                                    celltype = "B cell",
                                                    isotype = sample_combined[row, "isotype"],
                                                    VDJ_contig_id = contig_name,
                                                    VDJ_chain = "IGH",
                                                    VDJ_chain_count = 1,
                                                    VDJ_vgene = vgene,
                                                    VDJ_jgene = jgene,
                                                    VDJ_fwr1_nt = IgBLAST_annotations$fwr1,
                                                    VDJ_fwr1_aa = IgBLAST_annotations$fwr1_aa,
                                                    VDJ_cdr1_nt = IgBLAST_annotations$cdr1,
                                                    VDJ_cdr1_aa = IgBLAST_annotations$cdr1_aa,
                                                    VDJ_fwr2_nt = IgBLAST_annotations$fwr2,
                                                    VDJ_fwr2_aa = IgBLAST_annotations$fwr2_aa,
                                                    VDJ_cdr2_nt = IgBLAST_annotations$cdr2,
                                                    VDJ_cdr2_aa = IgBLAST_annotations$cdr2_aa,
                                                    VDJ_fwr3_nt = IgBLAST_annotations$fwr3,
                                                    VDJ_fwr3_aa = IgBLAST_annotations$fwr3_aa,
                                                    VDJ_cdr3_nt = IgBLAST_annotations$cdr3,
                                                    VDJ_cdr3_aa = IgBLAST_annotations$cdr3_aa,
                                                    VDJ_fwr4_nt = IgBLAST_annotations$fwr4,
                                                    VDJ_fwr4_aa = IgBLAST_annotations$fwr4_aa,
                                                    VDJ_sequence_nt_raw = raw_sequence,
                                                    VDJ_sequence_aa_raw = aa_sequence,
                                                    VDJ_sequence_nt_trimmed = VDJ_sequence_nt_trimmed,
                                                    VDJ_sequence_aa_trimmed = VDJ_sequence_aa_trimmed,
                                                    VDJ_consensus_nt_trimmed = NA,
                                                    VDJ_consensus_aa_trimmed = NA,
                                                    VDJ_germline_nt_raw = VDJ_germline_nt_raw,
                                                    VDJ_germline_aa_raw = VDJ_germline_aa_raw,
                                                    VDJ_germline_nt_trimmed = VDJ_germline_nt_trimmed,
                                                    VDJ_germline_aa_trimmed = VDJ_germline_aa_trimmed)
        # Return the combined VDJ dataframe row and annotations object
        return(VDJ_combined_unclonotyped_row)
      })
      # Remove null entries from the sample's combined data
      VDJ_combined_unclonotyped_unfiltered_sample <- VDJ_combined_unclonotyped_unfiltered_sample[!sapply(VDJ_combined_unclonotyped_unfiltered_sample, is.null)]
      # Combine the processed sample data into a single dataframe
      VDJ_combined_unclonotyped_unfiltered_sample <- do.call(rbind, VDJ_combined_unclonotyped_unfiltered_sample)
      # Return the combined sample data
      return(VDJ_combined_unclonotyped_unfiltered_sample)
    }))
  }

  # Function to calculate seqeunce identity between two sequences
  calc_seq_identity <- function(seq1, seq2){
    # Peform pairwise sequence alignment with default gap opening and extension costs
    alignment <- pwalign::pairwiseAlignment(pattern = seq1, subject = seq2)
    # Convert the aligned sequences (from pattern and subject) to character strings
    seq1_alignment <- as.character(alignment@pattern)
    seq2_alignment <- as.character(alignment@subject)
    # Split the aligned sequences into individual characters for comparison
    seq1_alignment <- strsplit(seq1_alignment, split = "")[[1]]
    seq2_alignment <- strsplit(seq2_alignment, split = "")[[1]]
    # Calculate the percentage of identical characters between the two aligned sequences (the number of matching characters divided by the length of the longer sequence, multiplied by 100)
    identity <- sum(seq1_alignment == seq2_alignment)/max(nchar(seq1), nchar(seq2))
    # Return the calculated sequence identity percentage
    return(identity)
  }

  # Function to assign each bulk sequence to clonotype from the single-cell data
  Bulk_Clonotyping <- function(VDJ_combined, tie.resolvement, seq.identity){

    # Process each sample separately
    bulkRNA_clonotyping <- lapply(unique(VDJ_combined$sample_id), function(sample){
    # Subset the VDJ dataframe for this sample
    VDJ_sample <- VDJ_combined[VDJ_combined$sample_id == sample,]

    # Extract unique bulk sequences for the current sample
    unique_bulkRNA_sequences <- unique(VDJ_sample[VDJ_sample$dataset == "bulk",
                                                    c("VDJ_sequence_nt_raw", "VDJ_vgene", "VDJ_jgene", "VDJ_cdr3_nt")])

    # Match bulk sequences to single-cell clonotypes
    clonotype_ids <- lapply(1:nrow(unique_bulkRNA_sequences), function(row){

      # Extract sequence and gene information for the current bulk sequence
      bulk_seq_raw <- unique_bulkRNA_sequences[row, "VDJ_sequence_nt_raw"]
      vgene <- unique_bulkRNA_sequences[row, "VDJ_vgene"]
      jgene <- unique_bulkRNA_sequences[row, "VDJ_jgene"]
      cdr3 <- unique_bulkRNA_sequences[row, "VDJ_cdr3_nt"]

      # Find matching single-cell sequences with the same V and J genes and CDR3 length
      matching_scRNA_sequences <- VDJ_sample[VDJ_sample$dataset == "single-cell" &
                                               VDJ_sample$VDJ_vgene == vgene &
                                               VDJ_sample$VDJ_jgene == jgene &
                                               nchar(VDJ_sample$VDJ_cdr3_nt) == nchar(cdr3),]

      # If no matching single-cell sequences are found, return NA
      if(nrow(matching_scRNA_sequences) == 0){return(NA)}
      # If a single matching single-cell sequence is found, return the clonotype ID
      else if(nrow(matching_scRNA_sequences) == 1){return(matching_scRNA_sequences$clonotype_id)}
      # If multiple matching single-cell sequences are found, proceed with clonotyping
      else if(nrow(matching_scRNA_sequences) > 1){
        # Extract clonotypes with identical CDR3 sequences
        clones_identical_cdr3 <- unique(matching_scRNA_sequences[matching_scRNA_sequences$VDJ_cdr3_nt == cdr3,"clonotype_id"])
        # If only one clonotype has an identical CDR3 sequence, return the clonotype ID
        if(length(clones_identical_cdr3) == 1){return(clones_identical_cdr3)}
        # If multiple clonotypes have identical CDR3 sequences, proceed with clonotyping
        else if(length(clones_identical_cdr3) > 1){
          # Select the clonotype(s) with the highest frequency
          clones_identical_cdr3_sizes <- sapply(clones_identical_cdr3, function(clone) unique(matching_scRNA_sequences[matching_scRNA_sequences$clonotype_id == clone, "clonotype_frequency"]))
          clones_identical_cdr3 <- clones_identical_cdr3[clones_identical_cdr3_sizes == max(clones_identical_cdr3_sizes)]
          # If only one clonotype has the highest frequency, return the clonotype ID
          if(length(clones_identical_cdr3) == 1){return(clones_identical_cdr3)}
          # If multiple clonotypes have the highest frequency, proceed with clonotyping
          else if(length(clones_identical_cdr3) > 1){
            # Select the clonotype(s) with sequences that have the minimum Levenshtein distance to the current bulk transcripts
            clones_identical_cdr3_dists <- sapply(clones_identical_cdr3, function(clone){
              clone_seqs <- VDJ_sample[VDJ_sample$dataset == "single-cell"& VDJ_sample$clonotype_id == clone, "VDJ_sequence_nt_raw"]
              mean_dist <- mean(sapply(clone_seqs, function(sc_seq_raw){stringdist::stringdist(bulk_seq_raw, sc_seq_raw, method = "lv")}))
            })
            clones_identical_cdr3 <- clones_identical_cdr3[clones_identical_cdr3_dists == min(clones_identical_cdr3_dists)]
            # If only one clonotype has the minimum Levenshtein distance, return the clonotype ID
            if(length(clones_identical_cdr3) == 1){return(clones_identical_cdr3)}
            # If multiple clonotypes remain, perform tie resolvement
            else if(length(clones_identical_cdr3) > 1){
              if(tie.resolvement == "all"){print(paste0("Multiple clonotypes match: ", bulk_seq_raw));return(paste(clones_identical_cdr3, collapse = "_"))}
              else if(tie.resolvement == "none"){return(NA)}
              else if(tie.resolvement == "random"){return(sample(clones_identical_cdr3, 1))}
          }
        }
        }
        # If no clonotypes have identical CDR3 sequences, proceed with identity-based clonotyping
        else if(length(clones_identical_cdr3) == 0){
          # Extract all unique clonotypes in the matching single-cell sequences
          all_clones <- unique(matching_scRNA_sequences$clonotype_id)

          # Initialize an empty vector for clonotypes with CDR3 sequences that have at least the threshold sequence identity
          clones_identity <- c()

          # For each clone, check sequence identity between bulk CDR3 and single-cell CDR3 sequences
          for(clone in all_clones){
            clones_identity <- c(clones_identity, clone)
            for(seq in matching_scRNA_sequences[matching_scRNA_sequences$clonotype_id == clone, "VDJ_cdr3_nt"]){
              identity_score <- calc_seq_identity(seq1 = cdr3, seq2 = seq)
              if(identity_score < seq.identity){clones_identity <- clones_identity[!clones_identity == clone]; break}
            }
          }

          # If only one clonotype meets the sequence identity threshold, return the clonotype ID
          if(length(clones_identity) == 1){return(clones_identity)}
          # If no clontype meets the sequence identitiy threshold, return NA
          else if(length(clones_identity) == 0){return(NA)}
          # If multiple clonotypes meet the sequence identity threshold, proceed with clonotyping
          else if(length(clones_identity) > 1){
            # Select the clonotype(s) with the highest frequency
            clones_identity_sizes <- sapply(clones_identity, function(clone) unique(matching_scRNA_sequences[matching_scRNA_sequences$clonotype_id == clone, "clonotype_frequency"]))
            clones_identity <- clones_identity[clones_identity_sizes == max(clones_identity_sizes)]
            # If only one clonotype has the highest frequency, return the clonotype ID
            if(length(clones_identity) == 1){return(clones_identity)}
            # If multiple clonotypes have the highest frequency, proceed with clonotyping
            else if(length(clones_identity) > 1){
              # Select the clonotype(s) with sequences that have the minimum Levenshtein distance to the current bulk transcripts
              clones_identity_dists <- sapply(clones_identity, function(clone){
                clone_seqs <- VDJ_sample[VDJ_sample$dataset == "single-cell" & VDJ_sample$clonotype_id == clone, "VDJ_sequence_nt_raw"]
                mean_dist <- mean(sapply(clone_seqs, function(sc_seq_raw){stringdist::stringdist(bulk_seq_raw, sc_seq_raw, method = "lv")}))
              })
              clones_identity <- clones_identity[clones_identity_dists == min(clones_identity_dists)]
              # If only one clonotype has the minimum Levenshtein distance, return the clonotype ID
              if(length(clones_identity) == 1){return(clones_identity)}
              # If multiple clonotypes remain, perform tie resolvement
              else if(length(clones_identity) > 1){
                if(tie.resolvement == "all"){print(paste0("Multiple clonotypes match: ", bulk_seq_raw));return(paste(clones_identity, collapse = "_"))}
                else if(tie.resolvement == "none"){return(NA)}
                else if(tie.resolvement == "random"){return(sample(clones_identity, 1))}
              }
            }
          }
        }
      }
    })

    print(paste0("No matching clonotypes were found for ", sum(is.na(clonotype_ids)), " out of ", nrow(unique_bulkRNA_sequences)," bulk sequences of sample ", sample))

    # Name the clonotype ID list by the bulk sequences
    names(clonotype_ids) <- unique_bulkRNA_sequences$VDJ_sequence_nt_raw
    # Make dataframe
    bulk_clonotypes <- data.frame(VDJ_sequence_nt_raw = unique_bulkRNA_sequences$VDJ_sequence_nt_raw,
                                  clonotype_id = unlist(clonotype_ids),
                                  clonotype_frequency = NA)
    # If multiple clonotypes were assigned (when tie.resolvement is "all") copy the sequence
    for (row in seq(1:nrow(bulk_clonotypes))){
      if(grepl(pattern = "_", bulk_clonotypes[row, "clonotype_id"])){
        clonotypes <- unlist(strsplit(bulk_clonotypes[row, "clonotype_id"], split = "_"))
        for(clonotype in clonotypes){
          bulk_clonotypes <- rbind(bulk_clonotypes, data.frame(VDJ_sequence_nt_raw = bulk_clonotypes[row, "VDJ_sequence_nt_raw"],
                                                               clonotype_id = clonotype,
                                                               clonotype_frequency = NA))
        }
      }
    }
    #Remove the rows with multiple clonotypes
    bulk_clonotypes <- bulk_clonotypes[grep(pattern = "_", bulk_clonotypes[, "clonotype_id"], invert = T),]

    # Add the clonotypes of the bulk sequences to the VDJ of this sample
    VDJ_sample_bulk <- VDJ_sample[VDJ_sample$dataset == "bulk",!(colnames(VDJ_sample) %in% c("clonotype_id", "clonotype_frequency"))]
    VDJ_sample_bulk <- dplyr::left_join(VDJ_sample_bulk, bulk_clonotypes, by = "VDJ_sequence_nt_raw")
    VDJ_sample <- rbind(VDJ_sample[VDJ_sample$dataset == "single-cell",], VDJ_sample_bulk)

    #remove NA clonotypes
    VDJ_sample <- VDJ_sample[!is.na(VDJ_sample$clonotype_id),]

    # Recalculate the clonotype frequency
    VDJ_sample %>% dplyr::group_by(clonotype_id) %>% dplyr::mutate(clonotype_frequency = dplyr::n()) -> VDJ_sample

    # Assign germline genes to the bulk sequences
    for (clonotype in unique(VDJ_sample$clonotype_id)){
      #Most abundant germline sequence
      nt_germline <- sort(table(unique(VDJ_sample[VDJ_sample$clonotype_id == clonotype, "VDJ_germline_nt_trimmed"])), decreasing = T)[1]
      aa_germline <- sort(table(unique(VDJ_sample[VDJ_sample$clonotype_id == clonotype, "VDJ_germline_aa_trimmed"])), decreasing = T)[1]
      VDJ_sample[VDJ_sample$clonotype_id == clonotype, "VDJ_germline_nt_trimmed"] <- names(nt_germline)
      VDJ_sample[VDJ_sample$clonotype_id == clonotype, "VDJ_germline_aa_trimmed"] <- names(aa_germline)
    }


    # Return the VDJ for the current sample
    return(VDJ_sample)
    })


    # Combine the clonotype assignments for all samples
    bulkRNA_clonotyping <- do.call(rbind, bulkRNA_clonotyping)

    return(bulkRNA_clonotyping)
  }


  #Read the bulk data table
  bulk.tsv <- as.data.frame(utils::read.table(bulk.tsv, header = TRUE, sep = "\t"))
  #Only keep bulk samples that are present in the single cell data
  bulk.tsv <- bulk.tsv[bulk.tsv[,bulk.tsv.sample.column] %in% sc.VDJ$sample_id,]
  #Add contig ID column
  bulk.tsv$VDJ_contig_id <- paste0("bulk_",bulk.tsv[,bulk.tsv.sample.column], "_", bulk.tsv[,bulk.tsv.barcode.column])

  #If no annotated bulk or sc data is provided, run IgBLAST on the unique sequences
  if(is.null(bulkRNA_seqs_annotations) | is.null(scRNA_seqs_annotations)){
    #1. Create fasta files of the unique single cell and bulk sequence
    #Create the temp directory
    temp_dir <- './temp_igblast'
    if(!dir.exists(temp_dir)) dir.create(temp_dir)
    #Write the unique single cell sequences to a fasta file
    sc_transcripts <- as.list(unique(sc.VDJ$VDJ_sequence_nt_raw))
    seqinr::write.fasta(sequences = sc_transcripts,
                        names = 1:length(sc_transcripts),
                        file.out = paste0(temp_dir,"/unique_scRNA_seqs.fasta"))
    #Write the unique bulk transcripts to a FASTA file
    bulk_transcripts <- as.list(unique(bulk.tsv[,bulk.tsv.sequence.column]))
    seqinr::write.fasta(sequences = bulk_transcripts,
                        names = 1:length(bulk_transcripts),
                        file.out = paste0(temp_dir,"/unique_bulkRNA_seqs.fasta"))


   # 2. Run IgBLAST to annotate the single cell and bulk sequences
    #For single-cell data
    print("Running IgBLAST on the single-cell sequences, this may take several hours for large datasets.")

    console_command <- paste0(os_prefix, 'AssignGenes.py igblast -s ', paste0(temp_dir,"/unique_scRNA_seqs.fasta"),
                                " -b ", igblast.dir," --organism ", organism, " --loci ig --format airr")
    system(console_command)
    #For bulk data
    print("Running IgBLAST on the bulk sequences, this may take several hours for large datasets.")
    console_command <- paste0(os_prefix, 'AssignGenes.py igblast -s ', paste0(temp_dir,"/unique_bulkRNA_seqs.fasta"),
                                " -b ", igblast.dir," --organism ", organism, " --loci ig --format airr")
    system(console_command)

    #3. Build combined VDJ dataframe with IgBLAST annotations
    # Read the IgBLAST annotations for single-cell and bulk sequences
    scRNA_seqs_annotations <- utils::read.table(paste0(temp_dir, "/unique_scRNA_seqs_igblast.tsv"), header = TRUE, sep = "\t")
    bulkRNA_seqs_annotations <- utils::read.table(paste0(temp_dir, "/unique_bulkRNA_seqs_igblast.tsv"), header = TRUE, sep = "\t")
    unlink(temp_dir, recursive = T)
  }
  else{
    # Read the IgBLAST annotations for single-cell and bulk sequences
    scRNA_seqs_annotations <- utils::read.table(scRNA_seqs_annotations, header = TRUE, sep = "\t")
    bulkRNA_seqs_annotations <- utils::read.table(bulkRNA_seqs_annotations, header = TRUE, sep = "\t")
  }


  # Combine single-cell and bulk IgBLAST annotations into a single dataframe
  annotations_combined <- unique(rbind(scRNA_seqs_annotations, bulkRNA_seqs_annotations))
  print(paste0("Start Transform_to_VDJ() ", Sys.time()))
  VDJ_combined <- Transform_to_VDJ(annotations_combined, sc.VDJ, bulk.tsv, bulk.tsv.sample.column, bulk.tsv.barcode.column,
                                   bulk.tsv.sequence.column, bulk.tsv.isotype.column, trim.FR1)
  print(paste0("End Transform_to_VDJ() ", Sys.time()))

  #4. Append the original clonotype IDs from the single-cell VDJ dataframe to the VDJ_OVA_combined_unfiltered dataframe
  VDJ_combined <- dplyr::left_join(VDJ_combined[,!(colnames(VDJ_combined) %in% c("clonotype_id", "clonotype_frequency"))],
                                   sc.VDJ[,c("sample_id", "barcode", "VDJ_contig_id", "clonotype_id", "clonotype_frequency")],
                                   by = c("sample_id", "barcode", "VDJ_contig_id"))

  #5. Assign the bulk transcript to the single-cell clonotypes
  print(paste0("Start Bulk_Clonotyping() ", Sys.time()))
  VDJ_clonotyped <- Bulk_Clonotyping(VDJ_combined, tie.resolvement, seq.identity)
  print(paste0("End Bulk_Clonotyping() ", Sys.time()))

  return(VDJ_clonotyped)
}

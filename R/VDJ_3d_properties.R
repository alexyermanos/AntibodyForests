#' Function to calculate 3D-structure propoperties such as the average charge, hydrophobicity, pKa shift and free energy of PDB files and add them to an AntibodyForests-object
#' @description Function to calculate protein 3D-structure properties of antibodies (or antibody-antigen complexes) and integrate them into an AntibodyForests-object.
#' @param VDJ a dataframe with V(D)J information such as the output of Platypus::VDJ_build(). Must contain columns sample_id, clonotype_id, barcode.
#' @param pdb.dir a directory containing PDB files.
#' @param file.df a dataframe of filenames (column file_name) to be used and sequence IDs (column sequence) corresponding to the the barcodes column of the VDJ dataframe.
#' @param properties a vector of properties to be calculated. Default is c("charge", "hydrophobicity").
#' - charge
#' - hydrophobicity
#' - pKa_shift: the acid dissociation constant shift upon binding of the antibody to the antigen (needs Propka output)
#' - free_energy: the free energy of binding of the antibody to the antigen at a certain pH (needs Propka output)
#' @param free_energy_pH the pH to be used to calculate the free energy of binding. Default is 7.
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
#' @param propka.dir a directory containing Propka output files. The propka filenames should be similar to the PDB filenames.


VDJ_3d_properties <- function(VDJ,
                              pdb.dir,
                              file.df,
                              properties,
                              sequence.region,
                              chain,
                              propka.dir,
                              free_energy_pH,
                              sub.sequence.column){
  
  #Check for missing or false input
  if(missing(VDJ)){stop("VDJ dataframe is missing")}
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
  if(missing(properties)){properties = c("charge", "hydrophobicity")}
  if(!all(properties %in% c("charge", "hydrophobicity", "pKa_shift", "free_energy"))){stop("properties must be one or more of 'charge', 'hydrophobicity', 'pKa_shift', 'free_energy'")}
  if(missing(propka.dir)){propka.dir = NULL} 
  if(is.null(propka.dir) & any(c("pKa_shift", "free_energy") %in% properties)){stop("propka.dir is missing")}
  if(!(is.null(propka.dir))){if(!dir.exists(propka.dir)){stop("propka.dir does not exist")}}
  if(missing(free_energy_pH)){free_energy_pH = 7}
  if(!is.numeric(free_energy_pH)){stop("free_energy_pH must be numeric")}
  if(missing(sub.sequence.column)){sub.sequence.column = NULL}
  if(!is.null(sub.sequence.column)){if(!sub.sequence.column %in% colnames(VDJ)){stop("sub.sequence.column does not exist in VDJ dataframe")}}
  if(sequence.region == "sub.sequence" & is.null(sub.sequence.column)){stop("sub.sequence.column is missing")}
    
  
  #Based on Lucas' VDJ_structure_analysis function
  calculate_charge <- function(pdb, chains, resnos){
    df_charge <- data.frame()
    pdb$atom$charge <- NULL
    #Calculate the charge of each residue per chain
    for(CHAIN in chains){
      df_chain <- pdb$atom |> dplyr::filter(chain == CHAIN) |> dplyr::distinct(resno, .keep_all = TRUE) |> dplyr::select(resno, chain)
      df_chain$charge <- pdb$atom |> dplyr::filter(chain == CHAIN) |> dplyr::distinct(resno, .keep_all = TRUE) |> dplyr::pull(resid) |> stringr::str_to_title() |> seqinr::a() |> Peptides::charge()
      df_charge <- rbind(df_charge, df_chain)
    }
    pdb$atom <- dplyr::left_join(pdb$atom, df_charge, by = c("resno","chain"))
    #Calculate the average charge depending on the residues supplied
    average_charge = mean(na.omit(pdb$atom[pdb$atom$resno %in% resnos,]$charge))
    return(average_charge)
  }
  
  #Based on Lucas' VDJ_structure_analysis function
  calculate_hydrophobicity <- function(pdb, chains, resnos){
    df_hydph <- data.frame()
    #Calculate the hydrophobicity of each residue per chain
    for(CHAIN in chains){
      df_chain <- pdb$atom |> dplyr::filter(chain == CHAIN) |> dplyr::distinct(resno, .keep_all = T) |> dplyr::select(resno,chain)
      df_chain$hydrophobicity <- pdb$atom |> dplyr::filter(chain == CHAIN) |> dplyr::distinct(resno, .keep_all = T) |> dplyr::pull(resid) |> stringr::str_to_title() |> seqinr::a() |> Peptides::hydrophobicity()
      df_hydph <- rbind(df_hydph,df_chain)
    }
    pdb$atom <- dplyr::left_join(pdb$atom, df_hydph, by = c("resno","chain"))
    #Calculate the average charge depending on the residues supplied
    average_hydrophobicity = mean(na.omit(pdb$atom[pdb$atom$resno %in% resnos,]$hydrophobicity))
    return(average_hydrophobicity)
  }
  
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
    dist_mat <- dist.xyz(coord_binder, coord_ligand)
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
    
    return(resnos)
  }
  
  #From Steropodon
  get_propka_df <- function(propka.dir, filename, chains, resnos){

    start <- grep('SUMMARY OF THIS PREDICTION', propka_out)
    end <- grep('Free energy of', propka_out)
    pka <- propka_out[(start+1):(end-3)]
    pka <- paste(pka, collapse = '\n')
    tc <- textConnection(pka)
    df <- utils::read.table(tc, as.is=T, fill=T, blank.lines.skip=F)
    colnames(df) <- c("aa", "position", "chain", "pKa", "model-pKa")
    df <- df[2:nrow(df),]
    df$position <- as.numeric(as.character(df$position))
    
    #Filter on chain and residues
    df <- df[df$chain %in% chains & df$position %in% resnos,]
    return(df)
  }
  
  calculate_pKa_shift <- function(propka_out, filename, chains, resnos){
    #Create propka dataframe
    start <- grep('SUMMARY OF THIS PREDICTION', propka_out)
    end <- grep('Free energy of', propka_out)
    pka <- propka_out[(start+1):(end-3)]
    pka <- paste(pka, collapse = '\n')
    tc <- textConnection(pka)
    df <- utils::read.table(tc, as.is=T, fill=T, blank.lines.skip=F)
    colnames(df) <- c("aa", "position", "chain", "pKa", "model-pKa")
    df <- df[2:nrow(df),]
    df$position <- as.numeric(as.character(df$position))
    
    #Filter on chain and residues
    df <- df[df$chain %in% chains & df$position %in% resnos,]
    
    #Calculate the pKa shift
    df[,c("pKa", "model.pKa")] <- lapply(df[,c("pKa", "model-pKa")], as.numeric)
    df$pKa_shift <- df$pKa - df$model.pKa
    average_pKa_shift <- mean(df$pKa_shift)
    return(average_pKa_shift)
  }
  
  calculate_free_energy <- function(propka_out, free_energy_pH){
    start <- grep('Free energy of', propka_out)
    end <- grep('The pH of optimum stability', propka_out)
    en <- propka_out[(start+1):(end-2)]
    en <- paste(en, collapse = '\n')
    tc <- textConnection(en)
    df <- utils::read.table(tc, as.is=T, fill=T, blank.lines.skip=F)
    close(tc)
    free_energy <- df[df$V1 == free_energy_pH, "V2"]
    return(free_energy)
  }

  
  #Get list of files if file.df is null
  if(is.null(file.df)){file.df = list.files(pdb.dir)}
  
  #Add the empty properties columns to the VDJ dataframe
  if("charge" %in% properties){VDJ$average_charge <- NA}
  if("hydrophobicity" %in% properties){VDJ$average_hydrophobicity <- NA}
  if("pKa_shift" %in% properties){VDJ$average_pKa_shift <- NA}
  if("free_energy" %in% properties){VDJ$free_energy <- NA}
  
  #Set the chains
  if(chain == "HC+LC"){chains = c("A", "B")}
  if(chain == "HC"){chains = "A"}
  if(chain == "LC"){chains = "B"}
  if(chain == "AG"){chains = "C"}
  if(chain == "whole.complex"){chains = c("A", "B", "C")}
  

  
  #Calculate properties for each structure
  for(row in seq(1:nrow(file.df))){
    file = file.df[row,]
    if(file$sequence %in% VDJ$barcode){
      #Read in the pdb file
      pdb <- bio3d::read.pdb(paste0("~/OneDrive - UMC Utrecht/Documenten/Structure_Daphne/PDBS_superimposed/", file$file_name))
      
      #Get the residues for calculations of average properties
      if(sequence.region == "full.sequence"){
        resnos <- unique(pdb$atom[pdb$atom$chain %in% chains, ]$resno)
      }
      if(sequence.region == "binding.residues"){
        resnos <- get_bs_res(pdb, chains_binder = c("A", "B"), chains_ligand = "C")
      }
      if(sequence.region == "sub.sequence"){
        #get the subsequence
        sub_seq <- VDJ[VDJ$barcode == file$sequence,sub.sequence]
        resnos <- get_subseq_res(pdb, sub_seq, chains)
      }
      
      if("charge" %in% properties){
        #Calculate average charge
        average_charge <- calculate_charge(pdb, chains, resnos)
        #Add to the VDJ dataframe
        VDJ[VDJ$barcode == file$sequence,]$average_charge <- average_charge
      }
      if("hydrophobicity" %in% properties){
        #Calculate average hydrophobicity
        average_hydrophobicity <- calculate_hydrophobicity(pdb, chains, resnos)
        #Add to the VDJ dataframe
        VDJ[VDJ$barcode == file$sequence,]$average_hydrophobicity <- average_hydrophobicity
      }
      if("pKa_shift" %in% properties | "free_energy" %in% properties){
        #Read in the propka output
        propka_out <- readLines(paste0(propka.dir, "/", substr(file$file_name, 1, nchar(file$file_name)-4),".pka"))
        
        if("pKa_shift" %in% properties){
          #Calculate the pKa shift
          average_pKa_shift <- calculate_pKa_shift(propka_out, file$file_name, chains, resnos)
          #Add to the VDJ dataframe
          VDJ[VDJ$barcode == file$sequence,]$average_pKa_shift <- average_pKa_shift
        }
        if("free_energy" %in% properties){
          #Calculate the free energy
          free_energy <- calculate_free_energy(propka_out, free_energy_pH)
          #Add to the VDJ dataframe
          VDJ[VDJ$barcode == file$sequence,]$free_energy <- free_energy
        }
      }


    } else {
      print(paste0("Sequence ", file$sequence, " not found in VDJ dataframe"))
      next
    }
    

    
  }
  return(VDJ)
}

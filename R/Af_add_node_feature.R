#' Function to add node features to an AntibodyForests-object
#' @description Function to add node features to an AntibodyForests-object
#' @param AntibodyForests_object AntibodyForests-object, output from Af_build()
#' @param feature.df Dataframe with features for each node. Must contain columns sample_id, clonotype_id, barcode and the features to be added.
#' @param feature.names Character vector with the names of the features to be added.
#' @return Returns an AntibodyForests-object with the features added to the nodes.
#' @export
#' @examples
#' af <- Af_add_node_feature(AntibodyForests::small_af,
#'                           feature.df = AntibodyForests::small_vdj,
#'                           feature.names = c("VDJ_dgene", "VDJ_jgene"))

Af_add_node_feature <- function(AntibodyForests_object,
                                             feature.df,
                                             feature.names){

  #Stop when no input is provided
  if(missing(AntibodyForests_object)){stop("Please provide a valid AntibodyForests_object input object.")}
  if(missing(feature.df)){stop("Please provide a valid feature.df argument.")}
  if(missing(feature.names)){stop("Please provide a valid feature.names argument.")}
  #Check if feature.names is in feature.df
  if(!all(feature.names %in% colnames(feature.df))){stop("Not all feature.names are present in the feature.df.")}


  out <- lapply(seq_along(AntibodyForests_object), function(sample){
    sample_name <- names(AntibodyForests_object)[[sample]]

    #Start list of clonotype names for this sample
    clonotype_name_list <- names(AntibodyForests_object[[sample_name]])
    out_list <- lapply(seq_along(AntibodyForests_object[[sample]]), function(clonotype){
      clonotype_name <- names(AntibodyForests_object[[sample_name]])[[clonotype]]

      #Subset the feature.df for this specific sample and clonotype
      feature.df_sub <- feature.df[feature.df$sample_id == sample_name & feature.df$clonotype_id == clonotype_name,]

      AntibodyForests_object[[sample]][[clonotype]]$nodes <- lapply(AntibodyForests_object[[sample]][[clonotype]]$nodes, function(node){
        for(feature in feature.names){
          #Subset the feature.df for this specific node
          feature.df_node_sub <- feature.df_sub[feature.df_sub$barcode %in% node$barcode,]

          #Match the order of barcodes
          if(length(unique(node$barcode)) > 1){
            feature.df_node_sub$barcode = factor(feature.df_node_sub$barcode, levels = node$barcode)
            feature.df_node_sub[order(feature.df_node_sub$barcode), ]
          }

          #Add node feature
          node[[feature]] <- feature.df_node_sub[,feature]
        }
        return(node)
      })
      return(AntibodyForests_object[[sample]][[clonotype]])
    })
    #Set the names of the clonotypes
    names(out_list) <- clonotype_name_list
    return(out_list)
  })
  #Set the names of the samples
  names(out) <- names(AntibodyForests_object)

  #Change class type into AntibodyForests-object
  antibodyforests_object <- base::structure(out, class = "AntibodyForests")

  return(antibodyforests_object)
}

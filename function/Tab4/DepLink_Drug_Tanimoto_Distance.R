#Calculate drug similarity

# fps.bits.idx1_list = pubchem_maccs_bits_idx #list format
# nDrug = nrow(drug)
# drugName = drug$broad_id

DepLink_Tanimoto_Distance_Matrix <- function(fps.bits.idx1_list, nDrug, drugName){
  
    simmat <- diag(nDrug) #create a matrix
    rownames( simmat) <- drugName
    colnames( simmat) <- drugName
    
  for (i in 1:nDrug) {
    for (j in i:nDrug) {
      fps1 <- fps.bits.idx1_list[[i]]
      fps2 <- fps.bits.idx1_list[[j]]
      tmp <- length(intersect(fps1,fps2))/length(union(fps1,fps2))
      simmat[i, j] <- tmp
      simmat[j, i] <- tmp
      
    }
    cat("col: ",j,"/ row: ", i, "\n")
  }
  
 
  return(simmat)
  
}
# 
# drug.sim.matrix <- DepLink_Tanimoto_Distance_Matrix(fps.bits.idx1_list = pubchem_maccs_bits_idx ,
#                        nDrug = length(pubchem_maccs_bits_idx),
#                        drugName = drug$broad_id)
# 
# saveRDS(drug.sim.matrix ,file = "data/primary-screen-replicate-collapsed-treatment-info_19Q4_4684drug_similarity_tanimoto.RDS")

#individual
# fps.bits.idx1 = pubchem_maccs_bits_idx[1:2]   #list format
# compared.idx2 = pubchem_maccs_bits_idx      #list format compared
# nDrug.input = 2 #how many drugs
# drugName.input = drug$broad_id[1:2]  #new drug names
# combineMatrix = drug.sim.matrix 

DepLink_Tanimoto_Distance_Individual <- function(fps.bits.idx1,
                             compared.idx2 ,
                             nDrug.input,
                             drugName.input,
                             combineMatrix) {
    simmat <- as.data.frame(matrix(data = NA, nrow = nrow(combineMatrix)+1, ncol = nDrug.input), stringsAsFactors = F)
    
    broad_id <- rownames(combineMatrix)
    rownames(simmat) <- c(broad_id,drugName.input)
    colnames(simmat) <- drugName.input
    
    
    if (!nDrug.input == 1) {
      for (i in 1:nDrug.input) {
        for (j in 1:length(compared.idx2)) {
          fps1 <- fps.bits.idx1[[i]]
          fps2 <- compared.idx2[[j]]
          tmp <- length(intersect(fps1, fps2)) / length(union(fps1, fps2))
          simmat[j, i] <- tmp
        }
      }
    } else{
      
      for (i in 1:length(compared.idx2)) {
        fps1 <- fps.bits.idx1[[1]]
        fps2 <- compared.idx2[[i]]
        tmp <- length(intersect(fps1, fps2)) / length(union(fps1, fps2))
        simmat[i, 1] <- tmp
      }
    }
    

    broad_id <- rownames(simmat)
    simmat <- cbind(broad_id, simmat)
    simmat <- simmat[which(!simmat$YourDrug == 1), ] #Remove correlation = 1
    return(simmat)
}

# 
#  newDrug <- DepLink_Tanimoto_Distance_Individual(fps.bits.idx1 = pubchem_maccs_bits_idx[1:2] ,
#                  compared.idx2 = pubchem_maccs_bits_idx ,
#                  nDrug.input = 2,
#                  drugName.input = drug$broad_id[1:2],
#                  combineMatrix = drug.sim.matrix )





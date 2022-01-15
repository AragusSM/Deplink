# Functions for DepLink Gene Query function.
# All plots and table generate functions are included in this document.
#    Project".: DepLink
#    Contributed by: Li-Ru Wang, Tapsya Nayak, Yu-Chiao Chiu, and Yidong Chen
#    Started from: October 2021
#    Affiliation: Greehey Children's Cancer Research Institute, 
#                 The University of Texas Health San Antonio
#


#######################################
#This function project the new drug in the PRISM TSNE object
#TSNE plot for new drug
DepLink_TSNEplot <- function(newDrug.similarity,drugTable,
                             drug.similarity,drug.tsne, target.drug.id = NULL){
  require(plotly)
  require(gplots)
  require(RColorBrewer)

 
  
  if (sum(newDrug.similarity$`Tanimoto similarity` == 1) == 1) {
    idx <- which(newDrug.similarity$`Tanimoto similarity`== 1)
    NewDrug.tsne.related <- drug.tsne[which(drug.tsne$broad_id %in% newDrug.similarity$broad_id[idx] ), ]
    similarDrug <- drug.tsne[which(drug.tsne$broad_id %in% newDrug.similarity$broad_id[2:4]), ]
  }else{
    similarDrug <- drug.tsne[which(drug.tsne$broad_id %in% newDrug.similarity$broad_id[1:3]), ]
    
    #Calculate middle point among three points
    NewDrug.tsne.related <- data.frame(X1 = mean(similarDrug$X1),
                                       X2 = mean(similarDrug$X2),
                                       X3 = mean(similarDrug$X3),
                                       drug_name = "Query drug")
}
    #Plot

    fig <-  plot_ly(data = drug.tsne[,1:4] ,x =  ~X1, y = ~X2, z = ~X3,
                    text = ~drug_name,
                    hoverinfo = "text", color = I("gray"),
                    name = "PRISM: 4684 Drug") %>% 
      add_markers(size = 6) %>%
      layout(scene = list(xaxis = list(title = 'tSNE1'),
                          yaxis = list(title = 'tSNE2'),
                          zaxis = list(title = 'tSNE3')))
           
    
    fig <- fig %>% add_trace(
      name = "Query drug",
      data = NewDrug.tsne.related,
      x = ~X1,
      y = ~X2,
      z = ~X3,
      text = ~drug_name,
      hoverinfo = "text",
      mode = "markers",
      type = "scatter3d",
      marker = list(size = 6, color = "red"))

    fig <- fig %>% add_trace(
      name = "Nearest neighbor",
      data = similarDrug,
      x = ~X1,
      y = ~X2,
      z = ~X3,
      text = ~drug_name,
      hoverinfo = "text",
      mode = "markers",
      type = "scatter3d",
      marker = list(size = 6, color = "green",symbol = "cross"))
    
    if (!is.null(target.drug.id)) {
      fig <- fig %>% add_trace(
        name = "Selected Drug",
        data = drug.tsne[which(drug.tsne$broad_id == target.drug.id),1:4],
        x = ~X1,
        y = ~X2,
        z = ~X3,
        text = ~drug_name,
        hoverinfo = "text",
        mode = "markers",
        type = "scatter3d",
        marker = list(size = 7, color = "yellow",symbol = "diamond",
                      line = list(width = 0.5, color = "black")))
    }
    
    return(fig) 
  
}


##################################
#This function calculate the tanimoto distance between new drug and PRISM drugs

DepLink_Tanimoto_Distance_Individual <- function(fps.bits.idx1,
                                                 compared.idx2 ,
                                                 nDrug.input,
                                                 drugName.input,
                                                 combineMatrix) {
  simmat <- as.data.frame(matrix(data = NA, nrow = nrow(combineMatrix), ncol = nDrug.input+3), stringsAsFactors = F)
  
  broad_id <- rownames(combineMatrix)
  rownames(simmat) <- broad_id

  colnames(simmat) <- c(drugName.input,"fps1","fps2","intersct")
  
  
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
      simmat[i, "fps1"] <- length(fps1)
      simmat[i, "fps2"] <- length(fps1)
      simmat[i, "intersct"] <- length(intersect(fps1, fps2))
    }
  }
  
  
  #broad_id <- rownames(simmat)
  simmat <- cbind(broad_id, simmat)
  
  
  #simmat <- simmat[which(!simmat$YourDrug == 1), ] #Remove correlation = 1
  
  return(simmat)
}





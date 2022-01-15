#TSNE plot for new drug

# newDrug <- data.frame()
# newDrug.similarity <- as.data.frame(drug.similarity[,100])
# newDrug.similarity <- cbind(rownames(NewDrug.similarity), NewDrug.similarity)
# colnames(newDrug.similarity) <- c("Broad_id","YourDrug")

DepLink_TSNEplot <- function(newDrug.similarity,drugTable,
                             drug.similarity,drug.tsne, target.drug.id = NULL){
  require(plotly)
  require(gplots)
  require(RColorBrewer)

  
  
  # if(is.null(newDrug.similarity)){
  #   source("function/DepLink_Drug_Distance.R")
  #   NewDrug.similarity <- fp.sim.taminoto(fps.bits.idx1 =  ,
  #                              compared.idx2 = pubchem_maccs_bits_idx ,
  #                              nDrug.input = 2,
  #                              drugName.input = drug$broad_id[1:2],
  #                              combineMatrix = drug.sim.matrix )
  # }
  # 
 
    NewDrug.similarity <- NewDrug.similarity#[order(NewDrug.similarity[,"YourDrug"],decreasing = T),]
    similarDrug <- drug.tsne[which(drug.tsne$broad_id %in% NewDrug.similarity$broad_id[1:3]), ]
    
    #Calculate middle point among three points
    NewDrug.tsne.related <- data.frame(X1 = mean(similarDrug$X1),
                                       X2 = mean(similarDrug$X2),
                                       X3 = mean(similarDrug$X3),
                                       drug_name = "YourDrug")
    
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
      name = "YourDrug",
      data = NewDrug.tsne.related,
      x = ~X1,
      y = ~X2,
      z = ~X3,
      text = ~drug_name,
      hoverinfo = "text",
      #color = "red",
      #color= ~group,
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
      #color = "red",
      #color= ~group,
      mode = "markers",
      type = "scatter3d",
      marker = list(size = 6, color = "green"))
    
    if (!is.null(target.drug.id)) {
      fig <- fig %>% add_trace(
        name = "Selected Drug",
        data = drug.tsne[which(drug.tsne$broad_id == target.drug.id),1:4],
        x = ~X1,
        y = ~X2,
        z = ~X3,
        text = ~drug_name,
        hoverinfo = "text",
        #color = "red",
        #color= ~group,
        mode = "markers",
        type = "scatter3d",
        marker = list(size = 7, color = "yellow",symbol = "star-open",
                      line = list(width = 0.5, color = "black")))
    }
    
    return(fig) 
  
}

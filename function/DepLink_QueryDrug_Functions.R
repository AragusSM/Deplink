# Functions for DepLink Drug Query function.
# All plots and table generate functions are included in this document.
#    Project".: DepLink
#    Contributed by: Li-Ru Wang, Tapsya Nayak, Yu-Chiao Chiu, and Yidong Chen
#    Started from: October 2021
#    Affiliation: Greehey Children's Cancer Research Institute, 
#                 The University of Texas Health San Antonio
#


################################################################################
### DELETE IT Take it from QueryGene
DepLink_PanCan2CcleID <- function( cancer_type , sampleInfo = sample_info_mc){
  
  ccle_indx = which( sampleInfo$primary_disease_fixed %in% cancer_type )
  ccle_id = sampleInfo$DepMap_ID[ccle_indx]
  
  return(ccle_id)
  
}


################### Histogram of Dependency Score ##############################
# The function takes following inputs: a drug dataset in drg_data from PRISM or GDSC, name of
# the name of data in dataset_name (= "PRISM" or "GDSC") and type of tumor in pancan_type.
# Both PRISM and Sanger GDSC datasets, have cell ids in row names and Broad IDs as column.

DepLink_Density_Plot_Drug <- function(drug.data.p,  drug.data.g, drug_check_table,
                                      drug_check_table_idx,
                                      cancer_type = "PanCan", sampleInfo = NULL) {
  
  drug_table = drug_check_table$`Data Source`[drug_check_table_idx]
  #Drug data source
  if (drug_table == "PRISM") {
    #PRISM
    if (cancer_type != "PanCan") {
      ccle_id = DepLink_PanCan2CcleID(cancer_type = cancer_type, sampleInfo = sampleInfo)
      drug.data.p = as.data.frame(drug.data.p[rownames(drug.data.p) %in% ccle_id, ])
      
    }
    
    drug.data = drug.data.p[,drug_check_table$idx[drug_check_table_idx]]
    #rownames(drug.data) = rownames(drug.data.p)
    #colnames(drug.data) = colnames(drug.data.p)[drug_check_table$idx[drug_check_table_idx]]
    drug.name =  drug_check_table$`Drug Name`[drug_check_table_idx]
    title = paste("Log fold change (PRISM) of ",  drug.name)
    color_n = "skyblue"
    fillcolor = "rgba(135, 206, 235, 0.5)"
  }else{
    #GDSC
    if (cancer_type != "PanCan") {
      ccle_id = DepLink_PanCan2CcleID(cancer_type = cancer_type, sampleInfo = sampleInfo)
      drug.data.g = as.data.frame(drug.data.g[rownames(drug.data.g) %in% ccle_id, ])
      
    }
    drug.data = drug.data.g[,drug_check_table$idx[drug_check_table_idx]]
    #rownames(drug.data) = rownames(drug.data.g)
    #colnames(drug.data) = colnames(drug.data.g)[drug_check_table$idx[drug_check_table_idx]]
    drug.name =  drug_check_table$`Drug Name`[drug_check_table_idx]
    title = paste("Log10(IC50) of",  drug.name,"(GDSC)")
    
    color_n = "peachpuff"
    fillcolor = "rgba(255, 218, 185, 0.5)"
  }
  
  
  
  
  
  # single drug query
  t = c(" ")
  fig = plot_ly(alpha = 0.3)
  
  fig <-
    plot_ly(
      x = round(density(drug.data[!is.na(drug.data)])$x, 3),
      y = round(density(drug.data[!is.na(drug.data)])$y, 3),
      type = "scatter",
      mode = "lines",
      fill = "tozeroy",
      fillcolor = fillcolor,
      line =  list(color = color_n),
      name = paste0(drug.name, " (n=",length(drug.data[!is.na(drug.data)]),")")
    )
  
  fig <-
    fig %>% layout(
      title = t,
      font = list(size = 10),
      yaxis = list(title = "Density of cell lines"),
      xaxis = list(
        title = title,
        zeroline = FALSE,
        showgrid = FALSE
      ),
      legend = list(x = 100, y = 0.9,
                    title = list(text = paste0("<b> ",cancer_type," </b>"))) ,
      showlegend = TRUE
    )
  
  return(fig)
}



###################### Box Plot of Dependency Scores ###########################
# The function takes following inputs: a drug dataset in drg_data from PRISM or GDSC, 
# the manifest sample_info file name, a index to drug from dataset in brd_idx, the 
# name of data in dataset_name (= "PRISM" or "GDSC") and type of tumor in pancan_type.
DepLink_BoxPlot_PanCan_drug <- function( drug.data.p.t, drug.data.g.t, drug_check_table,
                                         drug_check_table_idx, sampleInfo) {
  
  drug.name <- drug_check_table$`Drug Name`[drug_check_table_idx]
  
  drug_table = drug_check_table$`Data Source`[drug_check_table_idx]
  
  if (drug_table == "PRISM") {
    
    # Query vaue of target drug from prism
    drug.idx = drug_check_table$idx[drug_check_table_idx]
    box_data = gather(drug.data.p.t[drug.idx,], ACHID, scores) #Including NA value
    indx =  match(box_data$ACHID, sampleInfo$DepMap_ID)
    box_data["tumor"] = sampleInfo$primary_disease_fixed[indx]
    box_data["data"] = rep("PRISM",nrow(box_data))
    color_n = "skyblue"
    title = paste("Log fold change (PRISM) of",  drug.name)
    
  }else{
    # Query value of target drug from gdsc
    drug.idx = drug_check_table$idx[drug_check_table_idx]
    box_data = gather(drug.data.g.t[drug.idx,], ACHID, scores) #Including NA value
    indx =  match(box_data$ACHID, sampleInfo$DepMap_ID)
    box_data["tumor"] = sampleInfo$primary_disease_fixed[indx]
    box_data["data"] = rep("GDSC",nrow(box_data))
    color_n ="peachpuff"
    title = paste("Log10(IC50) of",  drug.name,"(GDSC)")
  }
  
  
  
  fig = plot_ly(box_data, x=~tumor, y=~scores, color = ~data, colors = color_n,type="box")
  fig = fig %>% layout(xaxis=list(titlefont=list(size=10)),
                       yaxis=list(title =  title, 
                                  titlefont=list(size=10)),
                       showlegend = FALSE)
  return(fig)
}



###################### Scatter Plot Drug vs. Drug ##############################
# The function takes two input: a drug dataset in drg_data from PRISM or GDSC, and
# two drug broad ids in brd_idx1 and brd_idx2, and dataset name (="PRISM","GDSC")
DepLink_Drug_ScatterPlot <- function( drug.data ,drug_1, drug_2 ,drugname1,drugname2,drug.source =c("PRISM","GDSC"),
                                      cancer_type,sampleInfo) {
  
  
  # find ccle ids corresponding to pancan type
  if ( cancer_type != "PanCan") { 
    ccle_id = DepLink_PanCan2CcleID(cancer_type = cancer_type, sampleInfo = sampleInfo)
    
    drug.data = drug.data[rownames(drug.data) %in% ccle_id, ]
    
  }
  
  if (drug.source == "PRISM") {
    if (nrow(drug.data) >= 3) {
      # prepare data for PRISM scatter plot of two genes.
      combined = as.data.frame(cbind(drug.data[,drug_1], drug.data[,drug_2]))
      combined = na.omit(combined)
      colnames(combined) = c("d1", "d2")
      combined["dt"]  = rep(paste0(drug.source," (n=",nrow(combined),")"), nrow(combined))
      
      fit = lm(d2 ~ d1, data = combined )
      
      t = paste(" ")
      fig <- plot_ly(data = combined, x = ~d1, y = ~d2, type="scatter", color= ~dt, mode = 'markers',
                     colors = "skyblue", alpha=0.5, showlegend=F)
      fig <- fig %>% add_trace(x = ~d1, y = fitted(fit), mode="lines", alpha=30, showlegend=T)
      fig <- fig %>% layout(xaxis = list( title = paste(drugname1, "(Log fold change)"), zeroline=FALSE),
                            yaxis = list( title = paste(drugname2, "(Log fold change)"), zeroline=FALSE ),
                            showlegend = TRUE)
      
      
    }
  }else{
    if (nrow(drug.data) >= 3) {
      # prepare data for GDSC scatter plot of two genes.
      combined = as.data.frame(cbind(drug.data[,drug_1], drug.data[,drug_2]))
      combined = na.omit(combined)
      colnames(combined) = c("d1", "d2")
      combined["dt"]  = rep(paste0(drug.source," (n=",nrow(combined),")"), nrow(combined))
      
      fit = lm(d2 ~ d1, data = combined )
      
      t = paste(" ")
      fig <- plot_ly(data = combined, x = ~d1, y = ~d2, type="scatter", color= ~dt, mode = 'markers',
                     colors = "peachpuff", alpha=0.5, showlegend=F)
      fig <- fig %>% add_trace(x = ~d1, y = fitted(fit), mode="lines", alpha=30, showlegend=T)
      fig <- fig %>% layout(xaxis = list( title = paste(drugname1, "(Log10(IC50))"), zeroline=FALSE),
                            yaxis = list( title = paste(drugname2, "(Log10(IC50))"), zeroline=FALSE ),
                            showlegend = TRUE)
      
    }
  }
  
  return(fig)
}


################# Correlation Coefficients Drug vs Drugs #######################
# The function takes two input: a drug dataset in drg_data from PRISM or GDSC, 
# a drug broad id in brd_idx, dataset name in dataset_name(="PRISM","GDSC"), cancer
# type in pancan_type and met i.e. method to be used for computing correlation coefficient.
DepLink_Codependency_drug <- function( drug.data,
                                       drug_check_table,
                                       drug_check_table_idx,
                                       cancer_type="PanCan",
                                       met = "pearson",
                                       sampleInfo  ){
  
  #Drug data source
  
  if (cancer_type != "PanCan") {
    ccle_id = DepLink_PanCan2CcleID(cancer_type = cancer_type, sampleInfo = sampleInfo)
    drug.data = as.data.frame(drug.data[rownames(drug.data) %in% ccle_id, ])
    
  }
  
  if ( nrow(drug.data) < 3 ){
    
    correlation_matrix = NA
  }else{
    d = as.data.frame( drug.data[,drug_check_table$idx[drug_check_table_idx]])
    d1 = as.data.frame( drug.data[,-drug_check_table$idx[drug_check_table_idx]])
    
    correlation_matrix <- cor(d,d1, use="pairwise.complete.obs", method = met)
    
    
  }
  
  correlation_matrix_update <- correlation_matrix
  
  for (i in 1:length(correlation_matrix)) {
    if (correlation_matrix[,i] %in% c(1, -1)) {
      #Exculde nrow < 3 after na.omit 
      df <- cbind(d, d1[,i])
      df.noNA <- na.omit(df)
      
      if (nrow(df.noNA) < 3) {
        correlation_matrix_update[,i] <- NA
      }else{
        correlation_matrix_update[,i] <- correlation_matrix[,i] 
      }
      
    }else{
      
      correlation_matrix_update[,i] <- correlation_matrix[,i] 
    }
    
  }
  correlation_matrix_update = as.numeric(round(correlation_matrix_update,3))
  
  return(correlation_matrix_update)
}


################## Histogram of Correlation Coefficient ########################
# The function takes following inputs: drug correlation cofficients in rho_1
# and drug broad id index in brd_id, and drug dataset name in data_name(="PRISM","GDSC")
DepLink_Density_Plot_CorCoeff_Drug <- function(rho, drug_name, drug.type = c("PRISM","GDSC")){
  
  rho = rho[!is.na(rho)]
  if (drug.type == "PRISM") {
    
    fig = plot_ly(alpha=0.3)
    fig <- fig %>% add_trace(x=round(density(rho)$x, 3), y=round(density(rho)$y, 3), type="scatter",mode="lines",
                             fill="tozeroy", line=list(color="skyblue"),
                             fillcolor = "rgba(135, 206, 235, 0.5)",
                             name=paste0("PRISM"))
  }else{
    fig = plot_ly(alpha=0.3)
    fig <- fig %>% add_trace(x=round(density(rho)$x, 3), y=round(density(rho)$y, 3), type="scatter",mode="lines",
                             fill="tozeroy", line=list(color="peachpuff"),
                             fillcolor = "rgba(255, 218, 185, 0.5)",
                             name=paste0("GDSC"))
  }
  
  fig <- fig %>% layout(title=t, font=list(size=10), yaxis=list(title = "Density of drugs"),
                        xaxis=list(title= paste0("Correlation coefficients with ", drug_name),zeroline=FALSE, showgrid=FALSE), legend = list(x = 100, y = 0.9),
                        showlegend = FALSE)
  
  return(fig)
}

################## Correlation Coefficients Drug vs Gene #######################
# The function takes the following inputs: one drug data such as PRISM/GDSC in drug_data,
# one gene data from DepMap/Sanger in gene_data, drug index in brd_idx.
DepLink_D2G_Codependency<- function( drug.data,
                                     gene.data,
                                     drug_id,
                                      cancer_type,
                                     met = "pearson",
                                     sampleInfo  ){
  
  
  #Drug data source
  
  if (cancer_type != "PanCan") {
    ccle_id = DepLink_PanCan2CcleID(cancer_type = cancer_type, sampleInfo = sampleInfo)
    drug.data = as.data.frame(drug.data[rownames(drug.data) %in% ccle_id, ])
    gene.data = as.data.frame(gene.data[rownames(gene.data) %in% ccle_id, ])
    
  }
  
  if ( nrow(drug.data) < 3 ){
    
    correlation_matrix_update = NA
  }else{
    drug.data = as.data.frame( drug.data[,drug_id])
    
   correlation_matrix <- cor(drug.data,gene.data, use="pairwise.complete.obs", method = met)
    
  correlation_matrix_update <- correlation_matrix
  
  for (i in 1:length(correlation_matrix)) {
    if (correlation_matrix[,i] %in% c(1, -1)) {
      #Exculde nrow < 3 after na.omit 
      df <- cbind(d, d1[,i])
      df.noNA <- na.omit(df)
      
      if (nrow(df.noNA) < 3) {
        correlation_matrix_update[,i] <- NA
      }else{
        correlation_matrix_update[,i] <- correlation_matrix[,i] 
      }
      
    }else{
      
      correlation_matrix_update[,i] <- correlation_matrix[,i] 
    }
    
  }
  correlation_matrix_update = as.numeric(round(correlation_matrix_update,3))
}
  return(correlation_matrix_update)
}


############## Histogram of Correlation Coefficients Drug vs Gene ##############
################### Two gene datasets - PRISM and Sanger #######################
# The function takes following inputs: drug data correlation with depMap and sanger dataset
# in depMap_rho and sanger_rho, drug BROAD ID in brd_id, and drug dataset name in data_name
# (="PRISM" or "GDSC")
DepLink_D2G_Density_Plot <- function( rho1,rho1_n, rho2, rho2_n, drug_name, drug_type){
  
  rho1 = rho1[!is.na(rho1)]
  rho2 = rho2[!is.na(rho2)]
  
  fig = plot_ly(alpha=0.3)
  if(length(rho1) >0){
    fig <- fig %>% add_trace(x=round(density(rho1)$x,3), y=round(density(rho1)$y,3), type="scatter", yaxis="y2", mode="lines",fill="tozeroy", line=list(color="skyblue"),
                             legendgroup="DepMap", name=paste0("Broad DepMap"),
                             showlegend = TRUE)
  }
  
  if (length(rho2) >0) {
    fig <- fig %>% add_trace(x=round(density(rho2)$x,3), y=round(density(rho2)$y,3), type="scatter", yaxis="y2", mode="lines",fill="tozeroy",
                             line=list(color="peachpuff"),fillcolor = "rgba(255, 218, 185, 0.5)",
                             legendgroup="Sanger", name=paste0("Sanger DepMap"),
                             showlegend = TRUE)
  }
  
  
  if (sum(length(rho1) > 0 & length(rho2) > 0 )== 1) {
    
    fig <- fig %>% layout(barmode="overlay", title=t, font=list(size=10), yaxis2=list(title = "Density of genes", overlaying="y"),
                          xaxis=list(title=paste0("Correlation coefficients with ", drug_name),zeroline=FALSE, showgrid=FALSE),
                          legend = list(x = 100, y = 0.9,
                                        title = list(text = paste0("<b> ",drug_type," vs. </b>"))
                                        ))
    
    
  }else{
    
    fig <- fig %>% layout( title=t, font=list(size=10), yaxis2=list(title = "Density of genes"),
                           xaxis=list(title=paste0("Correlation coefficients with ", drug_name),zeroline=FALSE, showgrid=FALSE),
                           legend = list(x = 100, y = 0.9,
                                         title = list(text = paste0("<b> ",drug_type," vs. </b>"))
                                         ), showlegend = TRUE)
  }
  return(fig)
}

######################## Scatter Plot Drug vs. Gene ############################
################### Two gene datasets - PRISM and Sanger #######################
# The function takes the following inputs: one drug data such as PRISM/GDSC in drug_data,
# DepMap and Sanger dataset in gene_data1 and gene_data2, drug index in brd_idx, gene name 
# in brd_idx and drug data name in data_name(="PRISM" and "GDSC")
DepLink_D2G_ScatterPlot <- function(drug_data, gene_data1,
                                    gene_data2,
                                    gene_name, broad_id, drug_name,
                                    sampleInfo, cancer_type,drug_source){
  
  drug_rownames = rownames(drug_data)
  gene1_rownames = rownames(gene_data1)
  gene2_rownames = rownames(gene_data2)
  
  
  interNames1 = intersect(drug_rownames, gene1_rownames)
  gene_data1 = gene_data1[interNames1,]
  drug_data1 = drug_data[interNames1,]
  
  interNames2 = intersect(drug_rownames, gene2_rownames)
  gene_data2 = gene_data2[interNames2,]
  drug_data2 = drug_data[interNames2,]
  
  
  if (cancer_type!="PanCan"){
    ccle_id = DepLink_PanCan2CcleID(cancer_type = cancer_type, sampleInfo = sampleInfo)
    gene_data1 = gene_data1[rownames(gene_data1) %in% ccle_id,]
    gene_data2 = gene_data2[rownames(gene_data2) %in% ccle_id,]
    drug_data1 = drug_data1[rownames(drug_data1) %in% ccle_id,]
    drug_data2 = drug_data2[rownames(drug_data2) %in% ccle_id,]
  }
  
  if(nrow(drug_data1) >= 3 & nrow(gene_data1) >=3){
    #gene_data1
    combined = as.data.frame(cbind(gene_data1[,gene_name], drug_data1[,broad_id]))
    colnames(combined) = c("g1","d1")
    combined = na.omit(combined)
    combined["dt"] = rep(paste0("Broad DepMap (n=",nrow(combined),")"), nrow(combined))
  }else{
    combined = data.frame(matrix(data = NA, nrow = 1,ncol = 3))
  }
  
  if (nrow(drug_data2) >=3 & nrow(gene_data2) >= 3) {
    #gene_data2
    combined2 = as.data.frame(cbind(gene_data2[,gene_name], drug_data2[,broad_id]))
    colnames(combined2) = c("g1","d1")
    combined2 = na.omit(combined2)
    combined2["dt"] = rep(paste0("Sanger DepMap (n=",nrow(combined2),")"), nrow(combined2))
  }else{
    combined2 = data.frame(matrix(data = NA, nrow = 1,ncol = 3))
  }
  
  if (nrow(combined) >= 3 & nrow(combined2) >=3) {
    
    merged = rbind(combined, combined2)
    fit = lm(g1 ~ d1*dt, data = merged )
    fig <- plot_ly(data = merged, x = ~d1, y = ~g1, type="scatter", color=~dt, alpha=0.5, mode="markers", 
                   colors=c("skyblue","peachpuff"), showlegend=F)
    fig <- fig %>% add_lines(x = ~d1, y = fitted(fit), mode="lines", alpha=30, showlegend=T)
  }
  
  if (nrow(combined) >= 3 & nrow(combined2) < 3) {
    
    merged = combined
    
    fit = lm(g1~d1, data = merged )
    fig <- plot_ly(data = merged, x = ~d1, y = ~g1, type="scatter", color=~dt, alpha=0.5, mode="markers", 
                   colors="skyblue", showlegend=F)
    fig <- fig %>% add_lines(x = ~d1, y = fitted(fit), mode="lines", alpha=30, showlegend=T)
  }
  
  if (nrow(combined) < 3 & nrow(combined2) >= 3) {
    merged = combined2
    fit = lm(g1 ~ d1, data = merged )
    fig <- plot_ly(data = merged, x = ~d1, y = ~g1, type="scatter", color=~dt, alpha=0.5, mode="markers", 
                   colors="peachpuff", showlegend=F)
    fig <- fig %>% add_lines(x = ~d1, y = fitted(fit), mode="lines", alpha=30, showlegend=T)
  }
  
  if (drug_source == "PRISM") {
    x.title = paste(drug_name,"(Log fold change)")
  }else{
    x.title = paste(drug_name,"(Log10 IC50)")
  }
  
  
  fig <- fig %>% layout(legend=list(title=list(text=paste0("<b> ",drug_source," vs. </b>")),font = list(size = 10)),
                        xaxis = list( title = x.title, zeroline=FALSE), 
                        yaxis = list( title = paste(gene_name,"(Chronos gene effect score)"), zeroline=FALSE), showlegend=T)
  return(fig)
  
}



# VisNetwork for Drug-Drug association. 
# input is the ddTable (see first table for details), drug (e.g. BAM7), and ndx (visual table entries in ddTable)
DepLink_D2D_ReturnNet <- function( ddTable, drug, ndx ) {
  #  ddTable = ddTable()
  
  nodes <- data.frame( id = 1:11, label = c(drug, ddTable$`Drug Name`[ndx]), shape = "diamond", color = c( "skyblue", rep( "orange", 10)) )
  edges <- data.frame( id = 1:10, from = rep(1, 10), to = 2:11, width = abs(10*ddTable$Correlation[ndx]), color = "skyblue" )
  edges <- mutate(edges, color = ifelse( ddTable$Correlation[ndx] > 0, "red", "blue"))
  
  ledges <- data.frame(color = c("red", "blue"),
                       label = c("r>0", "r<0"),
                       width = 4,
                       arrows.scaleFactor = 0 )
  
  network <- visNetwork(nodes, edges, label = nodes$label) %>%
    visNodes(size = 20 ) %>%
    visLegend(width = 0.1, addEdges = ledges, useGroups = FALSE) %>%
    visInteraction(multiselect = TRUE,navigationButtons = TRUE) %>%
    #visOptions(selectedBy = "label",  highlightNearest = TRUE,  nodesIdSelection = TRUE) %>%
    visPhysics(stabilization = FALSE)
  
  return(network)
  
}


# VisNetwork for Drug-Gene association. 
# input is the dgTable (see panel 1I for table detail), drug (e.g. BAM7), and ndx (visual table entries in dgTable)

DepLink_D2G_ReturnNet <- function( dgtable, drug, ndx ) {
  
  # create nodes and edges. Different from ggTable, we set drug with different shapes.
  nodes <- data.frame( id = 1:11, label = c(drug, dgtable$`Gene Name`[ndx]), shape = c("diamond", rep("circle", 10)), color = c( "skyblue", rep( "orange", 10)) )
  edges <- data.frame( id = 1:10, from = rep(1, 10), to = 2:11, width = abs(10*dgtable$`Corr in Broad DepMap`[ndx]), color = "skyblue" )
  edges <- mutate(edges, color = ifelse( dgtable$`Corr in Broad DepMap`[ndx] > 0, "red", "blue"))
  
  ledges <- data.frame(color = c("red", "blue"),
                       label = c("r>0", "r<0"),
                       width = 4,
                       arrows.scaleFactor = 0 )
  
  network <- visNetwork(nodes, edges, label = nodes$label) %>%
    visNodes(size = 20 ) %>%
    visLegend(width = 0.1, addEdges = ledges, useGroups = FALSE) %>%
    visInteraction(multiselect = TRUE, navigationButtons = TRUE) %>%
    #visOptions(selectedBy = "label",  highlightNearest = TRUE,  nodesIdSelection = TRUE) %>%
    visPhysics(stabilization = FALSE)
  
  return(network)
  
}
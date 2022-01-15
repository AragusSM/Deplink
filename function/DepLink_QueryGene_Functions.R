# Functions for DepLink Gene Query function.
# All plots and table generate functions are included in this document.
#    Project".: DepLink
#    Contributed by: Li-Ru Wang, Tapsya Nayak, Yu-Chiao Chiu, and Yidong Chen
#    Started from: October 2021
#    Affiliation: Greehey Children's Cancer Research Institute, 
#                 The University of Texas Health San Antonio
#


#######################################
# This function takes one input pancan type and returns corresponding CCLE ID.
# This function is mainly called by other functions.
# Change function name and revise function (add sampleInfo variable, add ifelse() and revise plotly code)
DepLink_PanCan2CcleID <- function( cancer_type , sampleInfo = sample_info_mc){
  
  ccle_indx = which( sampleInfo$primary_disease_fixed %in% cancer_type )
  ccle_id = sampleInfo$DepMap_ID[ccle_indx]
  
  return(ccle_id)

}


#######################################
# The function takes following inputs: DepMap (d_data) and Sanger (s_data) datasets,
# a genes name in gene_name, such as "A1CF" and a pancan type such "Breast" (default is  "ALL")
# Both DepMap and Sanger datasets have unique column names (genes), as well as row names (cell-lines)
DepLink_Density_Plot_DepScores <- function(d_data, s_data, gene_name, cancer_type = "PanCan", sampleInfo = NULL) {
    if (length(gene_name) == 1) {
      # find ccle ids corresponding to pancan type
      if (cancer_type != "PanCan") {
        ccle_id = DepLink_PanCan2CcleID(cancer_type = cancer_type, sampleInfo = sampleInfo)
        
        d_data = d_data[, colnames(d_data) %in% ccle_id]
        s_data = s_data[, colnames(s_data) %in% ccle_id]
      }
      
      # Single gene query
      
      t = c(" ")
      fig = plot_ly(alpha = 0.3)
      
      
      if (is.data.frame(d_data) &
          ncol(d_data) >= 3 & is.data.frame(s_data)) {
        if (ncol(s_data) >= 3) {
          d_data.num = as.numeric(d_data[gene_name, ])
          
          fig <-
            fig %>% add_trace(
              x = round(density(d_data.num[!is.na(d_data.num)])$x, 3),
              y = round(density(d_data.num[!is.na(d_data.num)])$y, 3),
              type = "scatter",
              mode = "lines",
              fill = "tozeroy",
              line = list(color = "skyblue"),
              yaxis = "y2",
              name = paste0("Broad DepMap (n=", ncol(d_data), ")")
            )
          
          s_data.num = as.numeric(s_data[gene_name, ])
          fig <-
            fig %>% add_trace(
              x = round(density(s_data.num[!is.na(s_data.num)])$x, 3),
              y = round(density(s_data.num[!is.na(s_data.num)])$y, 3),
              type = "scatter",
              yaxis = "y2",
              mode = "lines",
              fill = "tozeroy",
              line = list(color = "peachpuff"),
              #legendgroup="Sanger",
              name = paste0("Sanger DepMap (n=", ncol(s_data), ")")
            )
          
          fig <-
            fig %>% layout(
              barmode = "overlay",
              title = t,
              font = list(size = 10),
              yaxis2 = list(title = "Density of cell lines", overlaying =
                              "y"),
              xaxis = list(
                title = paste0("Chronos gene effect score of ", gene_name),
                zeroline = FALSE,
                showgrid = FALSE
              ),
              legend = list(x = 100, y = 0.9,title = list(text = paste0("<b> ",cancer_type," </b>")))
            )
        } else{
          d_data.num = as.numeric(d_data[gene_name, ])
          
          fig <-
            fig %>% add_trace(
              x = round(density(d_data.num[!is.na(d_data.num)])$x, 3),
              y = round(density(d_data.num[!is.na(d_data.num)])$y, 3),
              type = "scatter",
              mode = "lines",
              fill = "tozeroy",
              line = list(color = "skyblue"),
              name = paste0("Broad DepMap (n=", ncol(d_data), ")")
            )
          
          fig <-
            fig %>% layout(
              title = t,
              font = list(size = 10),
              yaxis = list(title = "Density of cell lines"),
              xaxis = list(
                title = paste0("Chronos gene effect score of ", gene_name),
                zeroline = FALSE,
                showgrid = FALSE
              ),
              legend = list(x = 100, y = 0.9) ,
              showlegend = TRUE
            )
        }
        
      } else{
        d_data.num = as.numeric(d_data[gene_name, ])
        
        fig <-
          fig %>% add_trace(
            x = round(density(d_data.num[!is.na(d_data.num)])$x, 3),
            y = round(density(d_data.num[!is.na(d_data.num)])$y, 3),
            type = "scatter",
            mode = "lines",
            fill = "tozeroy",
            line = list(color = "skyblue"),
            name = paste0("Broad DepMap (n=", ncol(d_data), ")")
          )
        
        fig <-
          fig %>% layout(
            title = t,
            font = list(size = 10),
            yaxis = list(title = "Density of cell lines"),
            xaxis = list(
              title = paste0("Chronos gene effect score of ", gene_name),
              zeroline = FALSE,
              showgrid = FALSE
            ),
            legend = list(x = 100, y = 0.9,title = list(text = paste0("<b> ",cancer_type," </b>"))) ,
            showlegend = TRUE
          )
      }
      
    } else {
      # Multiple gene queries. Only up to 10 genes.
      
      if (length(gene_name) <= 10) {
        fig = plot_ly(alpha = 0.3)
        #gene_name=as.character(gene_name)
        d_data = d_data[gene_name, ]
        color_seq = RColorBrewer::brewer.pal(n = 10, name = "Set3")
        
        # Display depMap genes density plots
        for (i in 1:length(gene_name)) {
        fig <-  fig %>% add_trace(
                                  x = round(density(as.numeric(d_data[i, ]))$x, 3),
                                  y = round(density(as.numeric(d_data[i, ]))$y, 3),
                                  type = "scatter",
                                  yaxis = "y2",
                                  mode = "lines",
                                  fill = "tozeroy",
                                  colors = color_seq[i],
                                  legendgroup = "Broad DepMap",
                                  name = gene_name[i]
                                )
          fig <- fig %>%  layout(
                                  barmode = "overlay",
                                  title = t,
                                  font = list(size = 10),
                                  yaxis2 = list(title = "Density of cell lines", overlaying =  "y"),
                                  xaxis = list(
                                    title = paste0("Chronos gene effect score of ", gene_name),
                                    zeroline = FALSE,
                                    showgrid = FALSE
                                  ),
                                  legend = list(
                                    x = 100,
                                    y = 0.9,
                                    size = 10,
                                    title = list(text = '<b> Broad DepMap </b>')
                                  ),
                                  showlegend = T
                                )
        }
        
      }
    }
    return(fig)
  }


#######################################
## Updated 6 December 2021
# The function takes following input: DepMap (d_data) and Sanger (s_data) datasets, and
# a gene name such as "A1CF".
# Both DepMap and Sanger datasets have unique row names (genes), as well as column names (cell-lines)
# Remove Cancer type variable and change function name
DepLink_BoxPlot_PanCan <- function( d_data, s_data, gene_name , sampleInfo) {
  
  # Query dependency scores of gene name from d_data
  dt_data = gather(d_data[gene_name,], ACHID, scores)
  indx =  match(dt_data$ACHID, sampleInfo$DepMap_ID)
  dt_data["tumor"] = sampleInfo$primary_disease_fixed[indx]
  dt_data["data"] = rep("Broad DepMap",nrow(dt_data))
  
  # Query dependency scores of gene name from s_data
  st_data = gather(s_data[gene_name,], ACHID, scores)
  indx =  match(st_data$ACHID, sampleInfo$DepMap_ID)
  st_data["tumor"] = sampleInfo$primary_disease_fixed[indx]
  st_data["data"] = rep("Sanger DepMap",nrow(st_data))
  
  # Merge the data frames for plotting
  merged = rbind(dt_data,st_data)
  
  fig = plot_ly(merged, x=~tumor, y=~scores, color=~data, type="box", colors = c("skyblue","peachpuff"), alpha = 0.5)
  fig = fig %>% layout(boxmode="group", xaxis=list(titlefont=list(size=10)),
                       yaxis=list(title = paste0("Chronos gene effect score of ", gene_name), titlefont=list(size=10)))
  return(fig)
}


#######################################
## Updated 6 December 2021
# The function takes two input: DepMap (d_data) and Sanger (s_data) datasets, and
# two genes names in x and y, such as "A1CF", "A4GALT".
# Both DepMap and Sanger datasets have unique column names (genes), as well as row names (cell-lines)
DepLink_Gene_ScatterPlot <- function( d_data,depMap_data, s_data, gene_1, gene_2 ,
                                      cancer_type,sampleInfo) {
  
  #library(data.table)
  nsValue <- as.data.frame(which(is.na(depMap_data), arr.ind = TRUE))
  d_data.clean <- d_data[which(!rownames(d_data) %in% colnames(depMap_data)[unique(nsValue$col)]), ]
  
  # find ccle ids corresponding to pancan type
  if ( cancer_type != "PanCan") { 
    ccle_id = DepLink_PanCan2CcleID(cancer_type = cancer_type, sampleInfo = sampleInfo)
    
    d_data.clean = d_data.clean[rownames(d_data.clean) %in% ccle_id, ]
    s_data = s_data[rownames(s_data) %in% ccle_id, ]
  }
  
  if (nrow(s_data) < 3) {
    # prepare data for DepMap scatter plot of two genes.
    combined = as.data.frame(cbind(d_data.clean[,gene_1], d_data.clean[,gene_2]))
    colnames(combined) = c("g1", "g2")
    combined["dt"]  = rep(paste0("Broad DepMap (n=",nrow(combined),")"), nrow(combined))
    
    fit = lm(g2 ~ g1, data = combined )
    
    t = paste(" ")
    fig <- plot_ly(data = combined, x = ~g1, y = ~g2, type="scatter", color= ~dt, mode = 'markers',
                   colors = "skyblue", alpha=0.5, showlegend=F)
    fig <- fig %>% add_trace(x = ~g1, y = fitted(fit), mode="lines", alpha=30, showlegend=T)
    fig <- fig %>% layout(xaxis = list( title = paste(gene_1, "(Chronos gene effect score)"), zeroline=FALSE),
                          yaxis = list( title = paste(gene_2, "(Chronos gene effect score)"), zeroline=FALSE ),
                          showlegend = TRUE)
    
    
  }else{
  # prepare data for DepMap scatter plot of two genes.
  combined = as.data.frame(cbind(d_data.clean[,gene_1], d_data.clean[,gene_2]))
  colnames(combined) = c("g1", "g2")
  combined["dt"]  = replicate(nrow(combined),paste0("Broad DepMap (n=",nrow(combined),")"))

  # prepare data for Sanger scatter plot of two genes.
  combined_s = as.data.frame(cbind(s_data[,gene_1], s_data[,gene_2]))
  colnames(combined_s) = c("g1", "g2")
  combined_s["dt"] = rep(paste0("Sanger DepMap (n=",nrow(combined_s),")"), nrow(combined_s))

  # merge two data.frames  
  merged = rbind(combined, combined_s)
  fit = lm(g2 ~ g1*dt, data = merged )
  
  t = paste(" ")
  fig <- plot_ly(data = merged, x = ~g1, y = ~g2, type="scatter", color= ~dt, mode = 'markers',
                 colors = c("skyblue","peachpuff"), alpha=0.5, showlegend= F)
  fig <- fig %>% add_trace(x = ~g1, y = fitted(fit), mode="lines", alpha=30, showlegend=T)
  fig <- fig %>% layout(xaxis = list( title = paste(gene_1, "(Chronos gene effect score)"), zeroline=FALSE), 
                        yaxis = list( title = paste(gene_2, "(Chronos gene effect score)"), zeroline=FALSE ),
                        showlegend = TRUE
                        )
  }
  return(fig)
}


#######################################
# The function takes following inputs: DepMap (d_data), pancan type (pancan_type) such as "Breast" (default="ALL")
# and gene name (gene_name) such as "A1CF". Default type of correlation is pearson.
DepLink_Codependency <- function( data, gene_name, cancer_type="PanCan", met = "pearson",sampleInfo  ){

   if ( cancer_type!="PanCan" ){
     ccle_id = DepLink_PanCan2CcleID(cancer_type = cancer_type, sampleInfo = sampleInfo)
     data = data[rownames(data) %in% ccle_id, ]
   }
  if ( nrow(data) < 3 ){
  
    correlation_matrix = NA
  }else{
    d = as.data.frame(data[,gene_name])
    d1 = data[,!(names(data) %in% gene_name)]
    correlation_matrix <- cor(d,d1, use="pairwise.complete.obs", method = met)
    correlation_matrix <- as.numeric(round(correlation_matrix,3))
    
  }
  return(correlation_matrix)
}


# The function takes following inputs: two vectors of correlation coefficient in rho_1 and rho_2
# and gene name (gene_name) such as "A1CF".
DepLink_Density_Plot_CorCoeff <- function(rho_1,rho_2, gene_name){
  rho_1 = rho_1[!is.na(rho_1)]
  rho_2 = rho_2[!is.na(rho_2)]
  
  if(sum(is.na(rho_2)) > 0){
    fig = plot_ly(alpha=0.3)
    fig <- fig %>% add_trace(x=round(density(rho_1)$x, 3), y=round(density(rho_1)$y, 3), type="scatter",mode="lines",
                             fill="tozeroy", line=list(color="skyblue"),
                             #legendgroup="DepMap",
                             name=paste("Broad DepMap"))
    fig <- fig %>% layout(title=t, font=list(size=10), yaxis=list(title = "Density of genes"),
                          xaxis=list(title= paste0("Correlation coefficients with ", gene_name),zeroline=FALSE, showgrid=FALSE), legend = list(x = 100, y = 0.9),
                          showlegend = TRUE)
  }else{
  fig = plot_ly(alpha=0.3)
  fig <- fig %>% add_trace(x=round(density(rho_1)$x, 3), y=round(density(rho_1)$y, 3), type="scatter", 
                           yaxis="y2", mode="lines",fill="tozeroy", line=list(color="skyblue"),
                           legendgroup="DepMap", name=paste("Broad DepMap"))
  fig <- fig %>% add_trace(x=round(density(rho_2)$x, 3), y=round(density(rho_2)$y, 3), type="scatter", 
                           yaxis="y2", mode="lines",fill="tozeroy", line=list(color="peachpuff"),
                           line= list(dash = 'dash'), legendgroup="Sanger", name=paste("Sanger DepMap"))
  fig <- fig %>% layout(barmode="overlay", title=t, font=list(size=10), yaxis2=list(title = "Density of genes", overlaying="y"),
                        xaxis=list(title=paste0("Correlation coefficients with ", gene_name),zeroline=FALSE,
                                   showgrid=FALSE), legend = list(x = 100, y = 0.9),
                        showlegend = TRUE)
  }
  return(fig)
}


#######################################
# This function takes the following inputs dependency data in g_data, drug data in d_data,
# gene name in (gene_name) such as "A1CF", pan can type in (pancan_type) such as "Breast"
DepLink_G2D_Codependency <- function( gene_data, drug_data,
                                      gene_name, cancer_type,
                                      met = "pearson",sampleInfo ){
  if (cancer_type!="PanCan"){
    ccle_id = DepLink_PanCan2CcleID(cancer_type = cancer_type, sampleInfo = sampleInfo)
    gene_data = gene_data[rownames(gene_data) %in% ccle_id,]
    drug_data = drug_data[rownames(drug_data) %in% ccle_id,]
  }
  if (nrow(gene_data) < 3 || nrow(drug_data) < 3){
    correlation_matrix = NA
    
  } else{
    gene_data = as.data.frame(gene_data[,gene_name])
    correlation_matrix = cor(gene_data, drug_data, use="pairwise.complete.obs", method = met)
    
   # merged <- cbind(gene_data,drug_data)
    correlation_matrix_update <- correlation_matrix
    
    for (i in 1:length(correlation_matrix)) {
      if (correlation_matrix[,i] %in% c(1, -1)) {
        #Exculde nrow < 3 after na.omit 
        df <- cbind(gene_data, drug_data[,i])
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
}

# This function takes the following inputs: vectors of correlation coefficient of gene and drug.
# variable rho1: takes correlation coefficients between data depMap and PRISM
# variable rho2: takes correlation coefficients between data Sanger and PRISM
# variable rho3: takes correlation coefficients between data depMap and GDSC
# variable rho4: takes correlation coefficients between data Sanger and GDSC
# Variable gene_name: take gene names such as "A1CF", for legend display
DepLink_G2D_Density_Plot <- function( DepMap_PRISM, Sanger_PRISM,
                                      DepMap_GDSC, Sanger_GDSC,
                                      gene_name ){
  fig = plot_ly(alpha=0.3)
  rho1 = DepMap_PRISM[!is.na(DepMap_PRISM)]
  rho2 = Sanger_PRISM[!is.na(Sanger_PRISM)]
  rho3 = DepMap_GDSC[!is.na(DepMap_GDSC)]
  rho4 = Sanger_GDSC[!is.na(Sanger_GDSC)]
  
  if(length(rho1) >0){
    fig <- fig %>% add_trace(x=round(density(rho1)$x,3), y=round(density(rho1)$y,3), type="scatter", yaxis="y2", mode="lines",fill="tozeroy", line=list(color="skyblue"),
                             legendgroup="DepMap", name=paste0("Broad-PRISM"), #paste0("Broad-PRISM: ", gene_name, "(n=",length(rho1),")")
                             showlegend = TRUE)
  }
  
  if (length(rho2) >0) {
    fig <- fig %>% add_trace(x=round(density(rho2)$x,3), y=round(density(rho2)$y,3), type="scatter", yaxis="y2", mode="lines",fill="tozeroy", line=list(color="peachpuff"),
                             legendgroup="Sanger", name=paste0("Sanger-PRISM"), #paste0("Sanger-PRISM: ", gene_name,"(n=",length(rho2),")")
                             showlegend = TRUE)
  }
  if (length(rho3) >0) {
    fig <- fig %>% add_trace(x=round(density(rho3)$x,3), y=round(density(rho3)$y,3), type="scatter", yaxis="y2", mode="lines",fill="tozeroy", line=list(color="skyblue", dash="dash", width=3),
                             legendgroup="DepMap", name=paste0("Broad-GDSC"), #name=paste0("Broad-GDSC: ", gene_name, "(n=",length(rho3),")")
                             showlegend = TRUE)
  }
  if (length(rho4) > 0) {
    fig <- fig %>% add_trace(x=round(density(rho4)$x,3), y=round(density(rho4)$y,3), type="scatter", yaxis="y2", mode="lines",fill="tozeroy", line=list(color="peachpuff", dash="dash", width=3),
                             legendgroup="Sanger", name=paste0("Sanger-GDSC"), #paste0("Sanger-GDSC: ", gene_name, "(n=",length(rho4),")")
                             showlegend = TRUE)
  }
  
  if (sum(length(rho1) > 0 |length(rho2) > 0|length(rho3) > 0|length(rho4) > 0) == 1) {
    
    fig <- fig %>% layout( title=t, font=list(size=10), yaxis2=list(title = "Density of drugs"),
                          xaxis=list(title=paste0("Correlation coefficients with ", gene_name),zeroline=FALSE, showgrid=FALSE),
                          legend = list(x = 100, y = 0.9,font = list(size =10),
                                        title=list(text=paste0("<b> Gene-Drug Databases </b>"))), showlegend = TRUE
                          )

  }else{ 
  
  fig <- fig %>% layout(barmode="overlay", title=t, font=list(size=10), yaxis2=list(title = "Density of drugs", overlaying="y"),
                        xaxis=list(title= paste0("Correlation coefficients with ", gene_name),zeroline=FALSE, showgrid=FALSE),
                        legend = list(
                                      # x = 100,
                                      # y = 0.9,
                                      font = list(size =10),
                                      title=list(text="<b> Gene-Drug Databases </b>"))
                        , showlegend = TRUE
                        # layout(legend=list(title=list(text=paste0("<b> ",drug_source,": </b>"))),
                        #        xaxis = list( title = gene_name, zeroline=FALSE), 
                        #        yaxis = list( title = broad_id, zeroline=FALSE), showlegend=T)
                        )
  }
  return(fig)
}


#######################################
# This function takes the following inputs dependency data in d_data (DepMap), s_data (Sanger)
# gene name in (gene_name) such as "A1CF", drug broad id in (brd_id)
DepLink_G2D_ScatterPlot<- function(drug_data, gene_data1,
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
    fit = lm(d1 ~ g1*dt, data = merged )
    fig <- plot_ly(data = merged, x = ~g1, y = ~d1, type="scatter", color=~dt, alpha=0.5, mode="markers", 
                   colors=c("skyblue","peachpuff"), showlegend=F)
    fig <- fig %>% add_lines(x = ~g1, y = fitted(fit), mode="lines", alpha=30, showlegend=T)
  }
  
  if (nrow(combined) >= 3 & nrow(combined2) < 3) {
    
    merged = combined
    
    fit = lm(d1 ~ g1, data = merged )
    fig <- plot_ly(data = merged, x = ~g1, y = ~d1, type="scatter", color=~dt, alpha=0.5, mode="markers", 
                   colors="skyblue", showlegend=F)
    fig <- fig %>% add_lines(x = ~g1, y = fitted(fit), mode="lines", alpha=30, showlegend=T)
  }
  
  if (nrow(combined) < 3 & nrow(combined2) >= 3) {
    merged = combined2
    fit = lm(d1 ~ g1, data = merged )
    fig <- plot_ly(data = merged, x = ~g1, y = ~d1, type="scatter", color=~dt, alpha=0.5, mode="markers", 
                   colors="peachpuff", showlegend=F)
    fig <- fig %>% add_lines(x = ~g1, y = fitted(fit), mode="lines", alpha=30, showlegend=T)
  }
  
  if (drug_source == "PRISM") {
    y.title = paste(drug_name,"(Log fold change)")
  }else{
    y.title = paste(drug_name,"(Log10 IC50)")
  }
  
  
  
  fig <- fig %>% layout(legend=list(title=list(text=paste0("<b> ",drug_source," vs. </b>"))),
                        xaxis = list( title = paste(gene_name, "(Chronos gene effect score)"), zeroline=FALSE), 
                        yaxis = list( title = y.title, zeroline=FALSE), showlegend=T)
  return(fig)
  
}

# VisNetwork for Gene-Gene association. 
# input is the ggTable (see first table for details), gene (e.g. CDK6), and ndx (visual table entries in ggTable)
DepLink_G2G_ReturnNet <- function( ggTable, gene, ndx ) {
  #  ggTable = ggTable()
  
  nodes <- data.frame( id = 1:11, label = c(gene, ggTable$`Gene Name`[ndx]), shape = "circle", color = c( "skyblue", rep( "orange", 10)) )
  edges <- data.frame( id = 1:10, from = rep(1, 10), to = 2:11, width = abs(10*ggTable$`Corr in Broad`[ndx]), color = "skyblue" )
  edges <- mutate(edges, color = ifelse( ggTable$`Corr in Broad`[ndx] > 0, "red", "blue"))
  
  ledges <- data.frame(color = c("red", "blue"),
                       label = c("r>0", "r<0"),
                       width = 4,
                       arrows.scaleFactor = 0 )
  
  network <- visNetwork(nodes, edges, label = nodes$label) %>%
    visNodes(size = 20 ) %>%
    visLegend(width = 0.1, addEdges = ledges, useGroups = FALSE) %>%
    visInteraction(multiselect = FALSE,navigationButtons = TRUE) %>%
    #visOptions(selectedBy = "label",  highlightNearest = TRUE,  nodesIdSelection = TRUE) %>%
    visPhysics(stabilization = FALSE)
  
  return(network)
  
}

# VisNetwork for Gene-Drug association. 
# input is the gdTable (see panel 1I for table detail), gene (e.g. CDK6), and ndx (visual table entries in gdTable)
DepLink_G2D_ReturnNet <- function( gdTable, gene, ndx ) {

  # create nodes and edges. Different from ggTable, we set drug with different shapes.
  nodes <- data.frame( id = 1:11, label = c(gene, gdTable$`Drug Name`[ndx]), shape = c("circle", rep("diamond", 10)), color = c( "skyblue", rep( "orange", 10)) )
  edges <- data.frame( id = 1:10, from = rep(1, 10), to = 2:11, width = abs(10*gdTable$`Corr with Broad`[ndx]), color = "skyblue" )
  edges <- mutate(edges, color = ifelse( gdTable$`Corr with Broad`[ndx] > 0, "red", "blue"))
  
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


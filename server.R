#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library("rintrojs")
library("shiny")
library("shinyjs")
library("shinyBS")
library("shinydashboardPlus")
library("shinycssloaders")
library("shinydashboard")
library("shinyWidgets")
library("plotly")
library("DT")
library("VennDiagram")
library("rcdk")
library("tidyr")
library("visNetwork")
library("stringr")
library("dplyr")
library("sqldf")
#setwd("~/Desktop/DepLink/")


#Load function
source("function/DepLink_QueryNewDrug_Functions.R")
source("function/DepLink_QueryGene_Functions.R")
source("function/DepLink_OueryMultiGenes_Funciton.R")
source("function/DepLink_QueryDrug_Functions.R")
source("function/DepLink_Plotly_Message_Function.R")

#load HEML
#Welcome page 
# welcome <- read.delim("welcom_page_DepLink.html", stringsAsFactors = F)
# welcome  <- paste(test$X.html., collapse = "")

#Data of Tab1
#Gene
depMap_data = readRDS("data/DepMap_21Q4/CRISPR_GeneEffect.rds")
sanger_data = readRDS("data/SangerCRISPR_Chronos21Q2/GeneEffect.rds")
geneNames <- rownames(depMap_data)
sample_info_mc <-
    read.csv(
        file = "data/sample_info_21Q3_PrimaryTypeFixed.csv",
        sep = ",",
        stringsAsFactors = FALSE,
        header = TRUE
    )

#Drug
load("data/PRISM_19Q4/primary-screen-replicate-collapsed-logfold-change_treatment-info.rdata")
colnames(prism_data) = prism_info$broad_id
load("data/GDSC_2019/gdsc_data_info.rdata")
load("data/Tab3/DrugSearchTable_prism_gdsc.rdata")


#Data of Tab2
drugs.tab2 <- readRDS(file = "data/Tab2/pertsig_19811drugs.rds")
fileDropList <- readRDS(file = "data/Tab2/pertsig_Dropbox_URL.rds")
overlapGs <- intersect(geneNames,fileDropList$Gene )

#Data of Tab4
load("data/Tab4/primary-screen-replicate-collapsed-treatment_19Q4_bits.RData") #idx of bits
drug <- readRDS("data/primary-screen-replicate-collapsed-treatment_19Q4_SampleInfo.RDS") #drug info (exclude no SMILE code drug: BRD-K28042756 & BRD-U08520523)
drug.sim.matrix <- readRDS("data/Tab4/primary-screen-replicate-collapsed-treatment_19Q4_Tanimoto_distance.RDS") #Pre-calculate tanimoto distance
drug.tsne <- readRDS("data/Tab4/primary-screen-replicate-collapsed-treatment_19Q4_TSNE_Table.RDS") #The 3 dimension position of each drug


CancerType <- unique(sample_info_mc$primary_disease_fixed)
Cancer_Type_tab3 <- readRDS("data/Tab3/cancer_type_list.rds")

# Define server logic 
shinyServer(function(input, output, session){
    

    output$welcome <- renderUI({
        includeHTML("welcom_page_DepLink.html")
    })

    output$int_feat <- renderUI({
        includeHTML("help.html")
    })


    # #About
observeEvent(input$About, {
    showModal(modalDialog(
        size = "l",
        #title = "Help",
        includeHTML("welcom_page_DepLink.html"),
        #"Introduction of DepLink",
        easyClose = TRUE,
        footer = tagList(
            actionButton(
                    inputId = "about",
                    label = "Close",
                    icon = icon("xmark")
                )
            )
        ))
    })

#     #Close about page
    observeEvent(input$about, {
        removeModal()
    })

    ######################################


    ##################################

    #Help
    observeEvent(input$Help, {
        showModal(modalDialog(
            size = "l",
            #title = "Help",
            includeHTML("help.html"),
            #"Introduction of DepLink",
            easyClose = TRUE,
            footer = tagList(
                actionButton(
                    inputId = "help.all",
                    label = "Close",
                    icon = icon("xmark")
                )
            )
        ))
    })
#     #Close help all page
    observeEvent(input$help.all, {
        removeModal()
    })

# #Tab1
#     # show help of tab1
    observeEvent(input$Help_tab1, {
        showModal(modalDialog(
            size = "l",
            #title = "Help",
            includeHTML("help_tab1.html"),
            #"Introduction of DepLink",
            easyClose = TRUE,
            footer = tagList(
                actionButton(
                    inputId = "help.tab1",
                    label = "Close",
                    icon = icon("xmark")
                )
            )
        ))
    })

#     #Close help all page
    observeEvent(input$help.tab1, {
        removeModal()
    })

    #input control-1
    gene = reactive({
        if (!input$refresh)
            return()

        gene = isolate(input$gene_id)
        gene
    })

    cancer_type = reactive({
        if (!input$refresh)
            return()

        cancer_type = isolate(input$Type)
        cancer_type
    })
    
    
    
    output$g2g_plot1 <- renderPlotly({
        if (!input$refresh)
            return()

        gene = gene()
        cancer_type = cancer_type()
        DepLink_Density_Plot_DepScores(d_data = depMap_data,
                                       s_data = sanger_data,
                                       gene_name = gene,
                                       cancer_type = cancer_type,
                                       sampleInfo = sample_info_mc)

    })

    output$g2g_plot2 <- renderPlotly({
        if (!input$refresh)
            return()

        gene = gene()

        DepLink_BoxPlot_PanCan(d_data = depMap_data,
                               s_data = sanger_data,
                               gene_name = gene,
                               sampleInfo = sample_info_mc)
    })


    #Tab 1-2
    #Gene-Gene correlation:
    output$g2g_plot3 <- renderPlotly({
        if (!input$refresh)
            return()

        gene = gene()
        cancer_type = cancer_type()
        d_data_t <- as.data.frame(t(as.matrix(depMap_data)))
        s_data_t <- as.data.frame(t(as.matrix(sanger_data)))



        # Compute correlations
        rho_d = DepLink_Codependency(data = d_data_t,
                                     gene_name = gene,
                                     cancer_type= cancer_type,
                                     met = "pearson",
                                     sampleInfo = sample_info_mc)
        rho_s = DepLink_Codependency(data = s_data_t,
                                     gene_name = gene,
                                     cancer_type= cancer_type,
                                     met = "pearson",
                                     sampleInfo = sample_info_mc)

        # Display correlation histogram
        DepLink_Density_Plot_CorCoeff(rho_1 = rho_d,rho_2 = rho_s,gene_name = gene)

    })

    #g-g table
    ggTable <- reactive({
        if (!input$refresh)
            return()

            gene = gene()
            cancer_type = cancer_type()

            d_data_t <- as.data.frame(t(as.matrix(depMap_data)))
            s_data_t <- as.data.frame(t(as.matrix(sanger_data)))

            # Compute correlations
            rho_d = DepLink_Codependency(data = d_data_t,
                                         gene_name =gene,
                                         cancer_type= cancer_type,
                                         met = "pearson",
                                         sampleInfo = sample_info_mc)
            rho_s = DepLink_Codependency(data = s_data_t,
                                         gene_name = gene,
                                         cancer_type= cancer_type,
                                         met = "pearson",
                                         sampleInfo = sample_info_mc)

            ggTable <- data.frame(matrix(data = NA, nrow = length(geneNames)-1, ncol = 5),
                                  stringsAsFactors = F)
            colnames(ggTable) <- c("Gene Name", "Corr in Broad", "Corr in Sanger", "Max Corr",
                                   "Avg Corr")
            ggTable$`Gene Name` <- geneNames[geneNames != gene]
            if (sum(is.na(rho_s))>0) {
                ggTable$`Corr in Broad` <- rho_d
                ggTable$`Max Corr` <- rho_d
                ggTable$`Avg Corr` <- rho_d


            }else{
                ggTable$`Corr in Broad` <- rho_d
                ggTable$`Corr in Sanger` <- rho_s
                ggTable$`Max Corr` <- apply(X = ggTable[,c("Corr in Broad", "Corr in Sanger")],MARGIN = 1,FUN =max)
                ggTable$`Avg Corr` <- apply(X = ggTable[,c("Corr in Broad", "Corr in Sanger")],MARGIN = 1,FUN =mean)

            }

            ggTable[order(ggTable$`Max Corr`, decreasing = TRUE), ]

    })



    output$g2g_table1 <- renderDataTable({


        DT::datatable({
            if (!input$refresh)
                return()

            ggTable = ggTable()

        },
        extensions = c('Buttons', 'FixedColumns'), rownames = FALSE,
        options = list(
            dom = 'Bfrtip',
            buttons = c('excel', 'csv'),
            scrollX = TRUE
        ),
        selection = list(mode = "single"))

    })


    output$g2g_network <- renderVisNetwork({
        if (!input$refresh)
            return()

        ggTable = ggTable()
        gene = gene()
        ndx <- input$g2g_table1_rows_current
        DepLink_G2G_ReturnNet( ggTable = ggTable(), gene = gene(), ndx = ndx )
    } )

    output$g2g_scatter1 <- renderPlotly({
        if (!input$refresh)
            return()

            gene = gene()
            cancer_type = cancer_type()

            d_data_t <- as.data.frame(t(as.matrix(depMap_data)))
            s_data_t <- as.data.frame(t(as.matrix(sanger_data)))

            ggTable = ggTable()
            idx =input$g2g_table1_rows_selected

            if (is.null(idx)) {
                DepLink_Plotly_Message(selected_type = "gene",position = "above")
                # DepLink_Gene_ScatterPlot(d_data = d_data_t,
                #                          depMap_data = depMap_data,
                #                          s_data = s_data_t,
                #                          gene_1 = gene,
                #                          gene_2 = ggTable$`Gene Name`[1],
                #                          cancer_type = cancer_type,
                #                          sampleInfo = sample_info_mc)
            }else{
                DepLink_Gene_ScatterPlot(d_data = d_data_t,
                                         depMap_data = depMap_data,
                                         s_data = s_data_t,
                                         gene_1 = gene,
                                         gene_2 = ggTable$`Gene Name`[idx],
                                         cancer_type = cancer_type,
                                         sampleInfo = sample_info_mc)
            }

    })

    #Tab 1-3


    output$g2d_plot1 <- renderPlotly({
        if (!input$refresh)
            return()

            gene = gene()
            cancer_type = cancer_type()

            d_data_t <- as.data.frame(t(as.matrix(depMap_data)))
            s_data_t <- as.data.frame(t(as.matrix(sanger_data)))
            d_rownames = rownames(d_data_t)
            s_rownames = rownames(s_data_t)

            p_rownames = rownames(prism_data)
            g_rownames = rownames(gdsc_IC50)

            ccle_id = DepLink_PanCan2CcleID(cancer_type = cancer_type, sampleInfo = sample_info_mc)

            # DepMap-prism
            cellNames = intersect(d_rownames,p_rownames)
            d_data = d_data_t[ cellNames,]
            p_data = prism_data[ cellNames,]

            if ( cancer_type != "PanCan") {
                n.DepMap_PRISM = sum(cellNames %in% ccle_id)
            }else{
                n.DepMap_PRISM = length(cellNames)
            }
            depMap_prism_rho1 = DepLink_G2D_Codependency(gene_data  = d_data,
                                                         drug_data = p_data,
                                                         gene_name = gene,
                                                         cancer_type=cancer_type,
                                                         sampleInfo = sample_info_mc)
            # Sanger-prism
            cellNames = intersect(s_rownames,p_rownames)
            s_data = s_data_t[ cellNames,]
            p_data = prism_data[ cellNames,]

            if ( cancer_type != "PanCan") {
                n.Sanger_PRISM = sum(cellNames %in% ccle_id)
            }else{
                n.Sanger_PRISM = length(cellNames)
            }
            sanger_prism_rho2 = DepLink_G2D_Codependency(gene_data = s_data,
                                                         drug_data = p_data,
                                                         gene_name = gene,
                                                         cancer_type=cancer_type,
                                                         sampleInfo = sample_info_mc)
            #DepMap-gdsc
            cellNames = intersect(d_rownames,g_rownames)
            d_data = d_data_t[ cellNames,]
            g_data = gdsc_IC50[ cellNames,]
            if ( input$Type != "PanCan") {
                n.DepMap_GDSC = sum(cellNames %in% ccle_id)
            }else{
                n.DepMap_GDSC = length(cellNames)
            }

            depMap_GDSC_rho3 = DepLink_G2D_Codependency(gene_data = d_data,
                                                        drug_data =  g_data,
                                                        gene_name = input$gene_id,
                                                        cancer_type= input$Type,
                                                        sampleInfo = sample_info_mc)

            #Sanger_gdsc
            cellNames = intersect(s_rownames,g_rownames)
            s_data = s_data_t[ cellNames,]
            g_data = gdsc_IC50[ cellNames,]
            if ( input$Type != "PanCan") {
                n.Sanger_GDSC = sum(cellNames %in% ccle_id)
            }else{
                n.Sanger_GDSC = length(cellNames)
            }
            sanger_GDSC_rho4 = DepLink_G2D_Codependency(gene_data = s_data,
                                                        drug_data =  g_data,
                                                        gene_name = input$gene_id,
                                                        cancer_type=input$Type,
                                                        sampleInfo = sample_info_mc)


            DepLink_G2D_Density_Plot(DepMap_PRISM = depMap_prism_rho1,
                                     Sanger_PRISM =  sanger_prism_rho2,
                                     DepMap_GDSC =  depMap_GDSC_rho3,
                                     Sanger_GDSC =  sanger_GDSC_rho4,
                                     gene_name= input$gene_id )


    })


    gd_table <- reactive({

        if (!input$refresh)
            return()

            gene = gene()
            cancer_type = cancer_type()

            d_data_t <- as.data.frame(t(as.matrix(depMap_data)))
            s_data_t <- as.data.frame(t(as.matrix(sanger_data)))
            d_rownames = rownames(d_data_t)
            s_rownames = rownames(s_data_t)

            p_rownames = rownames(prism_data)
            g_rownames = rownames(gdsc_IC50)

            #prism
            # DepMap-prism
            cellNames = intersect(d_rownames,p_rownames)
            d_data = d_data_t[ cellNames,]
            p_data = prism_data[ cellNames,]

            if ( cancer_type != "PanCan") {
                ccle_id = DepLink_PanCan2CcleID(cancer_type = cancer_type, sampleInfo = sample_info_mc)
            }
            depMap_prism_rho1 = DepLink_G2D_Codependency(gene_data  = d_data,
                                                         drug_data = p_data,
                                                         gene_name = gene,
                                                         cancer_type=cancer_type,
                                                         sampleInfo = sample_info_mc)
            # Sanger-prism
            cellNames = intersect(s_rownames,p_rownames)
            s_data = s_data_t[ cellNames,]
            p_data = prism_data[ cellNames,]


            sanger_prism_rho2 = DepLink_G2D_Codependency(gene_data = s_data,
                                                         drug_data = p_data,
                                                         gene_name = gene,
                                                         cancer_type=cancer_type,
                                                         sampleInfo = sample_info_mc)
            table_prism <- data.frame(rep("PRISM", length(depMap_prism_rho1)),
                                      prism_info$broad_id,
                                      prism_info$name,
                                      depMap_prism_rho1,
                                      sanger_prism_rho2)
            colnames(table_prism)<- c("Data Source","Broad ID", "Drug Name", "Corr with Broad", "Corr with Sanger"  )

            #gdsc
            #DepMap-gdsc
            cellNames = intersect(d_rownames,g_rownames)
            d_data = d_data_t[ cellNames,]
            g_data = gdsc_IC50[ cellNames,]


            depMap_GDSC_rho3 = DepLink_G2D_Codependency(gene_data = d_data,
                                                        drug_data =  g_data,
                                                        gene_name = gene,
                                                        cancer_type=cancer_type,
                                                        sampleInfo = sample_info_mc)

            #Sanger_gdsc
            cellNames = intersect(s_rownames,g_rownames)
            s_data = s_data_t[ cellNames,]
            g_data = gdsc_IC50[ cellNames,]

            sanger_GDSC_rho4 = DepLink_G2D_Codependency(gene_data = s_data,
                                                        drug_data =  g_data,
                                                        gene_name = gene,
                                                        cancer_type=cancer_type,
                                                        sampleInfo = sample_info_mc)

            table_gdsc <- data.frame(rep("GDSC", length(depMap_GDSC_rho3)),
                                     gdsc_info$Broad_ID,
                                     gdsc_info$Drug_Name,
                                     depMap_GDSC_rho3,
                                     sanger_GDSC_rho4)

            colnames(table_gdsc)<- c("Data Source","Broad ID", "Drug Name", "Corr with Broad", "Corr with Sanger" )
            table_output <- rbind(table_prism, table_gdsc)


            if (sum(is.na(table_output $`Sanger DepMap`)) > 0) {
                table_output $`Max Corr` <- table_output $`Corr with Broad`
                table_output $`Avg Corr` <- table_output $`Corr with Sanger`

            }else{
                table_output$`Max Corr` <- apply(X = table_output[, c("Corr with Broad", "Corr with Sanger")],MARGIN = 1,FUN = max)
                table_output$`Avg Corr` <- round(apply(X = table_output[, c("Corr with Broad", "Corr with Sanger")],MARGIN = 1,FUN = mean), 3)

            }


            table_output <- table_output[order(table_output$`Max Corr`,decreasing = TRUE), ]

    })



    output$g2d_table1 <- renderDataTable({
        width.t <- session$clientData$output_g2d_table1_width
        height.t <- session$clientData$output_g2d_table1_height
        DT::datatable({
            if (!input$refresh)
                return()

            table_output = gd_table()
            table_output

        },
        extensions = c('Buttons', 'FixedColumns'), rownames = FALSE,width = width.t,height = height.t,
        options = list(
            dom = 'Bfrtip',
            buttons = c('excel', 'csv'),
            scrollX = TRUE
        ),
        selection = list(mode = "single"))

    })

    #Gene to Drug network
    output$g2d_network <- renderVisNetwork({
        if (!input$refresh)
            return()

        ndx <- input$g2d_table1_rows_current
        gdTable = gd_table()
        gene = gene()
        DepLink_G2D_ReturnNet(gdTable = gdTable, gene = gene, ndx = ndx )
    } )



    output$g2d_scatter1 <- renderPlotly({
        if (!input$refresh)
            return()

            gene = gene()
            cancer_type = cancer_type()


            d_data_t <- as.data.frame(t(as.matrix(depMap_data)))
            s_data_t <- as.data.frame(t(as.matrix(sanger_data)))

            gd_table = gd_table()
            idx =input$g2d_table1_rows_selected

            if (is.null(idx)) {
                DepLink_Plotly_Message(selected_type = "drug",position = "above")

            }else{
                #PRISM
                if (gd_table$`Data Source`[idx] == "PRISM") {
                    DepLink_G2D_ScatterPlot(drug_data = prism_data,
                                            gene_data1 = d_data_t,
                                            gene_data2 = s_data_t,
                                            gene_name = gene,
                                            broad_id = gd_table$`Broad ID`[idx],
                                            drug_name = gd_table$`Drug Name`[idx],
                                            sampleInfo = sample_info_mc,
                                            cancer_type = cancer_type,
                                            drug_source = "PRISM")
                }else{
                    #GDSC
                    DepLink_G2D_ScatterPlot(drug_data = gdsc_IC50,
                                            gene_data1 = d_data_t,
                                            gene_data2 = s_data_t,
                                            gene_name = gene,
                                            broad_id = gd_table$`Broad ID`[idx],
                                            drug_name = gd_table$`Drug Name`[idx],
                                            sampleInfo = sample_info_mc,
                                            cancer_type = cancer_type,
                                            drug_source = "GDSC")
                }
            }



    })

    #Tab2:
    # show help of tab2
    observeEvent(input$Help_tab2, {
        showModal(modalDialog(
            size = "l",
            #title = "Help",
            includeHTML("help_tab2.html"),
            #"Introduction of DepLink",
            easyClose = TRUE,
            footer = tagList(
                actionButton(
                    inputId = "help.tab2",
                    label = "Close",
                    icon = icon("xmark")
                )
            )
        ))
    })

    #Close help tab2 page
    observeEvent(input$help.tab2, {
        removeModal()
    })

    #Example
    observeEvent(input$gene_list_ex, {
        updateTextInput(session, "genelist", value = "CDK4\nCDK6\nCCND1")

    })

    #Runnig massage
    observeEvent(input$submit,{
        GeneList <- unlist(strsplit(input$genelist, split = "\n"))


        if (length(GeneList) <= 10 & length(GeneList) > 1) {
            showModal(modalDialog(
                size = "l",
                h3("Job successfully submitted. Results may take a few minutes to appear...\nPlease click anywhere to dismiss this message and do not refresh the page."),
                #"Introduction of DepLink",
                easyClose = TRUE#,
                # footer =tagList(
                #     actionButton(
                #         inputId = "close_tab2",
                #         label = "Close",
                #         icon = icon("xmark")
                #     )
                # )
            ))

        }else if(length(GeneList) == 1){
            showModal(modalDialog(
                size = "l",
                h3("Please use Module 1: Query a Gene"),
                #"Introduction of DepLink",
                easyClose = TRUE,
                footer =tagList(
                    actionButton(
                        inputId = "close_tab2",
                        label = "Close",
                        icon = icon("xmark")
                    )
                )
            ))
        }else{
            showModal(modalDialog(
                size = "l",
                h3("Please limit to ≤10 genes!"),
                #"Introduction of DepLink",
                easyClose = TRUE,
                footer =tagList(
                    actionButton(
                        inputId = "close_tab2",
                        label = "Close",
                        icon = icon("xmark")
                    )
                )
            ))
        }

    })

    #close tab2 runging message
    observeEvent(input$close_tab2, {
        removeModal()
    })

    #Check how many gene symbol in our data
    output$list <- reactive({
        GeneList <- unlist(strsplit(input$genelist, split = "\n"))
        match_msg <- sum(toupper(overlapGs) %in% toupper(GeneList))


        if(length(GeneList)<=10 & length(GeneList) > 1 & length(GeneList) == match_msg){
            paste("All",match_msg, "genes mapped!")
        }else if (length(GeneList) >10) {
            paste("Please limit to ≤10 genes!")
        }else if(length(GeneList) == 1){
            paste("Please use Module 1: Query a Gene!")

        }else{
            paste(match_msg, "genes mapped! Unmapped genes:",
                  paste(GeneList[!toupper(GeneList) %in% toupper(overlapGs)], collapse = ","))
        }


    })

    GeneList <- reactive({
        if(!input$submit)
            return()

        GeneList <- unlist(strsplit(isolate(input$genelist), split = "\n"))
        GeneList <- GeneList[toupper(GeneList) %in% toupper(overlapGs)]

        if (length(GeneList) <=10 &length(GeneList) >1 ) {
            GeneList
        }


    })

    output$gs_plot1 <- renderPlotly({
        if(!input$submit)
            return()

        GeneList = GeneList()
        if(is.null(GeneList))
            return()


        DepLink_Density_Plot_DepScores(d_data = depMap_data,gene_name = GeneList)
    })


    geneList_object <- reactive({
        if(!input$submit)
            return()

        GeneList = GeneList()
        if(is.null(GeneList))
            return()

        if (length(GeneList)<=10) {
            Direction <- isolate(input$checkdir)

            permDistr <- read.table(file = paste0("data/Tab2/permutation/normdistr_mu_std_", length(GeneList), "genes.txt"), header = FALSE)

           DepLink_MultiGenes(geneList = GeneList,direction = Direction,
                               fileDropList = fileDropList,
                               drugs = drugs.tab2,
                               permDistr = permDistr)
        }

    })


    gsTable <- reactive({
        if(!input$submit)
            return()

        GeneList = GeneList()
        if(is.null(GeneList))
            return()

        geneList_object = geneList_object()

        # table
        table_out_pertSig_P <- as.data.frame(matrix(data = NA, nrow = dim(geneList_object$pertSigs)[1]*2, ncol = dim(geneList_object$pertSigs)[2]))
        # table_out_pertSig_P_header <- as.data.frame(matrix(data = NA, nrow = dim(geneList_object$pertSig)[1]*2, 1))
        for (i in 1:dim(geneList_object$pertSigs)[1]) {
            table_out_pertSig_P[(2*i-1):(2*i),] <- rbind(geneList_object$pertSigs[i,],geneList_object$pertSigs_p[i,])
            rownames(table_out_pertSig_P)[(2*i-1):(2*i)] <- c(rownames(geneList_object$pertSigs)[i], rownames(geneList_object$pertSigs_p)[i])
            #  table_out_pertSig_P_header[(2*i-1):(2*i),] <- c(geneList_object$,paste0(geneList_object$," pVal"))
        }
        # rownames(table_out_pertSig_P) <- table_out_pertSig_P_header

        table_out <- as.data.frame(t(rbind(geneList_object$drugs[1:2,],
                                           table_out_pertSig_P,
                                           geneList_object$meanPerDrug,
                                           geneList_object$pVal_t,
                                           geneList_object$pVal_perm,
                                           geneList_object$drugs[3:4,])))
        rownames(table_out) <- seq(1,nrow(table_out),)
        table_out[,3:(2*length(GeneList)+5)] <- apply(table_out[,3:(2*length(GeneList)+5)], MARGIN = 2, FUN = as.numeric)
        table_out[,seq(from = 3, by = 2, length.out = length(GeneList)+1)] <- round(table_out[,seq(from = 3, by = 2, length.out = length(GeneList)+1)] , 2)
        table_out[,seq(from = 4, by = 2, length.out = length(GeneList)+1)] <- round(table_out[,seq(from = 4, by = 2, length.out = length(GeneList)+1)] , 5)
        table_out[,2*length(GeneList)+5] <- round(table_out[,2*length(GeneList)+5] , 10)

        table_out[,2*length(GeneList)+7] <- paste0('<a href="https://pubchem.ncbi.nlm.nih.gov/compound/',table_out[,2*length(GeneList)+7],'" target= "_blank">PubChem</a>')
        table_out <- cbind(table_out, as.data.frame(t(geneList_object$drugs[5,])))
        table_out <- table_out[order(table_out$`Permutation pVal`,decreasing = TRUE),]

    })


    output$gs_table1 <- renderDataTable({
        # GeneList = GeneList()
        DT::datatable({

            if(!input$submit)
                return()

            gsTable =  gsTable()
            if(is.null(gsTable))
                return()

            gsTable

        },
        extensions = c('Buttons', 'FixedColumns'), rownames = FALSE, escape = FALSE,
        options = list(
            dom = 'Bfrtip',
            buttons = c('excel', 'csv'),
            scrollX = TRUE,
            rowCallback = JS(
                "function(col, data) {",
                "for (i = data.length-5; i < data.length-3; i++) {",
                "if (data[i] < 1){",
                "$('td:eq('+i+')', col).html(data[i].toExponential(1));",
                "}",
                "}",
                "}")
        ),
        selection = list(mode = "single")

        )})

    output$gs_plot2 <- renderPlotly({
        if(!input$submit)
            return()

        GeneList = GeneList()
        if(is.null(GeneList) & is.null(gsTable))
            return()

        gsTable =  gsTable()
        geneList_object = geneList_object()

        permDistr <- geneList_object$permDistr
        #idx <- input$gs_table1_rows_selected

        if (is.null(input$gs_table1_rows_selected)) {
           fig =  DepLink_Plotly_Message(selected_type = "drug",position = "above")

        }else{
            idx.row = input$gs_table1_rows_selected

        idx.drug <- gsTable[idx.row, "Name"]

        idx <- which(drugs.tab2[1,] ==  idx.drug)

        meanPerDrug <- gsTable$`Mean Perturbation Score`[idx.row]

        fig = plot_ly(alpha = 0.1)
        fig <- fig %>% add_trace(x=seq(min(-abs(meanPerDrug),permDistr[1, idx]-3*permDistr[2, idx]),max(abs(meanPerDrug),permDistr[1, idx]+3*permDistr[2, idx]),0.01),
                                 y=dnorm(seq(min(-abs(meanPerDrug),permDistr[1, idx]-3*permDistr[2, idx]),max(abs(meanPerDrug),permDistr[1, idx]+3*permDistr[2, idx]),0.01),
                                         permDistr[1, idx],permDistr[2, idx]),
                                 type="scatter",mode="lines", fill="tozeroy",
                                 line=list(color="skyblue"),
                                 name = idx.drug
        )

        fig <- fig %>% layout(  title=t, font=list(size=10),
                                yaxis=list(title = "Density of permutations"),
                                shapes = vline(x = meanPerDrug),
                                annotations = annotation(x = meanPerDrug,
                                                         y = 0.8*max(dnorm(seq(min(-abs(meanPerDrug),permDistr[1, idx]-3*permDistr[2, idx]),max(abs(meanPerDrug),permDistr[1, idx]+3*permDistr[2, idx]),0.01),
                                                                                        permDistr[1, idx],permDistr[2, idx]))) ,
                                xaxis=list(title= paste0("Mean perturbation score of ",length(GeneList) ,"genes by ",idx.drug),
                                           zeroline=FALSE, showgrid=FALSE),
                                legend = list(x = 100, y = 0.9) ,
                                showlegend = FALSE)
        fig <- fig %>% add_annotations(x = 0,
                                       y =  0.9*max(dnorm(seq(min(-abs(meanPerDrug),permDistr[1, idx]-3*permDistr[2, idx]),max(abs(meanPerDrug),permDistr[1, idx]+3*permDistr[2, idx]),0.01),
                                                          permDistr[1, idx],permDistr[2, idx])),
                                       text = "Random permutations",
                                       xref = "x",
                                       yref = "y",
                                       showarrow = TRUE,
                                       arrowhead = 4,
                                       arrowsize = .5,
                                       ax = 20,
                                       ay = -40)
        }
        fig

    })

    #Tab3
    # show help of tab3
    observeEvent(input$Help_tab3, {
        showModal(modalDialog(
            size = "l",
            #title = "Help",
            includeHTML("help_tab3.html"),
            #"Introduction of DepLink",
            easyClose = TRUE,
            footer = tagList(
                actionButton(
                    inputId = "help.tab3",
                    label = "Close",
                    icon = icon("xmark")
                )
            )

        ))
    })

    #Close help tab3 page
    observeEvent(input$help.tab3, {
        removeModal()
    })

    drug.row <- reactive({
        #drugname = "tamoxifen"
        #drugname = "BAM7" #duplicate  "synephrine"     "metaproterenol" "doxycycline"
        drugname = input$drug_id
        idx.row <- unique(which(drugname == drugsearch.matrix, arr.ind = TRUE)[,1])
        idx.row
    })



    drug.table <- reactive({
        #drug.row = idx.row
        drug.row = drug.row()

        idx.prism <- NA
        m = 1
        idx.gdsc <- NA
        n=1

        for (i in 1:length(drug.row)) {
            if (drug.row[i] <= ncol(prism_data)) {
                #PRISM idx
                idx.prism[m] <- drug.row[i]
                m=m+1
            }else{
                #GDSC idx
                idx.gdsc[n] <- drug.row[i] -ncol(prism_data)
                n=n+1
            }
        }

        #Both
        if (length(idx.prism[!is.na(idx.prism)])>=1 & length(idx.gdsc[!is.na(idx.gdsc)]) >= 1) {
            prism_part = data.frame(rep("PRISM", length(idx.prism)), prism_info[idx.prism, c("broad_id" ,"name" , "dose" ,"screen_id" )],idx.prism)
            colnames(prism_part) <- c("Data Source", "Broad id", "Drug Name", "Dose", "Screen id","idx")
            gdsc_part = data.frame(rep("GDSC", length(idx.gdsc)),gdsc_info[idx.gdsc,c(1,2) ],rep(NA, length(idx.gdsc)),rep(NA, length(idx.gdsc)),idx.gdsc)
            colnames(gdsc_part) <- c("Data Source", "Broad id", "Drug Name", "Dose", "Screen id","idx")
            table_out= rbind(prism_part, gdsc_part)
        }

        #PRISM only
        if (length(idx.prism[!is.na(idx.prism)])>=1 & sum(is.na(idx.gdsc)) >=1) {
            prism_part = data.frame(rep("PRISM", length(idx.prism)), prism_info[idx.prism, c("broad_id" ,"name" , "dose" ,"screen_id" )],idx.prism)
            colnames(prism_part) <- c("Data Source", "Broad id", "Drug Name", "Dose", "Screen id","idx")
            table_out= prism_part
        }

        #GDSC only
        if (sum(is.na(idx.prism)) >=1& length(idx.gdsc[!is.na(idx.gdsc)]) >= 1) {
            gdsc_part = data.frame(rep("GDSC", length(idx.gdsc)),gdsc_info[idx.gdsc, ],rep(NA, length(idx.gdsc)),rep(NA, length(idx.gdsc)),idx.gdsc)
            colnames(gdsc_part) <- c("Data Source", "Broad ID", "Drug Name", "Dose", "Screen id","idx")
            table_out =  gdsc_part
        }
        table_out
    })

    #Return drug result
    output$check_table <- renderDataTable({

        DT::datatable({
            #input$refresh
            drug.table = drug.table()
            drug.table[,-6]

        },#caption ="Select one to start analysis (default: the first drug)",
        rownames = FALSE,
        options = list(
            dom = 'Brti'
        ),
        selection = list(mode = "single"))

    })
    # #Drug name
    # drugname = reactive({
    #     idx = input$check_table_rows_selected
    #     if (is.null(idx))
    #         return()
    #
    #     drugname = isolate(input$drug_id)
    #     drugname
    #     })

    #Tab3-1
    output$d2d_plot1 <- renderPlotly({

        idx = input$check_table_rows_selected
        drug.table = drug.table()

        if (is.null(idx)) {
            return()

        }

        DepLink_Density_Plot_Drug(drug.data.p = prism_data,
                                  drug.data.g = gdsc_IC50,
                                  drug_check_table = drug.table,
                                  drug_check_table_idx =idx ,
                                  cancer_type = input$Type_tab3,
                                  sampleInfo = sample_info_mc)


    })


    output$d2d_plot2 <- renderPlotly({

        idx = input$check_table_rows_selected
        drug.table = drug.table()

        if (is.null(idx)) {
            return()
        }
        DepLink_BoxPlot_PanCan_drug(drug.data.p.t = as.data.frame(t(prism_data)),
                                    drug.data.g.t = as.data.frame(t(gdsc_IC50)),
                                    drug_check_table = drug.table,
                                    drug_check_table_idx = idx,
                                    sampleInfo = sample_info_mc)
    })


    #Tab3-2
    #Drug-drug correlation
    #Gene-Gene correlation:
    output$d2d_plot3 <- renderPlotly({

        #input$refresh
        #gene = gene()
        #cancer_type = cancer_type()
        #d_data_t <- as.data.frame(t(as.matrix(depMap_data)))
        #s_data_t <- as.data.frame(t(as.matrix(sanger_data)))

        idx = input$check_table_rows_selected
        drug.table = drug.table()

        if (is.null(idx)) {
            return()
        }else{

            idx = input$check_table_rows_selected
        }

        # Compute correlations
        if(drug.table$`Data Source`[idx] == "PRISM"){
            corr <- DepLink_Codependency_drug(drug.data =prism_data,
                                              drug_check_table = drug.table,
                                              drug_check_table_idx =idx,
                                              cancer_type = input$Type_tab3,
                                              sampleInfo = sample_info_mc,
                                              met = "pearson")
        }else{
            corr <- DepLink_Codependency_drug(drug.data =gdsc_IC50,
                                              drug_check_table = drug.table,
                                              drug_check_table_idx = idx,
                                              cancer_type = input$Type_tab3,
                                              sampleInfo = sample_info_mc,
                                              met = "pearson")
        }

        # Display correlation density plot


        DepLink_Density_Plot_CorCoeff_Drug(rho = corr,
                                           drug_name = drug.table$`Drug Name`[idx],
                                           drug.type = drug.table$`Data Source`[idx])


    })

    #drug-drug table

    ddTable <- reactive({

        idx = input$check_table_rows_selected
        drug.table = drug.table()

        if (is.null(idx)) {
            return()
        }else{

            idx = input$check_table_rows_selected
        }

        # Compute correlations
        if(drug.table$`Data Source`[idx] == "PRISM"){
            corr <- DepLink_Codependency_drug(drug.data =prism_data,
                                              drug_check_table = drug.table,
                                              drug_check_table_idx =idx,
                                              cancer_type = input$Type_tab3,
                                              sampleInfo = sample_info_mc,
                                              met = "pearson")
            ddTable <- data.frame(matrix(data = NA, nrow = ncol(prism_data)-1, ncol= 6), stringsAsFactors = F)
            colnames(ddTable) <- c("Drug Name", "Broad ID", "Correlation","Dose", "Screen ID", "Target Gene")

            ddTable$`Drug Name` <- prism_info$name[- drug.table$idx[idx]]
            ddTable$`Broad ID` <- prism_info$broad_id[- drug.table$idx[idx]]
            ddTable$Correlation <- corr
            ddTable$Dose  <- prism_info$dose[- drug.table$idx[idx]]
            ddTable$`Screen ID`  <- prism_info$screen_id[- drug.table$idx[idx]]
            ddTable$`Target Gene`  <- prism_info$target[- drug.table$idx[idx]]

        }else{
            corr <- DepLink_Codependency_drug(drug.data =gdsc_IC50,
                                              drug_check_table = drug.table,
                                              drug_check_table_idx = idx,
                                              cancer_type =  "PanCan",
                                              sampleInfo = sample_info_mc,
                                              met = "pearson")

            ddTable <- data.frame(matrix(data = NA, nrow  = ncol(gdsc_IC50)-1, ncol= 3), stringsAsFactors = F)

            colnames(ddTable) <- c("Drug Name", "Broad ID", "Correlation")

            ddTable$`Drug Name` <- gdsc_info$Drug_Name[-drug.table$idx[idx]]
            ddTable$`Broad ID` <-gdsc_info$Broad_ID[-drug.table$idx[idx]]
            ddTable$Correlation <- corr

        }

        ddTable[order(ddTable$Correlation, decreasing = TRUE), ]

    })


    output$d2d_table1 <- renderDataTable({
        width.t <- session$clientData$output_d2d_table1_width
        height.t <- session$clientData$output_d2d_table1_height
        DT::datatable({

            ddTable = ddTable()

        },
        extensions = c('Buttons', 'FixedColumns'), rownames = FALSE,width = width.t,height = height.t,
        options = list(
            dom = 'Bfrtip',
            buttons = c('excel', 'csv'),
            scrollX = TRUE
        ),
        selection = list(mode = "single"))

    },server = TRUE)


    output$d2d_network <- renderVisNetwork({
        #drugname = drugname()
        ddTable = ddTable()

        ndx <- input$d2d_table1_rows_current

        if (is.null(ndx))
            return()
       else{

        DepLink_D2D_ReturnNet(ddTable = ddTable,drug =  input$drug_id,ndx = ndx)

            }

    } )


    output$d2d_scatter1 <- renderPlotly({
        drug.table = drug.table()
        idx.1 = input$check_table_rows_selected

        if (is.null(idx.1)) {
           return()
        }

        ddTable = ddTable()
        idx.2 = input$d2d_table1_rows_selected

        # if (is.null(idx.2)) {
        #     DepLink_Plotly_Message(selected_type = "drug",position = "above")
        # }

        if (drug.table$`Data Source`[idx.1] == "PRISM") {

            if (is.null(idx.2)) {
                DepLink_Plotly_Message(selected_type = "drug",position = "above")
            }else{
            DepLink_Drug_ScatterPlot(drug.data = prism_data,
                                     drug_1 = drug.table$`Broad id`[idx.1],
                                     drugname1 = drug.table$`Drug Name`[idx.1],
                                     drug_2 = ddTable$`Broad ID`[idx.2],
                                     drugname2 = ddTable$`Drug Name`[idx.2],
                                     drug.source = "PRISM",
                                     cancer_type = input$Type_tab3,
                                     sampleInfo = sample_info_mc)}

        }else{
            if (is.null(idx.2)) {
                DepLink_Plotly_Message(selected_type = "drug",position = "above")
            }else{
            DepLink_Drug_ScatterPlot(drug.data = gdsc_IC50,
                                     drug_1 = drug.table$`Broad id`[idx.1],
                                     drugname1 = drug.table$`Drug Name`[idx.1],
                                     drug_2 = ddTable$`Broad ID`[idx.2],
                                     drugname2 = ddTable$`Drug Name`[idx.2],
                                     drug.source = "GDSC",
                                     cancer_type =  input$Type_tab3,
                                     sampleInfo = sample_info_mc)
            }
        }
    })

    #Tab3-3 : Drug to Genes (Sanger & Broad) correlation
    output$d2g_plot1 <- renderPlotly({

        drug.table= drug.table()
        idx = input$check_table_rows_selected
        if (is.null(idx)) {
            return()
        }

        d_data_t <- as.data.frame(t(as.matrix(depMap_data)))
        s_data_t <- as.data.frame(t(as.matrix(sanger_data)))
        d_rownames = rownames(d_data_t)
        s_rownames = rownames(s_data_t)

        p_rownames = rownames(prism_data)
        g_rownames = rownames(gdsc_IC50)

        ccle_id = DepLink_PanCan2CcleID(cancer_type = input$Type_tab3, sampleInfo = sample_info_mc)

        if(drug.table$`Data Source`[idx] == "PRISM"){
            # PRISM-Broad
            cellNames = intersect(d_rownames,p_rownames)
            if ( input$Type_tab3 != "PanCan") {
                rho1_n = sum(cellNames %in% ccle_id)
            }else{
                rho1_n = length(cellNames)
            }

            d_data = d_data_t[ cellNames,]
            p_data = prism_data[ cellNames,]

            rho1 <- DepLink_D2G_Codependency(drug.data = p_data,
                                             gene.data = d_data,
                                             drug_id = drug.table$`Broad id`[idx],
                                             cancer_type = input$Type_tab3,
                                             met = "pearson",
                                             sampleInfo = sample_info_mc)


            # PRISM-Sanger
            cellNames = intersect(s_rownames,p_rownames)
            if ( input$Type_tab3 != "PanCan") {
                rho2_n = sum(cellNames %in% ccle_id)
            }else{
                rho2_n = length(cellNames)
            }

            s_data = s_data_t[ cellNames,]
            p_data = prism_data[ cellNames,]

            rho2 <- DepLink_D2G_Codependency(drug.data = p_data,
                                             gene.data = s_data,
                                             drug_id = drug.table$`Broad id`[idx],
                                             cancer_type = input$Type_tab3,
                                             met = "pearson",
                                             sampleInfo = sample_info_mc)
        }


        if (drug.table$`Data Source`[idx] == "GDSC") {

            #GDSC-Broad
            cellNames = intersect(d_rownames,g_rownames)
            if ( input$Type_tab3 != "PanCan") {
                rho1_n = sum(cellNames %in% ccle_id)
            }else{
                rho1_n = length(cellNames)
            }

            d_data = d_data_t[ cellNames,]
            g_data = gdsc_IC50[ cellNames,]
            rho1 <- DepLink_D2G_Codependency(drug.data = g_data,
                                             gene.data = d_data,
                                             drug_id = drug.table$`Broad id`[idx],
                                             cancer_type = input$Type_tab3,
                                             met = "pearson",
                                             sampleInfo = sample_info_mc)


            #GDSC-Sanger
            cellNames = intersect(s_rownames,g_rownames)

            if ( input$Type_tab3 != "PanCan") {
                rho2_n = sum(cellNames %in% ccle_id)
            }else{
                rho2_n = length(cellNames)
            }
            s_data = s_data_t[ cellNames,]
            g_data = gdsc_IC50[ cellNames,]
            rho2 <- DepLink_D2G_Codependency(drug.data = g_data,
                                             gene.data = s_data,
                                             drug_id = drug.table$`Broad id`[idx],
                                             cancer_type = input$Type_tab3,
                                             met = "pearson",
                                             sampleInfo = sample_info_mc)

        }


        DepLink_D2G_Density_Plot(rho1 =rho1 ,rho1_n = rho1_n, rho2 = rho2, rho2_n =rho2_n ,drug_name =drug.table$`Drug Name`[idx],
                                 drug_type =drug.table$`Data Source`[idx] )


    })

    #Drug to gene table
    dgTable <- reactive({
        drug.table= drug.table()
        idx = input$check_table_rows_selected
        if (is.null(idx)) {
            return()
        }

        d_data_t <- as.data.frame(t(as.matrix(depMap_data)))
        s_data_t <- as.data.frame(t(as.matrix(sanger_data)))
        d_rownames = rownames(d_data_t)
        s_rownames = rownames(s_data_t)

        p_rownames = rownames(prism_data)
        g_rownames = rownames(gdsc_IC50)

        ccle_id = DepLink_PanCan2CcleID(cancer_type = input$Type_tab3, sampleInfo = sample_info_mc)

        if(drug.table$`Data Source`[idx] == "PRISM"){
            # PRISM-Broad
            cellNames = intersect(d_rownames,p_rownames)
            if ( input$Type_tab3 != "PanCan") {
                rho1_n = sum(cellNames %in% ccle_id)
            }else{
                rho1_n = length(cellNames)
            }

            d_data = d_data_t[ cellNames,]
            p_data = prism_data[ cellNames,]

            rho1 <- DepLink_D2G_Codependency(drug.data = p_data,
                                             gene.data = d_data,
                                             drug_id = drug.table$`Broad id`[idx],
                                             cancer_type = input$Type_tab3,
                                             met = "pearson",
                                             sampleInfo = sample_info_mc)


            # PRISM-Sanger
            cellNames = intersect(s_rownames,p_rownames)
            if ( input$Type_tab3 != "PanCan") {
                rho2_n = sum(cellNames %in% ccle_id)
            }else{
                rho2_n = length(cellNames)
            }

            s_data = s_data_t[ cellNames,]
            p_data = prism_data[ cellNames,]

            rho2 <- DepLink_D2G_Codependency(drug.data = p_data,
                                             gene.data = s_data,
                                             drug_id = drug.table$`Broad id`[idx],
                                             cancer_type = input$Type_tab3,
                                             met = "pearson",
                                             sampleInfo = sample_info_mc)


            dgtable_out <- data.frame(rep("PRISM", length(rho1)),
                                      colnames(d_data),
                                      rho1,
                                      rho2)
        }


        if (drug.table$`Data Source`[idx] == "GDSC") {

            #GDSC-Broad
            cellNames = intersect(d_rownames,g_rownames)
            if ( input$Type_tab3 != "PanCan") {
                rho1_n = sum(cellNames %in% ccle_id)
            }else{
                rho1_n = length(cellNames)
            }

            d_data = d_data_t[ cellNames,]
            g_data = gdsc_IC50[ cellNames,]
            rho1 <- DepLink_D2G_Codependency(drug.data = g_data,
                                             gene.data = d_data,
                                             drug_id = drug.table$`Broad id`[idx],
                                             cancer_type = input$Type_tab3,
                                             met = "pearson",
                                             sampleInfo = sample_info_mc)


            #GDSC-Sanger
            cellNames = intersect(s_rownames,g_rownames)

            if ( input$Type_tab3 != "PanCan") {
                rho2_n = sum(cellNames %in% ccle_id)
            }else{
                rho2_n = length(cellNames)
            }
            s_data = s_data_t[ cellNames,]
            g_data = gdsc_IC50[ cellNames,]
            rho2 <- DepLink_D2G_Codependency(drug.data = g_data,
                                             gene.data = s_data,
                                             drug_id = drug.table$`Broad id`[idx],
                                             cancer_type = input$Type_tab3,
                                             met = "pearson",
                                             sampleInfo = sample_info_mc)

            dgtable_out <- data.frame(rep("GDSC", length(rho1)),
                                      colnames(d_data),
                                      rho1,
                                      rho2)
        }

        colnames(dgtable_out) <- c("Drug Source", "Gene Name", "Corr in Broad DepMap", "Corr in Sanger DepMap")
        if (sum(is.na(rho2))>=1) {
            dgtable_out$`Max Corr` <- dgtable_out$`Corrin Broad DepMap`
            dgtable_out$`Avg Corr`<- dgtable_out$`Corrin Broad DepMap`

        }else{
            dgtable_out$`Max Corr` <- round(apply(X = dgtable_out[, c("Corr in Broad DepMap", "Corr in Sanger DepMap")],MARGIN = 1,FUN = max),3)
            dgtable_out$`Avg Corr` <- round(apply(X = dgtable_out[, c("Corr in Broad DepMap", "Corr in Sanger DepMap")],MARGIN = 1,FUN = mean), 3)
        }
        dgtable_out<- dgtable_out[order(dgtable_out$`Max Corr`,decreasing = TRUE), ]
        dgtable_out
    })


    output$d2g_table1 <- renderDataTable({

        DT::datatable({

            dgTable  = dgTable()
            dgTable = dgTable[, -1]

        },
        extensions = c('Buttons', 'FixedColumns'), rownames = FALSE,
        options = list(
            dom = 'Bfrtip',
            buttons = c('excel', 'csv'),
            scrollX = TRUE
        ),
        selection = list(mode = "single"))

    },server = TRUE)


    output$d2g_network <- renderVisNetwork({

        dgtable= dgTable()
        ndx <- input$d2g_table1_rows_current


        if (is.null(ndx) ) {
            return()
        }else{

        DepLink_D2G_ReturnNet(dgtable = dgtable,
                              drug = input$drug_id,
                              ndx = ndx)}

    } )

    output$d2g_scatter1 <- renderPlotly({

        idx.1 = input$check_table_rows_selected
        drug.table = drug.table()

        idx.2 = input$d2g_table1_row_selected

        if (is.null(idx.1)) {
            return()
        }

        dgtable= dgTable()
        idx.2 = input$d2g_table1_rows_selected

        # if (is.null(idx.2)) {
        #     DepLink_Plotly_Message(selected_type = "gene",position = "above")
        # }
        if (drug.table$`Data Source`[idx.1] == "PRISM") {
            if (!is.null(idx.2)) {
            DepLink_D2G_ScatterPlot(drug_data = prism_data,
                                    gene_data1 = as.data.frame(t(depMap_data)),
                                    gene_data2 = as.data.frame(t(sanger_data)),
                                    gene_name =  dgtable$`Gene Name`[idx.2],
                                    broad_id = drug.table$`Broad id`[idx.1],
                                    drug_name = drug.table$`Drug Name`[idx.1],
                                    sampleInfo =  sample_info_mc,
                                    cancer_type = input$Type_tab3,
                                    drug_source = "PRISM" )
            }else{
                DepLink_Plotly_Message(selected_type = "gene",position = "above")
            }

        }else{
            if (!is.null(idx.2)) {
            DepLink_D2G_ScatterPlot(drug_data = gdsc_IC50,
                                    gene_data1 = as.data.frame(t(depMap_data)),
                                    gene_data2 = as.data.frame(t(sanger_data)),
                                    gene_name =  dgtable$`Gene Name`[idx.2],
                                    broad_id = drug.table$`Broad id`[idx.1],
                                    drug_name = drug.table$`Drug Name`[idx.1],
                                    sampleInfo =  sample_info_mc,
                                    cancer_type = input$Type_tab3,
                                    drug_source = "GDSC" )
            }else{
                DepLink_Plotly_Message(selected_type = "gene",position = "above")
            }
        }})

    #Tab4
    # show help of tab3
    observeEvent(input$Help_tab4, {
        showModal(modalDialog(
            size = "l",
            #title = "Help",
            includeHTML("help_tab4.html"),
            #"Introduction of DepLink",
            easyClose = TRUE,
            footer = tagList(
                actionButton(
                    inputId = "help.tab4",
                    label = "Close",
                    icon = icon("xmark")
                )
            )

        ))
    })

    #Close help tab4 page
    observeEvent(input$help.tab4, {
        removeModal()
    })

    #Example
    observeEvent(input$new_drug_ex, {
        updateTextInput(session, "smile_code", value = "CN(C)C(=N)NC(N)=N")

    })


    #New drug similarity data processing
    NewDrug.similarity.output <- reactive({

        smiles.newdrug <- parse.smiles(input$smile_code)
        fps.pubchem <-lapply(X = smiles.newdrug,FUN = get.fingerprint,type=c('pubchem'))
        fps.maccs   <-lapply(smiles.newdrug,get.fingerprint,type=c('maccs'))

        NewDrug_bits_idx <- list(c(fps.pubchem[[1]]@bits, 881+fps.maccs[[1]]@bits))

        NewDrug.similarity <- DepLink_Tanimoto_Distance_Individual(fps.bits.idx1 = NewDrug_bits_idx,
                                                                   compared.idx2 = pubchem_maccs_bits_idx ,
                                                                   nDrug.input = 1,
                                                                   drugName.input = "Query drug",
                                                                   combineMatrix = drug.sim.matrix )

        NewDrug.similarity.output <- merge(x = drug[, c(2,3,6:9,12)],
                                           y = NewDrug.similarity,
                                           by.x= "broad_id",
                                           by.y= "broad_id", sort = F)
        colnames(NewDrug.similarity.output)[ncol(NewDrug.similarity.output)-4] <- "SMILE code"
        NewDrug.similarity.output <- NewDrug.similarity.output[ order(NewDrug.similarity.output$`Query drug`, decreasing = T),]
        colnames(NewDrug.similarity.output)[which(colnames(NewDrug.similarity.output) == "Query drug")] <- "Tanimoto similarity"
        NewDrug.similarity.output$`Tanimoto similarity` <- round(NewDrug.similarity.output$`Tanimoto similarity` ,digits = 3)
        #NewDrug.similarity.output$`Drug to Gene` <- "Link to Query 3"

        NewDrug.similarity.output
    })


    output$ng_plot1 <- renderImage({

        if (is.null(input$smile_code)) {
            return()
        }

        smile.code.object <- parse.smiles(input$smile_code)
        mypngwidth <- session$clientData$output_ng_plot1_width
        mypngheight <- session$clientData$output_ng_plot1_height

        # A temp file to save the output.
        # This file will be removed later by renderImage

        outfile <- tempfile(fileext='.png')

        # Generate the png
        png(outfile)
        depictor <- get.depictor( width=400, height=400,zoom = 2)
        img <- view.image.2d(molecule = smile.code.object[[1]],
                             depictor = depictor)

        plot(NA, xlim = c(0, 2), ylim = c(0, 4), type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",ann = F,bty = 'n')
        rasterImage(image = img,
                    xleft =  0,ybottom =  0,
                    xright = 2, ytop = 4)
        dev.off()

        # Return a list containing the filename
        list(src = outfile,
             contentType = 'image/png',
             width = mypngwidth,
             height = mypngheight
             #alt = "This is alternate text"
        )
    }, deleteFile = TRUE)
    output$sim_drug <- renderDataTable(
        DT::datatable({

            if (is.null(input$smile_code))
                return()
            

            NewDrug.similarity.output = NewDrug.similarity.output()
            OutputTable = NewDrug.similarity.output[,which(!colnames(NewDrug.similarity.output) %in%c("fps1" ,"fps2","intersct"))]
            colnames(OutputTable) <- c("Broad ID", "Drug Name", "MOA", "Target Genes",
                                       "Disease Area", "Indication", "SMILES Code", "Tanimoto Similarity")
            OutputTable
        },
        extensions = c('Buttons', 'FixedColumns'), rownames = FALSE,
        options = list(
            dom = 'Bfrtip',
            buttons = c('excel', 'csv'),
            scrollX = TRUE
        ),
        selection = list(mode = "single")

        ))

    output$tsne_plot <- renderPlotly({


        NewDrug.similarity.output = NewDrug.similarity.output()

        idx = input$sim_drug_rows_selected
        if(!is.null(idx)){

            DepLink_TSNEplot(newDrug.similarity =  NewDrug.similarity.output ,
                             drugTable = drug,
                             drug.tsne = drug.tsne,
                             drug.similarity = drug.sim.matrix,
                             target.drug.id = NewDrug.similarity.output[idx, "broad_id"])
        }else{

            DepLink_Plotly_Message(selected_type = "drug" ,position = "above")
            # DepLink_TSNEplot(newDrug.similarity =  NewDrug.similarity.output ,
            #                  drugTable = drug,
            #                  drug.tsne = drug.tsne,
            #                  drug.similarity = drug.sim.matrix )
        }


    })



    output$ng_plot2 <- renderImage({
        NewDrug.similarity.output = NewDrug.similarity.output()

        mypngwidth <- session$clientData$output_ng_plot2_width
        mypngheight <- session$clientData$output_ng_plot2_height

        idx =input$sim_drug_rows_selected

        if(is.null(idx)){
            img.file <- "data/Tab4/drug_image/Message_plot2.png"
            #img <- readPNG(img.file)

        }else{

            # Generate the png
            img.file <- paste0("data/Tab4/drug_image/",NewDrug.similarity.output$broad_id[idx],".png")
            #img <- readPNG(img.file)

        }

        # Return a list containing the filename
        list(src = img.file,
             width = mypngwidth,
             height = mypngheight ,
             alt = "Please select one drug from the table blow"
        )
    }, deleteFile = FALSE)


    output$vennDia_plot <- renderPlot({
        NewDrug.similarity.output = NewDrug.similarity.output()


        idx = input$sim_drug_rows_selected  # index of selected row

        if(is.null(idx)){

            plot(c(0, 4), c(0, 4), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
            text(x = 2, y = 2, paste("Please select one drug from the table above"),
                 cex = 1, col = "black")


        }else{
            drugNames <- NewDrug.similarity.output[idx,"name"]
            VennDiagram::draw.pairwise.venn(area1 = NewDrug.similarity.output$fps1[idx],
                                            area2 = NewDrug.similarity.output$fps2[idx],
                                            cross.area = NewDrug.similarity.output$intersct[idx],
                                            scaled = FALSE,
                                            category = c("Query drug", drugNames),
                                            lty = rep("solid", 2),
                                            fill = c("#ECF87F", "#FF8300"),
                                            alpha = 0.5,
                                            cat.pos = c(0,0),
                                            cex = 2,
                                            cat.cex = 1)
        }

    })

 })

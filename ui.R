#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#' ---
#' title: DepLink: A missing link between cancer gene dependency and drug response
#' author: Y-C Chiu, T. Nayak, L-J Wang, and Y. Chen
#' date: Novemeber 2021
#' Affiliation: Greehey Children's Cancer Research Institute,
#'              UT Health San Antonio, San Antonio, TX
#' ---

#' The shiny server main program including ui() and server()

# Nov. 2021, app.R created, L-J. Wang
# Nov. 19, 2021, Header added.  Y. Chen

# You may need to install following packages
#install.packages("rintrojs")
#install.packages("shinyBS")
#install.packages("visNetwork")
library("htmltools")
library("htmlwidgets")
library("rintrojs")
library("shiny")
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
library("shinycssloaders")
library("shinyjs")

#setwd("~/Desktop/Test_all/")

#Load function
source("function/DepLink_QueryNewDrug_Functions.R")
source("function/DepLink_QueryGene_Functions.R")
source("function/DepLink_OueryMultiGenes_Funciton.R")

#load Data
#Figure explain
explain <- read.delim("Figure_explain.csv", sep = ",")

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


#Data of Tab4
drug <- readRDS("data/primary-screen-replicate-collapsed-treatment_19Q4_SampleInfo.RDS") #drug info (exclude no SMILE code drug: BRD-K28042756 & BRD-U08520523)



CancerType <- unique(sample_info_mc$primary_disease_fixed)
Cancer_Type_tab3 <- readRDS("data/Tab3/cancer_type_list.rds")
#Header setting
header <- dashboardHeader(title = span(img(src = "DepLink_logo_v1.svg", height= 80,width = 250)
                                       #style = "color:red",strong("DepLink")
),

tags$li(class = "dropdown",
        actionBttn(inputId = "Help",label = "Help",size = "s",  icon = icon("circle-info"),style = "bordered", color = "primary")) 
# tags$li(class = "dropdown",
#         actionBttn(inputId = "About",label = "About",size = "s",  icon = icon("magnifying-glass"),style = "bordered", color = "primary"))

# tags$li(class = "dropdown",
#         tags$a(href = "https://github.com/chenlabgccri/Prep4DeepDEP",
#                tags$img(src = "github_logo.png", height = 20, width = 20))
#)
)


ui <- htmlTemplate("index.html",
                 geneid = selectizeInput(
                       'gene_id',
                       'Gene Symbol',
                       choices = geneNames,
                       selected = "CDK6"

                   ),
                 cancertype =  selectizeInput(
                       inputId = "Type",
                       label = "Cancer Type",
                       choices = c("PanCan",sort(CancerType)),
                       selected = "PanCan"
                   ),
                 submit1 = actionBttn(inputId = "refresh",label = "Submit",
                                      size = "md",  icon = icon("play-circle"),
                                      style = "material-flat", color = "royal"),
                 genetogeneplot1 = shinycssloaders::withSpinner( plotlyOutput(outputId = "g2g_plot1", width = "auto")),
                 genetogeneplot2 = shinycssloaders::withSpinner(plotlyOutput(outputId = "g2g_plot2", width = "auto")),
                 genetogeneplot3 = shinycssloaders::withSpinner(plotlyOutput(outputId = "g2g_plot3" , width = "auto")),
                 genetogenetable1 = shinycssloaders::withSpinner(DT::dataTableOutput("g2g_table1")),
                 genetogenenet1 =  shinycssloaders::withSpinner(visNetworkOutput("g2g_network",height = "390px" )),
                 genetogenescatter1 = shinycssloaders::withSpinner(plotlyOutput(outputId = "g2g_scatter1",height = "390px")),
                 genetodrugplot1 = shinycssloaders::withSpinner(plotlyOutput(outputId = "g2d_plot1" , width = "auto")),
                 genetodrugtable1 = shinycssloaders::withSpinner(DT::dataTableOutput("g2d_table1")),
                 genetodrugnet1 = shinycssloaders::withSpinner(visNetworkOutput("g2d_network",height = "390px" )),
                 genetodrugscatter1 = shinycssloaders::withSpinner(plotlyOutput(outputId = "g2d_scatter1",height = "390px")),
                 regulation = radioButtons(
                                       "checkdir",
                                       label = "Genes to be...",
                                       choices = list("Up-regulated" = 1, "Down-regulated"= 0),
                                       selected = 0
                                   ),
                 list1 = textAreaInput(
                                       inputId = "genelist",
                                       label =  "Gene List (up to 10 genes)",
                                       placeholder = "Paste gene symbols here as separate lines, for example:\nCDK4\nCDK6\nCCND1",
                                       width = "500px",
                                       height = "250px"
                                   ),
                 example1 =  actionBttn(inputId = "gene_list_ex",label = "Example",size = "xs", style = "simple", color = "primary"),
                 submit2 =  actionBttn(inputId = "submit",label = "Submit",size = "xs",  icon = icon("play-circle"),style = "simple", color = "warning"),
                 text1 = verbatimTextOutput("list"),
                 genesplot1 = shinycssloaders::withSpinner(plotlyOutput(outputId = "gs_plot1")),
                 genestable1 = shinycssloaders::withSpinner(DT::dataTableOutput("gs_table1")),
                 genesplot2 = shinycssloaders::withSpinner(plotlyOutput(outputId = "gs_plot2")),
                 drugid = selectizeInput(
                                      'drug_id',
                                      'Drug Name',
                                      choices = sort(drugSearch.lst),
                                      # options = list(
                                      #   placeholder = 'Type to search for drug name',
                                      #   onInitialize = I('function() { this.setValue(""); }')
                                      # )
                                      selected = "PALBOCICLIB"
                                      ),
                 cancertype2 = selectizeInput(
                                      inputId = "Type_tab3",
                                      label = "Cancer Type",
                                      choices = c("PanCan",sort(Cancer_Type_tab3)),
                                      # options = list(
                                      #   placeholder = 'Type to search for cancer type',
                                      #   onInitialize = I('function() { this.setValue(""); }')),
                                      selected = "PanCan"
                                      ),
                 selectsource = dataTableOutput(outputId = "check_table"),
                 drugtodrugplot1 = shinycssloaders::withSpinner(plotlyOutput(outputId = "d2d_plot1")),
                 drugtodrugplot2 =  shinycssloaders::withSpinner(plotlyOutput(outputId = "d2d_plot2")),
                 drugtodrugplot3 = shinycssloaders::withSpinner(plotlyOutput(outputId = "d2d_plot3" , width = "auto")),
                 drugtodrugtable1 = shinycssloaders::withSpinner(DT::dataTableOutput("d2d_table1")),
                 drugtodrugnet1 = shinycssloaders::withSpinner(visNetworkOutput("d2d_network",height = "390px" )),
                 drugtodrugscatter1 = shinycssloaders::withSpinner(plotlyOutput(outputId = "d2d_scatter1",height = "390px")),
                 drugtogeneplot1 = shinycssloaders::withSpinner(plotlyOutput(outputId = "d2g_plot1" , width = "auto")),
                 drugtogenetable1 = shinycssloaders::withSpinner(DT::dataTableOutput("d2g_table1")),
                 drugtogenenet1 = shinycssloaders::withSpinner(visNetworkOutput("d2g_network" ,height = "390px")),
                 drugtogenescatter1 = shinycssloaders::withSpinner(plotlyOutput(outputId = "d2g_scatter1",height = "390px")),
                 smilecode = textInput("smile_code", label =  "New Drug (SMILE code)",
                                        placeholder = "Enter SMILE code...",width = 600),
                 example2 = actionBttn(inputId = "new_drug_ex",label = "Example",size = "md", style = "simple", color = "primary"),
                 newdrugplot1 = shinycssloaders::withSpinner (imageOutput(outputId = "ng_plot1")),
                 newdrugplot2 = shinycssloaders::withSpinner(imageOutput(outputId = "ng_plot2")),
                 vdplot = shinycssloaders::withSpinner(plotOutput(outputId = "vennDia_plot")),
                 simdrugs = shinycssloaders::withSpinner(dataTableOutput(outputId = "sim_drug")),
                 tsne = shinycssloaders::withSpinner(plotlyOutput(outputId = "tsne_plot",width = "auto"))
                 )
    




# Library required packages

# for data processing and analysis
# devtools::install_github("PheWAS/PheWAS")
library(qqman)
library(PheWAS)
library(stringr)
library(plotly)
library(purrr)
library(scales)

library(readr)
library(lubridate)

# for Shiny
library(shiny); library(shinydashboard); library(shinyjs)
library(ggrepel)
library(DT)
library(reshape2)
library(crosstalk)

# for Word cloud generator
library(ECharts2Shiny)

###########
## Capstone project
# shiny app

######################
## Variable Initialization
bar_order <- 0

# original color used in the MAP manuscript
color_vis <- c("blue","darkcyan","brown","darkorange1","magenta","darkblue",
               "darkseagreen4","red","coral4","chartreuse4","black","royalblue4",
               "firebrick","darkolivegreen","mediumspringgreen","purple","gray50")

# convert color to rgb format, for use in DT table
colvis_rgb <- c()
for (ele in color_vis){
  tmp <- col2rgb(ele)
  colvis_rgb <- c(colvis_rgb,str_glue("rgb({tmp[1]},{tmp[2]},{tmp[3]})"))
}

# put a message in console or server log; note this happens only when the app is started!
cat("uiStub application started...\n")

ui <- uiOutput("uiStub")            # single-output stub ui

server <- function(input, output, session) {
  
  cat("Session started.\n")                               # this prints when a session starts
  onSessionEnded(function() {cat("Session ended.\n\n")})  # this prints when a session ends
  
  
  #Initialize the values
  pheinfo1 <- reactiveVal();
  df1 <- reactiveVal();
  vis_df_all1 <- reactiveVal(); groupinfo_df1 <- reactiveVal()
  ratio_df1 <- reactiveVal()
  three_ms1 <- reactiveVal(); three_ms_vd1 <- reactiveVal(); three_mss1 <- reactiveVal()
  
  ##########
  #### Read in the data
  # load("data/pheinfo.rda")  # 1814*5; phecode info
  # I. pheinfo1
  pheinfo1 <- reactive({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    da <- read_csv(inFile$datapath) %>% as.data.frame
    return(da)
    
  })
  
  # Patient 5 and 6 are high in MS  - new patient data
  # load("data/sample_PheMAP.Rdata")   # dataset updated 07/09/2019. name: df, 1864*8 
  # II. df1 
  df1 <- reactive({
    inFile <- input$file2
    if (is.null(inFile))
      return(NULL)
    da <- read_csv(inFile$datapath) %>% as.data.frame
    return(da)
  })
  
  # III. vis_df_all1 
  vis_df_all1 <- reactive({
    
    if(is.null(df1())) return(NULL)
    
    df <- df1()
    phe_man <- df$phecode   # 1864, with 2 already filtered by Luwan; 
    
    vis_df_all <- tibble()
    
    for (i_indi in 2:(ncol(df)-1)){   # starting from the second columns
      
      phe_man_indiv1 <- phe_man %>% tibble %>% mutate(map_prob = df[,i_indi])
      colnames(phe_man_indiv1)[1] <- "phecode"
      
      # convert from e.g. 611_1 to 611.1
      phe_man_indiv1 <- phe_man_indiv1 %>% mutate(phecode = phecode %>% str_replace('_','.') %>% as.character)
      # phe_man_indiv1$phecode = phe_man_indiv1$phecode %>% str_replace('_','.') %>% as.character
      
      # remove the "." at the end of each phecode if there is any
      for (i in 1:length(phe_man_indiv1$phecode)){
        # print(phe_man_indiv1$phecode[i] %>% str_sub(.,start=nchar(.)) == ".")
        if(phe_man_indiv1$phecode[i] %>% str_sub(.,start=nchar(.)) == ".") {
          phe_man_indiv1$phecode[i]  <- substr(phe_man_indiv1$phecode[i],1, nchar(phe_man_indiv1$phecode[i])-1)
        }
      }
      
      #dim(phe_man_indiv1)  # 1864    2
      
      ##############
      # Critical part
      # map to phecode information
      binom <- addPhecodeInfo(phe_man_indiv1,groupnums = TRUE,groupcolors = TRUE)
      colnames(binom)[1] = "phenotype"      #rename the `phecode` column to phenotype
      
      
      binom$pheno <- binom$description
      binom$description <- paste0("Phecode: ",binom$phenotype,
                                  "\nPhenotype: ",binom$description, 
                                  "\nProb: ",binom$map_prob %>% round(3) ,sep="")
      # dim(binom)     # only 1815 6
      
      
      # extract list of significant phecodes for each individual patient 
      # Assume that it would be common that there would be NAs in individual patient's data (phe_man_indiv1) ; 
      # but generally fixed number of NAs in MAPcutoff file (MAPcutoff), 2 NAs in our case, already filtered for downstream analysis (MAPcutoff_filteredNA)
      
      sig_list <- c()
      # [new]
      phe_man_indiv1 <- phe_man_indiv1 %>% mutate(cutoff=df$cut.MAP)
      # phe_man_indiv1$cutoff = df$cut.MAP
      
      for (i in 1:dim(phe_man_indiv1)[1]){
        if (phe_man_indiv1$map_prob[i] > phe_man_indiv1$cutoff[i]) {
          sig_list <- c(sig_list,phe_man_indiv1$phecode[i])
        }
      }
      
      vis_df <- binom %>% 
        arrange(groupnum) %>% 
        mutate(phenotypes = rownames(binom) %>% as.numeric,
               is_highlight=ifelse(phenotype %in% sig_list, "yes", "no"))
      
      vis_df <- merge(vis_df, phe_man_indiv1[,c(1,3)], by.x = "phenotype", by.y = "phecode") %>% arrange(phenotypes) 
      colnames(vis_df)[1] <- "phecode"
      ### append into a final tibble that contains all the info
      vis_df_all <- rbind(vis_df_all,vis_df)  
      print(dim(vis_df_all))
    }
    
    return(vis_df_all)
    
  })
  
  
  
  # PheWAS group info summary
  # IV. groupinfo_df1
  groupinfo_df1 <- reactive({
    
    if(is.null(df1()) | is.null(pheinfo1()) | is.null(vis_df_all1())) return(NULL)
    #read in the data
    df <- df1()
    
    pheinfo <- pheinfo1()
    
    vis_df_all <- vis_df_all1()
    
    # PheWAS group info summary
    groupinfo_df <- table(pheinfo$groupnum) %>% names %>% as.data.frame    # 17 groups (1-15, 17-18) 
    colnames(groupinfo_df) <- "Groupnum"
    
    groupinfo_df$num <- table(pheinfo$groupnum) %>% unname
    
    groupinfo_df$group <- vis_df_all[1:(nrow(vis_df_all)/(ncol(df)-2)),] %>% group_by(group) %>% summarize(groupnum = groupnum[1]) %>% 
      arrange(groupnum) %>% .$group
    
    
    return(groupinfo_df)
  })
  
  
  # # Prepare data for the barchart plot
  # # ratio_df stores the information of the proportion of the phecodes above threshold in each PheWAS group.
  # V. groupinfo_df1
  ratio_df1 <- reactive({
    if(is.null(groupinfo_df1()) | is.null(vis_df_all1()) | is.null(df1()) ) return(NULL)
    #read in the data
    groupinfo_df <- groupinfo_df1()
    vis_df_all <- vis_df_all1()
    df <- df1()
    
    vis_df <- vis_df_all[1:(nrow(vis_df_all)/(ncol(df)-2)),]
    
    tmp_df <- vis_df[which(vis_df$is_highlight == "yes"),]
    tmpp_df <- table(tmp_df$groupnum) %>% names %>% as.data.frame
    colnames(tmpp_df) <- "Groupnum"
    tmpp_df$num <- table(tmp_df$groupnum) %>% unname
    ratio_df <- groupinfo_df
    flag = 0
    for (j in 1:nrow(groupinfo_df)){
      if(is.na(tmpp_df$Groupnum[j-flag]) | (tmpp_df$Groupnum[j-flag] != groupinfo_df$Groupnum[j])){
        ratio_df$num[j] <- 0
        flag = flag + 1
      } else if (tmpp_df$Groupnum[j-flag] == groupinfo_df$Groupnum[j]){
        ratio_df$num[j] <- tmpp_df$num[j-flag]/groupinfo_df$num[j]
      }
    }
    colnames(ratio_df)[2] <- "Proportion_abv_thrh"
    
    tmp <- c()
    for (ele in ratio_df$Proportion_abv_thrh){
      tmp <- c(tmp,percent(ele))
    }
    
    ratio_df$description <- paste0("PheWAS group: ",ratio_df$group,
                                   "\nProportion above threshold: ",
                                   tmp, sep="")
    
    return(ratio_df)
  })
  
  #############
  # for MS data 
  # read in the MS patient files
  # VI. three_ms1
  three_ms1 <- reactive({
    inFile <- input$file3
    if (is.null(inFile))
      return(NULL)
    three_ms <- read_csv(inFile$datapath)
    # reformat the file
    three_ms <- three_ms %>% mutate(Encounter = 1)    # easier to plot
    three_ms$StartDate <- as.Date(three_ms$StartDate, "%m/%d/%Y")
    
    return(three_ms)
  })
  
  # read in the Vitamin D file
  # VII. three_ms_vd1
  three_ms_vd1 <- reactive({
    inFile <- input$file4
    if (is.null(inFile))
      return(NULL)
    three_ms_vd <- read_csv(inFile$datapath)
    three_ms_vd$StartDate <- as.Date(three_ms_vd$StartDate, "%m/%d/%Y")
    
    three_ms_vd$description <- paste0("Start Date: ",three_ms_vd$StartDate, 
                                      "\nPatient Number: ",three_ms_vd$PatientNum,
                                      "\nCategory: Vitamin D CUI",
                                      "\nVitamin D Value: ",three_ms_vd$Value)
    return(three_ms_vd)
  })
  
  
  # VIII. read in ms_CUI and return three_mss
  three_mss1 <- reactive({
    inFile <- input$file5
    if (is.null(inFile) | is.null(three_ms1()))
      return(NULL)
    
    ms_CUI <- read_csv(inFile$datapath)
    three_ms <- three_ms1()
    
    # merge the original data file and CUIs together
    three_mss <- inner_join(three_ms,ms_CUI,by=c("ConceptCd", "Category"))
    
    three_mss$Description <- paste0("Start Date: ",three_mss$StartDate, 
                                    "\nPatient Number: ",three_mss$PatientNum,
                                    "\nCategory: ",three_mss$Category,
                                    "\nDescription: ",three_mss$Description)
    
    # aggregate the encounters by month
    three_mss <- three_mss %>% group_by(Month=floor_date(StartDate, "month"))
    three_mss <- three_mss %>% group_by(Year=floor_date(StartDate, "year"))
    
    
    three_mss$color <- factor(three_mss$Category, labels = RColorBrewer::brewer.pal(length(unique(three_mss$Category)), name = "Set2"))
    
    ####! Important 
    choices <<- unique(three_mss$PatientNum)
    
    return(three_mss)
  })
  
  ######### Read in the data end.
  ###############################
  
  
  
  ############################################################################
  # build menu; only show all pages after users finish uploading all the files
  ## This the default page; 
  ## !!! In order to enable data uploading along the way, don't put the codes inside observe({})
  fname = "home.R"
  cat(paste0("Session filename: ", fname, ".\n"))      # print the URL for this session
  source(fname, local=TRUE)                            # load and run server code for this page

  output$uiStub <- renderUI(tagList(             # a single-output stub ui basically lets you
    dashboardPage(
      dashboardHeader( 
        disable = TRUE
        # title = HTML("<h2><a href='?home'>Home</a>","</h2>")
      ),
      dashboardSidebar(disable = TRUE),
      dashboardBody(
        uiOutput("pageStub")                     # loaded server code should render the rest of the page to this output$
      )
    )) )
  
  
  ## show all pages immediately after users finish uploading all the files
  observe({
    if( !is.null(df1()) & !is.null(vis_df_all1()) & !is.null(groupinfo_df1()) & !is.null(ratio_df1())
        & !is.null(three_mss1()) & !is.null(three_ms_vd1())) {
      
      ################################ !!Important
      ### Write in the "global" environment so that we can get it from map-explorer.R
      df11 <<- df1()
      vis_df_all11 <<- vis_df_all1()
      groupinfo_df11 <<- groupinfo_df1()
      ratio_df11 <<- ratio_df1()
      three_ms_vd11 <<- three_ms_vd1()
      three_mss11 <<- three_mss1()
      
      ###  
      fname = "map-explorer.R"
      cat(paste0("Session filename: ", fname, ".\n"))      # print the URL for this session
      source(fname, local=TRUE)                            # load and run server code for this page
      
      output$uiStub <- renderUI(tagList(             # a single-output stub ui basically lets you
        fluidPage(                                  #     move the ui into the server function
          fluidRow(
            # column(12,
            #        HTML("<h2><a href='?home'>Home</a> | ",
            #             "<a href='?map-explorer'>MAP Explorer</a>",
            #             "</h2>")
            # )
            column(12,
                   HTML("<h5><a href='?home'>Back to Data Upload</a> ","</h5>") )# for column
          ),
          uiOutput("pageStub")                     # loaded server code should render the rest of the page to this output$
        )
      ))
      
      # output$uiStub <- renderUI(tagList(             # a single-output stub ui basically lets you
      #   dashboardPage(
      #     dashboardHeader(
      #       title = HTML("<h5><a href='?home'>Back to Home</a> ","</h5>")   #"Back to Home"
      #     ),
      #     dashboardSidebar(disable = TRUE),
      #     dashboardBody(
      #       uiOutput("pageStub")                     # loaded server code should render the rest of the page to this output$
      #     )
      #   )) )
      
    }  # for the `if`
  })
  ############################################################################
  # finish building menu
  
}


shinyApp(ui, server)


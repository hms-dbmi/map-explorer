
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
library(shiny)
library(shinyjs)
library(DT)
library(reshape2)
library(crosstalk)

# for Word cloud generator
library(ECharts2Shiny)


### Generally speaking, only 'pheinfo.rda' & 'sample_PheMAP.Rdata' are sufficient.
# Read in the data
# `dat` contains the MAP probabilities for each individual patient across all diseases
setwd("~/Desktop/Capstone")
load("data/MAPmanhattan.Rdata")    # dataset name: dat, 4*1864; 
load("data/MAPcutoff.Rdata")       # 1866    2
load("data/pheinfo.rda")           # 1814    5            


# Patient 5 and 6 are high in MS
load("data/sample_PheMAP.Rdata")   # dataset updated 07/09/2019. name: df, 1864*8  
# format:  phecode pat1 pat2 pat3 pat4 pat5 pat6 cut.MAP





# use_new_data=TRUE
# if(use_new_data){...}


phe_man <- df$phecode   # 1864, with 2 already filtered by Luwan; 

# filter out phecodes with NA cutoff, with 2 NAs
MAPcutoff_filteredNA <- MAPcutoff[!is.na(MAPcutoff$cutoff),]     # 1864    2

# PheWAS group info summary
groupinfo_df <- table(pheinfo$groupnum) %>% names %>% as.data.frame    # 17 groups (1-15, 17-18) 
colnames(groupinfo_df) <- "Groupnum"
groupinfo_df$num <- table(pheinfo$groupnum) %>% unname

vis_df_all <- tibble()

for (i_indi in 2:(ncol(df)-1)){   # starting from the second columns
  
  phe_man_indiv1 <- phe_man %>% tibble %>% mutate(map_prob = df[,i_indi]) 
  colnames(phe_man_indiv1)[1] <- "phecode"
  # dim(phe_man_indiv1)      # check dimension
  
  # convert from e.g. 611_1 to 611.1
  phe_man_indiv1 <- phe_man_indiv1 %>% mutate(phecode = phecode %>% str_replace('_','.') %>% as.character)
  # remove the "." at the end of each phecode if there is any
  for (i in 1:length(phe_man_indiv1$phecode)){
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
  
  
  # print out not matched phenotype to make sure none of them are significant   - not true
  # these missing phecodes might not belong to any of the big disease groups.
  # phe_man_indiv1[which(phe_man_indiv1$phecode %in% not_matched_phe),]
  not_matched_phe <- setdiff(phe_man_indiv1$phecode,binom$phenotype) 
  # length(not_matched_phe)    # 52 (1864-1812) , same as what they got
  
  not_matched_phe <- phe_man_indiv1[which(phe_man_indiv1$phecode %in% not_matched_phe),]
  
  # create a new group for storing 52 unmatched phecodes through PheWAS
  # binom <- add_row(binom,
  #                  phenotype = not_matched_phe$phecode[i],
  #                  description = "Other",
  #                  groupnum = 19,
  #                  group = "Others",
  #                  color = "orange",
  #                  map_prob = not_matched_phe$map_prob[i])
  binom$pheno <- binom$description
  binom$description <- paste0("Phecode: ",binom$phenotype,
                              "\nPhenotype: ",binom$description, 
                              "\nProb: ",binom$map_prob %>% round(8) ,sep="")
  # dim(binom)     # only 1812 6(1814-2)
  
  
  # MAPcutoff_filteredNA: 1864  2
  # format the MAPcutoff_filteredNA file
  MAPcutoff_filteredNA <- MAPcutoff_filteredNA %>% mutate(phecode = phecode %>% str_replace('_','.') %>% as.character)
  for (i in 1:length(MAPcutoff_filteredNA$phecode)){
    if(MAPcutoff_filteredNA$phecode[i] %>% str_sub(.,start=nchar(.)) == ".") {
      MAPcutoff_filteredNA$phecode[i]  <- substr(MAPcutoff_filteredNA$phecode[i], 1, nchar(MAPcutoff_filteredNA$phecode[i])-1)
    }
  }
  
  # extract list of significant phecodes for each individual patient 
  # Assume that it would be common that there would be NAs in individual patient's data (phe_man_indiv1) ; 
  # but generally fixed number of NAs in MAPcutoff file (MAPcutoff), 2 NAs in our case, already filtered for downstream analysis (MAPcutoff_filteredNA)
  
  sig_list <- c()
  # [new]
  phe_man_indiv1 <- phe_man_indiv1 %>% mutate(cutoff=df$cut.MAP)
  
  for (i in 1:dim(phe_man_indiv1)[1]){
    if (phe_man_indiv1$map_prob[i] > phe_man_indiv1$cutoff[i]) {
      sig_list <- c(sig_list,phe_man_indiv1$phecode[i])
    }
  }
  
  ######################
  # Actual visualization
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
  
  vis_df <- binom %>% 
    arrange(groupnum) %>% 
    mutate(phenotypes = rownames(binom) %>% as.numeric,
           is_highlight=ifelse(phenotype %in% sig_list, "yes", "no"))
  
  
  groupinfo_df$group <- vis_df %>% group_by(group) %>% summarize(groupnum = groupnum[1]) %>% 
    arrange(groupnum) %>% .$group

  # # Prepare data for the barchart plot
  # # ratio_df stores the information of the proportion of the phecodes above threshold in each PheWAS group.
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
  
  # # Plot the barchart   - plot it in Shiny
  # x = factor(ratio_df$group)   
  # x = factor(x,levels(x)[c(8,12,5,7,10,13,16,1,15,4,6,14,3,11,2,17,9)])    #reorder factor levels
  # 
  # ratio_p <- ggplot(ratio_df, aes(x=x, y=Proportion_abv_thrh, text = description)) + 
  #   geom_bar(stat="identity",fill = color_vis) + 
  #   
  #   ggtitle ("Proportion of phecodes above threshold") +
  #   ylab("Proportion") +
  #   scale_y_continuous(expand = c(0, 0), limits = c(0,max(ratio_df$Proportion_abv_thrh)+0.05)) +
  #   scale_x_discrete(name = "",labels=ratio_df$group) +
  #   
  #   theme_bw() + 
  #   theme(
  #     legend.position="none",
  #     panel.border = element_blank(),
  #     panel.grid.major.x = element_blank(),
  #     panel.grid.minor.x = element_blank(),
  #     axis.title.y = element_text(family = "sans",size=8),
  #     axis.text.x = element_text(family = "sans", angle = 45, hjust = 1, size = 6.5),
  #     plot.title = element_text(family = "sans",size = 10)) 
  # 
  # assign(str_glue("ratio_id{i_indi}"),ratio_p) 
  
  # !!!! Indexed on June 15
  # axisdf <- vis_df %>% group_by(group) %>% summarize(center=( max(phenotypes) + min(phenotypes) ) / 2 )
  
  ## make the manhattan plot   - make it in Shiny 
  # p <- ggplot(vis_df, aes(x = phenotypes, y = map_prob,text = description)) +
  #   
  #   # Show all points
  #   geom_point( aes(color=as.factor(groupnum)), alpha=0.25, size=1.2) +
  #   
  #   # custom X axis:
  #   scale_x_continuous(label = axisdf$group, breaks = axisdf$center) +
  #   scale_y_continuous(name = "MAP Probabilities",expand = c(0, 0),limits = c(0,1.1), breaks = c(0,0.25,0.5,0.75,1)) +     # remove space between plot area and x axis
  #   
  #   # Add highlighted points
  #   geom_point(data=subset(vis_df, is_highlight=="yes"), aes(color=as.factor(groupnum))) +
  #   
  #   scale_color_manual(values = color_vis) +
  #   
  #   # Custom the theme:
  #   theme_bw() +
  #   theme( 
  #     legend.position="none",
  #     panel.border = element_blank(),
  #     panel.grid.major.x = element_blank(),
  #     panel.grid.minor.x = element_blank(),
  #     axis.text.x = element_text(angle = 45, hjust = 1)
  #   )
  # 
  # assign(str_glue("p{i_indi}"),p) 
  
  # ggplotly() function returns a plotly object
  # ggplotly(p, tooltip="text")
  
  vis_df <- merge(vis_df, phe_man_indiv1[,c(1,3)], by.x = "phenotype", by.y = "phecode") %>% arrange(phenotypes) 
  # vis_df <- merge(vis_df, MAPcutoff_filteredNA, by.x = "phenotype", by.y = "phecode") %>% arrange(phenotypes) # add cutoff threshold info into vis_df
  colnames(vis_df)[1] <- "phecode"
  ### append into a final tibble that contains all the info
  vis_df_all <- rbind(vis_df_all,vis_df)   
  print(dim(vis_df_all))
  
}

#############
# for MS data 
# read in the MS patient files
three_ms <- read_csv("data/3_MS_patients.csv")  

# reformat the file
three_ms <- three_ms %>% mutate(Encounter = 1)    # easier to plot
three_ms$StartDate <- as.Date(three_ms$StartDate, "%m/%d/%Y")

# # peek at the data
# unique(three_ms$PatientNum)
# # [1]  68286  99492 106579
# unique(three_ms$PatientID)
# # [1] BWH-487265 BWH-479730 BWH-513823
# # Levels: BWH-479730 BWH-487265 BWH-513823
# unique(three_ms$Category)
# # [1] MS ICD        MS CUI        Brain MRI CUI Relapse CUI   Brain MRI CPT Vitamin D CUI
# # Levels: Brain MRI CPT Brain MRI CUI MS CUI MS ICD Relapse CUI Vitamin D CUI
# 
# # Just 3 patients - pick 2 w lots of encounters, and 1 w fewer. Merge all brain mri. 
# three_ms %>% filter(PatientNum == "68286") %>% dim
# # [1] 270   5 (6)
# three_ms %>% filter(PatientNum == "99492") %>% dim
# # [1] 178   5 (6)
# three_ms %>% filter(PatientNum == "106579") %>% dim
# # [1] 71  5 (6)

# read in the Vitamin D file
three_ms_vd <- read_csv("data/3_MS_patients_VD.csv")   # 48 5
three_ms_vd$StartDate <- as.Date(three_ms_vd$StartDate, "%m/%d/%Y")

three_ms_vd$description <- paste0("Start Date: ",three_ms_vd$StartDate, 
                                  "\nPatient Number: ",three_ms_vd$PatientNum,
                                  "\nCategory: Vitamin D CUI",
                                  "\nVitamin D Value: ",three_ms_vd$Value)

# read in the CUIs lists
ms_CUI <- read_csv("data/MS_CUI_ICD_list.csv")

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

pat_encounter_1 <- which(three_mss$PatientNum == "68286")
pat_encounter_2 <- which(three_mss$PatientNum == "99492")
pat_encounter_3 <- which(three_mss$PatientNum == "106579")
choices <- unique(three_mss$PatientNum)


###########
## Capstone project
# shiny app
bar_order <- 0


ui <- fluidPage(
  titlePanel("MAP Explorer"),
  
  sidebarPanel(
    h4("You can explore the MAP data in this Shiny App!"),    
    selectInput(inputId="individual_id",
                label="Select Patient ID: ",
                choices=c(1:(ncol(df)-2)), 
                # choices=c(1:nrow(dat)), 
                selected = 1, selectize = F),
    plotlyOutput("bar"),
    
    h3("How to use"),
    p("The default visualization shows the results for Patient 1.
      You can also select any Patient ID from the dropdown list.
      Then it would output a barchart showing the information of the proportion of the phecodes above threshold in each PheWAS group of this individual patient.
      It would also output a Manhattan plot showing the MAP probabilities of each phecodes of this selected patient.
      You can hover over the points in Manhattan plot or Barchart to get more information.
      Enjoy playing around with it!"),
    
    br(),
    h4("You can also upload your own data files here!"), 
    fileInput(inputId = "file1"
              , label = "Chooose file for Phecode info:"
              , multiple = F
              , buttonLabel = "Browse..."
              , placeholder = "No file selected"),
    fileInput(inputId = "file2"
              , label = "Chooose file for Patients and MAP cutoff:"
              , multiple = F
              , buttonLabel = "Browse..."
              , placeholder = "No file selected")
    ),
  
  mainPanel(
    useShinyjs(),
    tabsetPanel(
      tabPanel("Main",
               h3("MAP Manhattan Plot"),
               fluidRow(
                 plotlyOutput("plot",height = 500))
               # br(),
               # textOutput("sig_tab"),
               # br(),
               # textOutput("cond_num"),     #report number of phecodes above threshold
               # textOutput("brush"),        #report number of phecodes both above threshold and selected
               # br(),
               # checkboxInput("details","More details (with MAP cutoff):",FALSE),
               # DT::DTOutput("panel"),
               # tags$head(tags$style("#sig_tab{color: black; font-size: 20px;}")),
               # tags$head(tags$style("#cond_num{color: black; font-size: 15px; font-style:italic;}")),  
               # tags$head(tags$style("#brush{color: black; font-size: 14px; font-style:italic;}")))
               
      ),
      tabPanel("Info Table",
               textOutput("sig_tab"),
               br(),
               textOutput("cond_num"),     #report number of phecodes above threshold
               textOutput("brush"),        #report number of phecodes both above threshold and selected
               br(),
               checkboxInput("details","More details (with MAP cutoff):",FALSE),
               DT::DTOutput("panel"),
               tags$head(tags$style("#sig_tab{color: black; font-size: 20px;}")),
               tags$head(tags$style("#cond_num{color: black; font-size: 15px; font-style:italic;}")),  
               tags$head(tags$style("#brush{color: black; font-size: 14px; font-style:italic;}"))
      ),
      
      tabPanel("Word Cloud",
               h3("Word Cloud of Phenotypes based on MAP Probabilities"),
               # for Word Cloud
               # MUST load the ECharts javascript library in advance
               loadEChartsLibrary(),
               tags$div(id="wordcloud", style="width:100%;height:500px;"),
               deliverChart(div_id = "wordcloud")
               
      ),
      tabPanel("Detailed Evidence",
               h3("MS Data Overview"),
               column(width = 6,
                      selectInput(inputId="patient_num",
                                  label="Select Patient Number: ",
                                  choices=choices, 
                                  selected = 1),
                      # get rid of the extra line between two checkboxes
                      tags$style(".shiny-input-container {margin-bottom: 0px} .checkbox { margin-top: 0px; margin-bottom: 0px }"),
                      checkboxInput("show_stackbar","Show Stacked Bar Charts ",FALSE),
                      checkboxInput("show_penc","Show percentage",FALSE),  # Abs encounters -> percentage
                      checkboxInput("show_comp","Show comparisons between adjacent years",FALSE)), 
               column(width = 6,
                      selectInput(inputId="enco_type",
                                  label="Select Encounter type(s): ",
                                  choices=unique(three_mss$Category), multiple = T,
                                  selected = unique(three_mss$Category))),
               br(), br(), br(), br(), br(), br(), br(),
               plotlyOutput("dat_year", height = 300),
               br(),
               plotlyOutput("dat_month", height = 300),
               br(),
               plotlyOutput("dat_daily", height = 300)
      ),
      tabPanel("Detailed Evidence I try",
               h3("MS Data Overview"),
               column(width = 6,
                      selectInput(inputId="patient_num2",
                                  label="Select Patient Number: ",
                                  choices=choices, 
                                  selected = 1),
                      # get rid of the extra line between two checkboxes
                      tags$style(".shiny-input-container {margin-bottom: 0px} .checkbox { margin-top: 0px; margin-bottom: 0px }"),
                      checkboxInput("show_stackbar2","Show Stacked Bar Charts ",FALSE),
                      checkboxInput("show_penc2","Show percentage",FALSE),  # Abs encounters -> percentage
                      checkboxInput("show_comp2","Show comparisons between adjacent years",FALSE)), 
               column(width = 6,
                      selectInput(inputId="enco_type2",
                                  label="Select Encounter type(s): ",
                                  choices=unique(three_mss$Category), multiple = T,
                                  selected = unique(three_mss$Category))),
               br(), br(), br(), br(), br(), br(), br(),
               plotlyOutput("dat_all", height = 500),
               uiOutput("back"),
               plotlyOutput("dat_comp", height = 500)
      ),
      tabPanel("Detailed Evidence II",
               h3("MS Data Overview"),
               column(width = 6,
                      selectInput(inputId="patient_vd_num",
                                  label="Select Patient Number: ",
                                  choices=choices, 
                                  selected = 1)), # modify the size of the input box
               column(width = 6,
                      selectInput(inputId="enco_vd_type",
                                  label="Select Encounter type(s): ",
                                  choices=unique(three_mss$Category), multiple = T,
                                  selected = unique(three_mss$Category))),
               br(), br(), br(), br(), br(), br(), br(),
               h4("Encounters by Day"),
               plotlyOutput("all_six", height = 600),
               br(),
               h4("Vitamin D Levels"),
               plotlyOutput("vitd"))
    ) # for tabsetPanel
    
  )) # ")" for mainPanel & fluidPage


server <- function(input, output, session) {
  
  # For `Detailed Evidence I try` tabset
  # for maintaining the state of drill-down variables
  
  ###  leverage `customdata` as a workaround
  
  
  # when select a new patient, only show the annual overview
  # observeEvent(input$patient_num2, {
  #   dat_all(NULL)
  # })
  # 
  # observeEvent(input$enco_type2,{
  #   
  #   output$dat_all <- renderPlotly({
  #     id <- match(input$patient_num2, choices)   # choices=c(68286,99492,106579)
  #     pat_encounter <- which(three_mss$PatientNum == choices[id])
  #     
  #     keep_category <- unique(three_mss$Category)[unique(three_mss$Category) %in% input$enco_type]
  #     
  #     select_df <- three_mss[pat_encounter,] %>% filter(Category %in% keep_category)
  #     
  #     years <- unique(select_df$Year)
  #     months <- unique(select_df$Month)
  #     days <- unique(select_df$StartDate)
  # 
  #     k <- select_df %>% count(Category,Year,color)
  #     k$Description <- paste0("Year: ",str_sub(k$Year,1,4), 
  #                             "\nCategory: ",k$Category,
  #                             "\nEncounter: ",k$n)
  #     
  #     #Default:group bar charts & abs # of encounters
  #     yaxi=list(title='Encounter', visible=T); my_barmode='group' ;title="Encounters Aggregated by Year"  
  #     
  #     title_style <- list(text=title,xanchor="left", yanchor="top",showarrow=F,xref = "paper",
  #                         yref = "paper", align = "center",x = -0.05, y = 1.15, font=list(size=16,color='black'))
  #     
  #     sd1 <- SharedData$new(k)
  #     
  #     pc <- sd1 %>%
  #       plot_ly(source = "dat_all", name =~Category, # name of the legend
  #               x = ~Year, y = ~n, color=~color, type="bar",    # ensure each Category has unique color
  #               text=~Description, hoverinfo="text") %>%  #,opacity=~opacity
  #       #add_bars(x = ~Year, y = ~Encounter, color=~color) %>% # ensure each Category has unique color
  #       layout(barmode=my_barmode, 
  #              # title = "Encounters Aggregated by Year",
  #              annotations=title_style,    # change the position of the title of plot_ly in r to the top left of the plot
  #              yaxis=yaxi,
  #              xaxis=list(title='Year', rangeslider=list(type="date"), visible=T))
  #     
  #     
  #     if (is.null(dat_all())) {
  #       pic_year <<- pc
  #       return(pic_year) 
  #     } 
  #     pc 
  #     
  #   })
  # 
  # }) 
  #   
  
  
  # For Detailed Evidence tabset
  # for maintaining the state of drill-down variables
  dat_year <- reactiveVal()
  dat_month <- reactiveVal()
  dat_daily <- reactiveVal()
  
  # when clicking on a bar chart, zoom into the next subcategory (year->month, month->day)
  observeEvent(event_data("plotly_click", source = "dat_year"), {
    dat_year(event_data("plotly_click", source = "dat_year")$x)
    dat_month(NULL)
    dat_daily(NULL)
  })
  
  observeEvent(event_data("plotly_click", source = "dat_month"), {
    dat_month(event_data("plotly_click", source = "dat_month")$x)
    dat_daily(NULL)
  })
  
  # when select a new patient, only show the annual overview
  observeEvent(input$patient_num, {
    dat_year(NULL) 
    dat_month(NULL)
    dat_daily(NULL)
  })
  
  
  observeEvent(input$enco_type,{
    
    output$dat_year <- renderPlotly({
      
      id <- match(input$patient_num, choices)
      pat_encounter <- which(three_mss$PatientNum == choices[id])
      
      keep_category <- unique(three_mss$Category)[unique(three_mss$Category) %in% input$enco_type]
      
      select_df <- three_mss[pat_encounter,] %>% filter(Category %in% keep_category)
      
      k <- select_df %>% count(Year,Category,color)
      k$Description <- paste0("Year: ",str_sub(k$Year,1,4), 
                              "\nCategory: ",k$Category,
                              "\nEncounter: ",k$n)
      
      #Default:group bar charts & abs # of encounters
      yaxi=list(title='Encounter', visible=T); my_barmode='group' ;title="Encounters Aggregated by Year"  
      
      
      # consider 8 (2*2*2) different scenarios: stacked/grouped bars; percentage/raw encounteres; show comparison or not
      if(input$show_stackbar == TRUE) my_barmode='stack'
      if(input$show_penc == TRUE) {
        # percent
        k1=select_df %>% count(Year)
        k_lj=left_join(k,k1,by="Year")
        k$n = k_lj$n.x/k_lj$n.y; 
        k$Description <- paste0("Year: ",str_sub(k$Year,1,4), 
                                "\nCategory: ",k$Category,
                                "\nPercentage: ",percent(k$n))
        
        yaxi=list(title='Percentage per Year', visible=T,tickformat = "%")
        # k_perc <<- k
      }
      
      ##############
      ##############
      ##Make Comparison - Year
      target=c('Brain MRI CPT','Brain MRI CUI','MS CUI','MS ICD','Relapse CUI','Vitamin D CUI')
      if(input$show_comp == TRUE){
        
        first_else=0; k$comp=0; flag=0 
        for(i in 1:nrow(k)){
          yr_new=k$Year[i]
          #1st if - only run in the 1st i 
          if(i==1) {yr=k$Year[i]; info_new=c(); info=c()} 
          
          #2nd if - only run in the last i
          if(i==nrow(k)) {
            yr=k$Year[i]-1  #make it bypass the third if and go to else
            if(flag==0){
              info[k$Category[i]]=k$n[i]
            }else{
              info_new[k$Category[i]]=k$n[i]
            }
          }
          
          #3rd big if  
          if(yr_new==yr) {
            if(flag==0){
              info[k$Category[i]]=k$n[i]
            }else{
              info_new[k$Category[i]]=k$n[i]
            }
            
          }else { # for making comparison (& if two yrs don't match, move to store in new yrs)
            yr=k$Year[i]
            if(flag==0) {  # 0->1
              info_odd=info; info_oddk=info     #as flag=0 finishes storing entries for that year
              info=c()                          # set it to default
              info_new[k$Category[i]]=k$n[i]; flag=1  #"even" yr - info_new (flag=1)
              
            }else{         # 1->0
              info_even=info_new; info_evenk=info_new     #as flag=1 finishes storing entries for that year
              info_new=c()                      # set it to default
              info[k$Category[i]]=k$n[i]; flag=0      #"odd" yr - info (flag=0)
            }
            
            if(first_else==1){  # starting from the 2nd time of unmatch of two yrs, run the following code chunk
              if(flag==0){   #actually is flag==1 as it is set to 0 right above; [even-odd]
                info_even[!(names(info_even) %in% names(info_odd))]=0
                info_comp = c(info_even[names(info_even) %in% names(info_odd)]-info_odd[names(info_odd) %in% names(info_even)],
                              info_even[!(names(info_even) %in% names(info_odd))])
                #Arrange the list in target order
                tmp=info_comp[match(target,names(info_comp))]; tmp=tmp[!is.na(tmp)]
                if(i!=nrow(k)){
                  k$comp[(i-length(info_comp)):(i-1)]=tmp  
                }else{  # account for the different scenario in the last yr
                  k$comp[(i-length(info_comp)+1):i]=tmp  
                }
                
                info_even=info_evenk  # undo the change made to info_even
              }else{    # [odd-even]
                info_odd[!(names(info_odd) %in% names(info_even))]=0
                info_comp = c(info_odd[names(info_odd) %in% names(info_even)]-info_even[names(info_even) %in% names(info_odd)],
                              info_odd[!(names(info_odd) %in% names(info_even))])
                #Arrange the list in target order
                tmp=info_comp[match(target,names(info_comp))]; tmp=tmp[!is.na(tmp)]
                if(i!=nrow(k)){
                  k$comp[(i-length(info_comp)):(i-1)]=tmp  
                }else{    # account for the different scenario in the last yr
                  k$comp[(i-length(info_comp)+1):i]=tmp  
                }
                info_odd=info_oddk   # undo the change made to info_odd
              }
            }
            first_else=1   
          } # a big else 
        }  # for loop
        
        k$n=k$comp  # override n to comp
        
        # show the text when hovered over
        if(input$show_penc == TRUE){
          k$Description <- paste0("Year: ",str_sub(k$Year,1,4), 
                                  "\nCategory: ",k$Category,
                                  "\nDifference: ",percent(k$n))
        }else{
          k$Description <- paste0("Year: ",str_sub(k$Year,1,4), 
                                  "\nCategory: ",k$Category,
                                  "\nDifference: ",k$n)
        }
        
        title='Difference on Encounters from Two Consecutively-Recorded Years'
      }  # the input if
      
      
      title_style <- list(text=title,xanchor="left", yanchor="top",showarrow=F,xref = "paper",
                          yref = "paper", align = "center",x = -0.05, y = 1.15, font=list(size=16,color='black'))
      
      sd1 <- SharedData$new(k)
      
      pc <- sd1 %>%
        plot_ly(source = "dat_year", name =~Category, # name of the legend
                x = ~Year, y = ~n, color=~color, type="bar",    # ensure each Category has unique color
                text=~Description, hoverinfo="text") %>%  #,opacity=~opacity
        #add_bars(x = ~Year, y = ~Encounter, color=~color) %>% # ensure each Category has unique color
        layout(barmode=my_barmode, 
               # title = "Encounters Aggregated by Year",
               annotations=title_style,    # change the position of the title of plot_ly in r to the top left of the plot
               yaxis=yaxi,
               xaxis=list(title='Year', rangeslider=list(type="date"), visible=T))
      
      
      if (is.null(dat_year())) {
        pic_year <<- pc
        return(pic_year) 
      } 
      pc 
      # bscols(
      #   widths=c(3,NA),
      #   list(
      #     filter_checkbox('category', 'Category', sd1, ~Category, inline=F)
      #   ),
      #   pc)
    })
    
    output$dat_month <- renderPlotly({
      if (is.null(dat_year())) return(NULL)
      
      id <- match(input$patient_num, choices)
      pat_encounter <- which(three_mss$PatientNum == choices[id])
      
      keep_category <- unique(three_mss$Category)[unique(three_mss$Category) %in% input$enco_type]
      
      sd <- three_mss[pat_encounter,] %>% 
        filter(Year == dat_year(), Category %in% keep_category)
      
      yyear <- sd$Year[1] %>% substr(1, 4)   # which year is clicked
      # sd2 <- SharedData$new(sd)
      
      k <- sd %>% count(Month,Category,color)
      
      k$Description <- paste0("Month: ",str_sub(k$Month,1,7), 
                              "\nCategory: ",k$Category,
                              "\nEncounter: ",k$n)
      
      
      #Default:group bar charts & abs # of encounters
      yaxi=list(title='Encounter', visible=T); my_barmode='group'; title=paste0("Encounters Aggregated by Month (Year ", yyear,")")
      
      
      # consider 8 (2*2*2) different scenarios: stacked/grouped bars; percentage/raw encounteres; show comparison or not
      if(input$show_stackbar == TRUE) my_barmode='stack'
      if(input$show_penc == TRUE) {
        # percent
        k1= sd %>% count(Month)
        k_lj=left_join(k,k1,by="Month")
        k$n = k_lj$n.x/k_lj$n.y; 
        k$Description <- paste0("Month: ",str_sub(k$Month,1,7), 
                                "\nCategory: ",k$Category,
                                "\nPercentage: ",percent(k$n))
        
        yaxi=list(title='Percentage per Month', visible=T,tickformat = "%")
        
        # k_ori <<- k
      }
      
      
      ##############
      ##############
      ##Make Comparison - month
      target=c('Brain MRI CPT','Brain MRI CUI','MS CUI','MS ICD','Relapse CUI','Vitamin D CUI')
      if(input$show_comp == TRUE){
        
        first_else=0; k$comp=0; flag=0 
        for(i in 1:nrow(k)){
          yr_new=k$Month[i]
          #1st if - only run in the 1st i 
          if(i==1) {yr=k$Month[i]; info_new=c(); info=c()} 
          
          #2nd if - only run in the last i
          if(i==nrow(k)) {
            yr=k$Month[i]-1  #make it bypass the third if and go to else
            if(flag==0){
              info[k$Category[i]]=k$n[i]
            }else{
              info_new[k$Category[i]]=k$n[i]
            }
          }
          
          #3rd big if  
          if(yr_new==yr) {
            if(flag==0){
              info[k$Category[i]]=k$n[i]
            }else{
              info_new[k$Category[i]]=k$n[i]
            }
            
          }else { # for making comparison (& if two months don't match, move to store in new yrs)
            yr=k$Month[i]
            if(flag==0) {  # 0->1
              info_odd=info; info_oddk=info     #as flag=0 finishes storing entries for that month
              info=c()                          # set it to default
              info_new[k$Category[i]]=k$n[i]; flag=1  #"even" month - info_new (flag=1)
              
            }else{         # 1->0
              info_even=info_new; info_evenk=info_new     #as flag=1 finishes storing entries for that month
              info_new=c()                      # set it to default
              info[k$Category[i]]=k$n[i]; flag=0      #"odd" month - info (flag=0)
            }
            
            if(first_else==1){  # starting from the 2nd time of unmatch of two months, run the following code chunk
              if(flag==0){   #actually is flag==1 as it is set to 0 right above; [even-odd]
                info_even[!(names(info_even) %in% names(info_odd))]=0
                info_comp = c(info_even[names(info_even) %in% names(info_odd)]-info_odd[names(info_odd) %in% names(info_even)],
                              info_even[!(names(info_even) %in% names(info_odd))])
                #Arrange the list in target order
                tmp=info_comp[match(target,names(info_comp))]; tmp=tmp[!is.na(tmp)]
                if(i!=nrow(k)){
                  k$comp[(i-length(info_comp)):(i-1)]=tmp  
                }else{  # account for the different scenario in the last month
                  k$comp[(i-length(info_comp)+1):i]=tmp  
                }
                
                info_even=info_evenk  # undo the change made to info_even
              }else{    # [odd-even]
                info_odd[!(names(info_odd) %in% names(info_even))]=0
                info_comp = c(info_odd[names(info_odd) %in% names(info_even)]-info_even[names(info_even) %in% names(info_odd)],
                              info_odd[!(names(info_odd) %in% names(info_even))])
                #Arrange the list in target order
                tmp=info_comp[match(target,names(info_comp))]; tmp=tmp[!is.na(tmp)]
                if(i!=nrow(k)){
                  k$comp[(i-length(info_comp)):(i-1)]=tmp  
                }else{    # account for the different scenario in the last month
                  k$comp[(i-length(info_comp)+1):i]=tmp  
                }
                info_odd=info_oddk   # undo the change made to info_odd
              }
            }
            first_else=1   
          } # a big else 
        }  # for loop
        
        k$n=k$comp  # override n to comp
        
        # show the text when hovered over
        if(input$show_penc == TRUE){
          k$Description <- paste0("Month: ",str_sub(k$Month,1,7), 
                                  "\nCategory: ",k$Category,
                                  "\nDifference: ",percent(k$n))
        }else{
          k$Description <- paste0("Month: ",str_sub(k$Month,1,7), 
                                  "\nCategory: ",k$Category,
                                  "\nDifference: ",k$n)
        }
        title=paste0("Difference on Encounters from Two Consecutively-Recorded Months (Year ", yyear,")")
        
        # k <<- k
        
      }  # the input if
      ## The end of the comparison
      
      title_style <- list(text=title,xanchor="left", yanchor="top",showarrow=F,xref = "paper",
                          yref = "paper", align = "center",x = -0.05, y = 1.15, font=list(size=16,color='black'))
      
      # when a year with entries only from one date, still show bars  
      if(length(unique(k$Month))==1) {nr=nrow(k); k[nr+1,] = k[nr,]; k$Month[nr+1]=k$Month[nr]+10; k$n[nr+1]=0}
      
      sd2 <- SharedData$new(k)
      
      pc <- sd2 %>%
        plot_ly(source = "dat_month", name =~Category, type='bar',
                x = ~Month, y = ~n, color=~color, 
                text=~Description, hoverinfo="text") %>% 
        layout(barmode=my_barmode, annotations=title_style,
               yaxis=yaxi, xaxis=list(range=c(paste0(yyear,"-01-01"),paste0(yyear,"-12-31")),title='Month', 
                                      rangeslider=list(type="date",range=c(paste0(yyear,"-01-01"),paste0(yyear,"-12-31"))), visible=T))
      
      if (is.null(dat_month())) {
        pic_month <<- pc
        return(pic_month)
      }
      pc
      
    })
    
    output$dat_daily <- renderPlotly({
      if (is.null(dat_month())) return(NULL)
      
      id <- match(input$patient_num, choices)
      pat_encounter <- which(three_mss$PatientNum == choices[id])
      
      keep_category <- unique(three_mss$Category)[unique(three_mss$Category) %in% input$enco_type]
      
      sd <- three_mss[pat_encounter,] %>% 
        filter(Month == dat_month(), Category %in% keep_category)
      mmonth <- sd$Month[1] %>% substr(1, 7)
      # sd3 <- SharedData$new(sd)
      
      k <- sd %>% count(StartDate,Category,color)
      k$Description <- paste0("Date: ",k$StartDate, 
                              "\nCategory: ",k$Category,
                              "\nEncounter: ",k$n)
      
      #Default:group bar charts & abs # of encounters
      yaxi=list(title='Encounter', visible=T); my_barmode='group'; title=paste0("Encounters by Day (Month ", mmonth,")")
      
      # k_ori <<- k
      
      # consider 8 (2*2*2) different scenarios: stacked/grouped bars; percentage/raw encounteres; show comparison or not
      if(input$show_stackbar == TRUE) my_barmode='stack'
      if(input$show_penc == TRUE) {
        # percent
        k1= sd %>% count(StartDate)
        k_lj=left_join(k,k1,by="StartDate")
        k$n = k_lj$n.x/k_lj$n.y; 
        k$Description <- paste0("Date: ",k$StartDate, 
                                "\nCategory: ",k$Category,
                                "\nPercentage: ",percent(k$n))
        
        yaxi=list(title='Percentage per Day', visible=T,tickformat = "%")
      }
      
      ##############
      ##############
      ##Make Comparison - Day
      target=c('Brain MRI CPT','Brain MRI CUI','MS CUI','MS ICD','Relapse CUI','Vitamin D CUI')
      if(input$show_comp == TRUE){
        
        first_else=0; k$comp=0; flag=0 
        for(i in 1:nrow(k)){
          yr_new=k$StartDate[i]
          #1st if - only run in the 1st i 
          if(i==1) {yr=k$StartDate[i]; info_new=c(); info=c()} 
          
          #2nd if - only run in the last i
          if(i==nrow(k)) {
            yr=k$StartDate[i]-1  #make it bypass the third if and go to else
            if(flag==0){
              info[k$Category[i]]=k$n[i]
            }else{
              info_new[k$Category[i]]=k$n[i]
            }
          }
          
          #3rd big if  
          if(yr_new==yr) {
            if(flag==0){
              info[k$Category[i]]=k$n[i]
            }else{
              info_new[k$Category[i]]=k$n[i]
            }
            
          }else { # for making comparison (& if two days don't match, move to store in new days)
            yr=k$StartDate[i]
            if(flag==0) {  # 0->1
              info_odd=info; info_oddk=info     #as flag=0 finishes storing entries for that day
              info=c()                          # set it to default
              info_new[k$Category[i]]=k$n[i]; flag=1  #"even" day - info_new (flag=1)
              
            }else{         # 1->0
              info_even=info_new; info_evenk=info_new     #as flag=1 finishes storing entries for that day
              info_new=c()                      # set it to default
              info[k$Category[i]]=k$n[i]; flag=0      #"odd" day - info (flag=0)
            }
            
            if(first_else==1){  # starting from the 2nd time of unmatch of two days, run the following code chunk
              if(flag==0){   #actually is flag==1 as it is set to 0 right above; [even-odd]
                info_even[!(names(info_even) %in% names(info_odd))]=0
                info_comp = c(info_even[names(info_even) %in% names(info_odd)]-info_odd[names(info_odd) %in% names(info_even)],
                              info_even[!(names(info_even) %in% names(info_odd))])
                #Arrange the list in target order
                tmp=info_comp[match(target,names(info_comp))]; tmp=tmp[!is.na(tmp)]
                if(i!=nrow(k)){
                  k$comp[(i-length(info_comp)):(i-1)]=tmp  
                }else{  # account for the different scenario in the last day
                  k$comp[(i-length(info_comp)+1):i]=tmp  
                }
                
                info_even=info_evenk  # undo the change made to info_even
              }else{    # [odd-even]
                info_odd[!(names(info_odd) %in% names(info_even))]=0
                info_comp = c(info_odd[names(info_odd) %in% names(info_even)]-info_even[names(info_even) %in% names(info_odd)],
                              info_odd[!(names(info_odd) %in% names(info_even))])
                #Arrange the list in target order
                tmp=info_comp[match(target,names(info_comp))]; tmp=tmp[!is.na(tmp)]
                if(i!=nrow(k)){
                  k$comp[(i-length(info_comp)):(i-1)]=tmp  
                }else{    # account for the different scenario in the last day
                  k$comp[(i-length(info_comp)+1):i]=tmp  
                }
                info_odd=info_oddk   # undo the change made to info_odd
              }
            }
            first_else=1   
          } # a big else 
        }  # for loop
        
        k$n=k$comp  # override n to comp
        
        # show the text when hovered over
        if(input$show_penc == TRUE){
          k$Description <- paste0("Date: ",k$StartDate,
                                  "\nCategory: ",k$Category,
                                  "\nDifference: ",percent(k$n))
        }else{
          k$Description <- paste0("Date: ",k$StartDate,
                                  "\nCategory: ",k$Category,
                                  "\nDifference: ",k$n)
        }
        title=paste0("Difference on Encounters from Two Consecutively-Recorded Days (Month ", mmonth,")")
        # k <<- k
        
      }  # the input if
      ## The end of the comparison
      
      title_style <- list(text=title,xanchor="left", yanchor="top",showarrow=F,xref = "paper",
                          yref = "paper", align = "center",x = -0.05, y = 1.15, font=list(size=16,color='black'))
      
      
      # when a month with entries only from one date, still show bars  
      if(length(unique(k$StartDate))==1) {nr=nrow(k); k[nr+1,] = k[nr,]; k$StartDate[nr+1]=k$StartDate[nr]+10; k$n[nr+1]=0}
      
      sd3 <- SharedData$new(k)
      
      # figure out the number of days in a month
      thirty_one=c(1,3,5,7,8,10,12); thirty=c(4,6,9,11); 
      if(str_sub(mmonth,6,7) %in% thirty_one) dat='-31'
      else if(str_sub(mmonth,6,7) %in% thirty) dat='-30'
      else dat='-28'
      
      pc <- sd3 %>%
        plot_ly(source = "dat_daily", name =~Category,
                x = ~StartDate, y = ~n, color=~color, type="bar",
                text=~Description, hoverinfo="text") %>% 
        layout(barmode=my_barmode, annotations=title_style ,
               yaxis=yaxi, xaxis=list(range=c(paste0(mmonth,"-01"),paste0(mmonth,dat)),title='Date', 
                                      rangeslider=list(type="date",range=c(paste0(mmonth,"-01"),paste0(mmonth,dat))), visible=TRUE))
      
      if (is.null(dat_daily())) {
        pic_daily <<- pc
        return(pic_daily)
      }
      pc
    })
    
  }) 
  
  
  ###############
  # the bar chart
  x = factor(ratio_df$group)   
  x = factor(x,levels(x)[c(8,12,5,7,10,13,16,1,15,4,6,14,3,11,2,17,9)])    #reorder factor levels
  
  output$bar <- renderPlotly({
    
    mat <- match(input$individual_id,1:(ncol(df)-2))
    # mat <- match(input$individual_id,1:nrow(dat)) 
    min_sub <- nrow(vis_df)*mat-nrow(vis_df)+1 
    max_sub <- nrow(vis_df)*mat  
    tmp_df <- subset(vis_df_all[min_sub:max_sub, ], is_highlight=="yes")
    
    # Prepare data for the barchart plot
    # ratio_df stores the information of the proportion of the phecodes above threshold in each PheWAS group.
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
    
    # Reordering the bar charts by proportion that is above threshold
    ratio_df <- ratio_df %>% mutate(color = color_vis) %>% arrange(desc(Proportion_abv_thrh)) 
    
    bar_order <<- ratio_df %>% .$Groupnum
    
    tmp <- ggplot(ratio_df, aes(x=x, y=Proportion_abv_thrh, text = description)) + 
      geom_bar(stat="identity",fill = ratio_df$color) +
      
      ggtitle ("Proportion of phecodes above threshold") +
      ylab("Proportion") +
      scale_y_continuous(expand = c(0, 0), limits = c(0,max(ratio_df$Proportion_abv_thrh)+0.05)) +
      scale_x_discrete(name = "", labels=ratio_df$group) +
      
      theme_bw() + 
      theme(
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size=8),   # family = "sans",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
        plot.title = element_text(size = 8),
        panel.background = element_rect(fill = "transparent",colour = NA),   # These two rows make the background of the barchart transparent
        plot.background = element_rect(fill = "transparent",colour = NA)) 
    
    ggplotly(tmp,tooltip="text",source="bar")
    # ggplotly(get(str_glue("ratio_id{input$individual_id}")), tooltip="text")
    
  })         
  
  
  
  # For MAIN tabset
  ####################
  # the Manhattan plot
  output$plot <- renderPlotly({
    
    phegrp_highlight <- event_data("plotly_click", source = "bar")$x
    
    mat <- match(input$individual_id,1:(ncol(df)-2))
    # mat <- match(input$individual_id,1:nrow(dat)) # `dat` contains the MAP probabilities for each individual patient across all diseases
    min_sub <- nrow(vis_df)*mat-nrow(vis_df)+1 
    max_sub <- nrow(vis_df)*mat  
    sub_df <- vis_df_all[min_sub:max_sub, ]
    
    man_df <- tibble()
    for (ele in bar_order){
      dff <- sub_df[sub_df$groupnum == ele,]
      man_df <- rbind(man_df,dff)  
    }
    
    man_df$phenotypes <- 1:nrow(man_df)
    subman_df <- subset(man_df, is_highlight=="yes")
    
    axisdf <- man_df %>% group_by(group) %>% summarize(center=( max(phenotypes) + min(phenotypes) ) / 2 )
    
    if(!is.null(phegrp_highlight)){
      subman_df <- subset(man_df,man_df$groupnum==bar_order[phegrp_highlight]) %>% subset(is_highlight=="yes")
    }
    
    tmp <- ggplot(man_df, aes(x = phenotypes, y = map_prob,text = description)) +
      
      # Show all points
      geom_point( aes(color=as.factor(groupnum)), alpha=0.25, size=1.2) +
      
      # custom X axis:
      scale_x_continuous(label = axisdf$group, breaks = axisdf$center) +
      scale_y_continuous(name = "MAP Probabilities",expand = c(0, 0),limits = c(0,1.1), breaks = c(0,0.25,0.5,0.75,1)) +     # remove space between plot area and x axis
      
      # Add highlighted points
      geom_point(data=subman_df, aes(color=as.factor(groupnum))) +
      
      scale_color_manual(values = color_vis) +
      
      # Custom the theme:
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    # make the plot interactive  
    # important to add `source="plot"` here, so that downstream plot can recall the results
    ggplotly(tmp, tooltip="text",source="plot")  
    
    # ggplotly(get(str_glue("p{input$individual_id}")), tooltip="text")
  })
  
  
  # For Info Table 
  ####################
  output$sig_tab <- renderText({
    paste("Significant Phecodes Table for Patient",input$individual_id)
  })
  
  # output the total number of phecodes that are above their corresponding threshold
  output$cond_num <- renderText({
    mat <- match(input$individual_id,1:(ncol(df)-2))
    # mat <- match(input$individual_id,1:nrow(dat)) # `dat` contains the MAP probabilities for each individual patient across all diseases
    min_sub <- nrow(vis_df)*mat-nrow(vis_df)+1 
    max_sub <- nrow(vis_df)*mat        
    
    paste("The total number of phecodes that are above their corresponding threshold is",
          subset(vis_df_all[min_sub:max_sub, ], is_highlight=="yes") %>% nrow,".")
  })
  
  
  ## !!
  # Manhattan plot selection will actively influence both info table & wordcloud 
  output$brush <- renderText({
    
    lasso <-event_data("plotly_relayout", source = "plot")
    phegrp_highlight <- event_data("plotly_click", source = "bar")$x
    
    mat <- match(input$individual_id,1:(ncol(df)-2))
    min_sub <- nrow(vis_df)*mat-nrow(vis_df)+1 
    max_sub <- nrow(vis_df)*mat  
    sub_df <- vis_df_all[min_sub:max_sub, ]
    
    man_df <- tibble()
    for (ele in bar_order){
      dff <- sub_df[sub_df$groupnum == ele,]
      man_df <- rbind(man_df,dff)  
    }
    
    man_df$phenotypes <- 1:nrow(man_df)
    subman_df <- subset(man_df, is_highlight=="yes")
    
    
    if (!is.null(lasso) & length(lasso)==4) {
      x.min <- lasso[1]; x.max <- lasso[2]
      y.min <- lasso[3]; y.max <- lasso[4]
      
      subman_df <- subman_df %>% filter(phenotypes >= x.min & phenotypes <= x.max &
                                          map_prob >= y.min & map_prob <= y.max)
    }
    
    if(!is.null(phegrp_highlight)){
      subman_df <- subset(subman_df,subman_df$groupnum==bar_order[phegrp_highlight]) 
    }
    
    num_count <- nrow(subman_df) %>% as.numeric
    
    paste0("The total number of phecodes that are both above their threshold and are selected is ", num_count,".")
    
  })
  
  
  # Show the result in the DT table                        
  output$panel <- DT::renderDT({
    
    lasso <-event_data("plotly_relayout", source = "plot")
    phegrp_highlight <- event_data("plotly_click", source = "bar")$x
    
    
    mat <- match(input$individual_id,1:(ncol(df)-2))
    min_sub <- nrow(vis_df)*mat-nrow(vis_df)+1 
    max_sub <- nrow(vis_df)*mat  
    sub_df <- vis_df_all[min_sub:max_sub, ]
    
    man_df <- tibble()
    for (ele in bar_order){
      dff <- sub_df[sub_df$groupnum == ele,]
      man_df <- rbind(man_df,dff)  
    }
    
    man_df$phenotypes <- 1:nrow(man_df)
    subman_df <- subset(man_df, is_highlight=="yes")
    
    if (!is.null(lasso) & length(lasso)==4) {
      x.min <- lasso[1]; x.max <- lasso[2]
      y.min <- lasso[3]; y.max <- lasso[4]
      
      subman_df <- subman_df %>% filter(phenotypes >= x.min & phenotypes <= x.max &
                                          map_prob >= y.min & map_prob <= y.max)
    }
    
    subman_df$phecode_pheno <- paste0(subman_df$pheno, " [", subman_df$phecode, "]")
    
    if(!is.null(phegrp_highlight)){
      subman_df <- subset(subman_df,subman_df$groupnum==bar_order[phegrp_highlight]) 
    }
    
    if(input$details == "TRUE"){
      subman_df <- data.frame(# Phecodes = sub_df$phecode, 
        Group = subman_df$group,
        cl = subman_df$groupnum,
        Phenotype = subman_df$phecode_pheno,
        MAP_prob = subman_df$map_prob %>% round(4),
        MAP_cutoff = as.numeric(subman_df$cutoff) %>% round(4))
      names(subman_df)[c(1,3,4,5)] <- c("PheWAS Group","Phenotype [Phecodes]","Map Probability","MAP Cutoff")
    } else{
      subman_df <- data.frame(# Phecodes = sub_df$phecode, 
        Group = subman_df$group,
        cl = subman_df$groupnum,
        Phenotype = subman_df$phecode_pheno,
        MAP_prob = subman_df$map_prob %>% round(4)) 
      names(subman_df)[c(1,3,4)] <- c("PheWAS Group","Phenotype [Phecodes]","Map Probability")
    }
    
    
    # Sort the table first by MAP probabilities then by category (only showing the “Yes” phenotypes)
    datatable(subman_df %>% arrange(desc(`Map Probability`),cl),
              options = list(pageLength = 20)) %>%    # each time shows only 20 rows in the output table
      formatStyle('cl',
                  backgroundColor = styleEqual(c(1:15,17:18), colvis_rgb))
    
    
  })
  
  
  ###############
  # Word Cloud  - only show what's selected in the Manhattan plot ; click thing.. not yet available here...
  observeEvent(input$individual_id, {
    
    phegrp_highlight <- event_data("plotly_click", source = "bar")$x
    
    # Try to fix a minor problem (has to move the window a little bit to make wordcloud display
    # just initially move it a little bit, then it's fine)
    if(is.null(event_data("plotly_relayout", source = "plot"))){
      
      mat <- match(input$individual_id,1:(ncol(df)-2))
      min_sub <- nrow(vis_df)*mat-nrow(vis_df)+1 
      max_sub <- nrow(vis_df)*mat  
      sub_df <- vis_df_all[min_sub:max_sub, ]
      
      man_df <- tibble()
      for (ele in bar_order){
        dff <- sub_df[sub_df$groupnum == ele,]
        man_df <- rbind(man_df,dff)  
      }
      
      man_df$phenotypes <- 1:nrow(man_df)
      subman_df <- subset(man_df, is_highlight=="yes")
      
      if(!is.null(phegrp_highlight)){
        subman_df <- subset(subman_df,subman_df$groupnum==bar_order[phegrp_highlight]) 
      }
      
      word_cl <- subman_df %>% .[,c(7,6)]  #extract columns map_prob and pheno
      colnames(word_cl) <- c("name","value")
      
      # the actual word cloud generating function
      renderWordcloud("wordcloud", data = word_cl,
                      shape = 'circle',
                      rotationRange = c(-50, 50),
                      grid_size = 5, sizeRange = c(25, 40))
    }
    
    observeEvent(event_data("plotly_relayout", source = "plot"), {
      lasso <-event_data("plotly_relayout", source = "plot")
      phegrp_highlight <- event_data("plotly_click", source = "bar")$x
      
      mat <- match(input$individual_id,1:(ncol(df)-2))
      min_sub <- nrow(vis_df)*mat-nrow(vis_df)+1 
      max_sub <- nrow(vis_df)*mat  
      sub_df <- vis_df_all[min_sub:max_sub, ]
      
      man_df <- tibble()
      for (ele in bar_order){
        dff <- sub_df[sub_df$groupnum == ele,]
        man_df <- rbind(man_df,dff)  
      }
      
      man_df$phenotypes <- 1:nrow(man_df)
      subman_df <- subset(man_df, is_highlight=="yes")
      
      if (!is.null(lasso) & length(lasso)==4) {
        x.min <- lasso[1]; x.max <- lasso[2]
        y.min <- lasso[3]; y.max <- lasso[4]
        
        subman_df <- subman_df %>% filter(phenotypes >= x.min & phenotypes <= x.max &
                                            map_prob >= y.min & map_prob <= y.max)
      }
      
      if(!is.null(phegrp_highlight)){
        subman_df <- subset(subman_df,subman_df$groupnum==bar_order[phegrp_highlight]) 
      }
      
      word_cl <- subman_df %>% .[,c(7,6)]  #extract columns map_prob and pheno
      colnames(word_cl) <- c("name","value")
      
      # the actual word cloud generating function
      renderWordcloud("wordcloud", data = word_cl,
                      shape = 'circle',
                      rotationRange = c(-90, 90),
                      grid_size = 5, sizeRange = c(25, 40))
    })
  })
  
  
  
  
  # For VD tabset
  ###############
  output$vitd <- renderPlotly({
    
    id <- match(input$patient_vd_num, choices)
    pat_encounter <- which(three_ms_vd$PatientNum == choices[id])
    
    sd1 <- SharedData$new(three_ms_vd[pat_encounter,])
    
    # Prepare: add vertical lines under each marker
    line_list <- list()
    for(i in 1:nrow(three_ms_vd[pat_encounter,])){ 
      line_color <- "#FC8D62"
      line_list[[i]] <- 
        list(type      = "line",
             fillcolor = line_color,
             line      = list(color = line_color),
             opacity   = 0.5,
             x0        = three_ms_vd[pat_encounter,]$StartDate[i],
             x1        = three_ms_vd[pat_encounter,]$StartDate[i],
             xref      = "x",
             y0        = 0, 
             y1        = three_ms_vd[pat_encounter,]$Value[i],
             yref      = "y")
    }
    
    pc <- sd1 %>%
      plot_ly(x = ~StartDate, y = ~Value, colors="#FC8D62",color="#FC8D62", #fix the `requested palette with 3 different levels` issue.
              type="scatter", mode="markers",   
              text=~description, hoverinfo="text") %>%  
      hide_legend() %>%
      layout(yaxis=list(title='Vitamin D values', visible=T), 
             xaxis=list(title='Year', rangeslider=list(type="date"), visible=T),
             shapes=line_list)  # add vertical lines under each marker
    
    pc 
    
  })
  
  observeEvent(input$enco_vd_type,{
    output$all_six <- renderPlotly({
      
      ay <- list(               # set up for the yaxis
        autotick = FALSE,
        ticks = "outside",
        tick0 = 0,
        dtick = 1,
        ticklen = 1,
        tickwidth = 1,
        tickcolor = toRGB("grey")
        #title='Encounters'
      )
      
      id <- match(input$patient_vd_num, choices)
      pat_encounter <- which(three_mss$PatientNum == choices[id])
      
      keep_category <- unique(three_mss$Category)[unique(three_mss$Category) %in% input$enco_vd_type]
      
      sd1 <-three_mss[pat_encounter,c(1,3,5,6,7,10)] %>%
        filter(Category %in% keep_category) %>%
        transform(id = as.integer(factor(Category))) %>%
        arrange(id) 
      
      pc <- sd1 %>%
        plot_ly(name =~Category,
                x = ~StartDate, y = ~Encounter, color=~color,
                text=~Description, hoverinfo="text",
                yaxis = ~paste0("y", id)) %>%
        # if you really do need explicit widths on a date axis, you can specify them as milliseconds.
        add_bars(width=1000*3600*30) %>%    # set consistent bar width
        layout(bargap = 0.05,   # set bar gap
               yaxis=ay,
               xaxis=list(title='Date',rangeslider=list(type="date", thickness=0.05), visible=T)) %>%  #add rangeslider
        subplot(nrows = 6, shareX = TRUE,
                margin = 0.03) 
      
      pc
    })
    
  })
  
  
}

shinyApp(ui, server)

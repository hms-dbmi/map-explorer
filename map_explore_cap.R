
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

# for plotting
library(shiny)
library(shinyjs)
library(DT)
library(reshape2)
library(crosstalk)



# Read in the data

# `dat` contains the MAP probabilities for each individual patient across all diseases
# setwd("~/Desktop/capstone")
load("data/MAPmanhattan.Rdata")    # dataset name: dat, 4*1864; 
load("data/MAPcutoff.Rdata")       # 1866    2
load("data/pheinfo.rda")           # 1814    5

phe_man <- colnames(dat)               # 1864, with 2 columns already filtered by Luwan; 
# the corresponding MAPcutoff for these two phecodes are NAs

# filter out phecodes with NA cutoff
MAPcutoff_filteredNA <- MAPcutoff[which(!is.na(MAPcutoff$cutoff)),]     # 1864    2
# PheWAS group info summary

groupinfo_df <- table(pheinfo$groupnum) %>% names %>% as.data.frame    # 17 groups (1-15, 17-18) 
colnames(groupinfo_df) <- "Groupnum"
groupinfo_df$num <- table(pheinfo$groupnum) %>% unname

vis_df_all <- tibble()

for (i_indi in 1:nrow(dat)){
  
  phe_man_indiv1 <- phe_man %>% tibble %>% mutate(map_prob = dat[i_indi,]) 
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
      MAPcutoff_filteredNA$phecode[i]  <- substr(MAPcutoff_filteredNA$phecode[i], 1,              nchar(MAPcutoff_filteredNA$phecode[i])-1)
    }
  }
  
  # extract list of significant phecodes for each individual patient 
  # Assume that it would be common that there would be NAs in individual patient's data (phe_man_indiv1) ; 
  # but generally fixed number of NAs in MAPcutoff file (MAPcutoff), 2 NAs in our case, already filtered for downstream analysis (MAPcutoff_filteredNA)
  sig_list <- c()
  for (i in 1:length(phe_man_indiv1$phecode)){
    if(is.na(phe_man_indiv1$map_prob[i])) next
    which_row <- match(phe_man_indiv1$phecode[i], MAPcutoff_filteredNA$phecode)
    if (phe_man_indiv1$map_prob[which_row] > MAPcutoff_filteredNA$cutoff[which_row]) {
      sig_list <- c(sig_list,phe_man_indiv1$phecode[which_row])
      # print(which_row)
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
  
  # Prepare data for the barchart plot
  # ratio_df stores the information of the proportion of the phecodes above threshold in each PheWAS group.
  tmp_df <- vis_df[which(vis_df$is_highlight == "yes"),]
  tmpp_df <- table(tmp_df$groupnum) %>% names %>% as.data.frame
  colnames(tmpp_df) <- "Groupnum"
  tmpp_df$num <- table(tmp_df$groupnum) %>% unname
  ratio_df <- groupinfo_df
  flag = 0
  for (j in 1:nrow(groupinfo_df)){
    if(tmpp_df$Groupnum[j-flag] != groupinfo_df$Groupnum[j]){
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
  
#   # Plot the barchart  - plot it in Shiny
#   x = factor(ratio_df$group)   
#   x = factor(x,levels(x)[c(8,12,5,7,10,13,16,1,15,4,6,14,3,11,2,17,9)])    #reorder factor levels
  
#   ratio_p <- ggplot(ratio_df, aes(x=x, y=Proportion_abv_thrh, text = description)) + 
#     geom_bar(stat="identity",fill = color_vis) + 
    
#     ggtitle ("Proportion of phecodes above threshold") +
#     ylab("Proportion") +
#     scale_y_continuous(expand = c(0, 0), limits = c(0,max(ratio_df$Proportion_abv_thrh)+0.05)) +
#     scale_x_discrete(name = "",labels=ratio_df$group) +
    
#     theme_bw() + 
#     theme(
#       legend.position="none",
#       panel.border = element_blank(),
#       panel.grid.major.x = element_blank(),
#       panel.grid.minor.x = element_blank(),
#       axis.title.y = element_text(family = "sans",size=8),
#       axis.text.x = element_text(family = "sans", angle = 45, hjust = 1, size = 6.5),
#       plot.title = element_text(family = "sans",size = 10)) 
  
#   assign(str_glue("ratio_id{i_indi}"),ratio_p) 
  
  axisdf <- vis_df %>% group_by(group) %>% summarize(center=( max(phenotypes) + min(phenotypes) ) / 2 )
  
#   ## make the manhattan plot  - make it in Shiny
#   p <- ggplot(vis_df, aes(x = phenotypes, y = map_prob,text = description)) +
    
#     # Show all points
#     geom_point( aes(color=as.factor(groupnum)), alpha=0.25, size=1.2) +
    
#     # custom X axis:
#     scale_x_continuous(label = axisdf$group, breaks = axisdf$center) +
#     scale_y_continuous(name = "MAP Probabilities",expand = c(0, 0),limits = c(0,1.1), breaks = c(0,0.25,0.5,0.75,1)) +     # remove space between plot area and x axis
    
#     # Add highlighted points
#     geom_point(data=subset(vis_df, is_highlight=="yes"), aes(color=as.factor(groupnum))) +
    
#     scale_color_manual(values = color_vis) +
    
#     # Custom the theme:
#     theme_bw() +
#     theme( 
#       legend.position="none",
#       panel.border = element_blank(),
#       panel.grid.major.x = element_blank(),
#       panel.grid.minor.x = element_blank(),
#       axis.text.x = element_text(angle = 45, hjust = 1)
#     )
  
#   assign(str_glue("p{i_indi}"),p) 
  
  # ggplotly() function returns a plotly object
  # ggplotly(p, tooltip="text")
  
  vis_df <- merge(vis_df, MAPcutoff_filteredNA, by.x = "phenotype", by.y = "phecode") %>% arrange(phenotypes) # add cutoff threshold info into vis_df
  colnames(vis_df)[1] <- "phecode"
  ### append into a final tibble that contains all the info
  vis_df_all <- rbind(vis_df_all,vis_df)   
  print(dim(vis_df_all))
  
}

#############
# for MS data 
# read in the MS patient files
three_ms <- read_csv("data/3_MS_patients.csv")  #,sep=",", header = T, stringsAsFactors = F)

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

three_ms_vd$description <- paste0("Patient Number: ",three_ms_vd$PatientNum,
                                  "\nStart Date: ",three_ms_vd$StartDate, 
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

# pat_encounter_1 <- which(three_mss$PatientNum == "68286")
# pat_encounter_2 <- which(three_mss$PatientNum == "99492")
# pat_encounter_3 <- which(three_mss$PatientNum == "106579")
choices <- c(68286, 99492, 106579)

###########
## Capstone project
# shiny app


ui <- fluidPage(
  titlePanel("MAP Explorer"),
  
  sidebarPanel(
    h4("You can explore the MAP data in this Shiny App!"),    
    selectInput(inputId="individual_id",
                label="Select Patient ID: ",
                choices=c(1:4), 
                selected = 1, selectize = F),
    plotlyOutput("bar"),
    
    h3("How to use"),
    p("The default visualization shows the results for Patient 1.
      You can also select any Patient ID from the dropdown list.
      Then it would output a barchart showing the information of the proportion of the phecodes above threshold in each PheWAS group of this individual patient.
      It would also output a Manhattan plot showing the MAP probabilities of each phecodes of this selected patient.
      You can hover over the points in Manhattan plot or Barchart to get more information.
      Enjoy playing around with it!")
    ),
  
  mainPanel(
    useShinyjs(),
    tabsetPanel(
      tabPanel("Main",
               h3("MAP Manhattan Plot"),
               fluidRow(
                 plotlyOutput("plot",height = 500),
                 br(),
                 h3("Word Cloud of Phenotypes based on MAP Probabilities"),
                 # for Word Cloud
                 # MUST load the ECharts javascript library in advance
                 loadEChartsLibrary(),
                 tags$div(id="wordcloud", style="width:100%;height:500px;"),
                 deliverChart(div_id = "wordcloud"))            
               
      ),
      
      tabPanel("Info Table",
               br(),
               textOutput("sig_tab"),
               br(),
               textOutput("cond_num"),     #report number of phecodes above threshold
               textOutput("brush"),        #report number of phecodes both above threshold and selected
               br(),
               radioButtons("details", "More details (with MAP cutoff):", choices=c("Yes","No"), selected = "No"),
               DT::DTOutput("panel"),
               tags$head(tags$style("#sig_tab{color: black; font-size: 20px;}")),
               tags$head(tags$style("#cond_num{color: black; font-size: 15px; font-style:italic;}")),  
               tags$head(tags$style("#brush{color: black; font-size: 14px; font-style:italic;}"))
               
      ),
      tabPanel("MS",
               h3("MS Data Overview"),
               selectInput(inputId="patient_num",
                           label="Select Patient Number: ",
                           choices=c(68286, 99492, 106579), 
                           selected = 1,
                           width = "30%"), # modify the size of the input box
               br(),
               plotlyOutput("dat_year", height = 300),
               br(),
               plotlyOutput("dat_month", height = 300),
               br(),
               plotlyOutput("dat_daily", height = 300)
      ),
      tabPanel("VD",
               h3("Vitamin D Levels"),
               selectInput(inputId="patient_vd_num",
                           label="Select Patient Number: ",
                           choices=c(68286, 99492, 106579), 
                           selected = 1,
                           width = "30%"), # modify the size of the input box
               br(),
               plotlyOutput("vitd"),
               br(),
               plotlyOutput("all_six", height = 600))
    ) # for tabsetPanel
    
  )) # ")" for mainPanel & fluidPage


server <- function(input, output, session) {
  
  # For MS tabset
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
  
  ####note: try to highlight only the clicked/selected bar; why with 'if' statement, it works.
  # there are some issue in 'if' statement, if keep clicking on the barchart on the same year, there would be mistakes
  
  output$dat_year <- renderPlotly({
    
    id <- match(input$patient_num, choices)
    pat_encounter <- which(three_mss$PatientNum == choices[id])
    
    if (!is.null(dat_year())) {
      selected_df <- three_mss[pat_encounter,] %>% mutate(opacity = ifelse(dat_year() == Year,1,0.2))
    } else{
      selected_df <- three_mss[pat_encounter,] %>% mutate(opacity = 1)
    }
    
    sd1 <- SharedData$new(selected_df)
    # sd1 <- SharedData$new(three_mss[pat_encounter,])
    
    pc <- sd1 %>%
      plot_ly(source = "dat_year", name =~Category, # name of the legend
              x = ~Year, y = ~Encounter, color=~color, type="bar",    # ensure each Category has unique color
              text=~Description, hoverinfo="text") %>%  #,opacity=~opacity
      #add_bars(x = ~Year, y = ~Encounter, color=~color) %>% # ensure each Category has unique color
      layout(barmode='stack', title = "Encounters Aggregated by Year",
             yaxis=list(title='Encounters', visible=T), 
             xaxis=list(title='Year', rangeslider=list(type="date"), visible=T)) 
    
    
    if (is.null(dat_year())) {
      pic_year <<- pc
      return(pic_year) 
    } #else{
    #   
    #   pc <- sd1 %>%
    #     plot_ly(source = "dat_year", name =~Category, # name of the legend
    #             x = ~Year, y = ~Encounter, color=~color, type="bar",
    #             text=~Description, hoverinfo="text", marker = list(opacity=ifelse(dat_year() == ~Year,1,0.6))) %>% 
    #     
    #     # add_bars(x = ~Year, y = ~Encounter, color=~color) %>% # ensure each Category has unique color
    #     layout(barmode='stack', title = "Encounters Aggregated by Year",
    #            yaxis=list(title='Encounters', visible=T), 
    #            xaxis=list(title='Year', rangeslider=list(type="date"), visible=T))
    # }
    
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
    
    sd <- three_mss[pat_encounter,] %>% 
      filter(Year == dat_year())
    yyear <- sd$Year[1] %>% substr(1, 4)   # which year is clicked
    sd2 <- SharedData$new(sd)
    pc <- sd2 %>%
      plot_ly(source = "dat_month", name =~Category,
              text=~Description, hoverinfo="text") %>% 
      suppressWarnings %>%
      add_bars(x = ~Month, y = ~Encounter, color=~color) %>%
      layout(barmode='stack', title = paste0("Encounters Aggregated by Month (Year ", yyear,")"),
             yaxis=list(title='Encounters', visible=TRUE), xaxis=list(title='Month', rangeslider=list(type="date"), visible=TRUE))
    
    if (is.null(dat_month())) {
      pic_month <<- pc
      return(pic_month)
    }
    pc
    # bscols(
    #   widths=c(3,NA),
    #   list(
    #     filter_checkbox('category', 'Category', sd2, ~Category, inline=F)
    #   ),
    #   pc)
  })
  
  output$dat_daily <- renderPlotly({
    if (is.null(dat_month())) return(NULL)
    
    id <- match(input$patient_num, choices)
    pat_encounter <- which(three_mss$PatientNum == choices[id])
    
    sd <- three_mss[pat_encounter,] %>% 
      filter(Month == dat_month())
    mmonth <- sd$Month[1] %>% substr(1, 7)
    sd3 <- SharedData$new(sd)
    pc <- sd3 %>%
      plot_ly(source = "dat_daily", name =~Category,
              text=~Description, hoverinfo="text") %>% 
      add_bars(x = ~StartDate, y = ~Encounter, color=~color) %>%
      layout(barmode='stack', title = paste0("Encounters Aggregated by Day (Month ", mmonth,")"),
             yaxis=list(title='Encounters', visible=TRUE), xaxis=list(title='Date', rangeslider=list(type="date"), visible=TRUE))
    
    if (is.null(dat_daily())) {
      pic_daily <<- pc
      return(pic_daily)
    }
    pc
    # bscols(
    #   widths=c(3,NA),
    #   list(
    #     filter_checkbox('category', 'Category', sd3, ~Category, inline=F)
    #   ),
    #   pc)
  })
  
  
  # For MAIN tabset
  ####################
  # the Manhattan plot
  output$plot <- renderPlotly({
    
    mat <- match(input$individual_id,1:nrow(dat)) # `dat` contains the MAP probabilities for each individual patient across all diseases
    min_sub <- nrow(vis_df)*mat-nrow(vis_df)+1 
    max_sub <- nrow(vis_df)*mat  
    sub_df <- vis_df_all[min_sub:max_sub, ]
    
    
    tmp <- ggplot(sub_df, aes(x = phenotypes, y = map_prob,text = description)) +
      
      # Show all points
      geom_point( aes(color=as.factor(groupnum)), alpha=0.25, size=1.2) +
      
      # custom X axis:
      scale_x_continuous(label = axisdf$group, breaks = axisdf$center) +
      scale_y_continuous(name = "MAP Probabilities",expand = c(0, 0),limits = c(0,1.1), breaks = c(0,0.25,0.5,0.75,1)) +     # remove space between plot area and x axis
      
      # Add highlighted points
      geom_point(data=subset(sub_df, is_highlight=="yes"), aes(color=as.factor(groupnum))) +
      
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
    ggplotly(tmp, tooltip="text")
    # ggplotly(get(str_glue("p{input$individual_id}")), tooltip="text")
  })
  
  
  ###############
  # Word Cloud
  observeEvent(input$individual_id, {
    mat <- match(input$individual_id,1:nrow(dat)) # `dat` contains the MAP probabilities for each individual patient across all diseases
    min_sub <- nrow(vis_df)*mat-nrow(vis_df)+1 
    max_sub <- nrow(vis_df)*mat  
    word_cl <- vis_df_all[min_sub:max_sub, ] %>% 
               filter(is_highlight == "yes") %>% 
               .[,c(7,6)]  #extract columns map_prob and pheno
    colnames(word_cl) <- c("name","value")
    
    renderWordcloud("wordcloud", data =word_cl,
                    shape = 'circle',
                    rotationRange = c(-90, 90),
                    grid_size = 5, sizeRange = c(25, 40))
  })
  
  
  ###############
  # the bar chart
  x = factor(ratio_df$group)   
  x = factor(x,levels(x)[c(8,12,5,7,10,13,16,1,15,4,6,14,3,11,2,17,9)])    #reorder factor levels
  
  output$bar <- renderPlotly({
    tmp <- ggplot(ratio_df, aes(x=x, y=Proportion_abv_thrh, text = description)) + 
      geom_bar(stat="identity",fill = color_vis) + 
      
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
        axis.title.y = element_text(family = "sans",size=8),
        axis.text.x = element_text(family = "sans", angle = 45, hjust = 1, size = 6.5),
        plot.title = element_text(family = "sans",size = 10)) 
    
    ggplotly(tmp,tooltip="text")
    # ggplotly(get(str_glue("ratio_id{input$individual_id}")), tooltip="text")
    
  })
  
  # For Info Table tabset
  ####################
  output$sig_tab <- renderText({
    paste("Significant Phecodes Table for Patient",input$individual_id)
  })
  
  # output the total number of phecodes that are above their corresponding threshold
  output$cond_num <- renderText({
    
    mat <- match(input$individual_id,1:nrow(dat)) # `dat` contains the MAP probabilities for each individual patient across all diseases
    min_sub <- nrow(vis_df)*mat-nrow(vis_df)+1 
    max_sub <- nrow(vis_df)*mat        
    
    paste("The total number of phecodes that are above their corresponding threshold is",
          subset(vis_df_all[min_sub:max_sub, ], is_highlight=="yes") %>% nrow,".")
  })
  
  
  
  # Show the result in the DT table                        
  output$panel <- DT::renderDT({
    
    lasso <- event_data("plotly_relayout")
    
    mat <- match(input$individual_id,1:nrow(dat)) # `dat` contains the MAP probabilities for each individual patient across all diseases
    min_sub <- nrow(vis_df)*mat-nrow(vis_df)+1 
    max_sub <- nrow(vis_df)*mat  
    sub_df <- subset(vis_df_all[min_sub:max_sub, ], is_highlight=="yes")
    
    if (!is.null(lasso) & length(lasso)==4) {
      x.min <- lasso[1]; x.max <- lasso[2]
      y.min <- lasso[3]; y.max <- lasso[4]
      
      sub_df <- sub_df %>% filter(phenotypes >= x.min & phenotypes <= x.max &
                                    map_prob >= y.min & map_prob <= y.max)
    }
    
    if(input$details == "Yes"){
      sub_df <- data.frame(Phecodes = sub_df$phecode, 
                           Group = sub_df$group,
                           cl = sub_df$groupnum,
                           Phenotype = sub_df$pheno,
                           MAP_prob = sub_df$map_prob %>% round(4),
                           MAP_cutoff = sub_df$cutoff %>% round(4))
    } else{
      sub_df <- data.frame(Phecodes = sub_df$phecode, 
                           Group = sub_df$group,
                           cl = sub_df$groupnum,
                           Phenotype = sub_df$pheno,
                           MAP_prob = sub_df$map_prob %>% round(4)) 
    }
    
    # Sort the table first by category then by MAP probabilities (only showing the “Yes” phenotypes)
    datatable(sub_df %>% arrange(cl,desc(MAP_prob)),
              options = list(pageLength = 10)) %>%    # each time shows only 10 rows in the output table
      formatStyle('cl',
                  backgroundColor = styleEqual(c(1:15,17:18), colvis_rgb))
    
    
  })
  
  output$brush <- renderText({
    lasso <- event_data("plotly_relayout")
    
    mat <- match(input$individual_id,1:nrow(dat)) # `dat` contains the MAP probabilities for each individual patient across all diseases
    min_sub <- nrow(vis_df)*mat-nrow(vis_df)+1 
    max_sub <- nrow(vis_df)*mat  
    sub_df <- subset(vis_df_all[min_sub:max_sub, ], is_highlight=="yes")  # create a subset for each individual patient
    
    if (!is.null(lasso) & length(lasso)==4) {
      x.min <- lasso[1]; x.max <- lasso[2]
      y.min <- lasso[3]; y.max <- lasso[4]
      
      sub_df <- sub_df %>% filter(phenotypes >= x.min & phenotypes <= x.max &
                                    map_prob >= y.min & map_prob <= y.max)
    }
    
    num_count <- nrow(sub_df) %>% as.numeric
    
    paste0("The total number of phecodes that are both above their threshold and are selected is ", num_count,".")
    
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
      line_color <- "#FFD92F"
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
      plot_ly(x = ~StartDate, y = ~Value, colors="#FFD92F",color="#FFD92F", #fix the `requested palette with 3 different levels` issue.
              type="scatter", mode="markers",   
              text=~description, hoverinfo="text") %>%  
      hide_legend() %>%
      layout(yaxis=list(title='Vitamin D values', visible=T), 
             xaxis=list(title='Year', rangeslider=list(type="date"), visible=T),
             shapes=line_list)  # add vertical lines under each marker
    
    pc 
    
  })
  
  
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
    
    
    sd1 <-three_mss[pat_encounter,c(1,3,5,6,7,10)] %>%
      transform(id = as.integer(factor(Category))) %>%
      arrange(id) #%>% mutate(width=0.01)
    
    pc <- sd1 %>%
      plot_ly(name =~Category,
              x = ~StartDate, y = ~Encounter, color=~color,
              text=~Description, hoverinfo="text",
              yaxis = ~paste0("y", id)) %>%
      # if you really do need explicit widths on a date axis, you can specify them as milliseconds.
      add_bars(width=1000*3600*30) %>%    # set consistent bar width
      layout(title = "Encounters Aggregated by Year",
             bargap = 0.05,   # set bar gap
             yaxis=ay,
             xaxis=list(title='Date',rangeslider=list(type="date", thickness=0.05), visible=T)) %>%  #add rangeslider
      subplot(nrows = 6, shareX = TRUE,
              margin = 0.03) 
    
    pc
  })
  
}

shinyApp(ui, server)

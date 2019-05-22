
# Library required packages

# for data processing and analysis
# devtools::install_github("PheWAS/PheWAS")
library(qqman)
library(PheWAS)
library(stringr)
library(plotly)
library(purrr)
library(scales)

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



###########
## Capstone project
# shiny app


library(shiny)
library(plotly)
library(shinyjs)

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
                 plotlyOutput("plot"),
                 br(),
                 textOutput("sig_tab"),
                 br(),
                 textOutput("cond_num"),     #report number of phecodes above threshold
                 textOutput("brush"),        #report number of phecodes both above threshold and selected
                 br(),
                 DT::DTOutput("panel") ,
                 tags$head(tags$style("#sig_tab{color: black; font-size: 20px;}")),
                 tags$head(tags$style("#cond_num{color: black; font-size: 15px; font-style:italic;}")),  
                 tags$head(tags$style("#brush{color: black; font-size: 14px; font-style:italic;}")))
      ),
      
      tabPanel("MS",
               h3("MS Data Overview"),
               plotlyOutput("ms_bar"),
               sliderInput(
                 "stdat", "Start date:",
                 min=min(three_ms$StartDate),
                 max=max(three_ms$StartDate),
                 value=c(min(three_ms$StartDate), max(three_ms$StartDate))
               ))
      
    )
    
  )) # ")" for mainPanel & fluidPage


server <- function(input, output) {
  
  # renderPlotly() also understands ggplot2 objects!
  
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
    
    sub_df <- data.frame(Phecodes = sub_df$phecode, 
                         Group = sub_df$group,
                         cl = sub_df$groupnum,
                         Phenotype = sub_df$pheno,
                         MAP_prob = sub_df$map_prob %>% round(4),
                         MAP_cutoff = sub_df$cutoff %>% round(4))
    
    datatable(sub_df,options = list(pageLength =5)) %>%    # each time shows only 5 rows in the output table
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
  
  
  # the second tabset
  output$ms_bar <- renderPlotly({
    mat <- match(input$individual_id,1:nrow(dat))
    if(mat == '1') {
      
    }
    p1 <- ggplot(three_ms[pat_encounter_1,], 
                 aes(x=StartDate, y= encounter, fill=Category,text = description)) + 
      geom_bar(stat="identity", position = "stack") +
      
      scale_x_date(date_breaks = "3 months") + 
      
      # Custom the theme:
      theme_bw() +
      theme(legend.position="none",
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
    
    ggplotly(p1,tooltip="text")
  })
  
  
}

shinyApp(ui, server)

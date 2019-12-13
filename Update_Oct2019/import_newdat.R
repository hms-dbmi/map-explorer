
setwd("~/Desktop/Capstone")

# Read in the encounter table with or without relapse
new_enco = read_csv('data/data Oct 2019/select_encounters_withrelapse.csv')
new_enco$StartDate <- mdy(new_enco$StartDate)
# dim(new_enco)
# [1] 815 143
unique(new_enco$PatientNum)
# [1]  11634  75468  98204  31860 136178

new_enco_noRelapse = read_csv('data/data Oct 2019/select_patients_norelapse.csv')
# new_enco_noRelapse$StartDate <- mdy(new_enco_noRelapse$StartDate)   # the data has already been correctly formatted
unique(new_enco_noRelapse$PatientNum)


# Read in the vitamin D table
new_vitd=read.table('data/data Oct 2019/25_OH_Vit_D_result_D2_excluded_20171018.csv',sep = "|")

# The patients in the Vitamin D file are an exhaustive list
# Check whether all of the five patients in the Encounter file are in the Vitamin D file. Yes!
all(unique(new_enco$PatientNum) %in% new_vitd$V1)  #TRUE

all(unique(new_enco_noRelapse$PatientNum) %in% new_vitd$V1)  #TRUE



# read in the lasso glmnet coefficient & description table
# BY Beta
lasso_coeff = read_csv('data/data Oct 2019/cv_glmnet_lasso.csv') %>% arrange(desc(beta))


### Solu-Medrol is a corticosteroid with anti-inflammatory action. It is used to treat acute exacerbations in patients with multiple sclerosis (MS). 

### Lesion of brain
### In multiple sclerosis (MS), the body mistakenly attacks the protective layer around the nerves in the brain and spinal cord (also known as myelin). 
### These damaged areas are called plaques or lesions. Everyone with MS will get lesions with varying severity. 
### However, lesions tend to happen more in people with !!!!!relapsing MS. Healthcare providers monitor lesions to track disease progression.


#  !!
lasso_coeff[c(1,4:5,7,8,10:12,14:15,18,121),]
c(1,4:5,7,8:15)      # interested i
# [1]  1  4  5  7  8  9 10 11 12 13 14 15

tp=lasso_coeff[c(1,4:5,7,8,10:12,14:15,18,121),] %>% .[,1] %>% unlist %>% unname
table(colnames(new_enco_noRelapse) %in% tp)   # all selected encounter types are in the new patient data with no relapse
# FALSE  TRUE 
# 194    12
# length(tp) #12

as.character(lasso_coeff[i,1])  # the i-th CUIs/CPT groups/Phecodes
as.character(lasso_coeff[i,3])  # the i-th description of the corresponding CUIs/CPT groups/Phecodes

# inspect the data
# new_enco[,colnames(new_enco)==as.character(lasso_coeff[1,1])] %>% head

five_mss_withRl <- new_enco[,c(2,4,which(colnames(new_enco) %in% as.character(unlist(lasso_coeff[c(1,4:5,7,8,10:12,14:15,18,121),1]))))] 
five_mss_withRl <- gather(five_mss_withRl,'Category','Encounter',-PatientNum,-StartDate)


seq_Category <- left_join(five_mss_withRl,lasso_coeff,by=c('Category'="term")) %>% arrange(desc(beta))


ori='Recurrent disease'; ord=1  # ord c(1:12)
for(i in 1:nrow(seq_Category)){
  if(seq_Category$desc[i] != ori){
    ori=seq_Category$desc[i]
    ord=ord+1
    seq_Category$desc[i] <- paste0(ord,' | ',seq_Category$desc[i])
  }else{
    seq_Category$desc[i] <- paste0(ord,' | ',seq_Category$desc[i])
  }
}

order_Category <- unique(seq_Category$desc)


five_mss_withRl <- seq_Category %>% arrange(StartDate,desc(beta)) %>% mutate(Category=desc)  %>% .[,-5]  %>%  # get rid of the actual CUIs/CPT groups/Phecodes
                        filter(Encounter==0)    # remove those rows with only 0 encounter
 
# five_mss_withRl <- five_mss_withRl %>% mutate(Category=as.character(lasso_coeff[1,3]),
#                                               Code=as.character(lasso_coeff[1,1]))
# colnames(five_mss_withRl)[3] <- 'Encounter'

five_mss_withRl$Description <- paste0("Start Date: ",five_mss_withRl$StartDate, 
                                      "\nPatient Number: ",five_mss_withRl$PatientNum,
                                      "\nCategory: ",five_mss_withRl$Category,
                                      "\nDescription: ",five_mss_withRl$Description)


# aggregate the encounters by month
five_mss_withRl <- five_mss_withRl %>% group_by(Month=floor_date(StartDate, "month"))
five_mss_withRl <- five_mss_withRl %>% group_by(Year=floor_date(StartDate, "year"))


# five_mss_withRl$color <- 'orange'
five_mss_withRl$color <- factor(five_mss_withRl$Category, labels = RColorBrewer::brewer.pal(length(unique(five_mss_withRl$Category)), name = "Set3"))


################################
#### For patient with no relapse

# lasso_coeff[c(1,4:5,7,8,10:12,14:15,18,121),1]

nfive_mss_withRl <- new_enco_noRelapse[,c(2,4,which(colnames(new_enco_noRelapse) %in% as.character(unlist(lasso_coeff[c(1,4:5,7,8,10:12,14:15,18,121),1]))))] 
nfive_mss_withRl <- gather(nfive_mss_withRl,'Category','Encounter',-PatientNum,-StartDate)


seq_Category <- left_join(nfive_mss_withRl,lasso_coeff,by=c('Category'="term")) %>% arrange(desc(beta))


ori='Recurrent disease'; ord=1  # ord c(1:12)
for(i in 1:nrow(seq_Category)){
  if(seq_Category$desc[i] != ori){
    ori=seq_Category$desc[i]
    ord=ord+1
    seq_Category$desc[i] <- paste0(ord,' | ',seq_Category$desc[i])
  }else{
    seq_Category$desc[i] <- paste0(ord,' | ',seq_Category$desc[i])
  }
}

n_order_Category <- unique(seq_Category$desc)


nfive_mss_withRl <- seq_Category %>% arrange(StartDate,desc(beta)) %>% mutate(Category=desc)  %>% .[,-5]  %>%  # get rid of the actual CUIs/CPT groups/Phecodes
  filter(Encounter==0)    # remove those rows with only 0 encounter

# five_mss_withRl <- five_mss_withRl %>% mutate(Category=as.character(lasso_coeff[1,3]),
#                                               Code=as.character(lasso_coeff[1,1]))
# colnames(five_mss_withRl)[3] <- 'Encounter'

nfive_mss_withRl$Description <- paste0("Start Date: ",nfive_mss_withRl$StartDate, 
                                      "\nPatient Number: ",nfive_mss_withRl$PatientNum,
                                      "\nCategory: ",nfive_mss_withRl$Category,
                                      "\nDescription: ",nfive_mss_withRl$Description)


# aggregate the encounters by month
nfive_mss_withRl <- nfive_mss_withRl %>% group_by(Month=floor_date(StartDate, "month"))
nfive_mss_withRl <- nfive_mss_withRl %>% group_by(Year=floor_date(StartDate, "year"))


# five_mss_withRl$color <- 'orange'
nfive_mss_withRl$color <- factor(nfive_mss_withRl$Category, labels = RColorBrewer::brewer.pal(length(unique(nfive_mss_withRl$Category)), name = "Set3"))

## we need `nfive_mss_withRl`, and `n_order_Category` in the web app




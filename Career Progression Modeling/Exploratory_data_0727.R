############################################
## Initialize the environment, load the data
## setwd, list.files, load data
############################################
getwd()
setwd('/usr/local/DiscoveryData')
list.files()

vec<-(list.files())
sapply(vec, file.info)

getwd()
list.files(pattern = "BD*")

## read data and create different datasets
BD_2013_06 <- read.csv("BD_Rep_2013_06_06-29-2013.csv", sep = ",", header = TRUE)
BD_2013_12 <- read.csv("BD_Rep_2013_12_12-28-2013.csv", sep = ",", header = TRUE)
BD_2014_06 <- read.csv("BD_Rep_2014_06_06-28-2014.csv", sep = ",", header = TRUE)
BD_2014_12 <- read.csv("BD_Rep_2014_12_01-02-2015.csv", sep = ",", header = TRUE)
BD_2015_06 <- read.csv("BD_Rep_2015_06_06-13-2015.csv", sep = ",", header = TRUE)
BD_2015_12 <- read.csv("BD_Rep_2015_12_12-26-2015.csv", sep = ",", header = TRUE)
BD_2016_06 <- read.csv("BD_Rep_2016_06_06-25-2016.csv", sep = ",", header = TRUE)
BD_2016_12 <- read.csv("BD_Rep_2016_12_12-31-2016.csv", sep = ",", header = TRUE)
BD_2017_06 <- read.csv("BD_Rep_2017_06_06-24-2017.csv", sep = ",", header = TRUE)
BD_2017_12 <- read.csv("BD_Rep_2017_12_12-30-2017.csv", sep = ",", header = TRUE)
BD_2018_06 <- read.csv("BD_Rep_2018_06_06-30-2018.csv", sep = ",", header = TRUE)
BD_2018_12 <- read.csv("BD_Rep_2018_12_12-29-2018.csv", sep = ",", header = TRUE)
BD_2019_06 <- read.csv("BD_Rep_2019_06_06-29-2019.csv", sep = ",", header = TRUE)
BD_2019_12 <- read.csv("BD_Rep_2019_12_12-28-2019.csv", sep = ",", header = TRUE)
BD_2020_06 <- read.csv("BD_Rep_2020_06_06-27-2020.csv", sep = ",", header = TRUE)
BD_2020_12 <- read.csv("BD_Rep_2020_12_12-26-2020.csv", sep = ",", header = TRUE)
BD_2021_06 <- read.csv("BD_Rep_2021_06_06-26-2021.csv", sep = ",", header = TRUE)
#after reading the data, we want to do descriptive analysis
dim(BD_2013_06);dim(BD_2013_12);dim(BD_2014_06)
##566107 190; 570786 190; 569532 190
##we can find that those datasets have different number of rows but same number of columns
##next, we want to see str of dataset BD_2013_06
str(BD_2013_06)
class(BD_2013_06)
##column names
colnames_1306 <- colnames(BD_2013_06)
colnames_1312 <- colnames(BD_2013_12)
##now we want to see whether all columns names are same
sum(colnames_1306 == colnames_1312) != 0 ##this shows all columns names are same, we can just focus on one colnames
rm(colnames_1312)
##find BDFirmName ofBD_2013_06
BD_FirmN_1306 <- BD_2013_06$BDFirmName
View(table(BD_FirmN_1306))
## choose "Merrill Lynch Professional Clearing Corp." and " Merrill Lynch, Pierce, Fenner & Smith Incorporate"
library(dplyr)
ML_2013_06 <- BD_2013_06 %>%
                      filter(BDFirmName == "Merrill Lynch Professional Clearing Corp."| BDFirmName == "Merrill Lynch, Pierce, Fenner & Smith Incorporated" )
rm(BD_2013_06)
ML_2013_12 <- BD_2013_12 %>%
  filter(BDFirmName == "Merrill Lynch Professional Clearing Corp."| BDFirmName == "Merrill Lynch, Pierce, Fenner & Smith Incorporated" )
rm(BD_2013_12)
ML_2014_06 <- BD_2014_06 %>%
  filter(BDFirmName == "Merrill Lynch Professional Clearing Corp."| BDFirmName == "Merrill Lynch, Pierce, Fenner & Smith Incorporated" )
rm(BD_2014_06)
ML_2014_12 <- BD_2014_12 %>%
  filter(BDFirmName == "Merrill Lynch Professional Clearing Corp."| BDFirmName == "Merrill Lynch, Pierce, Fenner & Smith Incorporated" )
rm(BD_2014_12)
ML_2015_06 <- BD_2015_06 %>%
  filter(BDFirmName == "Merrill Lynch Professional Clearing Corp."| BDFirmName == "Merrill Lynch, Pierce, Fenner & Smith Incorporated" )
rm(BD_2015_06)
ML_2015_12 <- BD_2015_12 %>%
  filter(BDFirmName == "Merrill Lynch Professional Clearing Corp."| BDFirmName == "Merrill Lynch, Pierce, Fenner & Smith Incorporated" )
rm(BD_2015_12)
ML_2016_06 <- BD_2016_06 %>%
  filter(BDFirmName == "Merrill Lynch Professional Clearing Corp."| BDFirmName == "Merrill Lynch, Pierce, Fenner & Smith Incorporated" )
rm(BD_2016_06)
ML_2016_12 <- BD_2016_12 %>%
  filter(BDFirmName == "Merrill Lynch Professional Clearing Corp."| BDFirmName == "Merrill Lynch, Pierce, Fenner & Smith Incorporated" )
rm(BD_2016_12)
ML_2017_06 <- BD_2017_06 %>%
  filter(BDFirmName == "Merrill Lynch Professional Clearing Corp."| BDFirmName == "Merrill Lynch, Pierce, Fenner & Smith Incorporated" )
rm(BD_2017_06)
ML_2017_12 <- BD_2017_12 %>%
  filter(BDFirmName == "Merrill Lynch Professional Clearing Corp."| BDFirmName == "Merrill Lynch, Pierce, Fenner & Smith Incorporated" )
rm(BD_2017_12)
ML_2018_06 <- BD_2018_06 %>%
  filter(BDFirmName == "Merrill Lynch Professional Clearing Corp."| BDFirmName == "Merrill Lynch, Pierce, Fenner & Smith Incorporated" )
rm(BD_2018_06)
ML_2018_12 <- BD_2018_12 %>%
  filter(BDFirmName == "Merrill Lynch Professional Clearing Corp."| BDFirmName == "Merrill Lynch, Pierce, Fenner & Smith Incorporated" )
rm(BD_2018_12)
ML_2019_06 <- BD_2019_06 %>%
  filter(BDFirmName == "Merrill Lynch Professional Clearing Corp."| BDFirmName == "Merrill Lynch, Pierce, Fenner & Smith Incorporated" )
rm(BD_2019_06)
ML_2019_12 <- BD_2019_12 %>%
  filter(BDFirmName == "Merrill Lynch Professional Clearing Corp."| BDFirmName == "Merrill Lynch, Pierce, Fenner & Smith Incorporated" )
rm(BD_2019_12)
ML_2020_06 <- BD_2020_06 %>%
  filter(BDFirmName == "Merrill Lynch Professional Clearing Corp."| BDFirmName == "Merrill Lynch, Pierce, Fenner & Smith Incorporated" )
rm(BD_2020_06)
ML_2020_12 <- BD_2020_12 %>%
  filter(BDFirmName == "Merrill Lynch Professional Clearing Corp."| BDFirmName == "Merrill Lynch, Pierce, Fenner & Smith Incorporated" )
rm(BD_2020_12)
ML_2021_06 <- BD_2021_06 %>%
  filter(BDFirmName == "Merrill Lynch Professional Clearing Corp."| BDFirmName == "Merrill Lynch, Pierce, Fenner & Smith Incorporated" )
rm(BD_2021_06)
ML_datalst <- c(ML_2013_06, ML_2013_12, ML_2014_06, ML_2014_12, ML_2015_06, ML_2015_12, ML_2016_06, ML_2016_12, ML_2017_06, ML_2017_12, ML_2018_06, ML_2018_12, ML_2019_06, ML_2019_12,
                ML_2020_06, ML_2020_12, ML_2021_06)
##we want to check for any nulls in Discovery ID
l = length(ML_datalst)
na_c <- c()
for (i in 1:l){
    append(na_c, sum(is.nan(ML_datalst[i]$DiscoveryContactID)))
}
print(sum(na_c))
##this shows there is no null in Discovery ID
##now check whether there is any duplicate Discovery ID in each dataset
dup_c <- c()
for (i in 1:l){
  append(dup_c, sum(duplicated(ML_datalst[i]$DiscoveryContactID)))
}
print(sum(dup_c)) ## this shows for each dataset, DiscoveryContactID is unique and can be used as key
##now I want to merge those datasets from different year by DiscoveryContactID
ML_2013 <- merge.data.frame(ML_2013_06, ML_2013_12, by = "DiscoveryContactID")
ML_2014 <- merge.data.frame(ML_2014_06, ML_2014_12, by = "DiscoveryContactID")
ML_2015 <- merge.data.frame(ML_2015_06, ML_2015_12, by = "DiscoveryContactID")
ML_2016 <- merge.data.frame(ML_2016_06, ML_2016_12, by = "DiscoveryContactID")
ML_2017 <- merge.data.frame(ML_2017_06, ML_2017_12, by = "DiscoveryContactID")
ML_2018 <- merge.data.frame(ML_2018_06, ML_2018_12, by = "DiscoveryContactID")
ML_2019 <- merge.data.frame(ML_2019_06, ML_2019_12, by = "DiscoveryContactID")
ML_2020 <- merge.data.frame(ML_2020_06, ML_2020_12, by = "DiscoveryContactID")
dim(ML_2013) ## 28496, 379
dim(ML_2013_06)[1] ; dim(ML_2013_12)[1] ## 30195, 190 ; 30077, 190 ## there are some data without same ID
View(table(ML_2013$Title.x))
View(table(ML_2013$Title.y))## found there are many missing values in Title column, so I will clean the dataset
##Inspect titles in the dataset
Title_2013_06 <- ML_2013[c("DiscoveryContactID","Title.x")]
Title_2013_12 <- ML_2013[c("DiscoveryContactID",  "Title.y")]
Title_2014_06 <- ML_2014[c("DiscoveryContactID","Title.x")]
Title_2014_12 <- ML_2014[c("DiscoveryContactID",  "Title.y")]
Title_2015_06 <- ML_2015[c("DiscoveryContactID","Title.x")]
Title_2015_12 <- ML_2015[c("DiscoveryContactID",  "Title.y")]
Title_2016_06 <- ML_2016[c("DiscoveryContactID","Title.x")]
Title_2016_12 <- ML_2016[c("DiscoveryContactID",  "Title.y")]
Title_2017_06 <- ML_2017[c("DiscoveryContactID","Title.x")]
Title_2017_12 <- ML_2017[c("DiscoveryContactID",  "Title.y")]
Title_2018_06 <- ML_2018[c("DiscoveryContactID","Title.x")]
Title_2018_12 <- ML_2018[c("DiscoveryContactID",  "Title.y")]
Title_2019_06 <- ML_2019[c("DiscoveryContactID","Title.x")]
Title_2019_12 <- ML_2019[c("DiscoveryContactID", "Title.y")]
Title_2020_06 <- ML_2020[c("DiscoveryContactID", "Title.x")]
Title_2020_12 <- ML_2020[c("DiscoveryContactID", "Title.y")]
Title_2021_06 <- ML_2021_06[c("DiscoveryContactID", "Title")]
library(ggplot2)
Title_all <- merge.data.frame(Title_2013_06, Title_2013_12, by = "DiscoveryContactID")
Title_all <- merge.data.frame(Title_all, Title_2014_06, by = "DiscoveryContactID")
Title_all <- merge.data.frame(Title_all, Title_2014_12, by = "DiscoveryContactID")
Title_all <- merge.data.frame(Title_all, Title_2015_06, by = "DiscoveryContactID")
Title_all <- merge.data.frame(Title_all, Title_2015_12, by = "DiscoveryContactID")
Title_all <- merge.data.frame(Title_all, Title_2016_06, by = "DiscoveryContactID")
Title_all <- merge.data.frame(Title_all, Title_2016_12, by = "DiscoveryContactID")
Title_all <- merge.data.frame(Title_all, Title_2017_06, by = "DiscoveryContactID")
Title_all <- merge.data.frame(Title_all, Title_2017_12, by = "DiscoveryContactID")
Title_all <- merge.data.frame(Title_all, Title_2018_06, by = "DiscoveryContactID")
Title_all <- merge.data.frame(Title_all, Title_2018_12, by = "DiscoveryContactID")
Title_all <- merge.data.frame(Title_all, Title_2019_06, by = "DiscoveryContactID")
Title_all <- merge.data.frame(Title_all, Title_2019_12, by = "DiscoveryContactID")
Title_all <- merge.data.frame(Title_all, Title_2020_06, by = "DiscoveryContactID")
Title_all <- merge.data.frame(Title_all, Title_2020_12, by = "DiscoveryContactID")
Title_all <- merge.data.frame(Title_all, Title_2021_06, by = "DiscoveryContactID")
colnames(Title_all) <- c("id", "2013_06", "2013_12", "2014_06", "2014_12", "2015_06","2015_12","2016_06","2016_12","2017_06","2017_12",
                         "2018_06","2018_12","2019_06","2019_12","2020_06","2020_12","2021_06")
dim(Title_all) #13143,18 
View(Title_all) ## we cannot simply drop missing values because some of them may have values in other columns
title_c <- c("2013_06", "2013_12", "2014_06", "2014_12", "2015_06","2015_12","2016_06","2016_12","2017_06","2017_12",
             "2018_06","2018_12","2019_06","2019_12","2020_06","2020_12","2021_06")
head_5 <- list()
for(t in 1:length(title_c)){
     element <- as.data.frame(sort(table(Title_all[title_c[t]]), decreasing = TRUE)[1:6])[,1]
     head_5[[t]] <- element
} ## we get the  highest 5 occurance of Titles each column
## now I want to find the number of occurance changing titles for each id
chan_list <- list()
nchan_list <- list()
for (i in 1:13143){
      row <- Title_all[i, 2:18]
      u_title <- unique(as.list(row[row != ""]))
      chan_list[[i]] <- u_title
      nchan_list[i] <- length(u_title)
}
nchan_c <- unlist(nchan_list)
length(nchan_c)
id_change <- data.frame(Title_all$id,  nchan_c)
length(chan_list[id_change$nchan_c >= 7]) ## check the career progression for for first id who has 7 chnages
# 7
length(chan_list[id_change$nchan_c == 6]) ## 88
length(chan_list[id_change$nchan_c == 5]) ## 535
length(chan_list[id_change$nchan_c == 4]) ## 2030
length(chan_list[id_change$nchan_c==3]) ## 3813
length(chan_list[id_change$nchan_c == 2]) ##3937
length(chan_list[id_change$nchan_c == 1]) ## 2027
length(chan_list[id_change$nchan_c == 0]) ## 706
id_7 <- (id_change %>% filter(nchan_c == 7))
id_7_n <- id_7$Title_all.id
id_1 <- (id_change %>% filter(nchan_c == 1))
id_1_n  <- id_1$Title_all.id
View(ML_2013%>% filter(DiscoveryContactID %in% id_7_n))
View(ML_2013%>% filter(DiscoveryContactID %in% id_1_n))
succe_2013_7 <- ML_2013 %>% filter(DiscoveryContactID %in% id_7_n)
succe_2013_1 <- ML_2013 %>% filter(DiscoveryContactID %in% id_1_n)
length(succe_2013_7$SuccessLikelihood.x)
length(succe_2013_7$SuccessLikelihood.x[succe_2013_7$SuccessLikelihood.x == "Higher"])
##0.7142857
length(succe_2013_1$SuccessLikelihood.x)
length(succe_2013_1$SuccessLikelihood.x[succe_2013_1$SuccessLikelihood.x == "Higher"])
##0.44
succe_2020_7 <- ML_2020 %>% filter(DiscoveryContactID %in% id_7_n)
length(succe_2020_7$SuccessLikelihood.x)
length(succe_2020_7$SuccessLikelihood.x[succe_2020_7$SuccessLikelihood.x == "Higher"])
id_5 <- (id_change %>% filter(nchan_c == 5))
id_5_n  <- id_5$Title_all.id
succe_2013_5 <- ML_2013 %>% filter(DiscoveryContactID %in% id_5_n)
length(succe_2013_5$SuccessLikelihood.x)
length(succe_2013_5$SuccessLikelihood.x[succe_2013_5$SuccessLikelihood.x == "Higher"])
##297 535 ##0.55
id_6 <- (id_change %>% filter(nchan_c == 6))
id_6_n  <- id_6$Title_all.id
succe_2013_6 <- ML_2013 %>% filter(DiscoveryContactID %in% id_6_n)
length(succe_2013_6$SuccessLikelihood.x)
length(succe_2013_6$SuccessLikelihood.x[succe_2013_6$SuccessLikelihood.x == "Higher"])
##40,88
##0.45
id_4 <- (id_change %>% filter(nchan_c == 4))
id_4_n  <- id_4$Title_all.id
succe_2013_4 <- ML_2013 %>% filter(DiscoveryContactID %in% id_4_n)
length(succe_2013_4$SuccessLikelihood.x)
length(succe_2013_4$SuccessLikelihood.x[succe_2013_4$SuccessLikelihood.x == "Higher"])
##1168; 2030
##0.575
id_3 <- (id_change %>% filter(nchan_c == 3))
id_3_n  <- id_3$Title_all.id
succe_2013_3 <- ML_2013 %>% filter(DiscoveryContactID %in% id_3_n)
length(succe_2013_3$SuccessLikelihood.x)
length(succe_2013_3$SuccessLikelihood.x[succe_2013_3$SuccessLikelihood.x == "Higher"])
##2328;3813
#0.61
id_2 <- (id_change %>% filter(nchan_c == 2))
id_2_n  <- id_2$Title_all.id
succe_2014_2 <- ML_2014 %>% filter(DiscoveryContactID %in% id_2_n)
length(succe_2014_2$SuccessLikelihood.x)
length(succe_2014_2$SuccessLikelihood.x[succe_2014_2$SuccessLikelihood.x == "Higher"])
##2026;3937
## 7-1: 0.7 0.45 0.55 0.575 0.61 0.50
## now I want to draw the histogram to show the distribution of data
tadata_col <- names(Title_all)[2:18]
for (i in 1:17){
  hist(log(as.data.frame(table(Title_all[tadata_col[i]]))$Freq))
}
## investigate the career progression for head5 title in 2013_06
head5_2013_06 <- head_5[1]
View(table(Title_all$`2014_06`[Title_all$`2013_06` == "Financial Advisor"]))
View(log(as.data.frame(table(Title_all$`2015_06`[Title_all$`2013_06` == "Financial Advisor"]))$Freq))
## we can see the distrbution of frequency of the titles
summary(as.data.frame(table(Title_all$`2015_06`[Title_all$`2013_06` == "Financial Advisor"]))$Freq)
###########################################################################
#first attempt to write transition matrix
#I will assume the five states of the matrix are head_5$'2013_06
five_states <- unlist(head_5[1])[2:6]
fif_tran_mat <- matrix(0,5,5)
library(data.table)
Title_all_1 <- copy(Title_all)
levels <- five_states
labels <- c(1,2,3,4,5)
for(i in 2:18){
  Title_all_1[,i] <- factor(Title_all_1[,i], levels = levels, labels = labels, ordered = TRUE)
  
}
for (i in 1:13143){
      null_index <- c()
  for (j in 2:17){
      cur_ <- Title_all_1[i,j]
      next_ <- Title_all_1[i,j+1]
    if ((!is.nan(cur_)) & (!is.nan(next_)) ){
      cur_v <- as.numeric(cur_)
      next_ <- as.numeric(next_)
      fif_tran_mat[cur_v,next_] <- fif_tran_mat[cur_v,next_] + 1
    }
    if ((!is.nan(cur_v)) & (is.nan(next_))){
      Title_all_1[i,j+1] <- Title_all_1[i,j]
      append(null_index, j+1)
    }
  }
      for(val in null_index){
        Title_all_1[i, val] <- null
      }
}
fif_tran_mat/rowSums(fif_tran_mat)

###
head_10 <- list()
for(t in 1:length(title_c)){
  element <- as.data.frame(sort(table(Title_all[title_c[t]]), decreasing = TRUE)[1:11])[,1]
  head_10[[t]] <- element
}
ten_states <- unlist(head_10[1])[2:11]
ten_tran_mat <- matrix(0,10,10)
Title_all_2 <- copy(Title_all)
levels <- ten_states
labels <- c(1,2,3,4,5,6,7,8,9,10)
for(i in 2:18){
  Title_all_2[,i] <- factor(Title_all_2[,i], levels = levels, labels = labels, ordered = TRUE)
  
}
for (i in 1:13143){
  null_index <- c()
  for (j in 2:17){
    cur_ <- Title_all_2[i,j]
    next_ <- Title_all_2[i,j+1]
    if ((!is.nan(cur_)) & (!is.nan(next_)) ){
      cur_v <- as.numeric(cur_)
      next_ <- as.numeric(next_)
      ten_tran_mat[cur_v,next_] <- ten_tran_mat[cur_v,next_] + 1
    }
    if ((!is.nan(cur_v)) & (is.nan(next_))){
      Title_all_2[i,j+1] <- Title_all_2[i,j]
      append(null_index, j+1)
    }
  }
  for(val in null_index){
    Title_all_2[i, val] <- null
  }
}
ten_tran_mat/rowSums(ten_tran_mat)
## next, I will find the 10 x 10 transition matrix for all data
## I will find top ten titles for all data
concat_all <-Title_all$`2013_06`
for (t in 2:(length(title_c)-1)){
  concat_all <- append(concat_all,Title_all[title_c[t]])
}

ten_state1 <- as.data.frame(sort(table(unlist(concat_all)),decreasing = TRUE))[2:11,1]
ten_tran_mat_1 <- matrix(0,10,10)
Title_all_3 <- copy(Title_all)[,1:17]
levels <- ten_state1
labels <- c(1,2,3,4,5,6,7,8,9,10)
for(i in 2:17){
  Title_all_3[,i] <- factor(Title_all_3[,i], levels = levels, labels = labels, ordered = TRUE)
  
}
for (i in 1:13143){
  null_index <- c()
  for (j in 2:16){
    cur_ <- Title_all_3[i,j]
    next_ <- Title_all_3[i,j+1]
    if ((!is.nan(cur_)) & (!is.nan(next_)) ){
      cur_v <- as.numeric(cur_)
      next_ <- as.numeric(next_)
      ten_tran_mat_1[cur_v,next_] <- ten_tran_mat_1[cur_v,next_] + 1
    }
    if ((!is.nan(cur_v)) & (is.nan(next_))){
      Title_all_3[i,j+1] <- Title_all_3[i,j]
      append(null_index, j+1)
    }
  }
  for(val in null_index){
    Title_all_3[i, val] <- null
  }
}
ten_tran_mat_1
tm_10  <- ten_tran_mat_1/rowSums(ten_tran_mat_1)
View(table(ten_state1))
## now we want to run Simple MC model
install.packages("markovchain")
library(markovchain)

MC_title <- new("markovchain", states = as.vector(ten_state1), byrow = T, transitionMatrix = tm_10, name = "TitleProgression")
MC_st <- steadyStates(MC_title) ##cannot use steady state
#now I want to estimate the value of ten_state1 in 2021_06
data_2021_06 <-  Title_all$`2021_06`[Title_all$`2021_06` %in% ten_state1]##actual value 
View(table(data_2021_06))
MCdata_2021_06 <- sum(as.data.frame(table(data_2021_06))$Freq)*MC_st
data_2020_12 <-  Title_all$`2020_12`[Title_all$`2020_12` %in% ten_state1]
MCMCdata_2021_06<- t(as.matrix(as.data.frame(table(data_2020_12))$Freq)) %*% tm_10
MCMCdata_2021_06 
##now we want to find the sqrt mean square error
MCMC_2021_06_diff <- (MCMCdata_2021_06 - as.data.frame(table(data_2021_06))$Freq)
for (i in 1: 10){
  MCMC_2021_06_diff[i] <- (MCMC_2021_06_diff[i])^2 ##squared sum
}
MCMC_2021_06_error <- sqrt(sum(MCMC_2021_06_diff)/10)
MCMC_2021_06_error #30.5876
##the prediction is quite good, I want to form a table to show the model prediction
titlenames <- as.vector(as.data.frame(table(data_2021_06))[,1])
data2021_06 <- as.vector(as.data.frame(table(data_2021_06))[,2])
MCMC2021_06 <- as.vector(MCMCdata_2021_06)
MCMC_diff <- data2021_06 - MCMC2021_06
MCMC2021_06_report <- data.frame(titlenames, data2021_06, MCMC2021_06, MCMC_diff)
View(MCMC2021_06_report)
##I will make the time period as three years, 2013-06, 2016-06, 2019-06
Title_2yearsp <- Title_all[c("2013_06", "2016_06", "2019_06")]
View(Title_2yearsp)
concat_all2 <-Title_2yearsp$`2013_06`
title_c2 <-c("2013_06", "2016_06","2019_06")
for (t in 2:length(title_c2)){
  concat_all2 <- append(concat_all2,Title_2yearsp[title_c2[t]])
}

ten_state2 <- as.data.frame(sort(table(unlist(concat_all2)),decreasing = TRUE))[2:11,1]
ten_state2 <-ten_state2[c(1,2,3,4,5,6,7,8,9,11)]
View(table(ten_state2))
ten_tran_mat_2 <- matrix(0,10,10)
Title_all_4 <- copy(Title_2yearsp)
levels <- ten_state2
labels <- c(1,2,3,4,5,6,7,8,9,10)
for(i in 1:3){
  Title_all_4[,i] <- factor(Title_all_4[,i], levels = levels, labels = labels, ordered = TRUE)
  
}
for (i in 1:13143){
  null_index <- c()
  for (j in 1:2){
    cur_ <- Title_all_4[i,j]
    next_ <- Title_all_4[i,j+1]
    if ((!is.nan(cur_)) & (!is.nan(next_)) ){
      cur_v <- as.numeric(cur_)
      next_ <- as.numeric(next_)
      ten_tran_mat_2[cur_v,next_] <- ten_tran_mat_2[cur_v,next_] + 1
    }
    if ((!is.nan(cur_v)) & (is.nan(next_))){
      Title_all_4[i,j+1] <- Title_all_4[i,j]
      append(null_index, j+1)
    }
  }
  for(val in null_index){
    Title_all_4[i, val] <- null
  }
}
ten_tran_mat_2
tm_10_2  <- ten_tran_mat_2/rowSums(ten_tran_mat_2)
View(table(ten_state2))
## now we want to run Simple MC model
install.packages("markovchain")
library(markovchain)

MC_title2 <- new("markovchain", states = as.vector(ten_state2), byrow = T, transitionMatrix = tm_10_2, name = "TitleProgression2")
MC_st2 <- steadyStates(MC_title2) ##cannot use steady state
#now I want to estimate the value of ten_state1 in 2021_06
data_2021_06 <-  Title_all$`2021_06`[Title_all$`2021_06` %in% ten_state2]##actual value 
View(table(data_2021_06))
MC2data_2021_06 <- sum(as.data.frame(table(data_2021_06))$Freq)*MC_st2
data_2019_06 <-  Title_2yearsp$`2019_06`[Title_2yearsp$`2019_06` %in% ten_state2]
MCMCdata2_2021_06<- t(as.matrix(as.data.frame(table(data_2019_06))$Freq)) %*% tm_10_2
MCMCdata2_2021_06
titlenames <- as.data.frame(table(data_2021_06))[,1]
data2021_06 <- as.vector(as.data.frame(table(data_2021_06))[,2])
MCMC2021_06 <- as.vector(MCMCdata2_2021_06)
#MCMC_diff <- data2021_06 - MCMC2021_06
data_2021_06df <- data.frame(title = titlenames, data2021_06)
MCMC2021_06df<- data.frame(title = ten_state2, MCMC2021_06)
MCMC2021_06_report2 <- merge(data_2021_06df, MCMC2021_06df, by = "title")
View(MCMC2021_06_report2)
MCMC_2021_06_diff2 <- MCMC2021_06_report2[,2] - MCMC2021_06_report2[,3]
MCMC2021_06_report2 <- cbind(MCMC2021_06_report2, MCMC_2021_06_diff2)
#  MCMC_2021_06_diff2[i] <- (MCMC_2021_06_diff2[i])^2 ##squared sum
#}
#MCMC_2021_06_error2 <- sqrt(sum(MCMC_2021_06_diff2)/10)
#MCMC_2021_06_error2#233.92
selected_06M_data <- Title_all[c("id", "2013_06","2014_06","2015_06", "2016_06","2017_06","2018_06","2019_06","2020_06","2021_06")]
concat_06M <-selected_06M_data$`2013_06`
for (t in 3:(length(selected_06M_data)-1)){
  concat_06M <- append(concat_06M,selected_06M_data[,t])
}
ten_state06M <- as.vector(as.data.frame(sort(table(unlist(concat_06M)),decreasing = TRUE))[2:11,1])
top10_selected_06 <- selected_06M_data %>% filter(`2013_06` %in% ten_state06M,
                                                  `2014_06`%in% ten_state06M,
                                                  `2015_06` %in% ten_state06M,
                                                  `2016_06` %in% ten_state06M,
                                                  `2017_06` %in% ten_state06M,
                                                  `2018_06` %in% ten_state06M,
                                                  `2019_06` %in% ten_state06M,
                                                  `2020_06` %in% ten_state06M)
                                                  #`2021_06`%in% ten_state06M)
View(top10_selected_06)
dim(top10_selected_06) #1307,10
ten_tran_mat_5 <- matrix(0,10,10)
Title_all_5 <- copy(top10_selected_06)
levels <- ten_state06M
labels <- c(1,2,3,4,5,6,7,8,9,10)
for(i in 2:9){
  Title_all_5[,i] <- factor(Title_all_5[,i], levels = levels, labels = labels, ordered = TRUE)
  
}
for (i in 1:1307){
  for (j in 2:8){
    cur_ <- Title_all_5[i,j]
    next_ <- Title_all_5[i,j+1]
    cur_v <- as.numeric(cur_)
    next_ <- as.numeric(next_)
    ten_tran_mat_5[cur_v,next_] <- ten_tran_mat_5[cur_v,next_] + 1
  }
}
ten_tran_mat_5
tm_10_5  <- ten_tran_mat_5/rowSums(ten_tran_mat_5)
tm_10_5
View(tm_10_5)
##we want to find the communication between state first
c_matrix <- tm_10_5 > 0
count <- 0
for(i in 1:9){
  for(j in (i+1):10){
    if(c_matrix[i,j]+c_matrix[j,i] == 2){ ##have communication
      print(c(i,j))
      count <- count +1
    }
  }
}
count
MC_title5 <- new("markovchain", states = ten_state06M, byrow = T, transitionMatrix = tm_10_5, name = "TitleProgression5")
MC_st5 <- steadyStates(MC_title5)
MC_st5
##title 12-21 1 class, and 2 recurrent states, 8 transient states
dim(MC_st5)

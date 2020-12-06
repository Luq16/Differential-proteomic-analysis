# read proteinGroups.txt file
raw_df = read.delim("proteinGroups.txt", stringsAsFactors = FALSE,
                    colClasses = "character")
# for data manipulation
library(dplyr)

# Filter out contaminants hits
raw_df = raw_df %>%
  filter(Potential.contaminant != "+") %>%
  filter(Reverse != "+") %>%
  filter(Only.identified.by.site != "+")


#select interested columns
df<-dplyr::select(raw_df,Majority.protein.IDs,Gene.names, Protein.names,Unique.peptides,
                  Unique.peptides.rep1_Bkd_ctrl,Unique.peptides.rep2_Bkd_ctrl,Unique.peptides.rep3_Bkd_ctrl,
                  Unique.peptides.rep1_ctr_5min,Unique.peptides.rep2_ctr_5min,Unique.peptides.rep3_ctr_5min,
                  Unique.peptides.rep1_ctr_10min,Unique.peptides.rep2_ctr_10min,Unique.peptides.rep3_ctr_10min,
                  Unique.peptides.rep1_ctr_15min,Unique.peptides.rep2_ctr_15min,Unique.peptides.rep3_ctr_15min,
                  Unique.peptides.rep1_NA_5min,Unique.peptides.rep2_NA_5min,Unique.peptides.rep3_NA_5min,
                  Unique.peptides.rep1_NA_10min,Unique.peptides.rep2_NA_10min,Unique.peptides.rep3_NA_10min,
                  Unique.peptides.rep1_NA_15min,Unique.peptides.rep2_NA_15min,Unique.peptides.rep3_NA_15min,
                  Unique.peptides.rep1_fab_5min,Unique.peptides.rep2_fab_5min,Unique.peptides.rep3_fab_5min,
                  Unique.peptides.rep1_fab_10min,Unique.peptides.rep2_fab_10min,Unique.peptides.rep3_fab_10min,
                  Unique.peptides.rep1_fab_15min,Unique.peptides.rep2_fab_15min,Unique.peptides.rep3_fab_15min,
                  Intensity.rep1_Bkd_ctrl,Intensity.rep2_Bkd_ctrl, Intensity.rep3_Bkd_ctrl,
                  Intensity.rep1_ctr_5min, Intensity.rep2_ctr_5min, Intensity.rep3_ctr_5min,
                  Intensity.rep1_ctr_10min, Intensity.rep2_ctr_10min, Intensity.rep3_ctr_10min,
                  Intensity.rep1_ctr_15min, Intensity.rep2_ctr_15min, Intensity.rep3_ctr_15min,
                  Intensity.rep1_NA_5min,Intensity.rep2_NA_5min,Intensity.rep3_NA_5min,
                  Intensity.rep1_NA_10min, Intensity.rep2_NA_10min, Intensity.rep3_NA_10min,
                  Intensity.rep1_NA_15min, Intensity.rep2_NA_15min, Intensity.rep3_NA_15min,
                  Intensity.rep1_fab_5min, Intensity.rep2_fab_5min, Intensity.rep3_fab_5min,
                  Intensity.rep1_fab_10min, Intensity.rep2_fab_10min, Intensity.rep3_fab_10min,
                  Intensity.rep1_fab_15min, Intensity.rep2_fab_15min, Intensity.rep3_fab_15min,
                  Unique.peptides, Reverse, Potential.contaminant, Mol..weight..kDa.)

#filter out proteins not identified in background control
dt_not_bckgrnd=df%>%filter(Intensity.rep1_Bkd_ctrl==0&Intensity.rep2_Bkd_ctrl==0&Intensity.rep3_Bkd_ctrl==0)
#filter out proteins identified in background control
dat_bckgrnd=df%>%filter(Intensity.rep1_Bkd_ctrl>0|Intensity.rep2_Bkd_ctrl>0|Intensity.rep3_Bkd_ctrl>0)

#convert columns to numeric
dat_bckgrnd[,4:64]<-sapply(dat_bckgrnd[,4:64], as.numeric)
dt_not_bckgrnd[,4:64]<-sapply(dt_not_bckgrnd[,4:64], as.numeric)

#make copies of variables
dat_4=dat_bckgrnd
dat_5=dt_not_bckgrnd

#log transform intensity columns
dat_4[,35:64]<-log2(dat_4[,35:64])
dat_5[,35:64]<-log2(dat_5[,35:64])

#convert infinite to NAs
dat_4[mapply(is.infinite, dat_4)] <- NA
dat_5[mapply(is.infinite, dat_5)] <- NA

#calculate mean of log intensity  identified in background control
dat_4$mean_bckd_int<-rowMeans(subset(dat_4, select = c(Intensity.rep1_Bkd_ctrl,Intensity.rep2_Bkd_ctrl,
                                                       Intensity.rep3_Bkd_ctrl)), na.rm = TRUE)

#calculate mean of log intensity  identified in NonAct
dat_4$mean_NA_int<-rowMeans(subset(dat_4, select = c(Intensity.rep1_NA_5min,Intensity.rep2_NA_5min,Intensity.rep3_NA_5min,
                                                     Intensity.rep1_NA_10min, Intensity.rep2_NA_10min, Intensity.rep3_NA_10min,
                                                     Intensity.rep1_NA_15min, Intensity.rep2_NA_15min, Intensity.rep3_NA_15min)), na.rm = TRUE)


#subtract mean intensity of background control from control (no H2O2)
dat_6=dat_4%>%mutate(fold=dat_4$mean_NA_int-dat_4$mean_bckd_int)

#filter out fold < 3
dat_7=filter(dat_6, dat_6$fold>=3)

#add extra column to be able to use rbind
dat_5$mean_bckd_int<-NA
dat_5$mean_NA_int<-NA
dat_5$fold<-NA

#bind dat_5 and Data_7:
data<-rbind(dat_5,dat_7)

data_3.1=dplyr::select(data,Majority.protein.IDs,Gene.names, Protein.names,Unique.peptides,fold,
                       Intensity.rep1_NA_5min,Intensity.rep2_NA_5min,Intensity.rep3_NA_5min,
                       Intensity.rep1_NA_10min, Intensity.rep2_NA_10min, Intensity.rep3_NA_10min,
                       Intensity.rep1_NA_15min, Intensity.rep2_NA_15min, Intensity.rep3_NA_15min,
                       Intensity.rep1_fab_5min, Intensity.rep2_fab_5min, Intensity.rep3_fab_5min,
                       Intensity.rep1_fab_10min, Intensity.rep2_fab_10min, Intensity.rep3_fab_10min,
                       Intensity.rep1_fab_15min, Intensity.rep2_fab_15min, Intensity.rep3_fab_15min,
                       Unique.peptides, Reverse, Potential.contaminant, Mol..weight..kDa.,mean_bckd_int,mean_NA_int)
data_3.3=subset(data_3.1,rowSums(is.na(data_3.1[6:14]))<8|rowSums(is.na(data_3.1[15:17]))<2|rowSums(is.na(data_3.1[18:20]))<2|rowSums(is.na(data_3.1[21:23]))<2)#at leat 2 in one condition

#filter >2 unique peptide
data_3.4<-filter(data_3.3, Unique.peptides>1)

data_3.4$Majority.protein.IDs<- as.character(data_3.4$Majority.protein.IDs)
data_3.4$Majority.protein.IDs<- sub(';.*',"",data_3.4$Majority.protein.IDs)

data_3.5=subset(data_3.4,rowSums(is.na(data_3.4[6:23]))<7)#2 out of 3

#####Imputation kNN
library(MSnbase)
library(data.table)
dat_impt=readMSnSet2(data_3.5, ecol= 6:23, fnames = "Majority.protein.IDs")#msnset mode

knn_imputation<- impute(dat_impt, "knn")#imputation
knn_imputation=data.frame(knn_imputation)%>%t()%>%data.frame()#transpose
knn_imputation$Majority.protein.IDs=rownames(knn_imputation)#make rowname columns
data_3.6=select(data_3.5, Majority.protein.IDs, Protein.names, Gene.names)#select from data_3.5. to prpare for merge
knn_imputation=merge(knn_imputation,data_3.6, by="Majority.protein.IDs")
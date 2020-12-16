<<<<<<< HEAD
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


###########################Differentia analysis
library(NormalyzerDE)
dataFp <- system.file(package="NormalyzerDE", "extdata", "knn_imputation2of3Impt.txt")
designFp <- system.file(package="NormalyzerDE", "extdata", "design2.txt")
outDir <- "/Users/luqmanawoniyi/Desktop/proteomicData/data_2019/normalizerDE/filterByR/filter2Of3_imput"
normalyzer(jobName="knn_imputation2of3Impt", noLogTransform= TRUE, designPath=designFp, dataPath=dataFp,
           outputDir=outDir)

#stat
normMatrixPath <- paste(outDir, "knn_imputation2of3Impt/Quantile-normalized.txt", sep="/")
normalyzerDE("knn_imputation2of3Impt_norm",designFp,normMatrixPath,outputDir=outDir,comparisons=c("4-1", "5-2", "6-3"),condCol="group")
###################################################################################/DE


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


###########################Differentia analysis
library(NormalyzerDE)
dataFp <- system.file(package="NormalyzerDE", "extdata", "knn_imputation2of3Impt.txt")
designFp <- system.file(package="NormalyzerDE", "extdata", "design2.txt")
outDir <- ""
normalyzer(jobName="knn_imputation2of3Impt", noLogTransform= TRUE, designPath=designFp, dataPath=dataFp,
           outputDir=outDir)

#stat
normMatrixPath <- paste(outDir, "knn_imputation2of3Impt/Quantile-normalized.txt", sep="/")
normalyzerDE("knn_imputation2of3Impt_norm",designFp,normMatrixPath,outputDir=outDir,comparisons=c("4-1", "5-2", "6-3"),condCol="group")
###################################################################################/DE

#####################################################significants and FC Filteration
setwd("/Users/luqmanawoniyi/Desktop/proteomicData/data_2019/normalizerDE/filterByR/filter2Of3_imput/knn_imputation2of3Impt_norm")
data_5=read.table(file = 'knn_imputation2of3Impt_norm_stats.csv', sep = ';', header = TRUE)
sig_all<-data_5%>%filter(X4.1_AdjPVal<=0.05|X5.2_AdjPVal<=0.05|X6.3_AdjPVal<=0.05)
#write.csv2(sig_all, "sigall.csv")
#sig in 5min FC not considered
sig_5min_noFC<-data_5%>%filter(X4.1_AdjPVal<=0.05)
#sig in 10min FC not considered
sig_10min_noFC<-data_5%>%filter(X5.2_AdjPVal<=0.05)
#sig in 15min FC not considered
sig_15min_noFC<-data_5%>%filter(X6.3_AdjPVal<=0.05)
##cal. squre root
sig_all$FC5min<-'^'(sig_all$X4.1_log2FoldChange,2)
sig_all$FC10min<-'^'(sig_all$X5.2_log2FoldChange,2)
sig_all$FC15min<-'^'(sig_all$X6.3_log2FoldChange,2)

#filter significant in all
sigFC_all_1_5<-sig_all%>%filter(FC5min>2.25|FC10min>2.25|FC15min>2.25)
sigFC_all_2<-sig_all%>%filter(FC5min>=4|FC10min>=4|FC15min>=4)
sig_5min<-sig_all%>%filter(FC5min>=2.25&X4.1_AdjPVal<=0.05)
sig_10min<-sig_all%>%filter(FC10min>=2.25&X5.2_AdjPVal<=0.05)
sig_15min<-sig_all%>%filter(FC15min>=2.25&X6.3_AdjPVal<=0.05)
write.csv(sig_10min, "sig_10min.csv")
write.csv(sig_5min, "sig_5min.csv")
write.csv(sig_15min, "sig_15min.csv")
write.csv(sigFC_all_1_5, "sigall_1.5FC.csv")
write.csv(sigFC_all_2, "sigall_2FC.csv")
write.csv(sig_all, "sig_all.csv")
write.csv(sig_10min_noFC, "sig_10min_noFc.csv")
write.csv(sig_5min_noFC, "sig_5min_noFc.csv")
write.csv(sig_15min_noFC, "sig_15min_noFc.csv")

######Up and down reg
dat<-read.csv("sig_5min_noFc.csv")
dat_2<-read.csv("sig_10min_noFc.csv")
dat_3<-read.csv("sig_15min_noFc.csv")
up_5min<-dat%>%filter(X4.1_log2FoldChange>0)
up_10min<-dat_2%>%filter(X5.2_log2FoldChange>0)
up_15min<-dat_3%>%filter(X6.3_log2FoldChange>0)

down_5min<-dat%>%filter(X4.1_log2FoldChange<0)
down_10min<-dat_2%>%filter(X5.2_log2FoldChange<0)
down_15min<-dat_3%>%filter(X6.3_log2FoldChange<0)
write.csv(up_5min, "up_5min.csv")
write.csv(up_10min, "up_10min.csv")
write.csv(up_15min, "up_15min.csv")

write.csv(down_5min, "down_5min.csv")
write.csv(down_10min, "down_10min.csv")
write.csv(down_15min, "down_15min.csv")

#####################################################/end of Filtering and significant filter

#######################volcano plot using enhanceed volcano
#make gene names rowname
names<- make.unique(as.character(data_5$Gene.names))
rownames(data_5) = make.names(names, unique=TRUE)
setwd("/Users/luqmanawoniyi/Desktop/proteomicData/data_2019/normalizerDE/filterByR/filter2Of3_imput/knn_imputation2of3Impt_norm")

library(EnhancedVolcano)
#5min
library(svglite)
svglite("volc5.svg")
EnhancedVolcano(data_5,
                lab = rownames(data_5),
                x = 'X4.1_log2FoldChange',
                y = 'X4.1_AdjPVal',
                xlim=c(-6,8.5),
                ylim = c(-0.2,7.5),
                title = 'non-activated versus 5min Fab2 activated',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                titleLabSize = 20,
                shape = 20,
                pointSize=3,
                colAlpha=1,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                #selectLab=c("IgK.V", "Ighm", "Myo1g","Irf5", "Tacc1" ,"Usp7", "Cybb", "Hspa14", "Tor1aip1",
                #         "Shq1", "Sssca1", "Polr2c", "Stk10", "Plcb3", "Mapk1", "Nfkbia", "Lyn","Lmtk2",
                #        "Vti1b","Smad2","Cd79a", "Cd79b"),
                selectLab=c("IgK.V","Ighm","Shq1", "Polr2c","Lyn","Pik3r1","Mark2"),
                colConnectors = 'grey10',
                legend=c('NS','Log2(fold-change)','P value',
                         'P value & Log2(fold-change'),
                subtitle="Differential expression",
                axisLabSize= 22,
                labSize=5,
                boxedLabels=TRUE
)
dev.off()

svglite("volc10.svg")
EnhancedVolcano(data_5,
                lab = rownames(data_5),
                x = 'X5.2_log2FoldChange',
                y = 'X5.2_AdjPVal',
                xlim=c(-6,8.5),
                ylim = c(-0.2,10),
                title = 'non-activated versus 10min Fab2 activated',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                titleLabSize = 20,
                shape = 20,
                pointSize=3,
                colAlpha=1,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                selectLab=c("IgK.V","Ighm","Shq1", "Polr2c","Lyn","Pik3r1","Mark2"),
                colConnectors = 'grey10',
                legend=c('NS','Log2(fold-change)','P value',
                         'P value & Log2(fold-change'),
                subtitle="Differential expression",
                axisLabSize= 22,
                labSize=5,
                boxedLabels=TRUE
)
dev.off()

svglite("volc15.svg")
EnhancedVolcano(data_5,
                lab = rownames(data_5),
                x = 'X6.3_log2FoldChange',
                y = 'X6.3_AdjPVal',
                xlim=c(-6,8.5),
                ylim = c(-0.2,7.5),
                title = 'non-activated versus 15min Fab2 activated',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                titleLabSize = 20,
                shape = 20,
                pointSize=3,
                colAlpha=1,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                selectLab=c("IgK.V","Ighm","Shq1", "Polr2c","Lyn","Pik3r1","Mark2"),
                colConnectors = 'grey10',
                legend=c('NS','Log2(fold-change)','P value',
                         'P value & Log2(fold-change'),
                subtitle="Differential expression",
                axisLabSize= 22,
                labSize=5,
                boxedLabels=TRUE
)
dev.off()
#############################################/end of volcano

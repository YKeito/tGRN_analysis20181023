#"~/Nakano_RNAseq/network_analysis/script/tGRN_analysis20181023/tGRN_masterTable.R"
#################################################################################################set##########################################################################################
allgenes <- read.table("~/Nakano_RNAseq/network_analysis/base/allgenes_attribute.txt", sep = "\t", header = T, stringsAsFactors = F)
allgenes <- allgenes$AGI
TAIR <- read.table("~/Nakano_RNAseq/network_analysis/base/Tair_annotation.txt", sep = "\t", fill = T, quote = "", header = T, stringsAsFactors = F)
TaleMine <- read.table("~/Nakano_RNAseq/network_analysis/base/TaleMine.csv", sep = ",", fill = T, quote = "", header = T, stringsAsFactors = F)

#AGI
Botrytis_cinerea <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/defense/pathogen infection/B.cinerea_newlist.txt", sep = "\t", header = T, stringsAsFactors = F)
Botrytis_cinerea$Gene.Locus <- toupper(Botrytis_cinerea$Gene.Locus)
colnames(Botrytis_cinerea) <- "AGI"
PstDC3000_hrp <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/defense/pathogen infection/PstDC3000_hrp-_DEGs.txt", sep = "\t", header = T, stringsAsFactors = F)
MeJA_DEGs <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/MeJA_DEGs.txt", sep = "\t", header = T, stringsAsFactors = F)
BTH <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/defense/BTH_affected_AGI.txt", sep = "\t", header = T, stringsAsFactors = F)
wrky18 <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/ChIPData/wrky18_2h_flg22_TargetGenes.txt", sep = "\t", header = T, stringsAsFactors = F)
wrky33 <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/ChIPData/wrky33_2h_flg22_TargetGenes.txt", sep = "\t", header = T, stringsAsFactors = F)
wrky40 <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/ChIPData/wrky40_2h_flg22_TargetGenes.txt", sep = "\t", header = T, stringsAsFactors = F)
colnames(wrky18) <- "AGI"
colnames(wrky33) <- "AGI"
colnames(wrky40) <- "AGI"


CY <- c("CY15", "CY16", "CY20")
time <- c("1h", "3h", "12h", "24h")
allCY <- c()
obnames <- c()
i <- 1
for(i in i:length(CY)){
  n <- 1
  for(n in n:length(time)){
    temp <- data.frame(AGI = rownames(allRNASeq[allRNASeq[, paste0(CY[i], "_", time[n], "_q_value")] < 0.05, ]), stringsAsFactors = F)
    T_data <- rep("No", times = length(allgenes))
    names(T_data) <- allgenes
    T_data[match(temp[, "AGI"], names(T_data))] <- "Yes"
    allCY <- cbind(allCY, T_data)
    obnames <- c(obnames, paste0(CY[i], "_", time[n]))
    n <- n+1
  }
  i <- i+1
}
colnames(allCY) <- paste0(obnames, "_FDR0.05")

temp <- list(Botrytis_cinerea, PstDC3000_hrp, MeJA_DEGs, BTH, wrky18, wrky33, wrky40)
obnames <- c("Botrytis_cinerea", "PstDC3000_hrp", "MeJA_DEGs", "BTH", "wrky18", "wrky33", "wrky40")
allsample <- c()
i <- 1
for(i in i:length(obnames)){
  T_data <- rep("No", times = length(allgenes))
  names(T_data) <- allgenes
  T_data[match(temp[[i]][, "AGI"], names(T_data))] <- "Yes"
  allsample <- cbind(allsample, T_data)
  print(i)
  i <- i+1
}
colnames(allsample) <- obnames

#AGI, TF
TF_family <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/TF/arabidopsis_TF_family.txt", sep = "\t", header = T, stringsAsFactors = F)
T_data <- rep("No", times = length(allgenes))
names(T_data) <- allgenes
T_data[match(TF_family$AGI, names(T_data))] <- TF_family$TF


tGRN <- read.table("~/bigdata/yasue/tGRN_Groping/NodeTable/inflation4.csv", header=T, sep=",", stringsAsFactors = F)
#MCL
MCLNum <- rep(NA, times = length(allgenes))
names(MCLNum) <- allgenes
MCLNum[match(tGRN$name, names(MCLNum))] <- tGRN$X__mclCluster

#Degree
Degree <- rep(0, times = length(allgenes))
names(Degree) <- allgenes
Degree[match(tGRN$name, names(Degree))] <- tGRN$Degree

#BetweennessCentrality
BC <- rep(0, times = length(allgenes))
names(BC) <- allgenes
BC[match(tGRN$name, names(BC))] <- tGRN$BetweennessCentrality

MasterTable <- data.frame(AGI = allgenes,
                          symbol = TaleMine$symbol[match(allgenes, TaleMine$AGI)],
                          TF = T_data,
                          MCLNum = MCLNum,
                          Degree = Degree,
                          BetweennessCentrality = BC,
                          allsample,
                          allCY,
                          annotation = TAIR$annotation[match(allgenes, TAIR$AGI)],
                          stringsAsFactors = F
)
rownames(MasterTable) <- c()
#save(MasterTable, file = "~/bigdata/yasue/tGRN_Groping/inflation4/.RData/New_MasterTable_inflation4.RData")

#write.table(MasterTable, file = "bigdata/yasue/tGRN_Groping/inflation4/MasterTable/New_MasterTable_inflation4.txt", append=F, quote = F, sep = "\t", row.names = F)
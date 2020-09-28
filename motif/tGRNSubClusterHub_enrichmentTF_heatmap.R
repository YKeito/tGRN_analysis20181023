#"~/Nakano_RNAseq/network_analysis/script/tGRN_analysis20181023/motif/tGRNSubClusterHub_enrichmentTF_heatmap.R"
TFname <- unique(tGRN_MotifIDTF$TF)
TFname <- TFname[TFname != ""]
T.data <- MasterTable[!is.na(MasterTable$MCLNum), ]
T_TF <- list(T.data %>% filter(TF != TFname[n], MCLNum < 243) %>% select(MCLNum) %>% unlist(., use.names = F))
n <- 1
for(n in n:length(TFname)){
  T_TF <- c(T_TF, list(T.data %>% filter(TF == TFname[n], MCLNum < 243) %>% select(MCLNum) %>% unlist(., use.names = F)))
}
names(T_TF) <- c("others", TFname)

temp <- unique(unlist(T_TF, use.names = F))
count <- c()
n <- 1
for(n in n:length(T_TF)){
  m <- 1
  for(m in m:length(temp)){
    count <- c(count, sum(unlist(T_TF) == temp[m]))
    m <- m+1
  }
  n <- n+1
}

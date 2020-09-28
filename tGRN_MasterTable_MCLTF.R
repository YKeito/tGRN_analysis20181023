#"~/Nakano_RNAseq/network_analysis/script/tGRN_analysis20181023/tGRN_MasterTable_MCLTF.R"
T_data <- MasterTable[!is.na(MasterTable$MCLNum), ]
temp <- unique(T_data$MCLNum)
test <- c()
i <- 1
for(i in i:length(temp)){
  test <- c(test, sum(T_data$MCLNum == i))
}


T_data <- T_data[T_data$MCLNum <= length(test[test >= 3]), ]
T_MCLNum <- c()
i <- 1
for(i in i:length(unique(T_data$MCLNum))){
  if(sum(T_data$TF[T_data$MCLNum == i] != "No") != 0){
    T_MCLNum <- c(T_MCLNum, i)
  }
  i <- i+1
}

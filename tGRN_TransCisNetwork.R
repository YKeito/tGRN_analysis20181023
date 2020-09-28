"~/Nakano_RNAseq/network_analysis/script/tGRN_analysis20181023/tGRN_TransCisNetwork.R"
library(dplyr)
#マニュアルで転写因子コア結合配列を分類したデータを読み込ませた
tGRN_MotifID_consensus
##Cis側----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#シス側:転写因子ファミリーの情報を基にクラスターを分類
TFname <- unique(unique(names(MotifList)))
allCis <- list()
i <- 1
for(i in i:length(TFname)){
  allCis <- c(allCis, list(unique(tGRN_MotifID_consensus$MCLNum[grep(TFname[i], tGRN_MotifID_consensus$TF)])))
  i <- i+1
}
names(allCis) <- TFname

#trans側---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#ノード数3以上で転写因子があるクラスター番号の抽出とその転写因子ファミリーの取得---------------------------------------------------------------------------------------------------------
T_data <- MasterTable[!is.na(MasterTable$MCLNum), ]
temp <- unique(T_data$MCLNum)
test <- c()
i <- 1
for(i in i:length(temp)){
  test <- c(test, sum(T_data$MCLNum == i))
}
T_data <- T_data %>% filter(MCLNum <= length(test[test >= 3]), TF != "No")
#Transのデータ整理------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#各クラスターの転写因子ファミリーをリスト化したが、今回はシス側で着目している転写因子ファミリーのクラスターしか必要ないのでシスでエンリッチしたファミリーがTransでもあるのかgrepで確認
#grep("A", list(c("A", "B", "C", "A"), c("A", "B"), c("A", "C")))
#[1] 1 2 3
#list型でgrepは複数ヒットしてもその階層にあるかないかしか返さない
CisName <- TFname
Trans <- list()
i <- 1
for(i in i:length(CisName)){
  Trans <- c(Trans, list(T_data %>% filter(TF == CisName[i]) %>% select(AGI, MCLNum, TF, symbol, annotation)))
  print(i)
  i <- i+1
}
names(Trans) <- CisName
#Trans-Cis---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#各クラスターの転写因子のファミリー(Trans)と各クラスターでエンリッチしたシス配列の転写因子ファミリーが一致していたとき結合する可能性があるという考えの基TransとCis側の紐付け
#Source側は当然Targetの数だけ必要なのでrep
#Target側はひたすら回す
allTransList <- list()
allCisList <- list()
m <- 1
for(m in m:length(Trans)){
  TransList <- c()
  CisList <- c()
  if(nrow(Trans[[m]]) != 0){
    n <- 1
    for(n in n:nrow(Trans[[m]])){
      TransList <- c(TransList, rep(Trans[[m]]$AGI[n], each = length(allCis[[m]])))
      CisList <- c(CisList, allCis[[m]])
      n <- n+1
    }
  }
  allTransList <- c(allTransList, list(TransList))
  names(allTransList[[m]]) <- rep(TFname[m], times = length(allTransList[[m]]))
  allCisList <- c(allCisList, list(CisList))
  names(allCisList[[m]]) <- rep(TFname[m], times = length(allCisList[[m]]))
  m <- m+1
}
TransCisNetwork <- data.frame(Trans = unlist(allTransList),
                              value = rep(1, times = length(unlist(allTransList))),
                              Cis = unlist(allCisList),
                              stringsAsFactors = F
)
T.AGI <- unique(TransCisNetwork$Trans)
names(T.AGI) <- paste0("MCLNum", MasterTable$MCLNum[match(T.AGI, MasterTable$AGI)])
temp <- unique(names(T.AGI))
T.MCLNum <- c()
TT.AGI <- c()
n <- 1
for(n in n:length(temp)){
  T.MCLNum <- c(T.MCLNum, rep(temp[n], times = sum(temp[n] == names(T.AGI))))
  TT.AGI <- c(TT.AGI, T.AGI[names(T.AGI) == temp[n]])
}

TransCisNetwork <- rbind(TransCisNetwork, data.frame(Trans = T.MCLNum, 
                                                     value = rep(1, times = length(T.MCLNum)), 
                                                     Cis = TT.AGI, 
                                                     stringsAsFactors = F)
)

#Transのattributeファイル作成----------------------------------------------------------------------------------------------------------------------------------------------------------------
attrTrans <- T_data %>% filter(TF == "bHLH" | TF == "WRKY" | TF == "bZIP" | TF == "NAC" | TF == "ERF" | TF == "RAV") %>% select(AGI, symbol, TF, MCLNum, annotation)
#Trans-Cisネットワーク, Cytoscape format--------------------------------------------------------------------------------------------------------------------------------------------------
write.table(TransCisNetwork, "bigdata/yasue/tGRN_Groping/inflation4/TransCis_network/TransCisNetwork.txt", sep = "\t", quote = F, row.names = F)
write.table(attrTrans, "bigdata/yasue/tGRN_Groping/inflation4/TransCis_network/attrTrans.txt", sep = "\t", quote = F, row.names = F)

#"~/Nakano_RNAseq/network_analysis/script/tGRN_analysis20181023/diffExpTransCisNetwork.R"
temp <- unique(CY15_ExpNetwork$sourcetime)
temp2 <- unique(CY15_ExpNetwork$targettime)
condition <- c("CY15", "CY16", "CY20")
time <- c("1h", "3h", "12h")
CY <- list(CY15_ExpNetwork, CY16_ExpNetwork, CY20_ExpNetwork)
TransCis_pair <- paste0(TransCisNetwork$Trans, TransCisNetwork$Cis)
e <- 1
for(e in e:length(CY)){
  i <- 1
  for(i in i:length(temp)){
    T_data <- CY[[e]][CY[[e]][, "sourcetime"] == temp[i] & CY[[e]][, "targettime"] == temp2[i], ]
    sink <- c()
    k <- 1
    for(k in k:length(TransCis_pair)){
      sink <- c(sink, grep(TransCis_pair[k], paste0(T_data$source, T_data$target)))
    }
    T_data <- T_data[sink, ]
    T_data <- T_data[!duplicated(T_data), ]
    
    assign(paste0(condition[e], "_", time[i], "_ExpNetwork"), T_data)
    title <- paste0("bigdata/yasue/tGRN_Groping/inflation4/", condition[e], "/diffExpNetwork/", time[i], "/", condition[e], "_", time[i], "_ExpNetwork.txt")
    write.table(T_data, title, sep = "\t", quote = F, row.names = F)
    i <- i+1
  }
  e <- e+1
}

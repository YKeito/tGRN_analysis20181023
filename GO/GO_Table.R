#"~/Nakano_RNAseq/network_analysis/script/tGRN_analysis20181023/GO/GO_Table.R"
t <- proc.time()
i <- 1
for(i in i:length(T_MCLNum)){
  T_GOID <- allMCL_enrichment_score[match(paste0("MCLNum", T_MCLNum[i]), names(allMCL_enrichment_score))]
  T_AGI <- AGI_list[[match(paste0("MCLNum", T_MCLNum[i]), names(AGI_list))]]
  T_GO_term <- c()
  n <- 1
  for(n in n:length(names(T_GOID[[1]]))){
    T_GO_term <- c(T_GO_term, unique(GO_TAIR$`GO term`[GO_TAIR$`GO ID` == names(T_GOID[[1]])[n]]))
    n <- n+1
  }
  GO_Table <- data.frame(GO_ID = names(T_GOID[[1]]),
                         pvalue = T_GOID[[1]],
                         qvalue = p.adjust(T_GOID[[1]], method = "BH"),
                         GO_term = T_GO_term,
                         T_AGI = paste(T_AGI, collapse = "|"),
                         stringsAsFactors = F
  )
  title <- paste0("bigdata/yasue/tGRN_Groping/inflation4/GO/", "MCLNum", T_MCLNum[i], "_", "GO_Table.txt")
  write.table(GO_Table, file = title, sep = "\t", quote = F, row.names = F)
  print(i)
  i <- i+1
}
t1 <- proc.time() - t
print(t1)
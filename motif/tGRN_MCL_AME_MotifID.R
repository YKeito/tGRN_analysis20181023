#"~/Nakano_RNAseq/network_analysis/script/tGRN_analysis20181023/tGRN_MCL_AME_MotifID.R"
T_filename <- list.files("~/bigdata/yasue/tGRN_Groping/inflation4/motif/motif_results", full.names = T)
filename <- c()
Motif_ID <- c()
obnames <- c()
consensus <- c()
motif_alt_ID <- c()
i <- 1
for(i in i:length(T_MCLNum)){
  filename <- T_filename[grep(paste0("MCLNum", T_MCLNum[i], "_"), T_filename)]
  title <- paste0(filename, "/ame.tsv")
  motif_ame_results <- try(read.table(title, sep = "\t", header = T, stringsAsFactors = F), silent = TRUE)
  if(class(motif_ame_results) != "try-error"){
    Motif_ID <- c(Motif_ID, motif_ame_results$motif_ID)
    consensus <- c(consensus, motif_ame_results$consensus)
    motif_alt_ID <- c(motif_alt_ID, motif_ame_results$motif_alt_ID)
    obnames <- c(obnames, rep(paste0("MCLNum", T_MCLNum[i]), times = length(motif_ame_results$motif_ID)))
  }
  i <- i+1
}

tGRN_MotifID_consensus <- data.frame(MCLNum = obnames,
                                     MotifID = Motif_ID,
                                     consensu = consensus,
                                     motif_alt_ID = motif_alt_ID,
                                     stringsAsFactors = F
                                     )

write.table(tGRN_MotifID_consensus, file = "bigdata/yasue/tGRN_Groping/inflation4/inflation4_MotifID_consensus.txt", append=F, quote = F, sep = "\t", row.names = F)
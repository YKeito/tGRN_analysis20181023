#"~/Nakano_RNAseq/network_analysis/script/tGRN_analysis20181023/motif/tGRN_MCL_AME_MotifID.R"
T_filename <- list.files("~/bigdata/yasue/tGRN_Groping/inflation4/motif/motif_results", full.names = T)
filename <- c()
Motif_ID <- c()
obnames <- c()
consensus <- c()
motif_alt_ID <- c()
i <- 1
for(i in i:length(T_MCLNum)){#T_MCLNum:"~/Nakano_RNAseq/network_analysis/script/tGRN_analysis20181023/SubCluster_Over3Nodes.R"
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
#AGRIS(https://agris-knowledgebase.org/AtcisDB/bindingsites.html)
#WRKY:TTGAC
#bHLH(G-box promoter motif):CACGTG
#ERF(GCC-box promoter motif):GCCGCC
#MYB binding site promoter:(A/C)ACC(A/T)A(A/C)C
#reference:Regulating the Regulators: The Control of Transcription Factors in Plant Defense Signaling
#bZIP:ACGT
#NAC:CGT[G/T], ACG
#RAV:CAACA
MotifList <- c("TTGAC", "GTCAA", "CACGTG", "GCCGCC", "GGCGGC", "ACGT", "CGT", "ACG", "CAACA", "TGTTG")
names(MotifList) <- c("WRKY", "WRKY", "bHLH", "ERF", "ERF", "bZIP", "NAC", "NAC", "RAV", "RAV")
CoreMotif <- rep(0, times = nrow(tGRN_MotifID_consensus))
TF <- rep(0, times = nrow(tGRN_MotifID_consensus))
i <- 1
for(i in i:length(MotifList)){
  if(length(CoreMotif[grep(MotifList[i], tGRN_MotifID_consensus$consensu)]) != 0){
    n <- 1
    for(n in n:length(CoreMotif[grep(MotifList[i], tGRN_MotifID_consensus$consensu)])){
      if(!is.na(CoreMotif[grep(MotifList[i], tGRN_MotifID_consensus$consensu)][n])){
        if(CoreMotif[grep(MotifList[i], tGRN_MotifID_consensus$consensu)][n] ==  0){
          CoreMotif[grep(MotifList[i], tGRN_MotifID_consensus$consensu)][n] <- MotifList[i]
          TF[grep(MotifList[i], tGRN_MotifID_consensus$consensu)][n] <- names(MotifList)[i]
        }else{
          CoreMotif[grep(MotifList[i], tGRN_MotifID_consensus$consensu)][n] <- paste0(CoreMotif[grep(MotifList[i], tGRN_MotifID_consensus$consensu)][n], "|", MotifList[i])
          TF[grep(MotifList[i], tGRN_MotifID_consensus$consensu)][n] <- paste0(TF[grep(MotifList[i], tGRN_MotifID_consensus$consensu)][n], "|", names(MotifList)[i])
        }
      }
      n <- n+1
    }
  }
  i <- i+1
}
tGRN_MotifID_consensus <-tGRN_MotifID_consensus %>% mutate(CoreMotif = CoreMotif,
                                                           TF = TF)
write.table(tGRN_MotifID_consensus, file = "bigdata/yasue/tGRN_Groping/inflation4/motif/inflation4_MotifID_consensus.txt", append=F, quote = F, sep = "\t", row.names = F)
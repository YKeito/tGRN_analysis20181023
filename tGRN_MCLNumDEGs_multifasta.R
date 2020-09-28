#"~/Nakano_RNAseq/network_analysis/script/tGRN_analysis20181023/motif/tGRN_MCLNumDEGs_multifasta.R"
###########################################
#upstream_500 <- read.table(file = "~/bigdata/yasue/motif/TAIR10_upstream_500_20101028", fill = T, sep = ",", stringsAsFactors = F)
#upstream_500 <- as.character(unlist(upstream_500))
temp <- grep("chr", upstream_500)


AGI <- c()
i <- 1
total <- length(temp)
for(i in i:total){
  AGI <- c(AGI, substr(upstream_500[temp[i]], 2, 10))
  print(i)
  i <- i+1
}

#最後の一周だけ自動化できなかった。
allsequence <- list()
presequence <- c()
sequence <- c()
i <- 1
total <- length(temp)
for(i in i:c(total-1)){
  presequence <- upstream_500[c(temp[i]+1):c(temp[i+1]-1)]
  n <- 1
  for(n in n:length(presequence)){
    sequence <- paste0(sequence, presequence[n])
    n <- n+1
  }
  allsequence <- c(allsequence, list(sequence))
  sequence <- c()
  print(i)
  i <- i+1
}
#残りの最後の一周を追加
presequence <- upstream_500[c(temp[total]+1):length(upstream_500)]
sequence <- paste0(sequence, presequence[n])
allsequence <- c(allsequence, list(sequence))

#data.frame
up_500bp <- data.frame(AGI = AGI,
                       sequence = unlist(allsequence),
                       stringsAsFactors = F
)
#####################対象の遺伝子群を引っ張ってくる#######################
load("~/bigdata/yasue/tGRN_Groping/inflation4/MasterTable_inflation4.RData")
t <- proc.time()
library("stringr")

temp <- MasterTable[!is.na(MasterTable$MCLNum), ]
T_AGI <- c()
T_data <- c()
T_control <- c()
T_control_data <- c()
total_i <- length(unique(temp$MCLNum[!is.na(temp$MCLNum)]))
i <- 1
for(i in i:total_i){
  T_AGI <- filter(temp, MCLNum == i)[, "AGI"]
  T_data <- up_500bp[match(T_AGI, up_500bp$AGI), ]
  T_control <- sample(up_500bp$AGI, length(T_AGI))
  T_control_data <- up_500bp[match(T_control, up_500bp$AGI), ]
  
  fasta <- c()
  allfasta <- c()
  cont_fasta <- c()
  cont_allfasta <- c()
  total_o <- nrow(T_data)
  o <- 1
  for(o in o:total_o){
    fasta <- rbind(str_sub(T_data[, "sequence"][o], start=1, end=80),
                   str_sub(T_data[, "sequence"][o], start=81, end=160),
                   str_sub(T_data[, "sequence"][o], start=161, end=240),
                   str_sub(T_data[, "sequence"][o], start=241, end=320),
                   str_sub(T_data[, "sequence"][o], start=321, end=400),
                   str_sub(T_data[, "sequence"][o], start=401, end=480),
                   str_sub(T_data[, "sequence"][o], start=481, end=500)
    )
    cont_fasta <- rbind(str_sub(T_control_data[, "sequence"][o], start=1, end=80),
                        str_sub(T_control_data[, "sequence"][o], start=81, end=160),
                        str_sub(T_control_data[, "sequence"][o], start=161, end=240),
                        str_sub(T_control_data[, "sequence"][o], start=241, end=320),
                        str_sub(T_control_data[, "sequence"][o], start=321, end=400),
                        str_sub(T_control_data[, "sequence"][o], start=401, end=480),
                        str_sub(T_control_data[, "sequence"][o], start=481, end=500)
    )
    data_fastaAGI <- paste0(">", T_data[, "AGI"][o])
    control_fastaAGI <- paste0(">", T_control_data[, "AGI"][o])
    allfasta <- c(allfasta, rbind(data_fastaAGI, fasta))
    cont_allfasta <- c(cont_allfasta, rbind(control_fastaAGI, cont_fasta))
    o <- o+1
  }
  target <- paste0("~/bigdata/yasue/tGRN_Groping/inflation4/motif/multi-fasta/target/", "tGRN_MCLNum", i, "_upstream500.fasta")
  control <- paste0("~/bigdata/yasue/tGRN_Groping/inflation4/motif/multi-fasta/control/", "control_tGRN_MCLNum", i, "_upstream500.fasta")
  write.table(allfasta, file = target, append = F, quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(cont_allfasta, file = control, append = F, quote = F, sep = "\t", row.names = F, col.names = F)
  print(i)
  i <- i+1
}

t1 <- proc.time() - t
print(t1)
#"~/Nakano_RNAseq/network_analysis/script/tGRN_analysis20181023/motif/tGRN_Subcluster_enrichmentTF_heatmap.R"
#どのクラスターにどの転写因子ファミリーのシス配列がエンリッチしたかヒートマップ化
tGRN_MotifID_consensus$TF[tGRN_MotifID_consensus$TF == 0] <- "others"
TFname <- c(unique(names(MotifList)), "others")
allMCLNum <- list()
k <- 1
for(k in k:length(TFname)){
  T_MotifID <- tGRN_MotifIDTF$MotifID[grep(TFname[k], tGRN_MotifID_consensus$TF)]
  T_MCLNum <- c()
  j <- 1
  for(j in j:length(T_MotifID)){
    T_MCLNum <- c(T_MCLNum, tGRN_MotifID_consensus$MCLNum[T_MotifID[j] == tGRN_MotifID_consensus$MotifID])
    j <- j+1
  }
  allMCLNum <- c(allMCLNum, list(T_MCLNum))
  k <- k+1
}
names(allMCLNum) <- TFname

UniMCLNum <- unique(tGRN_MotifID_consensus$MCLNum)
count <- list()
i <- 1
for(i in i:length(UniMCLNum)){
  total <- c()
  k <- 1
  for(k in k:length(allMCLNum)){
    total <- c(total, sum(allMCLNum[[k]] == UniMCLNum[i])/sum(unlist(allMCLNum) == UniMCLNum[i]))
    k <- k+1
  }
  names(total) <- TFname
  count <- c(count, list(total))
  i <- i+1
}
names(count) <- UniMCLNum
temp <- str_sub(UniMCLNum, start = 7, end = 9)

T_data <- data.frame(MCLNum = rep(paste0(str_sub("0000", start = 1, end  = nchar("000")-nchar(temp)), temp), each = length(TFname)),
                     value = as.vector(unlist(count)),
                     TF = rep(paste0("0", length(TFname):1, TFname), times = length(UniMCLNum)),
                     stringsAsFactors = F
)
library(ggplot2)
g <- ggplot(T_data, aes(x = TF, y = MCLNum, fill = value))
g <- g + geom_tile(color = "black", size = 0.1)
g <- g + scale_fill_gradient2(high = "red")
g <- g+theme_dark()
g <- g + theme_linedraw()
g <- g + coord_flip()
g <- g +　theme(axis.text=element_text(size=15))
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 0))
plot(g)
title <- paste0("bigdata/yasue/tGRN_Groping/inflation4/results_Fig/motif/", paste(str_split(Sys.Date(), pattern = "-", simplify = T), collapse = ""), "tGRN_Subcluster_MotifID_enrichment_heatmap.png")
ggsave(title, g, width = 8, height = 5)
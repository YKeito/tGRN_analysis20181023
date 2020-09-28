"~/Nakano_RNAseq/network_analysis/script/tGRN_analysis20181023/Fig/tGRN_MCL_AME_MotifID_heatmap.R"
T_data <- data.frame(sample = rep(unique(tGRN_MotifID_consensus$MCLNum), each = length(unique(tGRN_MotifID_consensus$MotifID))),
                     MotifID = rep(unique(tGRN_MotifID_consensus$MotifID), times = length(unique(tGRN_MotifID_consensus$MCLNum))),
                     value = rep(0, times = length(unique(tGRN_MotifID_consensus$MCLNum))*length(unique(tGRN_MotifID_consensus$MotifID))),
                     stringsAsFactors = F
)


target <- paste0(tGRN_MotifID_consensus$MCLNum, tGRN_MotifID_consensus$MotifID)
i <- 1
for(i in i:length(target)){
  T_data$value[which(paste0(T_data$sample, T_data$MotifID) == target[i])] <- 1
  i <- i+1
}

sample_sort <- c()
n <- 1
for(n in n:length(unique(T_data$sample))){
  sample_sort <- c(sample_sort, rep(c(length(unique(T_data$sample))+1-n), times = sum(T_data$sample == unique(T_data$sample)[n])))
}

df <- data.frame(T_data, 
                 sample_sort,
                 stringsAsFactors = F
)

g <- ggplot(df, aes(x = MotifID, y = reorder(sample, sample_sort), fill = value))
g <- g + geom_tile(color = "black", size = 0.1)
g <- g + scale_fill_gradient2(high = "red")
g <- g+theme_dark()
g <- g + theme_linedraw()
#g <- g + coord_flip()
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot(g)
ggsave("bigdata/yasue/tGRN_Groping/inflation4/results_Fig/motif/inflation4_MotifID_heatmap.png", g, width = 20, height = 10)
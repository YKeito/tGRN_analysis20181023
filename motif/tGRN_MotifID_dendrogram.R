#"~/Nakano_RNAseq/network_analysis/script/tGRN_analysis20181023/tGRN_MotifID_dendrogram.R"
library(stringdist)
uni_consensus <- unique(tGRN_MotifID_consensus$consensu)
uni_motif_alt_ID <- unique(tGRN_MotifID_consensus$motif_alt_ID)
names(uni_motif_alt_ID) <- unique(tGRN_MotifID_consensus$MotifID)
names(uni_consensus) <- unique(tGRN_MotifID_consensus$MotifID)

LCS_Table <- c()
n <- 1
for(n in n:length(uni_consensus)){
  test <- c()
  m <- 1
  for(m in m:length(uni_consensus)){
    test <- c(test, stringdist(uni_consensus[n], uni_consensus[m],  method = "lcs"))
    m <- m+1
  }
  LCS_Table <- cbind(LCS_Table, test)
  n <- n+1
}
rownames(LCS_Table) <- names(uni_consensus)
colnames(LCS_Table) <- names(uni_consensus)


#install.packages("factoextra")
#install.packages("dendextend")
#library(factoextra)
#library(ggplot2)
#library(igraph)

#k:クラスター数の指定
#k_colors, palette:色指定
#show_labels:横軸の名前を表示するか
#color_labels_by_k:クラスターのグループに応じて色を付けるか
#horiz縦軸と横軸を反転するかしないか
#rect:クラスターを囲む
#rect_fill:rectで囲んだ範囲を色塗する
#type:図の形式の指定

res.hc <- eclust(x = LCS_Table,
                 "hclust",
                 k = 5,
                 method = "euclidean",
                 graph = FALSE
)

g <- fviz_dend(res.hc,
               cex = 1.2,
               color_labels_by_k = TRUE,
               show_labels = TRUE,
               ggtheme = theme_bw(),
               horiz = FALSE,
               rect = TRUE,
               rect_fill = TRUE,
               type = "rectangle"
)
plot(g)

ggsave("bigdata/yasue/tGRN_Groping/inflation4/results_Fig/inflation4_MotifID_clustering.png", g, width = 20, height = 10)

T_data <- data.frame(ClusterNum = res.hc$cluster,
                     consensus = uni_consensus[match(names(res.hc$cluster), names(uni_consensus))],
                     motif_alt_ID = uni_motif_alt_ID[match(names(res.hc$cluster), names(uni_motif_alt_ID))],
                     stringsAsFactors = F
                     )
write.table(T_data, "bigdata/yasue/tGRN_Groping/inflation4/inflation4_Subcluster5_MotifID_summary.txt", sep = "\t", quote = F)
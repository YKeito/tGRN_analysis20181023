"~/Nakano_RNAseq/network_analysis/script/tGRN_analysis20181023/Fig/Comparison_CYtGRN.R"
#CY化合物毎のtGRNの構造を比較する
temp <- list(CY15_ExpNetwork, CY16_ExpNetwork, CY20_ExpNetwork)
object <- c("CY15_Edge", "CY16_Edge", "CY20_Edge")
CY_Edge <- list()
n <- 1
for(n in n:length(temp)){
  T_data <- temp[[n]]
  T_edge <- c()
  m <- 1
  for(m in m:nrow(T_data)){
    T_edge <- c(T_edge, paste(T_data[m, ], collapse = ","))
    m <- m+1
  }
  CY_Edge <- c(CY_Edge, list(T_edge))
  assign(object[n], T_edge)
  n <- n+1
}
names(CY_Edge) <- object

library(gplots)
CYall_Venn <- venn(CY_Edge)
names(attr(CYall_Venn,"intersections")) <- paste0(names(attr(CYall_Venn,"intersections")), ",")
T.data <- data.frame(CY = str_split(names(unlist(attr(CYall_Venn,"intersections"))), pattern = ",", simplify = T)[, 1],
                     str_split(unlist(attr(CYall_Venn,"intersections"), use.names = F), pattern = ",", simplify = T),
                     stringsAsFactors = F)
colnames(T.data) <- c("CY", "Source", "S_Time", "Target", "T_Time", "regulate", "attredge")
write.table(T.data, file = "~/bigdata/yasue/tGRN_Groping/inflation4/CYall_Venn.txt", sep = "\t", quote = F, row.names = F)
ggsave(filename = "~/bigdata/yasue/tGRN_Groping/inflation4/results_Fig/tGRN_Structure/CY_Edge.png", plot = plot(CYall_Venn), width = 6, height = 4)
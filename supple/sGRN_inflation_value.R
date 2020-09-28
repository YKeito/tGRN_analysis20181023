#"~/Nakano_RNAseq/network_analysis/script/tGRN_analysis20181023/supple/sGRN_inflation_value.R"
filename <- list.files("bigdata/yasue/tGRN_Groping/NodeTable", pattern = ".csv", full.names = T)
object_name <- c("inflation2.5", "inflation3", "inflation4", "inflation5")
n <- 1
for (n in 1:length(filename)){
  assign(object_name[n], read.table(filename[n],  header=T, sep=",", stringsAsFactors = F))
  n <- n+1
}
temp <- list(inflation2.5, inflation3, inflation4, inflation5)
MCLNum <- c()
Category <- c()
library(ggplot2)
n <- 1
for (n in 1:length(filename)){
  assign(object_name[n], read.table(filename[n],  header=T, sep=",", stringsAsFactors = F))
  MCLNum <- c(temp[[n]][, "X__mclCluster"][temp[[n]][, "X__mclCluster"] < 6][!is.na(temp[[n]][, "X__mclCluster"][temp[[n]][, "X__mclCluster"] < 6])])
  Category <- rep(object_name[n], times = length(MCLNum[[n]]))
  T_data <- data.frame(MCLNum = unlist(MCLNum),
                       Category = Category
  )
  g <- ggplot(T_data, aes(x = MCLNum))
  g <- g + geom_histogram(binwidth=0.5)
  g <- g + ggtitle(object_name[n])
  plot(g)
  title <- paste0("bigdata/yasue/tGRN_Groping/supple/", object_name[n], ".png")
  #ggsave(title, g)
  n <- n+1
}

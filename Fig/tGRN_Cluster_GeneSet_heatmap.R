#"~/Nakano_RNAseq/network_analysis/script/tGRN_analysis20181023/Fig/tGRN_Cluster_GeneSet_heatmap.R"
library(dplyr)
library(stringr)
library(ggplot2)
T_data <- MasterTable[!is.na(MasterTable$MCLNum), ]
temp <- unique(T_data$MCLNum)
test <- c()
test2 <- c()
i <- 1
for(i in i:length(temp)){
  test <- c(test, sum(T_data$MCLNum == i))
  test2 <- c(test2, sum(T_data$TF[T_data$MCLNum == i] != "No"))
}

T_data <- T_data[T_data$MCLNum <= length(test[test >= 3]), ]
T_data <- T_data[T_data$TF != "No", ]

MCL <- unique(T_data$MCLNum)
temp <- union(TransCisNetwork$Trans[grep("MCLNum", TransCisNetwork$Trans)], TransCisNetwork$Cis[grep("MCLNum", TransCisNetwork$Cis)])
temp <- str_split(temp, pattern = "Num", simplify = T)[, 2]
MCL <- intersect(MCL, temp)

sample <- c("Botrytis_cinerea", "PstDC3000_hrp", "MeJA_DEGs", "BTH")
allpvalue <- list()
k <- 1
for(k in k:length(sample)){
  N <- nrow(MasterTable)
  M <- sum(MasterTable[, sample[k]] == "Yes")
  pvalue <- c()
  T_MCLNum <- c()
  temp <- c()
  i <- 1
  for(i in i:length(MCL)){
    T_MCLNum <- c(T_MCLNum, MCL[i])
    temp <- c(temp, paste0(substr("0000", start = 1, stop = nchar("000")-nchar(T_MCLNum[i])), T_MCLNum[i]))
    T_data <- filter(MasterTable, MCLNum == T_MCLNum[i])
    n <- nrow(T_data)
    x <- sum(T_data[, sample[k]] == "Yes")
    pvalue <- c(pvalue, phyper(x-1, M, N-M, n, lower.tail = F))
    i <- i+1
  }
  names(pvalue) <- MCL
  allpvalue <- c(allpvalue, list(pvalue))
  k <- k+1
}

names(allpvalue) <- sample
T_data <- data.frame(MCLNum = rep(temp, times = 4),
                     pvalue = unlist(allpvalue),
                     enrichment = -log2(unlist(allpvalue)),
                     sample = rep(c("01B.cinerea", "02PstDC3000", "03MeJA", "04BTH"), each = length(MCL))
)
rownames(T_data) <- c()
T_data$pvalue[T_data$pvalue > 5e-2] <- 100
T_data$pvalue[T_data$pvalue <= 5e-2 & T_data$pvalue > 5e-5] <- 75
T_data$pvalue[T_data$pvalue <= 5e-5 & T_data$pvalue > 5e-10] <- 50
T_data$pvalue[T_data$pvalue <= 5e-10] <- 0



library(ggplot2)
g <- ggplot(data = T_data, aes(x = MCLNum, y = sample, fill = pvalue))
g <- g + geom_tile(color = "black", size = 0.5)
g <- g + scale_fill_gradient(low="seagreen",high="white")
g <- g + theme_linedraw()
g <- g + theme_set(theme_bw(base_size = 20))
g <- g +　theme(axis.text=element_text(size=24))
g <- g + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
g <- g + theme(legend.position="top")
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot(g)
ggsave("~/bigdata/yasue/tGRN_Groping/inflation4/results_Fig/GeneSet_heatmap.png", g, width = 24, height = 8)
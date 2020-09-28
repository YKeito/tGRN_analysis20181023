#"~/Nakano_RNAseq/network_analysis/script/tGRN_analysis20181023/Fig/tGRN_MCL_Node_graph.R"
library(ggplot2)
MCLNum <- union(TransCisNetwork$Trans, TransCisNetwork$Cis)
MCLNum <- MCLNum[grep("MCLNum", MCLNum)]
condition <- c("CY15", "CY16", "CY20")
e <- 1
for(e in e:length(condition)){
  all <- c()
  n <- 1
  for(n in n:length(MCLNum)){
    temp <- str_split(MCLNum, pattern = "MCLNum")[[n]][2]
    T_AGI <- MasterTable$AGI[MasterTable$MCLNum == temp]
    T_AGI <- T_AGI[!is.na(T_AGI)]
    T_data <- allRNASeq[match(T_AGI, rownames(allRNASeq)), ]
    T_data <- T_data[order(T_data[, paste0(condition[e], "_1h")], decreasing = T), ]
    sort <- 1:nrow(T_data)
    T_expression <- c(rep(0, times = nrow(T_data)),
                      T_data[, paste0(condition[e], "_1h")],
                      T_data[, paste0(condition[e], "_3h")],
                      T_data[, paste0(condition[e], "_12h")],
                      T_data[, paste0(condition[e], "_24h")]
    )
    
    df <- data.frame(AGI = rep(rownames(T_data), time = 5),
                     expression = T_expression,
                     Category = rep(c("00h", "01h", "03h", "12h", "24h"), each = length(T_AGI)),
                     time = rep(c(0, 1, 3, 12, 24), each = length(T_AGI)),
                     group = rep("", times = length(T_expression)),
                     AGI_sort = rep(sort, time = 5),
                     MCLNum = rep(paste0(substr("0000", start = 1, stop = nchar("000")-nchar(temp)), temp), times = nrow(T_data)*5),
                     stringsAsFactors = F
    )
    
    
    all <- rbind(all, df)
    ####geom_smooth####
    g <- ggplot(data = df, aes(x = time, y = expression))
    g <- g + geom_smooth(method = "loess", mapping = aes(x = time, y = expression), level = 0.95)
    g <- g + theme(legend.position="top")
    g <- g + ggtitle(paste0(condition[e], "_", MCLNum[n], "_NumNode:", length(T_AGI)))
    g <- g + theme_set(theme_bw(base_size = 10))
    g <- g +　theme(axis.text=element_text(size=18))
    plot(g)
    title <- paste0("bigdata/yasue/tGRN_Groping/inflation4/results_Fig/expression/smooth/", condition[e], "/", MCLNum[n], "_smooth", ".png")
    ggsave(title, g, width = 6, height = 4)
    
    ####geom_histo####
    T_data <- df
    T_data$expression[T_data$expression < -2] <- -2
    T_data$expression[T_data$expression > 2] <- 2
    g <- ggplot(data = T_data, aes(x = Category, y = reorder(AGI, AGI_sort), fill = expression))
    g <- g + geom_tile()
    g <- g + theme(legend.position="top")
    g <- g + scale_fill_gradient2(low = "blue", high = "red", midpoint = 0, limit = c(-2, 2))
    g <- g + ggtitle(paste0(condition[e], "_", MCLNum[n], "_NumNode:", length(T_AGI)))
    plot(g)
    title <- paste0("bigdata/yasue/tGRN_Groping/inflation4/results_Fig/expression/heatmap/", condition[e], "/", MCLNum[n], "_heatmap", ".png")
    ggsave(title, g, width = 8, height = 4)
    print(n)
    n <- n+1
  }
  
  g <- ggplot(data = all, aes(x = time, y = expression))
  g <- g + geom_smooth(method = "loess", mapping = aes(x = time, y = expression), level = 0.95)
  g <- g + facet_wrap(MCLNum ~ ., nrow=4, ncol = 8, scales = "free")
  #g <- g + theme_set(theme_bw(base_size = 12))
  g <- g + theme_set(theme_bw(base_size = 28))
  g <- g + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  #plot(g)
  title <- paste0("bigdata/yasue/tGRN_Groping/inflation4/results_Fig/expression/smooth/", condition[e], "_smooth.png")
  #ggsave(title, g, width = 12, height = 4, limitsize = F)
  ggsave(title, g, width = 24, height = 9, limitsize = F)
  e <- e+1
}

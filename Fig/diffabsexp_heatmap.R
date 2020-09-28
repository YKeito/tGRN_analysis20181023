#"~/Nakano_RNAseq/network_analysis/script/tGRN_analysis20181023/Fig/diffabsexp_heatmap.R"
condition <- c("CY15", "CY16", "CY20")
MCLNum <- attrTransCis$MCLNum
time <- c("1h", "3h", "12h")
e <- 1
for(e in e:length(condition)){
  diff <- c()
  T_MCLNum <- c()
  peak <- c()
  n <- 1
  for(n in n:length(MCLNum)){
    temp <- str_split(MCLNum, pattern = "MCLNum")[[n]][2]
    T_AGI <- filter(MasterTable, MCLNum == temp)[, "AGI"]
    T_data <- allRNASeq[match(T_AGI, rownames(allRNASeq)), ]
    AbsExpression <- c(mean(abs(T_data[, paste0(condition[e], "_1h")])),
                       mean(abs(T_data[, paste0(condition[e], "_3h")])),
                       mean(abs(T_data[, paste0(condition[e], "_12h")])),
                       mean(abs(T_data[, paste0(condition[e], "_24h")]))
    )
    diff00h_01h <- AbsExpression[1]-0
    diff03h_01h <- AbsExpression[2]-AbsExpression[1]
    diff12h_03h <- AbsExpression[3]-AbsExpression[2]
    diff24h_12h <- AbsExpression[4]-AbsExpression[3]
    diff00h_01h_sign <- sign(diff1h_3h)
    diff03h_01h_sign <- sign(diff03h_01h)
    diff12h_03h_sign <- sign(diff12h_03h)
    diff24h_12h_sign <- sign(diff24h_12h)
    peak <- c(peak, time[which.max(c(AbsExpression[1], AbsExpression[2], AbsExpression[3], AbsExpression[4]))])
    names(peak)[n] <- MCLNum[n]
    diff <- rbind(diff, cbind(diff00h_01h, diff03h_01h, diff12h_03h, diff24h_12h, diff00h_01h_sign, diff03h_01h_sign, diff12h_03h_sign, diff24h_12h_sign))
    T_MCLNum <- c(T_MCLNum, temp)
    print(n)
    n <- n+1
  }
  
  T_data <- data.frame(diffabs = c(diff[, 1], diff[, 2], diff[, 3], diff[, 4]),
                       Sample = rep(paste0("MCLNum", T_MCLNum), times = 4),
                       condition = rep(colnames(diff)[1:4], each = length(T_MCLNum)),
                       stringsAsFactors = F
  )
  
  attrTransCis_time <- data.frame(MCLNum = names(peak),
                                  peaktime = peak,
                                  stringsAsFactors = F
  )
  rownames(attrTransCis_time) <- c()
  title <- paste0("bigdata/yasue/tGRN_Groping/inflation4/TransCis_network/", condition[e], "/", condition[e], "_attrTransCis_time.txt")
  write.table(attrTransCis_time, title, sep = "\t", quote = F, row.names = F)
  
  library(ggplot2)
  g1 <- ggplot(T_data, aes(x = condition, y = Sample, fill = diffabs))
  g1 <- g1 + geom_tile()
  g1 <- g1 + theme_bw()
  g1 <- g1 + scale_fill_gradient2(low = "blue", high = "red", na.value = "white")
  g1 <- g1 + ggtitle(condition[e])
  plot(g1)
  
  exp_pattern <- data.frame(Sample = paste0("MCLNum", T_MCLNum),
                        diff,
                        stringsAsFactors = F
  )
  peak_all <- list()
  k <- 1
  for(k in k:nrow(exp_pattern)){
    check <- c()
    j <- 6
    for(j in c(6:8)){
      check <- c(check, sum(exp_pattern[k, j], exp_pattern[k, j+1]))
      j <- j+1
    }
    names(check) <- c("1h", "3h", "12h")
    peak_all <- c(peak_all, list(check))
    k <- k+1
  }
  names(peak_all) <- exp_pattern$Sample
  
  
  peak_time <- c() 
  k <- 1
  for(k in k:length(peak_all)){
    temp2 <- paste0(names(peak_all)[k], "_", names(which(peak_all[[k]] == 0)))
    temp2 <- str_split(temp2, pattern = "_")
    peak_time <- c(peak_time, temp2[[1]][2])
    names(peak_time)[k] <- temp2[[1]][1]
    k <- k+1
  }
  peak_time <- peak_time[peak_time != ""]
  
  
  exp_pattern[, 6:9][exp_pattern[, 6:9] == 1] <- "p"
  exp_pattern[, 6:9][exp_pattern[, 6:9] == -1] <- "n"
  difftime <- paste0(c("0h_01h", "03h_01h", "12h_03h", "24h_12h"), "_sign")
  check_max <- apply(exp_pattern[, 2:5], MARGIN = 1, FUN = which.max)
  check_max <- difftime[check_max]
  check_min <- apply(exp_pattern[, 2:5], MARGIN = 1, FUN = which.min)
  check_min <- difftime[check_min]
  
  m <- 1
  for(m in m:nrow(exp_pattern)){
    if(sum(exp_pattern[m, 6:9] == "p") != 0){
      exp_pattern[m, grep(check_max[m], colnames(exp_pattern))] <- "pp"
    }
    if(sum(exp_pattern[m, 6:9] == "n") != 0){
      exp_pattern[m, grep(check_min[m], colnames(exp_pattern))] <- "nn"
    }
    m <- m+1
  }
  
  ExpTable <- c()
  o <- 1
  for(o in o:length(time)){
    target <- list()
    regulate <- list()
    T_time <- list()
    source <- list()
    source_time <- list()
    q <- o+6
    for(q in q:ncol(exp_pattern)){
      target <- c(target, list(exp_pattern$Sample[which(exp_pattern[, q] == "pp" | exp_pattern[, q] == "nn")]))
      regulate <- c(regulate, list(exp_pattern[, q][which(exp_pattern[, q] == "pp" | exp_pattern[, q] == "nn")]))
      T_time <- c(T_time, list(rep(colnames(exp_pattern)[q], length(exp_pattern[, q][which(exp_pattern[, q] == "pp" | exp_pattern[, q] == "nn")]))))
      
      source <- c(source, list(names(peak_time[peak_time == time[o]])))
      source_time <- c(source_time, list(peak_time[peak_time == time[o]]))
      q <- q+1
    }
    allsource <- list()
    alltarget <- list()
    allregulate <- list()
    alltime <- list()
    allsource_time <- list()
    q <- 1
    for(q in q:length(target)){
      allsource <- c(allsource, list(rep(source[[q]], each = length(target[[q]]))))
      alltarget <- c(alltarget, list(rep(target[[q]], times = length(source[[q]]))))
      allregulate <- c(allregulate, list(rep(regulate[[q]], times = length(source[[q]]))))
      alltime <- c(alltime, list(rep(T_time[[q]], times = length(source[[q]]))))
      allsource_time <- c(allsource_time, list(rep(source_time[[q]], each = length(target[[q]]))))
    }
    attrTarget <- attrTransCis$attr[match(unlist(alltarget), attrTransCis$MCLNum)]
    attrsource <- attrTransCis$attr[match(unlist(allsource), attrTransCis$MCLNum)]
    
    ExpTable <- rbind(ExpTable, data.frame(source = unlist(allsource),
                                           sourcetime = unlist(allsource_time),
                                           attrsource = attrsource,
                                           target = unlist(alltarget),
                                           targettime = unlist(alltime),
                                           attrtarget = attrTarget,
                                           regulate = unlist(allregulate),
                                           interaction = rep(1, times = length(unlist(allregulate))),
                                           stringsAsFactors = F
    )
    )
    o <- o+1
  }
  attr_ExpTable <- data.frame(ExpTable = c(unique(ExpTable$source), unique(ExpTable$target)),
                              attr = attrTransCis$attr[match(c(unique(ExpTable$source), unique(ExpTable$target)), attrTransCis$MCLNum)],
                              stringsAsFactors = F)
  
  assign(paste0(condition[e], "_ExpTable"), ExpTable)
  assign(paste0(condition[e], "_attrExpTable"), attr_ExpTable)
  
  assign(paste0(condition[e], "_ExpNetwork"), ExpTable[grep("Trans", ExpTable$attrsource), ])

  title <- paste0("bigdata/yasue/tGRN_Groping/inflation4/TransCis_network/", condition[e], "/", condition[e], "_diffabs_table.txt")
  write.table(exp_pattern, title, sep = "\t", quote = F, row.names = F)
  
  title <- paste0("bigdata/yasue/tGRN_Groping/inflation4/TransCis_network/", condition[e], "/", condition[e], "_ExpTable.txt")
  write.table(ExpTable, title, sep = "\t", quote = F, row.names = F)
  
  title <- paste0("bigdata/yasue/tGRN_Groping/inflation4/TransCis_network/", condition[e], "/", condition[e], "_attrExpTable.txt")
  write.table(attr_ExpTable, title, sep = "\t", quote = F, row.names = F)
  
  title <- paste0("bigdata/yasue/tGRN_Groping/inflation4/TransCis_network/", condition[e], "/", condition[e], "_ExpNetwork.txt")
  write.table(ExpTable[grep("Trans", ExpTable$attrsource), ], title, sep = "\t", quote = F, row.names = F)
}
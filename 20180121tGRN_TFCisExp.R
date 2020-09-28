"~/Nakano_RNAseq/network_analysis/script/tGRN_analysis20181023/20180121tGRN_TFCisExp.R.R"
#やりたいことメモ-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#目的：CY化合物毎に動的な遺伝子制御ネットワークを作成するための.txtファイルを出力させる
#Trans-Cisネットワークのエッジの確からしさを高めるために、RNA-Seqで得られた実際の発現パターンを利用
#Trans側の発現ピークの後にCis側の劇的な発現変化が起きていればエッジを引くといった考え
#逆に発現のピークがずれているものや劇的に変化していないものはエッジは引かない。
#Cis-Trans Networkを作成したときのノードをMCLNumとしている。
#発現のピークはこの三つ。3h-1h, 12h-3h, 24h-12hで比較したときに差の符号が変わった時がピーク
#例えば3h-1h=+, 12h-3h=-であれば3hがピーク
#使ったパッケージ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(stringr)
library(dplyr)
library(purrr)
#自動化するための準備----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
condition <- c("CY15", "CY16", "CY20")
MCLNum <- unique(TransCisNetwork$Cis)
MCLNum <- MCLNum[grep("MCLNum", MCLNum)]
time <- c("1h", "3h", "12h")
TransCis_pair <- paste0(TransCisNetwork$Trans, TransCisNetwork$Cis)　#最終的にCis-Transを統合させるから
e <- 1
#CY毎のfor入ります----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
for(e in e:length(condition)){
  #各時間の発現の差の計算。1h-0h, 3h-1h, 12h-3h, 24h-12h--------------------------------------------------------------------------------------------------------------------------------------
  diff <- c()　#発現の差を格納
  allExp <- c()　#各MCLNumの平均発現値を格納
  temp <- substr(MCLNum, start = 7, stop = 9)　#MCLNum数字のみを一時格納
  n <- 1
  #Cisクラスター毎のfor入ります------------------------------------------------------------------------------------------------------------------------------------------------------------------
  for(n in n:length(temp)){
    T_data <- allRNASeq[match(MasterTable %>% filter(MCLNum == temp[n]) %>% select(AGI) %>% unlist(use.names = F), rownames(allRNASeq)), ]
    T_Expression <- T_data %>% select(ends_with("h")) %>% select(starts_with(condition[e])) %>% select(-ends_with("48h")) %>% map(mean) %>% unlist()#mapは各列で同様な処理をするときに良く使う
    #クラスターの平均発現値を前後の時間で引いた値が正なのか負なのか数値化(正なら+1、負なら-1)
    allExp <- rbind(allExp, data.frame(Time01h = T_Expression[1],
                                       Time03h = T_Expression[2],
                                       Time12h = T_Expression[3],
                                       Time24h = T_Expression[4])
    )
    diff <- rbind(diff, data.frame(diff01h_00h = allExp$Time01h[n]-0,
                                   diff03h_01h = allExp$Time03h[n]-allExp$Time01h[n],
                                   diff12h_03h = allExp$Time12h[n]-allExp$Time03h[n],
                                   diff24h_12h = allExp$Time24h[n]-allExp$Time12h[n])
    )
    n <- n+1
  }
  rownames(allExp) <- MCLNum
  diff <- diff %>% mutate(diff01h_00h_sign = if_else(diff01h_00h > 0, true = "p", false = "n"), 
                          diff03h_01h_sign = if_else(diff03h_01h > 0, true = "p", false = "n"),
                          diff12h_03h_sign = if_else(diff12h_03h > 0, true = "p", false = "n"),
                          diff24h_12h_sign = if_else(diff24h_12h > 0, true = "p", false = "n"))
  rownames(diff) <- MCLNum
  #トランス側の発現値-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  T_data <- allRNASeq[match(unique(TransCisNetwork$Trans)[grep("AT", unique(TransCisNetwork$Trans))], rownames(allRNASeq)), ] %>% select(ends_with("h")) %>% select(starts_with(condition[e])) %>% select(-ends_with("48h"))
  Trans_Exp <- data.frame(AGI = rownames(T_data),
                          T_data %>% select(1)-0,
                          T_data %>% select(2)-T_data %>% select(1),
                          T_data %>% select(3)-T_data %>% select(2),
                          T_data %>% select(4)-T_data %>% select(3),
                          stringsAsFactors = F
                          )
  colnames(Trans_Exp) <- c("AGI", "diff01h_00h", "diff03h_01h", "diff12h_03h", "diff24h_12h")
  Trans_Exp <- Trans_Exp %>% mutate(diff01h_00h_sign = if_else(diff01h_00h > 0, true = "p", false = "n"), 
                                    diff03h_01h_sign = if_else(diff03h_01h > 0, true = "p", false = "n"),
                                    diff12h_03h_sign = if_else(diff12h_03h > 0, true = "p", false = "n"),
                                    diff24h_12h_sign = if_else(diff24h_12h > 0, true = "p", false = "n"))
  
  #CY15の全クラスターの平均発現値の差をヒートマップで可視化----------------------------------------------------------------------------------------------------------------------------------
  #全クラスターの各時間の平均発現値の差をggplotでヒートマップ化するためのデータフレーム
  T_data <- data.frame(diffexp = c(diff[, 1], diff[, 2], diff[, 3], diff[, 4]),
                       Sample = rep(paste0("MCLNum", temp), times = 4),
                       condition = rep(colnames(diff)[1:4], each = length(temp)),
                       stringsAsFactors = F
  )
  g1 <- ggplot(T_data, aes(x = condition, y = Sample, fill = diffexp))
  g1 <- g1 + geom_tile()
  g1 <- g1 + theme_bw()
  g1 <- g1 + scale_fill_gradient2(low = "blue", high = "red", na.value = "white")
  g1 <- g1 + ggtitle(paste0(condition[e], "_CisCluster"))
  plot(g1)
  #Cis(制御される)側-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  #発現の変化が大きいものを検出
  difftime <- paste0(c("01h_00h", "03h_01h", "12h_03h", "24h_12h"), "_sign")
  check_max <- c()
  check_min <- c()
  k <- 1
  for(k in k:nrow(diff)){
    temp <- diff[k, 5:ncol(diff)]
    temp <- temp == "p"
    if(sum(temp) == 0){#nのみ
      check_max <- c(check_max, NA)
      check_min <- c(check_min, difftime[which.min(diff[k, 1:4])])
    }
    if(sum(temp) != 0 & sum(temp) != 4){#p, n共にあるとき
      check_max <- c(check_max, difftime[which.max(diff[k, 1:4])])
      check_min <- c(check_min, difftime[which.min(diff[k, 1:4])])
    }
    if(sum(temp) == 4){#pのみ
      check_max <- c(check_max, difftime[which.max(diff[k, 1:4])])
      check_min <- c(check_min, NA)
    }
    k <- k+1
  }
  #このforではcheck_max, check_minを基に、diffにpp, nnを上書きするためのもの
  m <- 1
  for(m in m:nrow(diff)){
    if(sum(diff[m, 5:8] == "p") != 0){
      diff[m, grep(check_max[m], colnames(diff))] <- "pp"
    }
    if(sum(diff[m, 5:8] == "n") != 0){
      diff[m, grep(check_min[m], colnames(diff))] <- "nn"
    }
    m <- m+1
  }
  
  n <- 1
  for(n in n:nrow(diff)){
    if(sum(diff[n, 5:8] == "p") != 0){
      diff[n, grep(check_max[n], colnames(diff))] <- "pp"
    }
    if(sum(diff[n, 5:8] == "n") != 0){
      diff[n, grep(check_min[n], colnames(diff))] <- "nn"
    }
    n <- n+1
  }
  #Trans(制御する)側-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  #発現の変化が大きいものを検出
  difftime <- paste0(c("01h_00h", "03h_01h", "12h_03h", "24h_12h"), "_sign")
  check_max <- c()
  check_min <- c()
  k <- 1
  for(k in k:nrow(Trans_Exp)){
    temp <- Trans_Exp[k, 6:ncol(Trans_Exp)]
    temp <- temp == "p"
    if(sum(temp) == 0){#nのみ
      check_max <- c(check_max, NA)
      check_min <- c(check_min, difftime[which.min(Trans_Exp[k, 2:5])])
    }
    if(sum(temp) != 0 & sum(temp) != 4){#p, n共にあるとき
      check_max <- c(check_max, difftime[which.max(Trans_Exp[k, 2:5])])
      check_min <- c(check_min, difftime[which.min(Trans_Exp[k, 2:5])])
    }
    if(sum(temp) == 4){#pのみ
      check_max <- c(check_max, difftime[which.max(Trans_Exp[k, 2:5])])
      check_min <- c(check_min, NA)
    }
    k <- k+1
  }
  #このforではcheck_max, check_minを基に、Trans_Expにpp, nnを上書きするためのもの
  m <- 1
  for(m in m:nrow(Trans_Exp)){
    if(sum(Trans_Exp[m, 5:8] == "p") != 0){
      Trans_Exp[m, grep(check_max[m], colnames(Trans_Exp))] <- "pp"
    }
    if(sum(Trans_Exp[m, 5:8] == "n") != 0){
      Trans_Exp[m, grep(check_min[m], colnames(Trans_Exp))] <- "nn"
    }
    m <- m+1
  }
  
  n <- 1
  for(n in n:nrow(Trans_Exp)){
    if(sum(Trans_Exp[n, 5:8] == "p") != 0){
      Trans_Exp[n, grep(check_max[n], colnames(Trans_Exp))] <- "pp"
    }
    if(sum(Trans_Exp[n, 5:8] == "n") != 0){
      Trans_Exp[n, grep(check_min[n], colnames(Trans_Exp))] <- "nn"
    }
    n <- n+1
  }
  
  #発現データのみでの発現制御予測-----------------------------------------------------------------------------------------------------------------------------------------------------------
  #Source側-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  S_regulate <- c()
  n <- 1
  for(n in n:nrow(Trans_Exp)){
    m <- 6
    for(m in m:c(ncol(Trans_Exp)-1)){
      if(str_detect(Trans_Exp[n, m], "p") && str_detect(Trans_Exp[n, m+1], "n")){
        S_regulate <- c(S_regulate, paste0("yama_", str_sub(colnames(Trans_Exp)[m], start = 5, end = 7)))
        names(S_regulate)[length(S_regulate)] <- Trans_Exp$AGI[n]
      }
      if(str_detect(Trans_Exp[n, m], "n") && str_detect(Trans_Exp[n, m+1], "p")){
        S_regulate <- c(S_regulate, paste0("tani_", str_sub(colnames(Trans_Exp)[m], start = 5, end = 7)))
        names(S_regulate)[length(S_regulate)] <- Trans_Exp$AGI[n]
      }
      m <- m+1
    }
    n <- n+1
  }
  #自己フィードバック(シスクラスター -> TF)を作るために
  n <- 1
  for(n in n:nrow(diff)){
    m <- 5
    for(m in m:c(ncol(diff)-1)){
      if(str_detect(diff[n, m], "p") && str_detect(diff[n, m+1], "n")){
        S_regulate <- c(S_regulate, paste0("yama_", str_sub(colnames(diff)[m], start = 5, end = 7)))
        names(S_regulate)[length(S_regulate)] <- rownames(diff)[n]
      }
      if(str_detect(diff[n, m], "n") && str_detect(diff[n, m+1], "p")){
        S_regulate <- c(S_regulate, paste0("tani_", str_sub(colnames(diff)[m], start = 5, end = 7)))
        names(S_regulate)[length(S_regulate)] <- rownames(diff)[n]
      }
      m <- m+1
    }
    n <- n+1
  }
  
  #Target側-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  T_regulate <- c()
  n <- 1
  for(n in n:nrow(diff)){
    if(sum(str_detect(diff[n, ], "pp")) == 1){
      TT_regulate <- diff[n, setdiff(c(which(str_detect(diff[n, ], "pp"))-1):which(str_detect(diff[n, ], "pp")), c(1:4))]
      T_regulate <- rbind(T_regulate, data.frame(Target = rownames(diff[n, ]),
                                                 regulate = paste(TT_regulate, collapse = "_"),
                                                 Time = str_sub(colnames(diff)[which(str_detect(diff[n, ], "pp"))], start = 5, end = 11),
                                                 stringsAsFactors = F
      ))
    }
    if(sum(str_detect(diff[n, ], "nn")) == 1){
      TT_regulate <- diff[n, setdiff(c(which(str_detect(diff[n, ], "nn"))-1):which(str_detect(diff[n, ], "nn")), c(1:4))]
      T_regulate <- rbind(T_regulate, data.frame(Target = rownames(diff[n, ]),
                                                 regulate = paste(TT_regulate, collapse = "_"),
                                                 Time = str_sub(colnames(diff)[which(str_detect(diff[n, ], "nn"))], start = 5, end = 11),
                                                 stringsAsFactors = F
      ))
    }
    n <- n+1
  }
  #自己フィードバック(シスクラスター -> TF)を作るために
  n <- 1
  for(n in n:nrow(Trans_Exp)){
    if(sum(str_detect(Trans_Exp[n, ], "pp")) == 1){
      TT_regulate <- Trans_Exp[n, setdiff(c(which(str_detect(Trans_Exp[n, ], "pp"))-1):which(str_detect(Trans_Exp[n, ], "pp")), c(2:5))]
      T_regulate <- rbind(T_regulate, data.frame(Target = Trans_Exp$AGI[n],
                                                 regulate = paste(TT_regulate, collapse = "_"),
                                                 Time = str_sub(colnames(Trans_Exp)[which(str_detect(Trans_Exp[n, ], "pp"))], start = 5, end = 11),
                                                 stringsAsFactors = F
      ))
    }
    if(sum(str_detect(Trans_Exp[n, ], "nn")) == 1){
      TT_regulate <- Trans_Exp[n, setdiff(c(which(str_detect(Trans_Exp[n, ], "nn"))-1):which(str_detect(Trans_Exp[n, ], "nn")), c(2:5))]
      T_regulate <- rbind(T_regulate, data.frame(Target = Trans_Exp$AGI[n],
                                                 regulate = paste(TT_regulate, collapse = "_"),
                                                 Time = str_sub(colnames(Trans_Exp)[which(str_detect(Trans_Exp[n, ], "nn"))], start = 5, end = 11),
                                                 stringsAsFactors = F
      ))
    }
    n <- n+1
  }
  #発現データのみでエッジの作成-----------------------------------------------------------------------------------------------------------------------------------------------------------
  T_regulate <- T_regulate[nchar(T_regulate$regulate) > 2, ]
  temp <- c("03h_01h", "12h_03h", "24h_12h")
  ExpNetwork <- c()
  n <- 1
  for(n in n:length(time)){
    T_source <- names(S_regulate[grep(time[n], S_regulate)])
    T_sourceInfo <- str_split(S_regulate[grep(time[n], S_regulate)], pattern = "_", simplify = T)
    colnames(T_sourceInfo) <- c("katachi", "time")
    
    T_target <- T_regulate %>% filter(Time == temp[n] | Time == temp[n+1] | Time == temp[n+2])
    ExpNetwork <- rbind(ExpNetwork, data.frame(Source = rep(T_source, times = length(T_target$Target)),
                                               S_Time = T_sourceInfo[, 2],
                                               Target = rep(T_target$Target, each = length(T_source)),
                                               T_Time = rep(T_target$Time, each = length(T_source)),
                                               regulate = paste0(T_sourceInfo[, 1], "_", rep(T_target$regulate, each = length(T_source))),
                                               stringsAsFactors = F
    )
    )
    n <- n+1
  }
  #edgeの定義------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ExpNetwork$regulate[ExpNetwork$regulate == "tani_n_pp" | ExpNetwork$regulate == "yama_pp_nn" | ExpNetwork$regulate == "tani_p_pp" | ExpNetwork$regulate == "yama_n_nn" | ExpNetwork$regulate == "yama_p_nn" | ExpNetwork$regulate == "tani_nn_pp"] <- "possitive"
  ExpNetwork$regulate[ExpNetwork$regulate == "yama_n_pp" | ExpNetwork$regulate == "tani_pp_nn" | ExpNetwork$regulate == "yama_p_pp" | ExpNetwork$regulate == "tani_n_nn" | ExpNetwork$regulate == "tani_p_nn" | ExpNetwork$regulate == "yama_nn_pp"] <- "negative"
  #Trans-Cisと組み合わせ-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Exp_pair <- paste0(ExpNetwork$Source, ExpNetwork$Target)
  temp <- c()
  n <- 1
  for(n in n:length(TransCis_pair)){
    temp <- c(temp, which(TransCis_pair[n] == Exp_pair))
    n <- n+1
  }
  T.data <- ExpNetwork[temp, ]
  T.data <- T.data %>% mutate(attrEdge = paste0(T.data$S_Time, T.data$T_Time, T.data$regulate))
  T.data <- T.data[!duplicated(T.data), ]
  
  T.rownum <- grep("MCLNum", T.data$Source)
  temp <- T.data[T.rownum, ]
  T.MCLNum <- unique(temp$Source)
  TT.Time <- sort(unique(paste0(T.data$S_Time, T.data$T_Time)))
  test3 <- c()
  n <- 1
  for(n in n:length(T.MCLNum)){
    #Source:MCLNum
    test1 <- match(sort(unique(substr(temp$attrEdge[temp$Source == T.MCLNum[n]], start = 1, stop = 10))), TT.Time)
    #Target:MCLNumがいつか
    test2 <- min(match(sort(unique(substr(T.data$attrEdge[T.data$Target == T.MCLNum[n]], start = 1, stop = 10))), TT.Time))
    if(sum(test1 < test2) >= 1){
      for(m in which(test1 < test2)){
        test3 <- c(test3, rownames(temp[grep(TT.Time[which(test1 < test2)[m]], temp$attrEdge), ]))
        m <- m+1
      }
    }
    n <- n+1
  }
  T.data <- T.data[setdiff(rownames(T.data), setdiff(T.rownum, test3)), ]
  #attributeファイル--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  attrNode <- list(intersect(T.data$Source, T.data$Target),
                   setdiff(T.data$Source, T.data$Target),
                   setdiff(T.data$Target, T.data$Source)
  )
  names(attrNode) <- c("S&T,", "S,", "T,")
  attrNode <- data.frame(MCLNum = unlist(attrNode),
                         NodeColour = str_split(names(unlist(attrNode)), pattern = ",", simplify = T)[, 1],
                         stringsAsFactors = F)
  attrNode <- cbind(attrNode, MasterTable[match(attrNode$MCLNum, MasterTable$AGI), ] %>% select(TF, symbol, annotation))
  
  attrNode$symbol[is.na(attrNode$symbol)] <- str_sub(attrNode$MCLNum[is.na(attrNode$symbol)], start = 7, end = 9)
  rownames(attrNode) <- c()
  temp <- attrNode$MCLNum
  temp <- temp[grep("MCLNum", temp)]
  n <- 1
  for(n in n:length(temp)){
    attrNode$TF[attrNode$MCLNum == temp[n]] <- paste0("Cis:", paste(unique(names(unlist(allCisList)[temp[n] == unlist(allCisList)])), collapse = "|"))
    n <- n+1
  }
  
  T_time <- c("01h", "03h", "12h")
  n <- 1
  for(n in n:length(T_time)){
    temp <- T.data[T.data$S_Time == T_time[n], ]
    STNode <- unique(intersect(temp$Source, temp$Target))
    attrSTNode <- c()
    if(length(STNode) != 0){
      m <- 1
      for(m in m:length(STNode)){
        #Target側だけではなくSource側も意識しないといけない
        #S&Tの時間を結合させる
        attrSTNode <- c(attrSTNode, paste0("S&T", T_time[n], "_", paste(sort(unique(temp$T_Time[temp$Source == STNode[m] | temp$Target == STNode[m]])), collapse = "|")))
        names(attrSTNode)[m] <- STNode[m]
        m <- m+1
      }
    }
    SNode <- setdiff(temp$Source, STNode)
    TNode <- setdiff(temp$Target, STNode)
    attrTNode <- c()
    m <- 1
    for(m in m:length(TNode)){
      attrTNode <- c(attrTNode, paste0("T", paste(sort(unique(temp$T_Time[temp$Target == TNode[m]])), collapse = "|")))
      names(attrTNode)[m] <- TNode[m]
      m <- m+1
    }
    
    m <- 1
    attrSNode <- c()
    for(m in m:length(SNode)){
      attrSNode <- c(attrSNode, paste0("S", T_time[n], "*", paste(sort(unique(temp$T_Time[temp$Source == SNode[m]])), collapse = "|")))
      names(attrSNode)[m] <- SNode[m]
      m <- m+1
    }
    
    test2 <- data.frame(Node = c(SNode, TNode, STNode),
                        attrNode = c(attrSNode, attrTNode, attrSTNode)
    )
    rownames(test2) <- c()
    colnames(test2) <- c("Node", paste0(T_time[n], "_", "attrNode"))
    title <- paste0("bigdata/yasue/tGRN_Groping/inflation4/", condition[e], "/diffExpNetwork/Network/", condition[e], "_", T_time[n], "_attrNodecolor.txt")
    #write.table(test2, file = title, sep = "\t", quote = F, row.names = F)
    n <- n+1
  }
  #tGRN出力-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  assign(paste0(condition[e], "_ExpNetwork"), T.data)
  title <- paste0("~/bigdata/yasue/tGRN_Groping/inflation4/", condition[e], "/diffExpNetwork/Network/", paste(str_split(Sys.Date(), pattern = "-", simplify = T), collapse = ""), condition[e], "_tGRN_TFCisExp.txt")
  #write.table(T.data, title, sep = "\t", quote = F, row.names = F)
  #tGRN attribute出力-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  title <- paste0("~/bigdata/yasue/tGRN_Groping/inflation4/", condition[e], "/diffExpNetwork/Network/", condition[e], "_tGRN_TFCisExp_attrNode.txt")
  #write.table(attrNode, title, sep = "\t", quote = F, row.names = F)
  #tGRN発現データ-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  title <- paste0("~/bigdata/yasue/tGRN_Groping/inflation4/", condition[e], "/diffExpNetwork/", condition[e], "_tGRN_AvrallExp.txt")
  write.table(allExp, title, sep = "\t", quote = F, row.names = T)
  e <- e+1
}
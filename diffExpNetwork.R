#"~/Nakano_RNAseq/network_analysis/script/tGRN_analysis20181023/diffExpNetwork.R"
#目的：CY化合物毎に動的な遺伝子制御ネットワークを作成するための.txtファイルを出力させる
library(ggplot2)
condition <- c("CY15", "CY16", "CY20")
MCLNum <- attrTransCis$MCLNum #Cis-Trans Networkを作成したときのノードをMCLNumとしている。
time <- c("1h", "3h", "12h")#発現のピークはこの三つ。3h-1h, 12h-3h, 24h-12hで比較したときに差の符号が変わった時がピーク
#例えば3h-1h=+, 12h-3h=-であれば3hがピーク
TransCis_pair <- paste0(TransCisNetwork$Trans, TransCisNetwork$Cis)　#最終的にCis-Transを統合させるから
e <- 1
for(e in e:length(condition)){
  diff <- c()　#発現の差を格納
  T_MCLNum <- c()　#MCLNumがいらないからそれを除いた数字のみを格納
  allExp <- c()　#各MCLNumの平均発現値を格納
  n <- 1
  for(n in n:length(MCLNum)){
    temp <- str_split(MCLNum, pattern = "MCLNum")[[n]][2]　#MCLNum数字のみを一時格納
    T_AGI <- MasterTable$AGI[MasterTable$MCLNum == temp]　#クラスターのAGIを格納
    T_AGI <- T_AGI[!is.na(T_AGI)]　#クラスターの発現テーブルを格納
    T_data <- allRNASeq[match(T_AGI, rownames(allRNASeq)), ]
    #クラスターの平均発現値を時間ごとに格納
    T_Expression <- c(mean(T_data[, paste0(condition[e], "_1h")]),
                      mean(T_data[, paste0(condition[e], "_3h")]),
                      mean(T_data[, paste0(condition[e], "_12h")]),
                      mean(T_data[, paste0(condition[e], "_24h")])
    )
    #クラスターの平均発現値を前後の時間で引く
    diff01h_00h <- T_Expression[1]-0
    diff03h_01h <- T_Expression[2]-T_Expression[1]
    diff12h_03h <- T_Expression[3]-T_Expression[2]
    diff24h_12h <- T_Expression[4]-T_Expression[3]
    #クラスターの平均発現値を前後の時間で引いた値が正なのか負なのか数値化(正なら+1、負なら-1)
    diff01h_00h_sign <- sign(diff01h_00h)
    diff03h_01h_sign <- sign(diff03h_01h)
    diff12h_03h_sign <- sign(diff12h_03h)
    diff24h_12h_sign <- sign(diff24h_12h)
    
    allExp <- rbind(allExp, data.frame(Time01h = T_Expression[1],
                                       Time03h = T_Expression[2],
                                       Time12h = T_Expression[3],
                                       Time24h = T_Expression[4])
    )
    TTTime <- c(0, 1, 3, 12, 24)
    TTExp <- c(0, T_Expression[1], T_Expression[2], T_Expression[3], T_Expression[4])
    
    #クラスターの各時間の平均発現値の差と正負の判定を格納
    diff <- rbind(diff, cbind(diff01h_00h, diff03h_01h, diff12h_03h, diff24h_12h, diff01h_00h_sign, diff03h_01h_sign, diff12h_03h_sign, diff24h_12h_sign))
    T_MCLNum <- c(T_MCLNum, temp)
    print(n)
    n <- n+1
  }
  rownames(allExp) <- MCLNum
  
  #全クラスターの各時間の平均発現値の差をggplotでヒートマップ化するためのデータフレーム
  T_data <- data.frame(diffexp = c(diff[, 1], diff[, 2], diff[, 3], diff[, 4]),
                       Sample = rep(paste0("MCLNum", T_MCLNum), times = 4),
                       condition = rep(colnames(diff)[1:4], each = length(T_MCLNum)),
                       stringsAsFactors = F
  )
  
  library(ggplot2)
  g1 <- ggplot(T_data, aes(x = condition, y = Sample, fill = diffexp))
  g1 <- g1 + geom_tile()
  g1 <- g1 + theme_bw()
  g1 <- g1 + scale_fill_gradient2(low = "blue", high = "red", na.value = "white")
  g1 <- g1 + ggtitle(condition[e])
  plot(g1)
  
  #diffには全クラスターの各時間の平均発現値の差と正負の判定が格納されている
  #exp_patternはdiffを見やすくしたかったのとMCLNumを追加したかったために作成
  exp_pattern <- data.frame(Sample = paste0("MCLNum", T_MCLNum),
                            diff,
                            stringsAsFactors = F
  )
  
  #ここのforは各クラスターの平均発現値の差を正なら+1, 負なら-1と数値化しているのを利用して、符号の変化を知るために作成
  #最終的目標は、+1+(-1) or -1+(+1)が0なら符号が変わっているとしたい
  #peak_allには+1+(-1) or -1+(+1)=0, +1+1=2 or -1+(-1) = -2を格納
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
  
  #ここのforはpeak_allには+1+(-1) or -1+(+1)=0となっているものを検出することが目的
  #例えば(1h-0h) + (3h-1h) = 0なら1hが発現が変わったタイミング、ピークであると判定
  peak_time <- c()
  obname <- c()
  k <- 1
  for(k in k:length(peak_all)){
    temp2 <- paste0(names(peak_all)[k], "_", names(which(peak_all[[k]] == 0)))
    temp2 <- str_split(temp2, pattern = "_")
    w <- 1
    for(w in w:length(temp2)){
      peak_time <- c(peak_time, temp2[[w]][2])
      obname <- c(obname, temp2[[w]][1])
    }
    k <- k+1
  }
  names(peak_time) <- obname
  #peakが現れないサブクラスターも当然パターンとして考えられるため、そういったクラスターは""が格納されている。それを除く
  peak_time <- peak_time[peak_time != ""]
  
  #ここでは制御される(Cis側)のスクリプトで、発現が最も高いものと低いものを検出
  #最終的に発現の差が最も低い-のものはnn,発現の差が最も高い+のものは++と判定するための前処理
  #exp_patternは全クラスターの番号と各時間の平均発現値の差と正負の判定(+1, -1)が格納
  #+1ならp, -1ならnを入れる。なんでこんな二度手間をしているかはもう覚えてない。直すのが面倒なのでこのまま
  
  exp_pattern[, 6:9][exp_pattern[, 6:9] == 1] <- "p"
  exp_pattern[, 6:9][exp_pattern[, 6:9] == -1] <- "n"
  difftime <- paste0(c("0h_01h", "03h_01h", "12h_03h", "24h_12h"), "_sign")
  check_max <- c() #
  check_min <- c()
  k <- 1
  for(k in k:nrow(exp_pattern)){
    temp <- exp_pattern[k, 6:ncol(exp_pattern)]
    temp <- temp == "p"
    #常に正のものや、負のクラスターがある可能性を考慮して、
    #pが0なら発現の差がすべて-ということで、--しかない, check_minに--の時間を格納
    #pが4なら発現の差がすべてpということで、ppしかない, check_maxに++の時間を格納
    #0でも4でもないときはnもpもあるのでnn, pp共にある, check_minに--, check_maxに++の時間を格納
    if(sum(temp) == 0){#nのみ
      check_max <- c(check_max, NA)
      check_min <- c(check_min, difftime[which.min(exp_pattern[k, 2:5])])
    }
    if(sum(temp) != 0 & sum(temp) != 4){#p, n共にあるとき
      check_max <- c(check_max, difftime[which.max(exp_pattern[k, 2:5])])
      check_min <- c(check_min, difftime[which.min(exp_pattern[k, 2:5])])
    }
    if(sum(temp) == 4){#pのみ
      check_max <- c(check_max, difftime[which.max(exp_pattern[k, 2:5])])
      check_min <- c(check_min, NA)
    }
    k <- k+1
  }
  
  #このforではcheck_max, check_minを基に、exp_patternにpp, nnを上書きするためのもの
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
  #-1, +1をn, pにした理由を思い出した。--, ++を識別できるようにするためだった
  
  
  #このforは時系列発現データのみで制御関係を予測
  ExpTable <- c()
  o <- 1
  for(o in o:length(time)){
    target <- list()
    regulate <- list()
    T_time <- list()
    source <- list()
    source_time <- list()
    q <- o+6
    #制御する側のピークが1hのものと制御される側(diff03h_01h_sign, diff12h_03h_sign, diff24h_12h_sign)の発現が++ or --なら制御関係あり
    #制御する側のピークの検出が1hであれば、diff03h_01h+diff01h_00h == 0というものなので、diff01h_00h_signは制御の可能性がないと考えた
    for(q in q:ncol(exp_pattern)){
      target <- c(target, list(exp_pattern$Sample[which(exp_pattern[, q] == "pp" | exp_pattern[, q] == "nn")]))
      regulate <- c(regulate, list(exp_pattern[, q][which(exp_pattern[, q] == "pp" | exp_pattern[, q] == "nn")]))
      T_time <- c(T_time, list(rep(colnames(exp_pattern)[q], length(exp_pattern[, q][which(exp_pattern[, q] == "pp" | exp_pattern[, q] == "nn")]))))
      source <- c(source, list(names(peak_time[peak_time == time[o]])))
      source_time <- c(source_time, list(peak_time[peak_time == time[o]]))
      q <- q+1
    }
    
    allsource <- list() #allsourceには制御する数だけ同じ制御する側を増やさないといけない。Target側の数だけ増やす
    alltarget <- list() #alltargetには制御される数だけ制御される側を一様に増やさないといけない。Source側の数だけ増やす
    allregulate <- list() #allregulateには制御される制御関係だけ制御される関係を一様に増やさないといけない。Source側の数だけ増やす
    alltime <- list() #alltimeには制御される制御される側の時間を一様に増やさないといけない。Source側の数だけ増やす
    allsource_time <- list()　#allsource_timeには同じ制御する側の時間を増やさないといけない。Target側の数だけ増やす
    q <- 1
    for(q in q:length(target)){
      allsource <- c(allsource, list(rep(source[[q]], each = length(target[[q]]))))
      alltarget <- c(alltarget, list(rep(target[[q]], times = length(source[[q]]))))
      allregulate <- c(allregulate, list(rep(regulate[[q]], times = length(source[[q]]))))
      alltime <- c(alltime, list(rep(T_time[[q]], times = length(source[[q]]))))
      allsource_time <- c(allsource_time, list(rep(source_time[[q]], each = length(target[[q]]))))
    }
    
    #最終的にCis-Transを組み合わせるため、各クラスターがどっちなのかの情報
    attrTarget <- attrTransCis$attr[match(unlist(alltarget), attrTransCis$MCLNum)]
    attrsource <- attrTransCis$attr[match(unlist(allsource), attrTransCis$MCLNum)]
    
    #時系列発現データのみの制御関係をまとめた表
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
  
  #edgeの制御をこの段階ではpeakの平均発現の符号を考慮していない。
  #そのためたとえ制御する側が-で制御される側が--でも--となっているこれを++にしたい
  #制御する側が+で制御される側が++, 制御する側が-で制御される側が--は++
  #制御する側が+で制御される側が--, 制御する側が-で制御される側が++は--
  
  #sourceの平均発現値を獲得したい
  attredge <- rep(0, times = nrow(ExpTable))
  names(attredge) <- paste0(ExpTable$source, "_", ExpTable$sourcetime)
  o <- 1
  for(o in o:length(MCLNum)){
    k <- 1
    for(k in k:length(time)){
      temp <- exp_pattern[which(MCLNum[o] == exp_pattern$Sample), c(grep(paste0(time[k], "_sign"), colnames(exp_pattern)) - 1)]
      if(length(grep("p", temp)) > 0){
        attredge[names(attredge) == paste0(MCLNum[o], "_", time[k])] <- "p"
      }else{
        attredge[names(attredge) == paste0(MCLNum[o], "_", time[k])] <- "n"
      }
      k <- k+1
    }
    o <- o+1
  }
  #attredgeにはsourceの平均発現値が正ならp, 負ならnが入っている
  #ExpTable$regulateが制御される側のエッジ情報
  ExpTable$regulate[which(attredge == "p" & ExpTable$regulate == "pp")] <- "pp"
  ExpTable$regulate[which(attredge == "p" & ExpTable$regulate == "nn")] <- "nn"
  ExpTable$regulate[which(attredge == "n" & ExpTable$regulate == "pp")] <- "nn"
  ExpTable$regulate[which(attredge == "n" & ExpTable$regulate == "nn")] <- "pp"
  
  attr_ExpTable <- data.frame(ExpTable = c(unique(ExpTable$source), unique(ExpTable$target)),
                              attr = attrTransCis$attr[match(c(unique(ExpTable$source), unique(ExpTable$target)), attrTransCis$MCLNum)],
                              stringsAsFactors = F)
  
  assign(paste0(condition[e], "_ExpTable"), ExpTable)
  assign(paste0(condition[e], "_attrExpTable"), attr_ExpTable)
  
  ExpNetwork <- ExpTable[intersect(grep("Trans", ExpTable$attrsource), grep("Cis", ExpTable$attrtarget)), ]
  k <- 1
  temp <- c()
  for(k in k:length(TransCis_pair)){
    temp <- c(temp, grep(TransCis_pair[k], paste0(ExpNetwork$source, ExpNetwork$target)))
    k <- k+1
  }
  ExpNetwork <- ExpNetwork[temp, ]
  ExpNetwork <- ExpNetwork[!duplicated(ExpNetwork), ]
  ExpNetwork$regulate <- paste0(ExpNetwork$sourcetime, "_", ExpNetwork$targettime, "_", ExpNetwork$regulate)
  
  k <- 1
  for(k in k:length(time)){
    T_data <- ExpNetwork[ExpNetwork$sourcetime == time[k], ]
    STNode <- unique(intersect(T_data$source, T_data$target))
    attrSTNode <- list()
    allSTNode <- c()
    if(length(STNode) != 0){
      m <- 1
      for(m in m:length(STNode)){
        #Target側だけではなくSource側も意識しないといけない
        attrSTNode <- c(attrSTNode, list(paste0(time[k], "_", T_data$targettime[T_data$target == STNode[m]])))
        allSTNode <- c(allSTNode, rep(STNode[m], times = length(attrSTNode[[m]])))
        m <- m+1
      }
    }
    SNode <- unique(setdiff(T_data$source, STNode))
    TNode <- unique(setdiff(T_data$target, STNode))
    attrTNode <- c()
    m <- 1
    for(m in m:length(TNode)){
      temp <- unique(T_data$targettime[T_data$target == TNode[m]])
      attrTNode <- c(attrTNode, paste(temp, collapse = "|"))
      m <- m+1
    }
    
    T_data <- data.frame(Node = c(SNode, TNode, allSTNode),
                         attrNode = c(rep(time[k], times = length(SNode)), attrTNode, unlist(attrSTNode))
    )
    colnames(T_data) <- c("Node", paste0(time[k], "_", "attrNode"))
    title <- paste0("bigdata/yasue/tGRN_Groping/inflation4/", condition[e], "/diffExpNetwork/Network/", condition[e], "_", time[k], "_attrNodecolor.txt")
    write.table(T_data, file = title, sep = "\t", quote = F, row.names = F)
    k <- k+1
  }
  
  T_data <- data.frame(MCLNum = union(ExpNetwork$source, ExpNetwork$target),
                       Num = str_sub(union(ExpNetwork$source, ExpNetwork$target), start = 7, end = 9)
  )
  
  title <- paste0("bigdata/yasue/tGRN_Groping/inflation4/", condition[e], "/diffExpNetwork/Network/", condition[e], "_attrNodeLabel.txt")
  write.table(T_data, file = title, sep = "\t", quote = F, row.names = F)
  
  
  assign(paste0(condition[e], "_ExpNetwork"), ExpNetwork)
  title <- paste0("bigdata/yasue/tGRN_Groping/inflation4/", condition[e], "/diffExpNetwork/Network/", condition[e], "_ExpNetwork.txt")
  write.table(ExpNetwork, title, sep = "\t", quote = F, row.names = F)
  
  title <- paste0("bigdata/yasue/tGRN_Groping/inflation4/", condition[e], "/diffExpNetwork/supple/", condition[e], "_diffexp_table.txt")
  write.table(exp_pattern, title, sep = "\t", quote = F, row.names = F)
  title <- paste0("bigdata/yasue/tGRN_Groping/inflation4/", condition[e], "/diffExpNetwork/supple/", condition[e], "_ExpTable.txt")
  write.table(ExpTable, title, sep = "\t", quote = F, row.names = F)
  title <- paste0("bigdata/yasue/tGRN_Groping/inflation4/", condition[e], "/diffExpNetwork/supple/", condition[e], "_attrExpTable.txt")
  write.table(attr_ExpTable, title, sep = "\t", quote = F, row.names = F)
  
  e <- e+1
}
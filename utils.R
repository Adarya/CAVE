##########################
# Utility functions for preprocessing
# and analyzing low burden TP53 variants
##########################

# exon calls +- 5 splice site bases
exons_xslx <- function(df) {
  
  exon2 <-filter(df, Position>7579832 & Position<7579946)
  Exon2 <- cbind(exon2, data.frame("Exon"= "exon2"))
  
  exon3 <-filter(df, Position>7579693 & Position<7579727)
  Exon3 <- cbind(exon3, data.frame("Exon"= "exon3"))
  
  exon4 <- filter(df, Position>7579305 & Position<7579596)
  Exon4 <- cbind(exon4, data.frame("Exon"= "exon4"))
  
  exon5 <-filter(df, Position>7578364 & Position<7578560)
  Exon5 <- cbind(exon5, data.frame("Exon"= "exon5"))
  
  exon6 <-filter(df, Position>7578170 & Position<7578295)
  Exon6 <- cbind(exon6, data.frame("Exon"= "exon6"))
  
  exon7 <-filter(df, Position>7577492 & Position<7577614)
  Exon7 <- cbind(exon7, data.frame("Exon"= "exon7"))
  
  exon8 <-filter(df, Position>7577012 & Position<7577161)
  Exon8 <- cbind(exon8, data.frame("Exon"= "exon8"))
  
  exon9 <-filter(df, Position>7576846 & Position<7576932)
  Exon9 <- cbind(exon9, data.frame("Exon"= "exon9"))
  
  exon10 <-filter(df, Position>7573920 & Position<7574039)
  Exon10 <- cbind(exon10, data.frame("Exon"= "exon10"))
  
  exon11 <-filter(df, Position>7572921 & Position<7573014)
  Exon11 <- cbind(exon11, data.frame("Exon"= "exon11"))
  
  
  vaf_by_exons <- rbind(Exon2, Exon3, Exon4 ,Exon5 ,Exon6 ,Exon7 ,Exon8 ,Exon9, Exon10, Exon11)
  vaf_by_exons$VAF <- as.numeric(vaf_by_exons$Frequency) / 100
  vaf_by_exons$Exon <- factor(vaf_by_exons$Exon, 
                              levels = tolower(c("exon2", "Exon3", "Exon4" ,"Exon5" ,"Exon6",
                                                 "Exon7" ,"Exon8" ,"Exon9", "Exon10", "Exon11")))
  
  return(vaf_by_exons)
}

#####-------------#########
# normalize calls per exon length
#####-------------#########

exon_sizes <- as.data.frame(t(data.frame(
  exon2 = abs(7579832 - 7579946),
  exon3 = abs(7579693 - 7579727),
  exon4 = abs(7579305 - 7579596),
  exon5 = abs(7578364 - 7578560),
  exon6 = abs(7578170 - 7578295),
  exon7 = abs(7577492 - 7577614),
  exon8 = abs(7577012 - 7577161),
  exon9 = abs(7576846 - 7576932),
  exon10 = abs(7573920 - 7574039),
  exon11 = abs(7572921 - 7573014)
)))

exon_sizes$exon <- row.names(exon_sizes)
colnames(exon_sizes)[1] <- "size"
exon_sizes$size <- exon_sizes$size - 2

exon_size_norm_vaf <- function(exon_sizes_df, calls_df, threshold){
  VAF <- threshold
  
  above <- calls_df[calls_df$VAF >= VAF,]
  below <- calls_df[calls_df$VAF < VAF,]
  temp <- NULL
  
  for(i in 1:nrow(exon_sizes_df)){
    all <- nrow(calls_df[calls_df$Exon == exon_sizes_df$exon[i],])
    all <- all / exon_sizes_df$size[i]
    
    
    above_i <- nrow(above[above$Exon == exon_sizes_df$exon[i],])
    above_i <- above_i / exon_sizes_df$size[i]
    
    below_i <- nrow(below[below$Exon == exon_sizes_df$exon[i],])
    below_i <- below_i / exon_sizes_df$size[i]
    
    temp <- cbind(temp, c(all,above_i, below_i))
  }
  
  colnames(temp) <- exon_sizes_df$exon
  row.names(temp) <- c(paste0("All_",VAF), paste0("Above_",VAF), paste0("Below",VAF))
  
  return(temp)
}

#######################################
# GEIE
# function related to the GENIE analysis 
#######################################
pvals <- function(df, Genie, max_VAF, min_VAF, length_out, TCGA=NULL){

    VAF_seq <- seq(min_VAF, max_VAF, length.out=length_out)

    #subsetting and calculating according to a selected VAF
    d <- data.frame()
    t <- data.frame()
    pvals <- numeric()
    pvals_tcga <- rep(NA,length(VAF_seq))
    
    #go through each threshold
    for (i in 1:length(VAF_seq)){
      
      #high vaf
      sus <- df[df$VAF >= VAF_seq[i],]
      sus <- sus[!duplicated(sus$Position),]
      susdf <- data.frame("Group" = "Above threshold", "Position" = sus$Position)
      
      #low vaf
      art <- df[df$VAF < VAF_seq[i],]
      art <- art[!duplicated(art$Position),]
      art <- subset(art, !(art$Position %in% sus$Position))
      artdf <- data.frame("Group" = "Below threshold", "Position" = art$Position)
        
      # calculating frequency in GENIE
      freq_sus <- numeric()
      for(k in 1:NROW(susdf)){
          freq_sus[k] <- sum(as.numeric(susdf$Position[k])==as.numeric(Genie$Start_Position))
      }
      
      freq_art <- numeric()
      for(j in 1:NROW(artdf)){
          freq_art[j] <- sum(as.numeric(artdf$Position[j])==as.numeric(Genie$Start_Position))
      }
      
      wil <- wilcox.test(freq_sus, freq_art) #two-sided
      pvals[i] <- wil$p.value
      
      d <- rbind(d, cbind(nrow(susdf), sum(freq_sus), mean(freq_sus),
                          nrow(artdf), sum(freq_art), mean(freq_art)))
    
      #TCGA data format
      if(!is.null(TCGA)){
        freq_sus_TCGA <- numeric()
        for(k in 1:NROW(susdf)){
          temp <- TCGA[as.numeric(TCGA$Start.Pos) == as.numeric(susdf$Position[k]),]
          freq_sus_TCGA[k] <- nrow(temp)
        }
      
        freq_art_TCGA <- numeric()
        for(j in 1:NROW(artdf)){
          temp <- TCGA[as.numeric(TCGA$Start.Pos) == as.numeric(artdf$Position[j]),]
          freq_art_TCGA[j] <- nrow(temp)
        }

        wil_tcga <- wilcox.test(freq_sus_TCGA, freq_art_TCGA)$p.value
        t <- rbind(t, cbind(sum(freq_sus_TCGA), mean(freq_sus_TCGA), 
                            sum(freq_art_TCGA), mean(freq_art_TCGA),
                            wil_tcga))
      }
    }
    
    pdata <- data.frame("VAF_threshold" = VAF_seq, "P_value_genie" = pvals)

    colnames(d) <- c("n_pos_above", "sum_above", "mean_above", "n_pos_below", "sum_below","mean_below")
    
    if(!is.null(TCGA)){
      colnames(t) <- c("sum_above_TCGA", "mean_above_TCGA", "sum_below_TCGA", "mean_below_TCGA", "Pval_TCGA")
    }
    num_pvals <- as.numeric(as.character((pdata$P_value)))

    pdata <- cbind(pdata,d,num_pvals)
    if(!is.null(TCGA)){
    pdata <- cbind(pdata, t)
    }

    return(pdata)
}

################----------------------#
# check abundance of mutations in CSD #
###############-----------------------#
CSD_fun <- function(df, set6, n = 183, max_VAF, min_VAF, breaks){
  

  VAF_seq <- seq(min_VAF, max_VAF, breaks)
  
  #subseting and calculating the data according to a selected VAF
  data <- data.frame("VAF_threshold" = numeric(), 
                     "n_pos_csd"= numeric(), 
                     "n_pos_cohort"=numeric(), 
                     "Percent_above_threshold" = numeric(),
                     "Percent_below" = numeric(),
                     "Percent_positions_above" = numeric(), 
                     "Percent_of_position_of_below"=numeric(),
                     "above_below_CSD" = numeric(), 
                     "diff" = numeric(), 
                     "CSD_by_above"=numeric(), 
                     "CSD_by_below"=numeric())
  
  for (i in 1:length(VAF_seq)){
    sus <- df[df$VAF >= VAF_seq[i],]
    sus <- sus[!duplicated(sus$Position),]

    art <- df[df$VAF < VAF_seq[i],]
    art <- art[!duplicated(art$Position),]
    art <- subset(art, !(art$Position %in% sus$Position))
    
    freq_sus <- 0
    for(k in 1:NROW(sus)){
      if(as.numeric(sus$Position[k]) %in% set6$HG19_Start){
        freq_sus <- freq_sus +1
      }
    }
    freq_art <- 0
    for(j in 1:NROW(art)){
      if(as.numeric(art$Position[j]) %in% set6$HG19_Start){
        freq_art <- freq_art +1
      }
    }
    data[i,] <- c(VAF_seq[i], freq_sus, NROW(sus), (freq_sus/n)*100, 
                  (freq_art/n)*100, (NROW(sus)/(NROW(sus)+NROW(art)))*100, 
                  (NROW(art)/(NROW(sus)+NROW(art)))*100, freq_sus/freq_art, 
                  (freq_sus/n)*100-(NROW(sus)/(NROW(sus)+NROW(art)))*100, 
                  freq_sus/NROW(sus)*100, freq_art/NROW(art)*100)
    
  }
  return(data)
}

################----------------------#
# TP53_PROF model analysis #
###############-----------------------#
ml_model_fun <- function(unitedDF, max_VAF, min_VAF, breaks){
  
  VAF <- seq(min_VAF, max_VAF, breaks)
  
  ratio_df_above <- data.frame()
  ratio_df_below <- data.frame()
  
  for(i in 1:length(VAF)){
    above <- unitedDF[unitedDF$VAF >= VAF[i],]
    
    ratio_df_above <- rbind(ratio_df_above, 
                            c(VAF[i],nrow(above),
                              nrow(above[above$final_label == "D",]),
                              nrow(above[above$final_label == "ND",]),
                              nrow(above[above$final_label == "D",])/nrow(above[above$final_label == "ND",])))
    
    below <- unitedDF[unitedDF$VAF < VAF[i],]
    
    ratio_df_below <- rbind(ratio_df_below, 
                            c(VAF[i],nrow(below),
                              nrow(below[below$final_label == "D",]),
                              nrow(below[below$final_label == "ND",]),
                              nrow(below[below$final_label == "D",])/nrow(below[below$final_label == "ND",])))
  }
  colnames(ratio_df_above)<- c("VAF_threshold", "nrow","Pathogenic","Not_pathogenic","D/ND_Ratio")
  colnames(ratio_df_below)<- c("VAF_threshold", "nrow","Pathogenic","Not_pathogenic","D/ND_Ratio")
  return(list(ratio_df_above,ratio_df_below))
}


#####################
# dndscv
#######################
dnd_fun_above <- function(df, vaf_values){
  dnddf <- data.frame(VAF_treshold = vaf_values)
  dndsloc <- data.frame()
  gmuts <- data.frame()
  
  for(v in vaf_values){
    print(v)
    m <- dndscv(df[df$VAF >= v,], gene_list = c("TP53","PTEN"),
                max_muts_per_gene_per_sample = Inf,
                max_coding_muts_per_sample = Inf)
    dndsloc <- rbind(dndsloc, m$sel_loc[1,])
    gmuts <- rbind(gmuts, m$genemuts[2,])
  }
  
  return(cbind(dnddf, dndsloc, gmuts))
}

dnd_fun_below <- function(df, vaf_values){
  dnddf <- data.frame(VAF_treshold = vaf_values)
  dndsloc <- data.frame()
  gmuts <- data.frame()
  
  for(v in vaf_values){
    print(v)
    m <- dndscv(df[df$VAF < v,], gene_list = c("TP53","PTEN"),
                max_muts_per_gene_per_sample = Inf,
                max_coding_muts_per_sample = Inf)
    dndsloc <- rbind(dndsloc, m$sel_loc[1,])
    gmuts <- rbind(gmuts, m$genemuts[2,])
  }
  
  return(cbind(dnddf, dndsloc, gmuts))
}


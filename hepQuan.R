hepQuan <- function(scan,iso,ppm,dataset2,minscan){
  #peak smoothing, if scan numb res<1000, skip
  if (nrow(scan)>100){
    peaks <- smth.gaussian(scan$tic,window=100)
  }else{
    peaks<-scan$tic
  }
  #plot(peaks)
  peaks[which(is.na(peaks))] <- 0
  peaks <- as.numeric(peaks)
  
  #peaks <- peaks[101:(length(peaks)-100)]
  
  #step2 find peak centriod
  peak_find <- findpeaks(peaks,nups=3,ndowns=3,threshold = 2)
  peak_find <- as.data.frame(peak_find)
  
  
  #peak filter, delet peaks which baseline/peak <0.9
  peak_find <- peak_find[which(scan$tic[peak_find$V2]*0.95>scan$tic[peak_find$V3]),]
  
  
  #add peak group
  peak_find$group <- c(1:nrow(peak_find))
  
  
  
  #add peak group information of each scan numbers
  iso$peak.No <-0
  for (i in c(1:nrow(peak_find))){
    iso$peak.No[which(iso$scan_num>=peak_find$V3[i]&iso$scan_num<=peak_find$V4[i])] <- peak_find$group[i] 
  }
  
  #delete scan numbers which are not in peak group
  iso<- iso[which(iso$peak.No!=0),]
  
  
  #Step3 combine compositions with same MW
  
  #round monoMW
  iso$round_mw <- round(iso$monoisotopic_mw,2) 
  
  #calculate intensity
  #iso$intensity <- iso$abundance
  
  #add each MW counts
  iso$mw_count <-1
  
  #add bpi and time of each scan number
  index <- iso$scan_num
  #iso$bpi <- scan$bpi[index]
  iso$time <- scan$scan_time[index]
  
  #delete unnessary infromation
  iso <- iso[,c("scan_num","charge","mz","round_mw","mw_count","abundance","peak.No","monoisotopic_mw","time")]
  
  

  
  #combine data based mw, peak group, charge
  iso <- iso %>%
    group_by(round_mw,peak.No,charge) %>%
    summarise(
      mz = mean(mz),
      mono_mw = mean(monoisotopic_mw),
      scan_range = max(scan_num)-min(scan_num)+1,
      abundance = sum(abundance),
      #score = sum(log((abundance+1)/(max(scan_num)-min(scan_num+1)))),
      scan_count = sum(mw_count),
      time = mean(time)
    )
  iso$score <- round(log((iso$abundance/iso$scan_range),10)*log(iso$scan_range,10),1)
  
  #calculate score in this step
  
  
  #delete peaks ehich couunt<10
  #iso <- iso[which(iso$scan_count>=minscan),]

  
  #????????????
  #load data
  
  structure_searchingx <- function(mass,ppm,dataset){
    data <- dataset[which(abs(dataset$Theoretical.MW-mass)<0.5),c("Theoretical.MW","Structure","Adductive","dp",
                                                                  "unsatUA","GlcA","SO3")]
    #data <- data[which(round(data$Theoretical.MW,-1)==round(mass,-1)),]
    data$ppm <- round(((data[,1]-mass)/data[,1]*1000000),3)
    filter <- which((abs(data[,1]-mass)/data[,1]*1000000)<ppm)
    if (length(filter)>0) {
      data <- data[filter,]
    } else {
      data <- as.data.frame(matrix(0,1,8))
      colnames(data) <- c("Theoretical.MW","Structure","Adductive","dp","ppm","unsatUA","GlcA","SO3")
    }
    return(data)
  }
  
  
  res <-as.data.frame(matrix(0,1,17))
  colnames(res)<- c("Theoretical.MW","Structure","Adductive","dp","ppm","peak.No","charge","mz","mono_mw","scan_range","abundance","scan_count","score","time",
                    "unsatUA","GlcA","SO3")
  
  for(i in c(1:nrow(iso))){
    res_temp<-structure_searchingx(iso$mono_mw[i],ppm,dataset2)
    res_temp$peak.No <- iso$peak.No[i]
    res_temp$charge <- iso$charge[i]
    res_temp$mz <- iso$mz[i]
    res_temp$mono_mw <- iso$mono_mw[i]
    res_temp$abundance <- iso$abundance[i]
    res_temp$scan_range <- iso$scan_range[i]
    res_temp$scan_count <- iso$scan_count[i]
    res_temp$score <- iso$score[i]
    res_temp$time <- iso$time[i]
    res_temp <- res_temp[which(res_temp$Theoretical.MW!=0),]
    #delete high adductive
    res_temp <- res_temp[which(res_temp$Adductive<(res_temp$SO3+res_temp$unsatUA+res_temp$GlcA+res_temp$charge-2)),]
    #if(nrow(res_temp) >1){
      #res_t2 <- res_temp[which(t(as.data.frame(strsplit(res_temp$Structure,split=",")))[,4]!=3),]
      #stru_temp<-""
      #add_temp <-""
      #dp_temp <-""
      #for(i in c(1:nrow(res_temp))){
        #stru_temp <- paste0(res_temp$Structure[i],";",stru_temp)
        #add_temp <- paste0(res_temp$Adductive[i],";",add_temp)
        #dp_temp <- paste0(res_temp$dp[i],";",dp_temp)
      #}
      #res_temp <- res_temp[1,]
     # res_temp$Structure <- stru_temp
      #res_temp$Adductive <- add_temp
     # res_temp$dp <- dp_temp
      
   # }
    res <- rbind(res,res_temp)
  }
  rm(res_temp)
  
  
  #delete not matched peaks
  res <- res[which(res$Theoretical.MW!=0),]
  
  res <- res[,c("Theoretical.MW","Structure","Adductive","dp","ppm","peak.No","charge","mz","mono_mw","scan_range","abundance","scan_count","score","time",
                "unsatUA","GlcA","SO3","scan_count")]
  
  
  #res is Null, return res, if not, combine same structure
  if(nrow(res)!=0){
    #res$temp_count<- 1
    #combine same structure result
    #res <- aggregate(res[,c("Theoretical.MW","Degree.of.Polymerization","mono_mw","charge","Inten","peak.No","exep_ms",
                            #"time","mw.count","temp_count" )],by=list(res$Structure,res$charge),sum)
    
    #res$Theoretical.MW <- res$Theoretical.MW/res$temp_count
    #res$Degree.of.Polymerization <- res$Degree.of.Polymerization/res$temp_count
    #res$monoisotopic_mw <- res$monoisotopic_mw/res$temp_count
    #res$charge <-res$charge/res$temp_count
    #res$peak.No <-round(res$peak.No/res$temp_count,0)
    #res$time <- res$time/res$temp_count
    
    res <- res %>%
      group_by(Structure,charge,Adductive) %>%
      summarise(
        Theoretical.MW = mean(Theoretical.MW),
        #Adductive = mean(Adductive),
        dp = mean(dp),
        mz =round(mean(mz),4),
        mono_mw = mean(mono_mw),
        abundance = sum(abundance),
        score = round(sum(score),4),
        time = mean(time),
        scan_count=sum(scan_count),
        scan_range=sum(scan_range)
        
      )
    if(nrow(res)>30){
      res <- res[which(res$scan_range>=minscan),]
      res$log_ms <- log2(res$Theoretical.MW)
      model <-  lm(log_ms~time,data = res)
      pred.int <- predict(model, interval = "prediction")
      res <- cbind(res, pred.int)
      
      #res$fit_ms <- res$time*model$coefficients[2]+model$coefficients[1]
      # peaks in 0.77-1.3py will be kept
      #res$fit_msh <- log2(1.2*(2^res$fit_ms))
      #res$fit_msl <- log2(0.83*(2^res$fit_ms))
      res <-res[which(res$log_ms>res$lwr&res$log_ms<res$upr),]
    }
    #Delete un-trusted results. calculate linear regression, delete peaks out of  the linear regression

    
    #plot(res$time,log10(res$Theoretical.MW))
    
    
    #split structure and adductive information
    
    res$Exep.isotopic.mz <- round((res$mono_mw-1.00782*res$charge)/res$charge,4)

    
    #res <- res[which(res$scan_range<2*res$scan_count),]
  }else{
    res <- res[which(res$scan_range>=minscan),]
    res <- res[which(res$scan_count>=1),]
    #res <- res[which(res$scan_range<2*res$scan_count),]
  }
  
  
  
  
  return(res[,c("Structure","Adductive","Exep.isotopic.mz","charge","Theoretical.MW","abundance","score","scan_count","scan_range")])
  
}


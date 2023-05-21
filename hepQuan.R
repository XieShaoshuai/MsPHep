hepQuan <- function(scan,iso,ppm,dataset2,minscan, start, end){
  
  # Step1: TIC grouping
  scan_starting = scan$scan_num[1]-1
  scan$scan_num <- scan$scan_num-scan_starting
  iso$scan_num <- iso$scan_num-scan_starting
  
  
  #peak smoothing, if scan num <1000, skip TIC grouping step
  if (nrow(scan)>100){
    peaks <- smth.gaussian(scan$tic,window=40)
    peaks[which(is.na(peaks))] <- 0
    peaks <- as.numeric(peaks)
    #TIC grouping
    tic_group <- findpeaks(peaks,nups=3,ndowns=3,threshold = 2)
    tic_group <- as.data.frame(tic_group)
    
    #TIC grouping filter
    tic_group <- tic_group[which(scan$tic[tic_group$V2]*0.95>scan$tic[tic_group$V3]),]
  }else{
    tic_group <-as.data.frame(matrix(0,1,1))
    tic$V2 <- min(scan$scan_num)
    tic_group$V3 <- max(scan$scan_num)
  }

  #add peak group
  tic_group$group <- c(1:nrow(tic_group))

  #add peak group information of each scan numbers
  iso$peak.No <-0
  for (i in c(1:nrow(tic_group))){
    iso$peak.No[which(iso$scan_num>=tic_group$V3[i]&iso$scan_num<=tic_group$V4[i])] <- tic_group$group[i] 
  }
  
  #delete scan numbers which are not in peak group
  iso<- iso[which(iso$peak.No!=0),]
  
  #Step2 mass grouping
  
  #round monoMW
  iso$round_mw <- round(iso$monoisotopic_mw,2) 
  
  #add each MW counts
  iso$mw_count <-1
  
  #add bpi and time of each scan number
  index <- iso$scan_num
  iso$time <- scan$scan_time[index]
  
  #peak from start elution time to end elution time for analysis
  iso <- iso[which(iso$time>=start),]
  iso <- iso[which(iso$time<=end),]
  
  #delete unnessary infromation
  iso <- iso[,c("scan_num","charge","mz","round_mw","mw_count","abundance","peak.No","monoisotopic_mw","time")]
  
  iso_raw <- iso
  
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
  
  
  
  #delete peaks with scan number
  iso <- iso[which(iso$scan_range>=minscan),]
  


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
      colnames(data) <- c("Theoretical.MW","Structure","Adduct","dp","unsatUA","GlcA","SO3","ppm")
    } else {
      data <- as.data.frame(matrix(0,1,8))
      colnames(data) <- c("Theoretical.MW","Structure","Adduct","dp","unsatUA","GlcA","SO3","ppm")
    }
    return(data)
  }
  
  
  res <-as.data.frame(matrix(0,1,16))
  colnames(res)<- c("Theoretical.MW","Structure","Adduct","dp","ppm","peak.No","charge","mz","mono_mw","scan_range","abundance","scan_count","time",
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
    res_temp$time <- iso$time[i]
    res_temp <- res_temp[which(res_temp$Theoretical.MW!=0),]
    #delete high adductive
    res_temp <- res_temp[which(res_temp$Adduct<(res_temp$SO3+res_temp$unsatUA+res_temp$GlcA+res_temp$charge-2)),]
    res <- rbind(res,res_temp)
  }
  rm(res_temp)
  
  
  #delete not matched peaks
  res <- res[which(res$Theoretical.MW!=0),]
  
  res <- res[,c("Theoretical.MW","Structure","Adduct","dp","ppm","peak.No","charge","mz","mono_mw","scan_range","abundance","scan_count","time")]
  
  res$gaussian <- 0.5
  
  
  #calculate Gaussian similarity
  length <- length(res$Theoretical.MW)
  
  for (i in c(1:length)){
    data <- iso_raw %>% 
      filter(round_mw == round(res$mono_mw[i],2))
    
    if(nrow(data) > 2){
      td <- data$time
      d <- data$abundance
      mu <- data$time[data$abundance==max(data$abundance)]
      num_peak_pts <- length(data$abundance)
      
      sigma <- max(data$time)-min(data$time)
      h <- max(data$abundance)
      
      fit <- try(nls(d ~ SSgauss(td, mu, sigma, h)), silent = TRUE)
      
      if(class(fit) != "try-error"){
        gaussPts <- as.matrix(fitted(fit))
        gaussPts_std <- (gaussPts-mean(gaussPts))/sd(gaussPts)
        gaussPts_scale <- gaussPts_std/norm(gaussPts_std, type="F")
        
        d <- as.matrix(d)
        peak_intensity_std <- (d-mean(d))/sd(d)
        peak_intensity_scale <- peak_intensity_std/norm(peak_intensity_std, type="F")
        
        gauss_similarity <- sum(gaussPts_scale*peak_intensity_scale)
        #res$gaussian[i] <- gauss_similarity
        
      }else{
        gauss_similarity <- 0.5
      }
    }else{
      gauss_similarity <- 0.5
    }
    res$gaussian[i] <- gauss_similarity
  }
  
  #mingaussian <- min(res$gaussian[!is.na(res$gaussian)])
  #res$gaussian[is.na(res$gaussian)] <- mingaussian
  


  
  if(nrow(res)!=0){
    res <- res %>%
      group_by(Structure,charge,Adduct) %>%
      summarise(
        Theoretical.MW = mean(Theoretical.MW),
        #Adductive = mean(Adductive),
        dp = mean(dp),
        mz =round(mean(mz),4),
        mono_mw = mean(mono_mw),
        abundance = sum(abundance),
        time = mean(time),
        gaussian = max(gaussian),
        scan_count=sum(scan_count),
        scan_range=sum(scan_range)
        
      )
    if(nrow(res)>30){
      res <- res[which(res$scan_range>=minscan),]
      res$log_ms <- log2(res$Theoretical.MW)
      model <-  lm(log_ms~time,data = res)
      pred.int <- predict(model, interval = "prediction")
      res <- cbind(res, pred.int)
      
      res <-res[which((res$log_ms>res$lwr&res$log_ms<res$upr)|res$dp<4),]
    }
    
    
    res$Exep.isotopic.mz <- round((res$mono_mw-1.00782*res$charge)/res$charge,4)

    
    #res <- res[which(res$scan_range<2*res$scan_count),]
  }else{
    res <- res[which(res$scan_range>=minscan),]
    res <- res[which(res$scan_count>=1),]
    #res <- res[which(res$scan_range<2*res$scan_count),]
  }
  
  
  #calculate score
  res$score <- round(res$gaussian*log10(res$abundance),2)
  res$time <- round(res$time,2)
  
  return(res[,c("Structure","Adduct","Exep.isotopic.mz","charge","Theoretical.MW","abundance","score","time")])
  
}


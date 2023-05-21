hepQuan <- function(scan, iso, ppm, dataset, minscan, mintime,maxtime) {
  dataset <- read.delim("./database/Eno_hilic.txt", sep=";")
  df <- as.data.frame(matrix(unlist(dataset)))

  # slice the TIC to different groups for further data processing
  # if scan number is less 100, skip the slide step and treat tht TIC as a group
  if (nrow(scan) > 100) {
    # do peak smooting first to get proper group
    peaks <- smth.gaussian(scan$tic, window = 100)
    peaks[is.na(peaks)] <- 0
    peaks <- as.numeric(peaks)
    group <- findpeaks(peaks, nups = 3, ndowns = 3, threshold = 1)
    group <- as.data.frame(group)
    #group filter, delete groupss which baseline/group <0.95
    #group <- group[scan$tic[group$V2] * 0.95 > scan$tic[group$V3], ]
    group$group <- c(1 : nrow(group))
  } else {
    peaks <- scan$tic
  }

  #add bpi and time of each scan number
  index <- iso$scan_num
  iso$time <- scan$scan_time[index]
  rm(index)

  # add peak group information of each scan numbers
  iso$group.No <- 0
  for (i in c(1:nrow(group))) {
    iso$group.No[iso$scan_num >= group$V3[i] & iso$scan_num <= group$V4[i]] <- group$group[i] # nolint: line_length_linter.
  }

  # delete scan numbers which are not in peak group
  iso <- iso[which(iso$group.No != 0), ]

  #delete nosie peak that scan number less than threshold, default 3
  ## round monoMW
  iso$round_mw <- round(iso$monoisotopic_mw, 2)
  iso$round_mz <- round(iso$mz, 2)
  iso$count <- 1

  iso_summary <- iso  %>%
    group_by(round_mz, round_mw, group.No, charge) %>%
    summarise(
      mz = mean(mz),
      scan_num = max(scan_num) - min(scan_num) + 1,
      abundance = sum(abundance),
      scan_count = sum(count),
      time = median(time)
    )

  iso_summary <- iso_summary %>% 
    filter(scan_count >= 3 & time > mintime & time < maxtime)
  

  #search data  
  res <- as.data.frame(matrix(0, 1, 10))
  colnames(res) <- c(
    "Theoretical.MW", "Structure0", "Structure", "ppm",
    "group.No", "mz", "mono_mw", "abundance", "scan_count", "time")

  for (i in c(1:nrow(iso_summary))){
    mass <- iso_summary$round_mw[i]
    res_temp <- dataset[(abs(dataset$Theoretical.MW - mass) < 0.5), c(1:3)]
    res_temp$ppm <- abs(res_temp$Theoretical.MW-mass)/res_temp$Theoretical.MW * 1000000
    res_temp <- res_temp  %>%
      filter(ppm < 15)
    if (nrow(res_temp)>0){
      res_temp$group.No <- iso_summary$group.No[i]
      res_temp$mz <- iso_summary$mz[i]
      res_temp$mono_mw <- iso_summary$round_mw[i]
      res_temp$abundance <- iso_summary$abundance[i]
      res_temp$scan_count <- iso_summary$scan_count[i]
      res_temp$time <- iso_summary$time[i]
      res <- rbind(res, res_temp)
    }
  }
  
  # delete not matched peaks
  res <- res[which(res$Theoretical.MW != 0), ]

  #calculate regression line first, then group result
  
  if (nrow(res) > 30){
    res <- res %>% 
      filter(time > mintime & time < maxtime)
    res$log <- log2(res$Theoretical.MW)
    lm.out <- lm(log ~ time, data = res)
    newx = seq(min(res$time),max(res$time),by = 0.05)
    conf_interval <- predict(lm.out,interval="confidence",
                             level = 0.90)
    conf_interval <- predict(lm.out, interval = "prediction")
    plot(res$time, res$log, xlab="x", ylab="y", main="Regression")+
      abline(lm.out, col="lightblue")+
      lines(res$time, conf_interval[,2], col="red", lty=2)+
      lines(res$time, conf_interval[,3], col="blue", lty=2)
    res2 <- cbind(res, pred.int)
    res2 <- res2[which(res$log_ms > res$lwr & res$log_ms < res$upr), ]
  }
  

  if (nrow(res) != 0) {
    res2 <- res2 %>%
      group_by(Structure, Structure0, mz) %>%
      summarise(
        Theoretical.MW = mean(Theoretical.MW),
        mz = round(mean(mz), 4),
        mono_mw = mean(mono_mw),
        abundance = sum(abundance),
        time = mean(time),
        scan_count = sum(scan_count),
      )
  }



}

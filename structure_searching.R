structure_searching <- function(mz,charge,ppm,dataset){
  data <- dataset
  res<- as.data.frame(matrix(0,1,15))
  colnames(res) <-c("Theoretical.MW","Experimental.MW","isotopic.peak","Structure0","Structure","unsatUA",
                    "GlcA","GlcN" ,"Man", "Ac" ,
                    "SO3","anhydro",
                    "Adductive","dp","PPM")
  data$Experimental.MW <-round((mz*charge+1.00782*charge),5)
  data$PPM <- round(((data[,"Theoretical.MW"]-mz*charge-1.00782*charge)/data[,1]*1000000),3)
  filter <- which(abs(data$PPM)<ppm)
  if (length(filter)>0) {
    data <- data[filter,]
    #data <- data[,c(4,1,2,3,5)]
    data<-data[which(data$Experimental.MW!=0),]
  } else {
    data <-res
    data<-data[which(res$Experimental.MW!=0),]
  }
  
  return(data)
}


structure_searching2 <- function(mz,charge,ppm,dataset){
  data <- dataset
  res<- as.data.frame(matrix(0,1,15))
  colnames(res) <-c("Theoretical.MW","Experimental.MW","isotopic.peak","Structure0","Structure","unsatUA",
                    "GlcA","GlcN" ,"Man", "Ac" ,
                    "SO3","anhydro",
                    "Adductive","dp","PPM")
  for(i in c(1:5)){
    data1<-data
    data1$Experimental.MW <-round((mz*charge+1.00782*charge-(i-1)*1.00335),5)
    data1$PPM <- round(((data[,1]-data1$Experimental.MW)/data1[,1]*1000000),3)
    data1$isotopic.peak<-i
    filter <- which((abs((data1[,"Theoretical.MW"]-data1$Experimental.MW)/data1[,1]*1000000))<ppm)
    if (length(filter)>0) {
      data1 <- data1[filter,]
      #data1 <- data1[,c("Theoretical.MW","Experimental.MW","isotopic.peak","Structure","Adductive","dp","PPM")]
      #colnames(data1) <- colnames(res)
      res<- rbind(res,data1)
    } else {
      res<-res
    }
  }
  res<-res[which(res$Experimental.MW!=0),]
  return(res)
}

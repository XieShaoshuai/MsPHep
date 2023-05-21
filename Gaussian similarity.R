


data <- data[data$Time>43&data$Time<47,]

td <- data[,1]
d <- data$Intensity2
mu <- data$Time[data$Intensity2==max(data$Intensity2)]
num_peak_pts <- length(data[,3])

sigma <- max(data$Time)-min(data$Time)
h <- max(data$Intensity2)

fit <- try(nls(d ~ SSgauss(td, mu, sigma, h)), silent = TRUE)



gaussPts <- as.matrix(fitted(fit))
gaussPts_std <- (gaussPts-mean(gaussPts))/sd(gaussPts)
gaussPts_scale <- gaussPts_std/norm(gaussPts_std, type="F")

d <- as.matrix(d)
peak_intensity_std <- (d-mean(d))/sd(d)
peak_intensity_scale <- peak_intensity_std/norm(peak_intensity_std, type="F")

gauss_similarity <- sum(gaussPts_scale*peak_intensity_scale)



#plot
library(ggplot2)


pd <- as.data.frame(d)
colnames(pd) <- "real"
pd$time <- data$Time
pd$gaussian <- gaussPts

ggplot(pd) +
  geom_segment( aes(x=time, xend=time, y=0, yend=real))+
  geom_line(aes(x=time,y=gaussian),color="blue",size=1)+
  theme_bw()

















calculateGaussianSimilarity <- function(peakData, pts){
  peakrange <- peakData[c("rtmin", "rtmax")]
  ptsidx <- pts[, 1] >= peakrange[1] & pts[, 1] <= peakrange[2]
  intPts <- pts[ptsidx, ]
  
  num_peak_pts <- length(intPts[,2])
  
  if(length(intPts) > 2){
    td <- intPts[,1]
    d <- intPts[,2]
    mu <- peakData["rt"]
    sigma <- peakData["rtmax"] - peakData["rtmin"]
    h <- peakData["maxo"]
    
    fit <- try(nls(d ~ SSgauss(td, mu, sigma, h)), silent = TRUE)
    
    if(class(fit) != "try-error"){
      gaussPts <- as.matrix(fitted(fit))
      gaussPts_std <- (gaussPts-mean(gaussPts))/sd(gaussPts)
      gaussPts_scale <- gaussPts_std/norm(gaussPts_std, type="F")
      
      d <- as.matrix(d)
      peak_intensity_std <- (d-mean(d))/sd(d)
      peak_intensity_scale <- peak_intensity_std/norm(peak_intensity_std, type="F")
      
      gauss_similarity <- sum(gaussPts_scale*peak_intensity_scale)
      
    }else{
      gauss_similarity <- NA
    }
  }else{
    gauss_similarity <- NA
  }
  return(gauss_similarity)
}
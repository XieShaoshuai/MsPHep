library(IsoSpecR)
library(ggplot2)

iso_prob <- function(data,charge){
  res <- IsoSpecify( molecule = c(C= 6*data$unsatUA+6*data$GlcA+6*data$GlcN+6*data$Man+2*data$Ac,
                                  H=6*data$unsatUA+8*data$GlcA+11*data$GlcN+9*data$Man+2*data$Ac+2+3*data$Man+3*data$Adductive-2*data$anhydro,
                                  O=5*data$unsatUA+6*data$GlcA+4*data$GlcN+3*data$Man+data$Ac+3*data$SO3+1-data$anhydro,
                                  N=data$GlcN + data$Man+data$Adductive,
                                  S=data$SO3), 
                     stopCondition = .9999 )
  res <- as.data.frame(res)
  #if(nrow(res)>50){
   # res <- res[c(1:50),]
  #}
  
  res$mass_round <- round(res$mass,1)
  res <- res %>%
    group_by(mass_round) %>%
    summarise(
      mass = mean(mass),
      prob =sum(prob)
      
    )
  res$mz <- round((res$mass-1.00782*charge)/charge,3)
  #plot
  options(repr.plot.width = 2, repr.plot.height =3)
  ggplot(res[,], aes(x=mz, y=prob)) +
    geom_segment( aes(x=mz, xend=mz, y=0, yend=prob), color="red") +
    theme_light() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    xlab("m/z") +
    ylab("Abundance")
    
}


get_summary <- function(data) {
  if(nrow(data)>1){
    data %>%
      group_by(Structure) %>%
      summarise(
        #dp=dp,
        abundance = round(sum(abundance),0),
        score = round(sum(score),2)
      )
  }else{
    data <-data
  }

}


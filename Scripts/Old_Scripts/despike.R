############################################################
# spike detection function
############################################################
library(univOutl)
spike_flag <- function (NEE, block = 30, z = 5){
  flag <- rep(0,length(NEE))
  # 1. divide period
  num_b <- ceil(length(NEE)/48/block)
  for(i in 1:num_b){
    if(i==1){
      st <- 1
    } else{
      st <- last+1
    }
    if(i==num_b){
      last <- length(NEE)
    } else{
      last <- i*48*block
    }
    
    NEE_use <- NEE[st:last]
    
    if (length(which(is.na(NEE_use))) != length(NEE_use)){
      despike<- LocScaleB(NEE_use,k=z/0.6745,method = "MAD",return.dataframe=FALSE,exclude=NA)
      vec <- despike$outliers + st - 1
      flag[vec] <- 1
    }
    
    
  }
  return(flag)
}

## Diacamma data converter
## TrackingConverter.R
## N. Mizumoto

# Note
# --------------------------------------------------------------------------
# This script reads output from FastTrack tracking result (.txt) and
# convert it to next masking process, using python code (TrimFocusedTandem.ipynb)
# output will be in .csv
# --------------------------------------------------------------------------

# Packages -----------------------------------------------------------------
{
  library(data.table)
  library(stringr)
}

# Functions ----------------------------------------------------------------
## assuming they move straightly and constantly during non detection duration,
## I fill the non detected frame datasets
non.detect.fill <- function(x){
  non.detect.begin <- which(is.na(x))[c(T, diff(which(is.na(x))) > 1)]
  non.detect.end <- which(is.na(x))[c(diff(which(is.na(x))) > 1, T)]
  non.detect.num <- non.detect.end-non.detect.begin+1
  for(i in 1:length(non.detect.begin)){
    if(non.detect.begin[i]==1){
      x[non.detect.begin[i]:non.detect.end[i]] <- x[non.detect.end[i]+1]
    } else {
      x[non.detect.begin[i]:non.detect.end[i]] <- 
        (x[non.detect.end[i]+1] - x[non.detect.begin[i]-1]) / (non.detect.num[i]+1) * (1:non.detect.num[i]) + x[non.detect.begin[i]-1]
    }
  }
  return(x)
}

# Setup -------------------------------------------------------------------
{
  rm(list = ls())
  
  today <- Sys.Date()
  
  ## file location
  input.file.loc <- file.path("data/input")
  output.file.loc <- file.path("data/output")
  
}

# Computation --------------------------------------------------------------
{
  rawdata <- list.files(input.file.loc, full.names = TRUE, pattern = ".txt")
  data.name <- list.files(input.file.loc, full.names = FALSE, pattern = ".txt")
  
  for(j in 1:length(rawdata)){
    print(paste(j, "/", length(rawdata), data.name[j]))
    d <- read.table(rawdata[j], header=T)
    
    df <- d[,c("xBody", "yBody", "imageNumber", "id")]
    df <- df[df$id<1,]
    
    position <- (1:max(df$imageNumber)-1)
    x0 = y0 <- rep(NA, length(position))
    for(i in position){
      sdf <- df[df$imageNumber==i, ]
      if(sum(sdf$id == 0) > 0){
        x0[i+1] = sdf[sdf$id == 0, 1]
        y0[i+1] = sdf[sdf$id == 0, 2]
      }
    }
    
    x0 <- non.detect.fill(x0)
    y0 <- non.detect.fill(y0)
    
    df.output <- data.frame(position, x0, y0)
    
    df.output <- na.omit(df.output)
    write.table(df.output, file.path(output.file.loc, str_replace(data.name[j], "txt", "csv")),
                row.names = F, col.names=FALSE, sep=",")
  }
}

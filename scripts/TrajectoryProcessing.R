## Diacamma Tandem analysis
## TrajectoryProcessing.R
## N. Mizumoto

# Note
#---------------------------------------------------------------------------#
# This script processes the raw data for tandem trajectories of four ant/termite species.
# -Diacamma       sp.        (data from this study)
# -Temnothorax    rugatulus  (data from Valentini et al. 2020,    10.7554/eLife.55395)
# -Coptotermes    formosanus (data from Mizumoto and Dobata 2019, 10.1126/sciadv.aau6108)
# -Reticulitermes speratus   (data from Mizumoto and Dobata 2019, 10.1126/sciadv.aau6108)
# The raw data are at data/trajectory/raw/
# The processed data will be stored at data/trajectory/processed/
# The processed data will be used in TrajectoryOutputs.R for displaying results.
#---------------------------------------------------------------------------#

rm(list = ls())
preprocessAll()

#---------------------------------------------------------------------------#
{
  library(stringr)
  library(data.table)
  library(rinform) # transfer entropy
  
  ## constants
  fps    <- 29.97
  NMIN <- 15
  NSTEPS <- 45 # max interval for ss
  datasets = c("diacamma", "temnothorax", 
               "coptotermes", "reticulitermes")
  bodylength = c(10.54032, 2.34, 8.89, 5.5) # mm
  names(bodylength) = datasets
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
preprocessAll <- function(){
  ## data preparation
  #convert.rda.diacamma()
  
  ## transfer entropy analysis
  convert.rda.diacamma()
  normalizeTrajectoryLength()
  computeSpatialMeasureSamples()
  discretizeTrajectoriesLoop()
  computeGlobalITMeasures()
  findParamatersWithMaxNTE()
  computeLocalITMeasures()
  averageLocalMeasuresOverDistance()
  statsGlobalTransferEntropy()
  computeTandemStatistics()
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# This function reads csv files for Diacamma, and
# 1) process tracking error
# 2) output raw trajectories in .png file
# 3) convert them into R data.frame saved in .rda files
#---------------------------------------------------------------------------#
convert.rda.diacamma <- function(){
  ## data
  {
    rawdata <- list.files(("data/trajectory/raw/Diacamma-csv"), full.names = TRUE, pattern = ".csv")
    dataname <- list.files(("data/trajectory/raw/Diacamma-csv"), full.names = F, pattern = ".csv")
    dataname <- str_split(dataname, "_extract", simplify = T)[,1]
    dataname <- str_remove(dataname, "_eve")
    dataname <- paste0( str_sub(dataname, 1, 2), formatC(as.numeric(str_sub(dataname, 3, 4)), width=2, flag="0") )
    d.tracking.failed <- data.frame(fread("data/trajectory/Tracking_failed.csv", header=T))[,1:4]
  }
  
  ## scaling for Diacamma experiments
  {
    arena.size.pic = c(738.6533333, 755.24, 1108.373333, 1090.64)
    names(arena.size.pic) = c("KC", "KH", "KN", "NA")
    arena.size.mm = 1000
  }
  
  df <- data.frame()
  
  for(v in 1:length(rawdata)){
    
    # file info
    d <- data.frame(fread(rawdata[v], header=T))
    d <- d[,1:5]
    colnames(d) <- c("Frame", "LeaderX", "LeaderY", "FollowerX", "FollowerY")
    print(paste(v, "/", length(rawdata), "->", dataname[v]))
    
    Frametype = 0
    ID = paste(dataname[v],v,sep="-")
    Video = dataname[v]
    Colony <- str_sub(dataname[v], 1, 2)
    d <- data.frame(d, Colony, Video, ID, Frametype)
    
    # process tracking error
    # some videos have blackout (~2seconds), where I infer the individual location by assuming ants moving linearly and in constant speed
    # sometimes ants go under the arena, and cannot observe, where again I infer the individual location by assuming ants moving linearly and in constant speed
    # also remove pre-recording time
    d.ref <- d.tracking.failed[d.tracking.failed$id == Video, ]
    if(dim(d.ref)[1] > 0){
      for(i in 1:dim(d.ref)[1]){
        if(d.ref[i,4] == "pre"){
          d[d.ref[i,2]:d.ref[i,3], "Frametype"] <- -1
          d <- d[d$Frametype > -1,]
          d$Frame <- d$Frame - d$Frame[1]
        } else {
          d[d.ref[i,2]:d.ref[i,3], "Frametype"] <- 1
          pre <- d[d.ref[i,2]-1,2:5]
          post <- d[d.ref[i,3]+1,2:5]
          gapl <- (d.ref[i,3]-d.ref[i,2])+3
          for(j in 1:4){
            d[(d.ref[i,2]-1):(d.ref[i,3]+1), 1+j] <- seq(as.numeric(pre[j]), as.numeric(post[j]), length=gapl)
          }
        }
      }
    }
    
    ## scale
    d[, 2:5] <- d[, 2:5] * arena.size.mm / arena.size.pic[d[1,"Colony"]]
    
    df <- rbind(df,d)
  }
  f.name <- ("data/trajectory/raw/diacamma.rda")
  save(df, file = f.name)
  
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# This function reads the raw trajectories extracted from the videos 
# (in .rda format) of all datasets and shorten them to the same length of
# <NMIN> minutes.
#------------------------------------------------------------------------------#
normalizeTrajectoryLength <- function () {
  iodir            <- "data/trajectory/raw/"  
  nframes          <- round(fps * 60 * NMIN)
  
  for (ds in datasets) {
    fname <- paste(iodir, ds, ".rda", sep = "")
    load(fname)
    
    runs <- unique(df$ID)
    ndf  <- data.frame()
    cat("Dataset", ds, "-", length(runs), "tandem runs found\n")
    for (r in runs) {
      dftmp <- subset(df, ID == r)
      if (dim(dftmp)[1] < nframes) {
        cat("Warning: run", r, "is too short for given value of <NMIN> -> removed from datasets\n")
        next;
      }
      
      if (ds == "temnothorax") {
        dftmp  <- subset(dftmp, Frame < nframes)
      } else {
        maxlen <- dim(dftmp)[1]
        dftmp  <- dftmp[(maxlen - nframes + 1):maxlen,]
      }
      
      if (dim(dftmp)[1] != nframes) {
        cat("Warning: run", r, "is too short for given value of <NMIN> -> removed from datasets\n")
        next;
      }
      dftmp$Frame <- 0:(nframes - 1)
      ndf   <- rbind(ndf, dftmp)
    }
    
    df <- ndf
    fname <- paste(iodir, ds, NMIN, "min.rda", sep = "")
    save(df, file = fname)    
  }
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Compute step size, velocity while varying the sampling period 
# for a maximum number of samples for each run and dataset.
#------------------------------------------------------------------------------#
computeSpatialMeasureSamples <- function() {
  idir             <- "data/trajectory/raw/"
  odir             <- "data/trajectory/processed/"  
  maxStepsInterval <- NSTEPS
  maxNumSamples    <- 30000
  percentiles      <- seq(0.05, 0.35, by = 0.01)
  
  dfpercentiles <- data.frame(stringsAsFactors = !T)
  
  for (ds in datasets) {
    cat("Working on dataset: ", ds, "\n", sep = "")
    dfData  <- data.frame(stringsAsFactors = !T)
    fname   <- paste(idir, ds, NMIN, "min.rda", sep = "")
    load(fname)
    runs  <- unique(df$ID)
    
    for (r in runs) {
      dftmp  <- subset(df, ID == r)
      numPos <- min(dim(dftmp)[1], maxNumSamples)
      L      <- matrix(0, nrow = numPos, ncol = 2)
      L[, 1]  <- dftmp$LeaderX[1:numPos]
      L[, 2]  <- dftmp$LeaderY[1:numPos]
      F      <- matrix(0, nrow = numPos, ncol = 2)
      F[, 1]  <- dftmp$FollowerX[1:numPos]
      F[, 2]  <- dftmp$FollowerY[1:numPos]
      
      for (ss in 1:maxStepsInterval) {
        print(paste("runs:", r,"/", max(runs), "ss:",ss, "/", maxStepsInterval))
        tsteps <- seq(1, numPos, by = ss)
        l              <- L[tsteps, ]
        ld             <- sqrt(diff(l[, 1]) ^ 2 + diff(l[, 2]) ^ 2)
        nentries       <- length(ld)
        dftmp          <-
          data.frame(SS = rep(ss, nentries),
                     stringsAsFactors = !T)
        dftmp$ID       <- r
        dftmp$Role     <- "Leader"
        dftmp$Distance <- ld
        dftmp$Velocity <- ld / (ss / fps)
        dftmp$Dataset  <- ds
        dfData         <- rbind(dfData, dftmp)
        
        f              <- F[tsteps, ]
        fd             <- sqrt(diff(f[, 1]) ^ 2 + diff(f[, 2]) ^ 2)
        nentries       <- length(fd)
        dftmp          <-
          data.frame(SS = rep(ss, nentries),
                     stringsAsFactors = !T)
        dftmp$ID       <- r
        dftmp$Role     <- "Follower"
        dftmp$Distance <- fd
        dftmp$Velocity <- fd / (ss / fps)
        dftmp$Dataset  <- ds
        dfData         <- rbind(dfData, dftmp)
        
      } # End: ss
    } # End: runs
  
    names(dfData) <- c("SS", "ID", "Role", "Distance", "Velocity", "Dataset")
    fname <- paste(odir, "spatialMeasureSamples-", ds, NMIN, ".rda", sep = "")
    save(dfData, file = fname)
    
    # Compute second decile
    for (ss in 1:maxStepsInterval) {
      x   <- subset(dfData, SS == ss)$Distance
      qpercentiles <- quantile(x, prob = percentiles)
      dfpercentiles <- rbind(dfpercentiles,
                             c(ds, ss, qpercentiles,
                               length(x)),
                             stringsAsFactors = !T)
    }
  } # End: df
  
  fname <- paste(odir, "spatialMeasureSamples-decils", NMIN, ".txt", sep = "")
  names(dfpercentiles) <- c("Dataset", "SS", paste("D", percentiles * 100, sep = ""), "NS")
  write.table(dfpercentiles, file = fname, quote = !T, row.names = !T)
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# This function loops over different parameter configurations and,
# for each of them, calls the function discretizeTrajectories.
#------------------------------------------------------------------------------#
discretizeTrajectoriesLoop <- function () {
  idir      <- "data/trajectory/"
  stepsizes <- c(1:NSTEPS)
  
  fname  <- paste(idir, "processed/spatialMeasureSamples-decils", NMIN, ".txt", sep = "")
  decils <- read.table(fname, header = T)
  
  for(ds in datasets){
    fname <- paste(idir, "raw/", ds, NMIN, "min.rda", sep = "")
    for (ss in stepsizes) {
      tresh <- subset(decils, Dataset == ds & SS == ss)$D10
      discretizeTrajectories(fname, ss, tresh, NULL)
    }
  }
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# This function reads a dataset in the form of a .rda file containing a
# data.frame and adds two columns to it representing the discretized
# trajectories of the leader and of the follower.
#------------------------------------------------------------------------------#
discretizeTrajectories <- function (infname, stepsize = 1, threshold, outfname = NULL) {
  odir      <- "data/trajectory/processed/encoded/"
  trackX  <- c("LeaderX", "FollowerX")
  trackY  <- c("LeaderY", "FollowerY")
  trackM  <- c("LeaderM", "FollowerM")
  trackR  <- c("LeaderR", "FollowerR")
  trackMR <- c("LeaderMR", "FollowerMR")
  trackS  <- c(trackM, trackR, trackMR)
  
  cat("Working on dataset: ", infname, " SS = ", stepsize, "\n", sep = "")  
  
  encodeMotion <- function(x, y, threshold) {
    X     <- matrix(c(x, y), 2, 2)
    d     <- dist(X)
    if (d <= threshold) {
      # No motion state: 0
      entry <- 0
    } else {
      # Motion state: 1
      entry <- 1
    }
    entry
  }
  
  encodeRotation <- function(x, y, previous) {
    v1       <- c(x[2] - x[1], y[2] - y[1])
    v2       <- c(x[3] - x[2], y[3] - y[2])
    rotation <- v1[2] * v2[1] - v1[1] * v2[2]
    
    if (rotation < 0.0) {
      # Turn clockwise state: 0
      entry <- 0
    } else if (rotation > 0.0) {
      # Turn counter-clockwise state: 1
      entry <- 1
    } else if (rotation == 0) {
      # Collinear, no rotation: use previous entry if valid or
      # random rotation if not
      if (previous != -1) {
        entry <- previous
      } else {
        entry <- rbinom(1, 1, 0.5)
      }
    }
    entry
  }
  
  encodeMotionAndRotation <- function(x, y, threshold, previous) {
    X     <- matrix(c(x[2:3], y[2:3]), 2, 2)
    d     <- dist(X)
    if (d <= threshold) {
      # No motion state: 0
      entry <- 0
    } else {
      # Motion:
      v1       <- c(x[2] - x[1], y[2] - y[1])
      v2       <- c(x[3] - x[2], y[3] - y[2])
      rotation <- v1[2] * v2[1] - v1[1] * v2[2]
      
      if (rotation < 0.0) {
        # Turn clockwise state: 1
        entry <- 1
      } else if (rotation > 0.0) {
        # Turn counter-clockwise state: 2
        entry <- 2
      } else if (rotation == 0) {
        # Collinear, no rotation: use previous entry if valid or
        # random rotation if not
        if (previous > 0) {
          entry <- previous
        } else {
          entry <- 1 + rbinom(1, 1, 0.5)
        }
      }
    }
    entry
  }
  
  load(infname)  
  runs  <- unique(df$ID)
  dfnew <- data.frame()
  
  # Loop over each run..
  for (run in runs) {
    # Subset and subsample dataset..
    dfw           <- subset(df, df$ID == run)
    frames        <- seq(0, dim(dfw)[1], by = stepsize)
    dfw           <- dfw[dfw$Frame %in% frames, ]
    dfw[, trackS] <- -1
    nsteps        <- dim(dfw)[1]
    cat("..run: ", run, " (", nsteps, " tsteps)\n", sep = "")
    
    for (t in 2:(nsteps - 1)) {
      for (a in 1:2) {      
        x         <- dfw[c(t - 1, t, t + 1), trackX[a]]
        y         <- dfw[c(t - 1, t, t + 1), trackY[a]]
        
        dfw[t, trackM[a]]  <- encodeMotion(x[2:3], y[2:3], threshold)
        dfw[t, trackR[a]]  <- encodeRotation(x, y, dfw[t - 1, trackR[a]])
        dfw[t, trackMR[a]] <- encodeMotionAndRotation(x, y, threshold, dfw[t - 1, trackMR[a]])	
      }
    }
    dfnew <- rbind(dfnew, dfw)
  }
  dfnew <- subset(dfnew, LeaderM  != -1)
  dfnew <- subset(dfnew, LeaderR  != -1)
  dfnew <- subset(dfnew, LeaderMR != -1)
  
  # Save data file..
  if (is.null(outfname)) {
    outfname <- paste(odir, sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(infname)),
                      "-SS", stepsize, ".rda", sep = "")
  }
  save(dfnew, file = outfname)
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# This function computes global information theoretic measures for each tandem
# run in each dataset. Measures: entropy rate, transfer entropy.
#------------------------------------------------------------------------------#
computeGlobalITMeasures <- function(verbose = T) {
  idir      <- paste0("data/trajectory/processed/encoded/")
  odir      <- paste0("data/trajectory/processed/it-measures/")
  ncsamples <- 50        # Number of control samples..
  stepsizes <- c(1:NSTEPS)
  histories <- 1:17
  encodings <- c("M", "R", "MR")
  
  df <- data.frame(stringsAsFactors = !T)

  for(ds in datasets){
    for (ss in stepsizes) {
      # Load encoded data..
      fname  <- paste(idir, ds, NMIN, "min-SS", ss, ".rda", sep = "")
      load(fname)
      dfw    <- dfnew
      
      runs   <- unique(dfw$ID)
      nruns  <- length(runs)
      nsteps <- dim(dfw)[1] / nruns
      
      # Loop over encodings..
      dftmp1 <- data.frame(stringsAsFactors = !T)
      for (ee in encodings) {
        if (verbose) { cat("..working on dataset", ds, "stepsize", ss, "encoding", ee, "\n") }
        # Format time series..
        L <- matrix(-1, nrow = nsteps, ncol = nruns)
        F <- matrix(-1, nrow = nsteps, ncol = nruns)
        for (r in 1:nruns) {
          dfrun  <- subset(dfw, ID == runs[r]) 
          L[, r] <- dfrun[, paste("Leader", ee, sep = "")]
          F[, r] <- dfrun[, paste("Follower", ee, sep = "")]
        }
        
        # Loop over histories..
        for (k in histories) {
          # Compute IT metrics..
          ERl  <- colMeans(entropy_rate(L, k, local = T))
          ERf  <- colMeans(entropy_rate(F, k, local = T))
          TElf <- colMeans(transfer_entropy(L, F, k = k, local = T))
          TEfl <- colMeans(transfer_entropy(F, L, k = k, local = T))
          
          TElfc <- matrix(-1, nrow = ncsamples, ncol = nruns)
          TEflc <- matrix(-1, nrow = ncsamples, ncol = nruns)
          for (cs in 1:ncsamples) {
            x           <- sample(1:nruns)
            while(sum(x == 1:nruns) > 0) { x <- sample(1:nruns) }
            Fc          <- F[ , x]
            TElfc[cs, ] <- colMeans(transfer_entropy(L, Fc, k = k, local = T))
            TEflc[cs, ] <- colMeans(transfer_entropy(Fc, L, k = k, local = T))
          }
          
          # Store results..
          dftmp2          <- data.frame(Dataset = rep(ds, times = nruns), 
                                        Encoding = rep(ee, times = nruns), stringsAsFactors = !T)
          dftmp2$SS       <- ss
          dftmp2$K        <- k
          dftmp2$ID       <- runs
          dftmp2$ERl      <- ERl
          dftmp2$ERf      <- ERf
          dftmp2$TElf     <- TElf
          dftmp2$TEfl     <- TEfl
          dftmp2$TElfc    <- colMeans(TElfc)
          dftmp2$TEflc    <- colMeans(TEflc)
          dftmp1          <- rbind(dftmp1, dftmp2, stringsAsFactors = !T)
        }
      }
      df <- rbind(df, dftmp1, stringsAsFactors = !T)
    }
  }

  names(df)   <- c("Dataset", "Encoding", "SS", "K", "ID",
                   "ERl", "ERf", "TElf", "TEfl", "TElfc", "TEflc")
  df$Dataset <- as.factor(df$Dataset)
  df$Encoding <- as.factor(df$Encoding)
  df$SS       <- as.numeric(df$SS)
  df$K        <- as.numeric(df$K)
  df$ID       <- as.factor(df$ID)
  df$ERl      <- as.numeric(df$ERl)
  df$ERf      <- as.numeric(df$ERf)
  df$TElf     <- as.numeric(df$TElf)
  df$TEfl     <- as.numeric(df$TEfl)
  df$TElfc    <- as.numeric(df$TElfc)
  df$TEflc    <- as.numeric(df$TEflc)
  
  fname <- paste(odir, "globalITMeasures.rda", sep = "")
  save(df, file = fname)
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Find parameters (ss and k) that maximizes NTE for each species and encoding
#------------------------------------------------------------------------------#
findParamatersWithMaxNTE <- function() {
  iodir     <- "data/trajectory/processed/it-measures/"
  fname     <- paste0(iodir, "globalITMeasures.rda")
  load(fname) 
  
  Datasets  <- unique(df$Dataset)
  encodings <- unique(df$Encoding)
  SS        <- unique(df$SS)
  K         <- unique(df$K)
  
  dfp <- data.frame()
  for(ds in Datasets){
    for (ee in encodings) {
      # Aggregate the results..
      data <- data.frame()
      for (k in K) {
        for (ss in SS) {
          dfw   <- subset(df, Dataset == ds & Encoding == ee & K == k & SS == ss)
          TElf  <- mean(dfw$TElf)
          TEfl  <- mean(dfw$TEfl)
          TElfc <- TElf - mean(dfw$TElfc)
          TEflc <- TEfl - mean(dfw$TEflc)
          data  <- rbind(data, c(ds, k, ss, TElfc, TEflc))        
        }
      } 
      
     
      names(data) <- c("Datasets", "K", "SS", "TElfc", "TEflc")
      # IF TE after correction is negative, we set it to 0
      data[data$TElfc < 0 , "TElfc"] <- 0
      data[data$TEflc < 0 , "TEflc"] <- 0
      data$NTEc       <- abs(as.numeric(data$TElfc) - as.numeric(data$TEflc))
      data$direction  <- sign(as.numeric(data$TElfc) - as.numeric(data$TEflc))
      
      # Net transfer entropy
      row <- c(ds, as.character(ee), 
               as.numeric(data[data$NTEc == max(data$NTEc), c("K", "SS", "direction")]))
      dfp <- rbind(dfp, row, stringsAsFactors = !T)
    }
  }
  
  names(dfp)   <- c("Datasets", "Encoding", "K", "SS","Direction")
  dfp$Datasets <- as.factor(dfp$Datasets)
  dfp$Encoding <- factor(dfp$Encoding, levels = c("M", "R", "MR"))
  dfp$K        <- as.numeric(dfp$K)
  dfp$SS       <- as.numeric(dfp$SS)
  dfp$Direction       <- as.numeric(dfp$Direction)
  
  fname <- paste(iodir, "paramatersWithMaxNTE.rda", sep = "")
  save(dfp, file = fname)
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# This function computes local information theoretic measures for each tandem
# run in each dataset. Measures: entropy rate, transfer entropy.
#------------------------------------------------------------------------------#
computeLocalITMeasures <- function(verbose = T) {
  idir      <- "data/trajectory/processed/encoded/"
  odir      <- "data/trajectory/processed/it-measures/"
  encodings <- c("M", "R", "MR")
  
  fname <- paste(odir, "paramatersWithMaxNTE.rda", sep = "")
  load(fname) # dfp
  
  for (i in 1:dim(dfp)[1]) {
    ds <- dfp$Datasets[i]
    ee <- dfp$Encoding[i]
    k  <- dfp$K[i]
    ss <- dfp$SS[i]
    
    if (verbose) {
      cat("..working on dataset", ds, "stepsize", ss, "encoding", ee, "history", k, "\n")
    }
    
    # Load encoded data..
    fname  <- paste(idir, ds, NMIN, "min-SS", ss, ".rda", sep = "")
    load(fname)
    dfw    <- dfnew
    
    runs   <- unique(dfw$ID)
    nruns  <- length(runs)
    nsteps <- dim(dfw)[1] / nruns
    
    # Format time series..
    L <- matrix(-1, nrow = nsteps, ncol = nruns)
    F <- matrix(-1, nrow = nsteps, ncol = nruns)
    for (r in 1:nruns) {
      dfrun  <- subset(dfw, ID == runs[r]) 
      L[, r] <- dfrun[, paste("Leader", ee, sep = "")]
      F[, r] <- dfrun[, paste("Follower", ee, sep = "")]
    }
    
    # Compute IT metrics..
    ERl  <- entropy_rate(L, k, local = T)
    ERf  <- entropy_rate(F, k, local = T)
    TElf <- transfer_entropy(L, F, k = k, local = T)
    TEfl <- transfer_entropy(F, L, k = k, local = T)
    
    # Store results..
    dfw$ERl  <- NA
    dfw$ERf  <- NA
    dfw$TElf <- NA
    dfw$TEfl <- NA
    for (r in 1:nruns) {
      idxs           <- which(dfw$ID == runs[r])[(k + 1):nsteps]
      dfw$ERl[idxs]  <- ERl[, r]
      dfw$ERf[idxs]  <- ERf[, r]
      dfw$TElf[idxs] <- TElf[, r]
      dfw$TEfl[idxs] <- TEfl[, r]
    }
    
    df    <- dfw
    fname <- paste(odir, "localITMeasures-", ds, "-SS", ss, "K", k, ee, ".rda", sep = "")
    save(df, file = fname)
  }
}
#------------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# Obtain statistics of global Transfer Entropy,
# selected parameter configurations
#---------------------------------------------------------------------------#
statsGlobalTransferEntropy <- function() {
  dir <- "data/trajectory/processed/it-measures/"
  
  # Read data and parameters
  fname <- paste(dir, "globalITMeasures.rda", sep = "")
  load(fname)
  
  fname <- paste(dir, "paramatersWithMaxNTE.rda", sep = "")
  load(fname)
  
  dfp$Role <- factor(rep("L -> F", dim(dfp)[1]), levels = c("L -> F", "F -> L"))
  dfp[dfp$Direction == -1,]$Role <- "F -> L"
  
  
  # Compute statistics
  for (i in 1:dim(dfp)[1]) {
    dfw <- subset(df, Dataset == dfp$Datasets[i] &
                    Encoding == dfp$Encoding[i] &
                    K == dfp$K[i] & SS == dfp$SS[i])
    
    dfp$ERlmean[i]     <- mean(dfw$ERl)
    dfp$ERlsd[i]       <- sd(dfw$ERl)
    dfp$ERfmean[i]     <- mean(dfw$ERf)
    dfp$ERfsd[i]       <- sd(dfw$ERf)
    dfp$TElfmean[i]    <- mean(dfw$TElf)
    dfp$TElfsd[i]      <- sd(dfw$TElf)
    dfp$TEflmean[i]    <- mean(dfw$TEfl)
    dfp$TEflsd[i]      <- sd(dfw$TEfl)    
    dfp$TElfCmean[i]   <- mean(dfw$TElfc)    
    dfp$TElfCsd[i]     <- sd(dfw$TElfc)    
    dfp$TEflCmean[i]   <- mean(dfw$TEflc)
    dfp$TEflCsd[i]     <- sd(dfw$TEflc)
    dfp$NTEmean[i]     <- mean(abs(dfw$TElf - dfw$TEfl))
    dfp$NTEsd[i]       <- sd(abs(dfw$TElf - dfw$TEfl))
    dfp$NormNTEmean[i] <- mean(abs(dfw$TElf / dfw$ERf - dfw$TEfl / dfw$ERl))
    dfp$NormNTEsd[i]   <- sd(abs(dfw$TElf / dfw$ERf - dfw$TEfl / dfw$ERl))
    
    if (dfp$Role[i] == "F -> L") {
      dfp$NormTERolemean[i] <- mean((dfw$TEfl - dfw$TEflc) / dfw$ERl)
      dfp$NormTERolesd[i]   <- sd((dfw$TEfl - dfw$TEflc) / dfw$ERl)
    } else {
      dfp$NormTERolemean[i] <- mean((dfw$TElf - dfw$TElfc) / dfw$ERf)
      dfp$NormTERolesd[i] <- sd((dfw$TElf - dfw$TElfc) / dfw$ERf)
    }
    
    if (dfp$Role[i] == "F -> L") {
      dfp$TERolemean[i] <- mean(dfw$TEfl - dfw$TEflc)
      dfp$TERolesd[i]   <- sd(dfw$TEfl - dfw$TEflc)
    } else {
      dfp$TERolemean[i] <- mean(dfw$TElf - dfw$TElfc)
      dfp$TERolesd[i] <- sd(dfw$TElf - dfw$TElfc)
    }
    
    x <- wilcox.test(dfw$TElf, dfw$TElfc, alternative = "greater")
    dfp$TElfpvalueControl[i]  <- x$p.value
    dfp$TElfpvalueControlW[i] <- x$statistic
    
    x <- wilcox.test(dfw$TEfl, dfw$TEflc, alternative = "greater")
    dfp$TEflpvalueControl[i]  <- x$p.value
    dfp$TEflpvalueControlW[i] <- x$statistic
    
    x <- wilcox.test(dfw$TElf, dfw$TEfl, alternative = "greater", paired = T)
    dfp$TElfpvalue[i]  <- x$p.value
    dfp$TElfpvalueV[i] <- x$statistic
    
    x <- wilcox.test(dfw$TEfl, dfw$TElf, alternative = "greater", paired = T)
    dfp$TEflpvalue[i]  <- x$p.value
    dfp$TEflpvalueV[i] <- x$statistic
  }
  
  # Save statistics
  df <- dfp
  fname <- paste(dir, "globalITMeasuresStatistics.rda", sep = "")
  save(df, file = fname)
}
#---------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# This function computes average of local information theoretic measures,
# for binned inter-individual distance, for each tandem run in each dataset.
#------------------------------------------------------------------------------#
averageLocalMeasuresOverDistance <- function(verbose = T) {
  iodir   <- "data/trajectory/processed/it-measures/"
  
  # Load study parameters..
  fname <- paste(iodir, "paramatersWithMaxNTE.rda", sep = "")
  load(fname)
  
  # Loop over configurations..
  dfdata <- data.frame(stringsAsFactors = !T)
  for (pc in 1:dim(dfp)[1]) {
    ds <- as.character(dfp$Datasets[pc])
    ee <- as.character(dfp$Encoding[pc])
    k  <- as.numeric(dfp$K[pc])
    ss <- as.numeric(dfp$SS[pc])
    
    if (verbose) {
      cat("..working on dataset", ds, "stepsize", ss, "encoding", ee, "history", k, "\n")
    }
    
    fname <- paste(iodir, "localITMeasures-", ds, "-SS", ss, "K", k, ee, ".rda", sep = "")
    load(fname)
    
    df$Dist    <- sqrt((df$LeaderX - df$FollowerX)**2 + (df$LeaderY - df$FollowerY)**2)
    df$DistInc <- NA
    runs       <- unique(df$ID)
    for (r in runs) {
      idx <- which(df$ID == r)
      df$DistInc[idx[2:length(idx)]] <- df$Dist[idx[1:(length(idx) - 1)]] > df$Dist[idx[2:length(idx)]]
    }
    df <- na.omit(df)
    
    binsAll <- c(0, quantile(df$Dist, probs = seq(0.01, 1.0, by = 0.001)))
    binsInc <- c(0, quantile(subset(df, DistInc == T)$Dist, probs = seq(0.01, 1.0, by = 0.001)))
    binsDec <- c(0, quantile(subset(df, DistInc != T)$Dist, probs = seq(0.01, 1.0, by = 0.001)))
    for (b in 1:(length(binsAll) - 1)) {
      
      idx    <- which(df$Dist >= binsAll[b] & df$Dist < binsAll[b + 1])
      idxInc <- which(df$Dist >= binsInc[b] & df$Dist < binsInc[b + 1] & df$DistInc == T)
      idxDec <- which(df$Dist >= binsDec[b] & df$Dist < binsDec[b + 1] & df$DistInc != T)
      
      TElf <- mean(df$TElf[idx])
      row    <- c(ds, ee, k, ss, mean(df$Dist[idx]), "TElf", TElf)
      dfdata <- rbind(dfdata, row, stringsAsFactors = !T)
      
      TEfl <- mean(df$TEfl[idx])
      row    <- c(ds, ee, k, ss, mean(df$Dist[idx]), "TEfl", TEfl)
      dfdata <- rbind(dfdata, row, stringsAsFactors = !T)
      
      TElfInc <- mean(df$TElf[idxInc])
      row    <- c(ds, ee, k, ss, mean(df$Dist[idxInc]), "TElfInc", TElfInc)
      dfdata <- rbind(dfdata, row, stringsAsFactors = !T)
      
      TEflInc <- mean(df$TEfl[idxInc])
      row    <- c(ds, ee, k, ss, mean(df$Dist[idxInc]), "TEflInc", TEflInc)
      dfdata <- rbind(dfdata, row, stringsAsFactors = !T)
      
      TElfDec <- mean(df$TElf[idxDec])
      row    <- c(ds, ee, k, ss, mean(df$Dist[idxDec]), "TElfDec", TElfDec)
      dfdata <- rbind(dfdata, row, stringsAsFactors = !T)
      
      TEflDec <- mean(df$TEfl[idxDec])
      row    <- c(ds, ee, k, ss, mean(df$Dist[idxDec]), "TEflDec", TEflDec)
      dfdata <- rbind(dfdata, row, stringsAsFactors = !T)
    }
  }
  
  df <- dfdata
  names(df)   <- c("Datasets", "Encoding", "K", "SS", "Dist", "Stat", "Val")
  df$Datasets <- factor(df$Datasets)
  df$Encoding <- factor(df$Encoding)
  df$Stat     <- factor(df$Stat)
  for (i in c(3:5, 7)) { df[, i]  <- as.numeric(df[, i]) }
  
  fname <- paste(iodir, "localITMeasuresAveragedOverDistance.rda", sep = "")
  save(df, file = fname)
}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# This function reads coordinate datasets, and
# computes movement speed / accerelation / tandem run duration
# for each speceis
#------------------------------------------------------------------------------#
computeTandemStatistics <- function(){
  idir <- "data/trajectory/raw/"  
  odir <- "data/trajectory/processed/"
  
  df.all.speed = df.ind.speed = df.tandem.duration = data.frame(stringsAsFactors = !T)
  for(i in 1:4){
    
    ## load trajectory data
    ds = datasets[i]
    load(paste(idir, ds, 15, "min.rda", sep = ""))
    
    ## down sample for 5 FPS
    df <- df[df$Frame %% 6 == 0,]
    
    ## compute each characters
    {
      df$Time = df$Frame/30
      
      df$FollowerSpeed = c( NA,  sqrt( diff(df$FollowerX)^2 + diff(df$FollowerY)^2 ) )
      df$LeaderSpeed   = c( NA,  sqrt( diff(df$LeaderX)^2 + diff(df$LeaderY)^2 ) )
      df$FollowerAcc   = c( diff(df$FollowerSpeed), NA )
      df$LeaderAcc     = c( diff(df$LeaderSpeed), NA )
      df$Distance      = sqrt( (df$FollowerX -df$LeaderX)^2 + (df$FollowerY -df$LeaderY)^2 )
      df$FollowerSpeed[df$Frame==0] <- NA
      df$LeaderSpeed[df$Frame==0] <- NA
      df$FollowerAcc[df$Frame==0] <- NA
      df$LeaderAcc[df$Frame==0] <- NA
      df$FollowerAcc[df$Frame==max(df$Frame)] <- NA
      df$LeaderAcc[df$Frame==max(df$Frame)] <- NA
      
      df$Distance      = df$Distance / bodylength[i]
      df$FollowerSpeed = df$FollowerSpeed / bodylength[i]
      df$LeaderSpeed   = df$LeaderSpeed / bodylength[i]
      df$FollowerAcc   = df$FollowerAcc / bodylength[i]
      df$LeaderAcc     = df$LeaderAcc / bodylength[i]
      
      df$Tandem = df$Distance < 2
      
      angle_cal <- function(X, Y, Length){
        Ax <- (X[3:Length-1] - X[3:Length-2])
        Bx <- (X[3:Length] - X[3:Length-1])
        Ay <- (Y[3:Length-1] - Y[3:Length-2])
        By <- (Y[3:Length] - Y[3:Length-1])
        hugo <- (Ax * By - Ay * Bx + 0.000001)/abs(Ax * By - Ay * Bx + 0.000001)
        cos <- round((Ax * Bx + Ay * By) / ((Ax^2 + Ay^2)*(Bx^2 + By^2))^0.5,14)
        return(acos(cos)*hugo)
      }
      
      df$LeaderAngle = c(NA, angle_cal(df$LeaderX, df$LeaderY, dim(df)[1]), NA)
      df$LeaderAngle[df$Frame==0] <- NA
      df$LeaderAngle[df$Frame==max(df$Frame)] <- NA
      
      ## for Diacamma, we removed the frames with recording errors
      if(i==1){ # diacamma
        df <- df[df$Frametype == 0,]
      }
      
      df$Species = ds
    }
    
    ## summarize for each individual
    FollowerSpeed = tapply(df$FollowerSpeed, df[,c("ID", "Tandem")], mean, na.rm=T)
    LeaderSpeed   = tapply(df$LeaderSpeed,   df[,c("ID", "Tandem")], mean, na.rm=T)
    
    ## compute tandem duration
    {
      tandem = df$Tandem
      tan.timing = which(tandem)[c(T, diff(which(tandem))>1)]
      tan.end    = which(tandem)[c(diff(which(tandem))>1,T)]
      frame.zero = which(df$Frame == 0)
      frame.end  = which(df$Frame == max(df$Frame))
      frame.zero = frame.zero[tandem[frame.zero]]
      frame.end  = frame.end[tandem[frame.end]]
      
      tan.timing = unique(sort(c(tan.timing, frame.zero)))
      tan.end    = unique(sort(c(tan.end, frame.end)))
      
      tandem.duration = tan.end-tan.timing+1
      tandem.censor   = tan.end %in% frame.end
      
      ID = rep(0, length(tandem.duration))
      frame.zero = which(df$Frame == 0)
      for(f in 1:length(frame.end)){
        ID[tan.timing >= frame.zero[f]] = f
      }
    }
    
    ## output
    {
      df.all.speed <- rbind(df.all.speed, df[,c("Frame", "Colony", "ID", "FollowerSpeed", "LeaderSpeed", "FollowerAcc", "LeaderAcc", "Distance", "Species", "Tandem", "LeaderAngle")])
      
      
      dftemp <- rbind(
        data.frame(ID = row.names(FollowerSpeed),
                   Role = "Follower",
                   Type = "Tandem",
                   Speed = FollowerSpeed[,2],
                   Duration = tapply(df$Tandem, df$ID, sum)),
        data.frame(ID = row.names(FollowerSpeed),
                   Role = "Leader",
                   Type = "Tandem",
                   Speed = LeaderSpeed[,2],
                   Duration = tapply(df$Tandem, df$ID, sum)),
        data.frame(ID = row.names(FollowerSpeed),
                   Role = "Follower",
                   Type = "Separation",
                   Speed = FollowerSpeed[,1],
                   Duration = tapply(!df$Tandem, df$ID, sum)),
        data.frame(ID = row.names(FollowerSpeed),
                   Role = "Leader",
                   Type = "Separation",
                   Speed = LeaderSpeed[,1],
                   Duration = tapply(!df$Tandem, df$ID, sum))
      )
      dftemp$Species = ds
      df.ind.speed <- rbind(df.ind.speed, dftemp)
      
      df.tandem.duration = rbind(df.tandem.duration, 
                                 data.frame(Duration = tandem.duration,
                                            Censor = !tandem.censor,
                                            Species = ds,
                                            ID
                                 )
      )
    }
  }
  
  df.all.speed$Species <- factor(df.all.speed$Species, levels = datasets)
  
  df.ind.speed$Type    <- factor(df.ind.speed$Type, levels = c("Tandem", "Separation"))
  df.ind.speed$Species <- factor(df.ind.speed$Species, levels = datasets)
  
  df.tandem.duration$Species <- factor(df.tandem.duration$Species, levels = datasets)
  
  fname <- paste(odir, "tandemStatistics.rda", sep = "")
  save(df.all.speed, df.ind.speed, df.tandem.duration, file = fname)
}
#------------------------------------------------------------------------------#


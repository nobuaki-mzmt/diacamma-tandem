## Diacamma Tandem analysis
## TrajectoryOutputs.R
## N. Mizumoto

# Note
#---------------------------------------------------------------------------#
# This script displays the results for trajectory analysis of tandem runs.
# The processed data at data/trajectory/processed/ will be used.
# All the results will be stored at img/trajectory/
#---------------------------------------------------------------------------#

rm(list = ls())

#---------------------------------------------------------------------------#
{
  ## packages
  library(ggplot2)
  library(ggridges)
  library(viridis)
  
  library(lme4)
  library(car)
  library(multcomp)
  
  library(survival)
  library("survminer")
  require(coxme)
  
  ## constants
  {
    FPS = 30
    NMIN = 15
    
    pdfHeight  <- 7 
    pdfWidth   <- 10
    fontSize3C <- 28
    
    Datasets = datasets = c("diacamma", "temnothorax", "coptotermes", "reticulitermes")
    
    bodylength = c(10.54032, 2.34, 8.89, 5.5) # mm
    names(bodylength) = Datasets
  }
  
  output.all()
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
output.all <- function(){
  DisplayTrajectoryDiacamma()
  figNetTransferEntropyHeatmap()
  figNormalizedTransferEntropyBarplot()
  statglobalITMeasuresStatistics()
  figStepSizeDistribution()
  figLTEAveragedOverDistance()
  tandemOutput()
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# Display all the trajectories for Diacamma tandem runs
#---------------------------------------------------------------------------#
DisplayTrajectoryDiacamma <- function(){
  idir <- "data/trajectory/processed/"  
  odir <- "img/trajectory/"  
  f.name <- ("data/trajectory/raw/diacamma.rda")
  load(f.name)
  
  colonies = c("KC", "KH", "KN", "NA")
  
  for(cc in colonies){
    dftemp <- subset(df, Colony==cc)
    ggplot(dftemp) +
      geom_path(aes(x=LeaderX, y=LeaderY), col=viridis(2)[1], alpha=0.75) +
      geom_path(aes(x=FollowerX, y=FollowerY), col=viridis(2)[2], alpha=0.75) +
      coord_fixed() +
      ylab("y (mm)") +
      xlab("x (mm)") +
      ylim(c(-50,1050)) +
      xlim(c(450,1550)) +
      theme_bw()+
      theme(aspect.ratio = 1, legend.position = "none") +
      theme(strip.text = element_text(size = 8, margin = margin()))+
      theme(strip.background = element_rect(colour="#00000000", fill="#00000000"))+
      facet_wrap(.~Video, ncol=4)
    fname <- paste0(odir,"raw-trajectories/diacamma-",cc,"-alltrajectories.pdf")
    ggsave(fname, height = pdfHeight, width = pdfWidth)
    
    fname <- paste0(odir,"raw-trajectories/diacamma-",cc,"-alltrajectories.png")
    ggsave(fname, height = pdfHeight, width = pdfWidth)
  }
}  
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# Display heatmap of net transfer entropy over step length and history length
#---------------------------------------------------------------------------#
figNetTransferEntropyHeatmap <- function() {
  idir     <- "data/trajectory/processed/it-measures/"
  odir     <- "img/trajectory/"
  
  # Read data..
  fname     <- paste(idir, "globalITMeasures.rda", sep = "")
  load(fname) 
  datasets  <- unique(df$Dataset)
  encodings <- unique(df$Encoding)
  SS        <- unique(df$SS)
  K         <- unique(df$K) 
  
  theme_set(theme_minimal())
  for(ds in datasets){
    for (ee in encodings) {
      # Aggregate the results..
      data <- data.frame()
      for (k in K) {
        for (ss in SS) {
          dfw   <- subset(df, Dataset==ds & Encoding == ee & K == k & SS == ss)
          ERl   <- mean(dfw$ERl)
          ERf   <- mean(dfw$ERf)
          TElf  <- mean(dfw$TElf)
          TEfl  <- mean(dfw$TEfl)
          TElfc <- TElf - mean(dfw$TElfc)
          TEflc <- TEfl - mean(dfw$TEflc)
          data  <- rbind(data, c(k, ss, ERl, ERf, TElf, TEfl, TElfc, TEflc))        
        }
      } 
      names(data) <- c("K", "SS", "ERl", "ERf", "TElf", "TEfl", "TElfc", "TEflc")
      # IF TE after correction is negative, we set it to 0
      data[data$TElfc < 0 , "TElfc"] <- 0
      data[data$TEflc < 0 , "TEflc"] <- 0
      tmp <- abs(data$TElfc - data$TEflc)
      datamax  <- data[tmp == max(tmp),]
      
      # Net transfer entropy
      ggplot(data = data, aes(x = K, y = SS / FPS))+
        geom_raster(aes(fill = TElfc - TEflc), interpolate = !TRUE) +
        geom_point(data = datamax, aes(x = K, y = SS / FPS), shape = 18, size = 5.5, color="black") +
        labs(x = "History length (time steps)", y = "Sampling period (s)") +
        scale_fill_gradient2(midpoint = 0, mid="white", high="#cf2623", low="#1a6091") +
        theme(legend.position  = "bottom") +
        theme(legend.key.width = unit(5.5, "line")) +
        theme(text = element_text(size = fontSize3C)) +
        theme(plot.title = element_text(face = "bold", size = fontSize3C)) +
        theme(legend.title = element_blank())
      
      pdfName <- paste(odir, "/netTransferEntropyHeatmap-", ds, "-", ee, ".pdf", sep = "")  
      ggsave(pdfName, height = pdfHeight, width = pdfHeight)
    }
  }
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# Normalized Transfer Entropy Barplot
#---------------------------------------------------------------------------#
figNormalizedTransferEntropyBarplot <- function() {
  idir     <- "data/trajectory/processed/it-measures/"
  odir     <- "img/trajectory/"
  
  # Read data
  fname <- paste(idir, "globalITMeasuresStatistics.rda", sep = "")
  load(fname)
  
  idxFL <- which(df$Role == "F -> L")
  df$NormTERolemean[idxFL] <- -1 * df$NormTERolemean[idxFL]
  df$NormTERolesd[idxFL]   <- -1 * df$NormTERolesd[idxFL]
  df$TERolemean[idxFL] <- -1 * df$TERolemean[idxFL]
  df$TERolesd[idxFL]   <- -1 * df$TERolesd[idxFL]
  
  df$DatasetR <- paste(df$Role)
  df$DatasetR <- factor(df$DatasetR)
  
  df$Datasets <- factor(df$Datasets, levels=c("diacamma", "temnothorax", "coptotermes", "reticulitermes"))
  
  dftemp <- subset(df, Encoding!="MR")
  theme_set(theme_minimal())
  ggplot(data = dftemp, aes(x = Datasets, y = NormTERolemean, 
                            fill = Encoding, color = Encoding)) +
    geom_hline(yintercept = 0, color = "gray50") +
    geom_bar(stat = "identity", position = position_dodge2()) +
    geom_errorbar(stat = "identity", aes(ymin = NormTERolemean - NormTERolesd,
                                         ymax = NormTERolemean + NormTERolesd),
                  width = 0.2, position=position_dodge(0.9), color = "white") +         
    geom_errorbar(stat = "identity", aes(ymin = NormTERolemean,
                                         ymax = NormTERolemean + NormTERolesd),
                  width = 0.2, position=position_dodge(0.9)) +     
    scale_fill_viridis(discrete = T, end=0.5, direction = -1) +
    scale_color_viridis(discrete = T, end=0.5, direction = -1) +
    coord_cartesian(ylim = c(-0.15, 0.37)) +
    labs(x = "", y = "Normalized transfer entropy") + 
    theme(legend.position  = c(0.25 , 0.85)) +
    theme(axis.text.y = element_text(angle = 90)) +    
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size = fontSize3C)) +
    theme(text = element_text(size = fontSize3C))
  
  pdfName <- paste(odir, "/normalizedTransferEntropyBarplot.pdf", sep = "")  
  ggsave(pdfName, height = pdfHeight, width = pdfWidth * 1.15)
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# Statistical analysis comparing transfer entropy
#---------------------------------------------------------------------------#
statglobalITMeasuresStatistics <- function(){
  idir     <- "data/trajectory/processed/it-measures/"
  odir     <- "img/trajectory/"
  
  # Read data and parameters
  fname <- paste(idir, "globalITMeasuresStatistics.rda", sep = "")
  load(fname)
  dfp <- df
  
  # Output results as text file
  fname <- paste(odir, "globalITMeasuresStatistics.txt", sep = "")
  sink(fname)
  
  cat("----------------------------------------------\n")
  cat("Tests against surrogate data\n")
  cat("----------------------------------------------\n")
  for (i in 1:dim(dfp)[1]) {
    cat(as.character(dfp$Datasets[i]), as.character(dfp$Encoding[i]),"K",dfp$K[i],"SS",dfp$SS[i], "\n")
    cat("    L->F:", dfp$TElfpvalueControl[i], "W:", dfp$TElfpvalueControlW[i], "\n")
    cat("    F->L:", dfp$TEflpvalueControl[i], "W:", dfp$TEflpvalueControlW[i], "\n")
    if (i %% 3 == 0) cat("----------------------------------------------\n")
  }
  cat("\n")
  
  
  cat("----------------------------------------------\n")
  cat("Tests between Leader and Follower\n")
  cat("----------------------------------------------\n")
  for (i in 1:dim(dfp)[1]) {
    cat(as.character(dfp$Datasets[i]), as.character(dfp$Encoding[i]),"K",dfp$K[i],"SS",dfp$SS[i], "\n")
    cat("    L->F:", dfp$TElfpvalue[i], "V:", dfp$TElfpvalueV[i], "\n")
    cat("    F->L:", dfp$TEflpvalue[i], "V:", dfp$TEflpvalueV[i], "\n")
    if (i %% 3 == 0) cat("----------------------------------------------\n")
  }
  cat("\n")
  sink()
}  
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# Step size distributions and move/pause threshold
#---------------------------------------------------------------------------#
figStepSizeDistribution <- function() {
  idir     <- "data/trajectory/processed/"
  odir     <- "img/trajectory/"
  stepsInterval    <- seq(1, 45, by = 4)
  ssbounds         <- list(c(0, 9), c(0, 60), c(0, 40))
  
  theme_set(theme_minimal())
  for (ds in Datasets) {
    fname <- paste(idir, "spatialMeasureSamples-", ds, NMIN, ".rda", sep = "")
    load(fname)
    dfData <- subset(dfData, SS %in% stepsInterval)
    #dfData <- subset(dfData, SS %in% stepsInterval & Distance <= ssbounds[[dsidx]])
    
    # Plot the results..
    ggplot(dfData, aes(x = Distance, y = SS / FPS, group = SS, fill = factor(..quantile..))) +
      stat_density_ridges(geom = "density_ridges_gradient",
                          quantiles = c(0.1, 1.0), quantile_lines = TRUE,
                          panel_scaling = !T, rel_min_height = 0.001,
                          scale = 10, alpha = .7, color = "white", vline_color = "white")     +
      xlab("Step size (mm)") +
      ylab("Sampling period (s)") +
      scale_fill_viridis(discrete = T, alpha = .7, direction = 1) +
      ggtitle(ds) +
      theme(legend.position  = "none") +
      theme(text = element_text(size = fontSize3C)) +
      theme(plot.title = element_text(face = "bold", size = fontSize3C))
    
    pdfName <- paste(odir, "stepSizeDistribution-", ds, ".pdf", sep = "")
    ggsave(pdfName, height = pdfHeight, width = pdfWidth)
    
    ggplot(dfData, aes(x = Distance, y = SS / FPS, group = interaction(Role,SS), fill=Role)) +
      stat_density_ridges(geom = "density_ridges_gradient",
                          panel_scaling = !T, rel_min_height = 0.001,
                          scale = 10, alpha = .4, color = "white", vline_color = "white")     +
      xlab("Step size (mm)") +
      ylab("Sampling period (s)") +
      scale_fill_viridis(discrete = T, alpha = .4, direction = 1) +
      ggtitle(ds) +
      theme(legend.position  = "none") +
      theme(text = element_text(size = fontSize3C)) +
      theme(plot.title = element_text(face = "bold", size = fontSize3C))
    
    pdfName <- paste(odir, "stepSizeDistribution-", ds, "-Role.pdf", sep = "")
    ggsave(pdfName, height = pdfHeight, width = pdfWidth)
  }
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# Local Transfer entropy averaged over distance (show all direction encoding)
#---------------------------------------------------------------------------#
figLTEAveragedOverDistance <- function() {
  idir     <- "data/trajectory/processed/"
  odir     <- "img/trajectory/"

  # Load data.. 
  fname <- paste(idir, "it-measures/localITMeasuresAveragedOverDistance.rda", sep = "")
  load(fname)
  df      <- subset(df, Encoding != "MR")
  df$Stat <- as.factor(paste(df$Encoding, df$Stat, sep = ""))
  
  # Convert Dist from mm to Body length..
  for (ds in Datasets) {
    df[df$Datasets == ds,]$Dist <- df[df$Datasets == ds,]$Dist / bodylength[ds]
  }
  
  for(ds in Datasets) {
    dsIdx <- which(Datasets == ds)
    dfw   <- subset(df, Datasets == ds)
    
    xxlim    <- c(0.5, 2.5)
    yylim    <- c(-0.1, 0.4)
    llpos    <- c(0.35, 0.9)
    
    dfw      <- subset(dfw, Stat %in% c("RTElf", "MTEfl", "MTElf", "RTEfl"))
    dfw     <- subset(dfw, Dist >= xxlim[1] & Dist <= xxlim[2])
    
    dfw$Direction <- "lf"
    dfw[dfw$Stat == "MTEfl" | dfw$Stat == "RTEfl",]$Direction <- "fl"
    
    dfw$Encoding <- factor(dfw$Encoding, level=c("R","M"))
    
    # Plot the results..
    theme_set(theme_minimal())    
    ggplot(dfw, aes(x = Dist, y = Val, 
                    color = Direction, fill = Direction, linetype = Encoding)) +
      geom_hline(yintercept = 0, linetype = 1, color = "gray70") +
      #geom_vline(xintercept = vline[dsIdx], linetype = 3, size = 1, color = "red") +
      geom_smooth(span = 0.3, se = T, alpha = .1) +
      coord_cartesian(xlim = xxlim, ylim = yylim) +
      scale_fill_viridis(discrete = T, end=0.5, direction = -1) +
      scale_color_viridis(discrete = T, end=0.5, direction = -1) +
      xlab("Distance (body length)") +
      ylab("Predictive information (bits)") +
      labs(x = "", y = "Normalized transfer entropy") + 
      theme(legend.position  = c(0.25 , 0.85)) +
      theme(axis.text.y = element_text(angle = 90)) +    
      theme(legend.title = element_blank()) +
      theme(legend.text = element_text(size = fontSize3C)) +
      theme(text = element_text(size = fontSize3C))
    
    pdfName <- paste(odir, "localTEAveragedOverDistance-", ds, ".pdf", sep = "")
    ggsave(pdfName, height = pdfHeight, width = pdfHeight)
  }
  
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# Tandem speed / acceleration / duration
#---------------------------------------------------------------------------#
tandemOutput <- function(){
  idir <- "data/trajectory/processed/"  
  odir <- "img/trajectory/"  
  load(paste0(idir, "tandemStatistics.rda"))
  
  ## Acceleration analysis
  {
    ## randomly pickup 5000 datapoints for plotting
    {
      set.seed(100)
      sample.data <- c( sample(1:65013, 5000),
                        65013 + sample(1:89920, 5000),
                        65013 + 89920 + sample(1:76432, 5000),
                        65013 + 89920 + 76432 + sample(1:89920, 5000))
    }
    df.plot.temp <- df.all.speed[sample.data,]
    
    theme_set(theme_minimal())
    ggplot(df.plot.temp[df.plot.temp$Tandem,]) +
      geom_point(aes(x=Distance, y=LeaderAcc), col=viridis(3)[1], alpha = 0.1) +
      geom_point(aes(x=Distance, y=FollowerAcc), col=viridis(3)[2], alpha=0.1) +
      stat_smooth(aes(x=Distance, y=LeaderAcc), col=viridis(3)[1], method = "lm") +
      stat_smooth(aes(x=Distance, y=FollowerAcc), col=viridis(3)[2], method = "lm") +
      facet_wrap(~Species) +
      coord_cartesian(xlim = c(0.5,2.2), ylim=c(-1,1)) +
      labs(x = "Distance (BL)", y = "Acceleration (BL/sec2)") + 
      theme(legend.title = element_blank()) +
      theme(legend.text = element_text(size = fontSize3C)) +
      theme(text = element_text(size = fontSize3C))
    fname <- paste0(odir,"accelerationPerDistance.pdf")
    ggsave(fname, height = pdfHeight, width = pdfWidth)
    
    
    fname <- paste0(odir,"accelerationPerDistance.txt")
    sink(fname)
    cat("----------------------------------------------\n")
    cat("LMM, Acc vs Distance, for each role and species\n")
    cat("----------------------------------------------\n")
    
    df.stat <- df.all.speed[df.all.speed$Tandem,]
    for(i in 1:4){
      ds = Datasets[i]
      cat(ds, "\n")
      df.stat.temp <- df.stat[df.stat$Species == ds,]
      r <- lmer(LeaderAcc ~ Distance +(1|ID), data=df.stat.temp)
      cat("Leader\n")
      res <- summary(r)
      cat("slop+-sd =", res$coefficients[2,1:2], "\n")
      res <- Anova(r)
      cat("Chisq =", res$Chisq, "Df =", res$Df, "P =", res$`Pr(>Chisq)`, "\n")
      
      r <- lmer(FollowerAcc ~ Distance +(1|ID), data=df.stat.temp)
      cat("Follower\n")
      res <- summary(r)
      cat("slop+-sd =", res$coefficients[2,1:2], "\n")
      res <- Anova(r)
      cat("Chisq =", res$Chisq, "Df =", res$Df, "P =", res$`Pr(>Chisq)`, "\n")
    }
    sink()
  }
  
  ## Mean speed
  {
    df.temp <- df.ind.speed[df.ind.speed$Type=="Tandem", ]
    ggplot(df.temp, aes(x = Species, y=Speed, color=Role)) + 
      geom_boxplot() + 
      geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
      scale_color_viridis(discrete = T, end=0, begin=0.5) +
      coord_cartesian( ylim=c(0,.8)) +
      labs(x = "", y = "Mean speed (BL/sec)") + 
      theme(legend.position =  c(0.2,0.8)) +
      theme(legend.title = element_blank()) +
      theme(legend.text = element_text(size = fontSize3C)) +
      theme(text = element_text(size = fontSize3C))
    fname <- paste0(odir,"meanSpeedComparison.pdf")
    ggsave(fname, height = pdfHeight, width = pdfWidth)
    
    fname <- paste0(odir,"meanSpeedComparison.txt")
    sink(fname)
    cat("----------------------------------------------\n")
    cat("Paired t-test for mean speed during tandem run\n")
    cat("----------------------------------------------\n")
    cat("Note that P > 0.05 in Shapiro test for each\n")
    for(i in 1:4){
      ds = datasets[i]
      cat(ds, "\n")
      df.test <- df.temp[df.temp$Species == ds,]
      r <- t.test(df.test[df.test$Role=="Leader", "Speed"],
                  df.test[df.test$Role=="Follower", "Speed"],
                  paired = T)
      cat("t =", r$statistic, "df =", r$parameter, "p =", r$p.value, "\n")
      cat("estimate =", r$estimate, "\n")
    }
    sink()
  }
  
  ## Tandem duration
  {
    theme_set(theme_bw())
    ggplot(df.ind.speed[df.ind.speed$Role=="Leader" & 
                          df.ind.speed$Type=="Separation",],
           aes(x = Species, y=Duration/5, fill=Species)) + 
      geom_boxplot(alpha=0.2) + 
      geom_dotplot(binaxis = "y", stackdir = "center", dotsize=0.75, binwidth = 5, alpha=1) + 
      scale_fill_viridis(discrete = T, begin=1, end=0) +
      coord_cartesian(ylim=c(0,360)) +
      labs(x = "", y = "Separation duration (sec)") + 
      theme(legend.position  = "none") +
      theme(legend.title = element_blank()) +
      theme(legend.text = element_text(size = fontSize3C)) +
      theme(text = element_text(size = fontSize3C))
    fname <- paste0(odir,"SeparationDuration.pdf")
    ggsave(fname, height = pdfHeight, width = pdfWidth)
    
    theme_set(theme_bw())
    ggsurvplot(
      fit = survfit(Surv(Duration*0.2/60, Censor) ~ Species, 
                    type = "kaplan-meier", 
                    data = df.tandem.duration),
      conf.int = F, conf.int.style = "ribbon",
      xlab = "Duration (min)", 
      ylab = "Tandem Prob",
      palette = viridis(4)[c(4,1,2,3)],
      legend = c(0.8,0.9),
    )
    fname <- paste0(odir,"TandemDuration.pdf")
    ggsave(fname, height = pdfHeight/2, width = pdfWidth/2)
    
    m <- coxme(Surv(Duration*0.2/60, Censor) ~ Species + (1|ID), data = df.tandem.duration)
    summary(m)
    Anova(m)
    
    multicomparison<-glht(m,linfct=mcp(Species="Tukey"))
    summary(multicomparison)
  }
}
#---------------------------------------------------------------------------#

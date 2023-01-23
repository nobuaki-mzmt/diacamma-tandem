## Diacamma Tandem analysis
## NetworkOutput.R
## N. Mizumoto

# Note
#---------------------------------------------------------------------------#
# This script displays the results for trajectory analysis of tandem runs.
# The processed data at data/network/processed/ will be used.
# All the results will be stored at img/network/
#---------------------------------------------------------------------------#

rm(list = ls())
output_all()

#---------------------------------------------------------------------------#
{
  library(MASS)
  
  library(igraph)
  
  library(ggplot2)
  library(viridis)
  
  library(lme4)
  library(car)
  library(multcomp)
  
  library(tableHTML)
  
  library(stringr)
  
  library(Rmisc)
  
  pdfHeight  <- 7 
  pdfWidth   <- 10
  fontSize3C <- 28
}

output_all = function(){
  plotNetworks()
  plotSpecificNetworks()
  tableStateTransition()
  resultPropTandem()
  resultDegreeDstr()
  resultMotif()
  resultNetStat()
}

#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# Plot all tandem/transport networks
#---------------------------------------------------------------------------#
plotNetworks <- function(){
  idir  <- "data/network/processed/" 
  odir  <- "img/network/"
  load(paste0(idir,"networkEdges.RData"))
  
  dfEdges$Event = paste(dfEdges$Species, dfEdges$Colony, dfEdges$Treatment, sep="_")
  event = unique(dfEdges$Event)
  
  species = unique(dfEdges$Species)
  col5 = viridis(5, option="magma")
  names(col5) = species
  
  for(tt in c("tandem", "transport")){
    for(ss in species){
      if(tt=="transport" && ss!="temnothorax_rugatulus"){next;}
      pdf(paste0(odir, ss, "_", tt, ".pdf" ), width =12, height=7)
      par(mfrow=c(4,6), mar = c(1, 1, 1, 1) )
      for(i in 1:length(event)){
        dftemp <- subset(dfEdges, Species == ss & Event==event[i] & Type==tt)
        if(dim(dftemp)[1] == 0){next;}
        dfnet <- graph.data.frame(dftemp[,c("Leader", "Follower")],
                                  #vertices =subset(dfNodes,Colony==dftemp$Colony[1])[,2],
                                  directed=T)
        plot(dfnet,
             vertex.label = NA,
             vertex.size=6,
             edge.arrow.size  = 0.2,
             edge.width       = 1,
             edge.curved      = 0.2,
             vertex.color=col5[ss],
             layout=layout.fruchterman.reingold)
        title(paste0(dftemp$Colony[1], "-",  dftemp$Treatment[1],
                    " N=", vcount(dfnet),
                    " E=", ecount(dfnet)))
        
      }
      dev.off()
    }
  }
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# Plot specific tandem/transport networks for figs
#---------------------------------------------------------------------------#
plotSpecificNetworks <- function(){
  idir  <- "data/network/processed/" 
  odir  <- "img/network/"
  load(paste0(idir,"networkEdges.RData"))
  
  dfEdges$Event = paste(dfEdges$Species, dfEdges$Colony, dfEdges$Treatment, sep="_")
  unique(dfEdges$Event)
  
  pdf(paste0(odir, "represent_net", ".pdf" ), width =12, height=7)
  par(mfrow=c(1,6), mar = c(1, 1, 1, 1) )
  species <- c("diacamma_indicum_DI-320_1", "diacamma_sp_SB_1", "temnothorax_albipennis_1_4",
    "temnothorax_nylanderi_1_1_2", "temnothorax_rugatulus_6_3")
  
  for(i in 1:5){
    dftemp <- subset(dfEdges, Event== species[i])
    if(i == 5){
      dftemp <- subset(dfEdges, Event== species[i] & Type=="tandem")
    }
    dfnet <- graph.data.frame(dftemp[,c("Leader", "Follower")],
                              directed=T)
    if(i < 3){col3 = viridis(3)[3]}
    if(i > 2){col3 = viridis(3)[1]}
    plot(dfnet,
         vertex.label = NA,
         vertex.size=6,
         edge.arrow.size  = 0.2,
         edge.width       = 1,
         edge.curved      = 0.2,
         vertex.color=col3,
         layout=layout.fruchterman.reingold)
    title(paste0(dftemp$Colony[1], "-",  dftemp$Treatment[1],
                 " N=", vcount(dfnet),
                 " E=", ecount(dfnet)))
    
    if(i == 5){
      dftemp <- subset(dfEdges, Event== species[i] & Type=="transport")
      dfnet <- graph.data.frame(dftemp[,c("Leader", "Follower")],
                                directed=T)
      
      plot(dfnet,
           vertex.label = NA,
           vertex.size=6,
           edge.arrow.size  = 0.2,
           edge.width       = 1,
           edge.curved      = 0.2,
           vertex.color=viridis(3)[2],
           layout=layout.fruchterman.reingold)
      title(paste0(dftemp$Colony[1], "-",  dftemp$Treatment[1],
                   " N=", vcount(dfnet),
                   " E=", ecount(dfnet)))
    }
  }
  dev.off()
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# Output table of state transitions
#---------------------------------------------------------------------------#
tableStateTransition <- function(){
  idir  <- "data/network/processed/" 
  odir  <- "img/network/"
  fname <- paste(idir, "transitionProb.RData", sep = "")
  load(fname)
  for(i in 3:7){
    dftransitionsum[,i] <- round(as.numeric(dftransitionsum[,i]),4)
  }
  
  fname <- paste(odir, "transitionProb.html", sep = "")
  write_tableHTML(tableHTML(dftransitionsum), file = fname)
  
  dftemp <- dftransitioncount
  dftemp <- dftemp[dftemp$Before == "tandemFollower" | dftemp$Before == "transportFollower",]
  for(i in 3:7){
    dftemp[,i] <- as.numeric(dftemp[,i])
    dftemp[is.na(dftemp[,i]),i] <- 0
  }
  dftemp$recruited = dftemp$tandemLeader + dftemp$transportLeader
  dftemp$notrecruited = dftemp$tandemFollower + dftemp$transportFollower + dftemp$end
  dftemp = dftemp[c(1:2,7,9, 3, 5),c(1,2,8,9)]
  
  ## statistical analysis
  {
    fname <- paste(odir, "/transitionProbStatistics.txt", sep = "")  
    sink(fname)
    cat("----------------------------------------------\n")
    cat("Tests using Fisher's test\n")
    cat("Comparing the probability to be recuruited\n")
    cat("----------------------------------------------\n")
    
    cat("temno  tandem    vs temno   transport\n")
    r <- fisher.test(dftemp[c(1,2),3:4]) # temno  tandem    vs temno   transport
    cat("P: ", r$p.value, "\n")
    
    cat("temno  tandem    vs dia ind tandem\n")
    r <- fisher.test(dftemp[c(1,3),3:4]) # temno  tandem    vs dia ind tandem
    cat("P: ", r$p.value, "\n")
    
    cat("temno  tandem    vs dia sp  tandem\n")
    r <- fisher.test(dftemp[c(1,4),3:4]) # temno  tandem    vs dia sp  tandem
    cat("P: ", r$p.value, "\n")
    
    cat("temno  transport vs dia ind tandem \n")
    r <- fisher.test(dftemp[c(2,3),3:4]) # temno  transport vs dia ind tandem 
    cat("P: ", r$p.value, "\n")
    
    cat("temno  transport vs dia sp  tandem \n")
    r <- fisher.test(dftemp[c(2,4),3:4]) # temno  transport vs dia sp  tandem 
    cat("P: ", r$p.value, "\n")
    
    cat("dia sp tandem    vs dia ind  tandem \n")
    r <- fisher.test(dftemp[c(3,4),3:4]) # dia sp tandem    vs dia ind  tandem 
    cat("P: ", r$p.value, "\n")
    
    cat("\n----------------------------------\nFig.S2\n")
    cat("tem alb tandem    vs tem rug  tandem \n")
    r <- fisher.test(dftemp[c(5,1),3:4]) # dia sp tandem    vs dia ind  tandem 
    cat("P: ", r$p.value, "\n")
    
    cat("tem alb tandem    vs tem nyl  tandem \n")
    r <- fisher.test(dftemp[c(5,6),3:4]) # dia sp tandem    vs dia ind  tandem 
    cat("P: ", r$p.value, "\n")
    
    cat("tem alb tandem    vs dia ind  tandem \n")
    r <- fisher.test(dftemp[c(5,3),3:4]) # dia sp tandem    vs dia ind  tandem 
    cat("P: ", r$p.value, "\n")
    
    cat("tem alb tandem    vs dia sp   tandem \n")
    r <- fisher.test(dftemp[c(5,4),3:4]) # dia sp tandem    vs dia ind  tandem 
    cat("P: ", r$p.value, "\n")
    
    cat("tem alb tandem    vs tem rug  transport \n")
    r <- fisher.test(dftemp[c(5,2),3:4]) # dia sp tandem    vs dia ind  tandem 
    cat("P: ", r$p.value, "\n")
    
    cat("\n")
    cat("tem nyl tandem    vs tem rug  tandem \n")
    r <- fisher.test(dftemp[c(6,1),3:4]) # dia sp tandem    vs dia ind  tandem 
    cat("P: ", r$p.value, "\n")
    
    cat("tem nyl tandem    vs dia ind  tandem \n")
    r <- fisher.test(dftemp[c(6,3),3:4]) # dia sp tandem    vs dia ind  tandem 
    cat("P: ", r$p.value, "\n")
    
    cat("tem nyl tandem    vs dia sp   tandem \n")
    r <- fisher.test(dftemp[c(6,4),3:4]) # dia sp tandem    vs dia ind  tandem 
    cat("P: ", r$p.value, "\n")
    
    cat("tem nyl tandem    vs tem rug  transport \n")
    r <- fisher.test(dftemp[c(6,2),3:4]) # dia sp tandem    vs dia ind  tandem 
    cat("P: ", r$p.value, "\n")
   sink()
  }
  
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# Comparison of prop of individuals involved in tandem/transport
#---------------------------------------------------------------------------#
resultPropTandem <- function(){
  idir  <- "data/network/processed/" 
  odir  <- "img/network/"
  
  load(paste0(idir,"tandemAggregate.Rdata"))
  
  ## Recruit Proportion output
  {
    df.prop.recruit <- data.frame(
      Species = c("diacamma_indicum", "diacamma_sp", "temnothorax_rugatulus"),
      TandemInvolve.Mean   = tapply(df.sum$TandemInvolved/df.sum$ColonySize, df.sum$Species, mean),
      TandemInvolve.SD     = tapply(df.sum$TandemInvolved/df.sum$ColonySize, df.sum$Species, sd),
      TandemLeader.Mean    = tapply(df.sum$TandemLeader/df.sum$ColonySize, df.sum$Species, mean),
      TandemLeader.SD      = tapply(df.sum$TandemLeader/df.sum$ColonySize, df.sum$Species, sd),
      TandemFollowerMean   = tapply(df.sum$TandemFollower/df.sum$ColonySize, df.sum$Species, mean),
      TandemFollowerSD     = tapply(df.sum$TandemFollower/df.sum$ColonySize, df.sum$Species, sd),
      CarryInvolved.Mean   = tapply(df.sum$CarryInvolved/df.sum$ColonySize, df.sum$Species, mean),
      CarryInvolved.SD     = tapply(df.sum$CarryInvolved/df.sum$ColonySize, df.sum$Species, sd),
      Carry.Mean           = tapply(df.sum$Carrier/df.sum$ColonySize, df.sum$Species, mean),
      Carry.SD             = tapply(df.sum$Carrier/df.sum$ColonySize, df.sum$Species, sd),
      Carried.Mean         = tapply(df.sum$BeCarried/df.sum$ColonySize, df.sum$Species, mean),
      Carried.SD           = tapply(df.sum$BeCarried/df.sum$ColonySize, df.sum$Species, sd)
    )
    write.table(df.prop.recruit, paste0(odir, "RecruitProp.txt"))
  }
  
  ## modification for plot
  df.temp = rbind(
    data.frame(
      df.sum[, c("Species", "Colony", "Treatment")],
      Role = "Leader",
      Type = "tandem",
      Prop = df.sum$TandemLeader/df.sum$ColonySize
    ),
    data.frame(
      df.sum[, c("Species", "Colony", "Treatment")],
      Role = "Follower",
      Type = "tandem",
      Prop = df.sum$TandemFollower/df.sum$ColonySize
    ),
    data.frame(
      df.sum[, c("Species", "Colony", "Treatment")],
      Role = "Leader",
      Type = "transport",
      Prop = df.sum$Carrier/df.sum$ColonySize
    ),
    data.frame(
      df.sum[, c("Species", "Colony", "Treatment")],
      Role = "Follower",
      Type = "transport",
      Prop = df.sum$BeCarried/df.sum$ColonySize
    )
  )
    
  
  df.temp <- subset(df.temp, Species=="temnothorax_rugatulus" | Type != "transport")
  
  ## plot results
  dodge_width = 0.5
  pdfHeight = 4
  pdfWidth = 8
  theme_set(theme_minimal())    
  
  df.temp$Species <- factor(df.temp$Species, levels = c("diacamma_sp", "diacamma_indicum", "temnothorax_rugatulus"))
  
  ggplot(df.temp,
         aes(x=interaction(Species,Type), y=Prop, fill=Role)) +
    geom_boxplot(width = .5, position=position_dodge(dodge_width)) +
    geom_dotplot(binaxis='y', stackdir='center',
                 dotsize=1, binwidth=0.01, position=position_dodge(dodge_width)) +
    scale_fill_viridis(discrete = T, option = "magma", begin=0)+
    labs(x = "", y = "Proportion of involved individuals") +
    lims(y = c(0,1)) +
    theme(legend.position  = c(0.65 , 0.85)) +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size = fontSize3C)) +
    theme(text = element_text(size = fontSize3C))
  
  pdfName <- paste(odir, "/tandemPropComp.pdf", sep = "")  
  ggsave(pdfName, height = pdfHeight, width = pdfWidth)
  
  ## statistical analysis
  fname <- paste(odir, "/tandemPropCompStatistics.txt", sep = "")  
  sink(fname)
  cat("----------------------------------------------\n")
  cat("Tests using GLMM with binomial error\n")
  cat("Response variable: involved or not for each worker\n")
  cat("Explanatory variable: species\n")
  cat("Random effect: emigration event\n")
  cat("----------------------------------------------\n")
  df.sum$Treatment = paste0(df.sum$Colony,df.sum$Treatment)
  
  species <- unique(df.sum$Species)
  for(ss in species){
    df.temp <- subset(df.sum, Species != ss)
    if(ss == "temnothorax_rugatulus"){
      
      cat(unique(df.temp$Species)[1], "(tandem)     ", "vs", unique(df.temp$Species)[2], "(tandem)", "\n")
    
      r <- glmer(cbind(TandemLeader, (ColonySize-TandemLeader) ) ~ 
                   Species + (1|Treatment), 
                 family = binomial, data=df.temp)
      res <- Anova(r)
      cat("    Leader  , P:", round(res$`Pr(>Chisq)`,4), "Chisq:", round(res$Chisq,4),  "\n")
      
      r <- glmer(cbind(TandemFollower, (ColonySize-TandemLeader) ) ~ 
                   Species + (1|Treatment), 
                 family = binomial, data=df.temp)
      res <- Anova(r)
      cat("    Follower, P:", round(res$`Pr(>Chisq)`,4), "Chisq:", round(res$Chisq,4),  "\n")
      
    } else {
      
      # tandem vs tandem
      cat(unique(df.temp$Species)[1], "(tandem)      ", "vs", unique(df.temp$Species)[2], "(tandem)", "\n")
      
      r <- glmer(cbind(TandemLeader, (ColonySize-TandemLeader) ) ~ 
                   Species + (1|Treatment), 
                 family = binomial, data=df.temp)
      res <- Anova(r)
      cat("    Leader  , P:", round(res$`Pr(>Chisq)`,4), "Chisq:", round(res$Chisq,4),  "\n")
      
      r <- glmer(cbind(TandemFollower, (ColonySize-TandemLeader) ) ~ 
                   Species + (1|Treatment), 
                 family = binomial, data=df.temp)
      res <- Anova(r)
      cat("    Follower, P:", round(res$`Pr(>Chisq)`,4), "Chisq:", round(res$Chisq,4),  "\n")
      
      # tandem vs transport
      df.temp$Comp1 <- c(df.temp[df.temp$Species == "temnothorax_rugatulus",]$Carrier,
                         df.temp[df.temp$Species != "temnothorax_rugatulus",]$TandemLeader)
      
      cat(unique(df.temp$Species)[1], "(transport)", "vs", unique(df.temp$Species)[2], "(tandem)" , "\n")
      
      r <- glmer(cbind(Comp1, (ColonySize-Comp1) ) ~ 
                   Species + (1|Treatment), 
                 family = binomial, data=df.temp)
      res <- Anova(r)
      cat("    Leader  , P:", round(res$`Pr(>Chisq)`,4), "Chisq:", round(res$Chisq,4),  "\n")
      
      df.temp$Comp2 <- c(df.temp[df.temp$Species == "temnothorax_rugatulus",]$BeCarried,
                         df.temp[df.temp$Species != "temnothorax_rugatulus",]$TandemFollower)
      r <- glmer(cbind(Comp2, (ColonySize-Comp2) ) ~ 
                   Species + (1|Treatment), 
                 family = binomial, data=df.temp)
      res <- Anova(r)
      cat("    Follower, P:", round(res$`Pr(>Chisq)`,4), "Chisq:", round(res$Chisq,4),  "\n")
      
      
      
    }
  }
  sink()
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# Degree distribution
#---------------------------------------------------------------------------#
resultDegreeDstr <- function(){
  idir  <- "data/network/processed/" 
  odir  <- "img/network/"
  pdfHeight = 3
  pdfWidth = 4.5
  
  fname = paste0(idir,"networkStat.RData")
  load(fname)
  
  # plot distribution
  {
    dfdegree$Genus <- dfdegree$Species
    dfdegree$Genus[str_detect(dfdegree$Genus, "diacamma")] <- "diacamma"
    dfdegree$Genus[str_detect(dfdegree$Genus, "temnothorax")] <- "temnothorax"
    dfdegree$Genus[dfdegree$Type == "transport"] <- "temnothorax_carrying"
    
    theme_set(theme_minimal())
    
    df.plot <- subset(dfdegree,nNodes>20 & degree<6)
    sumrepdat <- summarySE(df.plot, measurevar = "prob",
                           groupvars=c("Direction", "degree", "Genus"))
    
    df.plot.plot <- subset(df.plot, Direction == "in")
    sumrepdat.plot <- subset(sumrepdat, Direction == "in")
    df.plot.plot$Genus <- factor(df.plot.plot$Genus, levels = c("temnothorax", "temnothorax_carrying", "diacamma"))
    sumrepdat.plot$Genus <- factor(sumrepdat.plot$Genus, levels = c("temnothorax", "temnothorax_carrying", "diacamma"))
    ggplot(data = df.plot.plot, aes(x=as.factor(degree), y=prob, fill=Genus)) + 
      geom_point(aes(x = as.numeric(degree)+.95, colour = Genus),
                 position = position_jitter(width = .05), size = .25, shape = 20)+
      #geom_boxplot(aes(x = as.factor(degree), y = prob, fill = Direction),
      #             outlier.shape = NA, alpha = .5, width = .1, colour = "black") +
      geom_line(data = sumrepdat.plot, 
                aes(x = as.numeric(degree)+1.05, y = prob, 
                    group = Genus, colour = Genus), linetype = 1)+
      geom_point(data = sumrepdat.plot, 
                 aes(x = as.numeric(degree)+1.05, y = prob,
                     group = Genus, colour = Genus), shape = 18) +
      geom_errorbar(data = sumrepdat.plot, 
                    aes(x = as.numeric(degree)+1.05, y = prob,
                        group = Genus, colour = Genus,
                        ymin = prob-sd, ymax = prob+sd), width = .05)+
      scale_fill_viridis(discrete = T)+
      scale_color_viridis(discrete = T) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "Degree", y = "Proportion of nodes") +
      theme_bw()+
      theme(legend.position  = c(0.75 , 0.85)) +
      theme(legend.title = element_blank()) +
      theme(legend.text = element_text(size = fontSize3C/2)) +
      theme(text = element_text(size = fontSize3C/2))
    pdfName <- paste(odir, "/indegree.pdf", sep = "")  
    ggsave(pdfName, height = pdfHeight, width = pdfWidth)
      
    df.plot.plot <- subset(df.plot, Direction == "out")
    sumrepdat.plot <- subset(sumrepdat, Direction == "out")
    df.plot.plot$Genus <- factor(df.plot.plot$Genus, levels = c("temnothorax", "temnothorax_carrying", "diacamma"))
    sumrepdat.plot$Genus <- factor(sumrepdat.plot$Genus, levels = c("temnothorax", "temnothorax_carrying", "diacamma"))
    ggplot(data = df.plot.plot, aes(x=as.factor(degree), y=prob, fill=Genus)) + 
      geom_point(aes(x = as.numeric(degree)+.95, colour = Genus),
                 position = position_jitter(width = .05), size = .25, shape = 20)+
      #geom_boxplot(aes(x = as.factor(degree), y = prob, fill = Direction),
      #             outlier.shape = NA, alpha = .5, width = .1, colour = "black") +
      geom_line(data = sumrepdat.plot, 
                aes(x = as.numeric(degree)+1.05, y = prob, 
                    group = Genus, colour = Genus), linetype = 1)+
      geom_point(data = sumrepdat.plot, 
                 aes(x = as.numeric(degree)+1.05, y = prob,
                     group = Genus, colour = Genus), shape = 18) +
      geom_errorbar(data = sumrepdat.plot, 
                    aes(x = as.numeric(degree)+1.05, y = prob,
                        group = Genus, colour = Genus,
                        ymin = prob-sd, ymax = prob+sd), width = .05)+
      scale_fill_viridis(discrete = T)+
      scale_color_viridis(discrete = T) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "Degree", y = "Proportion of nodes") + 
      theme_bw()+
      theme(legend.position  = c(0.75 , 0.85)) +
      theme(legend.title = element_blank()) +
      theme(legend.text = element_text(size = fontSize3C/2)) +
      theme(text = element_text(size = fontSize3C/2))
    pdfName <- paste(odir, "/outdegree.pdf", sep = "")  
    ggsave(pdfName, height = pdfHeight, width = pdfWidth)
  }
  
  # comparison of prop of passive workers
  df.ind.degree$passive = df.ind.degree$outdegree == 0 & df.ind.degree$indegree == 1
  
  dfDegreeSum = dfNetStat[,1:6]
  dfDegreeSum$ID = paste(dfDegreeSum$Species, dfDegreeSum$Colony,
                         dfDegreeSum$Treatment, dfDegreeSum$Type, sep="-")
  dfDegreeSum$Event = paste(dfDegreeSum$Species, dfDegreeSum$Colony,
                              dfDegreeSum$Treatment, sep="-")
  dfDegreeSum = dfDegreeSum[order(dfDegreeSum$ID),]
  
  dfDegreeSum$passive = tapply(df.ind.degree$passive, df.ind.degree$ID, sum)
  dfDegreeSum$active = tapply(!df.ind.degree$passive, df.ind.degree$ID, sum)
  
  dfDegreeSum$Genus <- dfDegreeSum$Species
  dfDegreeSum$Genus[str_detect(dfDegreeSum$Genus, "diacamma")] <- "diacamma"
  dfDegreeSum$Genus[str_detect(dfDegreeSum$Genus, "temnothorax")] <- "temnothorax"
  dfDegreeSum$Genus[dfDegreeSum$Type == "transport"] <- "temnothorax_carrying"
  
  
  r <- glmer(cbind(passive, active) ~ Genus + (1|Event), 
             family = binomial, data=dfDegreeSum)
  res <- Anova(r)
  cat("Chisq =", res$Chisq, "; df =", res$Df, "; P =", res$`Pr(>Chisq)`, "\n")
  
  multicomparison<-glht(r,linfct=mcp(Genus="Tukey"))
  res <- summary(multicomparison)
  
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# Motif analysis
#---------------------------------------------------------------------------#
resultMotif <- function(){
  idir  <- "data/network/processed/" 
  odir  <- "img/network/"
  pdfHeight = 2.5
  pdfWidth = 8
  
  fname = paste0(idir,"networkStat.RData")
  load(fname)
  
  dfmotif$Genus <- dfmotif$Species
  dfmotif$Genus[str_detect(dfmotif$Genus, "diacamma")] <- "diacamma"
  dfmotif$Genus[str_detect(dfmotif$Genus, "temnothorax")] <- "temnothorax"
  dfmotif$Genus[dfmotif$Type == "transport"] <- "temnothorax_carrying"
  
  df.plot <- subset(dfmotif,nNodes>20 & (motif.id == 4 | motif.id == 6))
  df.plot$Genus <- factor(df.plot$Genus, levels = c("diacamma", "temnothorax", "temnothorax_carrying"))
  dodge_width = 0.75
  ggplot(df.plot, aes(x=as.factor(motif.id), y=motif.prop, fill=Genus))+
    geom_boxplot(alpha = 0.4)+
    geom_dotplot(binaxis='y', stackdir='center',
                 dotsize=2, binwidth=0.01, position=position_dodge(dodge_width)) +
    scale_fill_manual(values = viridis(3)[c(3,1,2)]) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Motif", y = "Proportion of subgraphs") + 
    theme_classic()+
    theme(legend.position  = c(0.2 , 0.85)) +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size = fontSize3C/2)) +
    theme(text = element_text(size = fontSize3C/2))
  pdfName <- paste(odir, "/motif.pdf", sep = "")  
  ggsave(pdfName, height = pdfHeight, width = pdfWidth)
  
  plot(graph.isocreate(3,2)) 
  plot(graph.isocreate(3,4)) 
  plot(graph.isocreate(3,6)) 
  
  
  # statistical analysis
  dfmotif$ID = paste(dfmotif$Species, dfmotif$Colony, dfmotif$Treatment, dfmotif$Type, sep="-")
  dfmotif$Event = paste(dfmotif$Species, dfmotif$Colony, dfmotif$Treatment, sep="-")
  dfmotif = dfmotif[order(dfmotif$ID),]
  dfmotif_stat <- dfmotif[dfmotif$motif.id==0, c(1:5,11:13)]
  dfmotif_stat$MotifTotal = tapply(dfmotif$motif.count, dfmotif$ID,  sum, na.rm=T)
  dfmotif_stat$Motif4 = dfmotif[dfmotif$motif.id == 4, "motif.count"]
  dfmotif_stat$Motif6 = dfmotif[dfmotif$motif.id == 6, "motif.count"]
  
  dfmotif_stat = dfmotif_stat[dfmotif_stat$nNodes > 20,]

  fname <- paste(odir, "/tandemMotifComparison.txt", sep = "")  
  sink(fname)
  cat("----------------------------------------------\n")
  cat("Tests using GLMM\n")
  cat("Response variable: Motif is target or not\n")
  cat("Explanatory variable: Recruitment types (Dia tandem vs Temno tandem vs Temno carrying)\n")
  cat("Random effect: Emigration event\n")
  cat("----------------------------------------------\n")
  
  cat("Motif id = 4 (A -> B -> C)\n")
  r <- glmer(cbind(Motif4, (MotifTotal-Motif4) ) ~ Genus + (1|Event), 
             family = binomial, data=dfmotif_stat)
  res <- Anova(r)
  cat("Chisq =", res$Chisq, "; df =", res$Df, "; P =", res$`Pr(>Chisq)`, "\n")
  
  multicomparison<-glht(r,linfct=mcp(Genus="Tukey"))
  res <- summary(multicomparison)
  
  for(i in 1:length(capture.output(res))){
    cat( capture.output(res)[i], "\n" )
  }
  
  cat("Motif id = 6 (B <- A -> C)\n")
  r <- glmer(cbind(Motif6, (MotifTotal-Motif6) ) ~ Genus + (1|Event), 
             family = binomial, data=dfmotif_stat)
  res <- Anova(r)
  cat("Chisq =", res$Chisq, "; df =", res$Df, "; P =", res$`Pr(>Chisq)`, "\n")
  
  multicomparison<-glht(r,linfct=mcp(Genus="Tukey"))
  res <- summary(multicomparison)
  for(i in 1:length(capture.output(res))){
    cat( capture.output(res)[i], "\n" )
  }
  sink()
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# Some network statistics
#---------------------------------------------------------------------------#
resultNetStat <- function(){
  idir  <- "data/network/processed/" 
  odir  <- "img/network/"
  
  fname = paste0(idir,"networkStat.RData")
  load(fname)
  
  ## density
  dfNetStat$Genus <- dfNetStat$Species
  dfNetStat$Genus[str_detect(dfNetStat$Genus, "diacamma")] <- "diacamma"
  dfNetStat$Genus[str_detect(dfNetStat$Genus, "temnothorax")] <- "temnothorax"
  dfNetStat$Genus[dfNetStat$Type == "transport"] <- "temnothorax_carrying"
  
  ggplot(dfNetStat, aes(x=Genus, y=Density))+ geom_boxplot()
  
  # statistical analysis
  
  fname <- paste(odir, "/tandemNetworkComparison.txt", sep = "")  
  sink(fname)
  cat("----------------------------------------------\n")
  cat("Tests using LMM\n")
  cat("Response variable: Network density\n")
  cat("Explanatory variable: Recruitment types (Dia tandem vs Temno tandem vs Temno carrying)\n")
  cat("Random effect: Colony nested within species\n")
  cat("----------------------------------------------\n")
  
  r <- lmer(Density~Genus + (1|Species/Colony), data=dfNetStat) 
  res <- Anova(r)
  cat("Chisq =", res$Chisq, "; df =", res$Df, "; P =", res$`Pr(>Chisq)`, "\n")
  multicomparison<-glht(r,linfct=mcp(Genus="Tukey"))
  res <- summary(multicomparison)
  for(i in 1:length(capture.output(res))){
    cat( capture.output(res)[i], "\n" )
  }
  sink()
}
#---------------------------------------------------------------------------#
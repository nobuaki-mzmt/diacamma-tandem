## Diacamma Tandem analysis
## NetworkProcessing.R
## N. Mizumoto

# Note
#---------------------------------------------------------------------------#
# This script processes the raw data for tandem network of five ant species:
# -Diacamma    sp.        (data from this study)
# -Diacamma    indicum    (data from Kolay and Annagiri 2015, 10.1098/rsos.150104)
# -Temnothorax rugatulus  (data from Valentini et al.   2020, 10.1098/rspb.2019.2950)
# -Temnothorax albipennis (data from Richardson et al.  2018, 10.1098/rspb.2017.2726)
# -Temnothorax nylanderi  (data from Richardson et al.  2021, 10.1038/s42003-021-02048-7).
# The raw data are at data/network/raw/
# The processed data will be stored at data/network/processed/
# The processed data will be used in NetworkOutputs.R for displaying results.
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
rm(list = ls())
library(igraph)
process_all()
#---------------------------------------------------------------------------#


#---------------------------------------------------------------------------#
process_all <- function(){
  extractNetworks()
  TandemAggreagte(OnlyOldGoodMovement = F)
  RoleTransition()
  networkStatistics()
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# This function combine all extractNetworks functions and make one dataframe
#---------------------------------------------------------------------------#
extractNetworks <- function(){
  iodir  <- "data/network/processed/" 

  dftemp1 = dftemp2 <- data.frame()
  
  extractNetworks_tem.rug.()
  load(paste0(iodir,"networkEdges_tem_rug.RData"))
  dftemp1 = rbind(dftemp1, dfEdges)
  dftemp2 = rbind(dftemp2, dfNodes)
  
  extractNetworks_tem.alb.()
  load(paste0(iodir,"networkEdges_tem_alb.RData"))
  dftemp1 = rbind(dftemp1, dfEdges)
  dftemp2 = rbind(dftemp2, dfNodes)
  
  extractNetworks_tem.nyl.()
  load(paste0(iodir,"networkEdges_tem_nyl.RData"))
  dftemp1 = rbind(dftemp1, dfEdges)
  dftemp2 = rbind(dftemp2, dfNodes)
  
  extractNetworks_dia.ind.()
  load(paste0(iodir,"networkEdges_dia_ind.RData"))
  dftemp1 = rbind(dftemp1, dfEdges)
  dftemp2 = rbind(dftemp2, dfNodes)
  
  extractNetworks_dia.sp.()
  load(paste0(iodir,"networkEdges_dia_sp.RData"))
  dftemp1 = rbind(dftemp1, dfEdges)
  dftemp2 = rbind(dftemp2, dfNodes)
  
  dfEdges = dftemp1
  dfNodes = dftemp2
  
  fname <- paste(iodir, "networkEdges.RData", sep = "")
  save(dfEdges, dfNodes, file = fname)
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# This function reads data of Temnothorax rugatulus from 
#   Valentini et al. 2020 (10.1098/rspb.2019.2950)
# Then, convert it into .Rdata file with Edge list and Node list for this study.
#---------------------------------------------------------------------------#
extractNetworks_tem.rug. <- function() {
  idir  <- "data/network/raw/temnothorax_rugatulus/"
  odir  <- "data/network/processed/" 
  colonies        <- c(6, 208,3004)
  treatments      <- list(c(1:5), c(2:5), c(1:5))
  
  # Load aggregated statistics.. 
  fname <- paste(idir, "aggregatedStatistics.RData",sep="")
  load(fname)
  
  dfEdges <- data.frame(stringsAsFactors = !T)

  for(cc in colonies) {
    ccIdx   <- which(colonies == cc)
    dfNodes <- subset(df, colony == cc)[, c(19:20)]
    
    for (tt in treatments[[ccIdx]]) {
      # Load raw data..
      fname <- paste(idir, "colony", cc, "treatment", tt, ".txt", sep = "")
      dfw <- read.table(fname, header = T, fill = T)
      dfw <- dfw[-1,] # Remove initial time entry
      dfw$Time <- as.character(dfw$Time)
      
      # Extract edges
      for (j in 1:(dim(dfw)[1] - 1)) {
        # Find edge
        x1 <- as.POSIXlt(dfw$Time[j], format = "%H:%M:%S")
        x2 <- as.POSIXlt(dfw$Time[j + 1], format = "%H:%M:%S")
        
        x1Id <- subset(dfNodes, colorId == dfw$Ant[j])$id
        x2Id <- subset(dfNodes, colorId == dfw$Ant[j + 1])$id

        if (x1 == x2) {
          if (x1Id == x2Id) { next }
          
          # Interaction 1: LB <--> FB
          if (dfw$Behavior[j] == "LB" && dfw$Behavior[j + 1] == "FB") {
            dfw$LocationOrigin[j]
            dfw$SourceDestination[j]
            ttype <- "tandem"
            
            dfEdges <-
              rbind(dfEdges,
                    c(cc, tt, x1Id, x2Id, ttype = "tandem", 
                      From = dfw$LocationOrigin[j],
                      Into = dfw$SourceDestination[j],
                      Time = unclass(as.POSIXct(x1))[1]),
                    stringsAsFactors = !T)

          } else if (dfw$Behavior[j] == "FB" &&
                     dfw$Behavior[j + 1] == "LB") {
            
            dfEdges <-
              rbind(dfEdges,
                    c(cc, tt, x2Id, x1Id, ttype = "tandem", 
                      From = dfw$LocationOrigin[j],
                      Into = dfw$SourceDestination[j],
                      Time = unclass(as.POSIXct(x1))[1]),
                    stringsAsFactors = !T)

          }
          
          # Interaction 2: TA <--> BC
          if (dfw$Behavior[j] == "TA" &&
              dfw$Behavior[j + 1] == "BC") {
            dfEdges <-
              rbind(dfEdges, c(cc, tt, x1Id, x2Id, ttype = "transport", 
                               From = dfw$SourceDestination[j],
                               Into = dfw$LocationOrigin[j],
                               Time = unclass(as.POSIXct(x1))[1]),
                    stringsAsFactors = !T)
          } else if (dfw$Behavior[j] == "BC" &&
                     dfw$Behavior[j + 1] == "TA") {
            dfEdges <-
              rbind(dfEdges, c(cc, tt, x2Id, x1Id,ttype = "transport", 
                               From = dfw$SourceDestination[j],
                               Into = dfw$LocationOrigin[j],
                               Time = unclass(as.POSIXct(x1))[1]),
                    stringsAsFactors = !T)
          }
          
          # Interaction 3: TAO) --> BCO)
          if (dfw$Behavior[j] == "TAO" &&
              dfw$Behavior[j + 1] == "BCO") {
            dfEdges <-
              rbind(dfEdges, c(cc, tt, x1Id, x2Id, ttype = "transport", 
                               From = dfw$LocationOrigin[j],
                               Into = dfw$SourceDestination[j],
                               Time = unclass(as.POSIXct(x1))[1]),
                    stringsAsFactors = !T)
          } else if (dfw$Behavior[j] == "BCO" &&
                     dfw$Behavior[j + 1] == "TAO") {
            dfEdges <-
              rbind(dfEdges, c(cc, tt, x2Id, x1Id, ttype = "transport", 
                               From = dfw$LocationOrigin[j],
                               Into = dfw$SourceDestination[j],
                               Time = unclass(as.POSIXct(x1))[1]),
                    stringsAsFactors = !T)
          }
        }
      }
    }
  }
  dfNodes <- df[,c("colony", "id")]
  names(dfNodes) <- c("Colony", "ID")
  
  dfNodes$Species = "temnothorax_rugatulus"
  
  names(dfEdges) <-
    c("Colony", "Treatment", "Leader", "Follower", "Type", "From", "Into", "Time")
  dfEdges$Species = "temnothorax_rugatulus"
  fname <- paste(odir, "networkEdges_tem_rug.RData", sep = "")
  save(dfEdges, dfNodes, file = fname)
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# This function reads data of Temnothorax albipennis from
#   Richardson et al. 2018 (10.1098/rspb.2017.2726).
# Then, convert it into .Rdata file with Edge list and Node list for this study.
# Note that
# 1) only info for tandems (no carrying)
# 2) no information about colony members who did not engage in tandems
#---------------------------------------------------------------------------#
extractNetworks_tem.alb. <- function() {
  idir  <- "data/network/raw/temnothorax_albipennis/"
  odir  <- "data/network/processed/" 
  
  # Load aggregated statistics.. 
  fname <- paste(idir, "rspb20172726_si_002.txt",sep="")
  dfw = read.table(fname, header = T, fill = T)

  # remove failed tandems
  dfw <- subset(dfw, Fail.Success == "Success")
  
  dfEdges = data.frame(
    dfw[,c("Colony_ID", "Replicate", "Leader", "Follower")],
    Type = "tandem",
    dfw[,c("Origin", "target_choice", "StartTime")],
    Species = "temnothorax_albipennis"
  )
  
  names(dfEdges) <-
    c("Colony", "Treatment", "Leader", "Follower", "Type", "From", "Into", "Time", "Species")

  dfNodes = data.frame()
  Colonies = unique(dfw$Colony_ID)
  for(cc in Colonies){
    dftemp = subset(dfw, Colony_ID == cc)
    ID = unique( c(dftemp$Leader, dftemp$Follower) )
    
    dfNodes = rbind(dfNodes, 
      dfNodesTemp = data.frame(
        cc,
        ID
      )
    )
  }
  names(dfNodes) <- c("Colony", "ID")
  rownames(dfNodes) <- NULL
  dfNodes$Species = "temnothorax_albipennis"
  
  fname <- paste(odir, "networkEdges_tem_alb.RData", sep = "")
  save(dfEdges, dfNodes, file = fname)
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# This function reads data of Temnothorax nylanderi from
#   Richardson et al. 2021 (10.1038/s42003-021-02048-7).
# Then, convert it into .Rdata file with Edge list and Node list for this study.
# Note that
# 1) only info for tandems (no carrying)
# 2) no information about colony members who did not engage in tandems
# 3) no distinction between forward tandem, reverse tandem, and failed tandem 
#    thus, included all of these for this dataset.
#---------------------------------------------------------------------------#
extractNetworks_tem.nyl. <- function() {
  idir  <- "data/network/raw/temnothorax_nylanderi/"
  odir  <- "data/network/processed/" 
  
  # Load aggregated statistics.. 
  fname <- paste(idir, "Tandem_Running,Emigrations_1-5.txt",sep="")
  dfw = read.table(fname, header = T, fill = T, sep=",")
  
  # remove 5th emigration as it removed some members
  dfw <- subset(dfw, Treatment == "Control")
  
  dfEdges = data.frame(
    dfw[,c("Colony", "Emigration", "Leader", "Follower")],
    Type = "tandem",
    dfw[,c("Leaves_From", "Arrives_at", "Time_Start")],
    Species = "temnothorax_nylanderi"
  )
  
  names(dfEdges) <-
    c("Colony", "Treatment", "Leader", "Follower", "Type", "From", "Into", "Time", "Species")
  
  ## Note that no information about other colony members
  
  dfNodes = data.frame()
  Colonies = unique(dfw$Colony)
  for(cc in Colonies){
    dftemp = subset(dfw, Colony == cc)
    ID = unique( c(dftemp$Leader, dftemp$Follower) )
    
    dfNodes = rbind(dfNodes, 
                    dfNodesTemp = data.frame(
                      cc,
                      ID
                    )
    )
  }
  names(dfNodes) <- c("Colony", "ID")
  rownames(dfNodes) <- NULL
  dfNodes$Species = "temnothorax_nylanderi"
  
  fname <- paste(odir, "networkEdges_tem_nyl.RData", sep = "")
  save(dfEdges, dfNodes, file = fname)
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# This function reads data of Diacamma indicum from
#   Kolay and Annagiri 2015 (10.1098/rsos.150104).
# Then, convert it into .Rdata file with Edge list and Node list for this study.
#---------------------------------------------------------------------------#
extractNetworks_dia.ind. <- function() {
  idir  <- "data/network/raw/diacamma_indicum/"
  odir  <- "data/network/processed/" 
  
  # Load aggregated statistics.. 
  fname <- paste(idir, "CR-only.csv",sep="")
  dfw = read.csv(fname, header = T)
  
  dfEdges = data.frame(
    dfw[,c("Colony")],
    Treatment = 1,
    dfw[,c("ï..Leader", "Follower")],
    Type = "tandem",
    From = "old",
    Into = "new",
    dfw[,"Time"],
    Species = "diacamma_indicum"
  )
  
  names(dfEdges) <-
    c("Colony", "Treatment", "Leader", "Follower", "Type", "From", "Into", "Time", "Species")
  
  dfNodes = data.frame()
  colonies = unique(dfw$Colony)
  for(cc in colonies){
    dftemp = subset(dfw, Colony==cc)
    tandem.involved = unique(c(dftemp$ï..Leader, dftemp$Follower))
    non.tandem.involved = 1:(dftemp$ColonySize[1] - length(tandem.involved))
    ID = c(tandem.involved, non.tandem.involved)
    dfNodes = rbind(dfNodes, data.frame(cc, ID, "diacamma_indicum"))
  }
  names(dfNodes) = c("Colony", "ID", "Species")
  
  ## Note that no information about other colony members
  
  fname <- paste(odir, "networkEdges_dia_ind.RData", sep = "")
  save(dfEdges, dfNodes, file = fname)
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# This function reads raw data of Diacamma experiments (in this study).
# Then, convert it into .Rdata file with Edge list and Node list.
#---------------------------------------------------------------------------#
extractNetworks_dia.sp. <- function() {
  idir  <- "data/network/raw/diacamma_sp/"
  odir  <- "data/network/processed/" 
  colonies        <- c("K4", "K8", "SA", "SB")
  colonysize <- c(115, 127, 111, 129)+1
  names(colonysize) = colonies
  
  dfEdges <- data.frame(stringsAsFactors = !T)
  dfNodes <- data.frame(stringsAsFactors = !T)
  for(cc in colonies) {
    fname <- paste(idir, "Tandem_New_", cc, "_fixed.csv", sep = "")
    dfw <- read.csv(fname)
    dfw$Colony = cc
    From = Into = rep("old", dim(dfw)[1])
    From[dfw$Direction == "Out"] = "new"
    Into[dfw$Direction == "In"] = "new"
    
    dfEdges <- rbind(dfEdges, 
                     data.frame(
                       dfw[, c("Colony")],
                       Treatment = 1,
                       dfw[, c("Leader_ID", "Follower_ID")],
                       Type = "tandem",
                       From, Into,
                       Time = dfw[, c("Time")]
                     )
    )
    dfNodes <- rbind(dfNodes, data.frame(Colony = cc, ID = 1:colonysize[cc]-1))
  }
  names(dfEdges) = c("Colony", "Treatment", "Leader", "Follower", "Type", "From", "Into", "Time")
  dfEdges$Species = "diacamma_sp"
  dfNodes$Species = "diacamma_sp"
  fname <- paste(odir, "networkEdges_dia_sp.RData", sep = "")
  save(dfEdges, dfNodes, file = fname)
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# This function reads networkEdges.RData,
# and then, count the number of tandem/carry events for each migration
#---------------------------------------------------------------------------#
TandemAggreagte <- function(OnlyOldGoodMovement = F){
  iodir  <- "data/network/processed/" 
  load(paste0(iodir,"networkEdges.RData"))
  
  # remove Temnothorax albipennis and Temnothorax nylanderi
  dfEdges <- subset(dfEdges, Species!="temnothorax_albipennis" & Species!="temnothorax_nylanderi")
  
  if(OnlyOldGoodMovement){
    dfEdges = subset(dfEdges, (From=="old" | From=="old?") & (Into=="new" | Into=="good") )
  }
  
  dfEdges$Event = paste(dfEdges$Species, dfEdges$Colony, dfEdges$Treatment, sep="_")
  event = unique(dfEdges$Event)
  
  df.sum <- data.frame()
  for(ee in event){
    dftemp1 <- subset(dfEdges, Event==ee)
    dftemp2 <- subset(dfNodes, Colony == dftemp1[1, "Colony"] & Species == dftemp1[1, "Species"])
    
    tandem.leader = (unique((dftemp1[dftemp1$Type=="tandem",]$Leader)))
    tandem.follower = (unique((dftemp1[dftemp1$Type=="tandem",]$Follower)))
    tandem.involved = unique(c(tandem.leader, tandem.follower))
    carrier = (unique((dftemp1[dftemp1$Type=="transport",]$Leader)))
    becarried = (unique((dftemp1[dftemp1$Type=="transport",]$Follower)))
    carry.involved = unique(c(carrier, becarried))
    
    df.sum = rbind(df.sum,
      data.frame(
        Event = ee,
        Species = dftemp1[1, "Species"],
        Colony = dftemp1[1, "Colony"],
        ColonySize = dim(dftemp2)[1],
        Treatment = dftemp1[1, "Treatment"],
        TandemNum = sum(dftemp1$Type == "tandem"),
        CarryNum = sum(dftemp1$Type == "transport"),
        TandemLeader = length(tandem.leader),
        TandemFollower = length(tandem.follower),
        TandemInvolved = length(tandem.involved),
        Carrier = length(carrier),
        BeCarried = length(becarried),
        CarryInvolved = length(carry.involved)
      )
    )
  }
  if(OnlyOldGoodMovement){
    fname <- paste(iodir, "tandemAggregate_sub.RData", sep = "")
  } else {
    fname <- paste(iodir, "tandemAggregate.RData", sep = "")
  }
  save(df.sum, file = fname)
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# Role transition probability
#---------------------------------------------------------------------------#
RoleTransition <- function(){
  iodir  <- "data/network/processed/" 
  load(paste0(iodir,"networkEdges.RData"))
  
  # remove Temnothorax albipennis and Temnothorax nylanderi
  # dfEdges <- subset(dfEdges, Species!="temnothorax_albipennis" & Species!="temnothorax_nylanderi")
  
  dfEdges$Event = paste(dfEdges$Species, dfEdges$Colony, dfEdges$Treatment, sep="_")
  event = unique(dfEdges$Event)
  
  dftransition = data.frame()
  for(ee in event){
    cat(paste("Analyze event:", ee, "\n"))
    
    dfEdgesTemp <- subset(dfEdges, Event==ee)
    dfNodesTemp <- subset(dfNodes, Colony == dfEdgesTemp[1, "Colony"] & Species == dfEdgesTemp[1, "Species"])
    
    for(ii in dfNodesTemp$ID){
      dftemp <- dfEdgesTemp[dfEdgesTemp$Leader==ii | dfEdgesTemp$Follower==ii,]
      if(dim(dftemp)[1] == 0){ cat(paste(ii, "no tandem or carry\n")); next; }
      Role = rep("Follower", dim(dftemp)[1])
      Role[dftemp$Leader == ii] <- "Leader"
      Role = paste(dftemp$Type, Role, sep="")
      
      if(length(Role)>1){ 
        for(rr in 2:length(Role)-1){
          Transition = paste(Role[rr], Role[rr+1], sep="-")
          dftransition <- rbind(dftransition, 
                                c(unlist(dftemp[1, c("Species", "Colony", "Treatment")]),
                                  ii, Role[rr], Role[rr+1], Transition))
        }
      }
      Transition = paste(Role[length(Role)], "end", sep="-")
      dftransition <- rbind(dftransition, 
                            c(unlist(dftemp[1, c("Species", "Colony", "Treatment")]),
                              ii, Role[length(Role)], "end", Transition))
    } 
  }
  names(dftransition) = c("Species", "Colony", "Treatment", "ID", "Before", "After", "Transition")
  
  
  RoleNames <- c("tandemLeader", "tandemFollower", "transportLeader", "transportFollower")
  Species <- unique(dftransition$Species)
  dftransitionsum <- data.frame()
  for(ss in Species){
    for(rr in RoleNames){
      dftemp = subset(dftransition, Species==ss & Before==rr)
      transitions <- table(dftemp$After)/(dim(dftemp)[1])
      transitions <- transitions[c(RoleNames, "end")]
      dftransitionsum <- rbind(dftransitionsum, c(ss, rr, transitions))
    }
  }
  names(dftransitionsum) = c("Species", "Before", c(RoleNames, "end"))
  
  RoleNames <- c("tandemLeader", "tandemFollower", "transportLeader", "transportFollower")
  Species <- unique(dftransition$Species)
  dftransitionsum <- data.frame()
  for(ss in Species){
    for(rr in RoleNames){
      dftemp = subset(dftransition, Species==ss & Before==rr)
      transitions <- table(dftemp$After)
      transitions <- transitions[c(RoleNames, "end")]
      dftransitionsum <- rbind(dftransitionsum, c(ss, rr, transitions))
    }
  }
  names(dftransitionsum) = c("Species", "Before", c(RoleNames, "end"))
  dftransitioncount <- dftransitionsum
  
  fname <- paste(iodir, "transitionProb.RData", sep = "")
  save(dftransitionsum, dftransitioncount, file = fname)
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# Network statitics
#---------------------------------------------------------------------------#
networkStatistics <- function(){
  iodir  <- "data/network/processed/" 
  load(paste0(iodir,"networkEdges.RData"))
  
  dfEdges$Event = paste(dfEdges$Species, dfEdges$Colony, dfEdges$Treatment, sep="_")
  event = unique(dfEdges$Event)
  
  types = sort(unique(dfEdges$Type))
  
  dfNetStat = dfmotif = dfdegree = df.ind.degree <- data.frame()
  for(ee in event){
    print(paste(which(event==ee), "/", length(event), ee))
    for(tt in types){
      dftemp <- subset(dfEdges, Type==tt & Event == ee)
      if(dim(dftemp)[1]==0){next;}
      dfnet <- graph.data.frame(dftemp[,c("Leader", "Follower")], directed=T)
      
      ## statistics per network
      {
        nEdges = ecount(dfnet)
        nNodes = vcount(dfnet)
        Density = graph.density(dfnet)
        Transitivity = transitivity(dfnet)
        Reciprocity = reciprocity(dfnet)
        Diameter = diameter(dfnet)
        dfNetStat <- rbind(dfNetStat,
                           data.frame(
                             dftemp[1,c("Species", "Colony", "Treatment","Type")],
                             nEdges,
                             nNodes,
                             Density,
                             Transitivity,
                             Reciprocity,
                             Diameter
                             
                           )
        )
      }
      
      ## degree distribution
      {
        outdegree = degree.distribution(dfnet, mode = "out", cumulative = T)
        indegree = degree.distribution(dfnet, mode = "in", cumulative = T)
        dfdegree = rbind(dfdegree, 
                            data.frame(
                              dftemp[1,c("Species", "Colony", "Treatment","Type")],
                              nNodes = nNodes,
                              Direction = "out",
                              degree = 1:length(outdegree)-1,
                              prob = outdegree
                            )
        )
        dfdegree = rbind(dfdegree, 
                         data.frame(
                           dftemp[1,c("Species", "Colony", "Treatment","Type")],
                           nNodes = nNodes,
                           Direction = "in",
                           degree = 1:length(indegree)-1,
                           prob = indegree
                         )
        )
      }
      
      # individual level data
      {
        outdegree = degree(dfnet,mode="out")
        indegree = degree(dfnet,mode="in")
        
        df.ind.degree = rbind(df.ind.degree,
                              data.frame(
                                dftemp[1,c("Species", "Colony", "Treatment","Type")],
                                name = names(indegree), outdegree, indegree
                              )
                              )
  
      }
      
      ## motif
      {
        c_real <- graph.motifs(dfnet,3)
        c_rand <- data.frame()
        outdegree = degree(dfnet,mode="out")
        indegree = degree(dfnet,mode="in")
        
        for(i in 1:1000){
          g_rand <- degree.sequence.game(outdegree,indegree)
          c_rand <- rbind(c_rand,graph.motifs(g_rand,3))
        }
        
        names(c_rand)<-0:15
        
        z.score = (c_real-colMeans(c_rand))/sqrt(colMeans(c_rand**2)-colMeans(c_rand)**2) 
        sp.score = z.score / sqrt(sum(z.score^2, na.rm=T))
        
        dfmotif <- rbind(dfmotif,
          data.frame( dftemp[1,c("Species", "Colony", "Treatment","Type")],
                      nNodes,
                      motif.id = 0:15,
                      motif.count = c_real,
                      motif.prop = c_real/sum(c_real, na.rm = T),
                      z.score, sp.score
          )
        )
      }
    }
  }
  fname = paste0(iodir,"networkStat.RData")
  save(dfNetStat, dfmotif, dfdegree, df.ind.degree, file = fname)
}
#---------------------------------------------------------------------------#

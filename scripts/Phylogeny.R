## Diacamma Tandem analysis
## Phylogeny.R
## N. Mizumoto

# Note
#---------------------------------------------------------------------------#
# This script is for phylogenetic comparative analysis of tandem behavior
# This does both processing data and output the results
# Phylogeny data comes from Nelsen et al. 2018 (10.1073/pnas.1719794115)
# Tandem data is from Reeves and Moreau 2019 (10.26049/ASP77-2-2019-10)
#  and add tamdem information of Diacamma and Hypoponera
#---------------------------------------------------------------------------#

rm(list = ls())
#---------------------------------------------------------------------------#
{
  library(phytools)
  library(stringr)
  library(extrafont)
  loadfonts()
  
  processPhylogeny()
  simplifiedPhylogeny()
  outputPhylogeny()
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# This function reads phylogeny and foraging information, and then
# output genus-level tree and its tip's tandem information
#---------------------------------------------------------------------------#
processPhylogeny <- function(){
  iodir <- "data/phylogeny/"  
  f.name <- ("Dryad_Supplementary_File_7_ML_TREE_treepl_185.tre")
  full.tree = read.tree(paste0(iodir, f.name))
  
  ## convert species-level tree into genus-level tree
  {
    label <- full.tree$tip.label
    genus = str_sub(label, 1, str_locate(label, "_")[,1]-1)
    
    remove.overwrap = rep(F, length(genus))
    for(i in 2:length(genus)){
      if(genus[i] %in% genus[2:i-1]){
        remove.overwrap[i] = T
      }
    }
    
    genus.tree = drop.tip(full.tree, label[remove.overwrap])
    label <- genus.tree$tip.label
    genus = str_sub(label, 1, str_locate(label, "_")[,1]-1)
    genus.tree$tip.label = genus
  }
  
  ## Remove tips without foraging information
  {
    f.name <- ("reevesmoreau-antforaging-asp2019-electronicsupplement-2.csv")
    d <- read.csv(paste0(iodir, f.name))
    names(d) = c("Genus", "Species", "Foraging")
    
    d = subset(d, Foraging != "?")
    
    genus.with.foraging.info = unique(d$Genus)
    label <- genus.tree$tip.label
    
    genus.tree.foraging = drop.tip(genus.tree, label[!(label %in% genus.with.foraging.info)])
    
    df <- data.frame(genus = genus.tree.foraging$tip.label)
    tandem.genus = d[d$Foraging %in% c("1", "0&1", "0&1&5"),"Genus"]
    df$tandem = (df$genus %in% tandem.genus)*1
    df[df$genus=="Diacamma",]$tandem = 1
    df[df$genus=="Hypoponera",]$tandem = 1
    
    tandem = df$tandem
    names(tandem) = df$genus
  }
  
  f.name <- ("genusLevelTandemTree.RData")
  save(genus.tree.foraging, tandem, file = paste0(iodir, f.name))
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# This function creates simplified phylogeny for presentation purposes
#---------------------------------------------------------------------------#
simplifiedPhylogeny <- function(){
  iodir <- "data/phylogeny/"  
  f.name <- ("genusLevelTandemTree.RData")
  load(paste0(iodir, f.name))

  ## integrate genera into tribe
  {
    Dorylinae <- c("Ooceraea", "Leptanilloides","Simopone","Cylindromyrmex","Acanthostichus","Zasphinctus","Parasyscia", "Lioponera", "Aenictus","Dorylus", "Cheliomyrmex",  "Labidus", "Nomamyrmex", "Eciton",  "Neivamyrmex")
    simple.tree <-  drop.tip(genus.tree.foraging, Dorylinae[Dorylinae!="Eciton"])
    simple.tree$tip.label[simple.tree$tip.label=="Eciton"] = "Dorylinae"
    
    Dolichoderinae <- c("Liometopum", "Tapinoma", "Technomyrmex", "Dolichoderus", "Azteca", "Forelius", "Dorymyrmex", "Leptomyrmex","Anonychomyrma","Philidris","Ochetellus", "Iridomyrmex","Linepithema")
    simple.tree <-  drop.tip(simple.tree, Dolichoderinae[Dolichoderinae!="Liometopum"])
    simple.tree$tip.label[simple.tree$tip.label=="Liometopum"] = "Dolichoderinae"
    
    Attini <- c("Lachnomyrmex", "Blepharidatta","Wasmannia","Allomerus", "Acanthognathus", "Daceton","Orectognathus","Apterostigma",
                "Myrmicocrypta","Cyatta", "Cyphomyrmex", "Mycetophylax", "Sericomyrmex", "Trachymyrmex", "Acromyrmex", "Atta",
                "Basiceros","Strumigenys", "Cephalotes", "Pheidole")
    simple.tree <-  drop.tip(simple.tree, Attini[Attini!="Pheidole"])
    simple.tree$tip.label[simple.tree$tip.label=="Pheidole"] = "Attini"
    
    Amblyoponini <- c("Onychomyrmex","Amblyopone","Prionopelta","Myopopone","Stigmatomma","Mystrium" )
    simple.tree <-  drop.tip(simple.tree, Amblyoponini[Amblyoponini!="Mystrium"])
    simple.tree$tip.label[simple.tree$tip.label=="Mystrium"] = "Amblyoponini"
    
    Solenopsidini = c("Stegomyrmex","Megalomyrmex","Monomorium","Solenopsis","Myrmicaria")
    simple.tree <-  drop.tip(simple.tree, Solenopsidini[Solenopsidini!="Solenopsis"])
    simple.tree$tip.label[simple.tree$tip.label=="Solenopsis"] = "Solenopsidini"
    
    Stenammini = c("Veromessor","Novomessor","Aphaenogaster","Goniomma","Messor")
    simple.tree <-  drop.tip(simple.tree, Stenammini[Stenammini!="Messor"])
    simple.tree$tip.label[simple.tree$tip.label=="Messor"] = "Stenammini"
    
    Lasiini = c("Cladomyrma","Lasius","Myrmecocystus","Zatania","Prenolepis","Euprenolepis","Paratrechina" )
    simple.tree <-  drop.tip(simple.tree, Lasiini[Lasiini!="Lasius"])
    simple.tree$tip.label[simple.tree$tip.label=="Lasius"] = "Lasiini"
    
    Melophorini = c("Notostigma","Melophorus","Prolasius")
    simple.tree <-  drop.tip(simple.tree, Melophorini[Melophorini!="Notostigma"])
    simple.tree$tip.label[simple.tree$tip.label=="Notostigma"] = "Melophorini"
    
    Formicini = c("Proformica","Cataglyphis","Rossomyrmex","Polyergus","Formica" )
    simple.tree <-  drop.tip(simple.tree, Formicini[Formicini!="Formica"])
    simple.tree$tip.label[simple.tree$tip.label=="Formica"] = "Formicini"
    
    Pseudomyrmecinae = c("Tetraponera", "Pseudomyrmex")
    simple.tree <-  drop.tip(simple.tree, Pseudomyrmecinae[Pseudomyrmecinae!="Pseudomyrmex"])
    simple.tree$tip.label[simple.tree$tip.label=="Pseudomyrmex"] = "Pseudomyrmecinae"
    
    Ectatommini = c("Gnamptogenys","Rhytidoponera","Ectatomma")
    simple.tree <-  drop.tip(simple.tree, Ectatommini[Ectatommini!="Ectatomma"])
    simple.tree$tip.label[simple.tree$tip.label=="Ectatomma"] = "Ectatommini"
    
    Pogonomyrmecini=c("Patagonomyrmex","Pogonomyrmex")
    simple.tree <-  drop.tip(simple.tree, Pogonomyrmecini[Pogonomyrmecini!="Pogonomyrmex"])
    simple.tree$tip.label[simple.tree$tip.label=="Pogonomyrmex"] = "Pogonomyrmecini"
    
    Proceratiini = c("Proceratium","Discothyrea", "Probolomyrmex") # Probolomyrmex is a different tribe but within a clade
    simple.tree <-  drop.tip(simple.tree, Proceratiini[Proceratiini!="Proceratium"])
    simple.tree$tip.label[simple.tree$tip.label=="Proceratium"] = "Proceratiini"
    
    simple.tree <-  drop.tip(simple.tree, c("Leptanilla","Proceratiini","Amblyoponini",
                                            "Dorylinae",  "Aneuretus", "Dolichoderinae", 
                                            "Nothomyrmecia","Myrmecia","Pseudomyrmecinae",
                                            "Ectatommini"))
    
    
    tandem.simple = tandem
    names(tandem.simple)[names(tandem.simple) == "Eciton"] = "Dorylinae"
    names(tandem.simple)[names(tandem.simple) == "Liometopum"] = "Dolichoderinae"
    names(tandem.simple)[names(tandem.simple) == "Pheidole"] = "Attini"
    names(tandem.simple)[names(tandem.simple) == "Mystrium"] = "Amblyoponini"
    names(tandem.simple)[names(tandem.simple) == "Solenopsis"] = "Solenopsidini"
    names(tandem.simple)[names(tandem.simple) == "Messor"] = "Stenammini"
    names(tandem.simple)[names(tandem.simple) == "Lasius"] = "Lasiini"
    names(tandem.simple)[names(tandem.simple) == "Notostigma"] = "Melophorini"
    names(tandem.simple)[names(tandem.simple) == "Formica"] = "Formicini"
    names(tandem.simple)[names(tandem.simple) == "Pseudomyrmex"] = "Pseudomyrmecinae"
    names(tandem.simple)[names(tandem.simple) == "Ectatomma"] = "Ectatommini"
    names(tandem.simple)[names(tandem.simple) == "Pogonomyrmex"] = "Pogonomyrmecini"
    names(tandem.simple)[names(tandem.simple) == "Proceratium"] = "Proceratiini"
    
    tandem.simple = tandem.simple[names(tandem.simple) %in% simple.tree$tip.label]
    
    simple.tree$tip.label[simple.tree$tip.label=="Brachymyrmex"] = "Myrmelachistini"
    names(tandem.simple)[names(tandem.simple) == "Brachymyrmex"] = "Myrmelachistini"
    
    simple.tree$tip.label[simple.tree$tip.label=="Oecophylla"] = "Oecophylini"
    names(tandem.simple)[names(tandem.simple) == "Oecophylla"] = "Oecophylini"
    
    simple.tree$tip.label[simple.tree$tip.label=="Anoplolepis"] = "Plagiolepidini"
    names(tandem.simple)[names(tandem.simple) == "Anoplolepis"] = "Plagiolepidini"
    
    simple.tree$tip.label[simple.tree$tip.label=="Gigantiops"] = "Gigantiopini"
    names(tandem.simple)[names(tandem.simple) == "Gigantiops"] = "Gigantiopini"
    
    simple.tree$tip.label[simple.tree$tip.label=="Myrmoteras"] = "Myrmoteratini"
    names(tandem.simple)[names(tandem.simple) == "Myrmoteras"] = "Myrmoteratini"
    
    simple.tree$tip.label[simple.tree$tip.label=="Myrmica"] = "Myrmicini"
    names(tandem.simple)[names(tandem.simple) == "Myrmica"] = "Myrmicini"
  }
  
  if(F){
    plot(simple.tree)
    tiplabels(pie = to.matrix(tandem.simple, sort(unique(tandem.simple))), piecol = c("blue", "red"), 
              cex = 0.5)
  }
  
  f.name <- paste0(iodir, "simplifiedTree.RData")
  save(simple.tree, tandem.simple, file=f.name)
}
#---------------------------------------------------------------------------#

#---------------------------------------------------------------------------#
# This function outputs the results
#---------------------------------------------------------------------------#
outputPhylogeny <- function(){
  idir <- "data/phylogeny/" 
  odir <- "img/phylogeny/"
  f.name <- ("genusLevelTandemTree.RData")
  load(paste0(idir, f.name))
  
  f.name <- ("simplifiedTree.RData")
  load(paste0(idir, f.name))
  
  ## Ancestral State Reconstruction
  fitER <- rerootingMethod(genus.tree.foraging, tandem, model = "ER")

  ## Output Genus level tree
  {
    genus.tree.foraging <- ladderize(genus.tree.foraging)
    f.name <- paste0(odir, "genusLevelTree.pdf")
    pdf(f.name, width=8, height=25, family = "PT Serif", paper = "a4")
    plot(genus.tree.foraging)
    nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
               piecol = c("white", "red"), cex = 0.3)
    tiplabels(pie = to.matrix(tandem, sort(unique(tandem))), piecol = c("white", "red"), 
              cex = 0.3)
    dev.off()
  }
  
  ## Output simplified tree
  {
    f.name <- paste0(odir, "simplifiedTree.pdf")
    pdf(f.name, width=5, height=9, family = "PT Serif", paper = "a4")
    dotTree(simple.tree, as.factor(tandem.simple),
            colors=setNames(c("white","red"),
                            c("0","1")),
            ftype="i",fsize=0.7)
    dev.off()
  }
}
#---------------------------------------------------------------------------#





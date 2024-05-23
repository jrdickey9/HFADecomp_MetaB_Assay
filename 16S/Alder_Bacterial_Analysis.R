#load packages
library(microbiome)
#library(vegetarian) #error
library(phyloseq); packageVersion("phyloseq") #1.42.0
library(ggplot2); packageVersion("ggplot2") #3.5.1
library(gdata)
library(ecodist)
library(vegan)
library("car") 
library(dplyr)
library(biomformat)
library("ape") #phylogenetic tools packages
library(phytools) #phylogenetic tools packages
library(castor)
library(doParallel) #used for UniFrac but meh not super necessary for small data sets.
library(lme4)

#set seed
set.seed(32)

#read in global environment
#save(list=ls(),file="bacterialHFAanalyses.rda")
#setwd("~/Desktop/Jonathan_Dickey_LabMac/HFA/Analyses/Bacteria")
#load(file="bacterialHFAanalyses.rda")

#Read in biom table from working directory
setwd("~/Desktop/Jonathan_Dickey_LabMac/HFA/bacterial_prelim/HFABac_paired_dadatableR1")
ASV_reads<-read_biom("HFABac-feature-table-FINAL.biom")
str(ASV_reads)
ASV_table<-as.data.frame(as.matrix(biom_data(ASV_reads)))
ASV_table[1:10,1:10]
otu_tab<-t(ASV_table)
dim(otu_tab) #428 samples by 30,989 
otu_tab[1:10, 1:10]

#Read in metadata file
setwd("/Users/lab/Desktop/Jonathan_Dickey_LabMac/HFA/Analyses/Bacteria")
meta_data<-read.csv("updated_meta_bacterial_data.csv",header=TRUE) #428 observations x 17 variables

#Read in taxonomic information. #First go into Excel, open tsv, and resave as csv.
setwd("/Users/lab/Desktop/Jonathan_Dickey_LabMac/HFA/bacterial_prelim/HFABac_paired_taxaR9")
bac.taxa<-read.csv(file="HFABac-Tax-FINAL.csv")
str(bac.taxa)
rownames(bac.taxa) #whole numbers
colnames(bac.taxa) # "Feature.ID" "Domain"     "Phylum"     "Class"      "Order"      "Family"     "Genus"      "species"
rownames(bac.taxa)<-bac.taxa[,1] #assign the feature IDs to the row names
bac.taxa<-bac.taxa[,2:8] #get rid of the first column. 
colnames(bac.taxa) # "Domain"  "Phylum"  "Class"   "Order"   "Family"  "Genus"   "species"

#Reading in phylogenetic tree
setwd("/Users/lab/Desktop/Jonathan_Dickey_LabMac/HFA/bacterial_prelim/HFABac-rooted-tree-paired-seqsR1")
dated.16Stree<-read.tree(file="HFABac-rooted-tree-FINAL.nwk") #reading in tree, uses ape package
is.rooted(dated.16Stree) #TRUE
sample_names(dated.16Stree) #NULL
dated.16Stree$tip.label #for example "b395acbe0db0ffc093d4f54a95fc7169" 

#Creating phyloseq object
str(bac.taxa) #data.frame 30,989 x 7
bac.taxa.m<-as.matrix(bac.taxa) #VERY NECESSARY TO DO, DON'T SKIP. 
str(bac.taxa.m)
colnames(bac.taxa.m)
str(otu_tab) #428 x 30,989
df.bacterial.OTU<-as.data.frame(otu_tab)
str(meta_data) #data.frame 428 x 17

#Matching row names
rownames(df.bacterial.OTU)<-as.character(meta_data[,2])
colnames(df.bacterial.OTU) #accession numbers
rownames(meta_data)<-as.character(meta_data[,2]) #sample names
rownames(bac.taxa.m)<-as.character(rownames(bac.taxa.m)) #these the accession numbers
samp.names<-as.character(meta_data[,2]) #for example, "HOK2_4_HOK1_0_JRD"

#Matching sample names (originally marked NULL)
sample_names(df.bacterial.OTU)<-samp.names
sample_names(meta_data)<-samp.names
sample_names(bac.taxa.m)<-samp.names
sample_names(dated.16Stree)<-samp.names

#Matching taxa names (originally marked NULL)
taxa_names(df.bacterial.OTU)<-colnames(df.bacterial.OTU)
taxa_names(dated.16Stree)<-colnames(df.bacterial.OTU)
taxa_names(meta_data)<-colnames(df.bacterial.OTU)
taxa_names(bac.taxa.m)<-colnames(df.bacterial.OTU)

#Here is the actual phyloseq object
Bacterial_phylo<-phyloseq(otu_table(df.bacterial.OTU, taxa_are_rows=FALSE), sample_data(meta_data), tax_table(bac.taxa.m), phy_tree(dated.16Stree))
Bacterial_phylo@otu_table[1:10,1:10]
dim(Bacterial_phylo@otu_table) #428 samples by 30,989 taxa
rowSums(Bacterial_phylo@otu_table)
colSums(Bacterial_phylo@otu_table)

#Examining taxonomic ranks to examine where chloroplasts and mitochondria are nested within
#102 archaea, 30,609 bacteria, 19 euks, 259 unassigned #i'm removing everything but bacteria
#table(tax_table(Bacterial_phylo)[, "Domain"], exclude = NULL)
#table(tax_table(Bacterial_phylo)[, "Phylum"], exclude = NULL) #2,885 unassigned
#table(tax_table(Bacterial_phylo)[, "Class"], exclude = NULL) #examine 
#table(tax_table(Bacterial_phylo)[, "Order"], exclude = NULL) #1081 reads assigned as chloroplast at this taxonomic rank
#table(tax_table(Bacterial_phylo)[, "Family"], exclude = NULL) #1920 reads assigned as mitochondria at this taxonomic rank
#table(tax_table(Bacterial_phylo)[, "Genus"], exclude = NULL) 
#table(tax_table(Bacterial_phylo)[, "species"], exclude = NULL) #2 unidentified

#Need to prune out ASVs that were marked as unknown by Silva-138 but blastn showed as a non-bacterial target with reads >= 1500 across all samples. See Filtering_Taxa.R
#to remove: #cbb2217c90d6c95e84f8f1fdd48db5a4, #90a94bd9fd9ecb49c346d4a57ba7ff8d, #87c919242ac555b5fbda1625eaf0922b, #77f8f8987d530682cdbf661e0dc7d381, #301335cde77a0d707c2f2a7638a41a1a, #8b9b07020402a01a66a91a91c73d7b3d, #40eab637f322d7219342605c49239de7, #02a1bd1f1d3321139eb7f402aad2b151, #60407e5fdc0519bb0151b05ee5f490b5, #f1a276064ed37471ec73a61b1a95a4e1, #2454f25042adf21a60dbb34305cae1e0, #dab85ceacce2c97933f6b663da1e1b44 

badTaxa<-c("cbb2217c90d6c95e84f8f1fdd48db5a4","90a94bd9fd9ecb49c346d4a57ba7ff8d","87c919242ac555b5fbda1625eaf0922b","77f8f8987d530682cdbf661e0dc7d381","301335cde77a0d707c2f2a7638a41a1a","8b9b07020402a01a66a91a91c73d7b3d","40eab637f322d7219342605c49239de7","02a1bd1f1d3321139eb7f402aad2b151","60407e5fdc0519bb0151b05ee5f490b5","f1a276064ed37471ec73a61b1a95a4e1","2454f25042adf21a60dbb34305cae1e0","dab85ceacce2c97933f6b663da1e1b44") #12 ASVs to remove

pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)] #requires dplyr (%in% syntax)
  return(prune_taxa(allTaxa, physeq))
}

dim(Bacterial_phylo@otu_table) #428 x 30,989
Bacterial_phylo_2 = pop_taxa(Bacterial_phylo, badTaxa)
dim(Bacterial_phylo_2@otu_table) #428 x 30,977
str(Bacterial_phylo_2@tax_table)
str(Bacterial_phylo_2@phy_tree) #works!

#Replace taxa names to short hand
taxa_names(Bacterial_phylo_2) <- paste0("ASV", seq(ntaxa(Bacterial_phylo_2)))

#remove other stuff
p1<-subset_taxa(Bacterial_phylo_2,  !Domain %in% "Unassigned") #requires dplyr (%in% syntax)
p2<-subset_taxa(p1,  !Domain %in% "Archaea") #be no more!
p3<-subset_taxa(p2,  !Domain %in% "Eukaryota") #sayonara euks!
p4<-subset_taxa(p3,  !Order %in% " Chloroplast") #Chlorplasts be gone! 
p5<-subset_taxa(p4,  !Family %in% " Mitochondria") #Adios amigos!

#Remove KatharoSeq samples for the time being and water column samples
HFABac_phylo<-subset_samples(p5, SHORT_SAMPLE_NAME != "JRD367" & SHORT_SAMPLE_NAME != "JRD368" & SHORT_SAMPLE_NAME != "JRD369" & SHORT_SAMPLE_NAME != "JRD370" & SHORT_SAMPLE_NAME != "JRD387" & SHORT_SAMPLE_NAME != "16S_JRD388" & SHORT_SAMPLE_NAME != "16S_JRD389" & SHORT_SAMPLE_NAME != "16S_JRD390" & SHORT_SAMPLE_NAME != "JRD399" & SHORT_SAMPLE_NAME != "JRD400" & SHORT_SAMPLE_NAME != "JRD401" & SHORT_SAMPLE_NAME != "JRD402" & SHORT_SAMPLE_NAME != "16S_JRD417" & SHORT_SAMPLE_NAME != "16S_JRD418" & SHORT_SAMPLE_NAME != "16S_JRD419" & SHORT_SAMPLE_NAME != "16S_JRD420" & SHORT_SAMPLE_NAME != "JRD397" & SHORT_SAMPLE_NAME != "JRD398" & SHORT_SAMPLE_NAME != "JRD403" & SHORT_SAMPLE_NAME != "16S_JRD404" & SHORT_SAMPLE_NAME != "JRD405" & SHORT_SAMPLE_NAME != "JRD406" & SHORT_SAMPLE_NAME != "JRD407" & SHORT_SAMPLE_NAME != "JRD408" & SHORT_SAMPLE_NAME != "JRD409" & SHORT_SAMPLE_NAME != "JRD410" & SHORT_SAMPLE_NAME != "JRD411" & SHORT_SAMPLE_NAME != "JRD412" & SHORT_SAMPLE_NAME != "JRD413" & SHORT_SAMPLE_NAME != "JRD414" & SHORT_SAMPLE_NAME != "JRD415" & SHORT_SAMPLE_NAME != "JRD416" & SHORT_SAMPLE_NAME != "JRD425" & SHORT_SAMPLE_NAME != "JRD426" & SHORT_SAMPLE_NAME != "JRD427" & SHORT_SAMPLE_NAME != "JRD428") #JRD397 and after are water column samples 

#Meta data file with new dimensions
dim(HFABac_phylo@sam_data) #392 x 17
HFABac_meta_pruned<-HFABac_phylo@sam_data

#ASV table with out zeros and singletons
HFABac_phylo.rm1<- prune_taxa(taxa_sums(HFABac_phylo) > 0, HFABac_phylo) #remove ASVs with zero abundance attributed to KatharoSeq and water columns
dim(HFABac_phylo.rm1@otu_table) #392 x 27,334
HFABac_phylo.rm2<- prune_taxa(taxa_sums(HFABac_phylo.rm1) > 1, HFABac_phylo.rm1) #remove singletons
dim(HFABac_phylo.rm2@otu_table) #392 x 27,021

HFABac_phylo.ASVtable.2<-HFABac_phylo.rm2@otu_table

#subset samples by decomposition day
# Subset5_20<-subset_samples(HFABac_phylo.rm2, DECOMP_DAY !="0")
# Subset5_20@sam_data$DECOMP_DAY
# min(rowSums(Subset5_20@otu_table))

#Rarefying (haulted due to low read depth)
min(rowSums(HFABac_phylo.ASVtable.2)) #480 Alright, so there was lots of host DNA in some of these samples. 
sort(x=rowSums(HFABac_phylo.ASVtable.2),decreasing=FALSE) #Consider samples to toss. It's a good amount below 10K. 

#remove low read depth samples to rarefy to 10K
HFABac_phylo_5K<-subset_samples(HFABac_phylo.rm2, LONGER_SAMPLE_NAME != "SEK2_5_HOK1_15_JRD" & LONGER_SAMPLE_NAME != "HOK2_3_HOK1_5_JRD" & LONGER_SAMPLE_NAME != "HOK1_5_HOK1_5_JRD" & LONGER_SAMPLE_NAME != "HOK2_2_HOK1_5_JRD" & LONGER_SAMPLE_NAME != "HOK1_4_HOK1_5_JRD" & LONGER_SAMPLE_NAME != "HOK1_4_HOK1_10_JRD" & LONGER_SAMPLE_NAME != "HOK1_5_HOK1_10_JRD" & LONGER_SAMPLE_NAME != "HOK1_2_HOK1_10_JRD" & LONGER_SAMPLE_NAME != "SEK2_1_SEK1_20_JRD" & LONGER_SAMPLE_NAME != "SEK2_5_HOK1_0_JRD" & LONGER_SAMPLE_NAME != "SEK1_4_HOK1_5_JRD" & LONGER_SAMPLE_NAME != "SEK1_4_SEK2_20_JRD" & LONGER_SAMPLE_NAME != "SEK1_4_SEK2_0_JRD" & LONGER_SAMPLE_NAME != "SEK2_5_SEK2_10_JRD" & LONGER_SAMPLE_NAME != "SEK2_4_HOK2_0_JRD" & LONGER_SAMPLE_NAME != "SEK1_3_SEK2_5_JRD" & LONGER_SAMPLE_NAME != "SEK1_2_SEK1_10_JRD" & LONGER_SAMPLE_NAME != "SEK1_5_SEK2_0_JRD"  & LONGER_SAMPLE_NAME != "HOK1_2_SEK1_0_JRD" & LONGER_SAMPLE_NAME != "SEK1_1_SEK2_0_JRD" & LONGER_SAMPLE_NAME != "SEK2_3_SEK1_5_JRD" & LONGER_SAMPLE_NAME != "SEK2_4_SEK2_10_JRD" & LONGER_SAMPLE_NAME != "SEK2_5_SEK1_5_JRD" & LONGER_SAMPLE_NAME != "HOK1_5_SEK1_0_JRD" & LONGER_SAMPLE_NAME != "SEK2_2_SEK2_0_JRD" & LONGER_SAMPLE_NAME != "SEK2_1_HOK2_20_JRD" & LONGER_SAMPLE_NAME != "SEK1_1_SEK2_5_JRD" & LONGER_SAMPLE_NAME != "SEK2_2_SEK1_5_JRD" & LONGER_SAMPLE_NAME != "SEK1_3_SEK2_0_JRD" & LONGER_SAMPLE_NAME != "SEK1_1_HOK2_0_JRD" & LONGER_SAMPLE_NAME != "SEK2_1_SEK2_10_JRD" & LONGER_SAMPLE_NAME != "HOK1_5_HOK2_15_JRD" & LONGER_SAMPLE_NAME != "SEK1_2_SEK2_5_JRD" & LONGER_SAMPLE_NAME != "HOK1_1_SEK1_0_JRD" & LONGER_SAMPLE_NAME != "HOK2_4_HOK1_5_JRD" & LONGER_SAMPLE_NAME != "HOK2_4_SEK2_5_JRD" & LONGER_SAMPLE_NAME != "SEK1_2_SEK2_10_JRD" & LONGER_SAMPLE_NAME != "SEK1_4_SEK2_5_JRD" & LONGER_SAMPLE_NAME != "HOK2_4_SEK1_0_JRD" & LONGER_SAMPLE_NAME != "HOK2_1_SEK2_10_JRD" & LONGER_SAMPLE_NAME != "HOK2_4_SEK2_10_JRD" & LONGER_SAMPLE_NAME != "HOK2_1_SEK1_5_JRD" & LONGER_SAMPLE_NAME != "SEK2_1_SEK2_0_JRD" & LONGER_SAMPLE_NAME != "SEK1_4_HOK2_15_JRD" & LONGER_SAMPLE_NAME != "HOK1_1_HOK2_5_JRD" & LONGER_SAMPLE_NAME != "HOK2_1_HOK2_15_JRD" & LONGER_SAMPLE_NAME != "SEK2_3_SEK2_10_JRD" & LONGER_SAMPLE_NAME != "SEK1_3_SEK2_10_JRD" & LONGER_SAMPLE_NAME != "SEK1_2_SEK2_0_JRD" & LONGER_SAMPLE_NAME != "HOK2_3_SEK1_0_JRD" & LONGER_SAMPLE_NAME != "HOK2_5_SEK1_0_JRD" & LONGER_SAMPLE_NAME != "HOK1_4_SEK1_0_JRD" & LONGER_SAMPLE_NAME != "HOK1_3_SEK1_0_JRD" & LONGER_SAMPLE_NAME != "HOK2_5_SEK2_10_JRD")




bac.tab.df<-as.data.frame(HFABac_phylo.ASVtable.2)
rdat<-rrarefy(bac.tab.df,487) #rarefy!

#Standardize abundances into proportions. 
rowSums(rdat) #Check if it worked
std.bac.tab<-decostand(bac.tab.df,"total") #replace object to rdat after rarefying. 
std.bac.tab[1:10,1:10] #Looks groovy! 

#Smush back together into single phyloseq object for w. UniFrac
bacterial.phylo.4analysis<-phyloseq(otu_table(std.bac.tab, taxa_are_rows=FALSE), sample_data(HFABac_meta_pruned), tax_table(HFABac_taxa_pruned), phy_tree(HFABac_phylo_pruned)) #this object contains non-rarefied data. Move to line 132 to execute rarefying lines and incorporate that ASV table into this object. 


#Create Quant Jaccard distance matrix
drdat<-vegdist(std.bac.tab,"jaccard")

#Revisiting meta data to build db-rda in an appropriate way
str(HFABac_meta_pruned) #392 x 15 
colnames(HFABac_meta_pruned)
HFABac_meta_pruned$ORIGIN_SITE<-as.factor(HFABac_meta_pruned$ORIGIN_SITE)
HFABac_meta_pruned$TREE_NUMBER<-as.factor(HFABac_meta_pruned$TREE_NUMBER)
HFABac_meta_pruned$DEPLOYMENT_SITE<-as.factor(HFABac_meta_pruned$DEPLOYMENT_SITE)
HFABac_meta_pruned$DECOMP_DAY<-as.factor(HFABac_meta_pruned$DECOMP_DAY)
HFABac_meta_pruned$HFA_2<-as.factor(HFABac_meta_pruned$HFA_2)
HFABac_meta_pruned$LOCATION_CODE<-as.factor(HFABac_meta_pruned$LOCATION_CODE)
HFABac_meta_pruned$HFA_CODE<-as.factor(HFABac_meta_pruned$HFA_CODE)
HFABac_meta_pruned$TERR_AQUA<-as.factor(HFABac_meta_pruned$TERR_AQUA)
HFABac_meta_pruned$LOCATION_CODEwDay0<-as.factor(HFABac_meta_pruned$LOCATION_CODEwDay0)
HFABac_meta_pruned$Deployed_River<-as.factor(HFABac_meta_pruned$Deployed_River)

basic.mod<-dbrda(drdat~HFABac_meta_pruned$ORIGIN_SITE+HFABac_meta_pruned$DECOMP_DAY+HFABac_meta_pruned$HFA_2)
h<-how(nperm=10000)
anova(basic.mod, permutations = h, by="margin") 
#                                   Df SumOfSqs      F    Pr(>F)
#  HFABac_meta_pruned$ORIGIN_SITE   3    3.019 3.1350 9.999e-05 ***
#  HFABac_meta_pruned$DECOMP_DAY    3    6.114 6.3488 9.999e-05 ***
#  HFABac_meta_pruned$HFA_2         1    0.224 0.6979    0.9091    
#  Residual                       383  122.950 

basic.mod.summary<-summary(basic.mod) 
basic.mod.summary$concont

#Removing objects in global environment that could help free up space and processing. 
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(p5)
rm(Bacterial_phylo)
rm(dated.16Stree)
rm(bac.taxa.m)
rm(ASV_reads)
rm(otu_tab)

#Trying a more complex model by accounting for variation among the 20 tree genotypes, this will help code information from origin site
mod2.noH2O<-dbrda(drdat~HFABac_meta_pruned$DECOMP_DAY+HFABac_meta_pruned$LOCATION_CODE+Condition(HFABac_meta_pruned$TREE_NUMBER))
anova(mod2.noH2O,permutations = h, by="margin")
#                                   Df SumOfSqs      F    Pr(>F)    
#HFABac_meta_pruned$DECOMP_DAY      4    8.057 6.2983 9.999e-05 ***
#HFABac_meta_pruned$LOCATION_CODE   4    1.216 0.9506    0.5891    
#Residual                         364  116.415                     

#I think from what I can gleam that an important model to run given Jackrel et al. 2019 "The origin, succession..." would be to investigate the fixed effects of decomposition day and deployment site. 
mod3.noH2O<-dbrda(drdat~HFABac_meta_pruned$DECOMP_DAY+HFABac_meta_pruned$DEPLOYMENT_SITE+Condition(HFABac_meta_pruned$TREE_NUMBER))
anova(mod3.noH2O,permutations = h, by="margin")
#                                    Df SumOfSqs      F    Pr(>F)  
#HFABac_meta_pruned$DECOMP_DAY        3    6.114 6.3672 9.999e-05 ***
#HFABac_meta_pruned$DEPLOYMENT_SITE   3    0.798 0.8306    0.8614    
#Residual                           365  116.834  

mod4.b.noH2O<-dbrda(drdat~HFABac_meta_pruned$DECOMP_DAY+HFABac_meta_pruned$DEPLOYMENT_SITE+HFABac_meta_pruned$LOCATION_CODEwDay0+Condition(HFABac_meta_pruned$TREE_NUMBER))
anova(mod4.b.noH2O,permutations = h, by="margin")
#                                        Df SumOfSqs      F    Pr(>F)
#HFABac_meta_pruned$DECOMP_DAY           3    4.169 4.2975 9.999e-05 *** #There are 5 groups. I'm not sure why I'm being penalized with n-2. 
#HFABac_meta_pruned$DEPLOYMENT_SITE      3    1.047 1.0791    0.2704  #I dont know why but this should be 4 DF, there's 5 groups.   
#HFABac_meta_pruned$LOCATION_CODEwDay0   4    1.191 0.9204    0.6839  #There are 6 groups, again not sure why I'm getting 4 df
#Residual           

#the goal is this next model is to reflect what was done in Jackrel et al. 2019 -- WHICH I think is to include interactions terms
mod5.noH20<-dbrda(drdat~HFABac_meta_pruned$DECOMP_DAY+HFABac_meta_pruned$DEPLOYMENT_SITE+HFABac_meta_pruned$LOCATION_CODEwDay0+HFABac_meta_pruned$DECOMP_DAY*HFABac_meta_pruned$DEPLOYMENT_SITE+HFABac_meta_pruned$DECOMP_DAY*HFABac_meta_pruned$LOCATION_CODEwDay0+HFABac_meta_pruned$DEPLOYMENT_SITE*HFABac_meta_pruned$LOCATION_CODEwDay0)
anova(mod5.noH20,permutations = h, by="margin")
#                                                                           Df SumOfSqs      F    Pr(>F) 
#HFABac_meta_pruned$DECOMP_DAY:HFABac_meta_pruned$DEPLOYMENT_SITE           9    2.791 0.9447    0.6927    
#HFABac_meta_pruned$DECOMP_DAY:HFABac_meta_pruned$LOCATION_CODEwDay0       12    4.715 1.1971    0.0179 *  
#HFABac_meta_pruned$DEPLOYMENT_SITE:HFABac_meta_pruned$LOCATION_CODEwDay0   4    2.442 1.8598 9.999e-05 ***
#Residual                                                                 355  116.525                     

mod6.noH20<-dbrda(drdat~HFABac_meta_pruned$DECOMP_DAY+HFABac_meta_pruned$DEPLOYMENT_SITE+HFABac_meta_pruned$LOCATION_CODEwDay0+HFABac_meta_pruned$DECOMP_DAY*HFABac_meta_pruned$DEPLOYMENT_SITE+HFABac_meta_pruned$DECOMP_DAY*HFABac_meta_pruned$LOCATION_CODEwDay0+HFABac_meta_pruned$DEPLOYMENT_SITE*HFABac_meta_pruned$LOCATION_CODEwDay0+Condition(HFABac_meta_pruned$TREE_NUMBER))
anova(mod6.noH20,permutations = h, by="margin")
#                                                                          Df SumOfSqs      F    Pr(>F) 
#HFABac_meta_pruned$DECOMP_DAY:HFABac_meta_pruned$DEPLOYMENT_SITE           9    2.791 0.9735 0.579542    
#HFABac_meta_pruned$DECOMP_DAY:HFABac_meta_pruned$LOCATION_CODEwDay0       12    4.715 1.2336 0.007499 ** 
#HFABac_meta_pruned$DEPLOYMENT_SITE:HFABac_meta_pruned$LOCATION_CODEwDay0   4    2.920 2.2918 0.000200 ***
#Residual                                                                 336  107.032                    

#write.csv(HFABac_meta_pruned,file="HFABac_meta_pruned.csv") #trying to break off the phyloseq formal class in case its interacting with my degrees of freedom.
HFABac_meta_pruned1<-read.csv(file="HFABac_meta_pruned.csv") #basic data frame, good. Make all the things factor again.
HFABac_meta_pruned1$ORIGIN_SITE<-as.factor(HFABac_meta_pruned1$ORIGIN_SITE)
HFABac_meta_pruned1$TREE_NUMBER<-as.factor(HFABac_meta_pruned1$TREE_NUMBER)
HFABac_meta_pruned1$DEPLOYMENT_SITE<-as.factor(HFABac_meta_pruned1$DEPLOYMENT_SITE)
HFABac_meta_pruned1$DECOMP_DAY<-as.factor(HFABac_meta_pruned1$DECOMP_DAY)
HFABac_meta_pruned1$HFA_2<-as.factor(HFABac_meta_pruned1$HFA_2)
HFABac_meta_pruned1$LOCATION_CODE<-as.factor(HFABac_meta_pruned1$LOCATION_CODE)
HFABac_meta_pruned1$HFA_CODE<-as.factor(HFABac_meta_pruned1$HFA_CODE)
HFABac_meta_pruned1$TERR_AQUA<-as.factor(HFABac_meta_pruned1$TERR_AQUA)
HFABac_meta_pruned1$LOCATION_CODEwDay0<-as.factor(HFABac_meta_pruned1$LOCATION_CODEwDay0)
HFABac_meta_pruned1$Deployed_River<-as.factor(HFABac_meta_pruned1$Deployed_River)

#Quant. Jaccard RDA with no H20 samples, all members of the bacterial community with decomposition day (0-20 at a 5 factor level), deployed river (either sekiu or hoko) and location code (home and away 5 point code w Day 0 as not deployed), while accounting for tree genotype as a random effect. This is about to annoy me, there is no home and away for fresh leaves. It would be nice to see them cluster independently though -- away from terrestrial leaves. 
this.model<-dbrda(drdat~HFABac_meta_pruned1$DECOMP_DAY+HFABac_meta_pruned1$Deployed_River+HFABac_meta_pruned1$LOCATION_CODEwDay0+Condition(HFABac_meta_pruned1$TREE_NUMBER))
anova(this.model,permutations=h, by="margin")
#                                        Df SumOfSqs      F    Pr(>F)
#HFABac_meta_pruned$DECOMP_DAY           3    4.169 4.2936 9.999e-05 ***
#HFABac_meta_pruned$Deployed_River       1    0.293 0.9052    0.5594    
#HFABac_meta_pruned$LOCATION_CODEwDay0   4    0.950 0.7340    0.9909  #still unsure why i'm being penalized by n-2 but OK. Rare species might have a driving effect here.


str(HFABac_meta_pruned)

#ORDINATION ANALYSES
bacterial.phylo.4analysis1<-phyloseq(otu_table(bac.tab.df, taxa_are_rows=FALSE), sample_data(HFABac_meta_pruned), tax_table(HFABac_taxa_pruned), phy_tree(HFABac_phylo_pruned)) 
str(bacterial.phylo.4analysis1@sam_data)
dim(bacterial.phylo.4analysis1@otu_table)
bacterial.phylo.4analysis1@sam_data$DECOMP_DAY[bacterial.phylo.4analysis1@sam_data$DECOMP_DAY=="0"]<-"Fresh Leaves"
bacterial.phylo.4analysis1@sam_data$DECOMP_DAY[bacterial.phylo.4analysis1@sam_data$DECOMP_DAY=="5"]<-"Day 5"
bacterial.phylo.4analysis1@sam_data$DECOMP_DAY[bacterial.phylo.4analysis1@sam_data$DECOMP_DAY=="10"]<-"Day 10"
bacterial.phylo.4analysis1@sam_data$DECOMP_DAY[bacterial.phylo.4analysis1@sam_data$DECOMP_DAY=="15"]<-"Day 15"
bacterial.phylo.4analysis1@sam_data$DECOMP_DAY[bacterial.phylo.4analysis1@sam_data$DECOMP_DAY=="20"]<-"Day 20"
bacterial.phylo.4analysis1@sam_data$DECOMP_DAY<-factor(bacterial.phylo.4analysis1@sam_data$DECOMP_DAY,levels=c("Fresh Leaves","Day 5","Day 10","Day 15","Day 20"))
bacterial.phylo.4analysis1@sam_data$Deployed_River_test<-factor(bacterial.phylo.4analysis1@sam_data$Deployed_River_test,levels=c("Hoko River","Sekiu River"))
bacterial.phylo.4analysis1@sam_data$ORIGIN_RIVER[bacterial.phylo.4analysis1@sam_data$ORIGIN_RIVER=="Hoko"]<-"Hoko River"
bacterial.phylo.4analysis1@sam_data$ORIGIN_RIVER[bacterial.phylo.4analysis1@sam_data$ORIGIN_RIVER=="Sekiu"]<-"Sekiu River"
bacterial.phylo.4analysis1@sam_data$ORIGIN_RIVER<-factor(bacterial.phylo.4analysis1@sam_data$ORIGIN_RIVER,levels=c("Hoko River","Sekiu River"))
bacterial.phylo.4analysis1@sam_data$TREE_NUMBER<-as.factor(bacterial.phylo.4analysis1@sam_data$TREE_NUMBER)
bacterial.phylo.4analysis1@sam_data$LOCATION_CODEwDay0<-as.factor(bacterial.phylo.4analysis1@sam_data$LOCATION_CODEwDay0)

#Unconstrained Ordinations #toggel between bacterial.phylo.4analysis and bacterial.phylo.4analysis1 for prop abund and non rarefied
hfa_bac_pcoa <- ordinate(physeq = bacterial.phylo.4analysis, method = "PCoA", distance = "jaccard", binary= FALSE) #includes all taxa, jaccard distance
wu_hfa_bac_pcoa <- ordinate(physeq = bacterial.phylo.4analysis, method = "PCoA", distance = "wunifrac") #includes all taxa, jaccard distance

#Plot 
plot_ordination(physeq = bacterial.phylo.4analysis, ordination = hfa_bac_pcoa, color = "DECOMP_DAY", shape = "Deployed_River_test", axes = c(1,2)) + 
  scale_color_manual(values = c("#117878","#55b5e8","#0371b0","#d95e00","#cb7ca8")) +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)

plot_ordination(physeq = bacterial.phylo.4analysis, ordination = wu_hfa_bac_pcoa, color = "DECOMP_DAY", shape = "Deployed_River_test", axes = c(1,2)) + 
  scale_color_manual(values = c("#117878","#55b5e8","#0371b0","#d95e00","#cb7ca8")) +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)

hfa_bac_pcoa1 <- ordinate(physeq = bacterial.phylo.4analysis, method = "PCoA", distance = "bray") #includes all taxa, bray distance

plot_ordination(physeq = bacterial.phylo.4analysis1, ordination = hfa_bac_pcoa1, color = "DECOMP_DAY", shape = "Deployed_River_test") + 
  scale_color_manual(values = c("#117878","#55b5e8","#0371b0","#d95e00","#cb7ca8")) +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)

#transform sample counts to relative abundance
bacterial.phylo.r<-transform_sample_counts(bacterial.phylo.4analysis1, function(x) x / sum(x))
dim(bacterial.phylo.r@otu_table) #392 x 27021
bacterial.phylo.ra<-filter_taxa(bacterial.phylo.r, function(x) sum(x) > 0.01, TRUE) #retain greater than ##% of the community using within sample relative abundance
dim(bacterial.phylo.ra@otu_table) #392 x 1853

hfa_bac_pcoa3 <- ordinate(physeq = bacterial.phylo.ra, method = "PCoA", distance = "jaccard", binary= FALSE) #includes all taxa, jaccard distance

#Plot 
plot_ordination(physeq = bacterial.phylo.ra, ordination = hfa_bac_pcoa3, color = "DECOMP_DAY", shape = "Deployed_River_test", axes = c(1,2)) + 
  scale_color_manual(values = c("#117878","#55b5e8","#0371b0","#d95e00","#cb7ca8")) +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)

#NMDS
set.seed(32)
HFA_Bac_nmds <- ordinate(physeq = bacterial.phylo.ra, method = "NMDS", distance = "jaccard", binary= FALSE)

plot_ordination(physeq = bacterial.phylo.ra, ordination = HFA_Bac_nmds, color = "DECOMP_DAY", shape = "Deployed_River_test") + 
  scale_color_manual(values = c("#117878","#55b5e8","#0371b0","#d95e00","#cb7ca8")) +
  scale_shape_manual(values = c(1,19))+
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 3)

#constrained ordination
jacc.dist <- phyloseq::distance(physeq = bacterial.phylo.ra, method = "jaccard", binary=FALSE)

# CAP ordinate
cap_ord <- ordinate(physeq = bacterial.phylo.ra, method = "CAP", distance = jacc.dist, formula = ~ DECOMP_DAY + Deployed_River_test + LOCATION_CODEwDay0 + Condition(TREE_NUMBER))

cap_ord1<- ordinate(physeq = bacterial.phylo.ra, method = "CAP", distance = jacc.dist, formula = ~ DECOMP_DAY * Deployed_River_test * LOCATION_CODEwDay0 + Condition(TREE_NUMBER))

plot_ordination(physeq = bacterial.phylo.ra, ordination = cap_ord, color = "DECOMP_DAY", shape = "Deployed_River_test", axes = c(1,2)) + 
  scale_color_manual(values = c("#117878","#55b5e8","#0371b0","#d95e00","#cb7ca8")) +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)

# RDA ordinate
rda_ord <- ordinate(physeq = bacterial.phylo.ra, method = "RDA", distance = jacc.dist, formula = ~ DECOMP_DAY + Deployed_River_test + LOCATION_CODEwDay0 + Condition(TREE_NUMBER))

rda_ord1 <- ordinate(physeq = bacterial.phylo.ra, method = "RDA", distance = jacc.dist, formula = ~ DECOMP_DAY * Deployed_River_test * LOCATION_CODEwDay0 + Condition(TREE_NUMBER))

plot_ordination(physeq = bacterial.phylo.ra, ordination = rda_ord1, color = "DECOMP_DAY", shape = "Deployed_River_test", axes = c(1,2)) + 
  scale_color_manual(values = c("#117878","#55b5e8","#0371b0","#d95e00","#cb7ca8")) +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 3, stroke=1)

myround <- function(x) { trunc(x + 0.5) }

scale_reads <- function(physeq, n = min(sample_sums(physeq)), round = "floor") {
  
  # transform counts to n
  physeq.scale <- transform_sample_counts(physeq, 
                                          function(x) {(n * x/sum(x))}
  )
  
  # Pick the rounding functions
  if (round == "floor"){
    otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
  } else if (round == "round"){
    otu_table(physeq.scale) <- myround(otu_table(physeq.scale))
  }
  
  # Prune taxa and return new phyloseq object
  physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
  return(physeq.scale)
}

hfa_scale <- bacterial.phylo.4analysis1 %>%
  scale_reads(round = "round") 

hfa_scale_pcoa <- ordinate(physeq = hfa_scale, method = "PCoA", distance = "jaccard", binary= FALSE)

plot_ordination(physeq = hfa_scale, ordination = hfa_scale_pcoa, color = "DECOMP_DAY", shape = "Deployed_River_test", axes = c(1,2)) + 
  scale_color_manual(values = c("#117878","#55b5e8","#0371b0","#d95e00","#cb7ca8")) +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)

jacc.dist.scale <- phyloseq::distance(physeq = hfa_scale, method = "jaccard", binary=FALSE)

# CAP ordinate
cap_ord_scale <- ordinate(physeq = hfa_scale, method = "CAP", distance = jacc.dist.scale, formula = ~ DECOMP_DAY + Deployed_River_test + LOCATION_CODEwDay0 + Condition(TREE_NUMBER))
plot_ordination(physeq = hfa_scale, ordination = cap_ord_scale, color = "DECOMP_DAY", shape = "Deployed_River_test", axes = c(1,2)) + 
  scale_color_manual(values = c("#117878","#55b5e8","#0371b0","#d95e00","#cb7ca8")) +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)

#save(list=ls(),file="bacterialHFAordinations.rda")

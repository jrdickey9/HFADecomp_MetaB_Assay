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
#save(list=ls(),file="bacterialHFAanalyses5kRarefiedTESTMay30.rda")
#setwd("~/Desktop/Jonathan_Dickey_LabMac/HFA/Analyses/Bacteria")
#load(file="bacterialHFAanalyses.rda")

#Read in biom table from working directory
setwd("~/Desktop/Jonathan_Dickey_LabMac/HFA/bacterial_prelim/HFABac_paired_dadatableR1")
ASV_reads<-read_biom("HFABac-feature-table-FINAL.biom")
ASV_table<-as.data.frame(as.matrix(biom_data(ASV_reads)))
ASV_table[1:10,1:10]
otu_tab<-t(ASV_table)
dim(otu_tab) #428 samples by 30,989 
otu_tab[1:10, 1:10]
rownames(otu_tab) #ordered with 16S_JRD10 as first entry

#Read in metadata file
setwd("/Users/lab/Desktop/Jonathan_Dickey_LabMac/HFA/Analyses/Bacteria")
meta_data<-read.csv("updated_meta_bacterial_data.csv",header=TRUE) #428 observations x 17 variables
head(meta_data)
str(meta_data)
meta_data<-meta_data[order(meta_data$SHORT_SAMPLE_NAME),]
dim(meta_data) #428 x 18

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
rownames(bac.taxa.m)
str(otu_tab) #428 x 30,989
df.bacterial.OTU<-as.data.frame(otu_tab)
str(meta_data) #data.frame 428 x 18

#Matching row names
rownames(df.bacterial.OTU)<-as.character(meta_data[,4])
colnames(df.bacterial.OTU) #accession numbers
rownames(meta_data)<-as.character(meta_data[,4]) #sample names
rownames(bac.taxa.m)<-as.character(rownames(bac.taxa.m)) #these the accession numbers
samp.names<-as.character(meta_data[,4]) #for example, "HOK2_4_HOK1_0_JRD"

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
dim(HFABac_phylo@sam_data) #392 x 18
HFABac_meta_pruned<-HFABac_phylo@sam_data

dim(HFABac_phylo@otu_table) #27,596

#ASV table with out zeros and singletons
HFABac_phylo.rm1<- prune_taxa(taxa_sums(HFABac_phylo) > 0, HFABac_phylo) #remove ASVs with zero abundance attributed to KatharoSeq and water columns
dim(HFABac_phylo.rm1@otu_table) #392 x 22,706
HFABac_phylo.rm2<- prune_taxa(taxa_sums(HFABac_phylo.rm1) > 1, HFABac_phylo.rm1) #remove singletons
dim(HFABac_phylo.rm2@otu_table) #392 x 22,454

HFABac_phylo.ASVtable.2<-HFABac_phylo.rm2@otu_table

#subset samples by decomposition day
# Subset5_20<-subset_samples(HFABac_phylo.rm2, DECOMP_DAY !="0")
# Subset5_20@sam_data$DECOMP_DAY
# min(rowSums(Subset5_20@otu_table))

#Rarefying (haulted due to low read depth)
min(rowSums(HFABac_phylo.ASVtable.2)) #480 Alright, so there was lots of host DNA in some of these samples. 
sort(x=rowSums(HFABac_phylo.ASVtable.2),decreasing=FALSE) #Consider samples to toss. It's a good amount below 10K. 

#remove low read depth samples to rarefy to 5K. #I NEED TO REVIIST THIS and COMPARE TO THE NON RAREFIED SCRIPT. # a difference of 12 samples when i made things match before this line. Also is there a better way to do this?
HFABac_phylo_5K<-subset_samples(HFABac_phylo.rm2, SHORT_SAMPLE_NAME != "16S_JRD8" & SHORT_SAMPLE_NAME != "16S_JRD9" & SHORT_SAMPLE_NAME != "16S_JRD396" & SHORT_SAMPLE_NAME != "16S_JRD98" & SHORT_SAMPLE_NAME != "16S_JRD4" & SHORT_SAMPLE_NAME != "16S_JRD13" & SHORT_SAMPLE_NAME != "16S_JRD12" & SHORT_SAMPLE_NAME != "16S_JRD15" & SHORT_SAMPLE_NAME != "JRD362" & SHORT_SAMPLE_NAME != "16S_JRD10" & SHORT_SAMPLE_NAME != "16S_JRD11" & SHORT_SAMPLE_NAME != "JRD6" & SHORT_SAMPLE_NAME != "JRD384" & SHORT_SAMPLE_NAME != "JRD377" & SHORT_SAMPLE_NAME != "JRD14" & SHORT_SAMPLE_NAME != "JRD374" & SHORT_SAMPLE_NAME != "JRD3" & SHORT_SAMPLE_NAME != "JRD381"  & SHORT_SAMPLE_NAME != "JRD358" & SHORT_SAMPLE_NAME != "JRD382" & SHORT_SAMPLE_NAME != "JRD92" & SHORT_SAMPLE_NAME != "JRD351" & SHORT_SAMPLE_NAME != "JRD378" & SHORT_SAMPLE_NAME != "JRD350" & SHORT_SAMPLE_NAME != "JRD361" & SHORT_SAMPLE_NAME != "JRD85" & SHORT_SAMPLE_NAME != "JRD363" & SHORT_SAMPLE_NAME != "JRD220" & SHORT_SAMPLE_NAME != "JRD375" & SHORT_SAMPLE_NAME != "JRD349" & SHORT_SAMPLE_NAME != "JRD383" & SHORT_SAMPLE_NAME != "JRD16" & SHORT_SAMPLE_NAME != "JRD379" & SHORT_SAMPLE_NAME != "JRD18" & SHORT_SAMPLE_NAME != "JRD265" & SHORT_SAMPLE_NAME != "JRD376" & SHORT_SAMPLE_NAME != "JRD94" & SHORT_SAMPLE_NAME != "JRD352" & SHORT_SAMPLE_NAME != "JRD1" & SHORT_SAMPLE_NAME != "JRD5" & SHORT_SAMPLE_NAME != "JRD87" & SHORT_SAMPLE_NAME != "JRD391" & SHORT_SAMPLE_NAME != "JRD95" & SHORT_SAMPLE_NAME != "JRD373" & SHORT_SAMPLE_NAME != "JRD355" & SHORT_SAMPLE_NAME != "JRD394" & SHORT_SAMPLE_NAME != "JRD395" & SHORT_SAMPLE_NAME != "JRD338" & SHORT_SAMPLE_NAME != "JRD366" & SHORT_SAMPLE_NAME != "JRD86" & SHORT_SAMPLE_NAME != "JRD19" & SHORT_SAMPLE_NAME != "JRD2" & SHORT_SAMPLE_NAME != "JRD254" & SHORT_SAMPLE_NAME != "JRD17" & SHORT_SAMPLE_NAME != "JRD380" & SHORT_SAMPLE_NAME != "JRD83" & SHORT_SAMPLE_NAME != "JRD392" & SHORT_SAMPLE_NAME != "JRD385" & SHORT_SAMPLE_NAME != "JRD353" & SHORT_SAMPLE_NAME != "JRD96" & SHORT_SAMPLE_NAME != "JRD84" & SHORT_SAMPLE_NAME != "JRD354" & SHORT_SAMPLE_NAME != "JRD356" & SHORT_SAMPLE_NAME != "JRD357" & SHORT_SAMPLE_NAME != "JRD91" & SHORT_SAMPLE_NAME != "JRD393")


HFABac_phylo.5K.rm1<- prune_taxa(taxa_sums(HFABac_phylo_5K) > 0, HFABac_phylo_5K) #remove zeros 
dim(HFABac_phylo.5K.rm1@otu_table) #326 x 22,328 

HFABac_phylo_5K_ASVtab<-HFABac_phylo.5K.rm1@otu_table
min(rowSums(HFABac_phylo_5K_ASVtab))
bac.tab.df<-as.data.frame(HFABac_phylo_5K_ASVtab)
rdat<-rrarefy(bac.tab.df,5222) #rarefy!

#Standardize abundances into proportions. 
rowSums(rdat) #Check if it worked
std.bac.tab<-decostand(rdat,"total") #replace object to rdat after rarefying. 
std.bac.tab[1:10,1:10] #Looks groovy! 

#Smush back together into single phyloseq object; this contains rarefied and standardized data
bacterial.phylo.4analysis<-phyloseq(otu_table(std.bac.tab, taxa_are_rows=FALSE), sample_data(HFABac_phylo.5K.rm1@sam_data), tax_table(HFABac_phylo.5K.rm1@tax_table), phy_tree(HFABac_phylo.5K.rm1@phy_tree))

#Create Quant Jaccard distance matrix
drdat<-vegdist(std.bac.tab,"jaccard")

HFABac_meta_pruned5k<-HFABac_phylo.5K.rm1@sam_data
#Revisiting meta data to build db-rda in an appropriate way
str(HFABac_meta_pruned5k) #326 x 18
colnames(HFABac_meta_pruned5k)
HFABac_meta_pruned5k$ORIGIN_SITE<-as.factor(HFABac_meta_pruned5k$ORIGIN_SITE)
HFABac_meta_pruned5k$TREE_NUMBER<-as.factor(HFABac_meta_pruned5k$TREE_NUMBER)
HFABac_meta_pruned5k$DEPLOYMENT_SITE<-as.factor(HFABac_meta_pruned5k$DEPLOYMENT_SITE)
HFABac_meta_pruned5k$DECOMP_DAY<-as.factor(HFABac_meta_pruned5k$DECOMP_DAY)
HFABac_meta_pruned5k$HFA_2<-as.factor(HFABac_meta_pruned5k$HFA_2)
HFABac_meta_pruned5k$LOCATION_CODE<-as.factor(HFABac_meta_pruned5k$LOCATION_CODE)
HFABac_meta_pruned5k$HFA_CODE<-as.factor(HFABac_meta_pruned5k$HFA_CODE)
HFABac_meta_pruned5k$TERR_AQUA<-as.factor(HFABac_meta_pruned5k$TERR_AQUA)
HFABac_meta_pruned5k$LOCATION_CODEwDay0<-as.factor(HFABac_meta_pruned5k$LOCATION_CODEwDay0)
HFABac_meta_pruned5k$Deployed_River<-as.factor(HFABac_meta_pruned5k$Deployed_River)
HFABac_meta_pruned5k$Argonne_Run<-as.factor(HFABac_meta_pruned5k$Argonne_Run)

basic.mod<-dbrda(drdat~HFABac_meta_pruned5k$ORIGIN_SITE+HFABac_meta_pruned5k$DECOMP_DAY+HFABac_meta_pruned5k$HFA_2)
h<-how(nperm=10000)
#anova(basic.mod, permutations = h, by="margin") 
#                                   Df SumOfSqs      F    Pr(>F)
#HFABac_meta_pruned5k$ORIGIN_SITE   3    3.054 3.1740 9.999e-05 ***
#HFABac_meta_pruned5k$DECOMP_DAY    3    3.461 3.5970 9.999e-05 ***
#HFABac_meta_pruned5k$HFA_2         1    0.369 1.1499    0.2248    
#Residual                         329  105.521     

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
mod2.noH2O<-dbrda(drdat~HFABac_meta_pruned5k$DECOMP_DAY+HFABac_meta_pruned5k$LOCATION_CODE+Condition(HFABac_meta_pruned5k$TREE_NUMBER))
#anova(mod2.noH2O,permutations = h, by="margin")
#                                    Df SumOfSqs      F    Pr(>F)    
#HFABac_meta_pruned5k$DECOMP_DAY      4    4.366 3.4024 9.999e-05 ***
#HFABac_meta_pruned5k$LOCATION_CODE   4    2.296 1.7896 9.999e-05 ***
#Residual                           310   99.437

#I think from what I can gleam that an important model to run given Jackrel et al. 2019 "The origin, succession..." would be to investigate the fixed effects of decomposition day and deployment site. 
mod3.noH2O<-dbrda(drdat~HFABac_meta_pruned5k$DECOMP_DAY+HFABac_meta_pruned5k$DEPLOYMENT_SITE+Condition(HFABac_meta_pruned5k$TREE_NUMBER))
#anova(mod3.noH2O,permutations = h, by="margin")
#                                      Df SumOfSqs      F    Pr(>F)  
#HFABac_meta_pruned5k$DECOMP_DAY        4    3.992 3.2454 9.999e-05 ***
#HFABac_meta_pruned5k$DEPLOYMENT_SITE   3    6.087 6.5974 9.999e-05 ***
#Residual                             311   95.646  

mod4.b.noH2O<-dbrda(drdat~HFABac_meta_pruned5k$DECOMP_DAY+HFABac_meta_pruned5k$DEPLOYMENT_SITE+HFABac_meta_pruned5k$LOCATION_CODEwDay0+Condition(HFABac_meta_pruned5k$TREE_NUMBER))
#anova(mod4.b.noH2O,permutations = h, by="margin")
#                                        Df SumOfSqs      F    Pr(>F)
#HFABac_meta_pruned5k$DECOMP_DAY           3    3.397 3.7046 9.999e-05 ***
#HFABac_meta_pruned5k$DEPLOYMENT_SITE      3    5.430 5.9217 9.999e-05 ***
#HFABac_meta_pruned5k$LOCATION_CODEwDay0   4    1.802 1.4739    0.0037 ** 
#Residual                                307   93.844           

#the goal is this next model is to reflect what was done in Jackrel et al. 2019 -- WHICH I think is to include interactions terms
mod5.noH20<-dbrda(drdat~HFABac_meta_pruned5k$DECOMP_DAY+HFABac_meta_pruned5k$DEPLOYMENT_SITE+HFABac_meta_pruned5k$LOCATION_CODEwDay0+HFABac_meta_pruned5k$DECOMP_DAY*HFABac_meta_pruned5k$DEPLOYMENT_SITE+HFABac_meta_pruned5k$DECOMP_DAY*HFABac_meta_pruned5k$LOCATION_CODEwDay0+HFABac_meta_pruned5k$DEPLOYMENT_SITE*HFABac_meta_pruned5k$LOCATION_CODEwDay0)
#anova(mod5.noH20,permutations = h, by="margin")
#                                                                               Df SumOfSqs      F    Pr(>F) 
#HFABac_meta_pruned5k$DECOMP_DAY:HFABac_meta_pruned5k$DEPLOYMENT_SITE           9    4.205 1.5897 9.999e-05 ***
#HFABac_meta_pruned5k$DECOMP_DAY:HFABac_meta_pruned5k$LOCATION_CODEwDay0       12    4.404 1.2488  0.005399 ** 
#HFABac_meta_pruned5k$DEPLOYMENT_SITE:HFABac_meta_pruned5k$LOCATION_CODEwDay0   4    2.641 2.2461 9.999e-05 ***
#Residual                  

mod6.noH20<-dbrda(drdat~HFABac_meta_pruned5k$DECOMP_DAY+HFABac_meta_pruned5k$DEPLOYMENT_SITE+HFABac_meta_pruned5k$LOCATION_CODEwDay0+HFABac_meta_pruned5k$DECOMP_DAY*HFABac_meta_pruned5k$DEPLOYMENT_SITE+HFABac_meta_pruned5k$DECOMP_DAY*HFABac_meta_pruned5k$LOCATION_CODEwDay0+HFABac_meta_pruned5k$DEPLOYMENT_SITE*HFABac_meta_pruned5k$LOCATION_CODEwDay0+Condition(HFABac_meta_pruned5k$TREE_NUMBER))
#anova(mod6.noH20,permutations = h, by="margin")
#                                                                               Df SumOfSqs      F    Pr(>F) 
#HFABac_meta_pruned5k$DECOMP_DAY:HFABac_meta_pruned5k$DEPLOYMENT_SITE           9    4.281 1.6302 9.999e-05 ***
#HFABac_meta_pruned5k$DECOMP_DAY:HFABac_meta_pruned5k$LOCATION_CODEwDay0       12    4.377 1.2500    0.0034 ** 
#HFABac_meta_pruned5k$DEPLOYMENT_SITE:HFABac_meta_pruned5k$LOCATION_CODEwDay0   4    2.406 2.0614 9.999e-05 ***
#Residual                                                                     279   81.408                    

#write.csv(HFABac_meta_pruned,file="HFABac_meta_pruned.csv") #trying to break off the phyloseq formal class in case its interacting with my degrees of freedom.
# HFABac_meta_pruned1<-read.csv(file="HFABac_meta_pruned.csv") #basic data frame, good. Make all the things factor again.
# HFABac_meta_pruned1$ORIGIN_SITE<-as.factor(HFABac_meta_pruned1$ORIGIN_SITE)
# HFABac_meta_pruned1$TREE_NUMBER<-as.factor(HFABac_meta_pruned1$TREE_NUMBER)
# HFABac_meta_pruned1$DEPLOYMENT_SITE<-as.factor(HFABac_meta_pruned1$DEPLOYMENT_SITE)
# HFABac_meta_pruned1$DECOMP_DAY<-as.factor(HFABac_meta_pruned1$DECOMP_DAY)
# HFABac_meta_pruned1$HFA_2<-as.factor(HFABac_meta_pruned1$HFA_2)
# HFABac_meta_pruned1$LOCATION_CODE<-as.factor(HFABac_meta_pruned1$LOCATION_CODE)
# HFABac_meta_pruned1$HFA_CODE<-as.factor(HFABac_meta_pruned1$HFA_CODE)
# HFABac_meta_pruned1$TERR_AQUA<-as.factor(HFABac_meta_pruned1$TERR_AQUA)
# HFABac_meta_pruned1$LOCATION_CODEwDay0<-as.factor(HFABac_meta_pruned1$LOCATION_CODEwDay0)
# HFABac_meta_pruned1$Deployed_River<-as.factor(HFABac_meta_pruned1$Deployed_River)

#Quant. Jaccard RDA with no H20 samples, all members of the bacterial community with decomposition day (0-20 at a 5 factor level), deployed river (either sekiu or hoko) and location code (home and away 5 point code w Day 0 as not deployed), while accounting for tree genotype as a random effect. This is about to annoy me, there is no home and away for fresh leaves. It would be nice to see them cluster independently though -- away from terrestrial leaves. 
this.model<-dbrda(drdat~HFABac_meta_pruned5k$DECOMP_DAY+HFABac_meta_pruned5k$Deployed_River+HFABac_meta_pruned5k$LOCATION_CODEwDay0+Condition(HFABac_meta_pruned5k$TREE_NUMBER))
#anova(this.model,permutations=h, by="margin")
#                                        Df SumOfSqs      F    Pr(>F)
#HFABac_meta_pruned5k$DECOMP_DAY           3    3.508 3.7072 9.999e-05 ***
#HFABac_meta_pruned5k$Deployed_River       1    1.805 5.7207 9.999e-05 ***
#HFABac_meta_pruned5k$LOCATION_CODEwDay0   4    2.514 1.9928 9.999e-05 ***
#Residual                                309   97.470

this.model1<-dbrda(drdat~HFABac_meta_pruned5k$DECOMP_DAY+HFABac_meta_pruned5k$Deployed_River+HFABac_meta_pruned5k$LOCATION_CODEwDay0+HFABac_meta_pruned5k$DECOMP_DAY*HFABac_meta_pruned5k$Deployed_River+HFABac_meta_pruned5k$DECOMP_DAY*HFABac_meta_pruned5k$LOCATION_CODEwDay0+HFABac_meta_pruned5k$Deployed_River*HFABac_meta_pruned5k$LOCATION_CODEwDay0+Condition(HFABac_meta_pruned5k$TREE_NUMBER))

anova(basic.mod, permutations = h, by="margin") 
anova(mod2.noH2O,permutations = h, by="margin")
anova(mod3.noH2O,permutations = h, by="margin")
anova(mod4.b.noH2O,permutations = h, by="margin")
anova(mod5.noH20,permutations = h, by="margin")
anova(mod6.noH20,permutations = h, by="margin")
anova(this.model,permutations=h, by="margin")
anova(this.model1,permutations=h, by="margin")

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#ORDINATION ANALYSES; don't need to run this phylo object, just here for continuity
#bacterial.phylo.4analysis<-phyloseq(otu_table(std.bac.tab, taxa_are_rows=FALSE), sample_data(HFABac_phylo.5K.rm1@sam_data), tax_table(HFABac_phylo.5K.rm1@tax_table), phy_tree(HFABac_phylo.5K.rm1@phy_tree))
str(bacterial.phylo.4analysis)
bacterial.phylo.4analysis@otu_table[1:10,1:10]
bacterial.phylo.4analysis@sam_data$DECOMP_DAY[bacterial.phylo.4analysis@sam_data$DECOMP_DAY=="0"]<-"Fresh Leaves"
bacterial.phylo.4analysis@sam_data$DECOMP_DAY[bacterial.phylo.4analysis@sam_data$DECOMP_DAY=="5"]<-"Day 5"
bacterial.phylo.4analysis@sam_data$DECOMP_DAY[bacterial.phylo.4analysis@sam_data$DECOMP_DAY=="10"]<-"Day 10"
bacterial.phylo.4analysis@sam_data$DECOMP_DAY[bacterial.phylo.4analysis@sam_data$DECOMP_DAY=="15"]<-"Day 15"
bacterial.phylo.4analysis@sam_data$DECOMP_DAY[bacterial.phylo.4analysis@sam_data$DECOMP_DAY=="20"]<-"Day 20"
bacterial.phylo.4analysis@sam_data$DECOMP_DAY<-factor(bacterial.phylo.4analysis@sam_data$DECOMP_DAY,levels=c("Fresh Leaves","Day 5","Day 10","Day 15","Day 20"))
bacterial.phylo.4analysis@sam_data$Deployed_River_test<-factor(bacterial.phylo.4analysis@sam_data$Deployed_River_test,levels=c("Hoko River","Sekiu River"))
bacterial.phylo.4analysis@sam_data$ORIGIN_RIVER[bacterial.phylo.4analysis@sam_data$ORIGIN_RIVER=="Hoko"]<-"Hoko River"
bacterial.phylo.4analysis@sam_data$ORIGIN_RIVER[bacterial.phylo.4analysis@sam_data$ORIGIN_RIVER=="Sekiu"]<-"Sekiu River"
bacterial.phylo.4analysis@sam_data$ORIGIN_RIVER<-factor(bacterial.phylo.4analysis@sam_data$ORIGIN_RIVER,levels=c("Hoko River","Sekiu River"))
bacterial.phylo.4analysis@sam_data$TREE_NUMBER<-as.factor(bacterial.phylo.4analysis@sam_data$TREE_NUMBER)
bacterial.phylo.4analysis@sam_data$LOCATION_CODEwDay0<-as.factor(bacterial.phylo.4analysis@sam_data$LOCATION_CODEwDay0)

#Unconstrained Ordinations 
hfa_bac_pcoa <- ordinate(physeq = bacterial.phylo.4analysis, method = "PCoA", distance = "bray") 
wu_hfa_bac_pcoa <- ordinate(physeq = bacterial.phylo.4analysis, method = "PCoA", distance = "wunifrac") 
uwu_hfa_bac_pcoa <- ordinate(physeq = bacterial.phylo.4analysis, method = "PCoA", distance = "unifrac")

#Plot bray curtis pcoa 
plot_ordination(physeq = bacterial.phylo.4analysis, ordination = hfa_bac_pcoa, color = "DECOMP_DAY", shape = "Deployed_River_test", axes = c(1,2)) + 
  scale_color_manual(values = c("#016403","#35CEF2","#03045e","#ffb703","#f03b20")) +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="River",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right") 

#plot weighted unifrac pcoa
plot_ordination(physeq = bacterial.phylo.4analysis, ordination = wu_hfa_bac_pcoa, color = "DECOMP_DAY", shape = "Deployed_River_test", axes = c(1,2)) + 
  scale_color_manual(values = c("#016403","#35CEF2","#03045e","#ffb703","#f03b20")) +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="River",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right") 

#plot unweighted unifrac pcoa
plot_ordination(physeq = bacterial.phylo.4analysis, ordination = uwu_hfa_bac_pcoa, color = "DECOMP_DAY", shape = "Deployed_River_test", axes = c(2,3)) + 
  scale_color_manual(values = c("#016403","#35CEF2","#03045e","#ffb703","#f03b20")) +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="River",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right") 

#NMDS
set.seed(32)
HFA_Bac_nmds <- ordinate(physeq = bacterial.phylo.4analysis, method = "NMDS", distance = "bray")

plot_ordination(physeq = bacterial.phylo.4analysis, ordination = HFA_Bac_nmds, color = "DECOMP_DAY", shape = "Deployed_River_test") + 
  scale_color_manual(values = c("#016403","#35CEF2","#03045e","#ffb703","#f03b20")) +
  scale_shape_manual(values = c(1,19))+
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4)+
  theme_test()+
  labs(shape="River",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right") 

#constrained ordination
bray.dist <- phyloseq::distance(physeq = bacterial.phylo.4analysis, method = "bray")

# CAP ordinate
cap_ord <- ordinate(physeq = bacterial.phylo.4analysis, method = "CAP", distance = bray.dist, formula = ~ DECOMP_DAY + Deployed_River_test + LOCATION_CODEwDay0 + Condition(TREE_NUMBER))
cap_ord1<- ordinate(physeq = bacterial.phylo.4analysis, method = "CAP", distance = bray.dist, formula = ~ DECOMP_DAY * Deployed_River_test * LOCATION_CODEwDay0 + Condition(TREE_NUMBER))

plot_ordination(physeq = bacterial.phylo.4analysis, ordination = cap_ord, color = "DECOMP_DAY", shape = "Deployed_River_test", axes = c(1,2)) + 
  scale_color_manual(values = c("#016403","#35CEF2","#03045e","#ffb703","#f03b20")) +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="River",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right") 

plot_ordination(physeq = bacterial.phylo.4analysis, ordination = cap_ord1, color = "DECOMP_DAY", shape = "Deployed_River_test", axes = c(1,2)) + 
  scale_color_manual(values = c("#016403","#35CEF2","#03045e","#ffb703","#f03b20")) +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="River",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right") 

# RDA ordinate
rda_ord <- ordinate(physeq = bacterial.phylo.4analysis, method = "RDA", distance = bray.dist, formula = ~ DECOMP_DAY + Deployed_River_test + LOCATION_CODEwDay0 + Condition(TREE_NUMBER))
rda_ord1 <- ordinate(physeq = bacterial.phylo.4analysis, method = "RDA", distance = bray.dist, formula = ~ DECOMP_DAY * Deployed_River_test * LOCATION_CODEwDay0 + Condition(TREE_NUMBER))

plot_ordination(physeq = bacterial.phylo.4analysis, ordination = rda_ord, color = "DECOMP_DAY", shape = "Deployed_River_test", axes = c(1,2)) + 
  scale_color_manual(values = c("#016403","#35CEF2","#03045e","#ffb703","#f03b20")) +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="River",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right") 

plot_ordination(physeq = bacterial.phylo.4analysis, ordination = rda_ord1, color = "DECOMP_DAY", shape = "Deployed_River_test", axes = c(1,2)) + 
  scale_color_manual(values = c("#016403","#35CEF2","#03045e","#ffb703","#f03b20")) +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="River",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right") 

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#Days 5 to 20 Only
Subset5_20<-subset_samples(bacterial.phylo.4analysis, DECOMP_DAY !="Fresh Leaves")
Subset5_20@sam_data$DECOMP_DAY
dim(Subset5_20@otu_table) #308 x 22,328
Subset5_20.5K.rm1<- prune_taxa(taxa_sums(Subset5_20) > 0, Subset5_20) #remove zeros 
dim(Subset5_20.5K.rm1@otu_table) #308 x 9943

sunset0<-subset_samples(bacterial.phylo.4analysis, DECOMP_DAY !="Day 5" & DECOMP_DAY !="Day 10" & DECOMP_DAY !="Day 15" & DECOMP_DAY !="Day 20")
dim(sunset0@otu_table)
sunset0.rm1<- prune_taxa(taxa_sums(sunset0) > 0, sunset0) #remove zeros 
dim(sunset0.rm1@otu_table)
d0_all <- ordinate(physeq = sunset0.rm1, method = "PCoA", distance = "bray") 
plot_ordination(physeq = sunset0.rm1, ordination = d0_all, color = "DECOMP_DAY", shape = "ORIGIN_SITE", axes = c(2,3)) + 
  scale_color_manual(values = c("#35CEF2")) +
  scale_shape_manual(values = c(0, 15, 1, 19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="Origin Site",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right") 


#unconstrained ordination
hfa_bac_pcoa <- ordinate(physeq = Subset5_20.5K.rm1, method = "PCoA", distance = "bray") 
plot_ordination(physeq = Subset5_20.5K.rm1, ordination = hfa_bac_pcoa, color = "DECOMP_DAY", shape = "Deployed_River_test", axes = c(1,2)) + 
  scale_color_manual(values = c("#35CEF2","#03045e","#ffb703","#f03b20")) +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="River",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right") 

Hoko.River<-subset_samples(Subset5_20.5K.rm1, Deployed_River_test !="Sekiu River")
Sekiu.River<-subset_samples(Subset5_20.5K.rm1, Deployed_River_test !="Hoko River")
Hoko.River<- prune_taxa(taxa_sums(Hoko.River) > 0, Hoko.River)
Sekiu.River<- prune_taxa(taxa_sums(Sekiu.River) > 0, Sekiu.River)

Hoko.River@sam_data$DEPLOYMENT_SITE[Hoko.River@sam_data$DEPLOYMENT_SITE=="HOK1"]<-"Downstream Hoko"
Hoko.River@sam_data$DEPLOYMENT_SITE[Hoko.River@sam_data$DEPLOYMENT_SITE=="HOK2"]<-"Upstream Hoko"
Hoko.River@sam_data$DEPLOYMENT_SITE<-factor(Hoko.River@sam_data$DEPLOYMENT_SITE,levels=c("Upstream Hoko","Downstream Hoko"))
Hoko.River@sam_data$ORIGIN_SITE[Hoko.River@sam_data$ORIGIN_SITE=="HOK1"]<-"Downstream Hoko"
Hoko.River@sam_data$ORIGIN_SITE[Hoko.River@sam_data$ORIGIN_SITE=="HOK2"]<-"Upstream Hoko"
Hoko.River@sam_data$ORIGIN_SITE[Hoko.River@sam_data$ORIGIN_SITE=="SEK1"]<-"Downstream Sekiu"
Hoko.River@sam_data$ORIGIN_SITE[Hoko.River@sam_data$ORIGIN_SITE=="SEK2"]<-"Upstream Sekiu"
Hoko.River@sam_data$ORIGIN_SITE<-factor(Hoko.River@sam_data$ORIGIN_SITE,levels=c("Upstream Hoko","Downstream Hoko","Upstream Sekiu","Downstream Sekiu"))

Sekiu.River@sam_data$DEPLOYMENT_SITE[Sekiu.River@sam_data$DEPLOYMENT_SITE=="SEK1"]<-"Downstream Sekiu"
Sekiu.River@sam_data$DEPLOYMENT_SITE[Sekiu.River@sam_data$DEPLOYMENT_SITE=="SEK2"]<-"Upstream Sekiu"
Sekiu.River@sam_data$DEPLOYMENT_SITE<-factor(Sekiu.River@sam_data$DEPLOYMENT_SITE,levels=c("Upstream Sekiu","Downstream Sekiu"))

hfa_bac_pcoa_hoko <- ordinate(physeq = Hoko.River, method = "PCoA", distance = "bray") 
hfa_bac_pcoa_sekiu <- ordinate(physeq = Sekiu.River, method = "PCoA", distance = "bray") 

Hoko.River@sam_data[1:10,]
plot_ordination(physeq = Hoko.River, ordination = hfa_bac_pcoa_hoko, color = "DECOMP_DAY", shape = "ORIGIN_SITE", axes = c(1,2)) + 
  scale_color_manual(values = c("#35CEF2","#03045e","#ffb703","#f03b20")) +
  scale_shape_manual(values = c(0,15,1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="Origin Site",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right") 

plot_ordination(physeq = Sekiu.River, ordination = hfa_bac_pcoa_sekiu, color = "DECOMP_DAY", shape = "DEPLOYMENT_SITE", axes = c(1,2)) + 
  scale_color_manual(values = c("#35CEF2","#03045e","#ffb703","#f03b20")) +
  scale_shape_manual(values = c(0,15)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="Deployment Site",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right") 
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

#Hoko River by Decomposition Day; use bacterial.phylo.4analysis to see fresh leaves too
Hoko.River1<-subset_samples(bacterial.phylo.4analysis, Deployed_River_test !="Sekiu River")
Sekiu.River1<-subset_samples(bacterial.phylo.4analysis, Deployed_River_test !="Hoko River")
Hoko.River1<- prune_taxa(taxa_sums(Hoko.River1) > 0, Hoko.River1)
Sekiu.River1<- prune_taxa(taxa_sums(Sekiu.River1) > 0, Sekiu.River1)
Hoko.River1@sam_data$DEPLOYMENT_SITE[Hoko.River1@sam_data$DEPLOYMENT_SITE=="HOK1"]<-"Downstream Hoko"
Hoko.River1@sam_data$DEPLOYMENT_SITE[Hoko.River1@sam_data$DEPLOYMENT_SITE=="HOK2"]<-"Upstream Hoko"
Hoko.River1@sam_data$DEPLOYMENT_SITE<-factor(Hoko.River1@sam_data$DEPLOYMENT_SITE,levels=c("Upstream Hoko","Downstream Hoko"))
Sekiu.River1@sam_data$DEPLOYMENT_SITE[Sekiu.River1@sam_data$DEPLOYMENT_SITE=="SEK1"]<-"Downstream Sekiu"
Sekiu.River1@sam_data$DEPLOYMENT_SITE[Sekiu.River1@sam_data$DEPLOYMENT_SITE=="SEK2"]<-"Upstream Sekiu"
Sekiu.River1@sam_data$DEPLOYMENT_SITE<-factor(Sekiu.River1@sam_data$DEPLOYMENT_SITE,levels=c("Upstream Sekiu","Downstream Sekiu"))

HokoDayZERO<-subset_samples(Hoko.River1, DECOMP_DAY !="Day 5" & DECOMP_DAY !="Day 10" & DECOMP_DAY !="Day 15" & DECOMP_DAY !="Day 20")
HokoDayZERO<- prune_taxa(taxa_sums(HokoDayZERO) > 0, HokoDayZERO)
dim(HokoDayZERO@otu_table) #10 x 353
hoko_pcoa_d0 <- ordinate(physeq = HokoDayZERO, method = "PCoA", distance = "bray") 

HokoDayFIVE<-subset_samples(Hoko.River1, DECOMP_DAY !="Fresh Leaves" & DECOMP_DAY !="Day 10" & DECOMP_DAY !="Day 15" & DECOMP_DAY !="Day 20")
HokoDayFIVE<- prune_taxa(taxa_sums(HokoDayFIVE) > 0, HokoDayFIVE)
dim(HokoDayFIVE@otu_table) #35 x 1435
hoko_pcoa_d5 <- ordinate(physeq = HokoDayFIVE, method = "PCoA", distance = "bray") 

HokoDayTEN<-subset_samples(Hoko.River1, DECOMP_DAY !="Fresh Leaves" & DECOMP_DAY !="Day 5" & DECOMP_DAY !="Day 15" & DECOMP_DAY !="Day 20")
HokoDayTEN<- prune_taxa(taxa_sums(HokoDayTEN) > 0, HokoDayTEN)
dim(HokoDayTEN@otu_table) #39 x 2318
hoko_pcoa_d10 <- ordinate(physeq = HokoDayTEN, method = "PCoA", distance = "bray") 

HokoDay15<-subset_samples(Hoko.River1, DECOMP_DAY !="Fresh Leaves" & DECOMP_DAY !="Day 5" & DECOMP_DAY !="Day 10" & DECOMP_DAY !="Day 20")
HokoDay15<- prune_taxa(taxa_sums(HokoDay15) > 0, HokoDay15)
dim(HokoDay15@otu_table) #39 x 3699
hoko_pcoa_d15 <- ordinate(physeq = HokoDay15, method = "PCoA", distance = "bray") 

HokoDay20<-subset_samples(Hoko.River1, DECOMP_DAY !="Fresh Leaves" & DECOMP_DAY !="Day 5" & DECOMP_DAY !="Day 10" & DECOMP_DAY !="Day 15")
HokoDay20<- prune_taxa(taxa_sums(HokoDay20) > 0, HokoDay20)
dim(HokoDay20@otu_table) #39 x 4462
hoko_pcoa_d20 <- ordinate(physeq = HokoDay20, method = "PCoA", distance = "bray") 

h0<-plot_ordination(physeq = HokoDayZERO, ordination = hoko_pcoa_d0, color = "DECOMP_DAY", shape = "DEPLOYMENT_SITE", axes = c(1,2)) + 
  scale_color_manual(values = "#016403") +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="Deployment Site",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "none")

h5<-plot_ordination(physeq = HokoDayFIVE, ordination = hoko_pcoa_d5, color = "DECOMP_DAY", shape = "DEPLOYMENT_SITE", axes = c(1,2)) + 
  scale_color_manual(values = "#35CEF2") +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="Deployment Site",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  labs(y="Axis 2 [11.4%]", x = "Axis 1 [41.3%]") +
  theme(legend.position = "none") 

h10<-plot_ordination(physeq = HokoDayTEN, ordination = hoko_pcoa_d10, color = "DECOMP_DAY", shape = "DEPLOYMENT_SITE", axes = c(1,2)) + 
  scale_color_manual(values = "#03045e") +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="Deployment Site",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  labs(y="Axis 2 [23.4%]", x = "Axis 1 [12.2%]") +
  theme(legend.position = "none") 

h15<-plot_ordination(physeq = HokoDay15, ordination = hoko_pcoa_d15, color = "DECOMP_DAY", shape = "DEPLOYMENT_SITE", axes = c(1,2)) + 
  scale_color_manual(values = "#ffb703") +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="Deployment Site",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  labs(y="Axis 2 [12.9%]", x = "Axis 1 [26.7%]") +
  theme(legend.position = "none") 

h20<-plot_ordination(physeq = HokoDay20, ordination = hoko_pcoa_d20, color = "DECOMP_DAY", shape = "DEPLOYMENT_SITE", axes = c(1,2)) + 
  scale_color_manual(values = "#f03b20") +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="Deployment Site",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  labs(y="Axis 2 [12.7%]", x = "Axis 1 [19.3%]") +
  theme(legend.position = "none") 

SekiuDayZERO<-subset_samples(Sekiu.River1, DECOMP_DAY !="Day 5" & DECOMP_DAY !="Day 10" & DECOMP_DAY !="Day 15" & DECOMP_DAY !="Day 20")
SekiuDayZERO<- prune_taxa(taxa_sums(SekiuDayZERO) > 0, SekiuDayZERO)
dim(SekiuDayZERO@otu_table) #8 x 336
Sekiu_pcoa_d0 <- ordinate(physeq = SekiuDayZERO, method = "PCoA", distance = "bray") 

SekiuDayFIVE<-subset_samples(Sekiu.River1, DECOMP_DAY !="Fresh Leaves" & DECOMP_DAY !="Day 10" & DECOMP_DAY !="Day 15" & DECOMP_DAY !="Day 20")
SekiuDayFIVE<- prune_taxa(taxa_sums(SekiuDayFIVE) > 0, SekiuDayFIVE)
dim(SekiuDayFIVE@otu_table) #39 x 1713
Sekiu_pcoa_d5 <- ordinate(physeq = SekiuDayFIVE, method = "PCoA", distance = "bray") 

SekiuDayFIVE_UP<-subset_samples(SekiuDayFIVE, DEPLOYMENT_SITE != "Downstream Sekiu")
SekiuDayFIVE_UP<- prune_taxa(taxa_sums(SekiuDayFIVE_UP) > 0, SekiuDayFIVE_UP)
dim(SekiuDayFIVE_UP@otu_table) #20 x 1371
Sekiu_pcoa_d5_up <- ordinate(physeq = SekiuDayFIVE_UP, method = "PCoA", distance = "bray") 

SekiuDayTEN<-subset_samples(Sekiu.River1, DECOMP_DAY !="Fresh Leaves" & DECOMP_DAY !="Day 5" & DECOMP_DAY !="Day 15" & DECOMP_DAY !="Day 20")
SekiuDayTEN<- prune_taxa(taxa_sums(SekiuDayTEN) > 0, SekiuDayTEN)
dim(SekiuDayTEN@otu_table) #39 x 2499
Sekiu_pcoa_d10 <- ordinate(physeq = SekiuDayTEN, method = "PCoA", distance = "bray") 

SekiuDay15<-subset_samples(Sekiu.River1, DECOMP_DAY !="Fresh Leaves" & DECOMP_DAY !="Day 5" & DECOMP_DAY !="Day 10" & DECOMP_DAY !="Day 20")
SekiuDay15<- prune_taxa(taxa_sums(SekiuDay15) > 0, SekiuDay15)
dim(SekiuDay15@otu_table) #39 x 3791
Sekiu_pcoa_d15 <- ordinate(physeq = SekiuDay15, method = "PCoA", distance = "bray") 

SekiuDay20<-subset_samples(Sekiu.River1, DECOMP_DAY !="Fresh Leaves" & DECOMP_DAY !="Day 5" & DECOMP_DAY !="Day 10" & DECOMP_DAY !="Day 15")
SekiuDay20<- prune_taxa(taxa_sums(SekiuDay20) > 0, SekiuDay20)
dim(SekiuDay20@otu_table) #39 x 4431
Sekiu_pcoa_d20 <- ordinate(physeq = SekiuDay20, method = "PCoA", distance = "bray") 

s0<-plot_ordination(physeq = SekiuDayZERO, ordination = Sekiu_pcoa_d0, color = "DECOMP_DAY", shape = "DEPLOYMENT_SITE", axes = c(1,2)) + 
  scale_color_manual(values = "#016403") +
  scale_shape_manual(values = c(0,15)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="Deployment Site",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "none")

s5<-plot_ordination(physeq = SekiuDayFIVE, ordination = Sekiu_pcoa_d5, color = "DECOMP_DAY", shape = "DEPLOYMENT_SITE", axes = c(1,2)) + 
  scale_color_manual(values = "#35CEF2") +
  scale_shape_manual(values = c(0,15)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="Deployment Site",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  labs(y="Axis 2 [19.0%]", x = "Axis 1 [26.4%]") +
  theme(legend.position = "none") 

s10<-plot_ordination(physeq = SekiuDayTEN, ordination = Sekiu_pcoa_d10, color = "DECOMP_DAY", shape = "DEPLOYMENT_SITE", axes = c(1,2)) + 
  scale_color_manual(values = "#03045e") +
  scale_shape_manual(values = c(0,15)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="Deployment Site",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  labs(y="Axis 2 [12.4%]", x = "Axis 1 [23.7%]") +
  theme(legend.position = "none") 

s15<-plot_ordination(physeq = SekiuDay15, ordination = Sekiu_pcoa_d15, color = "DECOMP_DAY", shape = "DEPLOYMENT_SITE", axes = c(1,2)) + 
  scale_color_manual(values = "#ffb703") +
  scale_shape_manual(values = c(0,15)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="Deployment Site",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  labs(y="Axis 2 [16.4%]", x = "Axis 1 [20.9%]") +
  theme(legend.position = "none") 

s20<-plot_ordination(physeq = SekiuDay20, ordination = Sekiu_pcoa_d20, color = "DECOMP_DAY", shape = "DEPLOYMENT_SITE", axes = c(1,2)) + 
  scale_color_manual(values = "#f03b20") +
  scale_shape_manual(values = c(0,15)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="Deployment Site",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  labs(y="Axis 2 [11.0%]", x = "Axis 1 [22.5%]") +
  theme(legend.position = "none") 

plot_ordination(physeq = SekiuDayFIVE_UP, ordination = Sekiu_pcoa_d5_up, color = "DECOMP_DAY", shape = "ORIGIN_SITE", axes = c(1,2)) + 
  scale_color_manual(values = "#f03b20") +
  scale_shape_manual(values = c(0,15,1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="Deployment Site",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  labs(y="Axis 2 [11.0%]", x = "Axis 1 [22.5%]") +
  theme(legend.position = "none") 

library(gridExtra)
library(grid)
library(cowplot)

plot_grid(h5,h10,h15,h20,s5,s10,s15,s20,nrow = 2, ncol=4, labels=letters[1:8], align="hv")

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#Non-rarefied PCoA
str(HFABac_phylo.rm2) #originates from removing water samples and KatharoSeq samples and then removing ASVs with zero abundance and singletons. 
str(HFABac_phylo.rm2@otu_table) #392 x 22,454 
str(HFABac_phylo.rm2@tax_table) #1 : 22,454
str(HFABac_phylo.rm2@sam_data) #392 x 18 variables -- everything is character
str(HFABac_phylo.rm2@phy_tree) #1 : 22,454 tips 

HFABac_phylo.rm2@sam_data$DECOMP_DAY[HFABac_phylo.rm2@sam_data$DECOMP_DAY=="0"]<-"Fresh Leaves"
HFABac_phylo.rm2@sam_data$DECOMP_DAY[HFABac_phylo.rm2@sam_data$DECOMP_DAY=="5"]<-"Day 5"
HFABac_phylo.rm2@sam_data$DECOMP_DAY[HFABac_phylo.rm2@sam_data$DECOMP_DAY=="10"]<-"Day 10"
HFABac_phylo.rm2@sam_data$DECOMP_DAY[HFABac_phylo.rm2@sam_data$DECOMP_DAY=="15"]<-"Day 15"
HFABac_phylo.rm2@sam_data$DECOMP_DAY[HFABac_phylo.rm2@sam_data$DECOMP_DAY=="20"]<-"Day 20"
HFABac_phylo.rm2@sam_data$DECOMP_DAY<-factor(HFABac_phylo.rm2@sam_data$DECOMP_DAY,levels=c("Fresh Leaves","Day 5","Day 10","Day 15","Day 20"))
HFABac_phylo.rm2@sam_data$Deployed_River_test<-factor(HFABac_phylo.rm2@sam_data$Deployed_River_test,levels=c("Hoko River","Sekiu River"))
HFABac_phylo.rm2@sam_data$ORIGIN_RIVER[HFABac_phylo.rm2@sam_data$ORIGIN_RIVER=="Hoko"]<-"Hoko River"
HFABac_phylo.rm2@sam_data$ORIGIN_RIVER[HFABac_phylo.rm2@sam_data$ORIGIN_RIVER=="Sekiu"]<-"Sekiu River"
HFABac_phylo.rm2@sam_data$ORIGIN_RIVER<-factor(HFABac_phylo.rm2@sam_data$ORIGIN_RIVER,levels=c("Hoko River","Sekiu River"))
HFABac_phylo.rm2@sam_data$TREE_NUMBER<-as.factor(HFABac_phylo.rm2@sam_data$TREE_NUMBER)
HFABac_phylo.rm2@sam_data$LOCATION_CODEwDay0<-as.factor(HFABac_phylo.rm2@sam_data$LOCATION_CODEwDay0)

#Unconstrained Ordinations NR = Nonrarefied
hfa_bac_pcoa_NR <- ordinate(physeq = HFABac_phylo.rm2, method = "PCoA", distance = "bray") 
wu_hfa_bac_pcoa_NR <- ordinate(physeq = HFABac_phylo.rm2, method = "PCoA", distance = "wunifrac") 
uwu_hfa_bac_pcoa_NR <- ordinate(physeq = HFABac_phylo.rm2, method = "PCoA", distance = "unifrac")

#Plot bray curtis pcoa 
plot_ordination(physeq = HFABac_phylo.rm2, ordination = hfa_bac_pcoa_NR, color = "DECOMP_DAY", shape = "Deployed_River_test", axes = c(1,2)) + 
  scale_color_manual(values = c("#016403","#35CEF2","#03045e","#ffb703","#f03b20")) +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="River",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right") 

#plot weighted unifrac pcoa
plot_ordination(physeq = HFABac_phylo.rm2, ordination = wu_hfa_bac_pcoa_NR, color = "DECOMP_DAY", shape = "Deployed_River_test", axes = c(1,2)) + 
  scale_color_manual(values = c("#016403","#35CEF2","#03045e","#ffb703","#f03b20")) +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="River",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right") 

#plot unweighted unifrac pcoa
plot_ordination(physeq = HFABac_phylo.rm2, ordination = uwu_hfa_bac_pcoa_NR, color = "DECOMP_DAY", shape = "Deployed_River_test", axes = c(2,3)) + 
  scale_color_manual(values = c("#016403","#35CEF2","#03045e","#ffb703","#f03b20")) +
  scale_shape_manual(values = c(1,19)) +
  geom_point(aes(color = DECOMP_DAY), alpha = 1, size = 4, stroke=1)+
  theme_test()+
  labs(shape="River",color="Decomposition Day")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right") 

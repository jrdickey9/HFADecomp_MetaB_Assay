custom_col42 = c("#781156","#A51876","#D21E96","#E43FAD","#EA6CC0","#F098D3","#114578","#185EA5","#1E78D2","#3F91E4","#6CABEA","#98C4F0","#117878","#18A5A5","#3FE4E4","#6CEAEA","#98F0F0", "#117845","#18A55E","#1ED278","#3FE491","#6CEAAB","#98F0C4","#787811","#A5A518","#D2D21E","#E4E43F","#EAEA6C","#F0F098","#F7F7C5","#784511","#A55E18","#D2781E","#E4913F","#EAAB6C","#F0C498","#781122","#A5182F","#D21E2C","#E43F5B","#EA6C81","#F098A7")
custom_col37 = c("#781156","#A51876","#D21E96","#E43FAD","#EA6CC0","#F098D3","#114578","#185EA5","#1E78D2","#3F91E4","#6CABEA","#98C4F0","#117878","#18A5A5","#3FE4E4","#6CEAEA","#98F0F0", "#117845","#18A55E","#1ED278","#3FE491","#6CEAAB","#98F0C4","#787811","#A5A518","#D2D21E","#E4E43F","#EAEA6C","#F0F098","#F7F7C5","#784511","#A55E18","#D2781E","#E4913F","#EAAB6C","#F0C498","#781122")
custom_col6<-c("#781156","#A51876","#D21E96","#E43FAD","#EA6CC0","#F098D3")
cus_col20<-c("#781156","#A51876","#D21E96","#E43FAD","#EA6CC0","#F098D3","#114578","#185EA5","#1E78D2","#3F91E4","#6CABEA","#98C4F0","#117878","#18A5A5","#3FE4E4","#6CEAEA","#98F0F0", "#117845","#18A55E","#1ED278")
cus_col9<-c("#781156","#98C4F0","#6CEAEA","#3FE491","#D2D21E","#F7F7C5","#D2781E","#781122","#EA6C81")
rainbow<-c("#E70000","#FF8C00","#FFEF00","#00811F","#0044FF", "#760089","#55CDFC","#FFFFFF","#F7A8B8")
cuz_col20<-c("#98F0C4","#787811","#A5A518","#D2D21E","#E4E43F","#EAEA6C","#F0F098","#F7F7C5","#784511","#A55E18","#D2781E","#E4913F","#EAAB6C","#F0C498","#781122","#A5182F","#D21E2C","#E43F5B","#EA6C81","#F098A7")
custom_col50 = c("#781156","#A51876","#D21E96","#E43FAD","#EA6CC0","#F098D3","#114578","#185EA5","#1E78D2","#3F91E4","#6CABEA","#98C4F0","#117878","#18A5A5","#3FE4E4","#6CEAEA","#98F0F0", "#117845","#18A55E","#1ED278","#3FE491","#6CEAAB","#98F0C4","#787811","#A5A518","#D2D21E","#E4E43F","#EAEA6C","#F0F098","#F7F7C5","#784511","#A55E18","#D2781E","#E4913F","#EAAB6C","#F0C498","#781122","#A5182F","#D21E2C","#E43F5B","#EA6C81","#F098A7","#E70000","#FF8C00","#FFEF00","#00811F","#0044FF", "#760089","#55CDFC","#F7A8B8")

#Filter out low abundance ASVs using abund_phylo. Prune out low abundance taxa and only include Family that contribute more than 2% of the relative abundance of each sample

HFA_Bac_family <- bacterial.phylo.4analysis %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Family)                                      # Sort data frame alphabetically by phylum

str(HFA_Bac_family) #'data.frame':	5562 obs. of  25 variables

#count number of fam to set color pallete
Count = length(unique(HFA_Bac_family$Family))
Count #89

#read the families
unique(HFA_Bac_family$Family)

HFA_Bac_family$Family <- factor(HFA_Bac_family$Family, levels = c(" Acetobacteraceae"," Acholeplasmataceae"," Acidobacteriaceae_(Subgroup_1)"," Aeromonadaceae"," Arcobacteraceae"," Armatimonadaceae"," Bacillaceae"," Bacteriovoracaceae"," Bdellovibrionaceae"," Beijerinckiaceae"," Blfdi19"," Brevibacillaceae"," Burkholderiaceae"," Caulobacteraceae"," Cellvibrionaceae"," Chitinophagaceae"," Chromobacteriaceae"," Comamonadaceae"," Corynebacteriaceae"," Crocinitomicaceae"," Cyclobacteriaceae"," Cytophagaceae"," Deinococcaceae"," Devosiaceae", " Enterococcaceae" ," Erwiniaceae" ," Erysipelotrichaceae"," Exiguobacteraceae", " Flavobacteriaceae" , " Gracilibacteria" ," Haliangiaceae"," Halomonadaceae"," Holophagaceae", " Hymenobacteraceae" ," Hyphomonadaceae" ," Lachnospiraceae" ," Leptotrichiaceae", " LiUU-11-161"," Marinococcaceae"," Methylophagaceae", " Methylophilaceae" ," Microbacteriaceae"," Micromonosporaceae"," Micropepsaceae" ," Microscillaceae" ," Moraxellaceae", " Morganellaceae" ," Myxococcaceae" ," NS11-12_marine_group" ," Nannocystaceae"," Nodosilineaceae" ," Nostocaceae", " Oxalobacteraceae" ," Paludibacteraceae" ," Pedosphaeraceae" ," Peptostreptococcales-Tissierellales"," Pleomorphomonadaceae"," Propionibacteriaceae"," Pseudomonadaceae" ," Rhizobiaceae" ," Rhizobiales_Incertae_Sedis"," Rhodobacteraceae"," Rhodocyclaceae" ," Rubritaleaceae"," SM2D12" ," Salisediminibacteriaceae"," Sandaracinaceae" ," Saprospiraceae" ," Silvanigrellaceae"," Solirubrobacteraceae"," Sphingobacteriaceae"," Sphingomonadaceae"," Spirochaetaceae"," Spirosomaceae"," Spongiibacteraceae"," Sporichthyaceae"," Staphylococcaceae"," Steroidobacteraceae"," Sutterellaceae"," Synechococcaceae"," Unknown_Family"," Verrucomicrobiaceae"," Weeksellaceae"," Xanthobacteraceae"," Xanthomonadaceae"," Yersiniaceae"," env.OPS_17"," uncultured","Unassigned"))

str(HFA_Bac_family)
HFA_Bac_family$DECOMP_DAY[HFA_Bac_family$DECOMP_DAY=="0"]<-"Fresh Leaves"
HFA_Bac_family$DECOMP_DAY[HFA_Bac_family$DECOMP_DAY=="5"]<-"Day 5"
HFA_Bac_family$DECOMP_DAY[HFA_Bac_family$DECOMP_DAY=="10"]<-"Day 10"
HFA_Bac_family$DECOMP_DAY[HFA_Bac_family$DECOMP_DAY=="15"]<-"Day 15"
HFA_Bac_family$DECOMP_DAY[HFA_Bac_family$DECOMP_DAY=="20"]<-"Day 20"
HFA_Bac_family$DECOMP_DAY<-factor(HFA_Bac_family$DECOMP_DAY,levels=c("Fresh Leaves","Day 5","Day 10","Day 15","Day 20"))
HFA_Bac_family$Deployed_River_test<-factor(HFA_Bac_family$Deployed_River_test,levels=c("Hoko River","Sekiu River"))
HFA_Bac_family$ORIGIN_RIVER[HFA_Bac_family$ORIGIN_RIVER=="Hoko"]<-"Hoko River"
HFA_Bac_family$ORIGIN_RIVER[HFA_Bac_family$ORIGIN_RIVER=="Sekiu"]<-"Sekiu River"
HFA_Bac_family$ORIGIN_RIVER<-factor(HFA_Bac_family$ORIGIN_RIVER,levels=c("Hoko River","Sekiu River"))

unique(HFA_Bac_family$Phylum) #choose class to facet
#length(which(HFA_Bac_family$Phylum==" Proteobacteria")) #3673
#length(which(HFA_Bac_family$Phylum==" Firmicutes")) #87
#length(which(HFA_Bac_family$Phylum==" Acidobacteriota")) #6
#length(which(HFA_Bac_family$Phylum==" Campilobacterota")) #64
#length(which(HFA_Bac_family$Phylum==" Armatimonadota")) #1
#length(which(HFA_Bac_family$Phylum==" Bdellovibrionota")) #16
#length(which(HFA_Bac_family$Phylum==" Myxococcota")) #216
#length(which(HFA_Bac_family$Phylum==" Bacteroidota")) #1074
#length(which(HFA_Bac_family$Phylum==" Actinobacteriota")) #113
#length(which(HFA_Bac_family$Phylum==" Deinococcota")) #1
#length(which(HFA_Bac_family$Phylum==" Patescibacteria")) #1
#length(which(HFA_Bac_family$Phylum==" Fusobacteriota")) #2
#length(which(HFA_Bac_family$Phylum==" Cyanobacteria")) #127
#length(which(HFA_Bac_family$Phylum==" Verrucomicrobiota")) #52
#length(which(HFA_Bac_family$Phylum==" Spirochaetota")) #6

# ggplot(HFA_Bac_family, aes(x = DECOMP_DAY, y = Abundance, fill = Family)) + 
#   facet_grid(Deployed_River_test~.) +
#   geom_bar(stat = "identity") +
#   #scale_fill_manual(values = custom_col42) +
#   scale_x_discrete(
#     breaks = c("Fresh Leaves","Day 5","Day 10","Day 15","Day 20"),
#     labels = c("Fresh Leaves","Day 5","Day 10","Day 15","Day 20"), 
#     drop = FALSE
#   ) +
#   # Remove x axis title
#   theme(axis.title.x = element_blank()) + 
#   #
#   guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
#   ylab("Relative Abundance (Families > 1%) \n") 

spatial_plot <- ggplot(subset(HFA_Bac_family, Phylum %in% c(" Proteobacteria"," Bacteroidota")), aes(x=DECOMP_DAY, y=Abundance, fill=Family))+ facet_wrap(~Phylum)
spatial_plot + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values =custom_col50) +
  theme_bw() + 
  theme(legend.position="right") + guides(fill=guide_legend(nrow=10)) +
  theme(text=element_text(size=14)) +
  theme(axis.text.x = element_text(angle=45,hjust = 1))+
  labs(y="Relative Abundance (Families > 1%) \n", x = "Decomposition Day")

spatial_plot1 <- ggplot(subset(HFA_Bac_family, Phylum %in% c(" Proteobacteria")), aes(x=DECOMP_DAY, y=Abundance, fill=Family))+ facet_wrap(~Deployed_River_test)
spatial_plot1 + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values =custom_col50) +
  theme_bw() + 
  theme(legend.position="right") + guides(fill=guide_legend(nrow=10)) +
  theme(text=element_text(size=14)) +
  #theme(axis.text.x = element_text(angle=45,hjust = 1))+
  labs(y="Relative Abundance (Families > 1%) \n", x = "Decomposition Day")

spatial_plot2 <- ggplot(subset(HFA_Bac_family, Phylum %in% c(" Proteobacteria")), aes(x=DECOMP_DAY, y=Abundance, fill=Family))+ facet_wrap(~ORIGIN_RIVER)
spatial_plot2 + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values =custom_col50) +
  theme_bw() + 
  theme(legend.position="right") + guides(fill=guide_legend(nrow=10)) +
  theme(text=element_text(size=14)) +
  #theme(axis.text.x = element_text(angle=45,hjust = 1))+
  labs(y="Relative Abundance (Families > 1%) \n", x = "Decomposition Day")

spatial_plot3 <- ggplot(subset(HFA_Bac_family, Phylum %in% c(" Bacteroidota")), aes(x=DECOMP_DAY, y=Abundance, fill=Family))+ facet_wrap(~Deployed_River_test)
spatial_plot3 + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values =custom_col50) +
  theme_bw() + 
  theme(legend.position="right") + guides(fill=guide_legend(nrow=10)) +
  theme(text=element_text(size=14)) +
  #theme(axis.text.x = element_text(angle=45,hjust = 1))+
  labs(y="Relative Abundance (Families > 1%) \n", x = "Decomposition Day")

spatial_plot4 <- ggplot(subset(HFA_Bac_family, Phylum %in% c(" Bacteroidota")), aes(x=DECOMP_DAY, y=Abundance, fill=Family))+ facet_wrap(~ORIGIN_RIVER)
spatial_plot4 + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values =custom_col50) +
  theme_bw() + 
  theme(legend.position="right") + guides(fill=guide_legend(nrow=10)) +
  theme(text=element_text(size=14)) +
  #theme(axis.text.x = element_text(angle=45,hjust = 1))+
  labs(y="Relative Abundance (Families > 1%) \n", x = "Decomposition Day")

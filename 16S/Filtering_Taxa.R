#This script is used to examine abundant ASVs that were not assigned at the Phylum level
poop<-phyloseq(otu_table(df.bacterial.OTU, taxa_are_rows=FALSE), sample_data(meta_data), tax_table(bac.taxa.m), phy_tree(dated.16Stree))
poop@otu_table[1:10,1:10] #column names are feature IDs. 

#Examining taxonomic ranks to examine where chloroplasts and mitochondria are nested within
#102 archaea, 30,609 bacteria, 19 euks, 259 unassigned #i'm removing everything but bacteria
table(tax_table(poop)[, "Domain"], exclude = NULL)
table(tax_table(poop)[, "Phylum"], exclude = NULL) #2,885 unassigned
table(tax_table(poop)[, "Class"], exclude = NULL) #examine 
table(tax_table(poop)[, "Order"], exclude = NULL) #1081 reads assigned as chloroplast at this taxonomic rank
table(tax_table(poop)[, "Family"], exclude = NULL) #1920 reads assigned as mitochondria at this taxonomic rank
table(tax_table(poop)[, "Genus"], exclude = NULL) 
table(tax_table(poop)[, "species"], exclude = NULL) #2 unidentified

p1<-subset_taxa(poop,  !Domain %in% "Unassigned") #requires dplyr (%in% syntax)
p2<-subset_taxa(p1,  !Domain %in% "Archaea") #be no more!
p3<-subset_taxa(p2,  !Domain %in% "Eukaryota") #sayonara euks!
p4<-subset_taxa(p3,  !Order %in% " Chloroplast") #Chlorplasts be gone! 
p5<-subset_taxa(p4,  !Family %in% " Mitochondria") #Adios amigos!

p5@otu_table[1:10,1:10]
table(tax_table(p5)[, "Domain"], exclude = NULL) #27,608 Bacteria
table(tax_table(p5)[, "Phylum"], exclude = NULL) #2,613 unassigned

Unassigned.phyla<-subset_taxa(p5, !Phylum %in% " Abditibacteriota") 
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Acidobacteriota") 
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Actinobacteriota") 
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Armatimonadota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Bacteroidota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Bdellovibrionota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Campilobacterota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Chloroflexi")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Cyanobacteria")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Deferribacterota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Deferrisomatota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Deinococcota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Dependentiae")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Desulfobacterota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Elusimicrobiota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " FCPU426")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Fibrobacterota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Firmicutes")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Fusobacteriota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Gemmatimonadota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Hydrogenedentes")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Latescibacterota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Margulisbacteria")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " MBNT15")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Methylomirabilota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Myxococcota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " NB1-j")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Nitrospinota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Nitrospirota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Patescibacteria")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Planctomycetota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Proteobacteria")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " RCP2-54")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Rs-K70_termite_group")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " SAR324_clade(Marine_group_B)")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Spirochaetota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Sumerlaeota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Sva0485")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Synergistota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Thermotogota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Verrucomicrobiota")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " WPS-2")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " WS2")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " WS4")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Zixibacteria")
Unassigned.phyla<-subset_taxa(Unassigned.phyla, !Phylum %in% " Cloacimonadota")

table(tax_table(Unassigned.phyla)[, "Phylum"], exclude = NULL)

Unassigned.phyla.asv<-Unassigned.phyla@otu_table
dim(Unassigned.phyla.asv) #428 x 2613
Unassigned.phyla.asv[1:10,1:10]

sort(x=colSums(Unassigned.phyla.asv),decreasing=TRUE) #these are the most abundant unknown taxa at the phylum level

#cd37e3b7102cf737f4e3715bfa2767b5
#Uncultured bacterium clone OTU_5110 16S ribosomal RNA gene, partial sequence

#019aefd085b3aa408f785243898a989f
#Uncultured bacterium clone HAV7D9G03BXD3M 16S ribosomal RNA gene, partial sequence
#Uncultured bacterium clone HAV7D9G03CIG4V 16S ribosomal RNA gene, partial sequence

#ae0c3bd703d23f2b32825b7617388acc
#Uncultured bacterium clone 4362 16S ribosomal RNA gene, partial sequence

#28d722ef40731a8076a2c573c0b5517f
#Uncultured bacterium clone FL11e11_17313 16S ribosomal RNA gene, partial sequence
#uncultured Alphaproteobacteria bacterium partial 16S rRNA gene
#Uncultured bacterium clone OTU_6530 16S ribosomal RNA gene, partial sequence

#c467d4b4b5f787571d0b7d72c1504ecd
#Select seq KX704720.1	Uncultured bacterium clone 2538 16S ribosomal RNA gene, partial sequence
#Select seq LC296049.1	Uncultured prokaryote gene for urease, partial sequence, 16S_OTU:8772
#Uncultured bacterium clone 3508 16S ribosomal RNA gene, partial sequence

#ae45ca8ac6cb92305fb72086f24eb409
#Uncultured bacterium clone FL11e11_17313 16S ribosomal RNA gene, partial sequence

#9605257f3f9190341ec6d92b2b2a7054
#Uncultured bacterium clone OTU_5110 16S ribosomal RNA gene, partial sequence
#Uncultured bacterium clone 22616 16S ribosomal RNA gene, partial sequence
#Uncultured bacterium ZOTU658 gene for 16S rRNA, partial sequence

#0db77bbeae7a423da2597f4cce96259a
#Uncultured bacterium clone 972 16S ribosomal RNA gene, partial sequence
#Uncultured prokaryote gene for urease, partial sequence, 16S_OTU:10882

#39c7c4b5b69817d0e7ddcede3a1a2438
#Uncultured organism clone SBXZ_6098 16S ribosomal RNA gene, partial sequence
#Uncultured organism clone SBXZ_5355 16S ribosomal RNA gene, partial sequence
#Uncultured organism clone SBXZ_2983 16S ribosomal RNA gene, partial sequence

#cbb2217c90d6c95e84f8f1fdd48db5a4
#Uncultured bacterium clone OTU_15793 16S ribosomal RNA gene, partial sequence
#Reclinomonas americana ATCC 50284 mitochondrion, complete genome
#Reclinomonas americana ATCC 50283 mitochondrion, complete genome
#Reclinomonas americana ATCC 50633 mitochondrion, complete genome
#Moramonas marocensis clone Mm_8 mitochondrion, partial sequence
#Uncultured bacterium clone denovo6034 16S ribosomal RNA gene, partial sequence

#cdb96cbe6159d7fa2346100e7a915b26
#Uncultured bacterium clone FL11e11_17313 16S ribosomal RNA gene, partial sequence
#uncultured Alphaproteobacteria bacterium partial 16S rRNA gene
#uncultured Alphaproteobacteria bacterium partial 16S rRNA gene

#0c35177aa1c5b142ca778120e76e4b22
#Uncultured bacterium clone 15687 16S ribosomal RNA gene, partial sequence
#Uncultured bacterium clone OTU_8258 16S ribosomal RNA gene, partial sequence
#Uncultured bacterium clone Unfilt_MembF07 16S ribosomal RNA gene, partial sequence

#90a94bd9fd9ecb49c346d4a57ba7ff8d
#Uncultured bacterium clone OTU_15793 16S ribosomal RNA gene, partial sequence
#Reclinomonas americana ATCC 50284 mitochondrion, complete genome
#Reclinomonas americana ATCC 50283 mitochondrion, complete genome
#Reclinomonas americana ATCC 50633 mitochondrion, complete genome

#0fd260df84f9bfa8ab5d33a1563c2949
#Uncultured bacterium clone NK2_514 16S ribosomal RNA gene, partial sequence
#Uncultured bacterium clone OTU_18217 16S ribosomal RNA gene, partial sequence
#Uncultured bacterium gene for 16S rRNA, partial sequence, clone: Sm4rk109U

#7c81533e05291e23754772eac60d312e
#Uncultured bacterium clone OTU_22622 16S ribosomal RNA gene, partial sequence
#Uncultured bacterium clone OTU_19736 16S ribosomal RNA gene, partial sequence
#Uncultured bacterium clone OTU_19965 16S ribosomal RNA gene, partial sequence

#785c07d7c0810c77225117d297a71112
#Uncultured bacterium clone OTU_1285 16S ribosomal RNA gene, partial sequence
#Uncultured bacterium clone RCP2-70 16S ribosomal RNA gene, partial sequence
#Uncultured bacterium clone 3725 16S ribosomal RNA gene, partial sequence

#da4cd6d3100d96be09f3525010c3e1d4
#Uncultured bacterium clone OTU_23164 16S ribosomal RNA gene, partial sequence
#Uncultured bacterium clone 972 16S ribosomal RNA gene, partial sequence

#5b3b2b92fe836ce4f453ea55fdefd61e
#Uncultured bacterium clone HAV7D9G03CIG4V 16S ribosomal RNA gene, partial sequence
#Uncultured bacterium clone HAV7D9G03BXD3M 16S ribosomal RNA gene, partial sequence

#669ba05a7f8d92a94124929ce3626c88
#Uncultured bacterium clone OTU_2151 16S ribosomal RNA gene, partial sequence
#Uncultured bacterium clone OTU_5991 16S ribosomal RNA gene, partial sequence

#1d888199feff1c9752059c874a2591e6
#Uncultured bacterium clone FL11e11_17313 16S ribosomal RNA gene, partial sequence

#0a38eecc462a378cb513f2e715daf743
#Uncultured prokaryote gene for urease, partial sequence, 16S_OTU:19835
#Uncultured bacterium clone OTU_22167 16S ribosomal RNA gene, partial sequence

#112946ab336062faa6b206d39f22790c
#Uncultured prokaryote gene for urease, partial sequence, 16S_OTU:19699
#Uncultured bacterium clone OTU_8258 16S ribosomal RNA gene, partial sequence
#Uncultured bacterium clone 15687 16S ribosomal RNA gene, partial sequence

#971f2913d861f0229663a4526ee804ab
#Uncultured bacterium clone OTU_22849 16S ribosomal RNA gene, partial sequence
#Uncultured eukaryote clone IIGE3NT01BXI5K_F32_OTU026 12S ribosomal RNA gene, partial sequence; mitochondrial
#Uncultured eukaryote clone ZH1408 16S ribosomal RNA gene, partial sequence; mitochondrial
#Uncultured eukaryote clone IIGE3NT01ADWKQ_F09_OTU026 12S ribosomal RNA gene, partial sequence; mitochondrial

#87c919242ac555b5fbda1625eaf0922b
#Scolecobasidium phaeophorum strain CBS 206.96 small subunit ribosomal RNA gene, partial sequence; mitochondrial
#Trichodelitschia munkii voucher Kruys 201 (UPS) small subunit ribosomal RNA gene, partial sequence; mitochondrial

#4eb3439357a6d353540b91cac3239972
#Uncultured organism clone OTU_1154 16S ribosomal RNA gene, partial sequence
#uncultured bacterium DNA containing 16S-23S intergenic spacer region, clone 43
#Uncultured bacterium clone 2802 16S ribosomal RNA gene, partial sequence

#458d7060d2f403c098c8bc52bc15e41a
#Uncultured prokaryote gene for urease, partial sequence, 16S_OTU:2334
#Uncultured bacterium clone HP_2_uncultured_bacterium_26217 16S ribosomal RNA gene, partial sequence
#Uncultured bacterium clone 3019 16S ribosomal RNA gene, partial sequence
#Uncultured bacterium clone 12487 16S ribosomal RNA gene, partial sequence

#bc29c0e86aef58867529cfd474da3d57
#Uncultured bacterium clone 26432 16S ribosomal RNA gene, partial sequence

#77f8f8987d530682cdbf661e0dc7d381
#Chytriomyces confervae culture CBS:675.73 strain CBS 675.73 mitochondrion, complete genome
#Uncultured bacterium clone denovo7164 16S ribosomal RNA gene, partial sequence
#Spizellomyces palustris culture CBS:455.65 strain CBS 455.65 mitochondrion, complete genome

#881fa8516a280f3e07714b55768a3338
#Uncultured bacterium clone OTU_2151 16S ribosomal RNA gene, partial sequence
#Uncultured bacterium clone OTU_5991 16S ribosomal RNA gene, partial sequence

#7374eb528821b8df58c0ba4c74a71353
#uncultured bacterium partial 16S rRNA gene

#78df6647f6b8e22aad1b1237d1bc107a
#Uncultured bacterium partial 16S rRNA gene, clone O218706E04
#Uncultured prokaryote gene for urease, partial sequence, 16S_OTU:8772

#423ee6d9bf45e99541e8f54a8e13affa
#Uncultured bacterium clone FL11e11_17313 16S ribosomal RNA gene, partial sequence
#uncultured Alphaproteobacteria bacterium partial 16S rRNA gene

#ad5d1e0c2510337ed6c4fab283d15dae
#uncultured bacterium partial 16S rRNA gene

#63a7079e7784147e3f1cff3574aea9a8
#Uncultured bacterium clone FL11e11_17313 16S ribosomal RNA gene, partial sequence

#44a73b043bae664e23c8beca3de04746
#Uncultured bacterium clone 22638 16S ribosomal RNA gene, partial sequence

#6ec96b74df47237dacc1d2abcb3612aa
#Uncultured bacterium clone OTU_5110 16S ribosomal RNA gene, partial sequence

#a9dcd5a98bd9af1880f977f9424b9666 
#Uncultured bacterium clone FW104-043 16S ribosomal RNA gene, partial sequence

#301335cde77a0d707c2f2a7638a41a1a 
#Scolecobasidium phaeophorum strain CBS 206.96 small subunit ribosomal RNA gene, partial sequence; mitochondrial

#8b9b07020402a01a66a91a91c73d7b3d 
#Scolecobasidium phaeophorum strain CBS 206.96 small subunit ribosomal RNA gene, partial sequence; mitochondrial

#8c3529f01f1b68de2f2d3958e0670e6a 
#Uncultured bacterium clone NK2_514 16S ribosomal RNA gene, partial sequence

#40eab637f322d7219342605c49239de7 
#Uncultured bacterium clone cafs939 16S ribosomal RNA gene, partial sequence
#Navicula veneta isolate SZCZR1826 mitochondrion, complete genome

#02a1bd1f1d3321139eb7f402aad2b151 
#Diphylleia rotans isolate NIES-3764 mitochondrion, complete genome

#343aac843d3d7b35b411507c34aec6b3 
#Uncultured bacterium clone HP_2_uncultured_bacterium_26997 16S ribosomal RNA gene, partial sequence

#7b28c1b211656447d0672599e68268b3 
#Uncultured prokaryote gene for urease, partial sequence, 16S_OTU:4072

#4f383882e690d2d87a6fab20eacec066
#Uncultured bacterium clone OTU_5800 16S ribosomal RNA gene, partial sequence

#25abaf02a822597d9f8a64261204524b 
#Uncultured bacterium clone OTU_5110 16S ribosomal RNA gene, partial sequence

#60407e5fdc0519bb0151b05ee5f490b5 
#Picea glauca clone GQ03511_I18 mRNA sequence
#Taphrina wiesneri JCM 22204 mitochondrial DNA, complete sequence

#f1a276064ed37471ec73a61b1a95a4e1 
#Vorticellides astyliformis small subunit ribosomal RNA gene, partial sequence

#2454f25042adf21a60dbb34305cae1e0 
#Scolecobasidium phaeophorum strain CBS 206.96 small subunit ribosomal RNA gene, partial sequence; mitochondrial

#b69aba68b1c94762d41705795590f2ea 
#Uncultured bacterium clone 4362 16S ribosomal RNA gene, partial sequence

#4af5afd0ddd39f37820e619abf4d9e59 
#Uncultured bacterium clone OTU_2151 16S ribosomal RNA gene, partial sequence

#9888bd0a6d06d5ea8fc41bbc6891b8e2 
#Uncultured bacterium clone 3225 16S ribosomal RNA gene, partial sequence

#cb837f65c09020328a050279404f1ff1 
#Uncultured bacterium clone OTU_23273 16S ribosomal RNA gene, partial sequence

#9d91ae29ccbe6738444d66c9dc90cf87 
#Uncultured bacterium clone OTU_2151 16S ribosomal RNA gene, partial sequence

#8e419581ec1d12a9afc0059921b4140b 
#Uncultured bacterium clone 1165 16S ribosomal RNA gene, partial sequence

#dab85ceacce2c97933f6b663da1e1b44 
#Edaphochlamys debaryana mitochondrion, complete genome

#dbaa739f1eb9257d1940679c64f4e173 
#Uncultured bacterium clone OTU_2151 16S ribosomal RNA gene, partial sequence

#3b57fe76dfb2e5a936480272d4eea5fc 
#Uncultured bacterium clone F5K2Q4C04H64I0 16S ribosomal RNA gene, partial sequence

#23ff8284e910ddbe4452310864eecb5f 
#Uncultured bacterium clone FL11e11_17313 16S ribosomal RNA gene, partial sequence

#0ad5534180d222fdf48b87012bf04daf
#Uncultured bacterium clone HAV7D9G03CIG4V 16S ribosomal RNA gene, partial sequence

#e00fad3364ed3f52b31a10da432800cc 
#Uncultured bacterium clone F5K2Q4C04IZ3BC 16S ribosomal RNA gene, partial sequence

#93ad2cb4166368ebe22a2450dcf274aa 
#uncultured Bacterium partial 16S rRNA gene

#dab8bb8c2def5e61fd3531ed09494fb3
#Uncultured bacterium clone 3296 16S ribosomal RNA gene, partial sequence


#df12ef07fcac0d98047f5c3d80794e9d 
#Uncultured bacterium clone FL11e11_17313 16S ribosomal RNA gene, partial sequence

#034f1876d93ee3be044105bba9f87b1c #1500 reads 
#Uncultured bacterium clone OTU23424_AL203_4851880 16S ribosomal RNA gene, partial sequence

#to remove: 
#cbb2217c90d6c95e84f8f1fdd48db5a4
#90a94bd9fd9ecb49c346d4a57ba7ff8d
#87c919242ac555b5fbda1625eaf0922b
#77f8f8987d530682cdbf661e0dc7d381
#301335cde77a0d707c2f2a7638a41a1a 
#8b9b07020402a01a66a91a91c73d7b3d 
#40eab637f322d7219342605c49239de7 
#02a1bd1f1d3321139eb7f402aad2b151 
#60407e5fdc0519bb0151b05ee5f490b5 
#f1a276064ed37471ec73a61b1a95a4e1 
#2454f25042adf21a60dbb34305cae1e0 
#dab85ceacce2c97933f6b663da1e1b44 




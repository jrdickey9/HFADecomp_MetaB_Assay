#install.packages("/Users/lab/vegetarian",repos = NULL, type="source")
library(vegetarian)
dim(std.bac.tab) #326 x 22,328 I think it's OK to do this with Day 0 included as it will just be a column and can be removed downstream if not interesting. At least the calculation of alpha can be done. 

#Alpha Diversity at q= 0, 1, and 2 using non-rarefied, proportional abundandance data with no water samples 
alpha0<-apply(std.bac.tab,1,vegetarian::d,q=0) #species richness
alpha1<-apply(std.bac.tab,1,vegetarian::d,q=1) #exponential Shannon's entropy
alpha2<-apply(std.bac.tab,1,vegetarian::d,q=2) #inverse Simpson's Index

length(alpha0) #326

#metadata
colnames(HFABac_meta_pruned5k) 
Alpha_div_meta<-cbind(HFABac_meta_pruned5k[,4:17],alpha0,alpha1,alpha2)
str(Alpha_div_meta) #Transition to factor for mixed models and boxplots; alpha div scores are numeric.
Alpha_div_meta$Deployed_River_test<-as.factor(Alpha_div_meta$Deployed_River_test)
Alpha_div_meta$DEPLOYMENT_SITE_TEST<-as.factor(Alpha_div_meta$DEPLOYMENT_SITE_TEST)
Alpha_div_meta$ORIGIN_RIVER<-as.factor(Alpha_div_meta$ORIGIN_RIVER)

#remove day zero for the lmer
Ordered_Alpha_div_meta<-Alpha_div_meta[order(Alpha_div_meta$DECOMP_DAY),]
dim(Ordered_Alpha_div_meta)
Ordered_Alpha_div_meta<-Ordered_Alpha_div_meta[19:326,]

Ordered_Alpha_div_meta$DECOMP_DAY<-as.character(Ordered_Alpha_div_meta$DECOMP_DAY)
Ordered_Alpha_div_meta$DECOMP_DAY[Ordered_Alpha_div_meta$DECOMP_DAY=="5"]<-"Day 5"
Ordered_Alpha_div_meta$DECOMP_DAY[Ordered_Alpha_div_meta$DECOMP_DAY=="10"]<-"Day 10"
Ordered_Alpha_div_meta$DECOMP_DAY[Ordered_Alpha_div_meta$DECOMP_DAY=="15"]<-"Day 15"
Ordered_Alpha_div_meta$DECOMP_DAY[Ordered_Alpha_div_meta$DECOMP_DAY=="20"]<-"Day 20"
Ordered_Alpha_div_meta$DECOMP_DAY<-factor(Ordered_Alpha_div_meta$DECOMP_DAY,levels=c("Day 5","Day 10","Day 15","Day 20"))
str(Ordered_Alpha_div_meta)

#examining normality of alpha diversity data
library(fBasics)
hist(Alpha_div_meta$alpha0)
hist(Alpha_div_meta$alpha1) #more skewed
hist(Alpha_div_meta$alpha2) #even more skewed
dagoTest(log(Alpha_div_meta$alpha0))
dagoTest(log(Alpha_div_meta$alpha1))
dagoTest(log(Alpha_div_meta$alpha2)) #normalized

#Alpha Diversity Models; main effects are HFA location code and decomposition day as fixed effects and tree number as random effect. 
alpha0.mod<-lmer(Ordered_Alpha_div_meta$alpha0~Alpha_div_meta$LOCATION_CODE+Alpha_div_meta$DECOMP_DAY+(1+Alpha_div_meta$DEPLOYMENT_SITE|Alpha_div_meta$TREE_NUMBER)) #species richness
alpha1.mod<-lmer(Ordered_Alpha_div_meta$alpha1~Alpha_div_meta$LOCATION_CODE+Alpha_div_meta$DECOMP_DAY+(1+Alpha_div_meta$DEPLOYMENT_SITE|Alpha_div_meta$TREE_NUMBER)) #exponential Shannon's entropy
alpha2.mod<-lmer(Ordered_Alpha_div_meta$alpha2~Alpha_div_meta$LOCATION_CODE+Alpha_div_meta$DECOMP_DAY+(1+Alpha_div_meta$DEPLOYMENT_SITE|Alpha_div_meta$TREE_NUMBER)) #inverse Simpson's Index

#Alpha Diversity Models; main effects are HFA location code and decomposition day as fixed effects and tree number as random effect. 
alpha0.mod1<-lmer(Ordered_Alpha_div_meta$alpha0~Ordered_Alpha_div_meta$LOCATION_CODE+Ordered_Alpha_div_meta$DECOMP_DAY+(1+Ordered_Alpha_div_meta$DEPLOYMENT_SITE|Ordered_Alpha_div_meta$TREE_NUMBER)) #species richness
alpha1.mod2<-lmer(Ordered_Alpha_div_meta$alpha1~Ordered_Alpha_div_meta$LOCATION_CODE+Ordered_Alpha_div_meta$DECOMP_DAY+(1+Ordered_Alpha_div_meta$DEPLOYMENT_SITE|Ordered_Alpha_div_meta$TREE_NUMBER)) #exponential Shannon's entropy
alpha2.mod3<-lmer(Ordered_Alpha_div_meta$alpha2~Ordered_Alpha_div_meta$LOCATION_CODE+Ordered_Alpha_div_meta$DECOMP_DAY+(1+Ordered_Alpha_div_meta$DEPLOYMENT_SITE|Ordered_Alpha_div_meta$TREE_NUMBER)) #inverse Simpson's Index

#analysis of variance
Anova(alpha0.mod) 
Anova(alpha1.mod)
Anova(alpha2.mod)

#analysis of variance
Anova(alpha0.mod1) 
Anova(alpha1.mod2)
Anova(alpha2.mod3)

#Plots
decomp.colors<-c("#016403","#35CEF2","#03045e","#ffb703","#f03b20") 
Hoko.colors<-c("#e97cfa","#7c8afa")
Sekiu.colors<-c("#e88d09","#0e92cb")
River_col<-c("#7c8afa","#e97cfa","#e88d09","#0e92cb")

str(Alpha_div_meta)
Alpha_div_meta$DECOMP_DAY<-as.character(Alpha_div_meta$DECOMP_DAY)
Alpha_div_meta$DEPLOYMENT_SITE<-as.character(Alpha_div_meta$DEPLOYMENT_SITE)
Alpha_div_meta$DECOMP_DAY[Alpha_div_meta$DECOMP_DAY=="0"]<-"Fresh Leaves"
Alpha_div_meta$DECOMP_DAY[Alpha_div_meta$DECOMP_DAY=="5"]<-"Day 5"
Alpha_div_meta$DECOMP_DAY[Alpha_div_meta$DECOMP_DAY=="10"]<-"Day 10"
Alpha_div_meta$DECOMP_DAY[Alpha_div_meta$DECOMP_DAY=="15"]<-"Day 15"
Alpha_div_meta$DECOMP_DAY[Alpha_div_meta$DECOMP_DAY=="20"]<-"Day 20"
Alpha_div_meta$DECOMP_DAY<-factor(Alpha_div_meta$DECOMP_DAY,levels=c("Fresh Leaves","Day 5","Day 10","Day 15","Day 20"))
Alpha_div_meta$DEPLOYMENT_SITE[Alpha_div_meta$DEPLOYMENT_SITE=="HOK1"]<-"Downstream Hoko"
Alpha_div_meta$DEPLOYMENT_SITE[Alpha_div_meta$DEPLOYMENT_SITE=="HOK2"]<-"Upstream Hoko"
Alpha_div_meta$DEPLOYMENT_SITE[Alpha_div_meta$DEPLOYMENT_SITE=="SEK2"]<-"Upstream Sekiu"
Alpha_div_meta$DEPLOYMENT_SITE[Alpha_div_meta$DEPLOYMENT_SITE=="SEK1"]<-"Downstream Sekiu"
Alpha_div_meta$DEPLOYMENT_SITE<-factor(Alpha_div_meta$DEPLOYMENT_SITE,levels=c("Upstream Hoko","Downstream Hoko","Upstream Sekiu","Downstream Sekiu"))

#Species richness across leaf litter decomposition per deployment site
a0.plot<-ggplot(data=Alpha_div_meta, aes(x=DECOMP_DAY,y=alpha0,fill=DECOMP_DAY))+
  geom_boxplot(aes(fill=DECOMP_DAY), alpha=0.5, outlier.size=0) +
  geom_point(pch=21, position=position_jitterdodge()) +
  facet_wrap(~ DEPLOYMENT_SITE)+
  scale_fill_manual(values=decomp.colors) + # Boxplot fill color
  theme_bw() +
  labs(y="Bacterial Species Richness",x="Decomposition Day") +
  theme(legend.position = "right")

a0.plot+ labs(fill="Decomposition Day")+ theme(text = element_text(size = 15),legend.text.align = 0) +theme(axis.text.x = element_text(angle=35,hjust = 1))

#Exponential Shannnon's Diversity across leaf litter decomposition per deployment site
a1.plot<-ggplot(data=Alpha_div_meta, aes(x=DECOMP_DAY,y=alpha1,fill=DECOMP_DAY))+
  geom_boxplot(aes(fill=DECOMP_DAY), alpha=0.5, outlier.size=0) +
  geom_point(pch=21, position=position_jitterdodge()) +
  facet_wrap(~ DEPLOYMENT_SITE)+
  scale_fill_manual(values=decomp.colors) + # Boxplot fill color
  theme_bw() +
  labs(y="Bacterial Shannon's Diversity",x="Decomposition Day") +
  theme(legend.position = "right")

a1.plot+ labs(fill="Decomposition Day")+ theme(text = element_text(size = 15),legend.text.align = 0) +theme(axis.text.x = element_text(angle=35,hjust = 1))

#Exponential Shannnon's Diversity across leaf litter decomposition per deployment site
a2.plot<-ggplot(data=Alpha_div_meta, aes(x=DECOMP_DAY,y=alpha2,fill=DECOMP_DAY))+
  geom_boxplot(aes(fill=DECOMP_DAY), alpha=0.5, outlier.size=0) +
  geom_point(pch=21, position=position_jitterdodge()) +
  facet_wrap(~ DEPLOYMENT_SITE)+
  scale_fill_manual(values=decomp.colors) + # Boxplot fill color
  theme_bw() +
  labs(y="Inverse Simpson's Index",x="Decomposition Day") +
  theme(legend.position = "right")

a2.plot+ labs(fill="Decomposition Day")+ theme(text = element_text(size = 15),legend.text.align = 0) +theme(axis.text.x = element_text(angle=35,hjust = 1))

#Species richness; home vs away per decomposition day for each deployment river
str(Ordered_Alpha_div_meta)
Ordered_Alpha_div_meta$HFA_CODE<-as.character(Ordered_Alpha_div_meta$HFA_CODE)
Ordered_Alpha_div_meta$HFA_CODE[Ordered_Alpha_div_meta$HFA_CODE=="HOME"]<-"Home"
Ordered_Alpha_div_meta$HFA_CODE[Ordered_Alpha_div_meta$HFA_CODE=="AWAY"]<-"Away"
Ordered_Alpha_div_meta$HFA_CODE<-factor(Ordered_Alpha_div_meta$HFA_CODE,levels=c("Home","Away"))

library(tidyr)
small<-cbind(Ordered_Alpha_div_meta[,9:10],Ordered_Alpha_div_meta$Deployed_River_test) 
dim(small)
colnames(small)<-c("CODE","HFA","RIVER")
df<- small %>%
  unite(fill_var, c("HFA","RIVER"))
dim(df)
fill_var_forcolor<-df$fill_var

Ordered_Alpha_div_meta1<-cbind(Ordered_Alpha_div_meta,fill_var_forcolor)
str(Ordered_Alpha_div_meta1)
Ordered_Alpha_div_meta1$fill_var_forcolor[Ordered_Alpha_div_meta1$fill_var_forcolor=="Home_Hoko River"]<-"Hoko River Home Location"
Ordered_Alpha_div_meta1$fill_var_forcolor[Ordered_Alpha_div_meta1$fill_var_forcolor=="Away_Hoko River"]<-"Hoko River Away Location"
Ordered_Alpha_div_meta1$fill_var_forcolor[Ordered_Alpha_div_meta1$fill_var_forcolor=="Away_Sekiu River"]<-"Sekiu River Away Location"
Ordered_Alpha_div_meta1$fill_var_forcolor[Ordered_Alpha_div_meta1$fill_var_forcolor=="Home_Sekiu River"]<-"Sekiu River Home Location"
Ordered_Alpha_div_meta1$fill_var_forcolor<-factor(Ordered_Alpha_div_meta1$fill_var_forcolor,levels=c("Hoko River Home Location","Hoko River Away Location","Sekiu River Home Location","Sekiu River Away Location"))
River_col<-c("#7c8afa","#e97cfa","#e88d09","#0e92cb")
River_1<-c("slateblue4","orchid2","darkgoldenrod2","dodgerblue2")

HFA.a0.plot<-ggplot(data=Ordered_Alpha_div_meta1, aes(x=Deployed_River_test,y=alpha0,fill=fill_var_forcolor))+
  geom_boxplot(aes(fill=fill_var_forcolor), alpha=0.5, outlier.size=0) +
  geom_point(pch=21, position=position_jitterdodge()) +
  facet_grid(~ DECOMP_DAY)+
  scale_fill_manual(values=River_col) + # Boxplot fill color
  theme_bw() +
  labs(y="Bacterial Species Richness", x="Deployment River") +
  theme(legend.position = "right")

HFA.a0.plot+ labs(fill="Home vs Away")+ theme(text = element_text(size = 15),legend.text.align = 0) +theme(axis.text.x = element_text(angle=35,hjust = 1))

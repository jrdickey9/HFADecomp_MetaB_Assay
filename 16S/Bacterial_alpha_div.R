#Alpha Diversity at q= 0, 1, and 2 using non-rarefied, proportional abundandance data with no water samples 
alpha0<-apply(std.bac.tab.df.noH20,1,vegetarian::d,q=0) #species richness
alpha1<-apply(std.bac.tab.df.noH20,1,vegetarian::d,q=1) #exponential Shannon's entropy
alpha2<-apply(std.bac.tab.df.noH20,1,vegetarian::d,q=2) #inverse Simpson's Index

#Species Richness
alpha0<-lmer()

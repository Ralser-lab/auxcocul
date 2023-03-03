#################
### Set Paths ###
#################

proj_dir <- "/Users/aulakhs/Documents/Ralser Lab/Collaborations/Auxotroph co culture screen"
inp_data_dir <- paste0(proj_dir,"/Data/input/trp2_trp4_metabolomics")

###############################################################
### Initialize paths, load functions and graphic parameters ###
###############################################################

source(paste0(proj_dir,"/Rscripts/auxcocul_functions_graphparams_SKA.R"))

####################################################################
### Read in Trp2-Trp4 growth in SM, SM+Trp, SM+AnthranilicAcid data ###
####################################################################

setwd(inp_data_dir)

## Growth of each Auxotroph ##

OD <- read_excel("trpbiosynth_growthcurves_annotation.xlsx",sheet = 5)
colnames(OD)[1]<-"Sample_Name"
OD <- OD %>%
  reshape2::melt(id.vars=c("Sample_Name"),value.name="OD",variable.name="Replicate")

lo <- read_excel("trpbiosynth_growthcurves_annotation.xlsx",sheet = 3,col_names = F)
colnames(lo) <- c("Sample_Name","R1","R2","R3")
lo <- lo %>%
  reshape2::melt(id.vars="Sample_Name",variable.name="Replicate",value.name = "loc_96well")

## Annotaion file

label_anno <-  read_excel("trpbiosynth_growthcurves_annotation.xlsx",sheet = 2,col_names = F)

colnames(label_anno) <- c("Sample_Name","R1","R2","R3")
label_anno<- label_anno %>%
  reshape2::melt(id.vars=c("Sample_Name"),value.name="label",variable.name= "Replicate")

## Merge OD and anno

ODlablo <- merge(OD,label_anno,by=c("Sample_Name","Replicate"))
ODlablo <- merge(ODlablo,lo,by=c("Sample_Name","Replicate"))


# Plot OD
OD2plot <- separate(ODlablo,label, into = c("Strain","Media","Replicate"),sep="_")%>%
  dplyr::select(Strain,Media,OD,Replicate)%>%
  filter(Strain != "Ctrl" )%>%
  unique()%>%
  mutate(Strain = gsub("HIS3","WT",Strain))%>%
  group_by(Media,Strain)%>%
  mutate(Mean_OD600 = mean(OD,na.rm=T),
         SD_OD600 = sd(OD,na.rm=T),
         SE_OD600 = SD_OD600/sqrt(length(OD)))

write.csv(OD2plot,paste0(output_data_dir,"/OD600_trp_biosynthesis_pathway.csv"))

OD2plot$Strain = factor( OD2plot$Strain,
                         levels = c("WT",
                                    "TRP2","TRP4",
                                    "TRP2+TRP4+1:5",
                                    "TRP2+TRP4+1:2",
                                    "TRP2+TRP4+1:1",
                                    "TRP2+TRP4+2:1",
                                    "TRP2+TRP4+5:1" ))
OD2plot$Strain_Media = factor(paste(OD2plot$Strain,OD2plot$Media),
                              levels = c("WT SM",
                                         "WT SM+Trp",
                                         "WT SM+Ant",
                                         "WT SM+Ind",
                                         
                                         "TRP2 SM",
                                         "TRP2 SM+Trp",
                                         "TRP2 SM+Ant",
                                         "TRP2 SM+Ind",
                                         
                                         "TRP4 SM",
                                         "TRP4 SM+Trp",
                                         "TRP4 SM+Ant",
                                         "TRP4 SM+Ind",
                                         
                                         "TRP2+TRP4+1:5 SM",
                                         "TRP2+TRP4+1:2 SM",
                                         "TRP2+TRP4+1:1 SM",
                                         "TRP2+TRP4+2:1 SM",
                                         "TRP2+TRP4+5:1 SM"
                                         
                              ))

setwd(main_fig_dir)

pdf("od600_allcultures_trpbiosynth.pdf",width=12,height=9)

ggplot()+
  geom_point(data = OD2plot,
             aes(x=Strain_Media,
                 y=OD,
                 group=paste(Strain_Media),
                 colour= Media),
             size=2,
             alpha=1,
             position = position_jitterdodge())+
  
  geom_bar(data = unique(OD2plot[,c("Strain","Media","Mean_OD600","SE_OD600","Strain_Media")]),
           aes(x=Strain_Media,
               y=Mean_OD600,
               group=paste(Strain_Media),
               fill= Media,
               colour= Media),
           width = 0.6,
           size=0.25,alpha=0.4,stat="identity" )+
  
  geom_errorbar(data = unique(OD2plot[,c("Strain","Media","Mean_OD600","SE_OD600","Strain_Media")]),
                aes(x=Strain_Media,
                    y=Mean_OD600,
                    ymin=Mean_OD600-SE_OD600,
                    ymax=Mean_OD600+SE_OD600),
                size = 0.25,
                width=0.4,colour="black")+
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_auxcocul(base_size = 20)+
  theme(axis.text.x=element_text(angle=90, hjust = 1),
        legend.position = "bottom")+
  labs(x="",y="OD600 at sampling")+
  
  
  scale_y_continuous(breaks = seq(0, 1.8, by = 0.2))

dev.off()

rm(label_anno)

########################################################
### Read in Trp2-Trp4 Metabolite quantification data ###
########################################################

setwd(inp_data_dir)

dat <- read_excel("trpbiosynth_metabolomics.xlsx",sheet = 4)[2:5]
colnames(dat)[1:4] <-c("loc_96well","Tryptophan_uM","AnthranilicAcid_uM","Indole_uM")


met_df<- merge(dat,ODlablo,by="loc_96well")%>%
  separate(label, into = c("Strain","Media","Replicate"),sep="_")%>%
  reshape2::melt(id.vars=c("Strain","Media","Replicate","OD","loc_96well","Sample_Name"),
                 variable.name = "Metabolite", value.name="Conc_in_uM")%>%
  ## Filter out Indole 
  # filter(Metabolite != "Indole_uM" & Media != "SM+Ind")%>%
  mutate( # Multiple Concentration with 20 --- dilution factor of media -- from Sreejith
    Conc_in_uM = 20*Conc_in_uM,
    Conc_ODnorm = Conc_in_uM / OD)

met_df_sum <- met_df%>%
  # Get stats 
  group_by(Strain,Media,Metabolite)%>%
  mutate(Mean_MetConc = mean(Conc_ODnorm, na.rm=T),
         SD_MetConc = sd(Conc_ODnorm,na.rm=T),
         SE_MetConc = SD_MetConc/sqrt(length(Conc_ODnorm)),
         
         Mean_MetConc_unnorm = mean(Conc_in_uM, na.rm=T),
         SD_MetConc_unnorm = sd(Conc_in_uM,na.rm=T),
         SE_MetConc_unnorm = SD_MetConc_unnorm/sqrt(length(Conc_in_uM)))%>%
  ungroup()%>%
  na.omit()%>%
  mutate(Strain = gsub("HIS3","WT",Strain),
         Media = gsub("SM+Ant","SM+AnthranilicAcid",Media),
         Media = gsub("SM+Trp","SM+Tryptophan",Media),
         Media = gsub("SM+Ind","SM+Indole",Media))

## Convert Strain to factor and fix positions of TRP ratios

met_df_sum$Strain = factor( met_df_sum$Strain,
                            levels = c("WT",
                                       "TRP2","TRP4",
                                       "TRP2+TRP4+1:5",
                                       "TRP2+TRP4+1:2",
                                       "TRP2+TRP4+1:1",
                                       "TRP2+TRP4+2:1",
                                       "TRP2+TRP4+5:1" ))

ODnorm_SM <- ggplot()+
  facet_wrap(c("Metabolite","Media"),scales="free")+
  # barplot
  geom_bar(data = unique(filter(met_df_sum,Media =="SM" & 
                                  !Strain %in% c("TRP2","TRP4"))[,c("Strain","Mean_MetConc","Metabolite","Media")]),
           aes(x = Strain,
               y = Mean_MetConc,
               colour  = as.character(Metabolite),
               fill = as.character(Metabolite)),
           stat="identity",
           size=0.25, 
           width = 0.6,
           alpha = 0.4)+
  # point plot
  
  geom_point(data = filter(met_df_sum,Media =="SM"  & 
                             !Strain %in% c("TRP2","TRP4")),
             aes(x = Strain,
                 y = Conc_ODnorm,
                 colour = Metabolite),
             size=2, position = position_jitterdodge(jitter.width = 0.7))+
  # error bars on top of bar plot and above points 
  
  geom_errorbar(data = unique(filter(met_df_sum,Media =="SM" &
                                       !Strain %in% c("TRP2","TRP4"))[,c("Strain","Mean_MetConc","SE_MetConc","Metabolite","Media")]),
                aes(x = Strain,
                    y = Mean_MetConc,
                    ymin = Mean_MetConc-SE_MetConc, 
                    ymax = Mean_MetConc+SE_MetConc),
                size = 0.25,
                width = 0.4 , colour= "black") +
  # colour scale
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Dark2")[2:4])+
  scale_colour_manual(values = RColorBrewer::brewer.pal(4, "Dark2")[2:4])+
  
  # theme adjustments 
  
   theme_auxcocul(base_size = 20)+
  theme(axis.text.x = element_text(angle=90, hjust =1))+
  # axis labels adjustments 
  
  labs(x="",y="uM of metabolite in supernatant \n (normalised by OD600)",
       title = "OD normalized data SM media only",
       colour = "metabolite",
       fill = "metabolite")

## Plot raw data of cells cultivated in SM 

raw_SM <- ggplot()+
  facet_wrap(c("Metabolite","Media"),scales="free")+
  # barplot
  geom_bar(data = unique(filter(met_df_sum,Media =="SM" & 
                                  !Strain %in% c("TRP2","TRP4"))[,c("Strain","Mean_MetConc_unnorm","Metabolite","Media")]),
           aes(x = Strain,
               y = Mean_MetConc_unnorm,
               colour  = as.character(Metabolite),
               fill = as.character(Metabolite)),
           stat="identity",
           size=0.25, 
           width = 0.6,
           alpha = 0.4)+
  # point plot
  
  geom_point(data = filter(met_df_sum,Media =="SM" & 
                             !Strain %in% c("TRP2","TRP4")),
             aes(x = Strain,
                 y = Conc_in_uM,
                 colour = Metabolite),
             size=2, position = position_jitterdodge(jitter.width = 0.7))+
  # error bars on top of bar plot and above points 
  
  geom_errorbar(data = unique(filter(met_df_sum,Media =="SM" & 
                                       !Strain %in% c("TRP2","TRP4"))[,c("Strain","Mean_MetConc_unnorm","SE_MetConc_unnorm","Metabolite","Media")]),
                aes(x = Strain,
                    y = Mean_MetConc_unnorm,
                    ymin = Mean_MetConc_unnorm-SE_MetConc_unnorm, 
                    ymax = Mean_MetConc_unnorm+SE_MetConc_unnorm),
                size = 0.25,
                width = 0.4 , colour= "black") +
  # colour scale
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Dark2")[2:4])+
  scale_colour_manual(values = RColorBrewer::brewer.pal(4, "Dark2")[2:4])+
  
  # theme adjustments 
  
  theme_auxcocul(base_size = 20)+
  theme(axis.text.x = element_text(angle=90, hjust =1))+
  # axis labels adjustments 
  
  labs(x="",y="uM of metabolite in supernatant \n unnormalised",
       title = "Raw data SM media only",
       colour = "metabolite",
       fill = "metabolite")


zoom_in_anthranilic_acid_ODnorm <- ggplot()+
  facet_wrap(c("Metabolite","Media"),scales="free")+
  # barplot
  geom_bar(data = unique(filter(met_df_sum,Media =="SM" & 
                                  Metabolite == "AnthranilicAcid_uM" &
                                  !Strain %in% c("TRP2","TRP4","WT","TRP2+TRP4+1:5"))[,c("Strain","Mean_MetConc","Metabolite","Media")]),
           aes(x = Strain,
               y = Mean_MetConc,
               colour  = as.character(Metabolite),
               fill = as.character(Metabolite)),
           stat="identity",
           size=0.25, 
           width = 0.6,
           alpha = 0.4)+
  # point plot
  
  geom_point(data = filter(met_df_sum,Media =="SM" & 
                             Metabolite == "AnthranilicAcid_uM" & 
                             !Strain %in% c("TRP2","TRP4","WT","TRP2+TRP4+1:5")),
             aes(x = Strain,
                 y = Conc_ODnorm,
                 colour = Metabolite),
             size=2, position = position_jitterdodge(jitter.width = 0.7))+
  # error bars on top of bar plot and above points 
  
  geom_errorbar(data = unique(filter(met_df_sum,Media =="SM"  & 
                                       Metabolite == "AnthranilicAcid_uM" &
                                       !Strain %in% c("TRP2","TRP4","WT","TRP2+TRP4+1:5"))[,c("Strain","Mean_MetConc","SE_MetConc","Metabolite","Media")]),
                aes(x = Strain,
                    y = Mean_MetConc,
                    ymin = Mean_MetConc-SE_MetConc, 
                    ymax = Mean_MetConc+SE_MetConc),
                size = 0.25,
                width = 0.4 , colour= "black") +
  # colour scale
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Dark2")[2:4])+
  scale_colour_manual(values = RColorBrewer::brewer.pal(4, "Dark2")[2:4])+
  
  # theme adjustments 
  
   theme_auxcocul(base_size = 20)+
  theme(axis.text.x = element_text(angle=90, hjust =1))+
  # axis labels adjustments 
  
  labs(x="",y="uM of metabolite in supernatant \n (normalised by OD600)",
       title = "OD normalized data SM media only",
       colour = "",
       fill = "")


zoom_in_anthranilic_acid_raw <-  ggplot()+
  facet_wrap(c("Metabolite","Media"),scales="free")+
  # barplot
  geom_bar(data = unique(filter(met_df_sum,Media =="SM" & 
                                  Metabolite == "AnthranilicAcid_uM" &
                                  !Strain %in% c("TRP2","TRP4","WT","TRP2+TRP4+1:5"))[,c("Strain","Mean_MetConc_unnorm","Metabolite","Media")]),
           aes(x = Strain,
               y = Mean_MetConc_unnorm,
               colour  = as.character(Metabolite),
               fill = as.character(Metabolite)),
           stat="identity",
           size=0.25, 
           width = 0.6,
           alpha = 0.4)+
  # point plot
  
  geom_point(data = filter(met_df_sum,Media =="SM" & 
                             Metabolite == "AnthranilicAcid_uM" & 
                             !Strain %in% c("TRP2","TRP4","WT","TRP2+TRP4+1:5")),
             aes(x = Strain,
                 y = Conc_in_uM,
                 colour = Metabolite),
             size=2, position = position_jitterdodge(jitter.width = 0.7))+
  # error bars on top of bar plot and above points 
  
  geom_errorbar(data = unique(filter(met_df_sum,Media =="SM"  & 
                                       Metabolite == "AnthranilicAcid_uM" & 
                                       !Strain %in% c("TRP2","TRP4","WT","TRP2+TRP4+1:5"))[,c("Strain","Mean_MetConc_unnorm","SE_MetConc_unnorm","Metabolite","Media")]),
                aes(x = Strain,
                    y = Mean_MetConc_unnorm,
                    ymin = Mean_MetConc_unnorm-SE_MetConc_unnorm, 
                    ymax = Mean_MetConc_unnorm+SE_MetConc_unnorm),
                size = 0.25,
                width = 0.4 , colour= "black") +
  # colour scale
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Dark2")[2])+
  scale_colour_manual(values = RColorBrewer::brewer.pal(4, "Dark2")[2])+
  
  # theme adjustments 
  
  theme_auxcocul(base_size = 20)
  theme(axis.text.x = element_text(angle=90, hjust =1))+
  # axis labels adjustments 
  
  labs(x="",y="uM of metabolite in supernatant \n unnormalised",
       title = "Raw data SM media only",
       colour = "",
       fill = "")


### anthranilic acid in wt, trp2 and trp 4 in media containing tryptophan

anthranilic_acid_ODnorm_inSMtrp <- ggplot()+
  facet_wrap(c("Metabolite","Media"),scales="free")+
  # barplot
  geom_bar(data = unique(filter(met_df_sum,Media =="SM+Trp" & 
                                  Metabolite == "AnthranilicAcid_uM" &
                                  Strain %in% c("TRP2","TRP4","WT"))[,c("Strain","Mean_MetConc","Metabolite","Media")]),
           aes(x = Strain,
               y = Mean_MetConc,
               colour  = as.character(Metabolite),
               fill = as.character(Metabolite)),
           stat="identity",
           size=0.25, 
           width = 0.6,
           alpha = 0.4)+
  # point plot
  
  geom_point(data = filter(met_df_sum,Media =="SM+Trp" & 
                             Metabolite == "AnthranilicAcid_uM" &
                             Strain %in% c("TRP2","TRP4","WT")),
             aes(x = Strain,
                 y = Conc_ODnorm,
                 colour = Metabolite),
             size=2, position = position_jitterdodge(jitter.width = 0.7))+
  # error bars on top of bar plot and above points 
  
  geom_errorbar(data = unique(filter(met_df_sum,Media =="SM+Trp" & 
                                       Metabolite == "AnthranilicAcid_uM" &
                                       Strain %in% c("TRP2","TRP4","WT"))[,c("Strain","Mean_MetConc","SE_MetConc","Metabolite","Media")]),
                aes(x = Strain,
                    y = Mean_MetConc,
                    ymin = Mean_MetConc-SE_MetConc, 
                    ymax = Mean_MetConc+SE_MetConc),
                size = 0.25,
                width = 0.4 , colour= "black") +
  # colour scale
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Dark2")[2:4])+
  scale_colour_manual(values = RColorBrewer::brewer.pal(4, "Dark2")[2:4])+
  
  # theme adjustments 
  
   theme_auxcocul(base_size = 20)+
  theme(axis.text.x = element_text(angle=90, hjust =1))+
  # axis labels adjustments 
  
  labs(x="",y="uM of metabolite in supernatant \n (normalised by OD600)",
       title = "OD normalized data SM media only",
       colour = "",
       fill = "")

## save figures as pdf

setwd(main_fig_dir)

pdf("metabolite_quantification_trp2trp4_in_SM_ODnorm.pdf",width=18,height=8)
ODnorm_SM
dev.off()

pdf("zoom_in_anthranilicacid_trp2trp4_in_SM_ODnorm.pdf",width=6,height=5)
zoom_in_anthranilic_acid_ODnorm
dev.off()

setwd(sup_fig_dir)

pdf("metabolite_quantification_trp2trp4_in_SM_raw.pdf", width=18,height=8)
raw_SM
dev.off()


pdf("zoom_in_anthranilicacid_trp2trp4_in_SM_raw.pdf",width=8,height=6)
zoom_in_anthranilic_acid_raw
dev.off()

pdf("anthranilicacid_in_SMtrp.pdf",width=8,height=6)
anthranilic_acid_ODnorm_inSMtrp
dev.off()

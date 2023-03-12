#################
### Set Paths ###
#################

proj_dir <- "/Users/aulakhs/Documents/Ralser Lab/Collaborations/Auxotroph co culture screen"
inp_data_dir <- paste0(proj_dir,"/data/input/coculture_validation")

###############################################################
### Initialize paths, load functions and graphic parameters ###
###############################################################

source(paste0(proj_dir,"/Rscripts/auxcocul_functions_graphparams_SKA.R"))

###############################################
### Read in OD data from validated strains ###
##############################################

setwd(inp_data_dir)

df_alldata <- vector()

for(i in 1:3){
  
  df_OD = read_excel(paste0("cocul_validation_ODbodrod_passage_",i,".xlsx"),sheet = 1)%>%
    reshape2::melt(id.vars="Time",variable.name = "Culture",value.name ="OD600")
  
  if(i==1){
    negcntrl = df_OD[grep("[_]NA[_]",df_OD$Culture),]
  }
  BOD = read_excel(paste0("cocul_validation_ODbodrod_passage_",i,".xlsx"),sheet = 2)%>%
    reshape2::melt(id.vars="Time",variable.name = "Culture",value.name ="BOD")
  
  
  ROD = read_excel(paste0("cocul_validation_ODbodrod_passage_",i,".xlsx"),sheet = 3)%>%
    reshape2::melt(id.vars="Time",variable.name = "Culture",value.name ="ROD")
  
  
  df = merge(df_OD, BOD, by = c("Time","Culture"))
  df = merge(df, ROD, by= c("Time","Culture"))
  df$Passage = i
  
  df_alldata = rbind(df_alldata,df)
}


df_time_Transform <- df_alldata %>%
  separate(Culture, into = c("Culture","Ratio","Replicate"),sep="_")%>%
  mutate(Replicate = as.numeric(gsub("rep ","",Replicate)),
         Culture = factor(Culture, 
                          levels = c("met3∆-his2∆","met3∆-met1∆","met14∆-met5∆",
                                     "met14∆-trp4∆","met1∆-arg2∆","lys12∆-trp4∆",
                                     "met3∆-trp4∆","trp2∆-trp4∆","met14∆-arg2∆",
                                     "met14∆-his2∆","RB"))
  )

P1_end = max(filter(df_time_Transform,Passage ==1)$Time)
P2_end = max(filter(df_time_Transform,Passage ==2)$Time)

df_time_Transform <- df_time_Transform%>%
  mutate(Time = ifelse(Passage == 1, Time,
                       ifelse(Passage ==2, Time+P1_end,
                              ifelse(Passage ==3, Time+P1_end+P2_end,Time))))

setwd(main_fig_dir)

pdf("stability_validated_cocultures.pdf",width=15,height=10)
ggplot(filter(df_time_Transform,!Culture %in% c("met14∆-his2∆") & Passage != 3),
       aes())+
  facet_wrap(c("Culture","Ratio"),ncol=3,scales="free")+
  geom_point(aes(x=Time, y= OD600),colour = "gray")+
  geom_point(aes(x=Time, y= BOD),colour = "#155289",alpha=0.3)+
  geom_point(aes(x=Time, y= ROD),colour = "#94221F",alpha=0.3)+
  theme_auxcocul()
dev.off()

## Plot negative controls

negcntrl <- negcntrl %>%
  separate(Culture, into = c("Culture","Ratio","Replicate"),sep="_")%>%
  mutate(Replicate = as.numeric(gsub("rep ","",Replicate)))

pdf("negcontrols_stability_validated_cocultures.pdf",width=8,height=5)
ggplot(negcntrl,
       aes( x= Time,
            y= OD600,
            colour = Culture))+
  geom_point(alpha= 0.5)+
  ylim(0,0.7)+
  scale_colour_manual(values = c(RColorBrewer::brewer.pal(8, "Dark2"),
                                 RColorBrewer::brewer.pal(3, "Set2"))) +
  theme_auxcocul()
dev.off()



## Plot positive controls

poscntrl <- filter(df_time_Transform,Culture == "RB" & Time < 95)

pdf("poscontrols_stability_validated_cocultures.pdf",width=8,height=3)
ggplot(poscntrl)+
  facet_wrap(c("Culture","Ratio"))+
  geom_point(aes(x=Time, y= OD600),colour = "gray")+
  geom_point(aes(x=Time, y= BOD),colour = "#155289",alpha=0.3)+
  geom_point(aes(x=Time, y= ROD),colour = "#94221F",alpha=0.3)+
  ylim(0,0.7)+
  theme_auxcocul()
dev.off()



#############################################
### Compute mean and sd of all replicates ###
#############################################

df_time_Transform_smry <- df_time_Transform%>%
  group_by(Culture, Ratio,Replicate,Time)%>%
  mutate(ROD_pcnt = ROD / (BOD+ROD),
         BOD_pcnt = BOD / (BOD+ROD))%>%
  ungroup()%>%
  group_by(Culture, Ratio, Time)%>%
  summarise(Mean_OD600 = mean(OD600,na.rm=T),
            SD_OD600 = sd(OD600,na.rm=T),
            
            Mean_ROD = mean(ROD, na.rm=T),
            SD_ROD = sd(ROD,na.rm=T),
            
            Mean_BOD = mean(BOD,na.rm=T),
            SD_BOD = sd(BOD,na.rm=T),
            
            Mean_ROD_pcnt = 100*mean(ROD_pcnt,na.rm=T),
            SD_ROD_pcnt = 100*sd(ROD_pcnt,na.rm=T),
            
            Mean_BOD_pcnt = 100*mean(BOD_pcnt,na.rm=T),
            SD_BOD_pcnt = 100*sd(BOD_pcnt,na.rm=T))%>%
  ungroup()



pdf("stability_validated_cocultures_summariseddata.pdf",width = 10,height=12)
ggplot(filter(df_time_Transform_smry,!Culture %in% c("met14∆-his2∆") & Time < 95))+
  geom_point(aes(x = Time,
                 y = Mean_ROD),
             alpha = 0.3,
             size = 1,
             colour = "#94221F")+
  geom_line(aes(x = Time,
                y = Mean_ROD),
            colour = "#94221F",
            alpha = 0.7,
            size = 0.5)+
  geom_point(aes(x = Time,
                 y = Mean_BOD),
             alpha = 0.3,
             size = 1,
             colour = "#155289")+
  geom_line(aes(x = Time,
                y = Mean_BOD),
            colour = "#155289",
            alpha = 0.7,
            size = 0.5)+
  geom_point(aes(x= Time,
                 y= Mean_OD600),
             size = 1,
             colour = "gray", 
             alpha = 0.3)+
  geom_line(aes( x = Time,
                 y = Mean_OD600),
            colour = "gray",
            alpha = 0.7,
            size = 0.5)+
  facet_wrap(c("Culture","Ratio"),ncol=3,scales = "free")+
  theme_auxcocul()

dev.off()

## trp2-trp4 3 passages 


pdf(paste0(sup_fig_dir,"/stability_validated_cocultures_summariseddata_trp2trp4_3passages.pdf"),width = 4,height=4)
ggplot(filter(df_time_Transform_smry,Culture == "trp2∆-trp4∆"))+
  geom_point(aes(x = Time,
                 y = Mean_ROD),
             alpha = 0.3,
             size = 1,
             colour = "#94221F")+
  geom_line(aes(x = Time,
                y = Mean_ROD),
            colour = "#94221F",
            alpha = 0.7,
            size = 1)+
  geom_point(aes(x = Time,
                 y = Mean_BOD),
             alpha = 0.3,
             size = 1,
             colour = "#155289")+
  geom_line(aes(x = Time,
                y = Mean_BOD),
            colour = "#155289",
            alpha = 0.7,
            size = 1)+
  geom_point(aes(x= Time,
                 y= Mean_OD600),
             size = 1,
             colour = "gray", 
             alpha = 0.3)+
  geom_line(aes( x = Time,
                 y = Mean_OD600),
            colour = "gray",
            alpha = 0.7,
            size = 1)+
  facet_wrap(c("Culture","Ratio"),ncol=3,scales = "free")+
  theme_auxcocul()

dev.off()
########################################
### Plot fraction at end poitn P1 P2 ###
########################################

df_time_Transform_endpoint <- df_time_Transform_smry%>%
                              filter(Time %in% c(47, 94) & Culture != "met14∆-his2∆")%>%
                              mutate(Culture_Ratio = factor(paste(Culture,Ratio),
                                                            levels = c("RB 5:1",
                                                                       "RB 1:1",
                                                                       "RB 1:5",
                                                                       "met14∆-arg2∆ 1:5", "trp2∆-trp4∆ 1:1",
                                                                       "met3∆-trp4∆ 2:1", "lys12∆-trp4∆ 2:1", 
                                                                       "met1∆-arg2∆ 1:10","met14∆-trp4∆ 2:1",
                                                                       "met14∆-met5∆ 1:5","met3∆-met1∆ 1:5",
                                                                       "met3∆-his2∆ 1:5"

                                                            )))




pdf(paste0(main_fig_dir,"/stability_validated_cocultures_endpoint_fractions.pdf"),width = 8, height=9.6)
ggplot(df_time_Transform_endpoint)+
  
      # RFP
      geom_point(aes(x = Mean_ROD_pcnt,
                     y = Culture_Ratio,
                     shape = factor(Time)),
                     size = 3,
                     colour = "#94221F")+
      geom_line(aes(x = Mean_ROD_pcnt,
                    y = Culture_Ratio),
                    colour = "#94221F")+
      # BFP
      geom_point(aes(x = Mean_BOD_pcnt,
                     y = Culture_Ratio,
                     shape = factor(Time)),
                     size = 3,
                     position = position_nudge(y = 0.15),
                     colour = "#155289")+
      geom_line(aes(x = Mean_BOD_pcnt,
                    y = Culture_Ratio),
                    colour = "#155289",
                    position = position_nudge(y = 0.15))+
  theme_auxcocul()
dev.off()


write.csv(df_time_Transform_smry,paste0(output_data_dir,"/stabilitycheck_validated_cocultures_alldata_summarised.csv"),row.names = T)


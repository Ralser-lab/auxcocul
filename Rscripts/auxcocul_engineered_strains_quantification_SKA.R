#################
### Set Paths ###
#################

proj_dir <- "/Users/aulakhs/Documents/Ralser Lab/Collaborations/Auxotroph co culture screen"
inp_data_dir <- paste0(proj_dir,"/data/input/engineered_strains_bioproduction")


###############################################################
### Initialize paths, load functions and graphic parameters ###
###############################################################

source(paste0(proj_dir,"/Rscripts/auxcocul_functions_graphparams_SKA.R"))

####################################################################
### Read in Trp2-Trp4 growth in SM, SM+Trp, SM+Anthranilate data ###
####################################################################

setwd(inp_data_dir)

# Sheet 3 - OD, 4- Glucose, 5- Beta Alanine 6- 3HP

ES_ODs <- read_excel("engineered_strains_glucose_balanine_msa_quantification_HPJune2022.xlsx",sheet = 3) %>%
          reshape2::melt(id.vars="Sample")
colnames(ES_ODs) <- c("Sample","Timepoint_Replicate","OD600")

ES_Glu <- read_excel("engineered_strains_glucose_balanine_msa_quantification_HPJune2022.xlsx",sheet = 4) %>%
  reshape2::melt(id.vars="Sample")
colnames(ES_Glu) <- c("Sample","Timepoint_Replicate","Glucose")

ES_bAla <- read_excel("engineered_strains_glucose_balanine_msa_quantification_HPJune2022.xlsx",sheet = 5) %>%
  reshape2::melt(id.vars="Sample")
colnames(ES_bAla) <- c("Sample","Timepoint_Replicate","beta Alanine")

ES_msa <- read_excel("engineered_strains_glucose_balanine_msa_quantification_HPJune2022.xlsx",sheet = 7) %>%
  reshape2::melt(id.vars="Sample")
colnames(ES_msa) <- c("Sample","Timepoint_Replicate","MSA")


ES_data <- merge(ES_ODs,ES_Glu, by=c("Sample","Timepoint_Replicate"))
ES_data <- merge(ES_data,ES_bAla, by=c("Sample","Timepoint_Replicate"))
ES_data <- merge(ES_data,ES_msa, by=c("Sample","Timepoint_Replicate"))%>%
         #  filter(OD600 > 0.3)%>%
           reshape2::melt(id.vars=c("Sample","Timepoint_Replicate"),variable.name="parameter")%>%
           separate(Timepoint_Replicate, into = c("Timepoint","Replicate"),sep="_")

ES_data_ODnorm <- ES_data

ES_data<- ES_data%>%
           group_by(Sample,Timepoint,parameter)%>%
           mutate(Mean = mean(value,na.rm=T),
                     SD = sd(value,na.rm=T),
                     SE = SD/sqrt(length(value)))%>%
           mutate(col_fac = ifelse(grepl("BY",Sample),"WT",
                                   ifelse(grepl("RLA",Sample),"WT",
                                          ifelse(grepl("sHP4",Sample),
                                                 "WT",
                                          ifelse(grepl("△HIS2-△MET3",Sample),
                                                 "His2 Met3",
                                                 ifelse(grepl("△MET14-△MET5",Sample),
                                                              "Met14 Met5",
                                                        ifelse(grepl("△MET14-△TRP4",Sample),
                                                               "Met14 Trp2",
                                                               ifelse(grepl("△MET3-△HIS2",Sample),
                                                                      "His2 Met3",
                                                                      ifelse(grepl("△MET5-△MET14",Sample),
                                                                             "Met14 Met5",
                                                                             ifelse(grepl("△TRP4-△MET14",Sample),
                                                                                    "Met14 Trp2",
                                                                                    ifelse(grepl("△TRP4-△TRP2",Sample),
                                                                                           "Trp2 Trp4",
                                                                                           ifelse(grepl("△TRP2-△TRP4",Sample),
                                                                                                  "Trp2 Trp4",NA))))))))))))%>%
  filter(!Sample %in% c("sHP413","△MET14-△MET5 20:1","△MET14-△MET5 1:5","△MET14-△MET5 1:20",
                        "△MET5-△MET14 20:1","△MET5-△MET14 5:1","△MET5-△MET14 1:10"))%>%
  mutate(Sample = factor(Sample,
                         levels = c("RLAs879","BY-BY 10:1","BY-BY 1:1","BY-BY 1:10",
                                    # His3 Met3
                                    "△HIS2-△MET3 20:1","△HIS2-△MET3 5:1","△HIS2-△MET3 1:10",
                                    "△MET3-△HIS2 1:20","△MET3-△HIS2 1:5","△MET3-△HIS2 10:1",

                                    #Met14 Trp4
                                    "△MET14-△TRP4 10:1","△MET14-△TRP4 1:2","△MET14-△TRP4 1:10",
                                    "△TRP4-△MET14 1:10", "△TRP4-△MET14 2:1","△TRP4-△MET14 10:1",
                                    
                                    # Trp2 Trp
                                    "△TRP2-△TRP4 10:1","△TRP2-△TRP4 1:1","△TRP2-△TRP4 1:10",  
                                    "△TRP4-△TRP2 1:10" ,"△TRP4-△TRP2 1:1", "△TRP4-△TRP2 10:1"
                                    )))

setwd(main_fig_dir)
## plot
pdf(paste0(sup_fig_dir,"/engineered_strains_gluc_bala_msa_quantification_raw.pdf"),width=25,height=15)
ggplot()+
  facet_grid(parameter~Timepoint,scales="free")+
  # barplot
  geom_bar(data = unique(filter(ES_data, !Sample %in% c("BY-BY 10:1","BY-BY 1:1","BY-BY 1:10"))[,c("Sample","Mean","col_fac","parameter","Timepoint")] ),
                   aes(x=Sample,
                       y=Mean,
                       fill=col_fac),
           stat = "identity",
           width = 0.6,
           alpha = 0.4,
           size=0.1)+
  # pointplot
  geom_point(data = filter(ES_data, !Sample %in% c("BY-BY 10:1","BY-BY 1:1","BY-BY 1:10")),
           aes(x=Sample,
               y=value,
               colour=col_fac),
           position = position_jitterdodge(),
           size = 2)+
  # errorbars
  geom_errorbar(data = unique(filter(ES_data, !Sample %in% c("BY-BY 10:1","BY-BY 1:1","BY-BY 1:10"))[,c("Sample","Mean","SE","col_fac","parameter","Timepoint")] ),
                aes(x=Sample,
                    y=Mean,
                    ymin = Mean-SE,
                    ymax = Mean+SE),
                width = 0.4)+
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_colour_brewer(palette = "Dark2", direction = -1)+
  theme_auxcocul()+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle=90,hjust = 1))
dev.off()

pdf(paste0(main_fig_dir,"/engineered_strains_48h_only_gluc_bala_msa_quantification_raw.pdf"),width=10,height=20)
ggplot()+
         facet_grid(parameter~Timepoint,scales="free")+
         # barplot
         geom_bar(data = unique(filter(ES_data, Timepoint=="48h" &
                                         !Sample %in% c("BY-BY 10:1","BY-BY 1:1","BY-BY 1:10"))[,c("Sample","Mean","col_fac","parameter","Timepoint")] ),
                  aes(x=Sample,
                      y=Mean,
                      fill=col_fac),
                  stat = "identity",
                  width = 0.6,
                  alpha = 0.4,
                  size=0.1)+
         # pointplot
         geom_point(data = filter(ES_data,Timepoint=="48h" &
                                    !Sample %in% c("BY-BY 10:1","BY-BY 1:1","BY-BY 1:10")),
                    aes(x=Sample,
                        y=value,
                        colour=col_fac),
                    position = position_jitterdodge(),
                    size = 2)+
         # errorbars
         geom_errorbar(data = unique(filter(ES_data, Timepoint=="48h" &
                                              !Sample %in% c("BY-BY 10:1","BY-BY 1:1","BY-BY 1:10"))[,c("Sample","Mean","SE","col_fac","parameter","Timepoint")] ),
                       aes(x=Sample,
                           y=Mean,
                           ymin = Mean-SE,
                           ymax = Mean+SE),
                       width = 0.4)+
         scale_fill_brewer(palette = "Dark2", direction = -1)+
         scale_colour_brewer(palette = "Dark2", direction = -1)+
         theme_auxcocul()+
         theme(legend.position = "bottom",
               axis.text.x = element_text(angle=90,hjust = 1))
dev.off()


pdf(paste0(sup_fig_dir,"/wt_strains_24h_only_glu_bala_msa_quantification_raw.pdf"),width=4,height=9)

ggplot()+
  facet_grid(parameter~Timepoint,scales="free")+
  # barplot
  geom_bar(data = unique(filter(ES_data, Timepoint=="24h" &
                                  Sample %in% c("RLAs879","BY-BY 10:1","BY-BY 1:1","BY-BY 1:10"))[,c("Sample","Mean","col_fac","parameter","Timepoint")] ),
           aes(x=Sample,
               y=Mean,
               fill=col_fac),
           stat = "identity",
           width = 0.6,
           alpha = 0.4,
           size=0.1)+
  # pointplot
  geom_point(data = filter(ES_data,Timepoint=="24h" &
                             Sample %in% c("RLAs879","BY-BY 10:1","BY-BY 1:1","BY-BY 1:10")),
             aes(x=Sample,
                 y=value,
                 colour=col_fac),
             position = position_jitterdodge(),
             size = 2)+
  # errorbars
  geom_errorbar(data = unique(filter(ES_data, Timepoint=="24h" &
                                       Sample %in% c("RLAs879","BY-BY 10:1","BY-BY 1:1","BY-BY 1:10"))[,c("Sample","Mean","SE","col_fac","parameter","Timepoint")] ),
                aes(x=Sample,
                    y=Mean,
                    ymin = Mean-SE,
                    ymax = Mean+SE),
                width = 0.4)+
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_colour_brewer(palette = "Dark2", direction = -1)+
  theme_auxcocul()+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle=90,hjust = 1))
dev.off()


##########################
### OD normalized data ###
##########################


ES_data_ODnorm<- ES_data_ODnorm%>%
                 reshape2::dcast(Sample+Timepoint+Replicate ~ parameter, value.var="value")%>%
                 mutate(Glucose_pOD = Glucose/OD600,
                        bAlanine_pOD =`beta Alanine`/OD600,
                        MSA_pOD = MSA/OD600
                        )%>%
                  dplyr::select(Sample,Timepoint,Replicate,Glucose_pOD, bAlanine_pOD,MSA_pOD)%>%
                  reshape2::melt(id.vars= c("Sample","Timepoint","Replicate"))%>%
                  group_by(variable,Sample,Timepoint)%>%
                  mutate(Mean_ODnorm = mean(value,na.rm=T),
                         SD_ODnorm = sd(value,na.rm=T),
                         SE_ODnorm = SD_ODnorm/sqrt(length(value)))%>%
                 ungroup()%>%
                 mutate(col_fac = ifelse(grepl("BY",Sample),"WT",
                                  ifelse(grepl("RLA",Sample),"WT",
                                         ifelse(grepl("sHP4",Sample),
                                                "WT",
                                                ifelse(grepl("△HIS2-△MET3",Sample),
                                                       "His2 Met3",
                                                       ifelse(grepl("△MET14-△MET5",Sample),
                                                              "Met14 Met5",
                                                              ifelse(grepl("△MET14-△TRP4",Sample),
                                                                     "Met14 Trp2",
                                                                     ifelse(grepl("△MET3-△HIS2",Sample),
                                                                            "His2 Met3",
                                                                            ifelse(grepl("△MET5-△MET14",Sample),
                                                                                   "Met14 Met5",
                                                                                   ifelse(grepl("△TRP4-△MET14",Sample),
                                                                                          "Met14 Trp2",
                                                                                          ifelse(grepl("△TRP4-△TRP2",Sample),
                                                                                                 "Trp2 Trp4",
                                                                                                 ifelse(grepl("△TRP2-△TRP4",Sample),
                                                                                                            "Trp2 Trp4",NA))))))))))))%>%
              filter(!Sample %in% c("sHP413","△MET14-△MET5 20:1","△MET14-△MET5 1:5","△MET14-△MET5 1:20",
                                    "△MET5-△MET14 20:1","△MET5-△MET14 5:1","△MET5-△MET14 1:10"))%>%
              mutate(Sample = factor(Sample,
                                     levels = c("RLAs879","BY-BY 10:1","BY-BY 1:1","BY-BY 1:10",
                                                # His3 Met3
                                                "△HIS2-△MET3 20:1","△HIS2-△MET3 5:1","△HIS2-△MET3 1:10",
                                                "△MET3-△HIS2 1:20","△MET3-△HIS2 1:5","△MET3-△HIS2 10:1",
                                                
                                                #Met14 Trp4
                                                "△MET14-△TRP4 10:1","△MET14-△TRP4 1:2","△MET14-△TRP4 1:10",
                                                "△TRP4-△MET14 1:10", "△TRP4-△MET14 2:1","△TRP4-△MET14 10:1",
                                                
                                                
                                                # Trp2 Trp
                                                "△TRP2-△TRP4 10:1","△TRP2-△TRP4 1:1","△TRP2-△TRP4 1:10",  
                                                "△TRP4-△TRP2 1:10" ,"△TRP4-△TRP2 1:1", "△TRP4-△TRP2 10:1"
                                     )))



pdf(paste0(main_fig_dir,"/engineered_strains_48h_only_gluc_bala_msa_quantification_ODnorm.pdf"),width=10,height=15.5)

ggplot()+
  facet_grid(variable~Timepoint,scales="free")+
  # barplot
  geom_bar(data = unique(filter(ES_data_ODnorm, Timepoint=="48h" &
                                  !Sample %in% c("BY-BY 10:1","BY-BY 1:1","BY-BY 1:10"))[,c("Sample","Mean_ODnorm","col_fac","variable","Timepoint")] ),
           aes(x=Sample,
               y=Mean_ODnorm,
               fill=col_fac),
           stat = "identity",
           width = 0.6,
           alpha = 0.4,
           size=0.1)+
  # pointplot
  geom_point(data = filter(ES_data_ODnorm,  Timepoint=="48h" & !Sample %in% c("BY-BY 10:1","BY-BY 1:1","BY-BY 1:10")),
             aes(x=Sample,
                 y=value,
                 colour=col_fac),
             position = position_jitterdodge(),
             size = 2)+
  # errorbars
  geom_errorbar(data = unique(filter(ES_data_ODnorm,  Timepoint=="48h" &
                                     !Sample %in% c("BY-BY 10:1","BY-BY 1:1","BY-BY 1:10"))[,c("Sample","Mean_ODnorm","SE_ODnorm","col_fac","variable","Timepoint")] ),
                aes(x=Sample,
                    y=Mean_ODnorm,
                    ymin = Mean_ODnorm - SE_ODnorm,
                    ymax = Mean_ODnorm + SE_ODnorm),
                width = 0.4)+
  scale_fill_brewer(palette = "Dark2", direction = -1)+
  scale_colour_brewer(palette = "Dark2", direction = -1)+
  theme_auxcocul()+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle=90,hjust = 1))

dev.off()


## Save summarised data

write.csv(filter(ES_data, Timepoint !="72h"),
          paste0(output_data_dir,"/engineered_strains_glucose_bala_msa_quantification.csv"))


write.csv(filter(ES_data_ODnorm, Timepoint !="72h"),
          paste0(output_data_dir,"/engineered_strains_glucose_bala_msa_quantification_ODnorm.csv"))

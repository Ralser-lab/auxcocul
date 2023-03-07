#################
### Set Paths ### 
#################

proj_dir <- "/Users/aulakhs/Documents/Ralser Lab/Collaborations/Auxotroph co culture screen"

inp_data_dir <- paste0(proj_dir,"/Data/input/coculture_validation")

###############################################################
### Initialize paths, load functions and graphic parameters ###
###############################################################

source(paste0(proj_dir,"/Rscripts/auxcocul_functions_graphparams_SKA.R"))

##############################
### Read in processed data ###
##############################


valcc_OD <- read_excel(paste0(inp_data_dir,"/cocul_validation_OD_laracocultures.xlsx"), sheet = 1)%>%
            data.frame()%>%
            reshape2::melt(id.vars = c("monoculture_or_coculture","group","cultype"), value.name = "OD600_48h")%>%
            filter(!grepl("COLLECTION",monoculture_or_coculture) & 
                     !(group %in% c("his2Δ-met14Δ","arg2Δ-met5Δ","arg2Δ-met3Δ")))%>%
            dplyr::group_by(monoculture_or_coculture)%>%
            dplyr::mutate(Mean_OD_48h = mean(OD600_48h,na.rm = T),
                   SD_OD_48h = sd(OD600_48h, na.rm=T))%>%
            dplyr::ungroup()%>%
            mutate(group_culture_reord = factor(paste(group,monoculture_or_coculture),
                                                levels = c("WT BY4741",
                                                           "WT BY4741 pHLUM",
                                                           "met3Δ-his2Δ met3Δ-his2Δ",
                                                           "met3Δ-his2Δ met3Δ","met3Δ-his2Δ his2Δ",
                                                           "met3Δ-met1Δ met3Δ-met1Δ",
                                                           "met3Δ-met1Δ met3Δ","met3Δ-met1Δ met1Δ",
                                                           "met14Δ-met5Δ met14Δ-met5Δ",
                                                           "met14Δ-met5Δ met14Δ","met14Δ-met5Δ met5Δ",
                                                           "met14Δ-trp4Δ met14Δ-trp4Δ",
                                                           "met14Δ-trp4Δ met14Δ","met14Δ-trp4Δ trp4Δ",
                                                           "met1Δ-arg2Δ met1Δ-arg2Δ",
                                                           "met1Δ-arg2Δ met1Δ","met1Δ-arg2Δ arg2Δ",
                                                           "lys12Δ-trp4Δ lys12Δ-trp4Δ",
                                                           "lys12Δ-trp4Δ lys12Δ","lys12Δ-trp4Δ trp4Δ",
                                                           "met3Δ-trp4Δ met3Δ-trp4Δ",
                                                           "met3Δ-trp4Δ met3Δ","met3Δ-trp4Δ trp4Δ",
                                                           "trp2Δ-trp4Δ trp2Δ-trp4Δ",
                                                           "trp2Δ-trp4Δ trp2Δ","trp2Δ-trp4Δ trp4Δ",
                                                           "met14Δ-arg2Δ met14Δ-arg2Δ",
                                                           "met14Δ-arg2Δ met14Δ", "met14Δ-arg2Δ arg2Δ"
                                                     
                                                    )))
           
          
            
pdf(paste0(sup_fig_dir,"/validated_cocultures_OD600_48h.pdf"),width = 10, height = 6)
ggplot()+
  ## point
  geom_point(data = valcc_OD,
             aes( x = group_culture_reord,
                  y = OD600_48h,
                  colour = cultype),
                  position = position_jitterdodge(jitter.width = 0.6 ))+
  ## bar 
  
  geom_bar(data = unique(valcc_OD[,c("monoculture_or_coculture","Mean_OD_48h","cultype","group_culture_reord")]),
             aes( x = group_culture_reord,
                  y = Mean_OD_48h,
                  colour  = cultype,
                  fill = cultype),
           stat = "identity",
           width = 0.6,
           size = 0.25,
           alpha = 0.4)+

  
  ## errorbar
  geom_errorbar(data = unique(valcc_OD[,c("monoculture_or_coculture","Mean_OD_48h",
                                           "SD_OD_48h",
                                           "cultype","group_culture_reord")]),
                aes(
                 x = group_culture_reord,
                 y = Mean_OD_48h,
                 ymin = Mean_OD_48h - SD_OD_48h,
                 ymax = Mean_OD_48h + SD_OD_48h
                ),
                colour = "black",
                width = 0.3,
                size = 0.2)+
  
   theme_auxcocul(base_size = 10)+
   scale_fill_manual(values = cultype_colkey)+
   scale_colour_manual(values = cultype_colkey)+
   theme(axis.text.x = element_text(angle=90,size=8, hjust = 0.99),
            legend.position = "bottom")+
   labs(x="",y="OD600 at 48h ",fill="")

dev.off()

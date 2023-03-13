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
            mutate(remade_or_library = ifelse(grepl("BY4741",monoculture_or_coculture),"library",
                                                     ifelse(grepl("COLLECTION",monoculture_or_coculture),"library","remade")),
                   monoculture_or_coculture = gsub(" COLLECTION","",monoculture_or_coculture))%>%
            filter(!(group %in% c("his2Δ-met14Δ","arg2Δ-met5Δ","arg2Δ-met3Δ")))%>%
            dplyr::group_by(monoculture_or_coculture, group, remade_or_library)%>%
            dplyr::mutate(Mean_OD_48h = mean(OD600_48h,na.rm = T),
                   SD_OD_48h = sd(OD600_48h, na.rm=T),
                   SE_OD_48h = sd(OD600_48h, na.rm=T)/sqrt(length(OD600_48h)))%>%
            dplyr::ungroup()%>%
            mutate(group_culture_reord = factor(paste(group,monoculture_or_coculture,remade_or_library),
                                                levels = c("WT BY4741 library",
                                                           "WT BY4741 pHLUM library",
                                                           
                                                           "met3Δ-his2Δ met3Δ library", "met3Δ-his2Δ met3Δ remade",
                                                           "met3Δ-his2Δ his2Δ library",    "met3Δ-his2Δ his2Δ remade",
                                                           "met3Δ-his2Δ met3Δ-his2Δ library","met3Δ-his2Δ met3Δ-his2Δ remade",
                                                          
                                                           "met3Δ-met1Δ met3Δ library", "met3Δ-met1Δ met3Δ remade",
                                                           "met3Δ-met1Δ met1Δ library", "met3Δ-met1Δ met1Δ remade",
                                                           "met3Δ-met1Δ met3Δ-met1Δ library", "met3Δ-met1Δ met3Δ-met1Δ remade",
                                                           
                                                           "met14Δ-met5Δ met14Δ library","met14Δ-met5Δ met14Δ remade",
                                                           "met14Δ-met5Δ met5Δ library","met14Δ-met5Δ met5Δ remade",
                                                           "met14Δ-met5Δ met14Δ-met5Δ library","met14Δ-met5Δ met14Δ-met5Δ remade",
                                                           
                                                           "met14Δ-trp4Δ met14Δ library","met14Δ-trp4Δ met14Δ remade",
                                                           "met14Δ-trp4Δ trp4Δ library",  "met14Δ-trp4Δ trp4Δ remade",
                                                           "met14Δ-trp4Δ met14Δ-trp4Δ library","met14Δ-trp4Δ met14Δ-trp4Δ remade",
                                                           
                                                           "met1Δ-arg2Δ met1Δ library","met1Δ-arg2Δ met1Δ remade",
                                                           "met1Δ-arg2Δ arg2Δ library","met1Δ-arg2Δ arg2Δ remade",
                                                           "met1Δ-arg2Δ met1Δ-arg2Δ library","met1Δ-arg2Δ met1Δ-arg2Δ remade",
                                                           
                                                           "lys12Δ-trp4Δ lys12Δ library","lys12Δ-trp4Δ lys12Δ remade",
                                                           "lys12Δ-trp4Δ trp4Δ library","lys12Δ-trp4Δ trp4Δ remade",
                                                           "lys12Δ-trp4Δ lys12Δ-trp4Δ library","lys12Δ-trp4Δ lys12Δ-trp4Δ remade",
                                                           
                                                           "met3Δ-trp4Δ met3Δ library", "met3Δ-trp4Δ met3Δ remade",
                                                           "met3Δ-trp4Δ trp4Δ library", "met3Δ-trp4Δ trp4Δ remade",
                                                           "met3Δ-trp4Δ met3Δ-trp4Δ library","met3Δ-trp4Δ met3Δ-trp4Δ remade",
                                                           
                                                           "trp2Δ-trp4Δ trp2Δ library", "trp2Δ-trp4Δ trp2Δ remade",
                                                           "trp2Δ-trp4Δ trp4Δ library", "trp2Δ-trp4Δ trp4Δ remade",
                                                           "trp2Δ-trp4Δ trp2Δ-trp4Δ library","trp2Δ-trp4Δ trp2Δ-trp4Δ remade",
                                                      
                                                           "met14Δ-arg2Δ met14Δ library", "met14Δ-arg2Δ met14Δ remade", 
                                                           "met14Δ-arg2Δ arg2Δ library", "met14Δ-arg2Δ arg2Δ remade",
                                                           "met14Δ-arg2Δ met14Δ-arg2Δ library","met14Δ-arg2Δ met14Δ-arg2Δ remade"
                                                     
                                                    )))
           
          
            
pdf(paste0(sup_fig_dir,"/validated_cocultures_OD600_48h.pdf"),width = 10, height = 6)
ggplot()+
  ## point
  geom_point(data = valcc_OD,
             aes( x = group_culture_reord,
                  y = OD600_48h,
                  colour = cultype),
                  size = 0.75,
             position = position_jitterdodge(jitter.width = 0.4 ))+
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
                                           "SE_OD_48h",
                                           "cultype","group_culture_reord")]),
                aes(
                 x = group_culture_reord,
                 y = Mean_OD_48h,
                 ymin = Mean_OD_48h - SE_OD_48h,
                 ymax = Mean_OD_48h + SE_OD_48h
                ),
                colour = "black",
                width = 0.3,
                size = 0.2)+
   scale_y_continuous(breaks = seq(0, 1, by = 0.2))+
   theme_auxcocul(base_size = 10)+
   scale_fill_manual(values = cultype_libremad_colkey)+
   scale_colour_manual(values = cultype_libremad_colkey)+
   theme(axis.text.x = element_text(angle=90,size=8, hjust = 0.99),
            legend.position = "bottom")+
   labs(x="",y="OD600 at 48h ",fill="")

dev.off()

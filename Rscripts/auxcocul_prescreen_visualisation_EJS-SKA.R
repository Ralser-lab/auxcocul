#################
### Set Paths ### 
#################

proj_dir <- "/Users/aulakhs/Documents/Ralser Lab/Collaborations/Auxotroph co culture screen"

inp_data_dir <- paste0(proj_dir,"/Data/input/pre_growthscreen")

###############################################################
### Initialize paths, load functions and graphic parameters ###
###############################################################

source(paste0(proj_dir,"/Rscripts/auxcocul_functions_graphparams_SKA.R"))

#########################################
### Functions specific to this script ###
#########################################

#### Read tecan.OD excel files, create a df, and convert values to OD.spec

convert_tecanODfile2ODspecdf <- function(file_path) {
  
  tecan <- read.xlsx(file_path,sheetIndex = 1)
  
  r0 <- grep('A',tecan[,1],ignore.case = F)
  rf <- grep('End Time:',tecan[,1],ignore.case = F)
  
  tecan <- data.frame(tecan[r0:(rf-1),])
  colnames(tecan) <- c("row",seq(1:(ncol(tecan)-1)))
  
  spec <- reshape2::melt(tecan, id.vars = "row",variable.name = "column", value.name = "tecanOD" )%>%
          # convert to OD in spectrophotometer based on calibration curve between tecan and spec
           dplyr::mutate(tecanOD = as.numeric(tecanOD),
                        specOD = ( tecanOD +  0.08311471132883352)/0.19844904301220034 )
  
  return(spec)
}

##############################
### Read in prescreen data ###
##############################

# 384 well plate map 

plate_map_384 <- read.csv(paste0(inp_data_dir,"/prescreen_growthcurve_384well_platemap.csv"), stringsAsFactors = F)


# read in OD600 data at T0 and T18

OD600_tp_files <- dir(inp_data_dir)[grepl("T",dir(inp_data_dir))]

prescreen_SMgrowth <- vector()

for(i in 1:length(OD600_tp_files)){
  
  # read in excel file output from tecan and convert ODtecan to ODspectrophotometer
  ODspec_df <- convert_tecanODfile2ODspecdf(paste0(inp_data_dir,"/",OD600_tp_files[i]))
  
  # merge with annotation in plate map 
  ODspec_df <- merge(ODspec_df, plate_map_384, by.x = c("row","column"), by.y = c("row_384","column_384"))
  
  # add timepoint annotation
  ODspec_df$timepoint = gsub("[.]xlsx","",unlist(strsplit(OD600_tp_files[i],split = "_"))[[4]])
  
  prescreen_SMgrowth <- rbind(prescreen_SMgrowth,ODspec_df)
  
}

## fix deletion names

prescreen_SMgrowth <- prescreen_SMgrowth%>%
                      mutate(deletion = paste0(tolower(Gene),"Δ"))

###################################################
### Calculate delta OD ( OD600 T18  - OD600 T0) ###
###################################################

prescreen_SMgrowth <- prescreen_SMgrowth%>%
                        filter(ORF != "YDL227C")%>% # filter out HO strain as its slightly different from the other positive controls -- mutation in mating locus not biosynthesis
                        dplyr::select(-tecanOD)%>%
                        reshape2::dcast(ORF+deletion+row+column ~timepoint, value.var = "specOD")%>%
                        dplyr::mutate(deltaOD600_T18 = T18-T0)%>%

###################################################
### Calculate stats ( mean and standard error ) ###
###################################################

                        dplyr::group_by(ORF,deletion)%>%
                        dplyr::mutate(mean_deltaOD600 = mean(deltaOD600_T18,na.rm = T),
                                      sd_deltaOD600 = sd(deltaOD600_T18,na.rm = T),
                                      se_deltaOD600 =  sd(deltaOD600_T18,na.rm = T)/sqrt(length(deltaOD600_T18)))%>%
                        dplyr::ungroup()%>%

####################################################################################################
### Add annotation for prototrophic strains (positive controls for normal growth) and auxotrophs ###
####################################################################################################
                      
                      dplyr::mutate(strain_type = ifelse(ORF %in% c("BY4741_pHLUM","YOR202W","YLR303W"), "pos_control",
                                                         ifelse(ORF == "BLANK","neg_control","test_auxotroph")))
                    

###########################################################################################
### Calculate threshold for considering strains auxotrophic -- 20% of growth of control ###
###########################################################################################

ODthresh_for_auxotrophy <- round(0.2*mean(unique(filter(prescreen_SMgrowth, strain_type == "pos_control")[,c("ORF","mean_deltaOD600")])$mean_deltaOD600),1 )
                     
###############################################################################################################################
### Add annotation for considering a deletion mutant an auxotroph -- based on OD threshold just calculated from prototrophs ###
###############################################################################################################################          

prescreen_SMgrowth <- prescreen_SMgrowth%>%
                      dplyr::mutate(consider_auxotroph = ifelse(mean_deltaOD600 < ODthresh_for_auxotrophy,T,F),
                                    
#########################
### add colour factor ###
######################### 
                                    prescreen_result = ifelse(strain_type == "pos_control","prototroph",
                                                              ifelse(strain_type =="neg_control","blank",
                                                                     ifelse(consider_auxotroph,"suitable_aux_KO",
                                                                               "unsuitable_aux_KO"))))
                      
############
### plot ###
############

pdf(paste0(sup_fig_dir,"/prescreen_growth_of_deletion_mutants.pdf"),width = 5,height = 18)
ggplot()+
      # point plot 
      geom_point(data = prescreen_SMgrowth,
                 aes(x = reorder(deletion,mean_deltaOD600),
                     y = deltaOD600_T18,
                     colour = prescreen_result),
                 position = position_jitterdodge(jitter.width = 0.6),
                 size = 0.75
      )+
    # bar plot
    geom_bar(data = unique(prescreen_SMgrowth[,c("deletion","mean_deltaOD600","se_deltaOD600","prescreen_result")]),
             aes(x = reorder(deletion,mean_deltaOD600),
                 y = mean_deltaOD600,
                 colour = prescreen_result,
                 fill = prescreen_result),
             stat = "identity",
             width = 0.6,
             size = 0.1,
             alpha = 0.4
             )+
  
    # errorbar
    geom_errorbar(data = unique(prescreen_SMgrowth[,c("deletion","mean_deltaOD600","se_deltaOD600","prescreen_result")]),
             aes(x = reorder(deletion,mean_deltaOD600),
                 y = mean_deltaOD600,
                 ymin = mean_deltaOD600 - se_deltaOD600,
                 ymax = mean_deltaOD600 + se_deltaOD600),
                 colour = "black",
                 width = 0.3,
                 size = 0.2
                 )+
    # add horizontal line indicating deltaOD600 threshold
  geom_hline(yintercept = ODthresh_for_auxotrophy,
            colour = "red", linetype = "dashed",
            size = 0.25)+
  
  # flip coordinates as x axis labels become too small on A4 size otherwise
  coord_flip()+
  scale_colour_manual(values = prescreen_colkey)+
  scale_fill_manual(values = prescreen_colkey)+
  theme_auxcocul()+
  theme(axis.text.y = element_text(size=8, hjust = 0.99),
        legend.position = "bottom")+
  labs(x = "",
       y = "mean(delta OD600) at 18h")
dev.off()




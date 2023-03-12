library(readxl)
library(reshape2)
library(tidyr)
library(dplyr)
library(tidyverse)
library(xlsx)
library(venn)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(ggraph)
library(igraph)
library(ggforce)


##############################
### Initalize common paths ###
##############################

output_data_dir <- paste0(proj_dir,"/data/output")

main_fig_dir <- paste0(proj_dir,"/figures/main_figures")
sup_fig_dir <- paste0(proj_dir,"/figures/supplementary_figures")

dir.create(main_fig_dir,recursive = T)
dir.create(sup_fig_dir,recursive = T)
dir.create(output_data_dir,recursive = T)


##############################
### Set graphic parameters ###
##############################

## Theme ##
theme_auxcocul <- function(base_size = 25) {
  # Starts with theme_grey and then modify some parts
  theme_minimal(base_size = base_size) %+replace%
    theme(
      
      strip.background = element_rect(colour="black",fill=NA, size = 0.1),
      strip.text.x = element_text(size = 20,
                                  margin = margin(.2,0,.2,0, "cm")),
      strip.text.y = element_text(size = 20,
                                  margin = margin(.2,0,.2,0, "cm")),
      strip.switch.pad.grid = unit(0.2,"cm"),
      strip.placement = "outside",
      
      axis.text.x = element_text(size=16),
      axis.text.y = element_text(size=16,hjust=1),
      axis.ticks =  element_line(colour = "black",size=0.25), 
      axis.title.x= element_text(size=18),
      axis.title.y= element_text(size=18,angle=90,vjust=3),
      
      panel.background = element_blank(),
      panel.border = element_rect(colour="black",fill=NA, size = 0.1), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      
      plot.background = element_rect(colour="black",fill=NA, size = 0.1),
      plot.margin = unit(c(0.75,  0.75, 0.75, 0.75), "lines"),
      
      axis.line.x = element_line(color="black", size = 0.1),
      axis.line.y = element_line(color="black", size = 0.1)
    )
}

### Colour key for auxotrophy categories ###

colkey_aux_dark <- c( arginine = "#B22222",
                      `aromatic amino acids` ="#786ea8ff",
                      histidine = "#6E8B3D", 
                      `isoleucine/valine` = "#00688B",
                      lysine = "#008B8B",
                      `methionine/sulphate`= "#8B4789",
                      `methionine/siroheme`= "#8B4789",
                      tryptophan = "#EE7621",
                      uracil = "#8B1C62",
                      adenine = "#CD950C",
                      `transmembrane transport` = "#8B8B7A", 
                      `DNA damage response`="#8B8B7A",
                      `DNA translocase activity`="#8B8B7A",
                      `proteasome`="#8B8B7A",
                      `protein maturation`="#8B8B7A",
                      `unknown`="#8B8B7A",
                      `V-ATPase assembly`="#8B8B7A",
                      ribosome ="#8B8B7A",
                      ergosterol = "#8B8B7A",
                      `stress response` = "#8B8B7A",
                      ceramide = "#8B8B7A",
                      `vitamin B6`="#8B8B7A",
                      cytosine = "#8B8B7A",
                      clathrin = "#8B8B7A",
                      `mRNA splicing`="#8B8B7A",
                      other="#8B8B7A"
                      
)

colkey_aux_light <- c( arginine = "#F08080",
                       `aromatic amino acids` ="#bebadaff",
                       histidine = "#A2CD5A", 
                       `isoleucine/valine` = "#33A1C9",
                       lysine = "#79CDCD",
                       `methionine/sulphate`= "#D8BFD8",
                       `methionine/siroheme`= "#D8BFD8",
                       tryptophan = "#FFA54F",
                       uracil = "#CD6090",
                       adenine = "#ffed6fff",
                       `transmembrane transport` = "#CCCCCC", 
                       `DNA damage response`="#CCCCCC",
                       `DNA translocase activity`="#CCCCCC",
                       `proteasome`="#CCCCCC",
                       `protein maturation`="#CCCCCC",
                       `unknown`="#CCCCCC",
                       `V-ATPase assembly`="#CCCCCC",
                       ribosome ="#CCCCCC",
                       ergosterol = "#CCCCCC",
                       `stress response` = "#CCCCCC",
                       ceramide = "#CCCCCC",
                       `vitamin B6`="#CCCCCC",
                       cytosine = "#CCCCCC",
                       clathrin = "#CCCCCC",
                       `mRNA splicing`="#CCCCCC",
                       other="#CCCCCC"
                       
)


colkey_geneonto_dark <- c( `purine nucleobase biosynthesis`= "#FFD700",
                           `arginine biosynthesis` = "#B22222",
                           `ergosterol biosynthesis` = "#8B8B7A",
                           `histidine biosynthesis` = "#6E8B3D",
                           `siroheme biosynthesis` = "#8B4789",
                           `lysine biosynthesis` = "#008B8B",
                           `sulphate assimilation` = "#8B4789",
                           `methionine biosynthesis` = "#8B4789",
                           `tryptophan biosynthesis` = "#EE7621",
                           `pyrimidine nucleobase biosynthesis` = "#8B1C62",
                           `citruline biosynthesis` = "#B22222",
                           `chorismate biosynthesis` = "#786ea8ff",
                           `tyrosine biosynthesis` = "#786ea8ff",
                           `isoleucine biosynthesis` = "#00688B",
                           `ceramide biosynthesis` = "#2F4F4F",
                           
                           `ornithine transporter`= "#8B8B7A",
                           `oxidative stress response` = "#8B8B7A",
                           `cytoplasmic translation`= "#8B8B7A",
                           `ribosome structure` = "#8B8B7A",
                           `bud site selection` = "#8B8B7A",
                           `nucleobase transport` = "#8B8B7A",
                           `peroxisome associated` = "#8B8B7A",
                           `DNA translocase activity` = "#8B8B7A",
                           `DNA damage response`= "#8B8B7A",
                           `transmembrane transport`= "#8B8B7A",
                           `golgi transport`= "#8B8B7A",
                           `mRNA splicing`= "#8B8B7A",
                           `vacuolar acidification`= "#8B8B7A",
                           `proteasome assembly`= "#8B8B7A",
                           `unknown`= "#8B8B7A",
                           `iron-sulpher cluster`= "#8B8B7A"
                           
                           
)



colkey_geneonto_light <- c( `purine nucleobase biosynthesis`= "#ffed6fff",
                            `arginine biosynthesis` = "#F08080",
                            `ergosterol biosynthesis` = "#CD8162",
                            `histidine biosynthesis` = "#A2CD5A",
                            `siroheme biosynthesis` = "#D8BFD8",
                            `lysine biosynthesis` = "#79CDCD",
                            `sulphate assimilation` = "#D8BFD8",
                            `methionine biosynthesis` = "#D8BFD8",
                            `tryptophan biosynthesis` = "#EE7621",
                            `pyrimidine nucleobase biosynthesis` = "#CD6090",
                            `citruline biosynthesis` = "#F08080",
                            `chorismate biosynthesis` = "#bebadaff",
                            `tyrosine biosynthesis` = "#bebadaff",
                            `isoleucine biosynthesis` = "#33A1C9",
                            `ceramide biosynthesis` = "#668B8B",
                            
                            `ornithine transporter`= "#CCCCCC",
                            `oxidative stress response` = "#CCCCCC",
                            `cytoplasmic translation`= "#CCCCCC",
                            `ribosome structure` = "#CCCCCC",
                            `bud site selection` = "#CCCCCC",
                            `nucleobase transport` = "#CCCCCC",
                            `peroxisome associated` = "#CCCCCC",
                            `DNA translocase activity` = "#CCCCCC",
                            `DNA damage response`= "#CCCCCC",
                            `transmembrane transport`= "#CCCCCC",
                            `golgi transport`= "#CCCCCC",
                            `mRNA splicing`= "#CCCCCC",
                            `vacuolar acidification`= "#CCCCCC",
                            `proteasome assembly`= "#CCCCCC",
                            `unknown`= "#CCCCCC",
                            `iron-sulpher cluster`= "#CCCCCC"
                            
                            
)


cultype_colkey <- c(monocul1= "#807dbaff",
                    monocul2 ="#bcbddcff",
                    cocul = "#6a51a3ff",
                    `WT`= "#74c476ff")


prescreen_colkey <- c(unsuitable_aux_KO= "#EEC591",
                       suitable_aux_KO = "#6a51a3ff",
                      `prototroph`= "#74c476ff",
                       blank = "gray")


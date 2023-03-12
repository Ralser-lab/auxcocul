#################
### Set Paths ### 
#################

proj_dir <- "/Users/aulakhs/Documents/Ralser Lab/Collaborations/Auxotroph co culture screen"

inp_data_dir <- paste0(proj_dir,"/Data/input/main_growthscreen")

###############################################################
### Initialize paths, load functions and graphic parameters ###
###############################################################

source(paste0(proj_dir,"/Rscripts/auxcocul_functions_graphparams_SKA.R"))

#######################################
### Choose pvalue and FC thresholds ###
#######################################

pvthresh <- 0.05
FCthresh <- 1.5

##############################
### Read in processed data ###
##############################

setwd(inp_data_dir)

## Manually flagged auxotrophs ###

filtoutman <- c("CCS1","FUN12","PRO2","PHA2","BAS1","HO","MET17","HIS3","BY4741","ORF , Dubious (2)")

### All monocultures ####

monocul_QC <- read.csv("all_monoculture_zfactors.csv",stringsAsFactors = F)%>%
            mutate( condition = gsub("ORF [,] Uncharacterized [(]1[)]","YBR209W",condition),
                    condition = gsub("ORF [,] Uncharacterized [(]2[)]","YHL005C",condition))%>%
             filter(!condition %in% filtoutman)

allcol_passQC <- read.csv("all_QCpassed_conditions.csv",stringsAsFactors = F)%>%
          mutate( condition = gsub("ORF [,] Uncharacterized [(]1[)]","YBR209W",condition),
                  condition = gsub("ORF [,] Uncharacterized [(]2[)]","YHL005C",condition))%>%
     #    filter(!grepl("ORF-ts",ORF))%>%
         filter(!condition %in% filtoutman)%>%
         filter(! ( grepl(filtoutman[1],condition) |
                      grepl(filtoutman[2],condition) |
                      grepl(filtoutman[3],condition) |
                      grepl(filtoutman[4],condition) |
                      grepl(filtoutman[5],condition) |
                      grepl(filtoutman[6],condition) |
                      grepl(filtoutman[7],condition) |
                      grepl(filtoutman[8],condition) |
                      grepl(filtoutman[9],condition) |
                      grepl("ORF [,] Dubious [(]2[)]",condition) 
                  )
                )%>%
    ## filter out strains in some batches that grew in that batch but not in any other --> possible contamin
      
      filter(!(grepl("ORT1",condition) & batch == 6 ) &
            !(grepl("BAS1",condition) & batch == 10 ) )

cocul_passQC <- allcol_passQC%>%
                filter(grepl("[-]",condition))%>%
                pull(condition)%>%
                unique()

monocul_incoculpassQC <- data.frame(cocul = cocul_passQC)%>%
                          separate(cocul, into=c("mon1","mon2"),sep="[-]",remove=F)%>%
                          reshape2::melt(id.vars="cocul")%>%
                          pull(value)%>%
                          unique()

write.csv(monocul_incoculpassQC,paste0(output_data_dir,"/monoculutres_in_combinations_tested_QCpass.csv"),row.names=F)

stats <- read.csv("all_coculture_stats.csv",stringsAsFactors = F)%>%
          ### Filter stats for cocultures that passed QC checks
          filter(Cocul %in% cocul_passQC)%>%
          mutate(p.value.adj = p.adjust(t.test_p.value,method="BH"),
                 Significant = ifelse(p.value.adj < pvthresh &
                                      fold.change.48 > FCthresh, 1,0))%>%
        mutate( A = gsub("ORF [,] Uncharacterized [(]1[)]","YBR209W",A),
                A = gsub("ORF [,] Uncharacterized [(]2[)]","YHL005C",A),
                
                B = gsub("ORF [,] Uncharacterized [(]1[)]","YBR209W",B),
                B = gsub("ORF [,] Uncharacterized [(]2[)]","YHL005C",B),
                
                nearest.monoculture = gsub("ORF [,] Uncharacterized [(]1[)]","YBR209W",nearest.monoculture),
                nearest.monoculture = gsub("ORF [,] Uncharacterized [(]2[)]","YHL005C",nearest.monoculture),
                
                Cocul = gsub("ORF [,] Uncharacterized [(]1[)]","YBR209W",Cocul),
                Cocul = gsub("ORF [,] Uncharacterized [(]2[)]","YHL005C",Cocul),
        )%>%
        filter(A %in% monocul_incoculpassQC | B %in% monocul_incoculpassQC)

write.csv(stats,paste0(output_data_dir,"/all_coculture_stats_post_QC_for_volcanoplot.csv"),row.names = F)
## Function of ORF

orfcontext <- read.csv("monos_QCpassed_SGD_descriptions.csv",stringsAsFactors = F)%>%
              mutate(biosynthesis = ifelse(grepl("biosynthesis",gene.ontology),1,0))


#####################
### Define "hits" ###
#####################

hits <- stats %>%
        filter(Significant ==1)%>%
        pull(Cocul)

write.csv(hits,paste0(output_data_dir,"/hits_adjpvthresh_",pvthresh,"_fcthresh_",FCthresh,".csv"),
          row.names =F)

## auxotrophs present in hits
monocul_inHits <- data.frame(hits = hits)%>%
  separate(hits, into=c("mon1","mon2"),sep="[-]",remove=F)%>%
  reshape2::melt(id.vars="hits")%>%
  pull(value)%>%
  unique()
write.csv(monocul_inHits,paste0(output_data_dir,"/monocultures_in_hits.csv"),row.names=F)

top20 <- stats %>%
  filter(Significant ==1)%>%
  top_n(20,fold.change.48)%>%
  pull(Cocul)

vp_sig_pvadj<-ggplot(stats,
               aes(x=log2(fold.change.48),
                   y=-log(p.value.adj),
                   colour = factor(Significant),
                   fill= factor(Significant)))+
  geom_point(aes(size=log2(fold.change.48)),alpha=0.7)+
  geom_text_repel(data=filter(stats,Cocul %in% top20),
                  aes(x=log2(fold.change.48),
                      y=-log(p.value.adj),
                      label=Cocul),
                  size=2,
                  colour= "black",
                  max.overlaps = 11)+
  geom_hline(yintercept = -log(pvthresh),linetype="dashed")+
  geom_vline(xintercept = log2(FCthresh),linetype="dashed")+
  scale_colour_manual(values=c("gray","darkred"))+
  scale_fill_manual(values=c("gray","darkred"))+
  labs(x="Log2(Fold Change)",
       y="-Log(adjusted pvalue)")+
    theme_auxcocul(base_size = 20)+
  theme(legend.position="none")

setwd(main_fig_dir)

pdf(paste0("volcanoplot_all_tested_auxotrophpairs_adjpv_",pvthresh,"_fc_",FCthresh,".pdf"),
    width=9,height=5)
vp_sig_pvadj
dev.off()

########################################
### Heirarchichal edge bundling plot ###
########################################

setwd(inp_data_dir)

t1 <- read.csv("origin.lvl.1.csv",stringsAsFactors = F)
t2 <- read.csv("lvl.1.lvl.2.csv",stringsAsFactors = F)


### HEB adjusted PValue hits


test <- stats%>%
        filter(Significant ==1)%>%
        dplyr::select(A,B,fold.change.48)

# Filter for annotated origin and level nodes which are hits -- pv 

t2.2 <- filter(t2, to %in% c(test$A,test$B))%>%
  unique()

t1.2 <- filter(t1, to %in% t2.2$from)%>%
  unique()%>%
  arrange(to)

edges <- rbind(t1.2, t2.2)

connections <- test %>% 
  dplyr::select('A', 'B', 'fold.change.48')
colnames(connections) <- c("from", "to", 'value')

vertices <- data.frame(name = unique(c(as.character(edges$from), as.character(edges$to))), 
                       value = seq(1, length(as.character(edges$from))+1)) %>% #90, 39, 21 
  arrange(name)
vertices$group  <-  edges$from[match(vertices$name, edges$to)]


# Create 'count' column within vertices dataframe to indicate extent of connectivity
temp1 <- as.data.frame(table(connections$from))
temp2 <- as.data.frame(table(connections$to))
vertices$count=NA

for (i in 1:nrow(vertices)) {
  if (as.character(vertices$name[i]) %in% as.character(temp1$Var1)) {
    vertices$count[i] <- temp1[which(temp1 == as.character(vertices$name[i])),2]
  }else if (as.character(vertices$name[i]) %in% as.character(temp2$Var1)) {
    vertices$count[i] <- temp2[which(temp2 == as.character(vertices$name[i])),2]
  }
}


# Preparation to draw labels properly
vertices$id=NA
myleaves=which(is.na( match(vertices$name, edges$from) ))
nleaves=length(myleaves)

vertices$id[ myleaves ] = seq(1:nleaves)

# Build a network object from this dataset:
mygraph <- graph_from_data_frame(edges, vertices = vertices)

# The connection object must refer to the ids of the leaves:
from <- match(connections$from, vertices$name)
to <- match(connections$to, vertices$name)

setwd(main_fig_dir)

# plot

pdf(paste0("heb_allhits_adjpv_",pvthresh,"_fc_",FCthresh,".pdf"),width=8,height=9)

ggraph(mygraph, layout = 'partition', circular = TRUE) + 
  geom_conn_bundle(data = get_con(from = from, to = to),
                   width = 1,
                   alpha = .25, 
                   tension = .92) +
  geom_node_arc_bar(aes(fill = group,
                        filter=leaf,
                        r0=1.90,
                        r=2.10,
                        colour=NA) ,
                        alpha=1,
                        show.legend = F)+
  geom_node_point(aes(filter = leaf,
                      x = x*1.02, 
                      y = y*1.02, 
                      colour=group,
                      fill=NA,
                      size = count),
                      alpha=1) +
  geom_node_text(aes(x = x*1.17, 
                     y = y*1.17,
                     filter = leaf,
                     label=name,
                     angle=-((-node_angle(x, y)+90)%%180)+90),
                 size=3,
                 hjust="outward",
                 alpha=1,
                 show.legend = F) +
  scale_colour_manual(values= colkey_aux_light) +
  scale_fill_manual(values= colkey_aux_dark) +
  scale_size_continuous(range = c(5,15)) +
  coord_fixed() +
  theme_void() +
  theme(
        plot.margin=unit(c(0,0,0,0),"cm")) +
  expand_limits(x = c(-2.6, 2.6), y = c(-2.6, 2.6))

dev.off()

#################################################
### Distribution of auxotrophs through screen ###
#################################################

setwd(sup_fig_dir)

monocul_tested_anno <- orfcontext%>%
                       group_by(gene.ontology)%>%
                       mutate(NumAux_inGO = length(condition))%>%
                       ungroup()%>%
                       group_by(associated.product)%>%
                       mutate(NumAux_inAP = length(condition))%>%
                       ungroup()

## Monocultures in tested strains ###

tested_go<-ggplot(monocul_tested_anno,
       aes(x=reorder(gene.ontology,NumAux_inGO),
           fill = gene.ontology))+
  geom_bar(position="stack",stat="count",width = 0.6,size=0.25)+
  coord_flip()+
  theme_auxcocul(base_size = 20)+
  scale_fill_manual(values = colkey_geneonto_light)+
  theme(legend.position = "none")+
  labs(x="Number of auxotrophs",y="")

tested_ap<-ggplot(monocul_tested_anno,
                  aes(x=reorder(associated.product,NumAux_inAP),
                      colour = associated.product,
                      fill = associated.product))+
  geom_bar(position="stack",stat="count",width = 0.6,size=0.25)+
  coord_flip()+
  scale_colour_manual(values = colkey_aux_dark)+
  scale_fill_manual(values = colkey_aux_light)+
    theme_auxcocul(base_size = 20)+
  theme(legend.position = "none")+
  labs(y="Number of auxotrophs",x="")


monocul_hits_anno <- orfcontext%>%
  filter(condition %in% monocul_inHits)%>%
  group_by(gene.ontology)%>%
  mutate(NumAux_inGO = length(condition))%>%
  ungroup()%>%
  group_by(associated.product)%>%
  mutate(NumAux_inAP = length(condition))%>%
  ungroup()

## Monocultures in tested strains ###

monohits_go <- ggplot(monocul_hits_anno,
       aes(x=reorder(gene.ontology,NumAux_inGO),
           fill = gene.ontology))+
  geom_bar(position="stack",stat="count",width = 0.6,size=0.25)+
  coord_flip()+
    theme_auxcocul(base_size = 20)+
  scale_fill_manual(values = colkey_geneonto_light)+
  theme(legend.position = "none")+
  labs(y="Number of auxotrophs",x="")

monohits_ap<-ggplot(monocul_hits_anno,
                  aes(x=reorder(associated.product,NumAux_inAP),
                      fill = associated.product))+
  geom_bar(position="stack",stat="count",width = 0.6,size=0.25)+
  coord_flip()+
  scale_fill_manual(values = colkey_aux_light)+
    theme_auxcocul(base_size = 20)+
  theme(legend.position = "none")+
  labs(y="Number of auxotrophs",x="")


### Coculture distribution -- can count same auxotroph twice if present in 2 coculutres ###

monoculcounts_inhits <- data.frame(hits=hits)%>%
                        separate(hits, into =c("mono1","mono2"),sep="[-]",remove=F)%>%
                        reshape2::melt(id.vars ="hits")%>%
                        group_by(value)%>%
                        mutate(Num_Cocul = length(hits))%>%
                        ungroup()

monoculcounts_inhits <- merge(monoculcounts_inhits,orfcontext,by.x="value",by.y="condition",all.x = T)

monoculcounts_inhits_anno <- monoculcounts_inhits%>%
                              group_by(gene.ontology)%>%
                              mutate(NumAux_inGO = length(value))%>%
                              ungroup()%>%
                              group_by(associated.product)%>%
                              mutate(NumAux_inAP = length(value))%>%
                              ungroup()
                             
monohitcounts_go <- ggplot(monoculcounts_inhits_anno,
       aes(x=reorder(gene.ontology,NumAux_inGO),
           fill = gene.ontology))+
  geom_bar(position="stack",stat="count",width = 0.6,size=0.25)+
  coord_flip()+
  scale_fill_manual(values = colkey_geneonto_light)+
    theme_auxcocul(base_size = 20)+
  theme(legend.position = "none")+
  labs(y="Number of synergistic interactions",x="")

monohitcounts_ap<-ggplot(monoculcounts_inhits_anno,
                aes(x=reorder(associated.product,NumAux_inAP),
                    fill = associated.product))+
  geom_bar(position="stack",stat="count",width = 0.6,size=0.25)+
  coord_flip()+
  scale_fill_manual(values = colkey_aux_light)+
    theme_auxcocul(base_size = 20)+
  theme(legend.position = "none")+
  labs(y="Number of synergistic interactions",x="")

pdf("auxotroph_counts_at_screen_stages_by_GOanno_and_auxtype.pdf",width=20,height=15)

tested_ap + monohits_ap + monohitcounts_ap 

tested_go + monohits_go + monohitcounts_go
dev.off()

#################################
### Fractional Auxotroph type ###
#################################

forfrac_monotested <- monocul_tested_anno%>%
  mutate(AuxType = ifelse(associated.product %in%
                            c("methionine/sulphate","histidine",
                              "tryptophan","arginine","adenine",
                              "aromatic amino acids","lysine",
                              "methionine/siroheme","uracil",
                              "isoleucine/valine"),
                          associated.product,"other"))%>%
  group_by(AuxType)%>%
  summarize(Num_Bio = length(condition))%>%
  ungroup()%>%
  mutate(type ="Monocultures tested")

forfrac_monotested$Frac = forfrac_monotested$Num_Bio/sum(forfrac_monotested$Num_Bio)

forfrac_monohits <- monocul_hits_anno%>%
  mutate(AuxType = ifelse(associated.product %in%
                            c("methionine/sulphate","histidine",
                              "tryptophan","arginine","adenine",
                              "aromatic amino acids","lysine",
                              "methionine/siroheme","uracil",
                              "isoleucine/valine"),
                          associated.product, "other"))%>%
  group_by(AuxType)%>%
  summarize(Num_Bio = length(condition))%>%
  ungroup()%>%
  mutate(type ="Unique monocultures in hits")

forfrac_monohits$Frac = forfrac_monohits$Num_Bio/sum(forfrac_monohits$Num_Bio)

forfrac_monohitcounts <- monoculcounts_inhits_anno%>%
  mutate(AuxType = ifelse(associated.product %in%
                            c("methionine/sulphate","histidine",
                              "tryptophan","arginine","adenine",
                              "aromatic amino acids","lysine",
                              "methionine/siroheme","uracil",
                              "isoleucine/valine"),
                          associated.product,"other"))%>%
  group_by(AuxType)%>%
  summarize(Num_AuxType = length(value))%>%
  ungroup()%>%
  mutate(type = "Monoculture counts in hits")

forfrac_monohitcounts$Frac = forfrac_monohitcounts$Num_AuxType/sum(forfrac_monohitcounts$Num_AuxType)

allfrac <- rbind(forfrac_monotested[,c("AuxType","type","Frac")],forfrac_monohits[,c("AuxType","type","Frac")],forfrac_monohitcounts[,c("AuxType","type","Frac")])

allfrac$type <- factor(allfrac$type,levels = c("Monocultures tested",
                                               "Unique monocultures in hits",
                                               "Monoculture counts in hits"
))

allfrac$AuxType <- factor(allfrac$AuxType,levels = c("methionine/sulphate","methionine/siroheme",
                                                     "histidine", "tryptophan","arginine","adenine",
                                                     "aromatic amino acids","lysine","uracil",
                                                     "isoleucine/valine","other"))
write.csv(allfrac,paste0(output_data_dir,"/auxotrophy_type_across_screen.csv"),row.names=T)

setwd(main_fig_dir)
pdf("auxotroph_type_across_screen.pdf",width=8,height=10)
ggplot(allfrac,
       aes( x = type,
            y = Frac,
            colour=factor(AuxType),
            fill=factor(AuxType)))+
  geom_bar(position="stack",stat="identity",width = 0.6,size=0.25)+
  scale_colour_manual(values=colkey_aux_dark)+
  scale_fill_manual(values=colkey_aux_light)+
    theme_auxcocul(base_size = 20)+
  theme(axis.title = element_blank(),
        legend.position="bottom",
        axis.text.x = element_text(angle=90))
dev.off()



###################################
### Bar plots Top 25 cocultures ###
###################################

## Get OD values for each group of monoculture-coculture in top25
               
forBP <- filter(allcol_passQC,condition %in% top20 |
                       condition %in% unique(unlist(strsplit(top20,"[-]")))
                     )%>%
              dplyr::select(`OD600..48h.`,condition,batch,group)%>%
              group_by(condition,batch)%>%
              mutate(Median_OD600_48hr = median(`OD600..48h.`,na.rm=T),
                        MAD_OD600_48hr = mad(`OD600..48h.`,na.rm=T),
                        SE_OD600_48hr = MAD_OD600_48hr/sqrt(length(`OD600..48h.`)))%>%
              ungroup()%>%
              data.frame()

forBP_monoco <- vector()         

for(i in 1:length(top20)){
  
  b = unique(filter(forBP,condition == top20[i])$batch)

  df= filter(forBP, batch == b )
  
  cocul <- filter(df,condition == top20[i])[,c("condition","Median_OD600_48hr","MAD_OD600_48hr","SE_OD600_48hr","OD600..48h.")]
  cocul$cultype = "cocul"
  
  monocul1 <- filter(df,condition == unlist(
    strsplit(top20[i],"[-]"))[1])[,c("condition","Median_OD600_48hr","MAD_OD600_48hr","SE_OD600_48hr","OD600..48h.")]
  monocul1$cultype = "monocul1"
  
  monocul2 <- filter(df,condition == unlist(
    strsplit(top20[i],"[-]"))[2])[,c("condition","Median_OD600_48hr","MAD_OD600_48hr","SE_OD600_48hr","OD600..48h.")]
  monocul2$cultype = "monocul2"
  
  condgroup <- rbind(cocul,monocul1,monocul2)
  condgroup$Conditional_Group = top20[i]
  
  forBP_monoco <- rbind(forBP_monoco,condgroup)
}

## Add wildtype growth ###

setwd(inp_data_dir)

wt <- read.csv("all_QCpassed_conditions.csv",stringsAsFactors = F)%>%
      filter(condition %in% c("HO","MET17","HIS3","BY4741"))%>%
      mutate(condition ="Wild Type")%>%
      group_by(condition)%>%
      mutate(Median_OD600_48hr = median(`OD600..48h.`,na.rm=T),
                 MAD_OD600_48hr = mad(`OD600..48h.`,na.rm=T),
                 SE_OD600_48hr = MAD_OD600_48hr/sqrt(length(`OD600..48h.`)))%>%
      ungroup()%>%
      dplyr::select(condition,Median_OD600_48hr,MAD_OD600_48hr,SE_OD600_48hr,`OD600..48h.`)%>%
      as.data.frame()

wt$cultype = "WT"
wt$Conditional_Group = "Wild Type Control"

forBP_monoco <- rbind(forBP_monoco,wt)

#####################
### Make Bar Plot ###
#####################

order_bp_forbar <- forBP_monoco%>%
            dplyr::select(condition,Median_OD600_48hr,SE_OD600_48hr,cultype,Conditional_Group)%>%
            unique()%>%
            group_by(Conditional_Group)%>%
            mutate(OD_cocul = max(Median_OD600_48hr))%>%
            arrange(-OD_cocul)%>%
            ungroup()%>%
            mutate(x=1,
                   ord = cumsum(x))

order_bp_forpoint <- merge(forBP_monoco,order_bp_forbar[,c("Conditional_Group","condition","ord")],
                           by = c("condition","Conditional_Group"))
                      

setwd(main_fig_dir)

pdf("top20_cocultures_bar_plot.pdf",width=10,height=6)
ggplot()+
  geom_point(data = filter(unique(order_bp_forpoint[,c("condition","ord","OD600..48h.","cultype","Conditional_Group")]),OD600..48h. > 0.05 ),
             aes(x = reorder(paste(Conditional_Group,condition),ord),
                 y = OD600..48h.,
                 colour = cultype),
             size = 1, alpha = 0.8,
             position = position_jitterdodge(jitter.width = 0.6))+
  
  geom_bar(data = unique(order_bp_forbar[,c("condition","ord","Median_OD600_48hr","cultype","Conditional_Group")]),
           aes(x = reorder(paste(Conditional_Group,condition),ord),
               y = Median_OD600_48hr,
               colour = cultype,
               fill = cultype),
           alpha= 0.4, 
           stat = "identity", 
           width = 0.6,
           size = 0.25)+
  
  geom_errorbar( data = unique(order_bp_forbar[,c("condition","ord","Median_OD600_48hr","SE_OD600_48hr","cultype","Conditional_Group")]),
                 aes(x = reorder(paste(Conditional_Group,condition),ord),
                     y = Median_OD600_48hr,
                     ymin = Median_OD600_48hr-SE_OD600_48hr, 
                     ymax = Median_OD600_48hr+SE_OD600_48hr),
                 width=0.7,
                 size=0.5)+
  theme_auxcocul(base_size = 10)+
  scale_fill_manual(values = cultype_colkey)+
  scale_colour_manual(values = cultype_colkey)+
  theme(axis.text.x = element_text(angle=90,size=8, hjust = 0.99),
        legend.position = "bottom")+
  labs(x="",y="OD600 at 48h ",fill="")

dev.off()

#################
## Pie Charts ###
#################

setwd(inp_data_dir)

monoano_manual<-  read.csv("all_monos_w_dirindir_hit_anno.csv",
                           stringsAsFactors = F)       

monopie<-monoano_manual %>%
  group_by(AA_or_Nucl_biosynth_manual_anno)%>%
  summarize(frac=100*length(condition)/62)%>%
  ungroup()%>%
  mutate(label_loc = 100-cumsum(frac) + frac/2)%>%
ggplot(
       aes(x=1,
           y=frac,
           fill = AA_or_Nucl_biosynth_manual_anno)
)+
  geom_bar(position="stack",
           stat="identity",
           colour="white",
           size=0.25)+
  coord_polar("y")+
  geom_label(aes(x=1.5,
                 y=label_loc,
                 label=paste(round(frac,2),"%") ),
             colour="black",
             fill="white")+
  scale_fill_viridis_d(option="A",begin=0.15,end=0.85)+
    theme_auxcocul(base_size = 20)+
  theme(legend.position = "bottom",
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.x  = element_blank(),
        axis.ticks.y  = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank())+
  labs(title = "monocul",
    fill = "Function in amino/nucleic acid biosynthesis")




monohitpie<-monoano_manual %>%
  filter(condition %in% monocul_inHits)%>%
  group_by(AA_or_Nucl_biosynth_manual_anno)%>%
  summarize(frac=100*length(condition)/length(monocul_inHits))%>%
  ungroup()%>%
  mutate(label_loc = 100-cumsum(frac) + frac/2)%>%
  ggplot(
    aes(x=1,
        y=frac,
        fill = AA_or_Nucl_biosynth_manual_anno)
  )+
  geom_bar(position="stack",
           stat="identity",
           colour="white",
           size=0.5)+
  coord_polar("y")+
  geom_label(aes(x=1.5,
                 y=label_loc,
                 label=paste(round(frac,2),"%") ),
             colour="black",
             fill="white")+
  scale_fill_viridis_d(option="A",begin=0.15,end=0.85)+
    theme_auxcocul(base_size = 20)+
  theme(legend.position = "bottom",
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.x  = element_blank(),
        axis.ticks.y  = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank())+
  labs(title = "monocul in hits",
       fill = "")


## hits ##

hits_BS_manual_anno <- read.csv("hits_biosynth_manual_anno.csv",stringsAsFactors = F)


cocul_direct<-hits_BS_manual_anno%>%
  group_by(Number_direct)%>%
  summarize(frac= 100*(length(Coculture) /49))%>%
  mutate(label_loc = 100-cumsum(frac) + frac/2)%>%
  ggplot( aes(
             x=1,
             y=frac,
             fill =factor(Number_direct)
       ))+
  geom_bar(position="stack",
           stat="identity",
           colour="white",
           size=0.5)+
  coord_polar("y")+
  geom_label(aes(
    y=label_loc,
    label=paste(round(frac,2),"%") ),
    colour="black",
    fill="white")+
  scale_fill_viridis_d(option="D",begin=0,end=0.7,direction=-1)+
    theme_auxcocul(base_size = 20)+
  theme(legend.position = "bottom",
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.x  = element_blank(),
        axis.ticks.y  = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank())+
  labs( title = "cocul direct",
    fill = "")

cocul_indirect<-hits_BS_manual_anno%>%
  group_by(Number_direct_or_indirect)%>%
  summarize(frac= 100*(length(Coculture) /49))%>%
  mutate(label_loc = 100-cumsum(frac) + frac/2)%>%
  ggplot( aes(
    x=1,
    y=frac,
    fill =factor(Number_direct_or_indirect)
  ))+
  geom_bar(position="stack",
           stat="identity",
           colour="white",
           size=0.5)+
  coord_polar("y")+
  geom_label(aes(x=1.5,
    y=label_loc,
    label=paste(round(frac,2),"%") ),
    colour="black",
    fill="white")+
  scale_fill_viridis_d(option="D",begin=0,end=0.7,direction=-1)+
    theme_auxcocul(base_size = 20)+
  theme(legend.position = "bottom",
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.x  = element_blank(),
        axis.ticks.y  = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank())+
  labs( title = "cocul indirect",
       fill = "")

setwd(main_fig_dir)
pdf("direct_indirect_biosynthesis_function_pie_chart.pdf",width=15,height=7)
monopie+monohitpie+cocul_direct+cocul_indirect+plot_layout(ncol=4)
dev.off()


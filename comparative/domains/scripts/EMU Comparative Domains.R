library(tidyverse)
library(fmsb)
library(grid)
library(gridExtra)
library(UpSetR)
library(ggforestplot)
library(plotly)
library(viridis)

####MEROPS analysis####
mer.key <- read.delim("../../Comparative_pipeline/lib/merops_lib.families.tab", header=FALSE)

mer.key2 = mer.key %>%
  rename(`MERNUM`="V1", `Family`="V2", `Accession`="V3") %>%
  mutate(Type=str_sub(Family, 1, 1), Category=str_sub(Family, 2), Subfamily=str_extract_all(Category, "[aA-zZ]+", )) %>%
  mutate(Subfamily=modify_if(Subfamily, ~ length(.) == 0, 
                             ~ "None"), Strict_family=ifelse(Subfamily=="None",
                                                             Category, str_sub(Category, 1, -2L))) %>%
  mutate(Strict_family=str_c(Type, Strict_family))

chosen=c("Entomophthora_muscae_UCB.v3", "Entomophaga_maimaiga_ARSEF_7190.v1", "Strongwellsea_castrans_DrSc_Trinity", "Zoophthora_radicans_ATCC_208865", "Pandora_formicae_Trinity", "Conidiobolus_thromboides__ARSEF_4968", "Conidiobolus_coronatus_NRRL_28638.Conco1.v1")

mer.files=data.frame(names=str_replace(list.files("MEROPS/"), ".blasttab", "")) %>%
  filter(names %in% chosen) %>%
  mutate(loc=paste("MEROPS/", names, ".blasttab", sep="")) %>%
  mutate(source=as.character(1:length(loc)), Genome=as.factor(c("CCO", "CTH", "EMA", "EMU", "PFA", "SCA", "ZRA")))

merops.raw=mer.files$loc %>% map_dfr(read.delim, .id="source", header=F) %>%
  left_join(mer.files %>% select(source, Genome)) %>%
  select(-source)

merops=merops.raw %>%
  select(V1, V2, Genome) %>%
  rename(Name=`V1`, MERNUM=`V2`) %>%
  separate(Name, into=c("Strain", "Accession"), sep="\\|")

merops.counts=merops %>%
  group_by(Genome, MERNUM) %>%
  summarize(Count=n_distinct(Accession)) %>%
  group_by(Genome) %>%
  mutate(n=sum(Count)) %>%
  group_by(MERNUM) %>%
  mutate(n_Genomes=length(unique(Genome)), Genomes=toString(unique(Genome))) %>%
  left_join(mer.key2)

merops.composition=merops.counts %>%
  group_by(Genome) %>%
  summarize(n=sum(n_Genomes), MEROPSs=list(unique(MERNUM)), Families=list(unique(Family)))

merops.unique=merops.counts %>%
  filter(n_Genomes==1) %>%
  group_by(Genome) %>%
  summarize(n=n_distinct(MERNUM), MEROPS=list(unique(MERNUM)), Families=list(unique(Family)))

merops.family=merops.counts %>%
  group_by()

merops.dat=merops.counts %>%
  filter(n_Genomes>1) %>%
  select(-n_Genomes)

mer.phylo=levels(as.factor(merops.dat$Genome))[c(1, 2, 5, 6, 7, 3, 4)]

upset.dat=merops.composition$MEROPSs
names(upset.dat)=merops.composition$Genome

mer.uplt=upset(fromList(upset.dat), sets=mer.phylo, mb.ratio = c(0.55, 0.45), order.by = "freq", keep.order = TRUE)

mer.uplt

grid.text("MEROPS UpSet Plot", x = 0.65, y=0.95, gp=gpar(fontsize=20))

#MEROPS enrichment
enrich = function(dat, var='DOMAIN', group="Genome", thresh=0.01){
  out=data.frame()
  vars=unique(dat[[var]])
  len=length(vars)
  for(i in 1:len){
    sub=vars[i]
    #print(sub)
    tmp=dat[dat[[var]]==sub,]
    Ct=tmp$Count
    num=tmp$n
    names(Ct)=tmp[[group]]
    test=broom::tidy(pairwise.fisher.test(Ct, num, p.adjust.method="bonferroni"))
    if(dim(test)[1]>0){
      test2=data.frame(test) %>%
        filter(p.value<thresh) %>%
        mutate(Domain=sub)}else(test2=NA)
    out=rbind(out, test2)
  }
  return(out)
} #Enrichment function

merops.res=enrich(merops.dat, var="MERNUM", group="Genome")

merops.res2 = merops.res %>%
  group_by(group1, group2) %>%
  summarize(Domains=toString(Domain)) %>%
  pivot_wider(id_cols=group1, names_from=group2, values_from=Domains)

merops.res3=merops.res2[match(rev(mer.phylo), merops.res2$group1),]

res4=merops.res3[,match(c("group1", intersect(rev(mer.phylo), colnames(merops.res3))), colnames(merops.res3)),]

merops.res.stat=data.frame(Genome=c(merops.res$group1, merops.res$group2), Domain=c(merops.res$Domain, merops.res$Domain)) %>%
  group_by(Genome, Domain) %>%
  summarize(sig_n=length(Domain))

merops.meta=merops.dat %>%
  select(MERNUM, Genome, Count) %>%
  rename(Domain="MERNUM") %>%
  group_by(Domain) %>%
  mutate(avg=median(Count), Fold=ifelse((Count+avg)>0, Count/avg, 0),
         Direction=ifelse(Fold>1 | Fold==Inf, "Up", ifelse(Fold==1 | Fold==0, ifelse(Count<avg, "Down", "Equal"), "Down"))) %>%
  mutate(Fold=ifelse(Direction=="Down", 1/Fold, Fold))

merops.res.stat.meta=left_join(merops.res.stat, merops.meta)

merops.res.stat.meta$Genome=as.factor(merops.res.stat.meta$Genome)

#MEROPS clustering
merops.dist=merops.res.stat.meta %>%
  select(Domain, Genome, sig_n) %>%
  pivot_wider(id_cols="Domain", values_from=sig_n, names_from=Genome) %>% 
  mutate(
    across(where(is.numeric), ~replace_na(.x, 0))
  )

merops.dist2=as.matrix(merops.dist[2:length(merops.dist)])
names(merops.dist2)=merops.dist$Domain

merops.hclust=hclust(d = dist(merops.dist2))

merops.hclust2 <- data.frame(Domain=names(merops.dist2),
                             cluster=cutree(merops.hclust, length(merops.dist$Domain)))

merops.clust=left_join(merops.res.stat.meta, merops.hclust2)

merops.clust$Genome=as.factor(merops.clust$Genome)
merops.clust$Genome=factor(merops.clust$Genome, levels=levels(merops.clust$Genome)[match(mer.phylo, levels(merops.clust$Genome))])

merops.list = unique(subset(merops.res.stat.meta, merops.res.stat.meta$Genome=="EMU")$Domain)

Emus.merops.stat=merops.clust %>%
  filter(Domain %in% merops.list) %>%
  filter(Genome=="EMU") %>%
  ungroup() %>%
  select(Domain, Direction) %>%
  distinct() %>%
  rename("Emus.status"="Direction") %>%
  mutate(Emus.status=as.factor(Emus.status))

Emus.merops.clust=merops.clust %>%
  left_join(Emus.merops.stat) %>%
  mutate(Emus.status=as.character(Emus.status)) %>%
  mutate(Emus.status=replace_na(Emus.status, "Missing")) %>%
  filter(Emus.status!="Equal" & Emus.status!="Missing") %>%
  rename(`MERNUM`="Domain") %>%
  left_join(mer.key2) %>%
  rename(`Domain`="MERNUM") %>% 
  mutate(Direction=as.factor(Direction)) %>%
  mutate(Direction=factor(Direction, levels=levels(Direction)[c(3, 2, 1)]))

#MEROPS visualization
theme = theme_bw()+theme(text = element_text(size=15), axis.title.x = element_text(size=20), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15), title = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=15), strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid"))

mer.dwn=ggplot(subset(Emus.merops.clust, Emus.status=="Down"), aes(y=reorder(Domain, cluster), x=Genome, color=Direction, size=as.factor(sig_n)))+geom_point(alpha=0.8)+
  theme+geom_stripes(inherit.aes=F, aes(y=reorder(Domain, cluster)), odd = "#33333333", even = "#00000000")+theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0.95))+
  theme(axis.text.y=element_text(size=15), axis.text.x=element_text(size=15))+theme(legend.position="right")+ylab("MERNUM")+facet_grid(~Type, scales="free_x", space="free_x")+scale_color_manual(values=c("#127852", "#FF0000"))+
  labs(color="Direction\nvs. median", size="Significant\nDifferences")+coord_flip()+guides(color = guide_legend(override.aes = list(size=3), order=1))

mer.dwn

mer.up=ggplot(subset(Emus.merops.clust, Emus.status=="Up"), aes(y=reorder(Domain, cluster), x=Genome, color=Direction, size=sig_n))+geom_point(alpha=0.8)+
  theme+geom_stripes(inherit.aes=F, aes(y=reorder(Domain, cluster)), odd = "#33333333", even = "#00000000")+theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0.95))+
  theme(axis.text.y=element_text(size=15), axis.text.x=element_text(size=10))+theme(legend.position="right")+ylab("MERNUM")+facet_grid(~Type, scales="free_x", space="free_x")+scale_color_manual(values=c("#127852", "#8B8588", "#FF0000"))+labs(color="Direction\nvs. median", size="Significant\nDifferences")+coord_flip()+guides(color = guide_legend(override.aes = list(size=3), order=1))

mer.up

mer.dwn2=ggplot(subset(Emus.merops.clust, Emus.status=="Down" & Genome=="EMU"), aes(x=Count, y=reorder(Domain, cluster), fill=Fold))+geom_col()+
  theme+geom_stripes(odd = "#33333333", even = "#00000000")+
  theme(axis.text.y=element_text(size=15), axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x = element_blank())+theme(legend.position="right")+xlab("Count for EMU only")+facet_grid(~Type, scales="free_x", space="free_x")+coord_flip()+labs(fill="Fold vs.\nmedian")+scale_fill_viridis_c()+ggtitle("MEROPS Down in EMU compared to median")+theme(plot.title = element_text(hjust = 0.5))

mer.dwn2

mer.up2=ggplot(subset(Emus.merops.clust, Emus.status=="Up" & Genome=="EMU"), aes(x=Count, y=reorder(Domain, cluster), fill=Fold))+geom_col()+
  theme+geom_stripes(odd = "#33333333", even = "#00000000")+
  theme(axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x = element_blank())+theme(legend.position="right")+xlab("Count for EMU only")+facet_grid(~Type, scales="free_x", space="free_x")+coord_flip()+labs(fill="Fold vs.\nmedian")+scale_fill_viridis_c()+ggtitle("MEROPS Up in EMU compared to median")+theme(plot.title = element_text(hjust = 0.5))

mer.up2

Emus.merops.missing=merops.clust %>%
  left_join(Emus.merops.stat) %>%
  mutate(Emus.status=as.character(Emus.status)) %>%
  mutate(Emus.status=replace_na(Emus.status, "Missing from EMU")) %>%
  filter(Emus.status=="Missing from EMU") %>%
  rename("MERNUM"="Domain") %>%
  left_join(mer.key2) %>%
  rename("Domain"="MERNUM") %>%
  mutate(Direction=as.factor(Direction)) %>%
  mutate(Direction=factor(Direction, levels=levels(Direction)[c(2, 3, 1)]))

mer.plt3=ggplot(Emus.merops.missing, aes(y=reorder(Domain, rev(cluster)), x=Genome, color=Direction, size=sig_n))+geom_point(alpha=0.8)+
  theme+theme(axis.text.x = element_text(size=15, angle = 90, vjust=0.5, hjust=0.95))+theme(axis.text.y=element_text(size=15))+theme(legend.position="right")+xlab("")+facet_grid(rows="Type", scales="free", space="free")+scale_color_manual(values=c("#127852", "#8B8588", "#FF0000"))+labs(size="Significant\nDifferences")+geom_stripes(inherit.aes=F, aes(y=reorder(Domain, rev(cluster))), odd ="#33333333", even = "#00000000")+scale_x_discrete(limits = rev(mer.phylo))+ylab("MERNUM")+ggtitle("Missing from EMU")+theme(plot.title = element_text(hjust = 0.5))+guides(color = guide_legend(override.aes = list(size=3), order=1))+labs(color="Direction vs.\nmedian")

mer.plt3

Emus.merops.unique=merops.counts %>%
  filter(MERNUM %in% unlist(merops.unique[merops.unique$Genome=="EMU",]$MEROPS))  %>%
  left_join(mer.key2) %>%
  rename("Domain"="MERNUM")

mer.plt4=ggplot(Emus.merops.unique, aes(x=Domain, y=Count, fill=Type))+geom_col()+scale_fill_viridis_d()+
  theme+theme(axis.text.x = element_text(angle = 90, vjust=0.5))+xlab("MERNUM")+scale_x_discrete(limits = rev(Emus.merops.unique$Domain))+coord_flip()+ggtitle("Unique to EMU")+theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")+labs(fill="Peptidase\ntype")

mer.plt4

mer.dwn.1=ggplotGrob(mer.dwn)
mer.dwn.2=ggplotGrob(mer.dwn2)

maxWidth <- unit.pmax(mer.dwn.1$widths, mer.dwn.2$widths)

mer.dwn.1$widths= maxWidth
mer.dwn.2$widths= maxWidth

layout1 <- rbind(c(2),
                 c(1))

grid.arrange(mer.dwn.1, mer.dwn.2, layout_matrix=layout1, heights=c(0.4, 0.6))

mer.up.1=ggplotGrob(mer.up)
mer.up.2=ggplotGrob(mer.up2)

maxWidth <- unit.pmax(mer.up.1$widths, mer.up.2$widths)

mer.up.1$widths= maxWidth
mer.up.2$widths= maxWidth

grid.arrange(mer.up.1, mer.up.2, layout_matrix=layout1, heights=c(0.4, 0.6))

mer3=ggplotGrob(mer.plt3)
mer4=ggplotGrob(mer.plt4)

layout2 <- rbind(c(1,2))

grid.arrange(mer4, mer3, layout_matrix=layout2, widths=c(0.4, 0.6))

####CAZY analysis####
cazy.files=data.frame(names=str_replace(list.files("CAZY/"), ".tsv", "")) %>%
  filter(names %in% chosen) %>%
  mutate(loc=paste("CAZY/", names, ".tsv", sep="")) %>%
  mutate(source=as.character(1:length(loc)), Genome=as.factor(c("CCO", "CTH", "EMA", "EMU", "PFA", "SCA", "ZRA")))

cazy.raw=cazy.files$loc %>% map_dfr(read.delim, .id="source", header=F) %>%
  left_join(cazy.files %>% select(source, Genome)) %>%
  select(-source)

cazy=cazy.raw %>%
  select(V3, V1, Genome) %>%
  rename(Name=`V3`, CAZY=`V1`) %>%
  separate(Name, into=c("Strain", "Accession"), sep="\\|") %>%
  mutate(CAZY=str_replace_all(CAZY, ".hmm", ""))

cazy.counts=cazy %>%
  group_by(Genome, CAZY) %>%
  summarize(Count=n_distinct(Accession)) %>%
  group_by(Genome) %>%
  mutate(n=sum(Count)) %>%
  group_by(CAZY) %>%
  mutate(n_Genomes=length(unique(Genome)), Genomes=toString(unique(Genome)))

cazy.composition=cazy.counts %>%
  group_by(Genome) %>%
  summarize(n=sum(n_Genomes), CAZYs=list(unique(CAZY)))

cazy.unique=cazy.counts %>%
  filter(n_Genomes==1) %>%
  group_by(Genome) %>%
  summarize(n=n_distinct(CAZY), cazy=list(unique(CAZY)))

cazy.dat=cazy.counts %>%
  filter(n_Genomes>1) %>%
  select(-n_Genomes) 

cazy.phylo=levels(as.factor(cazy.dat$Genome))[c(1, 2, 5, 6, 7, 3, 4)]

upset.dat=cazy.composition$CAZYs
names(upset.dat)=cazy.composition$Genome

cazy.uplt=upset(fromList(upset.dat), sets=cazy.phylo, mb.ratio = c(0.55, 0.45), order.by = "freq", keep.order = TRUE)

cazy.uplt
grid.text("CAZY UpSet Plot",x = 0.65, y=0.95, gp=gpar(fontsize=20))

#CAZY enrichment
cazy.res=enrich(cazy.dat, var="CAZY", group="Genome")

cazy.res2 = cazy.res %>%
  group_by(group1, group2) %>%
  summarize(Domains=toString(Domain)) %>%
  pivot_wider(id_cols=group1, names_from=group2, values_from=Domains)

cazy.res3=cazy.res2[match(rev(cazy.phylo), cazy.res2$group1),]

res4=cazy.res3[,match(c("group1", intersect(rev(cazy.phylo), colnames(cazy.res3))), colnames(cazy.res3)),]

cazy.res.stat=data.frame(Genome=c(cazy.res$group1, cazy.res$group2), Domain=c(cazy.res$Domain, cazy.res$Domain)) %>%
  group_by(Genome, Domain) %>%
  summarize(sig_n=length(Domain))

cazy.meta=cazy.dat %>%
  select(CAZY, Genome, Count) %>%
  rename(Domain="CAZY") %>%
  group_by(Domain) %>%
  mutate(avg=median(Count), Fold=ifelse((Count+avg)>0, Count/avg, 0),
         Direction=ifelse(Fold>1 | Fold==Inf, "Up", ifelse(Fold==1 | Fold==0, ifelse(Count<avg, "Down", "Equal"), "Down"))) %>%
  mutate(Fold=ifelse(Direction=="Down", 1/Fold, Fold))

cazy.res.stat.meta=inner_join(cazy.res.stat, cazy.meta)

cazy.res.stat.meta$Genome=as.factor(cazy.res.stat.meta$Genome)

#CAZY clustering
cazy.dist=cazy.res.stat.meta %>%
  select(Domain, Genome, sig_n) %>%
  pivot_wider(id_cols="Domain", values_from=sig_n, names_from=Genome) %>% mutate(
    across(where(is.numeric), ~replace_na(.x, 0))
  )

cazy.dist2=as.matrix(cazy.dist[2:length(cazy.dist)])
names(cazy.dist2)=cazy.dist$Domain

cazy.hclust=hclust(d = dist(cazy.dist2))

cazy.hclust2 <- data.frame(Domain=names(cazy.dist2),
                           cluster=cutree(cazy.hclust, length(cazy.dist$Domain)))

cazy.clust=left_join(cazy.res.stat.meta, cazy.hclust2)

cazy.clust$Genome=as.factor(cazy.clust$Genome)
cazy.clust$Genome=factor(cazy.clust$Genome, levels=levels(cazy.clust$Genome)[match(cazy.phylo, levels(cazy.clust$Genome))])

cazy.list = unique(subset(cazy.res.stat.meta, cazy.res.stat.meta$Genome=="EMU")$Domain)

Emus.cazy.stat=cazy.clust %>%
  filter(Domain %in% cazy.list) %>%
  filter(Genome=="EMU") %>%
  ungroup() %>%
  select(Domain, Direction) %>%
  distinct() %>%
  rename("Emus.status"="Direction") %>%
  mutate(Emus.status=recode(Emus.status, "Down"="Down in\nEMU", "Up"="Up in\nEMU"))

Emus.cazy.clust=cazy.clust %>%
  left_join(Emus.cazy.stat) %>%
  mutate(Emus.status=replace_na(Emus.status, "Missing")) %>%
  filter(Emus.status!="Equal" & Emus.status!="Missing") %>%
  mutate(Direction=as.factor(Direction)) %>%
  mutate(Direction=factor(Direction, levels=levels(Direction)[c(2, 1, 3)])) 

#CAZY visualization
cazy.plt=ggplot(Emus.cazy.clust, aes(y=reorder(Domain, cluster), x=Genome, color=Direction, size=sig_n))+geom_point(alpha=0.8)+
  theme+geom_stripes(inherit.aes=F, aes(y=reorder(Domain, cluster)), odd = "#33333333", even = "#00000000")+theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0.95))+
  theme(axis.text.y=element_text(size=15), axis.text.x=element_text(size=10))+theme(legend.position="right", strip.text.y = element_text(size = 12, color = "black", face = "bold"),strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid"))+ylab("CAZY")+facet_grid(~Emus.status, scales="free_x", space="free_x")+scale_color_manual(values=c("#127852", "#8B8588", "#FF0000"))+labs(color="Direction\nvs. median", size="Significant\nDifferences")+coord_flip()+guides(colour = guide_legend(override.aes = list(size=3), order=1))

cazy.plt

cazy.plt2=ggplot(subset(Emus.cazy.clust, Genome=="EMU"), aes(x=Count, y=reorder(Domain, cluster), fill=Fold))+geom_col()+
  theme+geom_stripes(odd = "#33333333", even = "#00000000")+
  theme(axis.text.y=element_text(size=15), axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x = element_blank(), strip.text.y = element_text(size = 12, color = "black", face = "bold"),strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid"))+theme(legend.position="right")+xlab("Count for EMU only")+facet_grid(~Emus.status, scales="free_x", space="free_x")+coord_flip()+labs(fill="Fold vs.\nmedian")+scale_fill_viridis_c()

cazy.plt2

Emus.cazy.missing=cazy.clust %>%
  left_join(Emus.cazy.stat) %>%
  mutate(Emus.status=as.character(Emus.status)) %>%
  mutate(Emus.status=replace_na(Emus.status, "Missing from EMU")) %>%
  filter(Emus.status=="Missing from EMU") %>%
  mutate(Direction=as.factor(Direction)) %>%
  mutate(Direction=factor(Direction, levels=levels(Direction)[c(2, 1, 3)]))

Emus.cazy.unique=cazy.counts %>%
  filter(CAZY %in% unlist(cazy.unique[cazy.unique$Genome=="EMU",]$cazy))

cazy.plt3=ggplot(Emus.cazy.missing, aes(y=reorder(Domain, rev(cluster)), x=Genome, color=Direction, size=as.factor(sig_n)))+geom_point(alpha=0.8)+
  theme+theme(axis.text.x = element_text(angle = 90, vjust=0.5))+
  theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=10))+theme(legend.position="right")+xlab("Genome")+scale_color_manual(values=c("#127852", "#8B8588", "#FF0000"))+labs(size="Significant\nDifferences")+geom_stripes(inherit.aes=F, aes(y=reorder(Domain, rev(cluster))), odd = "#33333333", even = "#00000000")+scale_y_discrete(limits = rev(unique(Emus.cazy.missing$Domain)))+scale_x_discrete(limits = rev(cazy.phylo))+ylab("CAZY")+ggtitle("Missing from EMU")+theme(plot.title = element_text(hjust = 0.5))+guides(colour = guide_legend(override.aes = list(size=3), order=1))

cazy.plt3

cazy.plt4=ggplot(Emus.cazy.unique, aes(x=CAZY, y=Count))+geom_col(fill="black")+
  theme+theme(axis.text.x = element_text(angle = 90, vjust=0.5))+xlab("CAZY")+scale_x_discrete(limits = rev(Emus.cazy.unique$CAZY))+scale_y_continuous(breaks=c(0:2))+coord_flip()+ggtitle("Unique to EMU")+theme(plot.title = element_text(hjust = 0.5))

cazy.plt4

cazy.1=ggplotGrob(cazy.plt)
cazy.2=ggplotGrob(cazy.plt2)

maxWidth <- unit.pmax(cazy.1$widths, cazy.2$widths)

cazy.1$widths= maxWidth
cazy.2$widths= maxWidth

grid.arrange(cazy.1, cazy.2, layout_matrix=layout1, heights=c(0.4, 0.6))

cazy.3=ggplotGrob(cazy.plt3)
cazy.4=ggplotGrob(cazy.plt4)

grid.arrange(cazy.4, cazy.3, layout_matrix=layout2, widths=c(0.4, 0.6))

####Pfam analysis####
sigp <- read.table("~/Documents/GitHub/E_muscae_berkeley/comparative/domains/CAZY/Entomophthora_muscae_UCB.v3.run_dbcan/signalp.out", quote="\"", comment.char="")

#Combining noTm/TM results
sigp2 = sigp %>%
  select(V1, V12) %>%
  rename(`Name`="V1", `Method`="V12") %>%
  separate(Name, into=c("Strain", "Accession"), sep="\\|") %>%
  select(-Strain, -Method) %>%
  mutate(Prediction="Secreted") %>%
  distinct()

pfam.files=data.frame(names=str_replace(list.files("Pfam/"), ".domtbl.gz", "")) %>%
  filter(names %in% chosen) %>%
  mutate(loc=paste("Pfam/", names, ".domtbl.gz", sep="")) %>%
  mutate(source=as.character(1:length(loc)), Genome=as.factor(c("CCO", "CTH", "EMA", "EMU", "PFA", "SCA", "ZRA")))

pfam.raw=pfam.files$loc %>% map_dfr(read_table2, .id="source", comment="#", col_names=F) %>%
  left_join(pfam.files %>% select(source, Genome)) %>%
  select(-source)

pfam=pfam.raw %>%
  select(X4, X1, Genome) %>%
  rename(Name=`X1`, Pfam=`X4`) %>%
  separate(Name, into=c("Strain", "Accession"), sep="\\|")

sigp.pfam=pfam %>%
  left_join(sigp2) %>%
  mutate(Prediction=replace_na(Prediction, "Not Secreted")) %>%
  group_by(Pfam, Prediction) %>%
  summarize(Accessions_n=length(unique(Accession)), Accessions=toString(unique(Accession))) %>%
  group_by(Pfam) %>%
  mutate(Prediction_n=n_distinct(Prediction), Predictions=toString(unique(Prediction), sep="\n")) %>%
  pivot_wider(names_from="Prediction", values_from=c("Accessions_n", "Accessions")) %>%
  rename(Secreted="Accessions_n_Secreted", `Not Secreted`="Accessions_n_Not Secreted") %>%
  mutate(Secreted=replace_na(Secreted, 0),
         `Not Secreted`=replace_na(`Not Secreted`, 0)) %>%
  group_by(Pfam) %>%
  mutate(Secretion=ifelse(Prediction_n==2, "Both", Predictions)) %>%
  #select(-n, -Prediction) %>%
  rename("Domain"="Pfam")

sigp.pfam2 = sigp.pfam %>%
  select(Domain, Secretion)

pfam.counts.acc=pfam %>%
  group_by(Genome, Pfam) %>%
  summarize(Count=n_distinct(Accession)) %>%
  group_by(Genome) %>%
  mutate(n=sum(Count)) %>%
  group_by(Pfam) %>%
  mutate(n_Genomes=n_distinct(Genome), Genomes=toString(unique(Genome))) %>%
  mutate(Method="Accessions")

pfam.counts.dom=pfam %>%
  group_by(Genome, Pfam) %>%
  summarize(Count=length(Accession)) %>%
  group_by(Genome) %>%
  mutate(n=sum(Count)) %>%
  group_by(Pfam) %>%
  mutate(n_Genomes=n_distinct(Genome), Genomes=toString(unique(Genome))) %>%
  mutate(Method="Domains")

p.fams.combined=bind_rows(pfam.counts.acc, pfam.counts.dom) %>%
  select(Genome, Pfam, Count, Method) %>%
  pivot_wider(names_from="Method", values_from="Count")

pfam.genome=pfam %>%
  group_by(Accession, Genome) %>%
  mutate(n=n_distinct(Pfam))

pfam.composition=pfam.counts.acc %>%
  group_by(Genome) %>%
  summarize(n=sum(n_Genomes), pfams=list(unique(Pfam)))

pfam.unique=pfam.counts.acc %>%
  filter(n_Genomes==1) %>%
  group_by(Genome) %>%
  summarize(n=n_distinct(Pfam), pfam=list(unique(Pfam)))

pfam.dat=pfam.counts.acc %>%
  filter(n_Genomes>1) %>%
  select(-n_Genomes) 

pfam.phylo=levels(as.factor(pfam.dat$Genome))[c(1, 2, 5, 6, 7, 3, 4)]

upset.dat=pfam.composition$pfams
names(upset.dat)=pfam.composition$Genome

pfam.uplt=upset(fromList(upset.dat), sets=pfam.phylo, mb.ratio = c(0.55, 0.45), order.by = "freq", keep.order = TRUE)

pfam.uplt
grid.text("Pfam UpSet Plot",x = 0.65, y=0.95, gp=gpar(fontsize=20))

#Pfam enrichment
pfam.res=enrich(pfam.dat, var="Pfam", group="Genome")

pfam.res2 = pfam.res %>%
  group_by(group1, group2) %>%
  summarize(Domains=toString(Domain)) %>%
  pivot_wider(id_cols=group1, names_from=group2, values_from=Domains)

pfam.res3=pfam.res2[match(rev(pfam.phylo), pfam.res2$group1),] %>%
  drop_na(group1)

pfam.res4=pfam.res3[,match(c("group1", intersect(rev(pfam.phylo), colnames(pfam.res3))), colnames(pfam.res3)),]

pfam.res.stat=data.frame(Genome=c(pfam.res$group1, pfam.res$group2), Domain=c(pfam.res$Domain, pfam.res$Domain)) %>%
  group_by(Genome, Domain) %>%
  summarize(sig_n=length(Domain))

pfam.meta=pfam.dat %>%
  select(Pfam, Genome, Count) %>%
  rename(Domain="Pfam") %>%
  group_by(Domain) %>%
  mutate(avg=median(Count), Fold=ifelse((Count+avg)>0, Count/avg, 0),
         Direction=ifelse(Fold>1 | Fold==Inf, "Up", ifelse(Fold==1 | Fold==0, ifelse(Count<avg, "Down", "Equal"), "Down"))) %>%
  mutate(Fold=ifelse(Direction=="Down", 1/Fold, Fold))

pfam.res.stat.meta=inner_join(pfam.res.stat, pfam.meta)

pfam.res.stat.meta$Genome=as.factor(pfam.res.stat.meta$Genome)

#Pfam clustering
pfam.dist=pfam.res.stat.meta %>%
  select(Domain, Genome, sig_n) %>%
  pivot_wider(id_cols="Domain", values_from=sig_n, names_from=Genome) %>% mutate(
    across(where(is.numeric), ~replace_na(.x, 0))
  )

pfam.dist2=as.matrix(pfam.dist[2:length(pfam.dist)])
names(pfam.dist2)=pfam.dist$Domain

pfam.hclust=hclust(d = dist(pfam.dist2))

pfam.hclust2 <- data.frame(Domain=names(pfam.dist2),
                           cluster=cutree(pfam.hclust, length(pfam.dist$Domain)))

pfam.clust=left_join(pfam.res.stat.meta, pfam.hclust2)

pfam.clust$Genome=as.factor(pfam.clust$Genome)
pfam.clust$Genome=factor(pfam.clust$Genome, levels=levels(pfam.clust$Genome)[match(pfam.phylo, levels(pfam.clust$Genome))])

pfam.list = unique(subset(pfam.res.stat.meta, pfam.res.stat.meta$Genome=="EMU")$Domain)

Emus.pfam.stat=pfam.clust %>%
  filter(Domain %in% pfam.list) %>%
  filter(Genome=="EMU") %>%
  ungroup() %>%
  select(Domain, Direction) %>%
  distinct() %>%
  rename("Emus.status"="Direction")

Emus.pfam.clust=pfam.clust %>%
  left_join(Emus.pfam.stat) %>%
  mutate(Emus.status=replace_na(Emus.status, "Missing")) %>%
  filter(Emus.status!="Equal" & Emus.status!="Missing") %>%
  mutate(Direction=as.factor(Direction)) %>%
  mutate(Direction=factor(Direction, levels=levels(Direction)[c(3, 2, 1)])) %>%
  left_join(sigp.pfam2)

#Pfam visualization
pfam.dwn=ggplot(subset(Emus.pfam.clust, Emus.status=="Down"), aes(y=reorder(Domain, cluster), x=Genome, color=Direction, size=sig_n))+geom_point(alpha=0.8)+
  theme+geom_stripes(inherit.aes=F, aes(y=reorder(Domain, cluster)), odd = "#33333333", even = "#00000000")+theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0.95))+
  theme(axis.text.y=element_text(size=15), axis.text.x=element_text(size=15))+theme(legend.position="right")+ylab("Pfam")+facet_grid(~Secretion, scales="free_x", space="free_x")+scale_color_manual(values=c("#127852", "#8B8588", "#FF0000"))+labs(color="Direction\nvs. median", size="Significant\nDifferences")+coord_flip()+guides(color = guide_legend(override.aes = list(size=3), order=1))

pfam.dwn

pfam.up=ggplot(subset(Emus.pfam.clust, Emus.status=="Up"), aes(y=reorder(Domain, cluster), x=Genome, color=Direction, size=sig_n))+geom_point(alpha=0.8)+
  theme+geom_stripes(inherit.aes=F, aes(y=reorder(Domain, cluster)), odd = "#33333333", even = "#00000000")+theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0.95))+
  theme(axis.text.y=element_text(size=15), axis.text.x=element_text(size=8))+theme(legend.position="right")+ylab("Pfam")+facet_grid(~Secretion, scales="free_x", space="free_x")+scale_color_manual(values=c("#127852", "#8B8588", "#FF0000"))+labs(color="Direction\nvs. median", size="Significant\nDifferences")+coord_flip()+guides(color = guide_legend(override.aes = list(size=3), order=1))

pfam.up

pfam.dwn2=ggplot(subset(Emus.pfam.clust, Emus.status=="Down" & Genome=="EMU"), aes(x=Count, y=reorder(Domain, cluster), fill=Fold))+geom_col()+
  theme+geom_stripes(odd = "#33333333", even = "#00000000")+
  theme(axis.text.y=element_text(size=15), axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x = element_blank())+theme(legend.position="right")+xlab("Count for EMU only")+facet_grid(~Secretion, scales="free_x", space="free_x")+coord_flip()+labs(fill="Fold vs.\nmedian")+scale_fill_viridis_c()+ggtitle("Pfam Down in EMU compared to median")+theme(plot.title = element_text(hjust = 0.5))

pfam.dwn2

pfam.up2=ggplot(subset(Emus.pfam.clust, Emus.status=="Up" & Genome=="EMU"), aes(x=Count, y=reorder(Domain, cluster), fill=Fold))+geom_col()+
  theme+geom_stripes(odd = "#33333333", even = "#00000000")+
  theme(axis.text.y=element_text(size=10), axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x = element_blank())+theme(legend.position="right")+xlab("Count for EMU only")+facet_grid(~Secretion, scales="free_x", space="free_x")+coord_flip()+labs(fill="Fold vs.\nmedian")+scale_fill_viridis_c()+ggtitle("Pfam Up in EMU compared to median")+theme(plot.title = element_text(hjust = 0.5))

pfam.up2

Emus.pfam.missing=pfam.clust %>%
  left_join(Emus.pfam.stat) %>%
  mutate(Emus.status=as.character(Emus.status)) %>%
  mutate(Emus.status=replace_na(Emus.status, "Missing from EMU")) %>%
  filter(Emus.status=="Missing from EMU") %>%
  mutate(Direction=as.factor(Direction)) %>%
  mutate(Direction=factor(Direction, levels=levels(Direction)[c(3, 2, 1)])) %>%
  left_join(sigp.pfam2)

pfam.plt3=ggplot(Emus.pfam.missing, aes(y=reorder(Domain, rev(cluster)), x=Genome, color=Direction, size=as.factor(sig_n)))+geom_point(alpha=0.8)+
  theme+theme(axis.text.x = element_text(size=15, angle = 90, vjust=0.5, hjust=0.95))+theme(axis.text.y=element_text(size=9))+theme(legend.position="right")+xlab("Genome")+scale_color_manual(values=c("#127852", "#8B8588", "#FF0000"))+labs(size="Significant\nDifferences")+geom_stripes(inherit.aes=F, aes(y=reorder(Domain, rev(cluster))), odd ="#33333333", even = "#00000000")+scale_x_discrete(limits = rev(pfam.phylo))+ylab("Pfam")+ggtitle("Missing from EMU")+theme(plot.title = element_text(hjust = 0.5))+guides(color = guide_legend(override.aes = list(size=3), order=1))+labs(color="Direction vs.\nmedian")

pfam.plt3

Emus.pfam.unique=pfam.counts %>%
  filter(Pfam %in% unlist(pfam.unique[pfam.unique$Genome=="EMU",]$pfam))  %>%
  rename("Domain"="Pfam") %>%
  left_join(sigp.pfam2)

pfam.plt4=ggplot(Emus.pfam.unique, aes(x=Domain, y=Count, fill=Secretion))+geom_col()+scale_fill_viridis_d()+
  theme+theme(axis.text.x = element_text(angle = 90, vjust=0.5))+xlab("Pfam")+scale_x_discrete(limits = rev(Emus.pfam.unique$Domain))+coord_flip()+ggtitle("Unique to EMU")+theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")+labs(fill="Prediction")

pfam.plt4

pfam.dwn.1=ggplotGrob(pfam.dwn)
pfam.dwn.2=ggplotGrob(pfam.dwn2)

maxWidth <- unit.pmax(pfam.dwn.1$widths, pfam.dwn.2$widths)

pfam.dwn.1$widths= maxWidth
pfam.dwn.2$widths= maxWidth

layout1 <- rbind(c(2),
                 c(1))

grid.arrange(pfam.dwn.1, pfam.dwn.2, layout_matrix=layout1, heights=c(0.4, 0.6))

pfam.up.1=ggplotGrob(pfam.up)
pfam.up.2=ggplotGrob(pfam.up2)

maxWidth <- unit.pmax(pfam.up.1$widths, pfam.up.2$widths)

pfam.up.1$widths= maxWidth
pfam.up.2$widths= maxWidth

grid.arrange(pfam.up.1, pfam.up.2, layout_matrix=layout1, heights=c(0.4, 0.6))

pfam3=ggplotGrob(pfam.plt3)
pfam4=ggplotGrob(pfam.plt4)

layout2 <- rbind(c(1,2))

grid.arrange(pfam4, pfam3, layout_matrix=layout2, widths=c(0.5, 0.5))


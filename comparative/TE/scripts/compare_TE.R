library(tidyverse)
library(purrr)
library(readr)
library(fs)
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(viridis)

data_dir <- "data"
index_dir = "indexes"
infiles <- list.files(path=data_dir,pattern = "fasta.out.bz2$", recursive = TRUE,full.names=TRUE)
tefiles = gsub("data/","",infiles)
tefiles = gsub("(.Nanopore10X|.v\\d+)\\S+$","",tefiles,perl=TRUE)
tefiles <- tibble::rowid_to_column(data.frame(tefiles),"source")
indexes = list.files(path=index_dir,pattern = "fai$", recursive = TRUE,full.names=TRUE)
repland = list.files(path=data_dir,pattern = "divsum.tbl$", recursive = TRUE,full.names=TRUE)
genomenames = gsub("indexes/","",indexes)
genomenames = gsub("(.Nanopore10X|.v\\d+|.scaffolds)\\S+$","",genomenames,perl=TRUE) 

genomenames <- tibble::rowid_to_column(data.frame(genomenames),"source")

genome_sizes <- indexes %>% map_dfr(read_tsv, .id="source",col_names = c("Scaffold","length","offset","linebases","linewidth"), col_types=c("ciii")) %>%
    group_by(source) %>% summarize(totallen = sum(length)) %>% type_convert(col_types="ii")

genome_sizes <- genome_sizes %>% left_join(genomenames,by="source") %>% select(genomenames,totallen) %>% rename(genomename=genomenames,length=totallen) 

tefiles2=tefiles %>%
  mutate(species=tefiles,
         species=str_replace(species, "Entomophaga_maimaiga_var_ARSEF_7190", "EMA"),
         species=str_replace(species, "Entomophthora_muscae_UCB", "EMU"),
         species=str_replace(species, "Massospora_cicadina_MCPNR19", "MCI"),
         species=str_replace(species, "Zoophthora_radicans_ATCC_208865", "ZRA"),
         species=str_replace(species, "Conidiobolus_thromboides_FSU_785.Conth1", "CTH"),
         species=factor(species, levels=c("EMU", "EMA", "MCI", "ZRA", "CTH")))

tbl <- infiles %>% map_dfr(read_table,.id = "source",skip=3,comment="#",col_names=c("Score","div","del","ins","query",
                                                                                  "qstart","qend","qleft","strand","match_name",
                                                                                  "family", "hstart", "hend", "hleft", "ID", 'Overlaps_another'
                                                                                  ),na = ".") %>% separate(family, c("superfamily","subfamily"),sep="/",
                                                                                                           extra = "drop", fill = "right") %>% 
    type_convert(col_types="iidddciiccccccicic") %>% left_join(tefiles2,by="source") %>% select(-source, ) %>% rename(genomename=tefiles)

Unknown <- tbl %>% filter(superfamily == "Unknown")
SimpleRepeat <- tbl %>% filter(superfamily =="Simple_repeat")
rRNA <- tbl %>% filter(superfamily %in% c("rRNA"))
helitron <- tbl %>% filter(superfamily %in% c("RC"))

KnownTbl <- tbl %>% filter( !superfamily  %in% c("Unknown","Simple_repeat","Satellite","rRNA","Satellite","Low_complexity") ) %>% mutate(flen=abs(qend-qstart)+1)
DNAElem <- KnownTbl %>% filter(superfamily %in% c("DNA"))

DNAElem %>% mutate(subfamily=gsub("hAT-\\S+","hAT",subfamily,perl=TRUE)) %>% group_by(genomename,subfamily) %>% summarise( total_TE_len = sum(flen) ) %>% left_join(genome_sizes,by="genomename") %>% mutate(percent = 100 * (total_TE_len / length)) 

unique(KnownTbl$superfamily)
TEcounts <- KnownTbl %>% group_by(genomename,superfamily) %>% summarise( total_TE_len = sum(flen) ) %>% left_join(genome_sizes,by="genomename") %>% 
  mutate(percent = 100 * (total_TE_len / length)) %>% select(genomename,superfamily,total_TE_len,length,percent) %>% rename(species=genomename,genomelength=length) %>% 
  mutate_at("superfamily", str_replace, "RC", "Helitron")

TEcounts

TEcounts2 <- KnownTbl %>% mutate(subfamily=replace_na(subfamily,"")) %>%
  separate(subfamily, into=c("subfamily", "type"), sep="-") %>%
  mutate(subfamily=ifelse(subfamily=="", superfamily, subfamily)) %>%
  group_by(genomename,subfamily, superfamily) %>% 
  summarise(total_TE_len = sum(flen) ) %>% left_join(genome_sizes,by="genomename") %>% 
  mutate(percent = total_TE_len / length) %>% select(genomename,superfamily,subfamily,total_TE_len,length,percent) %>% rename(species=genomename,genomelen=length)
TEcounts2
#write_csv(TEcounts2,"TE_subfamily_counts.csv")

TEcounts3 <- TEcounts2 %>%
  mutate(species=str_replace(species, "Entomophaga_maimaiga_var_ARSEF_7190", "EMA"),
         species=str_replace(species, "Entomophthora_muscae_UCB", "EMU"),
         species=str_replace(species, "Massospora_cicadina_MCPNR19", "MCI"),
         species=str_replace(species, "Zoophthora_radicans_ATCC_208865", "ZRA"),
         species=str_replace(species, "Conidiobolus_thromboides_FSU_785.Conth1", "CTH"),
         species=factor(species, levels=c("EMU", "EMA", "MCI", "ZRA", "CTH"))) 

theme = theme_bw()+theme(text = element_text(size=15), axis.title.x = element_text(size=20), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15), title = element_text(size=20), legend.title = element_text(size=20), legend.text = element_text(size=15), strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid"))

p<-ggplot(TEcounts3,aes(x=species,y=percent,fill=superfamily)) + theme + geom_bar(position="dodge",stat="identity") + scale_fill_viridis_d() + scale_y_continuous(labels=percent)+ ylab("Percent")+ xlab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
#ggsave("TE_species.pdf",p,width=10,height=12)

p<-ggplot(TEcounts3 %>% filter(superfamily=="DNA",percent>0.001),aes(x=species,y=percent,fill=subfamily))+ theme+ geom_bar(position="dodge",stat="identity") + scale_fill_viridis_d() + scale_y_continuous(labels=percent)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Percent")+ xlab("")+facet_wrap(~superfamily)
p
#ggsave("TE_DNA_species.pdf",p,width=10,height=12)

p<-ggplot(TEcounts3 %>% filter(superfamily=="LTR"),aes(x=species,y=percent,fill=subfamily)) + theme + geom_bar(position="dodge",stat="identity") + scale_fill_viridis_d() + scale_y_continuous(labels=percent)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Percent")+ xlab("")+facet_wrap(~superfamily)
p
#ggsave("TE_LTR_species.pdf",p,width=10,height=12)

p<-ggplot(TEcounts3 %>% filter(superfamily=="LINE"),aes(x=species,y=percent,fill=subfamily)) + theme+ geom_bar(position="dodge",stat="identity") + 
  scale_fill_viridis_d() + scale_y_continuous(labels=percent)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Percent")+ xlab("")+facet_wrap(~superfamily)
p
#ggsave("TE_LINE_species.pdf",p,width=10,height=12)

#write_csv(tbl,"RM_scores.csv.bz2")

####Landcape analyses####
class_tbl <- repland %>% map_dfr(read.csv, .id = "source",sep="") %>% mutate(source=as.integer(source)) %>% left_join(tefiles2,by="source") %>% select(-source, -tefiles) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  pivot_longer(cols=-c(Div, species), names_to="Class", values_to="Count") %>%
  group_by(Class, species) %>%
  dplyr::mutate(Diff = Count-lag(Count, default = 0, order_by = -Div)) %>%
  separate(Class, into=c("Type"), extra="drop", remove=F) %>%
  rename("Genome"=species)

thresh=0.01

DNA_comb=class_tbl %>%
  filter(Type=="DNA" & Count>0) %>%
  arrange(Type, Class) %>%
  group_by(Class, Genome) %>%
  mutate(sm=sum(Count)) %>%
  ungroup() %>%
  mutate(tot=sum(Count), per=sm/tot) %>%
  mutate(Class=ifelse(per<thresh, "Other", Class))

plt1=ggplot(DNA_comb, aes(x=Div, y=Count, color=Class))+geom_point(size=1)+geom_line(size=1)+theme+scale_color_viridis_d()+facet_wrap(~Genome)+theme(axis.text.x=element_text(angle=90, size=15))+xlab("Div")

plt1

plt2=ggplot(DNA_comb, aes(x=Div, y=Diff, color=Class))+geom_line(size=1)+theme+scale_color_viridis_d()+facet_wrap(~Genome)+theme(axis.text.x=element_text(angle=90, size=15))+xlab("Div")

plt2

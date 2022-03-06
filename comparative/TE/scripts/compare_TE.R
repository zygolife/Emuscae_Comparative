library(tidyverse)
library(purrr)
library(readr)
library(fs)
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
data_dir <- "data"
index_dir = "indexes"
infiles <- list.files(path=data_dir,pattern = "fasta.out.bz2$", recursive = TRUE,full.names=TRUE)
tefiles = gsub("data/","",infiles)
tefiles = gsub("(.Nanopore10X|.v\\d+)\\S+$","",tefiles,perl=TRUE)
tefiles <- tibble::rowid_to_column(data.frame(tefiles),"source")
indexes = list.files(path=index_dir,pattern = "fai$", recursive = TRUE,full.names=TRUE)
genomenames = gsub("indexes/","",indexes)
genomenames = gsub("(.Nanopore10X|.v\\d+|.scaffolds)\\S+$","",genomenames,perl=TRUE) 

genomenames <- tibble::rowid_to_column(data.frame(genomenames),"source")
genome_sizes <- indexes %>% map_dfr(read_tsv, .id="source",col_names = c("Scaffold","length","offset","linebases","linewidth"), col_types=c("ciii")) %>%
    group_by(source) %>% summarize(totallen = sum(length)) %>% type_convert(col_types="ii")


genome_sizes <- genome_sizes %>% left_join(genomenames,by="source") %>% select(genomenames,totallen) %>% rename(genomename=genomenames,length=totallen) 

tbl <- infiles %>% map_dfr(read_table,.id = "source",skip=3,comment="#",col_names=c("Score","div","del","ins","query",
                                                                                  "qstart","qend","qleft","strand","match_name",
                                                                                  "family", "hstart", "hend", "hleft", "ID", 'Overlaps_another'
                                                                                  ),na = ".") %>% separate(family, c("superfamily","subfamily"),sep="/",
                                                                                                           extra = "drop", fill = "right") %>% 
    type_convert(col_types="iidddciiccccccicic") %>% left_join(tefiles,by="source") %>% select(-source) %>% rename(genomename=tefiles)

Unknown <- tbl %>% filter(superfamily == "Unknown")
SimpleRepeat <- tbl %>% filter(superfamily =="Simple_repeat")
rRNA <- tbl %>% filter(superfamily %in% c("rRNA"))
helitron <- tbl %>% filter(superfamily %in% c("RC"))
KnownTbl <- tbl %>% filter( !superfamily  %in% c("Unknown","Simple_repeat","Satellite","rRNA","Satellite","Low_complexity") ) %>% mutate(flen=abs(qend-qstart)+1)
unique(KnownTbl$superfamily)
TEcounts <- KnownTbl %>% group_by(genomename,superfamily) %>% summarise( total_TE_len = sum(flen) ) %>% left_join(genome_sizes,by="genomename") %>% 
  mutate(percent = 100 * (total_TE_len / length)) %>% select(genomename,superfamily,percent) %>% rename(species=genomename) %>% mutate_at("superfamily", str_replace, "RC", "Helitron")
TEcounts
p<-ggplot(TEcounts,aes(x=species,y=percent,fill=superfamily)) + geom_bar(position="dodge",stat="identity") + scale_fill_brewer(palette = "Set1") + scale_colour_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("TE_species.pdf",p,width=10,height=12)
TEcounts2 <- KnownTbl %>% group_by(genomename,superfamily,subfamily) %>% summarise( total_TE_len = sum(flen) ) %>% left_join(genome_sizes,by="genomename") %>% 
  mutate(percent = 100 * (total_TE_len / length)) %>% select(genomename,superfamily,subfamily,percent) %>% rename(species=genomename) %>% mutate_at("superfamily", str_replace, "RC", "Helitron")
TEcounts2

p<-ggplot(TEcounts2 %>% filter(superfamily=="DNA"),aes(x=species,y=percent,fill=subfamily)) + geom_bar(position="dodge",stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("TE_DNA_species.pdf",p,width=10,height=12)

p<-ggplot(TEcounts2 %>% filter(superfamily=="LTR"),aes(x=species,y=percent,fill=subfamily)) + geom_bar(position="dodge",stat="identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("TE_LTR_species.pdf",p,width=10,height=12)

#write_csv(tbl,"RM_scores.csv.bz2")

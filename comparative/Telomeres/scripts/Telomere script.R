library(tidyverse)

data_dir <- "reports"
infiles <- list.files(path=data_dir, pattern = ".txt$", recursive = TRUE,full.names=TRUE)
telofiles = gsub("reports/", "", infiles)
telofiles_key <- tibble::rowid_to_column(data.frame(telofiles),"source") %>%
  mutate(Genome=as.factor(c("APS", "EMA", "EMU", "MCI", "MPL", "NCR", "ZRA")))

teloresults <- infiles %>% map_dfr(read_tsv, .id="source",col_names = c("Scaffold","Direction","Sequence"), col_types=c("ccc")) %>%
  group_by(source) %>% filter(!grepl("Telomere", Scaffold)) %>%
  mutate(source=as.integer(source), Scaffold=gsub(" .*", "", Scaffold)) %>%
  left_join(telofiles_key %>% select(source, Genome)) %>% ungroup() %>%
  filter(Genome!="NCR" & Genome!="MPL")

telo_motif1="TAA[C]+"
telo_motif2="[G]+TTA"

GTGTTTTGTTCCCCCTTGATCCACTCTTTCTCTATGTTGCAAGCTTTGGCGCCATCAACCTTGGAGGGGTTAGTTACAAGAGGAGGTTTTGTGGGGACAT

teloresults2 <- teloresults %>%
  select(-source) %>%
  mutate(Sequence=toupper(Sequence)) %>%
  mutate(len=str_length(Sequence), count1=str_count(Sequence, telo_motif1), count2=str_count(Sequence, telo_motif2), sum_motifs=count1+count2, both_motifs=(count1*count2)>0) %>%
  pivot_wider(names_from="Direction", values_from="Direction", values_fn=length) %>%
  mutate(across(c("forward", "reverse"), ~replace_na(.x, 0)),
         both=ifelse(forward==1 & reverse==1, 1, 0))

teloresult3 <- teloresults2 %>%
  group_by(Genome) %>%
  summarize(FWD=sum(forward), REV=sum(reverse), Both=sum(both), Total_w_telo=FWD+REV-Both,
            min.motifs=min(sum_motifs), med.motifs=median(sum_motifs), max.motifs=max(sum_motifs), sd.motifs=sd(sum_motifs))




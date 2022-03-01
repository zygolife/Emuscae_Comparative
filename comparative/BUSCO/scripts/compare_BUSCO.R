library(tidyverse)
library(purrr)
library(readr)
library(fs)
library(UpSetR)
library(rJava)
library(venneuler)
library(grid)
data_dir <- "BUSCO_genome"
infiles <- list.files(path=data_dir,pattern = "full_table.tsv$", recursive = TRUE,full.names=TRUE)
names <- tibble(Taxon=gsub("\\.(masked|scaffolds\\.fa)","",basename(dirname(dirname(infiles))))) %>% rowid_to_column("source")

tbl <- infiles %>% map_dfr(read_tsv,.id = "source",skip=3,comment="#",col_names=c("BUSCO_ID","Status","Sequence","Gene Start","Gene End",
                                                "Strand","Score","Length","OrthoDB url","Description"),na = ".") %>% select(source,BUSCO_ID,Status) %>% type_convert(col_types="icc")


busco_status <- tbl %>% left_join(names) %>% select(-source) %>% mutate(Found = factor(case_when(Status == "Complete" ~ TRUE, 
                                                                                       Status == "Fragmented" ~ TRUE,
                                                                                       Status == "Duplicated" ~ TRUE,
                                                                                        Status == "Missing" ~ FALSE))) %>% select(-Status)
#busco_status %>% pivot_wider(names_from="BUSCO_ID",values_from=Found)
#busco_status <- tbl %>% left_join(names) %>% select(-source) %>% mutate(Found = case_when(Status == "Complete" ~ 1, Status == "Missing" ~ 0)) %>% select(-Status)

busco_fungi <- busco_status %>%
  mutate(Taxon="Fungi Buscos") %>%
  distinct() %>%
  group_by(Taxon) %>%
  summarize(Buscos=list(unique(BUSCO_ID)))

busco_genome <- busco_status %>%
  distinct() %>%
  filter(Found==T) %>%
  mutate(Taxon=str_sub(Taxon, end=20)) %>%
  mutate(Taxon=paste(Taxon, "Genome", sep="")) %>%
  distinct() %>%
  group_by(Taxon) %>%
  summarize(Buscos=list(unique(BUSCO_ID)))

data_dir <- "BUSCO_pep"
infiles <- list.files(path=data_dir,pattern = "full_table.tsv$", recursive = TRUE,full.names=TRUE)
names <- tibble(Taxon=gsub("\\.(masked|scaffolds\\.fa)","",basename(dirname(dirname(infiles))))) %>% rowid_to_column("source")

tbl <- infiles %>% map_dfr(read_tsv,.id = "source",skip=3,comment="#",col_names=c("BUSCO_ID","Status","Sequence","Gene Start","Gene End",
                                                                                  "Strand","Score","Length","OrthoDB url","Description"),na = ".") %>% select(source,BUSCO_ID,Status) %>% type_convert(col_types="icc")


busco_status <- tbl %>% left_join(names) %>% select(-source) %>% mutate(Found = factor(case_when(Status == "Complete" ~ TRUE, 
                                                                                                 Status == "Fragmented" ~ TRUE,
                                                                                                 Status == "Duplicated" ~ TRUE,
                                                                                                 Status == "Missing" ~ FALSE))) %>% select(-Status)

busco_pep <- busco_status %>%
  distinct() %>%
  filter(Found==T) %>%
  mutate(Taxon=str_sub(Taxon, end=20)) %>%
  mutate(Taxon=paste(Taxon, Taxon="Pep", sep="")) %>%
  distinct() %>%
  group_by(Taxon) %>%
  summarize(Buscos=list(unique(BUSCO_ID)))

data_dir <- "BUSCO_genome_euk"
infiles <- list.files(path=data_dir,pattern = "full_table.tsv$", recursive = TRUE,full.names=TRUE)
names <- tibble(Taxon=gsub("\\.(masked|scaffolds\\.fa)","",basename(dirname(dirname(infiles))))) %>% rowid_to_column("source")

tbl <- infiles %>% map_dfr(read_tsv,.id = "source",skip=3,comment="#",col_names=c("BUSCO_ID","Status","Sequence","Gene Start","Gene End",
                                                                                  "Strand","Score","Length","OrthoDB url","Description"),na = ".") %>% select(source,BUSCO_ID,Status) %>% type_convert(col_types="icc")


busco_status <- tbl %>% left_join(names) %>% select(-source) %>% mutate(Found = factor(case_when(Status == "Complete" ~ TRUE, 
                                                                                                 Status == "Fragmented" ~ TRUE,
                                                                                                 Status == "Duplicated" ~ TRUE,
                                                                                                 Status == "Missing" ~ FALSE))) %>% select(-Status)

busco_euk <- busco_status %>%
  mutate(Taxon="Euk Buscos") %>%
  distinct() %>%
  group_by(Taxon) %>%
  summarize(Buscos=list(unique(BUSCO_ID)))

busco_genome_euk <- busco_status %>%
  distinct() %>%
  filter(Found==T) %>%
  mutate(Taxon=str_sub(Taxon, end=20)) %>%
  mutate(Taxon=paste(Taxon, "Genome_euk", sep="")) %>%
  distinct() %>%
  group_by(Taxon) %>%
  summarize(Buscos=list(unique(BUSCO_ID)))

data_dir <- "BUSCO_pep_euk"
infiles <- list.files(path=data_dir,pattern = "full_table.tsv$", recursive = TRUE,full.names=TRUE)
names <- tibble(Taxon=gsub("\\.(masked|scaffolds\\.fa)","",basename(dirname(dirname(infiles))))) %>% rowid_to_column("source")

tbl <- infiles %>% map_dfr(read_tsv,.id = "source",skip=3,comment="#",col_names=c("BUSCO_ID","Status","Sequence","Gene Start","Gene End",
                                                                                  "Strand","Score","Length","OrthoDB url","Description"),na = ".") %>% select(source,BUSCO_ID,Status) %>% type_convert(col_types="icc")

busco_status <- tbl %>% left_join(names) %>% select(-source) %>% mutate(Found = factor(case_when(Status == "Complete" ~ TRUE, 
                                                                                                 Status == "Fragmented" ~ TRUE,
                                                                                                 Status == "Duplicated" ~ TRUE,
                                                                                                 Status == "Missing" ~ FALSE))) %>% select(-Status)

busco_pep_euk <- busco_status %>%
  distinct() %>%
  filter(Found==T) %>%
  mutate(Taxon=str_sub(Taxon, end=20)) %>%
  mutate(Taxon=paste(Taxon, "Pep_euk", sep="")) %>%
  distinct() %>%
  group_by(Taxon) %>%
  summarize(Buscos=list(unique(BUSCO_ID)))

busco_status2 <- bind_rows(busco_fungi,
                           busco_pep,
                           busco_genome)

list <- busco_status2$Buscos
names(list) <- busco_status2$Taxon

upset(fromList(list), order.by = "freq", keep.order=T, sets=names(list),
      nsets = 11,  text.scale=1.5, point.size = 3.5, line.size = 2,
      mainbar.y.label = "BUSCO marker recovery in Entos", sets.x.label = "Species Found"
)

grid.text("Fungi UpSet Plot", x = 0.65, y=0.95, gp=gpar(fontsize=15))

busco_found <- unique(unlist(busco_status2[busco_status2$Taxon!="Fungi Buscos",]$Buscos))

busco_fungi <- unique(unlist(busco_status2[busco_status2$Taxon=="Fungi Buscos",]$Buscos))

missing_fungi=setdiff(busco_fungi, busco_found)

missing_fungi

busco_status3 <- bind_rows(busco_euk,
                           busco_pep_euk,
                           busco_genome_euk)

list2=busco_status3$Buscos
names(list2)=busco_status3$Taxon

upset(fromList(list2), order.by = "freq", keep.order = TRUE,
      nsets = 11,  text.scale=1.5, point.size = 3.5, line.size = 2,
      mainbar.y.label = "BUSCO marker recovery in Entos", sets.x.label = "Species Found"
)

grid.text("Euk UpSet Plot", x = 0.65, y=0.95, gp=gpar(fontsize=15))

busco_found_euk <- unique(unlist(busco_status3[busco_status3$Taxon!="Euk Buscos",]$Buscos))

busco_euk <- unique(unlist(busco_status3[busco_status3$Taxon=="Euk Buscos",]$Buscos))

missing_euk=setdiff(busco_euk, busco_found_euk)

missing_euk

grid.text(
  "@hyphaltip",
  x = 0.90,
  y = 0.05,
  gp = gpar(
    fontsize = 10,
    fontface = 3
  )
)

# Venn drawing
foundset <- busco_statusVenn %>% select(-Status,-Found)
v <- venneuler(data.frame(foundset))

par(cex = 0.7) 
plot(v, main = "BUSCO set", cex.main = 1.5)
grid.text(
  "@hyphaltip",
  x = 0.52,
  y = 0.15,
  gp = gpar(
    fontsize = 10,
    fontface = 3
  )
)

# upset(foundsets,
#       query.legend = "bottom", nsets = 10, number.angles = 30, point.size = 3.5, line.size = 2,
#       mainbar.y.label = "BUSCO marker", sets.x.label = "Species Foun", 
#       queries = list(
#         list(
#           query = elements,
#           params = list("physicalActivityPerMonth", 0,4),
#           color = "#Df5286", 
#           active = T,
#           query.name = "Physically Active 1x/Week or Less"
#         )
#       ), 
#       attribute.plots = list(gridrows = 50, 
#                              plots = list(list(plot = histogram, x = "volunteerPerMonth", queries = T), 
#                                           list(plot = histogram, x = "minAgeRange", queries = T), 
#                                           list(plot = scatter_plot, x = "minAgeRange", y="volunteerPerMonth", queries = F)
#                              ), 
#                              ncols = 3
#       ) 
# )
# grid.text(
#   "@littlemissdata",
#   x = 0.9,
#   y = 0.02,
#   gp = gpar(
#     fontsize = 10,
#     fontface = 3
#   )
# )
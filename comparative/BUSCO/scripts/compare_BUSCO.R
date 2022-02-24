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

movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
                   header = T, sep = ";")

upset(busco_status,
      nsets = 5, number.angles = 30, point.size = 3.5, line.size = 2,
      mainbar.y.label = "BUSCO marker recovery in Entos", sets.x.label = "Species Found"
)
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
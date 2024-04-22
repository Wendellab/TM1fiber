library(RColorBrewer)
library(data.table)
library(tidyverse)
library(ggrepel)

options(scipen = 999)
set.seed(8675309)
setDTthreads(18) # setting the data.table threads
options(datatable.fread.datatable=FALSE)


load(file = "04.R-01-dataInput.filtered.RData")

PStable <- fread('PCWSCW.tsv')


exprNorm <- t(datExpr) %>% 
  as.data.frame() %>% 
  rownames_to_column("homoeolog") %>% 
  right_join(.,PStable)

CesAexpr <- exprNorm %>%
  pivot_longer(!(c(homoeolog,CesA,stage)),names_to="Sample", values_to="expression") %>% 
  left_join(.,coldata[,c("id","DPA")],by=c("Sample"="id")) %>%
  mutate(DPA = as.numeric(sub("DPA","",DPA))) %>% 
  mutate(label=ifelse(Sample=="AD1_TM1_plant10_24DPA",CesA,NA)) %>% 
  filter(!(CesA=="GhCESA6-E-At" | CesA=="GhCESA6-E-Dt")) %>% 
  add_row(homoeolog="Gorai", CesA="GhCESA6-Fa-Dt", stage="PCW", Sample=NULL,expression=-999,DPA=11,label=NULL) %>%
  mutate(CWgroup = stage) %>% 
  mutate(CWgroup = ifelse(CesA == "GhCESA8-A-Dt", "PCW", CWgroup)) %>%
  mutate(CWgroup = ifelse(CesA == "GhCESA8-A-At", "PCW", CWgroup))


CesAcolorstmp <- c("GhCESA7-A-At"="#FAEBD7",
                   "GhCESA7-A-Dt"= "#7FFFD4",
                   "GhCESA4-A-At"= "#FFE4C4",
                   "GhCESA4-A-Dt"= "#F0FFF0",
                   "GhCESA6-D-At"= "#F0FBF2",
                   "GhCESA6-D-Dt"= "#FF1493",
                   "GhCESA6-A-At"= "#2F4F4F",
                   "GhCESA6-A-Dt"= "#8FBC8F",
                   "GhCESA3-C-At"= "#9932CC",
                   "GhCESA3-C-Dt"= "#F0FFFF",
                   "GhCESA4-B-At"= "#FF8C00",
                   "GhCESA4-B-Dt"= "#556B2F",
                   "GhCESA3-A-At"= "#B8860B",
                   "GhCESA3-A-Dt"= "#00FFFF",
                   "GhCESA3-B-At"= "#FFF8DC",
                   "GhCESA3-B-Dt"= "#FF7F50",
                   "GhCESA6-Fb-At"= "#00FF00",
                   "GhCESA6-Fb-Dt"="#DAA520",
                   "GhCESA6-Fa-At"= "#FFD700",
                   "GhCESA1-B-At"= "#B22222",
                   "GhCESA1-B-Dt"= "#1E90FF",
                   "GhCESA7-B-At"= "#00BFFF",
                   "GhCESA7-B-Dt"= "#7FFF00",
                   "GhCESA1-A-At"= "#D2691E",
                   "GhCESA1-A-Dt"= "#5F9EA0",
                   "GhCESA8-A-At"= "#DEB887",
                   "GhCESA8-A-Dt"="#A52A2A",
                   "GhCESA6-B-At"= "#0000FF",
                   "GhCESA6-B-Dt"= "#E0FFFF",
                   "GhCESA6-C-At"= "#FFF0F5",
                   "GhCESA6-C-Dt"= "#F0E68C",
                   "GhCESA8-B-At"= "#FFFFF0",
                   "GhCESA8-B-Dt"= "#CD5C5C",
                   "GhCESA6-Fa-Dt"="#FF69B4")


CesAplot1 <- ggplot(CesAexpr, aes(x=DPA,y=expression,group=CesA)) +
  geom_smooth(aes(color = CesA),
              alpha = 0.2) +
  scale_color_manual(values=CesAcolorstmp)

smoothed <- ggplot_build(CesAplot1)$data[[1]] %>%
  group_by(group) %>%
  slice(n()) %>%
  ungroup() %>%
  select(colour, group, y) %>%
  left_join(.,as.data.frame(CesAcolorstmp) %>% rownames_to_column(var="CesA"), by=c("colour" = "CesAcolorstmp")) %>%
  select(CesA,y) %>%
  mutate(expression = y, DPA=24) 
  
  
labels <- filter(CesAexpr, DPA==max(DPA)) %>%
  arrange(CWgroup) %>%
  filter(!is.na(label)) %>%
  select(-expression) %>% 
  left_join(.,smoothed) 

CesAplot <- ggplot(CesAexpr, aes(x=DPA,y=expression,group=CesA)) +
  geom_smooth(aes(color = stage),
              alpha = 0.2) +
  xlim(5,40) +
  ylim(0,15) +  
  geom_text_repel(aes(label = CesA, color=stage), labels, direction="y", nudge_x=10) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size=18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 18))  +
  facet_wrap(~CWgroup) + 
  scale_color_manual(values=c("brown", "darkgreen")) +
  scale_x_continuous(breaks = c(5, 10, 15, 20, 25)) +
  theme(legend.position = "none")

ggsave("out-12.CesAExpressionHomoeolog.jpg", plot=CesAplot, width = 15, height = 10, units = "in")

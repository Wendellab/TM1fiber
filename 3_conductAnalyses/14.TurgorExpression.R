library(RColorBrewer)
library(data.table)
library(tidyverse)

options(scipen = 999)
set.seed(8675309)
setDTthreads(18) # setting the data.table threads
options(datatable.fread.datatable=FALSE)


load(file = "04.R-01-dataInput.filtered.RData")

turgorTable <- fread('Turgor.tsv')


exprNorm <- t(datExpr) %>% 
  as.data.frame() %>% 
  rownames_to_column("shared_name") %>% 
  right_join(.,turgorTable)

Turgorexpr <- exprNorm %>%
  pivot_longer(!(c(shared_name,name,DEtrend,involvement,cluster)),names_to="Sample", values_to="expression") %>% 
  left_join(.,coldata[,c("id","DPA")],by=c("Sample"="id")) %>%
  mutate(DPA = as.numeric(sub("DPA","",DPA))) %>%
  filter(! DEtrend=="flat")
  

keep <- c("Gorai.003G074400.1.D","Gorai.002G198900.1.D","Gorai.006G181300.1.A","Gorai.003G074400.1.A","Gorai.010G030700.1.A","Gorai.010G147100.1.A")

Turgorexpr <- Turgorexpr %>%
  filter(shared_name %in% keep)

Turgorplot <- ggplot(Turgorexpr, aes(x=DPA,y=expression,group=name)) +
  geom_smooth(aes(color = name),
              alpha = 0.2) +
  xlim(5,25) +
  ylim(0,15) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~shared_name) 


ggsave("out-14.TurgorExpressionHomoeolognamed.jpg", plot=Turgorplot, width = 20, height = 20, units = "in")


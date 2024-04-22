library(data.table)
library(tidyverse)
library(dplyr)

options(scipen = 999)
set.seed(8675309)
setDTthreads(12) 
options(datatable.fread.datatable=FALSE)



# get TPM
tpm <- fread("tpmMean.Purdue.pseudoAD1.clean.tsv") %>% 
  rename(shared_name = target_id) %>%
  select(-AD1_TM1_25DPA_tpm) %>%
  rename_all(~gsub("AD1_TM1_(\\d{2})DPA_tpm", "DPA\\1", .))


# get impulse category
id <- scan("out-03.transient_down.list", what="character")
iu <- scan("out-03.transient_up.list", what="character")
td <- scan("out-03.transition_down.list", what="character")
tu <- scan("out-03.transition_up.list", what="character")


DEtrend <- bind_rows(
  data.frame(shared_name=id, DEtrend = "transient down"), 
  data.frame(shared_name=iu, DEtrend = "transient up"),
  data.frame(shared_name=td, DEtrend = "transition down"), 
  data.frame(shared_name=tu, DEtrend = "transition up")) 

# get module
me <- fread("out-07.s3.module_annotation.txt") %>%
  select(homoeolog,ME_rld,tair10.defline) %>%
  rename(shared_name=homoeolog, ME=ME_rld, name=tair10.defline) %>% 
  mutate(name = na_if(name,"")) %>%
  mutate(name = ifelse(str_detect(name, "^Protein of unknown"), NA, name)) %>%
  mutate(name = ifelse(str_detect(name, "^Domain of unknown"), NA, name)) %>%
  mutate(name = ifelse(str_detect(name, "^Eukaryotic protein of unknown"), NA, name)) %>%
  mutate(name = ifelse(str_detect(name, "^Expressed protein"), NA, name)) %>%
  mutate(name = ifelse(str_detect(name, "^DOMAIN OF UNKNOWN"), NA, name)) %>%
  mutate(name = ifelse(str_detect(name, "^Arabidopsis thaliana protein of unknown"), NA, name)) %>%
  mutate(name = ifelse(str_detect(name, "^unknown"), NA, name)) %>%
  mutate(name = ifelse(str_detect(name, "^Uncharacterised"), NA, name)) %>%
  mutate(name = gsub("\\\\\\\\","",name)) %>%
  mutate(name = gsub("$","-like",name))

# get clusters
cluster <- fread("ME.Louvain.Infomap.groups.tsv") %>%
  select(-ME_rld) %>%
  rename(shared_name=homoeolog, cluster=group)

# get TFs
tfa <- fread("GoraiTFs.txt") %>%
  mutate(shared_name=gsub("00$","00.1.A",Gene_ID)) %>%
  select(shared_name,Family)

tfd <- tfa %>%
  mutate(shared_name=gsub(".1.A",".1.D",shared_name))

tf <- bind_rows(tfa,tfd) %>%
  mutate(Family=gsub("$"," (TF)",Family))

# wall genes from Siva
walla <- fread("WallGenesSiva.txt") %>%
  mutate(shared_name=gsub("00$","00.1.A",shared_name)) 

walld <- walla %>%
  mutate(shared_name=gsub(".1.A",".1.D",shared_name))

wall <- bind_rows(walla,walld) %>%
  rename(wallname=name)

# genes from Alex
csi <- fread("CSIandMTgrowthAlex.txt") %>%
  rename(csiname=name) %>%
  mutate(csiname=gsub("$","-like",csiname))

# names from previous list
other <- fread("otherNames.txt") %>%
  rename(other=name)

# merge df
df <- Reduce(function(x, y) merge(x, y, by = "shared_name", all = TRUE), list(tpm, DEtrend, me, cluster,tf,wall,csi,other)) %>%
  mutate(DEtrend = replace_na(DEtrend, "flat")) %>%
  mutate(involvement=coalesce(involvement.x, involvement.y)) %>%
  mutate(name = coalesce(Family,wallname,csiname,name,other,shared_name)) %>%
  select(shared_name,name,DEtrend,involvement,cluster,ME,Louvain,Infomap,everything()) %>%
  select(-c(wallname,Family,csiname,involvement.x,involvement.y,other))  %>%
  replace_na(list(involvement="unknown",cluster="none",ME=99999,Louvain=99999,Infomap=99999))
  



fwrite(df,file="NodeNames.tsv",quote=F,sep="\t")

# Assuming your dataframe is named 'your_dataframe'
# and it has 75000 rows

# Create a vector of indices to split the dataframe
indices <- seq(1, nrow(df), by = 10000)

# Split the dataframe into a list of smaller dataframes
nodechunks <- split(df, findInterval(1:nrow(df), indices, rightmost.closed = TRUE))

# Save each chunk as a separate file
for (i in seq_along(nodechunks)) {
  filename <- paste0("NodeNameschunk_", i, ".tsv")
  write.table(nodechunks[[i]], file = filename, row.names = FALSE, sep="\t", quote=F)
}

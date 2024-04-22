library(data.table)
library(tidyverse)

setDTthreads(20)
options(datatable.fread.datatable=F)
set.seed(8675309)

x <- fread("aggregated.bb.p05.tsv",header=F)
names(x) <- c("source","target","directed","aracne","clr","genie3","llr","mi","narromi","pcor","pearson","plsnet","spearman","svm","tigress","wgcna","IRP","NC_Score")

#separate scores from ranks
y <- x %>%
	separate_wider_delim(aracne, ";", names = c("aracne_score", "aracne_rank")) %>%
	separate_wider_delim(clr, ";", names = c("clr_score", "clr_rank")) %>%
	separate_wider_delim(genie3, ";", names = c("genie3_score", "genie3_rank")) %>%
	separate_wider_delim(llr, ";", names = c("llr_score", "llr_rank")) %>%
	separate_wider_delim(mi, ";", names = c("mi_score", "mi_rank")) %>%
	separate_wider_delim(narromi, ";", names = c("narromi_score", "narromi_rank")) %>%
	separate_wider_delim(pcor, ";", names = c("pcor_score", "pcor_rank")) %>%
	separate_wider_delim(pearson, ";", names = c("pearson_score", "pearson_rank")) %>%
	separate_wider_delim(plsnet, ";", names = c("plsnet_score", "plsnet_rank")) %>%
	separate_wider_delim(spearman, ";", names = c("spearman_score", "spearman_rank")) %>%
	separate_wider_delim(svm, ";", names = c("svm_score", "svm_rank")) %>%
	separate_wider_delim(tigress, ";", names = c("tigress_score", "tigress_rank")) %>%
	separate_wider_delim(wgcna, ";", names = c("wgcna_score", "wgcna_rank")) %>%
	separate_wider_delim(IRP, ";", names = c("IRP_score", "IRP_rank")) %>%
	separate_wider_delim(NC_Score, ";", names = c("NC_score", "NC_sdev"))

# number of total edges that are top 10%
nrow(y) * 0.1
# [1] 222490

# number of total Directed or Undirected edges
y %>% filter(directed == "Directed") %>% nrow
# [1] 721101
y %>% filter(directed == "Undirected") %>% nrow
# [1] 1503799

# number of Directed or Undirected edges in the top 10%
y %>% arrange(desc(IRP_rank)) %>% head(222490) %>% filter(directed=="Directed") %>% nrow
# [1] 72545 (32.6%)
y %>% arrange(desc(IRP_rank)) %>% head(222490) %>% filter(directed=="Undirected") %>% nrow
# [1] 149945 (67.4%)

# dataframes for the top 10% of directed edges, undirected edges, or both edges
d <- y %>% arrange(desc(IRP_rank)) %>% filter(directed=="Directed") %>% head(72110) 
u <- y %>% arrange(desc(IRP_rank)) %>% filter(directed=="Undirected") %>% head(150380)
b <- y %>% arrange(desc(IRP_rank)) %>% head(222490)

write.table(d, file="out-11.directed.top.10pct.tsv", row.names=F, sep="\t", quote=F)
write.table(u, file="out-11.undirected.top.10pct.tsv", row.names=F, sep="\t", quote=F)
write.table(b, file="out-11.both_directed.top.10pct.tsv", row.names=F, sep="\t", quote=F)


save(x, y, d, u, b, file = "11.seidr.bb.p05.RData")



# do this again, but exclude homoeolog - homoeolog edges
# how many edges are not, or are, between homoeologs
y  %>%
	filter(substring(source, 1, nchar(source) - 2) != substring(target, 1, nchar(target) - 2)) %>% 
	nrow
# [1] 2214139

y  %>%
	filter(substring(source, 1, nchar(source) - 2) == substring(target, 1, nchar(target) - 2)) %>% 
	nrow
# [1] 10761 (5%)

# dataframe excluding homoeolog-homoeolog edges
z <- y  %>%
	filter(substring(source, 1, nchar(source) - 2) != substring(target, 1, nchar(target) - 2))

# number of total Directed or Undirected edges
z %>% filter(directed == "Directed") %>% nrow
# [1] 715679 (a reduction from 721101; 99.3%)
z %>% filter(directed == "Undirected") %>% nrow
# [1] 1498460 (a reduction from 1503799; 99.6%)

# number of Directed or Undirected edges in the top 10%
z %>% arrange(desc(IRP_rank)) %>% head(221414) %>% filter(directed=="Directed") %>% nrow
# [1] 71971 (32.5%)
z %>% arrange(desc(IRP_rank)) %>% head(221414) %>% filter(directed=="Undirected") %>% nrow
# [1] 149443 (67.5%)

# dataframes for the top 10% of directed edges, undirected edges, or both edges
dz <- z %>% arrange(desc(IRP_rank)) %>% filter(directed=="Directed") %>% head(71568) 
uz <- z %>% arrange(desc(IRP_rank)) %>% filter(directed=="Undirected") %>% head(149846)
bz <- z %>% arrange(desc(IRP_rank)) %>% head(221414)

write.table(dz, file="out-11.directed.top.10pct.nohomoeo.tsv", row.names=F, sep="\t", quote=F)
write.table(uz, file="out-11.undirected.top.10pct.nohomoeo.tsv", row.names=F, sep="\t", quote=F)
write.table(bz, file="out-11.both_directed.top.10pct.nohomoeo.tsv", row.names=F, sep="\t", quote=F)


library(igraph)

# Louvain Method: The Louvain algorithm is a popular method for community detection 
# that optimizes modularity. It iteratively assigns nodes to communities to maximize 
# the modularity score, which quantifies the strength of the community structure in the network.

# Create an igraph graph object from the edge list
edges <- y %>% select(c(source,target,IRP_score))
graph <- graph_from_data_frame(edges, directed = FALSE)

# Perform community detection using the Louvain algorithm
communities <- cluster_louvain(graph, weights = E(graph)$IRP_score)

# Get the membership vector indicating which community each node belongs to
membership <- membership(communities)

table(membership)
# membership
#    1     2     3     4     5     6     7     8     9    10    11    12    13 
#17576  6595  6588  8628   862  5379  2827  4868  7463  2115   115   676  1166 
#   14    15    16    17    18    19    20    21    22    23    24    25    26 
#    2     6     2   171     4    61     4   338    72     2     2    21     2 
#   27    28    29    30    31    32    33    34    35    36    37    38    39 
#    2     5     2     2     5    13     2     2     2     2     2     2     4 
#   40    41    42    43    44    45    46    47    48    49    50    51    52 
#   11     2     2     2     2     3     2     2     2     3     2     2     2 
#   53    54    55    56    57    58    59    60    61    62    63    64    65 
#    2     2     2     2     3     2     2     2     2     2     2     2     2 
#   66    67    68    69    70    71    72    73    74    75    76    77    78 
#    2     2     3     3     2     5     2     3     2     5     2    10     4 
#   79    80    81    82    83    84    85    86    87    88    89    90    91 
#    2     2     2     2     2     2     2     2     2     7     2     2     2 
#   92    93    94    95    96    97    98    99   100   101   102   103   104 
#    2     2     2     2     2     2     2     2     2     2     2     2     2 
#  105   106   107   108   109   110   111   112   113   114   115   116   117 
#    2     2     2     2     2     2     2     2     2     4     2     2     2 
#  118   119   120   121   122   123   124   125   126   127   128   129   130 
#    2     2     2     2     2     2     4     2     2     2     2     2     2 
#  131   132   133   134   135   136   137   138   139   140   141   142   143 
#    2     2     2     2     2     2     2     2     2     2     2     2     2 
#  144   145   146   147   148   149   150   151   152   153   154   155   156 
#    2     2     2     2     2     2     2     2     2     2     2     2     2 
#  157   158   159   160   161   162   163   164   165   166   167   168   169 
#    2     3     5     2     3     2     2     3     4     2     2     2     2 
#  170   171   172   173   174   175   176   177   178   179   180   181   182 
#    2     2     3     2     2     2     2     2     2     2     2     2     2 
#  183   184   185   186   187   188 



# Infomap: Infomap is an algorithm that uses information theory principles 
# to detect communities in networks. It considers network structure and 
# optimizes the description length of random walks in the network to identify distinct modules.


# Perform community detection using the Infomap algorithm
icommunities <- cluster_infomap(graph, e.weights = E(graph)$weight)

# Get the membership vector indicating which community each node belongs to
imembership <- membership(icommunities)

# Print the membership vector
table(imembership)

# imembership
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16
#    6 1796 4834 3350  479 2149 1438  132   35 1523 2435 1149 1413 1983 1078 2199
#   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32
#   10 1461   37 1822 1736  935    5   49   24   26  177   51  934    8  485  107
#   33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  147   52 1018   21  952  155   42   22   10   23  362    8   98    5  368    9
#   49   50   51   52   53   54   55   56   57   58   59   60   61   62   63   64
#  951    4  714   17  141  391   90   53   14  202   23   82   23    6   34   46
#   65   66   67   68   69   70   71   72   73   74   75   76   77   78   79   80
#  223   28  187    8   35   36   37   10   11   16  101   29   20  168   11   31
#   81   82   83   84   85   86   87   88   89   90   91   92   93   94   95   96
#   10   16   16   18  172   12   10  309   15   26   14   23  542   25   10  134
#   97   98   99  100  101  102  103  104  105  106  107  108  109  110  111  112
#   17   25   31   21   10   69  701   21    7   86   20    9    8   21   10    7
#  113  114  115  116  117  118  119  120  121  122  123  124  125  126  127  128
#   12   80   13   15   20    8   11   20  252   10    9   10   18   14   25  341
#  129  130  131  132  133  134  135  136  137  138  139  140  141  142  143  144
#   80   18   26   10  216   41  143    9   24   66    7   96    5   18   11    5
#  145  146  147  148  149  150  151  152  153  154  155  156  157  158  159  160
#   99   19    9    5   21   21   12    9   11   19    8   14   19    9   14   32
#  161  162  163  164  165  166  167  168  169  170  171  172  173  174  175  176
#   19   30   14   17   20   12   40   14   11   25   52   17   18   11   13   10
#  177  178  179  180  181  182  183  184  185  186  187  188  189  190  191  192
#   11   48   15   28    5   47   12    9  106   13   48   46   54   12    9   22
#  193  194  195  196  197  198  199  200  201  202  203  204  205  206  207  208
#   27   12   18   14    8   15  120   14   14   14   10   28  178   96    9   85
#  209  210  211  212  213  214  215  216  217  218  219  220  221  222  223  224
#   16   13   13   18   11   19   41   11    8   10   61   85   20   32    8   36
#  225  226  227  228  229  230  231  232  233  234  235  236  237  238  239  240
#    2   16    5    8    3   22   39   17   16   89   20   14   26   85   26   23
#  241  242  243  244  245  246  247  248  249  250  251  252  253  254  255  256
#   21    2    9   32   83    9   23   12   11    8   26    7   17   19   12   13
#  257  258  259  260  261  262  263  264  265  266  267  268  269  270  271  272
#   11    2   84   12   23   17   14   51    8   16   50    7   26   20    9   12
#  273  274  275  276  277  278  279  280  281  282  283  284  285  286  287  288
#   18   10    5   14   11   12   15   11   22   29   32    4   11   13    7  110
#  289  290  291  292  293  294  295  296  297  298  299  300  301  302  303  304
#   98   44   18   10    6    5   12   20    5   12   19   13   58   49   12  123
#  305  306  307  308  309  310  311  312  313  314  315  316  317  318  319  320
#   66   17   13   23   19   16    9   14   19    4   27    6   14   53   18   11
#  321  322  323  324  325  326  327  328  329  330  331  332  333  334  335  336
#   68    8   11   15   10   12   44   13   29   12   12    9    4   11   12   12
#  337  338  339  340  341  342  343  344  345  346  347  348  349  350  351  352
#   16   26    7   15   19   18   13   61    5   11   16   11   19   24   19   11
# ...
# 1953 1954 1955 1956 1957 1958 1959 1960 1961 1962 1963 1964 1965 1966 1967 1968
#    2    2    2    2    2    2    2    2    2    3    3    2    2    2    2    3
# 1969 1970 1971 1972 1973 1974 1975 1976
#    2    2    2    2    2    2    2    2

save.image("clustered.edges.R")

# read in wgcna gene-module associations
wgcna <- fread("../s3.module&annotation.txt", header=T)

# get module-gene only
modules <- wgcna %>% select(homoeolog,ME_rld)

comm <- cbind(communities$names,communities$membership) %>% 
	as.data.frame %>%
	rename(homoeolog=V1, Louvain=V2)

icomm <- cbind(icommunities$names,icommunities$membership) %>% 
	as.data.frame %>%
	rename(homoeolog=V1, Infomap=V2)


joind <- left_join(modules, comm, by="homoeolog") %>% 
	left_join(., icomm, by="homoeolog") %>% 
	filter(ME_rld != 0) %>%
	na.omit()

uniqueN(joind, by=c("ME_rld","Louvain", "Infomap"))
# [1] 6519

groupings <- joind %>% 
	mutate(group = factor(paste0(ME_rld,",",Louvain,",",Infomap)))


sumgroup <- groupings %>% 
	group_by(group) %>% 
	summarise(membership = n() ) %>%
	arrange(desc(membership))

write.table(groupings, file="out-11.ME.Louvain.Infomap.groups.tsv",quote=F, row.names=F, sep="\t")


TF <- fread("GoraiTFs.txt")

TFs <- rbind(TF %>% mutate(Gene_ID = gsub("$",".1.A",Gene_ID)),TF %>% mutate(Gene_ID = gsub("$",".1.D",Gene_ID))) %>%
	rename(homoeolog = Gene_ID)

TFbb <- edges %>% 
	filter(source %in% TFs$homoeolog | target %in% TFs$homoeolog)

TFgroupings <- groupings %>%
	filter(homoeolog %in% TFs$homoeolog)

TFsumgroup <- TFgroupings %>% 
	group_by(group) %>% 
	summarise(membership = n() ) %>%
	arrange(desc(membership))

summaryTable <- data.frame(group=character(0),module=character(0),nodes=numeric(0),edges=numeric(0),TFanodes=numeric(0),TFaedges=numeric(0),TFnodes=numeric(0),TFedges=numeric(0))


for (group.no in unique(groupings$group)) {
	# get homoeologs in group
	lx <- groupings %>%
		filter(group==group.no) %>%
		pull(homoeolog)

	# get module number from group
	ly <- gsub("-.*","",group.no)

	# get TF homoeologs and links
	tx <- groupings %>%
		filter(group==group.no) %>%
		filter(homoeolog %in% TFbb$source | homoeolog %in% TFbb$target) %>%
		pull(homoeolog)

	# get TF-TF links
	ttx <- lx[lx %in% TFs$homoeolog]

	# if more than 9 nodes, write file
	if (length(lx) > 1) {
		group.e <- edges %>%
			filter(source %in% lx) %>%
			filter(target %in% lx)
		tgroup.e <- TFbb %>%
			filter(source %in% tx & target %in% tx) 
		ttgroup.e <- tgroup.e %>%
			filter(source %in% ttx & target %in% ttx)

		nedges <- nrow(group.e)
		tnedges <- nrow(tgroup.e)
		ttnedges <- nrow(ttgroup.e)

		nnodes <- length(unique(c(group.e$source,group.e$target)))
		tnnodes <- length(unique(c(tgroup.e$source,tgroup.e$target)))
		ttnnodes <- length(unique(c(ttgroup.e$source,ttgroup.e$target)))

		summaryTable <- summaryTable %>%
			add_row(group=group.no,module=paste0("ME",ly),nodes=nnodes,edges=nedges,TFanodes=tnnodes,TFaedges=tnedges,TFnodes=ttnnodes,TFedges=ttnedges)
		
		if (nnodes > 1 & nedges > 0) {
			write.table(group.e, file=paste0("grouped.edges/network.from.group",group.no,".from.module",ly,".with.",nnodes,"nodes.and.",nedges,"edges.tsv"), sep="\t",quote=F, row.names=F)
		}
		if (tnnodes > 1 & tnedges > 0) {
			write.table(tgroup.e, file=paste0("grouped.edges.TF/network.from.group",group.no,".from.module",ly,".with.",tnnodes,"nodes.and.",tnedges,"edges.tsv"), sep="\t",quote=F, row.names=F)
		}
		if (ttnnodes > 1 & ttnedges > 0) {
			write.table(ttgroup.e, file=paste0("grouped.edges.TF.only/network.from.group",group.no,".from.module",ly,".with.",ttnnodes,"nodes.and.",ttnedges,"edges.tsv"), sep="\t",quote=F, row.names=F)
		}
	}
}



write.table(summaryTable, file="ME.Louvain.Infomap.group.summary.tsv",quote=F, row.names=F, sep="\t")


############# top 10% all directed edges ##############


xd <- fread("aggregated.directed.edges.top10pct.tsv",header=F)
names(xd) <- c("source","target","IRP_score", "IRP_rank")

# number of total edges that are top 10%
nrow(xd) * 0.1
# [1] 849518.3


# dataframes for the top 10% of directed edges, undirected edges, or both edges
dd <- xd %>% arrange(desc(IRP_rank)) %>% head(849518) 

write.table(dd, file="directed.allDirected.top.10pct.tsv", row.names=F, sep="\t", quote=F)

# do this again, but exclude homoeolog - homoeolog edges
# how many edges are not, or are, between homoeologs
xd  %>%
	filter(substring(source, 1, nchar(source) - 2) != substring(target, 1, nchar(target) - 2)) %>% 
	nrow
# [1] 8487853

xd  %>%
	filter(substring(source, 1, nchar(source) - 2) == substring(target, 1, nchar(target) - 2)) %>% 
	nrow
# [1] 7330 (0.09%)

# dataframe excluding homoeolog-homoeolog edges
zd <- xd  %>%
	filter(substring(source, 1, nchar(source) - 2) != substring(target, 1, nchar(target) - 2))


# dataframes for the top 10% of directed edges, undirected edges, or both edges
dzd <- zd %>% arrange(desc(IRP_rank)) %>% head(84879) 

write.table(dzd, file="directedOnly.top.10pct.nohomoeo.tsv", row.names=F, sep="\t", quote=F)


library(igraph)

# Louvain Method: The Louvain algorithm is a popular method for community detection 
# that optimizes modularity. It iteratively assigns nodes to communities to maximize 
# the modularity score, which quantifies the strength of the community structure in the network.

# Create an igraph graph object from the edge list
edgesD <- xd %>% select(c(source,target,IRP_score))
graphD <- graph_from_data_frame(edgesD, directed = FALSE)

# Perform community detection using the Louvain algorithm
communitiesD <- cluster_louvain(graphD, weights = E(graphD)$IRP_score)

# Get the membership vector indicating which community each node belongs to
membershipD <- membership(communitiesD)

table(membershipD)
membershipD
#    1     2     3     4     5
#35619   483 32735  5315   252




# Infomap: Infomap is an algorithm that uses information theory principles 
# to detect communities in networks. It considers network structure and 
# optimizes the description length of random walks in the network to identify distinct modules.


# Perform community detection using the Infomap algorithm
icommunitiesD <- cluster_infomap(graphD, e.weights = E(graphD)$weight)

# Get the membership vector indicating which community each node belongs to
imembershipD <- membership(icommunitiesD)

# Print the membership vector
table(imembershipD)

imembershipD
#    1     2     3     4     5     6     7     8     9    10    11    12    13
#11547   807    10   807  4003     6   513  1964    21     9  1341    15    13
#   14    15    16    17    18    19    20    21    22    23    24    25    26
#    6  1964    10    12    43    12     8     9     2    14    10   890     9
#   27    28    29    30    31    32    33    34    35    36    37    38    39
# 1857    10    12  6650    14  3130   973   111    14   759    13  3538   570
#  599   600   601   602   603   604   605   606   607   608   609   610   611
#    6     5     5     5     8     6     3     5     8     2     2     2     8
#  612   613   614   615   616   617   618   619   620   621   622   623   624
#    3     2     4     4     2     5     5     2     3     7    12     3     2
#  625   626


save.image("clustered.directed.edges.Rdata")

# read in wgcna gene-module associations

# get module-gene only
modules <- wgcna %>% select(homoeolog,ME_rld)

commD <- cbind(communitiesD$names,communitiesD$membership) %>% 
	as.data.frame %>%
	rename(homoeolog=V1, Louvain=V2)

icommD <- cbind(icommunitiesD$names,icommunitiesD$membership) %>% 
	as.data.frame %>%
	rename(homoeolog=V1, Infomap=V2)


joindD <- left_join(modules, commD, by="homoeolog") %>% 
	left_join(., icommD, by="homoeolog") %>% 
	filter(ME_rld != 0) %>%
	na.omit()

uniqueN(joindD, by=c("ME_rld","Louvain", "Infomap"))
# [1] 2206

groupingsD <- joindD %>% 
	mutate(group = factor(paste0(ME_rld,"-",Louvain,"-",Infomap)))


sumgroupD <- groupingsD %>% 
	group_by(group) %>% 
	summarise(membership = n() ) %>%
	arrange(desc(membership))

write.table(groupingsD, file="ME.Louvain.Infomap.groups.Directed.tsv",quote=F, row.names=F, sep="\t")


TFD <- edgesD %>% 
	filter(source %in% TFs$homoeolog | target %in% TFs$homoeolog)

TFgroupingsD <- groupingsD %>%
	filter(homoeolog %in% TFs$homoeolog)

TFsumgroupD <- TFgroupingsD %>% 
	group_by(group) %>% 
	summarise(membership = n() ) %>%
	arrange(desc(membership))

summaryTableD <- data.frame(group=character(0),module=character(0),nodes=numeric(0),edges=numeric(0),TFanodes=numeric(0),TFaedges=numeric(0),TFnodes=numeric(0),TFedges=numeric(0))


for (group.no in unique(groupingsD$group)) {
	# get homoeologs in group
	lx <- groupingsD %>%
		filter(group==group.no) %>%
		pull(homoeolog)

	# get module number from group
	ly <- gsub("-.*","",group.no)

	# get TF homoeologs and links
	tx <- groupingsD %>%
		filter(group==group.no) %>%
		filter(homoeolog %in% TFD$source | homoeolog %in% TFD$target) %>%
		pull(homoeolog)

	# get TF-TF links
	ttx <- lx[lx %in% TFs$homoeolog]

	# if more than 9 nodes, write file
	if (length(lx) > 1) {
		group.e <- edgesD %>%
			filter(source %in% lx) %>%
			filter(target %in% lx)
		tgroup.e <- TFD %>%
			filter(source %in% tx & target %in% tx) 
		ttgroup.e <- tgroup.e %>%
			filter(source %in% ttx & target %in% ttx)

		nedges <- nrow(group.e)
		tnedges <- nrow(tgroup.e)
		ttnedges <- nrow(ttgroup.e)

		nnodes <- length(unique(c(group.e$source,group.e$target)))
		tnnodes <- length(unique(c(tgroup.e$source,tgroup.e$target)))
		ttnnodes <- length(unique(c(ttgroup.e$source,ttgroup.e$target)))

		summaryTableD <- summaryTableD %>%
			add_row(group=group.no,module=paste0("ME",ly),nodes=nnodes,edges=nedges,TFanodes=tnnodes,TFaedges=tnedges,TFnodes=ttnnodes,TFedges=ttnedges)
		
		if (nnodes > 1 & nedges > 0) {
			write.table(group.e, file=paste0("grouped.edges/network.from.group",group.no,".from.module",ly,".with.",nnodes,"nodes.and.",nedges,"edges.tsv"), sep="\t",quote=F, row.names=F)
		}
		if (tnnodes > 1 & tnedges > 0) {
			write.table(tgroup.e, file=paste0("grouped.edges.TF/network.from.group",group.no,".from.module",ly,".with.",tnnodes,"nodes.and.",tnedges,"edges.tsv"), sep="\t",quote=F, row.names=F)
		}
		if (ttnnodes > 1 & ttnedges > 0) {
			write.table(ttgroup.e, file=paste0("grouped.edges.TF.only/network.from.group",group.no,".from.module",ly,".with.",ttnnodes,"nodes.and.",ttnedges,"edges.tsv"), sep="\t",quote=F, row.names=F)
		}
	}
}



write.table(summaryTable, file="ME.Louvain.Infomap.group.summary.tsv",quote=F, row.names=F, sep="\t")



############### Directed bb edges #######################

library(igraph)

# Louvain Method: The Louvain algorithm is a popular method for community detection 
# that optimizes modularity. It iteratively assigns nodes to communities to maximize 
# the modularity score, which quantifies the strength of the community structure in the network.

# Create an igraph graph object from the edge list
edgesDB <- y %>% 
	filter(directed == "Directed") %>% 
	select(c(source,target,IRP_score))

graphDB <- graph_from_data_frame(edgesDB, directed = FALSE)

# Perform community detection using the Louvain algorithm
communitiesDB <- cluster_louvain(graphDB, weights = E(graphDB)$IRP_score)

# Get the membership vector indicating which community each node belongs to
membershipDB <- membership(communitiesDB)

table(membershipDB)
#membershipDB
#    1     2     3     4     5     6     7     8     9    10    11    12    13
# 3352  6165  1279  7265  5619  2405  6421  5705 12981  5080  1297   762     3
#   14    15    16    17    18    19    20    21    22    23    24    25    26
# 1608   162     2     2     2     2   280   216     3    60    95     7    54
#   27    28    29    30    31    32    33    34    35    36    37    38    39
#  445    58     2     2     2     2     5     2     2     2     8     2     2
#   40    41    42    43    44    45    46    47    48    49    50    51    52
#    2     3     2     2     2     2     3     2     8     2     2     6     2
#   53    54    55    56    57    58    59    60    61    62    63    64    65
#    2     6     6     3     2     2     2     2     2     2     2     2     2
#  274   275   276   277   278   279   280   281   282   283   284   285   286
#    4     2     2     2     3     2     2     2     2     2     2     2     2
#  287   288   289   290   291   292   293   294   295   296   297   298   299
#    2     2     2     2     2     2     2     3     2     2     2     5     3
#  300   301   302   303   304   305   306   307   308   309   310   311   312
#    4     2     3     2     2     3     2     2     2     2     2     2     2
#  313   314   315   316   317   318   319   320   321   322   323   324   325
#    2     2     2     2     2     3     2     2     2     2     4     2     2
#  326   327   328   329   330   331   332   333   334   335   336   337   338
#    2     2     2     2     2     2     2     2     2     2     2     2     2
#  339   340   341   342   343   344   345   346   347   348   349   350   351
#    2     2     3     2     4     2     3     3     2     2     2     3     2
#  352   353   354   355   356   357   358   359   360   361   362   363   364
#    2     2     2     2     2     2     2     2     2     2     2     3     2
#  365   366   367   368   369
#    2     2     2     2     2





# Infomap: Infomap is an algorithm that uses information theory principles 
# to detect communities in networks. It considers network structure and 
# optimizes the description length of random walks in the network to identify distinct modules.


# Perform community detection using the Infomap algorithm
icommunitiesDB <- cluster_infomap(graphDB, e.weights = E(graphDB)$weight)

# Get the membership vector indicating which community each node belongs to
imembershipDB <- membership(icommunitiesDB)

# Print the membership vector
table(imembershipDB)

# imembership
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16
#    6 1796 4834 3350  479 2149 1438  132   35 1523 2435 1149 1413 1983 1078 2199
#   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32
#   10 1461   37 1822 1736  935    5   49   24   26  177   51  934    8  485  107
#   33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  147   52 1018   21  952  155   42   22   10   23  362    8   98    5  368    9
#   49   50   51   52   53   54   55   56   57   58   59   60   61   62   63   64
#  951    4  714   17  141  391   90   53   14  202   23   82   23    6   34   46
#   65   66   67   68   69   70   71   72   73   74   75   76   77   78   79   80
#  223   28  187    8   35   36   37   10   11   16  101   29   20  168   11   31
#   81   82   83   84   85   86   87   88   89   90   91   92   93   94   95   96
#   10   16   16   18  172   12   10  309   15   26   14   23  542   25   10  134
#   97   98   99  100  101  102  103  104  105  106  107  108  109  110  111  112
#   17   25   31   21   10   69  701   21    7   86   20    9    8   21   10    7
#  113  114  115  116  117  118  119  120  121  122  123  124  125  126  127  128
#   12   80   13   15   20    8   11   20  252   10    9   10   18   14   25  341
#  129  130  131  132  133  134  135  136  137  138  139  140  141  142  143  144
#   80   18   26   10  216   41  143    9   24   66    7   96    5   18   11    5
#  145  146  147  148  149  150  151  152  153  154  155  156  157  158  159  160
#   99   19    9    5   21   21   12    9   11   19    8   14   19    9   14   32
#  161  162  163  164  165  166  167  168  169  170  171  172  173  174  175  176
#   19   30   14   17   20   12   40   14   11   25   52   17   18   11   13   10
#  177  178  179  180  181  182  183  184  185  186  187  188  189  190  191  192
#   11   48   15   28    5   47   12    9  106   13   48   46   54   12    9   22
#  193  194  195  196  197  198  199  200  201  202  203  204  205  206  207  208
#   27   12   18   14    8   15  120   14   14   14   10   28  178   96    9   85
#  209  210  211  212  213  214  215  216  217  218  219  220  221  222  223  224
#   16   13   13   18   11   19   41   11    8   10   61   85   20   32    8   36
#  225  226  227  228  229  230  231  232  233  234  235  236  237  238  239  240
#    2   16    5    8    3   22   39   17   16   89   20   14   26   85   26   23
#  241  242  243  244  245  246  247  248  249  250  251  252  253  254  255  256
#   21    2    9   32   83    9   23   12   11    8   26    7   17   19   12   13
#  257  258  259  260  261  262  263  264  265  266  267  268  269  270  271  272
#   11    2   84   12   23   17   14   51    8   16   50    7   26   20    9   12
#  273  274  275  276  277  278  279  280  281  282  283  284  285  286  287  288
#   18   10    5   14   11   12   15   11   22   29   32    4   11   13    7  110
#  289  290  291  292  293  294  295  296  297  298  299  300  301  302  303  304
#   98   44   18   10    6    5   12   20    5   12   19   13   58   49   12  123
#  305  306  307  308  309  310  311  312  313  314  315  316  317  318  319  320
#   66   17   13   23   19   16    9   14   19    4   27    6   14   53   18   11
#  321  322  323  324  325  326  327  328  329  330  331  332  333  334  335  336
#   68    8   11   15   10   12   44   13   29   12   12    9    4   11   12   12
#  337  338  339  340  341  342  343  344  345  346  347  348  349  350  351  352
#   16   26    7   15   19   18   13   61    5   11   16   11   19   24   19   11
# ...
# 1953 1954 1955 1956 1957 1958 1959 1960 1961 1962 1963 1964 1965 1966 1967 1968
#    2    2    2    2    2    2    2    2    2    3    3    2    2    2    2    3
# 1969 1970 1971 1972 1973 1974 1975 1976
#    2    2    2    2    2    2    2    2

save.image("clustered.edges.R")

# read in wgcna gene-module associations
wgcna <- fread("../s3.module&annotation.txt", header=T)

# get module-gene only
modules <- wgcna %>% select(homoeolog,ME_rld)

comm <- cbind(communities$names,communities$membership) %>% 
	as.data.frame %>%
	rename(homoeolog=V1, Louvain=V2)

icomm <- cbind(icommunities$names,icommunities$membership) %>% 
	as.data.frame %>%
	rename(homoeolog=V1, Infomap=V2)


joind <- left_join(modules, comm, by="homoeolog") %>% 
	left_join(., icomm, by="homoeolog") %>% 
	filter(ME_rld != 0) %>%
	na.omit()

uniqueN(joind, by=c("ME_rld","Louvain", "Infomap"))
# [1] 6519

groupings <- joind %>% 
	mutate(group = factor(paste0(ME_rld,"-",Louvain,"-",Infomap)))


sumgroup <- groupings %>% 
	group_by(group) %>% 
	summarise(membership = n() ) %>%
	arrange(desc(membership))

write.table(groupings, file="ME.Louvain.Infomap.groups.tsv",quote=F, row.names=F, sep="\t")


TF <- fread("../GoraiTFs.txt")

TFs <- rbind(TF %>% mutate(Gene_ID = gsub("$",".1.A",Gene_ID)),TF %>% mutate(Gene_ID = gsub("$",".1.D",Gene_ID))) %>%
	rename(homoeolog = Gene_ID)

TFbb <- edges %>% 
	filter(source %in% TFs$homoeolog | target %in% TFs$homoeolog)

TFgroupings <- groupings %>%
	filter(homoeolog %in% TFs$homoeolog)

TFsumgroup <- TFgroupings %>% 
	group_by(group) %>% 
	summarise(membership = n() ) %>%
	arrange(desc(membership))

summaryTable <- data.frame(group=character(0),module=character(0),nodes=numeric(0),edges=numeric(0),TFanodes=numeric(0),TFaedges=numeric(0),TFnodes=numeric(0),TFedges=numeric(0))


for (group.no in unique(groupings$group)) {
	# get homoeologs in group
	lx <- groupings %>%
		filter(group==group.no) %>%
		pull(homoeolog)

	# get module number from group
	ly <- gsub("-.*","",group.no)

	# get TF homoeologs and links
	tx <- groupings %>%
		filter(group==group.no) %>%
		filter(homoeolog %in% TFbb$source | homoeolog %in% TFbb$target) %>%
		pull(homoeolog)

	# get TF-TF links
	ttx <- lx[lx %in% TFs$homoeolog]

	# if more than 9 nodes, write file
	if (length(lx) > 1) {
		group.e <- edges %>%
			filter(source %in% lx) %>%
			filter(target %in% lx)
		tgroup.e <- TFbb %>%
			filter(source %in% tx & target %in% tx) 
		ttgroup.e <- tgroup.e %>%
			filter(source %in% ttx & target %in% ttx)

		nedges <- nrow(group.e)
		tnedges <- nrow(tgroup.e)
		ttnedges <- nrow(ttgroup.e)

		nnodes <- length(unique(c(group.e$source,group.e$target)))
		tnnodes <- length(unique(c(tgroup.e$source,tgroup.e$target)))
		ttnnodes <- length(unique(c(ttgroup.e$source,ttgroup.e$target)))

		summaryTable <- summaryTable %>%
			add_row(group=group.no,module=paste0("ME",ly),nodes=nnodes,edges=nedges,TFanodes=tnnodes,TFaedges=tnedges,TFnodes=ttnnodes,TFedges=ttnedges)
		
		if (nnodes > 1 & nedges > 0) {
			write.table(group.e, file=paste0("grouped.edges/network.from.group",group.no,".from.module",ly,".with.",nnodes,"nodes.and.",nedges,"edges.tsv"), sep="\t",quote=F, row.names=F)
		}
		if (tnnodes > 1 & tnedges > 0) {
			write.table(tgroup.e, file=paste0("grouped.edges.TF/network.from.group",group.no,".from.module",ly,".with.",tnnodes,"nodes.and.",tnedges,"edges.tsv"), sep="\t",quote=F, row.names=F)
		}
		if (ttnnodes > 1 & ttnedges > 0) {
			write.table(ttgroup.e, file=paste0("grouped.edges.TF.only/network.from.group",group.no,".from.module",ly,".with.",ttnnodes,"nodes.and.",ttnedges,"edges.tsv"), sep="\t",quote=F, row.names=F)
		}
	}
}



write.table(summaryTable, file="ME.Louvain.Infomap.group.summary.tsv",quote=F, row.names=F, sep="\t")







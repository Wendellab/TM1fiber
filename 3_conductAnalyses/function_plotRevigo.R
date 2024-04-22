library(treemapify)
library(scales)
library(ggplot2)

RevigoTreeMap <- function (TreeMapdf, me, type) {
  if (type=="BP") {
    what="Biological Processes, "
  } else if (type=="MF") {
    what="Molecular Function, "
  } else {
    what="error, "
  }
  
  if (me<10) {
    ME=paste0("ME0",me)
  } else if (me>9) {
    ME=paste0("ME0",me)
  } else {
    ME="error"
  }
  
  
  
  myplot <- ggplot(TreeMapdf, aes(area=value, fill=representative,label=description,subgroup=representative)) +
    ggtitle(paste0("     ", ME, ", ", what, length(universe[net$colors==me])," genes")) +
    geom_treemap(layout="squarified") + 
    geom_treemap_text(place = "centre",size = 12,reflow=T) 

  ggsave(filename=paste0("topGO/", ME, "_", type, "_wordmap.jpg"), device="jpeg", plot=myplot, width=12, height=12)
}
  
  
RevigoScatter <- function (scatterdf, me, type) {
  if (type=="BP") {
    what="Biological Processes, "
  } else if (type=="MF") {
    what="Molecular Function, "
  } else {
    what="error, "
  }

  if (me<10) {
    ME=paste0("ME0",me)
  } else if (me>9) {
    ME=paste0("ME0",me)
  } else {
    ME="error"
  }
  
  myplot <- ggplot( data = scatterdf ) + 
    ggtitle(paste0("     ", ME, ", ", what, length(universe[net$colors==me])," genes")) +
    geom_point( aes( plot_X, plot_Y, colour = value, size = log_size), alpha = I(0.6) ) +
    scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(scatterdf$value), 0) ) + 
    geom_point( aes(plot_X, plot_Y, size = log_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + 
    scale_size( range=c(5, 30)) + 
    theme_bw() + 
    geom_text( data = scatterdf [ scatterdf$dispensability < 0.15, ], aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 ) +
    labs (y = "semantic space x", x = "semantic space y") +
    theme(legend.key = element_blank()) +
    xlim(min(scatterdf$plot_X)-(max(scatterdf$plot_X) - min(scatterdf$plot_X))/10,max(scatterdf$plot_X)+(max(scatterdf$plot_X) - min(scatterdf$plot_X))/10) + 
    ylim(min(scatterdf$plot_Y)-(max(scatterdf$plot_Y) - min(scatterdf$plot_Y))/10,max(scatterdf$plot_Y)+(max(scatterdf$plot_Y) - min(scatterdf$plot_Y))/10) 

  ggsave(filename=paste0("topGO/", ME, "_", type, "_scatter.jpg"), device="jpeg", plot=myplot, width=12, height=12)
  
}  
###plotting qualtiy scores

library(ggplot2)
#install.packages("gridExtra")
library(gridExtra)
#install.packages("readr")
library(readr)

#######

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

split = set

# Function for making density plots of a single annotation
makeDensityPlot <- function(dataframe, xvar, split, xmin=min(dataframe[xvar], na.rm=TRUE), xmax=max(dataframe[xvar], na.rm=TRUE), alpha=0.5) {
  
  if(missing(split)) {
    return(ggplot(data=dataframe, aes_string(x=xvar)) + xlim(xmin,xmax) + geom_density() )
  }
  else {
    return(ggplot(data=dataframe, aes_string(x=xvar, fill=split)) + xlim(xmin,xmax) + geom_density(alpha=alpha) )
  }
}

# Function for making scatter plots of two annotations
makeScatterPlot <- function(dataframe, xvar, yvar, split, xmin=min(dataframe[xvar], na.rm=TRUE), xmax=max(dataframe[xvar], na.rm=TRUE), ymin=min(dataframe[yvar], na.rm=TRUE), ymax=max(dataframe[yvar], na.rm=TRUE), ptSize=1, alpha=0.6) {
  if(missing(split)) {
    return(ggplot(data=dataframe) + aes_string(x=xvar, y=yvar) + xlim(xmin,xmax) + ylim(ymin,ymax) + geom_point(size=ptSize, alpha=alpha) )
  }
  else {
    return(ggplot(data=dataframe) + aes_string(x=xvar, y=yvar) + aes_string(color=split) + xlim(xmin,xmax) + ylim(ymin,ymax) + geom_point(size=ptSize, alpha=alpha) )
  }
}

# Function for making scatter plots of two annotations with marginal density plots of each
makeScatterPlotWithMarginalDensity <- function(dataframe, xvar, yvar, split, xmin=min(dataframe[xvar], na.rm=TRUE), xmax=max(dataframe[xvar], na.rm=TRUE), ymin=min(dataframe[yvar], na.rm=TRUE), ymax=max(dataframe[yvar], na.rm=TRUE), ptSize=1, ptAlpha=0.6, fillAlpha=0.5) {
  empty <- ggplot()+geom_point(aes(1,1), colour="white") +
    theme(
      plot.background = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
    )
  
  if(missing(split)){
    scatter <- ggplot(data=dataframe) + aes_string(x=xvar, y=yvar) + geom_point(size=ptSize, alpha=ptAlpha) + xlim(xmin,xmax) + ylim(ymin,ymax) 
    plot_top <- ggplot(data=dataframe, aes_string(x=xvar)) + geom_density(alpha=fillAlpha) + theme(legend.position="none") + xlim(xmin,xmax) 
    plot_right <- ggplot(data=dataframe, aes_string(x=yvar)) + geom_density(alpha=fillAlpha) + coord_flip() + theme(legend.position="none") + xlim(ymin,ymax) 
  } 
  else{
    scatter <- ggplot(data=dataframe) + aes_string(x=xvar, y=yvar) + geom_point(size=ptSize, alpha=ptAlpha, aes_string(color=split)) + xlim(xmin,xmax) + ylim(ymin,ymax) 
    plot_top <- ggplot(data=dataframe, aes_string(x=xvar, fill=split)) + geom_density(alpha=fillAlpha) + theme(legend.position="none") + xlim(xmin,xmax) 
    plot_right <- ggplot(data=dataframe, aes_string(x=yvar, fill=split)) + geom_density(alpha=fillAlpha) + coord_flip() + theme(legend.position="none") + xlim(ymin,ymax) 
  }
  legend <- get_legend(scatter)
  scatter <- scatter + theme(legend.position="none")
  temp <- grid.arrange(plot_top, legend, scatter, plot_right, ncol=2, nrow=2, widths=c(4,1), heights=c(1,4))
  return(temp)
}



####my data


motherSNP.giab <- read_delim("motherSNP.giab.table","\t",escape_double = FALSE, col_types = cols(giab.callsets = col_character()),trim_ws = TRUE)

names(motherSNP.giab)[names(motherSNP.giab) == 'giab.callsets'] <- 'set'

motherSNP.giab$set <- as.character(motherSNP.giab$set)

qual = makeDensityPlot(motherSNP.giab, "QUAL")
qual
qual = makeDensityPlot(motherSNP.giab, "QUAL", xmax=10000, split="set")
qual

QD = makeDensityPlot(motherSNP.giab, "QD", xmax=100)
QD
QD = makeDensityPlot(motherSNP.giab, "QD", xmax=100,split="set")
QD

QD_DP = makeScatterPlot(motherSNP.giab, "QD", "DP", split = "set", ymax=1000)
QD_DP

QD_DP_md = makeScatterPlotWithMarginalDensity(motherSNP.giab, "QD", "DP", split = "set", ymax=1000)
QD_DP_md

MQrs= makeDensityPlot(motherSNP.giab, "QD", xmax=100)
MQrs

MQrs_set= makeDensityPlot(motherSNP.giab, "QD", xmax=100_split='set')
MQrs_set

MQrs= makeDensityPlot(motherSNP.giab, "MQRankSum", xmax=100)
MQrs

MQrs_set= makeDensityPlot(motherSNP.giab, "MQRankSum", xmax=20, split='set')
MQrs_set

StrandOddsRatio

SOR= makeDensityPlot(motherSNP.giab, "SOR", xmax=20)
SOR

SOR_set= makeDensityPlot(motherSNP.giab, "SOR", xmax=20, split='set')
SOR_set

MQrs_SOR_md = makeScatterPlotWithMarginalDensity(motherSNP.giab, "SOR", "MQRankSum", split = "set", ymax=10)
MQrs_SOR_md

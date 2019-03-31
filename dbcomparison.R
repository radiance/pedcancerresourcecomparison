# See ped-db-comparison.csv for 
# Downloaded data from 4 selected pediatric databases 
# median and standard deviation computation and plotting

csvdata <- read.csv('~/Documents/ped-db-comparison.csv', header=TRUE, sep=";", quote = "", stringsAsFactors=FALSE) #as.is = 1 
numdf = as.data.frame.matrix(csvdata[2:5])

# create list of medians and standard deviation
medianlist = list()
stdlist = list()

for (row in 1:nrow(numdf)) {
  x <- (numdf[row,])
  medianlist[row] <- median(unlist(x), na.rm = TRUE)
  stdlist[row] <- sd(unlist(x), na.rm = TRUE)
}

# plot medians and standard deviation as barplot
barplot(unlist(medianlist), main="Median Distribution across genes for selected dbs", xlab="median")
barplot(unlist(stdlist), main="Standard Deviation Distribution across genes for selected dbs", xlab="standard deviation")

cbiomed = list()
pedcbiomed = list()
pecanmed = list()
icgcmed = list()

for (row in 1:nrow(numdf)) {
  y <- (numdf[row,])
  
  cbiomed[row] <- abs(unlist(y[2]) - medianlist[[row]])
  pedcbiomed[row] <- abs(unlist(y[2]) - medianlist[[row]])
  icgcmed[row] <- abs(unlist(y[3]) - medianlist[[row]])
  pecanmed[row] <- abs(unlist(y[3]) - medianlist[[row]])
}

# plot each db's median deviation
xi <- 1:100
plot(xi, cbiomed, 
     xlabel="per gene", ylabel="deviation from median mutation", type="line", main = "Plot of median deviation against index of gene")
lines(xi, pedcbiomed, col="green")
lines(xi, icgcmed, col="red")
lines(xi, pecanmed, col="blue")


# Other Plotting examples
#plot(cbiomed,pedcbiomed, main="Scatterplot of cbio vs. pedcbio") 
#heatplot(csvdata)
#barplot(t(as.matrix(numdf)), beside=TRUE, main="Mutation allocation across dbs for top 100 genes", xlab="# Mutations")
#library(ggplot2)
#library(ggplotify)
#ggplot(data=numdf, x = Variable, y = Count)
#plot(numdf)
#matplot(numdf, numdf[], type="l")
#
#install.packages("circlize")
#devtools::install_github("jokergoo/circlize")
#library(circlize)
#chordDiagram(head(numdf))
#
#install.packages("plotly")
#library(plotly)
#packageVersion('plotly')
#plot_ly(csvdata, x = ~genes, y = ~mutation count, #text = ~gene, type = 'scatter', mode = 'markers' , marker = medianlist(opacity = 0.5, color = 'rgb(255, 65, 54)'), layout(title = 'Comparison across dbs',xaxis = list(showgrid = FALSE), yaxis = list(showgrid = FALSE)))
#devtools::install_github("ndphillips/yarrr", build_vignettes = TRUE)
#library("yarrr") 
#pirateplot(data=comparemedians)
#install.packages("cowplot")
#library("cowplot")
#cowplot(comparemedians)

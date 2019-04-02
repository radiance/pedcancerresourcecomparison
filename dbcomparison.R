# See ped-db-comparison.csv for 
# Downloaded data from 4 selected pediatric databases 
# median and standard deviation computation and plotting

# if you would like to use a path for windows: csvdata <- read.csv('file:///C:/Documents/ped-db-100mutations-export.csv', header=TRUE, sep=";", quote = "", stringsAsFactors=FALSE) #as.is = 1 
csvdata <- read.csv('~/Documents/ped-db-100mutations-export.csv', header=TRUE, sep=";", quote = "", stringsAsFactors=FALSE) #as.is = 1 
numdf = as.data.frame.matrix(csvdata[2:5])

# create list of medians, standard deviation and average
medianlist = list()
stdlist = list()
averagelist = list()

for (row in 1:nrow(numdf)) {
  x <- (numdf[row,])
  medianlist[row] <- median(unlist(x), na.rm = TRUE)
  stdlist[row] <- sd(unlist(x), na.rm = TRUE)
  averagelist[row] <- mean(unlist(x), na.rm = TRUE)
}


# plot avarage/mean, median and standard deviation as barplot
genelabels <- csvdata[1]
plot(unlist(stdlist), xlab="gene symbol", ylab="median, standard deviation and average from mutation frequencies", type="b", pch=2, col="gray8", frame.plot=FALSE, xaxt="n") #main="Median Distribution across genes for selected dbs"
axis(1, at=1:100, labels=unlist(genelabels), las=2, cex.axis=.5)
lines(unlist(medianlist), col="gray32", type="b", pch=5) #main="Standard Deviation Distribution across genes for selected dbs"
lines(unlist(averagelist), col="black", type="b", pch=9) #main="Average/mean across genes for selected dbs"
legend("topright", legend=c("median", "standard deviation ", "mean"), col=c(col="gray8", col="gray32", col="black"), pch=c(2,5,9), lty=2:2, cex=1, bty="n")

# compare db to median
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

#
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
#
#etc.
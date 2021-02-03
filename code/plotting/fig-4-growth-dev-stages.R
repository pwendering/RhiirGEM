# plot predicted growth at the three developmental stages ERM, IRM and ARB
library(wesanderson)

writeToFile = T

# save default graphical parameters
originalPar = par()

# read predicted growth rates from file
inFile <- paste(topDir, "results/developmental-stages/growth-rates-dev-stages.csv", sep = "")
data <- read.table(inFile, header = T, row.names = 1)

# predicted growth rates
growth = data[1,]

# growth rates predicted using +- standard deviation of total protein content
yPos = data[2:nrow(data),]

if (writeToFile) {png(file=paste(topDir,"analysis/figures/growth-dev-stages.png", sep=""),
                      units = "cm", width = 12,height = 10,res = 600, pointsize = 8)}

# define graphical parameters
pal = wes_palette("Cavalcanti1", ncol(growth), type = "discrete")
par(xpd = TRUE,mfrow=c(1,1),mar=c(4,5,2,0)+0.1)
cex.axis = 1.2
cex.lab = 1.5

# create a barplot
xPos <- barplot(as.vector(as.matrix(growth)),
        ylim = c(0,1.4*max(growth)),
        xlim = c(0,4),
        ylab = "",
        angle = 45, density = 50, border = NA,
        space = 0.15,
        col = pal,
        cex.axis = cex.axis,
        cex.lab = cex.lab
)

# add error bars
errorBarWidth=0.1
for (i in 1:ncol(data)) {
  # vertical line
  segments(xPos[i],yPos[1,i],xPos[i],yPos[2,i])
  # horizontal lines
  segments(xPos[i]-errorBarWidth/2,yPos[1,i],xPos[i]+errorBarWidth/2,yPos[1,i])
  segments(xPos[i]-errorBarWidth/2,yPos[2,i],xPos[i]+errorBarWidth/2,yPos[2,i])
}

# add x-axis labels
text(x = xPos, y = -0.1*max(growth), labels = colnames(growth), cex = 1.3, font = 2)

# add y-axis label
mtext(text = expression("predicted growth"*"  ["*h^-1*"]"), side = 2, line = 3, cex = 1.5, font = 3)

if (writeToFile) dev.off()

# reset graphical parameters
par(originalPar)
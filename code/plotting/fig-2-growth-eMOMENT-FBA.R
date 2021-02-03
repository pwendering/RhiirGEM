# plot predicted growth from FBA and eMOMENT against measured hyphal dry weight

library(plotrix)
library(scales)

writeToFile = T

# load general plotting functions and variables
load("code/plotting/rhiir_gem_plotting.rdata")

# save default graphical parameters
originalPar = par()


# read data from FBA and eMOMENT predictions
enzymeFBAFile <- "results/carbon-sources/growth.csv"
eFBA <- as.matrix(read.table(enzymeFBAFile, header = T, row.names = 1))

FBAFile <- "results/carbon-sources/growth-fba.csv"
FBA <- as.matrix(read.table(FBAFile, header = T, row.names = 1))

# define carbon source names
cSourceNames = c("D-Glucose", "D-Fructose", "Raffinose", "Melibiose")

if (writeToFile) {png(file = "results/figures/Figure-2.png",
                      units = "cm", width = 18,height = 12,res = 600, pointsize = 8)
}

# define graphical parameters
par(mgp = c(3,1,0), xpd = F, mar=c(5,5,4,2)+0.1)
pointChar = c(21:24)
cex.legend = 1.3
cex.lab = 1.5
cex.points = 2
cex.axis = 1.2
col.symBorder = "grey40"
colScaleSteps = c(0.05,0.3,1)

# define colors
colorFill = unlist(lapply(rep(c("darkblue","red"),each=nrow(eFBA)),function(x) {
  return(c(alpha(x,colScaleSteps[1]),alpha(x,colScaleSteps[2]),alpha(x,colScaleSteps[3])))}))
colorFill = matrix(colorFill,ncol=ncol(eFBA),byrow = T)

# create scatter plot with gapped y-axis (separating two scales):
eFBA = 10*eFBA
transformYdata <- function(x,from){pmin(x,from) + 0.05*pmax(x-from,0)}
from = 1.5*max(eFBA)
plot(x = rbind(totalMass,totalMass),
     y = transformYdata(rbind(FBA,eFBA),from),
     xlab = expression("hyphae dry weight [mg]"),
     ylab = expression("predicted growth rate ["*h^-1*"]"),
     yaxt = "n",
     cex.lab = cex.lab, cex=cex.points, cex.axis=cex.axis,
     col = col.symBorder, pch = pointChar,
     bg = colorFill)

# add y-axis
yTickPos = c(0,0.001, 0.005,0.01,0.05,0.1,0.15)
yTickLabels = c(0,c(0.001,0.005,0.01)/10,0.05,0.1,0.15)
axis(side=2,at=round(transformYdata(yTickPos,from),4),
     labels = yTickLabels)

box()

# add axis break
break.width = .01
axis.break(2, from,brw = break.width, breakcol = "white", style = "gap")
axis.break(4, from,brw = break.width, breakcol = "white", style = "gap")
axis.break(2, 1.008*from,brw = break.width, breakcol="black", style= "slash")
axis.break(4, 1.008*from,brw = break.width, breakcol="black", style= "slash")

# add a legend
usr = par("usr")
par(usr = c(0,1,0,1))
# add method labels
text(x = 0.03, y = c(.8,.3),
     labels = c("FBA", "eMOMENT"),
     cex = cex.legend, font = 2, pos = 4
)

# add correlations
cor.FBA = round(cor(x = as.vector(FBA),y=as.vector(totalMass),method = "spearman"),2)
text(x = .5,y = .8,labels = bquote(rho[S] ~ "=" ~ .(as.character(cor.FBA))),
     cex = 1.1)
cor.eFBA = round(cor(x = as.vector(eFBA),y=as.vector(totalMass),method = "spearman"),2)
text(x = .5,y = .17,labels = bquote(rho[S] ~ "=" ~ .(as.character(cor.eFBA))),
     cex = 1.1)

# add legends for carbon sources and concentrations
legend(y=.28, x=.68, legend = cSourceNames, pch = pointChar, bty = "n",
       col = "grey40", cex = cex.legend,y.intersp = 1.2, xjust = 0,pt.cex = 2)
legend(y=.28, x=.83, legend = c("0.01 mM","0.1 mM", "1M"),bty = "n",
       fill = c(alpha("black",colScaleSteps[1]),
                alpha("black",colScaleSteps[2]),
                alpha("black",colScaleSteps[3])),
       border = col.symBorder, y.intersp = 1.2,cex = cex.legend, xjust = 0)

par(usr = usr)

if (writeToFile) dev.off()

# reset graphical parameters
par(originalPar)
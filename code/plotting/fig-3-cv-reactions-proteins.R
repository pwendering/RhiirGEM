# plot the coefficients of variation (CV) of protein concentration and reaction fluxes
# across all media conditions

library(scales)

writeToFile = T

# save default graphical parameters
orgininalPar = par()

# load general plotting functions and variables
load("code/plotting/rhiir_gem_plotting.rdata")

# read list of subsystems
subSystemListFile <- "results/subsystems.lst"
subSystemList <- read.table(subSystemListFile, header = F, sep = "\t")
subSystemList = unlist(subSystemList)
uniqueSubsystems = unique(subSystemList)

# convert to factor for individual coloring of points
subsystemFactor = unlist(lapply(subSystemList,function(x)which(uniqueSubsystems==x)))

if (writeToFile) {
  png("results/figures/Figure-3.png",
      units = "cm", width = 15, height = 15, res = 600, pointsize = 7)
}

# define graphical parameters
coords.letter = c(-0.1,.98)
font.letter = 2
cex.letter = 3
cex.axis = 1.5
cex.lab = 2

par(mfrow=c(2,2),
    mgp = c(2,.7,0),
    mar = c(4,4.5,1,1)+0.1,
    xpd = T,
    cex = 0.6,
    cex.lab = 1.5,
    pch = 19,
    mex=1
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CV for proteins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# read data from file
proteinCVFile <- "results/carbon-sources/protein-variances-per-subsystem.csv"
protCV <- read.table(proteinCVFile, header = T, sep = "\t",na.strings = "NaN")

# reorder table based on the order or unique subsystems
newOrder = unlist(lapply(uniqueSubsystems,
                         function(x) which(gsub(x=colnames(protCV),
                                                pattern="_",replacement=" ")==x)))
protCV = protCV[,newOrder]

# remove subsystems with all zero entries
zeroCols = apply(protCV,2,function(x){sum(x,na.rm = T)})==0
protCV = protCV[!zeroCols]

# replace zeros with NA
protCV[protCV<=0] = NA

# convert data to matrix
protCV = as.matrix(protCV)

# logarithmize data 
protCV = log10(protCV)
outlier = protCV<(-5)
print("Number of protein CV outliers:")
print(length(which(outlier)))
protCV[outlier]=NA
# sort columns by median
newOrder = order(apply(protCV,2,function(x) median(x,na.rm = T)),decreasing = T)
protCV=protCV[,newOrder]

# adapt color coding for subsystems
box.colors = color.palette[!zeroCols]
box.colors = box.colors[newOrder]

# create a boxplot
boxplot(x = protCV,
        ylab = expression(log[10]*" CV protein abundance"),
        xaxt = "n",
        outline = FALSE,
        col = box.colors,
        cex.axis = cex.axis,
        cex.lab = cex.lab
)

# add x-axis title
title(xlab="subsystem",line = 1.5,cex.lab = cex.lab)

# add a swarm of points
points(x = rep(seq(1,ncol(protCV)),
               each=nrow(protCV))-0.1*rnorm(nrow(protCV)*ncol(protCV)),
       y = protCV, col = alpha("grey40",0.5))

# draw box
box(lwd = 2)

# add panel letter
usr = par("usr")
par(usr = c(0, 1, 0 ,1))
text(coords.letter[1],coords.letter[2],cex = cex.letter,
     "A", font = font.letter)
par(usr = usr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CV for reactions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# read data from file
rxnCVFile <- "results/carbon-sources/reaction-variances-per-subsystem.csv"
rxnCV <- read.table(rxnCVFile, header = T, sep = "\t",na.strings = "NaN")

# reorder table based on the order or unique subsystems
newOrder = unlist(lapply(uniqueSubsystems,
                         function(x) which(gsub(x=colnames(rxnCV),pattern="_",replacement=" ")==x)))
rxnCV = rxnCV[,newOrder]

# remove subsystems with all zero entries
zeroCols = apply(rxnCV,2,function(x){sum(x,na.rm = T)})==0
rxnCV = rxnCV[,!zeroCols]

# replace zeros with NA
rxnCV[rxnCV<=0] = NA

# convert data to matrix
rxnCV = as.matrix(rxnCV)

# logarithmize data 
rxnCV = log10(rxnCV)
outlier = rxnCV<(-5)
print("Number of reaction CV outliers:")
print(length(which(outlier)))
rxnCV[outlier] = NA
# sort columns by median
newOrder = order(apply(rxnCV,2,function(x) median(x,na.rm = T)),decreasing = T)
rxnCV=rxnCV[,newOrder]

# adapt color coding for subsystems
box.colors = color.palette[!zeroCols]
box.colors = box.colors[newOrder]

# create a boxplot
boxplot(x = rxnCV,
        ylab = expression(log[10]*" CV reaction flux"),
        xaxt = "n",
        outline = FALSE,
        col = box.colors,
        cex.axis = cex.axis,
        cex.lab = cex.lab
)

# add x-axis title
title(xlab="subsystem", line = 1.5, cex.lab = cex.lab)

# add a swarm of points
points(x = rep(seq(1,ncol(rxnCV)),each=nrow(rxnCV))-0.1*rnorm(nrow(rxnCV)*ncol(rxnCV)),
       y=rxnCV, col = alpha("grey40",0.5), pch = 19)
# draw box
box(lwd = 2)

# add panel letter
usr = par("usr")
par(usr = c(0, 1, 0 ,1))
text(coords.letter[1], coords.letter[2], cex = cex.letter,
     "B",font = font.letter)
par(usr = usr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ flux CV vs concentration CV ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# read flux CV data
rxnProtCVFile <- "results/carbon-sources/rxn-prot-match-cv.csv"
rxnProtCVMatch <- read.table(rxnProtCVFile, header = F, sep = "\t",na.strings = "NaN")
rxnProtCVMatch = as.matrix(rxnProtCVMatch)

# read concentration CV data
protCVFile <- "results/carbon-sources/prot-cv.csv"
protCV <- read.table(protCVFile, header = F, sep = "\t",na.strings = "NaN")
protCV = unlist(protCV)

# logarithmize data
protCV=log10(protCV)
outlier = protCV<(-5)
protCV[outlier]=NA

rxnProtCVMatch=log10(rxnProtCVMatch)

# create scatterplot
matplot(protCV,rxnProtCVMatch,pch=19,
      xlab = "",
      ylab = expression(log[10]*" CV reaction flux"),
      col = color.palette[subsystemFactor],
      cex.axis = cex.axis,
      cex.lab = cex.lab
)

# add x axis title
title(xlab = expression(log[10]*" CV protein abundance"), line = 3,
      cex.lab = cex.lab)

# add box
box(lwd = 2)

# add panel letter
usr = par("usr")
par(usr = c(0, 1, 0 ,1))
text(coords.letter[1],coords.letter[2],cex = cex.letter,
     "C",font = font.letter)
par(usr = usr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~ legend for subsystem color coding ~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# create an empty plot
plot(NA,type="n",xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
# add legend
legend(x = -.15, y = 1, legend = uniqueSubsystems, fill = color.palette,
       cex = 2.1, bty = "n")

if (writeToFile) dev.off()

# reset graphical parameters
par(orgininalPar)
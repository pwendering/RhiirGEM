# Generate plots for model stats

library("wesanderson")

writeToFile = T
os = .Platform$OS.type

# adapt input and output files based on the operating system
if (os == "unix") {
  topDir = "/stud/wendering/Documents/Rhizophagus_irregularis_GEM/"
} else if (os == "windows") {
  topDir = "C:\\Users/wende/MobaXterm/home/rhizophagus_irregularis_gem/"
}

# load general plotting functions and variables
load(paste(topDir,"R_scripts/rhiir_gem_plotting.rdata", sep = ""))

inputDir <- paste(topDir, "analysis/stats/", sep = "")
outputDir <- paste(topDir, "analysis/figures/", sep = "")

if (writeToFile) {
  png(paste(outputDir, "model-stats.png", sep = ""), units = "cm",
      width = 20, height = 20, pointsize = 10, res = 300)
}

# save default graphical parameters
originalPar = par();

# figure setting: 2 by 2 panels
par(mfrow=c(2,2))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ read data from file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# read list of subsystems
subSystemListFile <- paste(topDir, "analysis/growth-simulation/subsystems.lst", sep = "")
subSystemList <- read.table(subSystemListFile, header = F, sep = "\t")
subSystemList = unlist(subSystemList)
uniqueSubsystems = unique(subSystemList)

# create a list with all model statistics
stat = list()
stat$general <- read_stat_file(paste(inputDir, "model_stats.txt", sep = ""))
stat$thermodynamics <- read_stat_file(paste(inputDir, "thermodynamics.txt", sep = ""))
stat$balancing <- read_stat_file(paste(inputDir, "balancing.txt", sep = ""))
stat$generic_reactions <- read_stat_file(paste(inputDir, "generic_reactions.txt", sep = ""))
stat$co_lco <- read_stat_file(paste(inputDir, "co_lco.txt", sep = ""))
stat$lipid <- read_stat_file(paste(inputDir, "lipid_stats.txt", sep = ""))
stat$maps = read_stat_file(paste(inputDir, "kegg_map_rxns.txt", sep = ""))
stat$subsystems <- read_stat_file(paste(inputDir, "subsystem_rxns.txt", sep = ""))
stat$subsystems = stat$subsystems[unlist(lapply(uniqueSubsystems,
  function(x) which(gsub(names(stat$subsystems),pattern="\\.",replacement=" ")==x)))]
stat$kcats <- read_stat_file(paste(inputDir, "kcats.txt", sep = ""))
stat$metBRITE <- read_stat_file(paste(inputDir, "met-brite.txt", sep = ""))
stat$metBRITEOther <- read_stat_file(paste(inputDir, "met-brite-other.txt", sep = ""))
stat$metBRITELipids <- read_stat_file(paste(inputDir, "met-brite-lipids.txt", sep = ""))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ subsystems ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# combine subsystems with low numbers of reactions into "Other"
stat$subsystems$Other = 
  stat$subsystems$Other +
  stat$subsystems$Inositol.phosphate.metabolism + 
  stat$subsystems$Terpenoid.backbone.biosynthesis +
  stat$subsystems$Sink +
  stat$subsystems$Medium +
  stat$subsystems$Amino.sugar.and.nucleotide.sugar.metabolism
keepIdx = which(
  !(names(stat$subsystems) %in% c("Inositol.phosphate.metabolism",
                                  "Terpenoid.backbone.biosynthesis",
                                  "Sink",
                                  "Medium",
                                  "Amino.sugar.and.nucleotide.sugar.metabolism")))
stat$subsystems = stat$subsystems[keepIdx]

# adapt colors
pie.col = color.palette[keepIdx]


# order subsystems and colors based on the number of associated reactions
idx = order(stat$subsystems, decreasing = T)
pie.col = pie.col[idx]
stat$subsystems = stat$subsystems[idx]

# replace dots in subsystem names with spaces
annotation = gsub("\\.", " ", names(stat$subsystems))

# create pie plot
par(mar = c(3,4,10,4)+0.1)
pie(as.numeric(stat$subsystems),
    labels = NA,
    col = pie.col,
    cex=1.2,
    init.angle=90,
    lwd = 1.5,
    clockwise = T
)

# add a legend
usr = par("usr")
par(usr = c(0,1,0,1), xpd = T)
legend(x=0,y=1.5, legend = annotation, fill = pie.col, cex = 1.1,
       bty = "n")
par(usr = usr)


# ~~~~~~~~~~~~~~~~~~~~~ BRITE classification of metabolites ~~~~~~~~~~~~~~~~~~~~~~ #
# combine classes with low numbers of metabolites into "Other"
stat$metBRITE$Other = 
  stat$metBRITE$Other +
  stat$metBRITE$Pseudometabolites +
  stat$metBRITE$Terpenoids + 
  stat$metBRITE$Lipochitooligosaccharides + 
  stat$metBRITE$Minerals +
  stat$metBRITE$Glycans
keepIdx = which(
  !(names(stat$metBRITE) %in% c("Pseudometabolites",
                                "Terpenoids",
                                "Lipochitooligosaccharides",
                                "Minerals",
                                "Glycans")))
stat$metBRITE = stat$metBRITE[keepIdx]

# order classes based on the number of associated metabolites
idx = order(stat$metBRITE, decreasing = T)
stat$metBRITE = stat$metBRITE[idx]

# replace dots with spaces in BRITE classifications
annotation = gsub("\\.", " ", names(stat$metBRITE))

# define colors
pie.col = wes_palette("Royal1", length(stat$metBRITE),
                      type = "continuous")

# create a pie plot
par(mar = c(3,4,10,4)+0.1, xpd = T)
pie(as.numeric(stat$metBRITE),
    labels = NA,
    col = pie.col,
    cex=1.2,
    init.angle=0,
    clockwise = T
)

# add a legend
usr = par("usr")
par(usr = c(0,1,0,1), xpd = T)
legend(x=0,y=1.5, legend = annotation, fill = pie.col, cex = 1.1,
       bty = "n")
par(usr = usr)

# ~~~~~~~~~~~~~~~~~~~~~~~~ general reaction characteristics ~~~~~~~~~~~~~~~~~~~~~~~~~ #
# number of reactions in the model
n_rxns = stat$general$reactions

# create data matrix
rxn_general <- matrix(c(stat$general$reactions_with_EC, n_rxns-stat$general$reactions_with_EC,
                        stat$general$transport, n_rxns-stat$general$transport,
                        stat$general$gpr, n_rxns-stat$general$gpr,
                        n_rxns-stat$general$mass_imbalanced, stat$general$mass_imbalanced,
                        n_rxns-stat$general$charge_imbalanced, stat$general$charge_imbalanced,
                        n_rxns-stat$general$blocked, stat$general$blocked,
                        1035, 537),
                      nrow = 2)

# define row annotation for barplot
annotation = c("E.C. assigned", "transport", "GPR rule", "mass-balanced", "charge-balanced",
               "unblocked", expression(paste(Delta, G, degree, " assigned")),
               expression("k"[cat]*" assigned"))

# define colors for barplot
bar.col = wes_palette("Royal1", ncol(rxn_general), type = "continuous")

# create barplot
par(mar = c(5, 10, 2, 2) + 0.1)
rev_order = ncol(rxn_general):1
x_lim = c(0,n_rxns)
x <- barplot(rxn_general[,rev_order],
             beside = F,
             horiz = T,
             xaxt = "n",
             col = bar.col,
             xlim = x_lim,
             space = .1
)

# add row annotation
text(x = -55, y = x, labels = annotation[rev_order],
     xpd = T, srt = 0, cex = 1.3, adj = 1
)
# add numbers for each statistic
text(x = rxn_general[1,rev_order], y = x, labels = as.character(rxn_general[1,rev_order]),
     cex = 1.1, font = 4, pos = 2)
# add sum of reactions
text(x=n_rxns, y=x[length(x)], bquote(Sigma==.(n_rxns)), cex = 1.3, pos = 2)

if (writeToFile) dev.off()

# reset graphical parameters
par(originalPar)
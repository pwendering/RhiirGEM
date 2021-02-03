# comparison to other fungal models
library(pheatmap)
library(wesanderson)
library(ape)

writeToFile = T

filename <- paste(topDir, "results/fungal-models/jd_ec_class.txt", sep = "")

# NJ / UPGMA phylogeny run on MAFFT server (MAFFT 7; 31.8.2020)
treeFile <- paste(topDir, "results/fungal-models/18S-rRNA/18S_archaeopteryx_js.tre", sep = "")

output_dir = paste(topDir, "results/figures/", sep = "")

data <- read.table(filename, header = TRUE)
rownames(data) <- data[,1]
data = data[,-1]

# update rownames
rownames(data) <- gsub("_", ". ", rownames(data))

# generate a tree from the NJ tree 
tree <- read.tree(file = treeFile)

toDrop <- c("Eremothecium_gossypii", 
            "Lachancea_kluyveri")
idxDrop = unlist(lapply(toDrop,function(x){grep(x,tree$tip.label)}))
tree = drop.tip(tree, tree$tip.label[idxDrop])

# re-write the tree tips to species names as formatted in data
tree$tip.label = gsub("^[0-9]+_", "", tree$tip.label)
genusLetter <- substr(tree$tip.label[-grep("Rhizophagus_irregularis",tree$tip.label)], 1,1)
speciesNames = gsub("^[A-Za-z]+_", "", tree$tip.label[-grep("Rhizophagus_irregularis",tree$tip.label)])
clustLabels <- paste(genusLetter, ". ", speciesNames, sep = "")
clustLabels = gsub("_.*$","",clustLabels)
clustLabels <- clustLabels[length(clustLabels):1]

# find the indices of the rows in data corresponding to the guide tree labels
idxMatch <- match(clustLabels, rownames(data))
data = data[idxMatch,]

# remove rows and columns with only ones
data = data[,-which(apply(data, 2, sum)==nrow(data))]

# append row of ones for R. irregularis
data = rbind(data, rep(0,ncol(data)))
rownames(data)[nrow(data)] = "R. irregularis"
head(data)

# format species names for plotting
newnames <- as.expression(lapply(
  rownames(data),
  function(x) bquote(italic(.(x)))))

# color palette
pal_length = 10
pal = colorRampPalette(colors = c("light blue", "firebrick4"))(pal_length)


if (writeToFile) pdf(paste(output_dir, "comparison_to_fungal_models.pdf", sep = ""),
                     pointsize = 14)

par(mar=c(5,6,4,4)+0.1)

# --- plot the phylogenetic tree ---
plot(tree, type = "phylogram", use.edge.length = F,show.tip.label = F,edge.width = 2)

# --- plot the heatmap ---

# define breaks
br = c(seq(from = 0,to = 1, length.out =  pal_length))
pheatmap(as.matrix(data),
         cluster_rows = F,
         cluster_cols = F,
         color = pal,
         breaks = br,
         labels_row = as.expression(newnames),
         fontsize = 18,cellheight = 25,cellwidth=50
)

if(writeToFile) dev.off()



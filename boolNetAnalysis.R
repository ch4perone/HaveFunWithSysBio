
install.packages("BoolNet")
library("BoolNet")

setwd("/home/chaperone/projects/computationalSystemBiology")

network = loadSBML("MODEL1606020000.xml")

microenvironment = c("IL12e", "IFNGe", "IL2e", "IL4e", "IL21e", "TGFBe", "IL10e", "IL27e", "INSULIN")
pro_Th0 = c(0,0,0,0,0,0,0,0,0)
pro_Th1 = c(1,1,0,0,0,0,0,0,0)
pro_Th2 = c(0,0,1,1,0,0,0,0,0)
pro_Th17 =c(0,0,0,0,1,1,0,0,0)
pro_iTreg=c(0,0,1,0,0,1,0,0,0)
pro_Tr1 = c(0,0,0,0,0,0,1,1,0)

?perturbNetwork()
pert



network = fixGenes(network, fixIndices = microenvironment, values = pro_iTreg)
#fixGenes(network, "INSULIN", 0)

STG = getAttractors(network)

plotStateGraph(STG)


#
# Perturb
#

library(devtools)
install_github("mar-esther23/boolnet-perturb")
library(BoolNetPerturb)

#rules = data.frame()

#labelsTh17Treg = load("./labelsTh17Treg.rda")


labels = c()
rules = c()

labels = c("Th0", "Th1", "TBET+", "Th1R", "TH2", "GATA3+" , "Th2R", "Th17" , "RORGT+" , "iTreg", "IL10+", "TGFB+")
rules = c("!(TBET | GATA3 | RORGT | FOXP3 | IL10 | TGFB)", "(TBET & IFNG) & !(IL10 | TGFB | FOXP3)", " TBET & !(IFNG | IL10 | TGFB | FOXP3) ", " TBET & (IL10 | TGFB | FOXP3) ", "(GATA3 & IL4) & ! (IL10 | TGFB | FOXP3) ",
         " GATA3 & ! (IL4 | IL10 | TGFB| FOXP3)" , " GATA3 & (IL10 | TGFB | FOXP3)" , "RORGT & IL21 & ! IL10" ,
         "RORGT & ! (IL21 | IL10)" , "FOXP3 & TGFB & ! (TBET | GATA3 | RORGT)" , "IL10 & ! (TBET | GATA3 | FOXP3 | RORGT)",
         "TGFB & ! (TBET | GATA3 | FOXP3 | RORGT)")

df.rules = data.frame(labdels, rules)



BoolNetPerturb::labelAttractors(STG, df.rules)

#
# Label attractors of STG
#

#
# Get cell fate map
# >> computes for all attractors, every transient perturbation for explicit "genes"
# >> lists from initial state, which cell fate is attained after each transient flip
# >> make plot out of this


map = BoolNetPerturb::cellFateMap(network, genes = c("TBET", "IFNG", "GATA3", "IL2", "IL4", "RORGT", "IL21", "FOXP3", "TGFB", "IL10"), label.rules = df.rules)

#
# Compute adjancency matrix & plot states
#

plot_labels = c(labels, "IL10+TGFB+", "IL10+TGFB+/iTreg", "Th1R/Th1R", "RORGT+/TGFB+")
Adj = matrix(0, nrow = length(plot_labels), ncol = length(plot_labels))
colnames(Adj) = plot_labels
rownames(Adj) = plot_labels
             


for(cellType in plot_labels) {
  finalStates = map$final[which(map$initial == cellType)]  
  for (state in finalStates) {
    Adj[cellType, state] = Adj[cellType, state] + 1
  }
}

#
# Plot using igraph
#

library(igraph)
#import the sample_dw_adj.csv file:
#dat=read.csv(file.choose(),header=TRUE,row.names=1,check.names=FALSE) # read .csv file
#m=as.matrix(dat)
net=graph.adjacency(Adj,mode="directed",weighted=TRUE,diag=FALSE) #the only difference between this and the weighted network code is that mode="directed"
#net[from=V(net), to = V(net)] = 1
for(cellType in plot_labels) {
  if (Adj[cellType, cellType] > 0) {
    net[from = cellType, to = cellType ] = Adj[cellType, cellType]  
  }
}
layout.self = matrix(ncol = 2, nrow = length(plot_labels))
rownames(layout.self) = plot_labels
colnames(layout.self) = c("x", "y")

plot.igraph(net, vertex.size = 20, vertex.label=V(net)$name, layout=layout.circle, vertex.label.color="black",edge.color="black",edge.width=E(net)$weight/2, edge.arrow.size=1.5)



plot.igraph(net, vertex.size = 20, vertex.label=V(net)$name, layout=layout.circle,vertex.label.color="black", edge.color="black")


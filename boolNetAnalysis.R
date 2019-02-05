library(reshape2)
library(ggplot2)
library(scales)

#install.packages("BoolNet")
library("BoolNet")
#library(devtools)
#install_github("mar-esther23/boolnet-perturb")
library(BoolNetPerturb)

#
# Setup and Definition
#

setwd("/home/aliki/Documents/rechnergbioinformatik")
network = loadSBML("MODEL1606020000.xml")

microenvironment = c("IL12e", "IFNGe", "IL2e", "IL4e", "IL21e", "TGFBe", "IL10e", "IL27e", "INSULIN")
pro_Th0 = c(0,0,0,0,0,0,0,0,0)
pro_Th1 = c(1,1,0,0,0,0,0,0,0)
pro_Th2 = c(0,0,1,1,0,0,0,0,0)
pro_Th17 =c(0,0,0,0,1,1,0,0,0)
pro_iTreg=c(0,0,1,0,0,1,0,0,0)
pro_Tr1 = c(0,0,0,0,0,0,1,1,0)

add_pro_Th9 = c(0,0,0,1,0,1,0,0,0)
add_pro_Tfh = c(1,0,0,0,1,0,0,0,0) # maybe Il2e instead of IL12e

MicroEnv = matrix(c(pro_Th0,pro_Th1,pro_Th2,pro_Th17,pro_iTreg,pro_Tr1), ncol = length(microenvironment), byrow = TRUE)
colnames(MicroEnv) = microenvironment
rownames(MicroEnv) = c("pro_Th0", "pro_Th1", "pro_Th2", "pro_Th17", "pro_iTreg", "pro_Tr1")


network = fixGenes(network, fixIndices = microenvironment, values = add_pro_Th9)

STG = getAttractors(network)
plotStateGraph(STG)

#
# Attractor Rules for Labeling
#

labels = c("Th0", "Th1", "TBET+", "Th1R", "Th2", "GATA3+" , "Th2R", "Th17" , "RORGT+" , "iTreg", "IL10+", "TGFB+")
rules = c("!(TBET | GATA3 | RORGT | FOXP3 | IL10 | TGFB)", "(TBET & IFNG) & !(IL10 | TGFB | FOXP3)", " TBET & !(IFNG | IL10 | TGFB | FOXP3) ", " TBET & (IL10 | TGFB | FOXP3) ", "(GATA3 & IL4) & ! (IL10 | TGFB | FOXP3) ",
         " GATA3 & ! (IL4 | IL10 | TGFB| FOXP3)" , " GATA3 & (IL10 | TGFB | FOXP3)" , "RORGT & IL21 & ! IL10" ,
         "RORGT & ! (IL21 | IL10)" , "FOXP3 & TGFB & ! (TBET | GATA3 | RORGT)" , "IL10 & ! (TBET | GATA3 | FOXP3 | RORGT)",
         "TGFB & ! (TBET | GATA3 | FOXP3 | RORGT)")

hybrid_label = "Tbet+Gata3+Foxp3+"
hybrid_rule = "TBET & GATA3 & FOXP3"

df.rules = data.frame(labels, rules)

#
# Label attractors of STG
#


BoolNetPerturb::labelAttractors(STG, df.rules)

#
# Determin Cell Differentiation 
# > based upon microenvironment (without perturbations)
# > display basin of attraction size for different celltypes

getCellDifferentiationBasinSizes = function(network, micro_env, micro_val, insulin = 0, knockout=vector(), knockin=vector(), attractorLabels, label.rules) {
  basinSizes = rep(0, length(attractorLabels)) 
  names(basinSizes) = attractorLabels
  
  #fix microenvironment | knock-out genes | and insulin
  micro_val[length(micro_val)] = insulin

  for(gene in knockin) {
    micro_env = c(micro_env, gene)
    micro_val = c(micro_val, 1)
  }
  for(gene in knockout) {
    micro_env = c(micro_env, gene)
    micro_val = c(micro_val, 0)
  }
  network = fixGenes(network, fixIndices = micro_env, values = micro_val)
  attr = getAttractors(network, type = "synchronous")
  attrLabels = BoolNetPerturb::labelAttractors(attr, label.rules = label.rules)
  
  for (i in 1:length(attr$attractors)) {
    if (attrLabels[i] %in% attractorLabels) {
      basinSizes[attrLabels[i]] = basinSizes[attrLabels[i]] + attr$attractors[[i]]$basinSize
    } else {
      print(attrLabels[i])
    }
  }
  
  return(basinSizes)  
}

#getCellDifferentiationBasinSizes(network, micro_env = microenvironment, micro_val = pro_Th0, knockout = "IL10", insulin = 0, attractorLabels = c(labels, "IL10+TGFB+"), label.rules = df.rules)

# Compute entire matrix
M = vector()
for(i in 1:nrow(MicroEnv)) {
  pro_mic = MicroEnv[i,]
  M = c(M, getCellDifferentiationBasinSizes(network, micro_env = microenvironment, micro_val = pro_mic, insulin = 0, attractorLabels = c(labels, "IL10+TGFB+"), label.rules = df.rules))
}
M = matrix(M, nrow = nrow(MicroEnv), byrow = TRUE)
colnames(M) = c(labels, "IL10+TGFB+")
rownames(M) = rownames(MicroEnv)




M = M[,c(-3,-9)]
y_ord = rownames(M)
jpeg('basalHeatmap.jpg',units = "cm" , width = 25.0, height = 15.0, res = 600)
ggplot(melt(M), aes(Var2, ordered(Var1, rev(y_ord)) , fill=value)) +
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient(na.value = "white", low = "thistle1", high = "deeppink3", trans = log10_trans(), limits = c(NA, 1024)) +
  xlab("Cell type") +
  ylab("Micro-environment")
dev.off()

# Compute entire INS matrix
M_INS = vector()
for(i in 1:nrow(MicroEnv)) {
  pro_mic = MicroEnv[i,]
  M_INS = c(M_INS, getCellDifferentiationBasinSizes(network, micro_env = microenvironment, micro_val = pro_mic, insulin = 1, attractorLabels = c(labels, "IL10+TGFB+"), label.rules = df.rules))
}
M_INS = matrix(M_INS, nrow = nrow(MicroEnv), byrow = TRUE)
colnames(M_INS) = c(labels, "IL10+TGFB+")
rownames(M_INS) = rownames(MicroEnv)


M_INS = M_INS[,c(-3,-9)]
y_ord = rownames(M_INS)
jpeg('hyperHeatmap.jpg',units = "cm" , width = 25.0, height = 15.0, res = 600)
ggplot(melt(M_INS), aes(Var2, ordered(Var1, rev(y_ord)) , fill=value)) +
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient(na.value = "white", low = "thistle1", high = "deeppink3", trans = log10_trans(), limits = c(NA, 1024)) +
  xlab("Cell type") +
  ylab("Micro-environment")
dev.off()

# statistics

effectorSum = sum(M[, c("Th1", "Th2", "Th17")])
regulatorySum = sum(M[, c("Th1R", "Th2R", "iTreg", "IL10+" , "IL10+TGFB+")])

effectorSum / regulatorySum

effectorSumINS = sum(M_INS[, c("Th1", "Th2", "Th17")])
regulatorySumINS = sum(M_INS[, c("Th1R", "Th2R", "iTreg", "IL10+", "IL10+TGFB+")])

effectorSumINS / regulatorySumINS


# to compensate for Insulin: Use GATA3 (gene 3) => 0.6 eff/reg  (basically deletes Th2 (effector) cell types)


# choose colors from http://sape.inf.usi.ch/quick-reference/ggplot2/colour

# Question: IL10+TGFB+ how is this determined (shouldn't this be the intersection between IL10 & TGFB+ which it is clearly not)

#
# Testing asynchronous update
#   for specific initial states
#

compareUpdateMethod = function(network, init_state, micro_env, micro_val, label.rules) {
  fixedNet = fixGenes(network, fixIndices = micro_env, values = micro_val)
  attr_sync = getAttractors(fixedNet, type = "synchronous", startStates = list(c(init_state, micro_val)))
  attr_async = getAttractors(fixedNet, type = "asynchronous", startStates = list(c(init_state, micro_val)))
  
  sync = BoolNetPerturb::labelAttractors(attr_sync, label.rules = label.rules)
  async = BoolNetPerturb::labelAttractors(attr_async, label.rules = label.rules)
  
  return(data.frame(sync, async))
}

compareUpdateMethod(network, init_state = rep(0, 10), micro_env = microenvironment, micro_val = pro_iTreg, label.rules = df.rules)


for(i in 1:nrow(MicroEnv)) {
  pro_mic = MicroEnv[i,]
  cat('\n',rownames(MicroEnv)[i], ":\n")
  result = compareUpdateMethod(network, init_state = rep(1, 10), micro_env = microenvironment, micro_val = pro_mic, label.rules = df.rules)
  print(result$sync)
  print(result$async)
  cat(as.character(result$sync) == as.character(result$async), "\n")
}

# For all possible initial states
all_states = expand.grid(rep(list(0:1), 10))
Y = vector()
GROUP = vector()
globalAgreementVector = vector()

for (k in 1:5) {
  counter = 0
  counter_vector = rep(0, 11) # errors in vector according number of 1's, starting with 0 1's
  full_vector = rep(0, 11) 
  for (j in 1:nrow(all_states) ) { # might take awhile
    init = all_states[j,]
    num_ones = sum(init)
    
    for(i in 1:nrow(MicroEnv)) {
      pro_mic = MicroEnv[i,]
      #cat('\n',rownames(MicroEnv)[i], ":\n")
      result = compareUpdateMethod(network, init_state = init, micro_env = microenvironment, micro_val = pro_mic, label.rules = df.rules)
      #cat(unlist(result), "\n")
      
      
      full_vector[num_ones + 1] = full_vector[num_ones + 1] + 1
      if(as.character(result$sync) == as.character(result$async)) {
        counter = counter + 1
        counter_vector[num_ones + 1] = counter_vector[num_ones + 1] + 1
      }
      
      if(length(as.character(result$async)) > 1) {
        print(2)
      }
      
    }
  }

  cat("iteration: ", k, " done\n")  
  cat("total: ", counter, " percent", counter / (1024 * nrow(MicroEnv)) * 100, "\n")
  full_vector
  counter_vector
  globalAgreementVector = c(globalAgreementVector, counter / (1024 * nrow(MicroEnv)))
  
  Y = c(Y, counter_vector / full_vector)
  GROUP = c(GROUP, seq(0, 10))
} 


df = data.frame(Y, GROUP)

ggplot(melt(globalAgreementVector), aes(y=value, "all states")) +
  ylab("agreement of update methods") +
  geom_boxplot(fill = "lightblue")

ggplot(df, aes(x=GROUP, y=Y)) +
  geom_boxplot(aes(x = GROUP),  fill = "lightblue") +
  ylab("agreement of update methods") +
  xlab("initial number of genes active") 

#
# Construct Cell Fate Map
# >> computes for all attractors, every transient perturbation for explicit "genes"
# >> lists from initial state, which cell fate is attained after each transient flip
# >> make plot out of this


createCellFateMap = function(network, micro_env, micro_val, insulin = 0, df.rules, plot_labels) {
  micro_val[length(micro_val)] = insulin
  network = fixGenes(network, fixIndices = micro_env, values = micro_val)
  map = BoolNetPerturb::cellFateMap(network, genes = c("TBET", "IFNG", "GATA3", "IL2", "IL4", "RORGT", "IL21", "FOXP3", "TGFB", "IL10"), label.rules = df.rules)
  
  # Compute adjancency matrix & plot states
  Adj = matrix(0, nrow = length(plot_labels), ncol = length(plot_labels))
  colnames(Adj) = plot_labels
  rownames(Adj) = plot_labels
  
  for(cellType in plot_labels) {
    finalStates = map$final[which(map$initial == cellType)]  
    for (state in finalStates) {
      if (state %in% plot_labels) {
        Adj[cellType, state] = Adj[cellType, state] + 1  
      }
    }
  }
  
  # Plot using igraph
  library(igraph)
  #Adj["Th0", "Th1"] = Adj["Th0", "Th1"] + 1
  net=graph.adjacency(Adj,mode="directed",weighted=TRUE,diag=FALSE) #the only difference between this and the weighted network code is that mode="directed"
  
  for(cellType in plot_labels) {
    if (Adj[cellType, cellType] > 0) {
      print(cellType)
      net[from = cellType, to = cellType ] = Adj[cellType, cellType]  
    }
  }
  #net[from = "Th0", to = "Th1" ] = Adj["Th0", "Th1"] - 1
  layout.self = matrix(-2, ncol = 2, nrow = length(plot_labels))
  rownames(layout.self) = plot_labels
  colnames(layout.self) = c("x", "y")
  
  layout.self['Th0', ] = c(0,0)
  layout.self['Th1', ] = c(0.3,1) *3
  layout.self['Th1R', ] =c(0.7,0.5) *3
  layout.self['Th17', ] =c(1,0) *3
  layout.self['iTreg', ]=c(0.7,-0.5) *3
  #layout.self['RORGT+', ]=c(0.3,-1)
  layout.self["IL10+TGFB+", ] = c(-0.3,-1) *3
  layout.self['IL10+',] = c(-0.7, -0.5)  *3
  layout.self['Th2R',] = c(-1, 0) *3
  layout.self['Th2',] = c(-0.7, 0.5) *3
  #layout.self['TGFB+',] = c(-0.3, 1)

  layout.self['TBET+',] = c(0.3, -1) *3
  layout.self['GATA3+',] = c(-0.3, 1) *3
  
  print(layout.self)
  regCellTypes = c("iTreg", "Th1R", "Th2R", "Tr1", "IL10+", "TGFB+", "IL10+TGFB+")
  effCellTypes = vector() #c("Th1", "Th2", "Th17")
  colors = c("orange", "skyblue", "green", "red", "grey")[1 + V(net)$name %in% c("") + 2 * V(net)$name %in% regCellTypes + 3 * V(net)$name %in% effCellTypes  + 4 * V(net)$name %in% c("GATA3+", "TBET+", "Th0")]
  print(Adj)
  plt = plot.igraph(net, edge.arrow.size=0.75, vertex.color = colors, vertex.size = 32, vertex.label=V(net)$name, layout=layout.self,  vertex.label.color="black",edge.color="black",edge.width=E(net)$weight/3, edge.arrow.size=1.5)
  #plot.igraph(net, vertex.size = 20, vertex.label=V(net)$name, layout=layout.circle,vertex.label.color="black", edge.color="black")
  return(Adj)
}
plot.labels = c(labels[c(-12,-9)], "IL10+TGFB+")
adj = createCellFateMap(network, micro_env = microenvironment, micro_val = pro_Tr1, insulin = 0,  df.rules = df.rules, plot_labels = plot.labels)

#
# Why 12 states (13), but randomly only 10 (11) are used? TGFB+ & RORGT+ dropped? Why?
# Example iTreg. Why oversimplified .. compared to real network
# 
# What is Tr1? It is not defined in the table
# Why Tr1 showing up in the heatmap (11), not so in the cell map fate (alos 11 states)
# Paper is highly inconsistent
#
#
# Clear Definitions of regulatory (iTreg, Th1R, Th2R, Tr1??) vs effector cells (Th1, Th2, Th17)? What about the rest??
#

#
# Find ideal microenvironments
#
# Computes all permutations of mictoenvironments, fixes the network & computes basin of attraction for each target label (cell type)

findIdealMicroenvironment = function(network, insulin = 0, label.rules, target.label) {
  
}





#"IL10+TGFB+/iTreg", "Th1R/Th1R", "RORGT+/TGFB+"
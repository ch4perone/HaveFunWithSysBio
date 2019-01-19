
install.packages("BoolNet")
library("BoolNet")

setwd("/home/chaperone/projects/computationalSystemBiology")

network = loadSBML("MODEL1606020000.xml")

microenvironment = c("IL12e", "IFNGe", "IL2e", "IL4e", "IL21e", "TGFBe", "IL10e", "IL27e", "INSULIN")
pro_Th0 = c(0,0,0,0,0,0,0,0,0)
pro_Th1 = c(1,1,0,0,0,0,0,0,0)




network = fixGenes(network, fixIndices = microenvironment, values = pro_Th0)
#fixGenes(network, "INSULIN", 0)

STG = getAttractors(network)

plotStateGraph(STG)


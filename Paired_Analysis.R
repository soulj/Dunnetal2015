#script for analysis of paired OA damaged vs undamaged data

#Define the working directory
setwd("~/Testing123/PairedOAAnalysis/")

#Do you want to visualise the networks in cytoscape?
visualise = FALSE


#load the required libraries
library("DESeq2")
library("gdata")
library("data.table")
library("lattice")
library("RColorBrewer")
library("gridExtra")
library("RUVSeq")
library("EDASeq")
library("gplots")
library("lumi")
library("lumiHumanAll.db")
library("annotate")
library("limma")
library("Matrix")
library("igraph")
library("data.table")
library("RCytoscape") # also librarys cytoscape v2.8 to be open with the Cytoscape RPC plugin active
library("annotate")
library("org.Hs.eg.db")
library("VennDiagram")
library("reshape2")
library("ggplot2")
library("GO.db")
require("goseq")
library("reactome.db")
library("KEGG.db")


#load the pre-processed data to save time
load("./Ensembl2Genes.RData")
load("./PreprocessedCountsMatrix.RData")


#################################
#Differential Expression Analysis
#################################

#caculate fold change and pvalues with DESeq2
colData=data.frame(Patient=as.factor(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8)),Tissue=as.factor(c("PLC","DMC","PLC","DMC","PLC","DMC","PLC","DMC","PLC","DMC","PLC","DMC","PLC","DMC","PLC","DMC")))
dds<-DESeqDataSetFromMatrix(countData= countsMatrix.complete, colData= colData,design=~Patient + Tissue)

#caclulate the fold changes
dds<-DESeq(dds)

#control for the non-uniform variance in the read counts.
vsd <- varianceStabilizingTransformation(dds)
library("vsn")
colours=c(brewer.pal(12, "Paired"),"lightgrey", "black","#996633", "#663300")
batch=c(18,17,17,16,17,16,17,17,11,11,15,15,15,15,15,15)

ntop=500
intgroup =c("Patient","Tissue")
rv = apply(assay(vsd), 1, var)
select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca = prcomp(t(assay(vsd)[select, ]))
fac = factor(apply(as.data.frame(colData(vsd)[, intgroup, drop = FALSE]), 
                   1, paste, collapse = " : "))
xyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$x), cex = 2, aspect = "iso", col = colours, pch=batch,main = draw.key(key = list(rect = list(col = colours), text = list(levels(fac)), rep = FALSE)))


# Strong Batch effect - will correct with RUVSeq
countsMatrix2<-  sapply(countsMatrix.complete,as.numeric)
rownames(countsMatrix2)=rownames(countsMatrix.complete)
filter<-  apply(countsMatrix.complete, 1 , function(x)length(x[x>10]) >=8)
filtered<-countsMatrix.complete[filter,]
filtered=as.matrix(filtered)
dds<-DESeqDataSetFromMatrix(countData=filtered, colData= colData,design=~Patient + Tissue)


#caclulate the fold changes
dds<-DESeq(dds)
res<-as.data.frame(results(dds,contrast=c("Tissue","DMC","PLC")))
geneRank=res[order(abs(res$pvalue)),]
#take all but the top 5000 genes ranked by differential expression
empirical<-  rownames(geneRank)[  which  ( !(rownames(geneRank)%in%rownames (geneRank)[1:5000]))]

#Set the number of surrogate variables to be 2 - based on the resulting PCA
set=RUVg(filtered,empirical,k=2)

#Look at the distribution of the batch effect corrected samples

dds<-DESeqDataSetFromMatrix(countData=set$normalizedCounts, colData= colData,design=~Patient + Tissue)

#caclulate the fold changes
dds<-DESeq(dds)
vsd <- DESeq2::varianceStabilizingTransformation(dds)
par(mfrow=c(1,3))
notAllZero <- (rowSums(counts(dds))>0)
ntop=500
intgroup =c("Patient","Tissue")
rv = apply(assay(vsd), 1, var)
select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca = prcomp(t(assay(vsd)[select, ]))
fac = factor(apply(as.data.frame(colData(vsd)[, intgroup, drop = FALSE]), 
                   1, paste, collapse = " : "))
xyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$x), cex = 2, aspect = "iso", col = colours, pch=batch,main = draw.key(key = list(rect = list(col = colours), text = list(levels(fac)), rep = FALSE)))

distsRL<- dist(t    ( assay(vsd)))
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(Patient, Tissue, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))


# Use the calulated covariantes from RUVSeq to identify the differentially expressed genes
colData=data.frame(Patient=as.factor(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8)),W_1=set$W[,1],W_2=set$W[,2],Tissue=c("PLC","DMC","PLC","DMC","PLC","DMC","PLC","DMC","PLC","DMC","PLC","DMC","PLC","DMC","PLC","DMC"))
dds<-DESeqDataSetFromMatrix(countData=filtered, colData= colData,design=~Patient + W_1+ W_2 + Tissue)

#caclulate the fold changes
dds<-DESeq(dds)
res<-as.data.frame(results(dds,contrast=c("Tissue","DMC","PLC")))
geneRank=res[order(-abs(res$log2FoldChange)),]
geneRank=merge(geneRank,Genes,by.x="row.names",by.y="ID")
geneRank=geneRank[order(-abs(geneRank$log2FoldChange)),]
normalisedCounts=counts(dds, normalized =TRUE)
geneRank=merge(geneRank,normalisedCounts,by.y="row.names",by.x="Row.names")
geneRank=geneRank[order(-abs(geneRank$log2FoldChange)),]
geneRank$FoldChange<-2^geneRank$log2FoldChange
geneRank$FoldChange<-ifelse(geneRank$FoldChange < 1, -1/geneRank$FoldChange,geneRank$FoldChange)
write.table(geneRank,file="./Results/DESeq2Results.tab",row.names=F,col.names=T,sep="\t",quote=F)


#Just the differentially expressed genes

DEGs=geneRank[abs(geneRank$log2FoldChange)>=0.58,]
DEGs=DEGs[DEGs$padj<=0.1,]
DEGs<-na.omit(DEGs)
DEGs<-DEGs[,c(1,8,25,6,7)]
DEGs<-DEGs[ !duplicated(DEGs$Gene.name),]
colnames(DEGs)[1]<-"Ensembl Gene ID"
write.table(DEGs,file="./Results/DESeq2ResultsDEGs.tab",row.names=F,col.names=T,sep="\t",quote=F)


#split the ID
transcriptIDs<- do.call('rbind', strsplit(as.character(geneRank$Row.names),'.',fixed=TRUE))

########################################
#Pathway enrichment analysis with GOSeq
#########################################

genes=ifelse(geneRank$Gene.name %in% DEGs$Gene.name,1,0)
names(genes)=transcriptIDs[,1]


pwf=nullp(genes,"hg19","ensGene")

# find KEGG pathways

KEGG=goseq(pwf,'hg19','ensGene',test.cats="KEGG")
KEGG$padj=p.adjust(KEGG$over_represented_pvalue,method="BH")
KEGG.sig=KEGG[KEGG$padj<=0.1,]

en2eg=as.list(org.Hs.egENSEMBL2EG)
eg2kegg=as.list(org.Hs.egPATH)
grepKEGG=function(id,mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
kegg=lapply(en2eg,grepKEGG,eg2kegg)


resKEGG=within(geneRank, Row.names<-data.frame(do.call('rbind', strsplit(as.character(Row.names), '.', fixed=TRUE))))
resKEGG$Row.names=resKEGG$Row.names$X1

KEGGResults=list()
sigGenes=genes[genes==1]
for ( i in 1:length(KEGG.sig$category)) {
  #search kegg for the kegg term of interest and filter by differentially expressed genes
  test_term=KEGG.sig$category[i]
  index=sapply(kegg,function(x) test_term %in% x)
  termIDs=names(index[index=="TRUE"])
  
  sig=termIDs[termIDs %in% names(sigGenes)]
  sig=resKEGG[resKEGG$Row.names  %in% sig ,]
  KEGGResults[[test_term]]=sig$Gene.name
}

xx <- as.list(KEGGPATHID2NAME)
names(KEGGResults)=apply(KEGG.sig,1,function(x) xx[[unlist(x[1])]])
KEGGResults=lapply(KEGGResults,function(x) paste(x, sep="", collapse=" ") )
KEGGResults=data.frame(Term= names(KEGGResults),Genes = unlist(KEGGResults),Adj.pvalue=KEGG.sig$padj)
colnames(KEGGResults)[3]=c("Adjusted p-value")
write.table(KEGGResults,file="./Results/KEGGPathways.tab",col.names=T,row.names=F,sep="\t",quote=F)


#map from ensembl to REACTOME
en2eg=as.list(org.Hs.egENSEMBL2EG)
eg2reactome=as.list(reactomeEXTID2PATHID)
grepREACTOME=function(id,mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
reactome=lapply(en2eg,grepREACTOME,eg2reactome)
REACTOME=goseq(pwf,gene2cat=reactome)
REACTOME$padj=p.adjust(REACTOME$over_represented_pvalue,method="BH")
xx <- as.list(reactomePATHID2NAME)
REACTOME$Term=apply(REACTOME,1,function(x) xx[[unlist(x[1])]])
REACTOME.sig=REACTOME[REACTOME$padj<=0.1,]

#work out which proteins are in each catergory
#reactome holds ENSG to reactome ID

#fix the IDs
resReactome=within(geneRank, Row.names<-data.frame(do.call('rbind', strsplit(as.character(Row.names), '.', fixed=TRUE))))
resReactome$Row.names=resReactome$Row.names$X1

reactomeResults=list()
sigGenes=genes[genes==1]

for ( i in 1:length(REACTOME.sig$category)) {
  
  #search reactome for the reactome term of interest and filter by differentially expressed genes
  test_term=REACTOME.sig$category[i]
  index=sapply(reactome,function(x) test_term %in% x)
  termIDs=names(index[index=="TRUE"])
  
  sig=termIDs[termIDs %in% names(sigGenes)]
  sig=resReactome[resReactome$Row.names  %in% sig ,]
  reactomeResults[[test_term]]=sig$Gene.name
}
names(reactomeResults)=REACTOME.sig$Term

reactomeResults=lapply(reactomeResults,function(x) paste(x, sep="", collapse=" ") )
reactomeResults=data.frame(Term= names(reactomeResults),Genes = unlist(reactomeResults),Adj.pvalue=REACTOME.sig$padj)
#reactomeResults$Genes = sapply(lapply(reactomeResults$Genes, strwrap, width=40), paste, collapse="\n")
colnames(reactomeResults)[3]=c("Adjusted p-value")
write.table(reactomeResults,file="./Results/ReactomePathways.tab",col.names=T,row.names=F,sep="\t",quote=F)


##run phenome Express

#source the methods
source("./PhenomeExpress/src/HeterogeneousNetwork.R")
source("./PhenomeExpress/src/RHWN.R")
source("./PhenomeExpress/src/runGIGA.R")
source("./PhenomeExpress/src/runPhenoExpress.R")


#map from ensembl gene to entrez gene to uniprot

#remove the version number on the transcript ID

res=within(res, Row.names<-data.frame(do.call('rbind', strsplit(as.character(rownames(res)), '.', fixed=TRUE))))
colnames(res)[7]="ensembl_gene_id"
res$ensembl_gene_id=res$ensembl_gene_id$X1
egIDs <- stack(mget(na.omit(as.character(res$ensembl_gene_id)), org.Hs.egENSEMBL2EG, ifnotfound = NA))
res=merge(res,egIDs,by.x="ensembl_gene_id",by.y="ind")

#collpase the duplicated entrez IDs based on the highest foldchange
res <- data.table(res)
res=res[, .SD[which.max(abs(log2FoldChange)),], by=values]
res=as.data.frame(res)

#use uniprot and DAVID mapping to get the SwissProt IDs


#Use the David and Uniprot ID maps to match the EntrezID to Swissprot for the PPI network
Young_EnteztoSwiss_via_Uniprot <- read.delim("./GenenamesEntreztoUniprot_via_Uniprot_for8.txt")
Young_EnteztoSwiss_via_Uniprot=Young_EnteztoSwiss_via_Uniprot[,c(1,3)]
colnames(Young_EnteztoSwiss_via_Uniprot)=c("From","To")
Young_EnteztoSwiss_via_David <- read.delim("./GenenamesEntreztoUniprot_via_David_for8.txt", dec=",")
Young_EnteztoSwiss_via_David=Young_EnteztoSwiss_via_David[,1:2]
Young_EnteztoSwiss=rbind(Young_EnteztoSwiss_via_David,Young_EnteztoSwiss_via_Uniprot)
Young_EnteztoSwiss=Young_EnteztoSwiss[!duplicated(Young_EnteztoSwiss),]

dt=merge(res,Young_EnteztoSwiss,by.x="values",by.y="From")
colnames(dt)[9]="name"

#calculate the Pi value for scoring the nodes
dt$absFC=abs(dt$log2FoldChange)
dt$logAdjPval=-log10(dt$padj)
dt$Pi=(dt$absFC*dt$logAdjPval)
dt=na.omit(dt)
#load the HumanConsensusDB PPI network
load("./PhenomeExpress/Networks/ConsensusDB_graph.RData")

#filter the network based on expressed genes
presentList=na.omit(match(dt$name,V(ConsensusDB_graph)$name))
OA.network=induced.subgraph(ConsensusDB_graph,presentList)
OA.network=decompose.graph(OA.network)[[1]]

#filter expression table based on genes in the network
presentList=na.omit(match(V(OA.network)$name,dt$name))
dt=dt[presentList,]
#useful to assign the node with the entrez ID as well - for downstream analysis in cytoscape i.e mapping to genenames or functional annotation
V(OA.network)$EntrezID=as.character(dt$values)

#merge all vertices with the same entrezID but with different protein names
IDs=V(OA.network)$EntrezID
names(IDs)=1:length(IDs)
IDs2=match(IDs,IDs)
IDs=IDs[IDs2]
igraph.options(vertex.attr.comb=list(name="first"))
OA.network2=igraph::contract.vertices(OA.network,mapping=as.integer(names(IDs)))
OA.network2=decompose.graph(OA.network2)[[1]]
OA.network3=get.edgelist(OA.network2)
OA.network3=graph.data.frame(OA.network3,directed=F)
E(OA.network3)$Confidence=E(OA.network2)$Confidence
OA.network=OA.network3

presentList=na.omit(match(V(OA.network)$name,dt$name))
dt=dt[presentList,]




#select the phenotypes from the UberPheno ontology - the Phenomiser tool and manual searching of the ontolgy by relevent keywords is helpful for this
Phenotypes=c("HP:0005086","HP:0001387","MP:0003724","MP:0003436")

set.seed(123)

#run Phenome Express - set inital subnetwork number to 7 to give reasonably sized consensus sub-networks
OAResults=runPhenomeExpress(OA.network,dt,Phenotypes,"Human",max_number=7,sampleSize=10000)

#retrieve the significant sub-networks
subnetworks=OAResults[[1]]

#get the gene ontolgy term for each netwrok
TopGO=list()
for (subnetwork in seq_along(subnetworks)) {
  
  genes=ifelse(dt$name %in% V(subnetworks[[subnetwork]])$name ,1,0)
  names(genes)=dt$ensembl_gene_id
  genes=genes[unique(names(genes))]
  pwf=nullp(genes,"hg19","ensGene",plot.fit=F)
  GO.BP=goseq(pwf,"hg19","ensGene",test.cats=c("GO:BP"),method="Hypergeometric")
  GO.BP=GO.BP[1, "term"]
  TopGO[[subnetwork]]=GO.BP
}

#retrieve the table of p-values
sigTable=OAResults[[2]]

sigTable$GO=unlist(TopGO)
names(subnetworks)<-unlist(TopGO)
save(subnetworks,file="./Results/PhenomeExpress_subnetworks.RData")
NetworkSize=sapply(subnetworks, vcount)
sigTable$NetworkSize=NetworkSize
sigTable$Number=1:nrow(sigTable)
colnames(sigTable)=c("Network Number","Empirical p-value","Top GO Biological Process","Network size")
sigTable=sigTable[,c(1,4,2,3)]

write.table(sigTable,file="./Results/PhenomeExpressSummary.tab",row.names=F,col.names=T,sep="\t",quote=F)


z=getHeterogeneousNetwork(OA.network,"Human")[["genePheno"]] # note contains all proteins - including ones not present in network
phenoAnnotated=z[rownames(z) %in% Phenotypes,]
phenoAnnotated=phenoAnnotated[,colSums(phenoAnnotated)>0]
phenoAnnotated=colnames(phenoAnnotated)

#send all the sub-networks from PhenomeExpress to cytoscape
#colours the nodes according to the fold change
#black border if directly annotated to seed phenotype 


  
  V(OA.network)$EntrezID=as.character(dt$values)
  
  Genes$ID<- do.call('rbind', strsplit(as.character(Genes$ID),'.',fixed=TRUE))[,1]
  dt2=merge(dt,Genes,by.x="ensembl_gene_id",by.y="ID")
  Swiss2GeneSymbol<-dt2[,c(9,13)]
#  save(Swiss2GeneSymbol,file="./Results/Swiss2GeneSymbol.RData")
  

if (visualise == TRUE) {

  for(i in 1:length(subnetworks)) {
    presentList=na.omit(match(V(subnetworks[[i]])$name,V(OA.network)$name))
    tempGraph=induced.subgraph(OA.network,presentList)
    FC=dt2[na.omit(match(V(tempGraph)$name,dt2$name)),]
    V(tempGraph)$logFC=FC$log2FoldChange
    V(tempGraph)$GeneSymbol=as.character(FC$Gene.symbol)
    V(tempGraph)$GeneName=as.character(FC$Gene.Name)
    seedAnnotatedGenes=ifelse(V(tempGraph)$name %in% phenoAnnotated,1,0)
    V(tempGraph)$Seed=seedAnnotatedGenes
    
    
    #do the network creation stuff
    
    #convert the igraph object to a graphNEL object and intialise the attributes
    tempGraph.NEL=igraph.to.graphNEL(tempGraph)
    tempGraph.NEL=initEdgeAttribute(tempGraph.NEL,"Confidence","numeric",0)
    tempGraph.NEL=initEdgeAttribute(tempGraph.NEL,"weight","numeric",0)
    tempGraph.NEL=initNodeAttribute(tempGraph.NEL,"logFC","numeric",0)
    tempGraph.NEL=initNodeAttribute(tempGraph.NEL,"Seed","numeric",0)
    tempGraph.NEL=initNodeAttribute(tempGraph.NEL,"EntrezID","char",0)
    tempGraph.NEL=initNodeAttribute(tempGraph.NEL,"GeneSymbol","char",0)
    tempGraph.NEL=initNodeAttribute(tempGraph.NEL,"GeneName","char",0)
    
    nodeDataDefaults(tempGraph.NEL, "label") <- "name"
    nodeData(tempGraph.NEL,V(tempGraph)$name,"label") = V(tempGraph)$name
    tempGraph.NEL=initNodeAttribute(tempGraph.NEL,"label","char","name")
    
    
    #Open the cytoscape window and send the graph
    cw1 <- new.CytoscapeWindow (paste("PhenoExpressadded",as.character(i),sep=""), graph=tempGraph.NEL)
    #display the graph
    displayGraph (cw1)
    #select the layout
    layoutNetwork (cw1, layout.name='force-directed')
    
    #colour according to the logFC
    control.points <- c(-5,0,5)
    node.colors <- c ("#00AA00", "#00FF00", "#FFFFFF", "#FF0000", "#AA0000")
    setNodeColorRule (cw1, node.attribute.name='logFC', control.points, node.colors, mode='interpolate')
    
    
    setDefaultBackgroundColor (cw1, '#FFFFFF')
    
    #set the nodeborder to correspond to the seed phenotype annotated genes
    data.values <- c ("1", "0")
    line.widths = c ("15","1")
    setNodeBorderWidthRule (cw1, 'Seed', data.values, line.widths)
  }
}

##################################
#Comparison with existing studies
##################################

#compare our results with Ramos et al's results

#load the raw data as a GEO gset
load("./GSE57218.RData")
lumi.N <- lumiN(gset,method="rsn")
detection <- read.delim("./GSE57218_Non-normalized_data.txt")
rownames(detection)=detection$ID_REF
detection_index=grep("detection",colnames(detection))
detection=detection[detection_index]
detection=ifelse(detection<=0.05,1,0)
detection=detection[rowSums(detection)>32,]
IDs=rownames(detection)
chipVersion=getChipInfo(lumi.N)[[1]]
IDs=probeID2nuID(IDs,lib.mapping="lumiHumanIDMapping",species="Human",chipVersion=chipVersion)
expressed=IDs[,7]


#convert the probe ID to NUIDs
lumi.N = addNuID2lumi(lumi.N,lib.mapping="lumiHumanIDMapping")
#get the probe ID and intensities for each sample
dataMatrix <- exprs(lumi.N)
dim(dataMatrix)
dataMatrix=dataMatrix[rownames(dataMatrix) %in% expressed,]

#filter out unannotated probes
dataMatrix=dataMatrix[!is.na(getSYMBOL(rownames(dataMatrix), 'lumiHumanAll.db')),]
dataMatrix.complete=dataMatrix
#remove the healthy samples
dataMatrix=dataMatrix[,-c(1,2)]
dataMatrix=dataMatrix[,-c(68:71)]
dataMatrix=dataMatrix[,-c(29)]

#define the experimental conditions (factors)
sampleType=as.factor(c(rep(c("OA","Perserved"),33)))
sampleType=relevel(sampleType,"Perserved")
patient=as.factor(rep(1:33, each=2))

design <- model.matrix(~patient + sampleType)
#fit the linear model
fit <- lmFit(dataMatrix, design)

#calculate pvalues
fit <- eBayes(fit)

#calculate BH correction p values and store results table
results=topTable(fit,coef="sampleTypeOA",number=Inf)

#Annotate the probes with Entrez gene IDs
genes=as.data.frame(getEG(rownames(results), 'lumiHumanAll.db' ))
colnames(genes)=c("EntrezID")
results=merge(results,genes,by="row.names")

#collpase the duplicated entrez IDs based on the highest foldchange
results <- data.table(results)
results=results[, .SD[which.max(abs(logFC)),], by=EntrezID]
results=as.data.frame(results)
results=results[order(-abs(results$logFC)),]

#Annotate with the genename
genes=as.data.frame(getSYMBOL(as.character(results$Row.names), 'lumiHumanAll.db' ))
colnames(genes)=c("GeneName")
results=merge(results,genes,by.y="row.names",by.x="Row.names")
results=results[order(-abs(results$logFC)),]


resultsRamos<-results
save(resultsRamos,file="./Results/resultsRamos.RData")


#analyse the snelling et al data


targets <- read.delim("./microarray_results/targets.txt", header=T)
setwd("./microarray_results")
images <- read.maimages(targets,source="agilent")
setwd("..")

images.processed <- backgroundCorrect(images, method="normexp", offset=50)
images.processed=normalizeWithinArrays(images.processed,method="loess")
images.processed <- normalizeBetweenArrays(images.processed, method="quantile")
images.processed <- avereps(images.processed, ID=images.processed$genes$GeneName)

xx <- as.list(org.Hs.egSYMBOL)
images.processed=images.processed[images.processed$genes$ControlType!=1,]
images.processed=images.processed[images.processed$genes$GeneName %in% xx,]


design=modelMatrix(targets,ref="undamaged")
fit <- lmFit(images.processed, design)
fit <- eBayes(fit)
results.Snelling=topTable(fit,coef="damaged",number=Inf)



#get the sig genes
sigGenesRamos=results[abs(results$logFC)>=0.58,]
sigGenesRamos=na.omit(sigGenesRamos[sigGenesRamos$adj.P.Val<=0.1,])

Oursig=geneRank[abs(geneRank$log2FoldChange)>=0.58,]
Oursig=na.omit(Oursig[Oursig$padj<=0.1,])

Snellingsig=results.Snelling[abs(results.Snelling$logFC)>=0.58,]
Snellingsig=na.omit(Snellingsig[Snellingsig$adj.P.Val<=0.1,])


vennList=list(ours=Oursig$Gene.name,Ramos=sigGenesRamos$GeneName,Snelling=Snellingsig$GeneName)


venn.diag=venn.diagram(vennList,fill = c("red", "green","lightskyblue"),alpha = c(0.5, 0.5,0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3,  filename="./Results/Venn.jpg")

#get intersectection genes and extract the data from each table to allow plotting

chosenGenes=intersect(sigGenesRamos$GeneName,intersect(Snellingsig$GeneName,Oursig$Gene.name))
chosenGeneTable=merge(Oursig,Snellingsig,by.x="Gene.name",by.y="GeneName")
chosenGeneTable=merge(chosenGeneTable,sigGenesRamos,by.x="Gene.name",by.y="GeneName")

chosenGeneTable <- data.table(chosenGeneTable)

chosenGeneTable=as.data.frame(chosenGeneTable)
chosenGeneTable=data.frame(GeneName=chosenGeneTable$Gene.name,Ours=chosenGeneTable$log2FoldChange,Ramos=chosenGeneTable$logFC.y,Snelling=chosenGeneTable$logFC.x)

chosenGeneTable.long=melt(chosenGeneTable,value.name="logFC")

graph.genes<- ggplot(chosenGeneTable.long,aes(GeneName,logFC,fill=as.factor(variable)))+ geom_bar(position="dodge",stat="identity") + facet_wrap(~GeneName,nrow=3,scales = "free_x") + theme(axis.ticks = element_blank(), axis.title=element_text(size=18,face="bold"), axis.text.x = element_blank(),strip.text.x = element_text(size = 14,face="bold"),legend.position="bottom",legend.text=element_text(size=16,face="bold"),legend.title = element_text(colour="red", size = 18, face = "bold")) +  scale_y_continuous(name="log2 Fold Change" ,breaks=seq(-4, 4, 0.5)) + scale_x_discrete(name="Gene Name")
graph.genes<- graph.genes + scale_fill_manual(values= c("red", "green","lightskyblue"),name="",
                                              breaks=c("Ours", "Ramos", "Snelling"),
                                              labels=c("Present Study", "Ramos et al", "Snelling et al")) + theme(panel.background = element_rect(fill = "white",colour = "grey"), plot.background = element_rect(fill = "white",colour = "grey"))


ggsave("./Results/comparison.jpg",width=12,height=9.3)



chosen2=intersect(Oursig$Gene.name,sigGenesRamos$GeneName)
chosen2=setdiff(chosen2,chosenGenes)
chosenGeneTable=merge(Oursig,sigGenesRamos,by.x="Gene.name",by.y="GeneName")
chosenGeneTable=data.frame(GeneName=chosenGeneTable$Gene.name,Ours=chosenGeneTable$log2FoldChange,Ramos=chosenGeneTable$logFC)
chosenGeneTable=chosenGeneTable[ chosenGeneTable$GeneName %in% chosen2,]
chosenGeneTable.long=melt(chosenGeneTable,value.name="logFC")

graph.genes<- ggplot(chosenGeneTable.long,aes(GeneName,logFC,fill=as.factor(variable)))+ geom_bar(position="dodge",stat="identity") + facet_wrap(~GeneName,nrow=6,scales = "free_x") + theme(axis.ticks = element_blank(), axis.title=element_text(size=18,face="bold"), axis.text.x = element_blank(),strip.text.x = element_text(size = 14,face="bold"),legend.position="bottom",legend.text=element_text(size=16,face="bold"),legend.title = element_text(colour="red", size = 18, face = "bold")) +  scale_y_continuous(name="log2 Fold Change" ,breaks=seq(-4, 4, 0.5)) + scale_x_discrete(name="Gene Name")
graph.genes<- graph.genes + scale_fill_manual(values=c("red", "green"),name="",
                                              breaks=c("Ours", "Ramos"),
                                              labels=c("Present Study", "Ramos et al")) + theme(panel.background = element_rect(fill = "white",colour = "grey"), plot.background = element_rect(fill = "white",colour = "grey"))

ggsave("./Results/comparisonRamos.jpg",width=12,height=9.3)


chosen3=intersect(Oursig$Gene.name,Snellingsig$GeneName)
chosen3=setdiff(chosen3,chosenGenes)
chosenGeneTable=merge(Oursig,Snellingsig,by.x="Gene.name",by.y="GeneName")
chosenGeneTable=data.frame(GeneName=chosenGeneTable$Gene.name,Ours=chosenGeneTable$log2FoldChange,Snelling=chosenGeneTable$logFC)
chosenGeneTable=chosenGeneTable[ chosenGeneTable$GeneName %in% chosen3,]
chosenGeneTable.long=melt(chosenGeneTable,value.name="logFC")

graph.genes<- ggplot(chosenGeneTable.long,aes(GeneName,logFC,fill=as.factor(variable)))+ geom_bar(position="dodge",stat="identity") + facet_wrap(~GeneName,nrow=6,scales = "free_x") + theme(axis.ticks = element_blank(), axis.title=element_text(size=18,face="bold"), axis.text.x = element_blank(),strip.text.x = element_text(size = 14,face="bold"),legend.position="bottom",legend.text=element_text(size=16,face="bold"),legend.title = element_text(colour="red", size = 18, face = "bold")) +  scale_y_continuous(name="log2 Fold Change" ,breaks=seq(-4, 4, 0.5)) + scale_x_discrete(name="Gene Name")
graph.genes<- graph.genes + scale_fill_manual(values=c("red", "lightskyblue"),name="",
                                              breaks=c("Ours", "Snelling"),
                                              labels=c("Present Study", "Snelling et al")) + theme(panel.background = element_rect(fill = "white",colour = "grey"), plot.background = element_rect(fill = "white",colour = "grey"))

ggsave("./Results/comparisonSnelling.jpg",width=14,height=9.3)


#intersection pathway analysis

OursandRamos<-intersect(Oursig$Gene.name,sigGenesRamos$GeneName)
OursandSnelling<-intersect(Oursig$Gene.name,Snellingsig$GeneName)
RamosandSnelling<-intersect(sigGenesRamos$GeneName,Snellingsig$GeneName)

intersectionGenes<-c(OursandRamos,OursandSnelling,RamosandSnelling)

genes=ifelse(geneRank$Gene.name %in% intersectionGenes,1,0)
names(genes)=transcriptIDs[,1]
pwf=nullp(genes,"hg19","ensGene")




#map from ensembl to REACTOME
en2eg=as.list(org.Hs.egENSEMBL2EG)
eg2reactome=as.list(reactomeEXTID2PATHID)
grepREACTOME=function(id,mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
reactome=lapply(en2eg,grepREACTOME,eg2reactome)
REACTOME=goseq(pwf,gene2cat=reactome)
REACTOME$padj=p.adjust(REACTOME$over_represented_pvalue,method="BH")
xx <- as.list(reactomePATHID2NAME)
REACTOME$Term=apply(REACTOME,1,function(x) xx[[unlist(x[1])]])
REACTOME.sig=REACTOME[REACTOME$padj<=0.1,]

#work out which proteins are in each catergory
#reactome holds ENSG to reactome ID

#fix the IDs
resReactome=within(geneRank, Row.names<-data.frame(do.call('rbind', strsplit(as.character(Row.names), '.', fixed=TRUE))))
resReactome$Row.names=resReactome$Row.names$X1

reactomeResults=list()
sigGenes=genes[genes==1]

for ( i in 1:length(REACTOME.sig$category)) {
  
  #search reactome for the reactome term of interest and filter by differentially expressed genes
  test_term=REACTOME.sig$category[i]
  index=sapply(reactome,function(x) test_term %in% x)
  termIDs=names(index[index=="TRUE"])
  
  sig=termIDs[termIDs %in% names(sigGenes)]
  sig=resReactome[resReactome$Row.names  %in% sig ,]
  reactomeResults[[test_term]]=sig$Gene.name
}
names(reactomeResults)=REACTOME.sig$Term

reactomeResults=lapply(reactomeResults,function(x) paste(x, sep="", collapse=" ") )
reactomeResults=data.frame(Term= names(reactomeResults),Genes = unlist(reactomeResults),Adj.pvalue=REACTOME.sig$padj)
#reactomeResults$Genes = sapply(lapply(reactomeResults$Genes, strwrap, width=40), paste, collapse="\n")
colnames(reactomeResults)[3]=c("Adjusted p-value")
write.table(reactomeResults,file="./Results/IntersectionReactomePathways.tab",col.names=T,row.names=F,sep="\t",quote=F)

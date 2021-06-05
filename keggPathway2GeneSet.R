###############################################################################
# title: keggPathway2GeneSet
# author: Xiaowei
# time: Jan.7 2021
# function: kegg.download, path.to.geneSet, kegg.geneSet
###############################################################################


###############################################################################
# install packages 
###############################################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", dependencies = TRUE)

requiredPackages <- c("KEGGREST", "GSEABase")			  
newPackages <- requiredPackages[!(requiredPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) BiocManager::install(newPackages, ask = TRUE)



suppressPackageStartupMessages(library(KEGGREST))
suppressPackageStartupMessages(library(GSEABase))


###############################################################################
# Function -- kegg.download 
# description: download pathway from kegg
# input: One or more KEGG identifiers
# output: A list wrapping a KEGG flat file.
###############################################################################


kegg.download <- function(pathway.list){
  #因为keggGet()一次最大查询10条，所以这里先将pathway.list分成10*n个，每10个放在pathway.x中
  #确定pathway.x长度
  if (length(pathway.list)%%10 == 0){
    pathway.x.len <- length(pathway.list)%/%10
  }else{
    pathway.x.len <- length(pathway.list)%/%10 + 1
  }
  #pathway.x,每个中包含10个pathway
  pathway.x <- vector(mode = "list", length = pathway.x.len)
  for (i in 1:length(pathway.x)){
    min.x <- 10*(i -1)+1
    max.x <- 10*i
    pathway.x[[i]] <- names(pathway.list)[min.x:max.x]
    rm(min.x,max.x)
  }
  #下载n次，每次下载10个pathway
  pathway.kegg <- list()
  for (i in 1:pathway.x.len){
    pathway.kegg.x <- keggGet(pathway.x[[i]])
    pathway.kegg <- append(pathway.kegg, pathway.kegg.x)
    rm(pathway.kegg.x)
  }
  
  names(pathway.kegg) <- unlist(lapply(names(pathway.list), function(x){trimws(strsplit(x, ':', fixed = TRUE)[[1]][2])} ))
  
  return(pathway.kegg)
  
}

###############################################################################
# Function -- path.to.GeneSet
# description: make the result of kegg.download to GeneSet
# input: the result of kegg.download
# Output: A GeneSet object 
###############################################################################
path.to.geneSet <- function(kegg.path){
  genes <- kegg.path$GENE
  
  if(!is.null(genes)){
    genelist_entrez <- genes[1:length(genes)%%2 ==1]  #entrez
    
    gs <- GeneSet(geneIds = as.character(genelist_entrez), 
                  geneIdType = EntrezIdentifier(), 
                  #organism = 'hsa', 
                  collectionType = KEGGCollection(),
                  #longDescription = kegg.path$DESCRIPTION,
                  #shortDescription = names(kegg.path$PATHWAY_MAP),
                  setName = kegg.path$PATHWAY_MAP, 
    )
  }else{gs = NULL;print("Not had genes.")}
  
  return(gs)
}

###############################################################################
# Function -- kegg.geneSet
# description: download pathway of one organism from KEGG and make it as GeneSetCollection 
# input: 
# organism: a KEGG organism code (list via keggList("organism"))
# outputfile: name of GMT format file
# output: A GeneSetCollection object or/and GMT file, gene ID is Entrez 
###############################################################################

kegg.geneSet <- function(organism = "hsa", outputfile = NULL){
  pathway.list <- keggList("pathway", organism) #获取所有pathway的名称和kegg标识符
  
  kegg.path <- kegg.download(pathway.list = pathway.list) #download pathway
  
  kegg.geneSetCollection <- mapply(path.to.geneSet, kegg.path) #pathway to GeneSet
  
  null.index <- unlist(lapply(kegg.geneSetCollection, is.null)) #remove NULL
  kegg.geneSetCollection <- GeneSetCollection(kegg.geneSetCollection[!null.index])
  
  #导出为gmt文件
  if (!is.null(outputfile)){toGmt(x=kegg.geneSetCollection, con = outputfile)}
  
  return(kegg.geneSetCollection)
}

###############################################################################
# test
###############################################################################
# pathway.list <- keggList("pathway", organism = "hsa")
# kegg.path <- kegg.download(pathway.list = pathway.list[1:11])
# kegg.geneSet <- path.to.geneSet(kegg.path[[1]])
# 
# hsa.geneSet <- kegg.geneSet("hsa") # Homo sapiens (human)
# ptr.geneSet <- kegg.geneSet("ptr", outputfile = "ptr.gmt") # Pan troglodytes (chimpanzee)
# ggo.geneSet <- kegg.geneSet("ggo") # Gorilla gorilla gorilla (western lowland gorilla)
# pon.geneSet <- kegg.geneSet("pon") # Nomascus leucogenys (northern white-cheeked gibbon)

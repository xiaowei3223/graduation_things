###############################################################################
# title: msigdb.geneSet
# author: Xiaowei
# time: Jan.8 2021
# function: msigdb.geneSet
###############################################################################

###############################################################################
# install packages 
###############################################################################
if (!requireNamespace("msigdbr", quietly = TRUE)){install.packages("msigdbr", dependencies = TRUE)}
if (!requireNamespace("GSEABase", quietly = TRUE)){BiocManager::install("GSEABase", dependencies = TRUE)}

###############################################################################
# Function -- msigdbr.geneSet
# description: download pathway of one organism from msigdb and make it as GeneSetCollection 
# input: 
# species: Species name, such as Homo sapiens or Mus musculus,See more  via msigdbr_species()
# category: MSigDB collection abbreviation, such as H or C1. See more via msigdbr_collections()
# subcategory: MSigDB sub-collection abbreviation, such as CGP or BP. See more via msigdbr_collections()
# geneIdType: Default as "entrez". one of "entrez" and "symbol"
# outputfile: name of GMT format file
# output: A GeneSetCollection object or/and GMT file
###############################################################################
msigdb.geneSet <- function(species, category = NULL, subcategory = NULL, geneIdType ="entrez", outputfile = NULL){
  suppressPackageStartupMessages(library(msigdbr))
  suppressPackageStartupMessages(library(GSEABase))
  
  # 下载基因集数据框
  gs = msigdbr(species = species, category = category, subcategory = subcategory)
  
  if (geneIdType == "entrez"){
    genes = gs$entrez_gene
    geneIdsType = EntrezIdentifier()
  }else if (geneIdType == "symbol"){
    genes = gs$gene_symbol
    geneIdsType = SymbolIdentifier()
  }
  # 根据GeneSet name 转变成list
  gs.names <- unique(gs$gs_name)
  gs.list <- lapply(gs.names, function(x){
    unique(genes[which(gs$gs_name == x)])
  })
  names(gs.list) = gs.names
  
  # 变成GeneSetCollection
  gs.geneSetCollection <- GeneSetCollection(mapply(function(x,y){
    genes = as.character(unlist(x))
    GeneSet(geneIds = genes,
            geneIdType = geneIdsType,
            collectionType = NullCollection(),
            setName = y
    )
  }, gs.list, gs.names))
  
  #导出为gmt文件
  if (!is.null(outputfile)){toGmt(x=gs.geneSetCollection, con = outputfile)}
  
  return(gs.geneSetCollection)
}



###############################################################################
# test
###############################################################################
# hsa.msigdb <- msigdb.geneSet(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
# hsa.msigdb <- msigdb.geneSet(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG", outputfile = "msigdb.human.kegg.gmt")
# hsa.msigdb <- msigdb.geneSet(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG", geneIdType = "symbol")


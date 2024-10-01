## ----loadBiobase----------------------------------------------------------------
library(Biobase)


## ----ExpressionSet, fig.cap="Structure of the <tt>ExpressionSet</tt> class, showing its slots and their meaning. Reproduced from Klaus, B., & Reisenauer, S. (2018)", echo=FALSE----
knitr::include_graphics("images/Structure-of-Bioconductors-ExpressionSet-class.png")


## ----simulateData---------------------------------------------------------------
expressionValues <- matrix (rnorm (300), nrow=30)
colnames(expressionValues) <- paste0("sample",1:10)
head(expressionValues)


## ----simulateCovariates---------------------------------------------------------
targets <- data.frame(sampleNames = paste0("sample",1:10),
                      group=c(paste0("CTL",1:5),paste0("TR",1:5)),
                      age = rpois(10, 30), 
                      sex=as.factor(sample(c("Male", "Female"),10,replace=TRUE)),
                      row.names=1)
head(targets, n=10)


## ----simulateGeneInfo-----------------------------------------------------------
myGenes <-  paste0("gene",1:30)


## ----simulateInfo---------------------------------------------------------------
myInfo=list(myName="Alex Sanchez", 
            myLab="Bioinformatics Lab",
            myContact="alex@somemail.com", 
            myTitle="Practical Exercise on ExpressionSets")
show(myInfo)


## ----PCA1-----------------------------------------------------------------------
pcs <- prcomp(expressionValues)
names(pcs)
barplot(pcs$sdev)
plot(pcs$rotation[,1], pcs$rotation[,2], 
     main="Representation of first two principal components")
text(pcs$rotation[,1], pcs$rotation[,2], targets$group, cex=0.8, pos=3)


## ----sortGenesByVar-------------------------------------------------------------
variab <- apply(expressionValues, 1, sd)
orderedGenes <- myGenes[order(variab, decreasing=TRUE)]
head(variab[order(variab, decreasing=TRUE)])
head(orderedGenes)


## ----subsetExpressions----------------------------------------------------------
newExpress<- expressionValues[,-9]
newTargets <- targets[-9,]
wrongNewTargets <- targets [-10,]


## ----creaExpressionSet1---------------------------------------------------------
myEset <- ExpressionSet(expressionValues)
class(myEset)
show(myEset)


## ----AnnotatedDataFrame2--------------------------------------------------------
columnDesc <-  data.frame(labelDescription= c("Treatment/Control", 
                                                "Age at disease onset", 
                                                "Sex of patient (Male/Female"))
myAnnotDF <- new("AnnotatedDataFrame", data=targets, varMetadata= columnDesc)
show(myAnnotDF)


## ----addAnnotatedDataFrame------------------------------------------------------
phenoData(myEset) <- myAnnotDF


## ----creaEset2------------------------------------------------------------------
myEset <- ExpressionSet(assayData=expressionValues, phenoData=myAnnotDF)
show(myEset)


## ----creaEset3------------------------------------------------------------------
myEset <- ExpressionSet(assayData=expressionValues,
                        phenoData=myAnnotDF,
                        featureNames =myGenes)
# show(myEset)


## ----label=MIAME----------------------------------------------------------------
myDesc <- new("MIAME", name= myInfo[["myName"]],
            lab= myInfo[["myLab"]],
            contact= myInfo[["myContact"]] ,
            title=myInfo[["myTitle"]])
print(myDesc)


## ----creaEset4------------------------------------------------------------------
myEset <- ExpressionSet(assayData=expressionValues,
                        phenoData=myAnnotDF,
                        fetureNames =myGenes,
                        experimentData = myDesc)
# show(myEset)


## ----usingExpressionSets--------------------------------------------------------
dim(exprs(myEset))
class(phenoData(myEset))
class(pData(phenoData(myEset)))
head(pData(phenoData(myEset)))
head(pData(myEset))


## ----smallEset------------------------------------------------------------------
smallEset <- myEset[1:15,c(1:3,6:8)]
dim(exprs(smallEset))
dim(pData(smallEset))
head(pData(smallEset))
all(colnames(exprs(smallEset))==rownames(pData(smallEset)))


## ----creaEset5------------------------------------------------------------------
youngEset <- myEset[,pData(myEset)$age<30]
dim(exprs(youngEset))
head(pData(youngEset))


## ----GEOquery-------------------------------------------------------------------
if (!require(GEOquery)) {
  BiocManager::install("GEOquery")
}
require(GEOquery)
gse <- getGEO("GSE27174", GSEMatrix=TRUE, AnnotGPL=TRUE)


## ----gseObject------------------------------------------------------------------
class(gse)
names(gse)
length(gse)
gse[[1]]
esetFromGEO <- gse[[1]]


## ----esetFromGEO----------------------------------------------------------------
head(exprs(esetFromGEO))


## ----esetfromGEO2---------------------------------------------------------------
colnames(pData(esetFromGEO))
pData(esetFromGEO)[,39:40]


## ----getGSD---------------------------------------------------------------------
gds <- getGEO("GDS4155")


## ----gdsObject------------------------------------------------------------------
class(gds)
slotNames(gds)


## ----gdsMetaData----------------------------------------------------------------
head(Meta(gds))


## ----esetFromGDS----------------------------------------------------------------
eset <- GDS2eSet(gds,do.log2=FALSE)
eset


## ----sessionInfo----------------------------------------------------------------
sessionInfo()


## ----creaIndex------------------------------------------------------------------
# An "index.html" file is created to allow visualitzation in the web using github pages
file.copy(from="Introduction_2_Bioc_classes_4_tabular_data.html", to="index.html", overwrite=TRUE)


## ----eval=FALSE-----------------------------------------------------------------
## # The R code for the document can be extracted from the document with the
## # knitr::purl() command
## # knitr::purl("Introduction_2_Bioc_classes_4_tabular_data.qmd")


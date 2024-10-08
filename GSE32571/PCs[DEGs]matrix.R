library(readr)
library(GEOquery)
library(hgu133plus2.db)
library(hgu133acdf)
library(limma)
library(hgu133plus2cdf)
library(GSEABase)
library(GOstats)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(curl)
library(RCurl)
library(affy)
library(readr)
library(hgu133a.db)
library(genefilter)
library(multtest)
library(affyPLM)
library(pheatmap)
library(pacman)
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%

my_id <- "GSE32571"
gse <- getGEO(my_id)

## check how many platforms used
length(gse)

gse2 <- gse[[1]]
gse2

pData(gse2) ## print the sample information
fData(gse2) ## print the gene annotation
exprs(gse2) ## print the expression data
pData(gse2)$data_processing[1] #Check the normalisation

## exprs get the expression levels as a data frame and get the distribution
summary(exprs(gse2))
boxplot(exprs(gse2),outline=FALSE,main = "Box Plot [GSE32571]",col=seq(1:23))

expdata <- exprs(gse2)

{
expdata <-log2(exprs(gse2))
boxplot(expdata,outline=FALSE,main = "Box Plot log2[GSE30994]",col=seq(1:23))
summary(expdata)
head(expdata)
head(exprs(gse2))
}

{
  #**"Depend on the data"
  #From this output we clearly see that the values go beyond 16,
  #so we will need to perform a log2 transformation
  colnames(exprs(gse2))
  colnames(exprs(gse2)[,1:24])
  expdata <- bind_cols(log2(exprs(gse2)[,1:24]), exprs(gse2)[,25:dim(exprs(gse2))[2]])
  
  boxplot(expdata,outline=FALSE,main = "Box Plot log(GSE62431)",col=seq(1:23))
  summary(expdata)
  
  
}


#====================[Explore The Data]=========================================


#***For your own data, you will have to decide which columns will be useful in the analysis
sampleInfo <- pData(gse2)
colnames(sampleInfo)

## source_name_ch1 seems to contain factors we might need for the analysis. Let's pick just those columns
sampleInfo <- select(sampleInfo, "disease stage:ch1")
head(sampleInfo)

## Optionally, rename to more convenient column names
sampleInfo <- rename(sampleInfo, Condition = "disease stage:ch1")


#***Unsupervised analysis is a good way to get an understanding of the sources of variation in the data.
#It can also identify potential outlier samples.

corMatrix <- cor(expdata,use="c") ## argument use="c" stops an error if there are any missing data points
rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix,
         annotation_col=sampleInfo)    

#**** A complementary approach is to use Principal Components Analysis (PCA).
## MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX
pca <- prcomp(t(exprs(gse2)),)
## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=Condition,label=paste("Condition", Condition))) + geom_point()



{
  #*** Exporting the data
  
  features <- fData(gse2)
  View(features)
  colnames(features)
  ### Look at the features data frame and decide the names of the columns you want to keep
  features <- select(features,"Gene Symbol","Gene Title")
  full_output <- cbind(features,exprs(gse2))
  write_csv(full_output, path="gse_full_output.csv")
  
  }


#==============================[limma]=========================================

#By far the most-popular package for performing differential expression is limma.
#The design matrix is a matrix of 0 and 1s; one row for each sample and one column for each sample group.
#A 1 in a particular row and column indicates that a given sample (the row) belongs to a given group (column).

design <- model.matrix(~0+sampleInfo$Condition)
design

## the column names are a bit ugly, so we will rename
colnames(design)
colnames(design) <-  c("Healthy","Cancer")


#***Coping with outliers
#It is tempting to discard any arrays which seem to be outliers prior to differential expressions. 
#However, this is done at the expense of sample-size which could be an issue for small experiments. 
#A compromise, which has been shown to work well is to calculate **weights to define the reliability of each sample.
aw <- arrayWeights(expdata,design)
aw

#***The lmFit function is used to fit the model to the data. 
#The result of which is to estimate the expression level in each of the groups that we specified.


fit <- lmFit(expdata, design,
             weights = aw)
head(fit$coefficients)

##In order to perform the differential analysis, we have to define the contrast that we are interested in. In our case we only have two groups and one contrast of interest.
##Multiple contrasts can be defined in the makeContrasts function

contrasts <- makeContrasts(Cancer - Healthy,
                           Healthy - Cancer, levels=design)

## can define multiple contrasts
## e.g. makeContrasts(Group1 - Group2, Group2 - Group3,....levels=design)

fit2 <- contrasts.fit(fit, contrasts)

##Finally, apply the empirical Bayes’ step 
##to get our differential expression statistics and p-values.

fit2 <- eBayes(fit2)
head(fit2$coefficients)

topTable(fit2, coef=1) ## topTable(fit2, coef=2)
### to see the results of the second contrast (if it exists)

#If we want to know how many genes are differentially-expressed overall
#we can use the decideTests function.
table(decideTests(fit2))

## to perform adjust.method="fdr" FDR means ***Benjamini & Hochberg (1995) ("BH")***
top1 <- limma::topTable(fit2, coef=1, number= nrow(fit2), adjust.method="fdr",sort.by="none")
top2 <- limma::topTable(fit2, coef=2, number=nrow(fit2), adjust.method="BH",sort.by="none")


#========================[Annotating the data]=============================

anno <- fData(gse2)
colnames(anno)
head(anno$ID)
top1[7] <- select(anno,"Entrez_Gene_ID") 
top1[8] <- select(anno,"Symbol")

#========================[Manipulating with Probes]=============================

# [1] select the Comparison items
full_results1 <- as.data.frame(top1) #Hep3B_ATAD2 - Hep3B_control
colnames(full_results1)
{
  
  full_results2 <- topTable(fit2, coef=2, number=Inf) #PLC_ATAD2 - PLC_control
  full_results3 <- topTable(fit2, coef=3, number=Inf) #"SNU449_ATAD2","SNU449_control"
  full_results4 <- topTable(fit2, coef=4, number=Inf)
  head(full_results1)
  
}

# [2] Save the the comparison statistic results
write.table(full_results1, file="statistic_results.txt", sep="\t", row.names=T, col.names=NA)
stat_data <- read_delim("statistic_results.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# [3] do i need the probe id again?  no, then let's drop it
stat_data=stat_data[,-1]

# [4] remove nulls : some probes were not mapped to any gene symbol
full_results1 = stat_data[ ! is.na(stat_data$Symbol),]

#[5] check duplciation of of gene symbols?  
x=duplicated(full_results1$Symbol)  
sum(x)

#[6] ### yes .. why ? probesets?  solutions : aggregation
exp.data=full_results1[-dim(full_results1)[2]]
exp.data=apply(full_results1,2, as.numeric)

# [7] #### remove duplication
exp.data.agg= aggregate(exp.data, list(full_results1$Symbol),FUN=mean)
names(exp.data.agg)

# [8] rename the row name -> Gene Symbol 
rownames(exp.data.agg)=exp.data.agg$Group.1
exp.data.agg=exp.data.agg[,c(-1,-9)]

# [9]
processed_data <- exp.data.agg
write.table(processed_data,file = "GSE69223_stats.txt",row.names = T,col.names = T,quote = F)


#=======================[Exporting & Visualization Of DE results]================================

#[1] Exporting

degs.res= processed_data[abs(processed_data$logFC) >= 1  & processed_data$adj.P.Val <= 0.05, ]
degs.res_UP=processed_data[processed_data$logFC >= 1   & processed_data$adj.P.Val <= 0.05,]  
degs.res_DOWN=processed_data[processed_data$logFC  <= -1   & processed_data$adj.P.Val <= 0.05,] 

write.table(degs.res,file = "DEGs[ID&names].txt",row.names = T,col.names = T,quote = F)
write.table(rownames(degs.res),file = "DEGs.txt",row.names = F,col.names = F,quote = F)
write.table(rownames(degs.res_UP),file = "DEGs_UP.txt",row.names = F,col.names = F,quote = F)
write.table(rownames(degs.res_DOWN),file = "DEGs_DOWN.txt",row.names = F,col.names = F,quote = F)

#[2] Visualization
# [A] ----------- Volcano (UP/Down)-----------  

res2 <- processed_data %>%
  mutate(gene_type = case_when(logFC >= 1  & adj.P.Val <= 0.05 ~ "up",
                               logFC <= -1  & adj.P.Val <= 0.05 ~ "down",
                               TRUE ~ "ns"))   
cols <- c("up" = "#FF0000", "down" = "#00FF00", "ns" = "grey") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

res2 %>%
  ggplot(aes(x = (logFC),
             y = -log10(adj.P.Val),
             fill = gene_type,    
             size = gene_type,
             alpha = gene_type)) + 
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black") + 
  geom_hline(yintercept = (0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-4, 4, 2)),       
                     limits = c(-4, 4)) 

#[2] Visualization
# [B] ----------- HeatMap(DEGs)-----------  

# [ Selecting the set]
full_results <- topTable(fit2, coef=1, number=Inf)
full_results <- tibble::rownames_to_column(full_results,"ID")
# [1]
my_genes <- rownames(degs.res)
# [2] # take care of Gene.Symbol maybe Symbol!!! depends on Annotation
ids_of_interest <-  filter(full_results,Gene.Symbol %in% my_genes) %>% 
  pull(ID)
# [3] # take care of Gene.Symbol maybe Symbol!!! depends on Annotation
gene_names <-  filter(full_results,Gene.Symbol %in% my_genes) %>% 
  pull(Gene.Symbol)
# [4]
ex <- exprs(gse2)[,1:6]
gene_matrix <- ex[ids_of_interest,]
# [5]
ctrl.indecies= c(1:3) #### or  seq(from=9,to=dim(exp)[2])
case.indecies=c(4:6)
ctrl.vector=rep("Ctrl", length(ctrl.indecies))
case.vector=rep("Case", length(case.indecies))
Sample=c( ctrl.vector, case.vector)
annotation=as.data.frame(Sample)
rownames(annotation)= colnames(ex)
annotation

# Specify colors
Sample = c("lightgreen", "navy")
names(Sample) = c("Ctrl", "Case")
ann_colors = list(Sample = Sample)

pheatmap(gene_matrix,labels_row = ids_of_interest, annotation = annotation, annotation_colors = ann_colors, scale="row",fontsize_row = 4,fontsize_col = 5)
pheatmap(gene_matrix,labels_row = gene_names, annotation = annotation, annotation_colors = ann_colors, scale="row",fontsize_row = 4,fontsize_col = 5)

heatmap.2(gene_matrix,labels_row = gene_names, cexCol=0.7,col = rev(redblue(26)), scale = "row")

#=======================[Questions]=======================================

## Get the results for particular gene of interest
filter(full_results1, Gene.Symbol == "SMOX")
## Get results for genes with TP53 in the name
filter(full_results1, grepl("TP53", Gene.Symbol))




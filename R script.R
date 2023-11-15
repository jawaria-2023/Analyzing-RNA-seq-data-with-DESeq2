library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)

# Read the data into R
seqdata <- read.delim("data/GSE60450_Lactation-GenewiseCounts.txt", stringsAsFactors = FALSE)
# Read the sample information into R
sampleinfo <- read.delim("data/SampleInfo.txt", stringsAsFactors = TRUE)

#############
head(seqdata)
dim(seqdata)

# Format the data
Let’s create a new data object, countdata, that contains only the counts for the 12 samples.

# Remove first two columns from seqdata
countdata <- seqdata[,-(1:2)]
# Look at the output
head(countdata)
# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]
#Take a look at the output
head(countdata)
#Now take a look at the column names
colnames(countdata)
# using substr, you extract the characters starting at position 1 and stopping at position 7 of the colnames
colnames(countdata) <- substr(colnames(countdata),start=1,stop=7)
#Take a look at the output.
head(countdata)

table(colnames(countdata)==sampleinfo$SampleName)

##########Convert counts to DGEList object
y <- DGEList(countdata)
# have a look at y
y
# See what slots are stored in y
names(y)
# Library size information is stored in the samples slot
y$samples

We can also store the groups for the samples in the DGEList object.

group <- paste(sampleinfo$CellType,sampleinfo$Status,sep=".")
# Take a look
group

# Convert to factor
group <- factor(group)
# Take another look.
group

# Add the group information into the DGEList
y$samples$group <- group
y$samples

#####Adding annotation
columns(org.Mm.eg.db)

ann <- select(org.Mm.eg.db,keys=rownames(y$counts),columns=c("ENTREZID","SYMBOL","GENENAME"))
'select()' returned 1:1 mapping between keys and columns
# Have a look at the annotation
head(ann)

##Let’s double check that the ENTREZID column matches exactly to our y$counts rownames.
table(ann$ENTREZID==rownames(y$counts))

###Filtering lowly expressed genes
# Obtain CPMs
myCPM <- cpm(countdata)
# Have a look at the output
head(myCPM)
# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)

# Summary of how many TRUEs there are in each row
# There are 11433 genes that have TRUEs in all 12 samples.
table(rowSums(thresh))
# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
# We will look at the first sample
plot(myCPM[,1],countdata[,1])


# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5)

#Library size and distribution plots
First, we can check how many reads we have for each sample in the y.
y$samples$lib.size
We can also plot the library sizes as a barplot to see whether there are any major discrepancies between the samples more easily.

# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(y$samples$lib.size,names=colnames(y),las=2)
# Add a title to the plot
title("Barplot of library sizes")

# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

#######Multidimensional scaling plots
plotMDS(y)

# We specify the option to let us plot two plots side-by-sde
par(mfrow=c(1,2))
# Let's set up colour schemes for CellType
# How many cell types and in what order are they stored?
levels(sampleinfo$CellType)

## Let's choose purple for basal and orange for luminal
col.cell <- c("purple","orange")[sampleinfo$CellType]
data.frame(sampleinfo$CellType,col.cell)

# Redo the MDS with cell type colouring
plotMDS(y,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend("topleft",fill=c("purple","orange"),legend=levels(sampleinfo$CellType))
# Add a title
title("Cell type")

# Similarly for status
levels(sampleinfo$Status)
plotMDS(y,col=col.status)
legend("topleft",fill=c("blue","red","black"),legend=levels(sampleinfo$Status),cex=0.8)
title("Status")

#### Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[sampleinfo$CellType]

# Plot the heatmap
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")

mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
aheatmap(highly_variable_lcpm,col=rev(morecols(50)),main="Top 500 most variable genes across samples",annCol=sampleinfo[, 3:4],labCol=group, scale="row")


































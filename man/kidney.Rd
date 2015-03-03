\name{kidney}
\alias{kidney}
\docType{data}
\title{
Kidney Renal Clear Cell Carcinoma [KIRC] RNA-Seq data
}
\description{
The KIRC RNA-seq dataset from The Cancer Genome Atlas containing 20531 genes and 72 paired columns of data with rows corresponding to genes and columns corresponding to replicates; replic vector specifies replicates and treatment vector specifies non-tumor and tumor group samples respectively within replicate.

The maximum possible sample size that can be simulated with this dataset is 36 replicates in each of two treatment groups. However, it is recommended to use a sample size of 20 or lower for simulation studies to ensure each simulated dataset is sufficiently different from the other simulated datasets.

Disclaimer: 
The version of the KIRC dataset provided is:
unc.edu_KIRC.IlluminaHiSeq_RNASeqV2.Level_3.1.5.0. 

The Cancer Genome Atlas updates its datasets periodically. The latest version of the KIRC dataset can always be downloaded from:
https://tcga-data.nci.nih.gov/tcga/. 

Please appropriately reference the source below when using this dataset. The source code used to assemble this dataset is provided below.
}
\usage{data(kidney)}
\format{
List containing:
\itemize{
\item counts: matrix of RNA-seq data for 20531 sampled genes and 72 paired columns from individuals with Kidney Renal Clear Cell Carcinoma.
\item replic: vector detailing which column in counts matrix belongs to each individual.
\item treatment: vector detailing whether each column in counts matrix is a non-tumor or tumor sample.
}
}
\source{
https://tcga-data.nci.nih.gov/tcga/

The Cancer Genome Atlas Research Network (2013). Comprehensive molecular characterization of clear cell renal cell carcinoma. Nature, 499(7456), 43-49.
}
\examples{
data(kidney)

\dontrun{
### Source code used to assemble KIRC dataset

### htmlToText function used to scrape barcode data from uuid
htmlToText <- function(input, ...) {
  ###---PACKAGES ---###
  require(RCurl)
  require(XML) 
  
  ###--- LOCAL FUNCTIONS ---###
  # Determine how to grab html for a single input element
  evaluate_input <- function(input) {    
    # if input is a .html file
    if(file.exists(input)) {
      char.vec <- readLines(input, warn = FALSE)
      return(paste(char.vec, collapse = ""))
    }
    
    # if input is html text
    if(grepl("</html>", input, fixed = TRUE)) return(input)
    
    # if input is a URL, probably should use a regex here instead?
    if(!grepl(" ", input)) {
      # downolad SSL certificate in case of https problem
      if(!file.exists("cacert.perm")) {
        download.file(url = "http://curl.haxx.se/ca/cacert.pem", destfile = "cacert.perm")
      }
      return(getURL(input, followlocation = TRUE, cainfo = "cacert.perm"))
    }
    
    # return NULL if none of the conditions above apply
    return(NULL)
  }
  
  # convert HTML to plain text
  convert_html_to_text <- function(html) {
    doc <- htmlParse(html, asText = TRUE)
    text <- xpathSApply(doc, paste0("//text()",
                   "[not(ancestor::script)][not(ancestor::style)]",
                   "[not(ancestor::noscript)][not(ancestor::form)]"), xmlValue)
    return(text)
  }
  
  # format text vector into one character string
  collapse_text <- function(txt) {
    return(paste(txt, collapse = " "))
  }
  
  ###--- MAIN ---###
  # STEP 1: Evaluate input
  html.list <- lapply(input, evaluate_input)
  
  # STEP 2: Extract text from HTML
  text.list <- lapply(html.list, convert_html_to_text)
  
  # STEP 3: Return text
  text.vector <- sapply(text.list, collapse_text)
  return(text.vector)
}

### Specify path name for folder containing raw counts for each sample
mainDir <- getwd()
folder.path <- "unc.edu_KIRC.IlluminaHiSeq_RNASeqV2.Level_3.1.5.0"

### Determine list of files containing summarized raw counts
file.list <- dir(file.path(mainDir, folder.path))
keep <- grepl("genes.results", file.list)
file.list <- file.list[keep]

### Create summarized count matrix.
### Get n.row and n.col for summarized count matrix number of genes in first
### sample and number of total samples from file.list

file.temp <- file.path(mainDir, folder.path, file.list[1])
n.row <- nrow(read.table(file = file.temp, header = TRUE))
n.col <- length(file.list)

### initialize counts matrix
counts <- matrix(NA, nrow = n.row, ncol = n.col)

### read in raw read counts from file.list
for(i in 1:n.col){
  file.temp <- file.path(mainDir, folder.path, file.list[i])
  counts[, i] <- read.table(file.temp, header = TRUE)$raw_count
}

### Data was summarized using RSEM software which produces non_integer 
### counts for ambiguous reads. Counts are rounded as a preprocessing
### step.
counts <- round(counts)

### Get uuid's for each sample
uuid <- substr(file.list, start = 9, stop = 44)

### Create urls from uuid list
urls <- paste(rep("https://tcga-data.nci.nih.gov/uuid/uuidws/mapping/xml/uuid/",
                  length(uuid)), uuid, sep = "")

### Scrape barcodes from urls
l <- length(urls)
barcodes <- vector("character", l)
for(i in 1:l){
  barcodes[i] <- htmlToText(urls[i])
}


barcodes <- substr(barcodes, start = 1, stop = 28)

### Get metadata on which samples were taken from each individual, 
### tumor type of sample, etc. from barcodes for each sample
metadata <- data.frame(barcodes)

### Study Participant
metadata$participant <- substr(barcodes, start = 9, stop = 12)

### Sample type code. See:
### https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm?codeTable=Sample\%20type
### for full list of codes and details on TCGA barcodes.
### 01: Primary Solid Tumor
### 02: Recurrent Solid Tumor
### 05: Additional New Primary
### 06: Metastatic Tumor
### 11: Solid Tissue Normal
metadata$type <- substr(barcodes, start = 14, stop = 15)

### Only keep Primary Solid Tumor and Solid Tissue Normal
keep.metadata <- metadata$type == "01" | metadata$type == "11"
metadata <- metadata[keep.metadata, ]
counts <- counts[, keep.metadata]

### Code from 01 to Tumor and 11 to Non-Tumor for easy identifiability
metadata$tumor <- "Non-Tumor"
metadata$tumor[metadata$type == "01"] <- "Tumor"

### tag participant, type, and tumor as factors
metadata$participant <- as.factor(metadata$participant)
metadata$type <-as.factor(metadata$type)
metadata$tumor <- as.factor(metadata$tumor)

### Sort and subset down to paired data
sorting <- 
  SortData(counts, treatment = metadata$tumor, 
           replic = metadata$participant, sort.method = "paired")$sorting

counts <- counts[, sorting]
metadata <- metadata[sorting, ]
metadata$participant <- factor(metadata$participant)
metadata$barcodes <- factor(metadata$barcodes)

kidney <- vector("list", 3)
kidney[[1]] <- counts
kidney[[2]] <- metadata$participant
kidney[[3]] <- metadata$tumor
names(kidney) <- c("counts", "replic", "treatment")

###Save file
save(kidney, "kidney.rda")
}
}
\keyword{datasets}

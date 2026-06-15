#install.packages('e1071')
#install.packages('parallel')

#if (!requireNamespace("BiocManager", quietly = TRUE))
   # install.packages("BiocManager")
#BiocManager::install("preprocessCore", version = "3.8")


setwd("G:\\CIBERSORT")
source("CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=100, QN=TRUE)
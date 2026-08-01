if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(version="3.22", c("DESeq2","edgeR","fgsea"), ask=FALSE)
install.packages(c("data.table","ggplot2","cowplot","readr","dplyr"))

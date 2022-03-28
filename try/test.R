setwd('D:\\zhanglab\\DMS_opt-master')
library(Biostrings)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library (stringdist)
library(tibble)
library(tidyr)
library("magrittr")

pKas <- read_csv("data/pKas.csv",col_names = FALSE)%>%column_to_rownames(var = "X1")
Eisenberg <<- read_csv("data/Eisenberg.csv", col_names = FALSE) %>%column_to_rownames(var = "X1")

#net_charge +test of one AA sequence
AA_counts <- sapply(rownames(pKas), 
                    function(x) str_count("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEW", x))
print(AA_counts)

dod <- 1/(10^(pKas - pH) + 1)
dod[, 1] <- -1*dod[, 1]
dod[, 2] <- 1 - dod[, 2]
dod[c(2,3,4,20), 3] <- -1*dod[c(2,3,4,20), 3]
dod[c(7,9,15), 3] <- 1 - dod[c(7,9,15), 3]

net_charge <- sum(AA_counts*rowSums(dod, na.rm = TRUE))
print(net_charge)

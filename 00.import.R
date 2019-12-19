library(Seurat)
library(dplyr)
library(pheatmap)
library(monocle)
library(velocyto.R)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(fgsea)
library(biomaRt)
library(knitr)
library(xtable)
library(gameofthrones)
library(cowplot)

setwd(".")

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
col=ggplotColours(n = 9)

eil <- read.delim("ressources/EIL.txt", sep="\t", header=FALSE)
colnames(eil) <- c("name","stade")
genes.early <- eil[which(eil$stade == "E"),]$name
genes.intermediate <- eil[which(eil$stade == "I"),]$name
genes.late <- eil[which(eil$stade == "L"),]$name
genes.vv <- eil$name

cc.genes <- readLines(con = "ressources/regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

cell_type_color <- c("SCRG1" = "#466cb9",
                     "COL1A2" = "#FCCC0A",
                     "COL1A2 cycling" = "#ffc107",
                     "SCRG1 cycling" = "#03a9f4",
                     "KRT14" = "#c997e7",
                     "SPP1" = "#A31E22",
                     "MET" = "#aea9ce",
                     "MDH1" = "#2da9d2",
                     "CA2" = "#a7d6a7",
                     "late" = "#ff6600",
                     "sain" = "grey",
                     "-virus" = "#9b9b98",
                     "s2.sain" = "grey",
                     "s3.sain" = "grey",
                     "infected" = "#ff9000",
                     "+virus" = "#F3766E",
                     "s2.infected" = "#ff6600",
                     "s3.infected" = "#ff6600",
                     "bystander" = "#006600",
                     "naive" = "#BDBDBD")

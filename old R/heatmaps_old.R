## ---------------------------
##
## Script Name: OP Analysis heatmaps
##
## Purpose of Script: Summarize results of OP analysis in heatmap graphics.
##
## Author: Zachary Rowson
##
## Date Created: 2022-02-17
##
## Email: Rowson.Zachary@epa.gov
##
## Working Directory: "/ccte/home2/zrowson/Desktop/Stephanie Projects/OP Analysis"
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------


library(data.table)
library(pheatmap)
library(viridis)


# heatmaps - Developmental Exposures ---------------------------------------------


# Load tcplfits for developmental Exposure
load("~/Desktop/Stephanie Projects/OP Analysis/pipelining/pipelined data/Developmental Exposure/Padilla_OP_Dev_tcplfits.rda")


##
## Directional hitcall heatmap ##
##

## gather directional hitcalls for each chemical endpoint pair
hits.dir <- lapply(tcplfits.dev_n, function(chm) lapply(chm, function(endp) endp[["summary"]]$hitcall*sign(endp[["summary"]]$top)))
temp <- lapply(hits.dir, function(chm) unlist(chm))
hits.dir <- do.call("rbind", temp); rm(temp)

## Order columns
ac <- c("AUC_L","avgS_L", "avgA_L", "avgJ_L","hbt_L",
        "strtlA","strtlAavg", "strtlF",
        "AUC_D", "avgS_D","avgA_D", "avgJ_D", "hbt_D",
        "AUC_T","avgS_T","AUC_r")
hits.dir <- hits.dir[, ac]

## Create annotations
ann.df <- data.frame(Phase = c(rep("Light",5), rep("Transition",3), rep("Dark",5), rep("Other",3)))
row.names(ann.df) <- ac
ann_colors <- c(Light = "white", Transition = "grey", Dark = "Black", Other = "red")

## Specify heatmap key breaks
breaks <- seq(from=-1, to=1, by=0.25)
colors <- rev(RColorBrewer::brewer.pal(length(breaks)-1, "RdBu"))
dir.hitsmap.dev <- pheatmap::pheatmap(hits.dir,
                                  main = "Directional Hitcall Heatmap for Padilla OP Chemicals \n Developmental Exposures",
                                  clustering_method = "ward.D2",
                                  cluster_rows = T, cluster_cols = F,
                                  cutree_rows = 7,
                                  breaks = breaks,
                                  color = colors,
                                  annotation_col = ann.df, annotation_colors = list(Phase=ann_colors))

## extract chemical order after clustering
order <- dir.hitsmap.dev$tree_row$order



# Gather BMC's
bmd <- lapply(tcplfits.dev_n, function(chm) unlist(lapply(chm, function(fit) {
  hitcall <- fit[["summary"]]$hitcall
  bmd <- fit[["summary"]]$bmd
  if (!(hitcall>0.8 & !is.na(bmd) & bmd != 0)) {
    bmd <- 10000
  }
  bmd
})))
bmd <- do.call("rbind", bmd)
bmd <- log10(bmd[, ac])

## fit heat map with BMD values by chemical-endpoint pair
labs <- c(-2:3,paste0("4 log(","\U03BC","M)"))
bmdmap.dev <- pheatmap::pheatmap(bmd[order,],
                             main = "BMC Heatmap for Padilla OP Chemicals \n Developmental Exposures",
                             cluster_rows = F, cluster_cols = F,
                             color = rev(viridis(8)),
                             legend_breaks = -2:4,
                             legend_labels = labs,
                             annotation_col = ann.df, annotation_colors = list(Phase=ann_colors))

##
## Save heatmap Graphics ##
##

setwd("~/Desktop/Stephanie Projects/OP Analysis/graphics/heatmaps")

png("Padilla_OP_Dev_Directional hitsmap.png",
    width = 720,
    height = 720)
print(dir.hitsmap.dev)
dev.off()

png("Padilla_OP_Dev_BMC heatmap.png",
    width = 720,
    height = 720)
print(bmdmap.dev)
dev.off()


rm(list=ls())

# heatmaps - Acute Exposures ---------------------------------------------


# Load tcplfits for acute Exposure
load("~/Desktop/Stephanie Projects/OP Analysis/pipelining/pipelined data/Acute Exposure/Padilla_OP_Acute_tcplfits.rda")


##
## Directional hitcall heatmap ##
##

## gather directional hitcalls for each chemical endpoint pair
hits.dir <- lapply(tcplfits.act_n, function(chm) lapply(chm, function(endp) endp[["summary"]]$hitcall*sign(endp[["summary"]]$top)))
temp <- lapply(hits.dir, function(chm) unlist(chm))
hits.dir <- do.call("rbind", temp); rm(temp)

## Order columns
ac <- c("AUC_L","avgS_L", "avgA_L", "avgJ_L","hbt_L",
        "strtlA","strtlAavg", "strtlF",
        "AUC_D", "avgS_D","avgA_D", "avgJ_D", "hbt_D",
        "AUC_T","avgS_T","AUC_r")
hits.dir <- hits.dir[, ac]

## Create annotations
ann.df <- data.frame(Phase = c(rep("Light",5), rep("Transition",3), rep("Dark",5), rep("Other",3)))
row.names(ann.df) <- ac
ann_colors <- c(Light = "white", Transition = "grey", Dark = "Black", Other = "red")

## Specify heatmap key breaks
breaks <- seq(from=-1, to=1, by=0.25)
colors <- rev(RColorBrewer::brewer.pal(length(breaks)-1, "RdBu"))
dir.hitsmap.act <- pheatmap::pheatmap(hits.dir,
                                      main = "Directional Hitcall Heatmap for Padilla OP Chemicals \n Acute Exposures",
                                      clustering_method = "ward.D2",
                                      cluster_rows = T, cluster_cols = F,
                                      cutree_rows = 7,
                                      breaks = breaks,
                                      color = colors,
                                      annotation_col = ann.df, annotation_colors = list(Phase=ann_colors))

## extract chemical order after clustering
order <- dir.hitsmap.act$tree_row$order


##
## BMD heat map ##
##

## Gather BMC's
bmd <- lapply(tcplfits.act_n, function(chm) unlist(lapply(chm, function(fit) {
  hitcall <- fit[["summary"]]$hitcall
  bmd <- fit[["summary"]]$bmd
  if (!(hitcall>0.8 & !is.na(bmd) & bmd != 0)) {
    bmd <- 10000
  }
  bmd
})))
bmd <- do.call("rbind", bmd)
bmd <- log10(bmd[, ac])

## fit heat map with BMD values by chemical-endpoint pair
labs <- c(-2:3,paste0("4 log(","\U03BC","M)"))
bmdmap.act <- pheatmap::pheatmap(bmd[order,],
                                 main = "BMC Heatmap for Padilla OP Chemicals \n Acute Exposures",
                                 cluster_rows = F, cluster_cols = F,
                                 color = rev(viridis(8)),
                                 legend_breaks = -2:4,
                                 legend_labels = labs,
                                 annotation_col = ann.df, annotation_colors = list(Phase=ann_colors))

##
## Save heatmap Graphics ##
##

setwd("~/Desktop/Stephanie Projects/OP Analysis/graphics/heatmaps")

png("Padilla_OP_Acute_Directional hitsmap.png",
    width = 720,
    height = 720)
print(dir.hitsmap.act)
dev.off()

png("Padilla_OP_Acute_BMC heatmap.png",
    width = 720,
    height = 720)
print(bmdmap.act)
dev.off()


rm(list=ls())

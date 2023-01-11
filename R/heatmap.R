## ---------------------------
##
## Script Name: OP Analysis Heatmaps
##
## Purpose of Script:
##
## Author: Zachary Rowson
##
## Date Created: 2022-09-16
##
## Email: Rowson.Zachary@epa.gov
##
## Working Directory: "/ccte/home2/zrowson/Desktop/StephanieProjects/OPAnalysis/pipelining"
##
## ---------------------------
##
## Notes: New version of heatmaps for OP Analysis
##
##
## ---------------------------


library(data.table)
library(ComplexHeatmap)
library(viridis)

load("~/Desktop/StephanieProjects/OPAnalysis/pipelining/pipelined data/Developmental Exposure/Padilla_OP_Dev_tcplfits.rda")
load("~/Desktop/StephanieProjects/OPAnalysis/pipelining/pipelined data/Acute Exposure/Padilla_OP_Acute_tcplfits.rda")

setwd("~/Desktop/StephanieProjects/OPAnalysis/graphics/heatmaps")

# Developmental Heatmap ---------------------------------------------------


# Gather BMD data
ac <- c("avgS_L", "avgA_L", "avgJ_L","hbt_L",
        "strtlA","strtlAavg", "strtlF",
        "avgS_D","avgA_D", "avgJ_D", "hbt_D",
        "avgS_T","AUC_r")

## Gather BMC's
bmd.dev <- lapply(tcplfits.dev_n, function(chm) unlist(lapply(chm, function(fit) {
  hitcall <- fit[["summary"]]$hitcall
  bmd <- fit[["summary"]]$bmd
  if (!(hitcall>0.8 & !is.na(bmd) & bmd != 0)) {
    bmd <- 10000
  }
  bmd
})))
bmd.dev <- do.call("rbind", bmd.dev)
bmd.dev <- log10(bmd.dev[, ac])
# bmd <- bmd[,-c(1,9,14)]

# Identify active chemicals
actives.dev <- rownames(bmd.dev)[apply(bmd.dev, 1, function(row) any(row < 4))]

# Isolate data of interest
to.fit <- bmd.dev

# Layer matrix
# Create a matrix specifying if an up, down, up and down, or down and up arrow should be printed in cells.
to.fit.dev <- lapply(tcplfits.dev_n, function(chm) chm[-c(1,9,14)])
layer.mat.dev <- lapply(to.fit.dev, function(chm) unlist(lapply(chm, function(fit) {
  hitcall <- fit[["summary"]]$hitcall
  fit_method <- fit[["summary"]]$fit_method
  dir <- sign(fit[["summary"]]$top)
  if (hitcall>0.8 & fit_method=="gnls" & dir==1) {
    layer <- 0
  } else if (hitcall>0.8 & fit_method=="gnls" & dir==-1) {
    layer <- 1
  } else if (hitcall>0.8 & dir==1) {
    layer <- 2
  } else if (hitcall>0.8 & dir==-1) {
    layer <- 3
  } else layer <- 4
})))

layer.mat.dev <- do.call("rbind", layer.mat.dev)
layer.mat.dev1 <- layer.mat.dev[row.names(layer.mat.dev) %in% actives.dev,]

# Labels
# Column labels
col_labels <- c(expression("Average Speed in Light"^1), "Average Acceleration in Light", "Average Jerk in Light", expression("Habituation in Light"^3),
                "Startle Acceleration", "Startle Relative to Avg. Speed in Light", "Startle Fold-Change",
                expression("Average Speed in Dark"^1), "Average Acceleration in Dark", "Average Jerk in Dark", expression("Habituation in Dark"^3),
                "Average Speed in Both Phases", expression("AUC in Dark / AUC in Light Ratio"^2)) # Superscripts notate references in poster
# Legends

# Custom heat legend.
f2 = circlize::colorRamp2(seq(min(bmd.dev), max(bmd.dev), length = 8), rev(viridis(8)), space = "sRGB")
heat_lgd = ComplexHeatmap::Legend(col_fun = f2,
                  title = paste0("BMC log(","\U03BC","M)"),
                  title_position = "lefttop",
                  legend_width = unit(4,"cm"),
                  direction = "horizontal")
# Column annotation by phase legend.
ann_lgd = ComplexHeatmap::Legend(labels = c("Light","Transition","Dark","Light+Dark"),
                 title = "Phase",
                 title_position = "leftcenter",
                 labels_gp = gpar(fontsize=8),
                 title_gp = gpar(fontsize=8),
                 legend_gp = grid::gpar(fill = c("white","grey","black","red")),
                 border = TRUE,
                 nrow = 1,
                 column_gap = unit(5, 'mm'))
# Legend for cell signal arrows.
dir_lgd = ComplexHeatmap::Legend(labels = c(paste("\U2191","Gain"),
                            paste("\U2193","Loss"),
                            paste("\U21C5","GainLoss"),
                            paste("\U21F5","LossGain")),
                 title = "Signal Direction",
                 labels_gp = gpar(fontsize=8),
                 title_gp = gpar(fontsize=8),
                 title_position = "leftcenter",
                 nrow = 1,
                 column_gap = unit(0, 'mm')) # Will produce warnings, don't worry.
lgd_list <- packLegend(ann_lgd, dir_lgd)

# Create column annotation indicating the phase of the LMR that is described.
column_ha <- ComplexHeatmap::HeatmapAnnotation(Phase = factor(c(rep("Light",4), rep("Transition",3),rep("Dark",4),rep("Light+Dark",2)),
                                                              levels=c("Light","Transition","Dark","Light+Dark")),
                                               border = TRUE,
                                               simple_anno_size = unit(0.25, 'cm'),
                                               col = list(Phase=c("Light"="white",
                                                                  "Transition"="grey",
                                                                  "Dark"="black",
                                                                  "Light+Dark"="red")),
                                               annotation_legend_param = list(nrow = 1),
                                               show_annotation_name = FALSE,
                                               show_legend = FALSE)

# Create function to add arrows indicating signal direction in cells.
cell_fun.dev <- function(j, i, x, y, width, height, fill) {
  if (layer.mat.dev1[i,j] == 0) {
    grid.text("\U21C5", x, y, gp=gpar(fontsize=14))
  } else if (layer.mat.dev1[i,j] == 1) {
    grid.text("\U21F5", x, y, gp=gpar(fontsize=14))
  } else if (layer.mat.dev1[i,j] == 2) {
    grid.text("\U2191", x, y, gp=gpar(fontsize=14))
  } else if (layer.mat.dev1[i,j] == 3) {
    grid.text("\U2193", x, y, gp=gpar(fontsize=14))
  }
}

# Heatmap
# Create main heat map.
htlist.dev <- ComplexHeatmap::Heatmap(bmd.dev[row.names(bmd.dev) %in% actives.dev,],

                                  # Specify some parameters for heat legend.
                                  name = paste0("BMC log(","\U03BC","M)"),
                                  col = f2,
                                  border_gp = grid::gpar(col="black",lwd=1),
                                  rect_gp=grid::gpar(col="grey"),
                                  show_heatmap_legend = TRUE,
                                  heatmap_legend_param = list(legend_height = unit(3.5,"cm"),
                                                              direction = "vertical",
                                                              title_gp = grid::gpar(fontsize=10),
                                                              labels_gp = grid::gpar(fontsize=10)),

                                  # Heatmap width and height
                                  # width = unit(16,"in"),
                                  # height = unit(18,"in"),

                                  # Column label parameters
                                  column_labels = col_labels,
                                  column_names_rot = 45,


                                  # Split columns by phase.
                                  column_split = factor(c(rep("Light",4), rep("Transition",3),rep("Dark",4),rep("Light+Dark",2)),
                                                        levels=c("Light","Transition","Dark","Light+Dark")),

                                  # Specify some column parameters.
                                  column_title = NULL,
                                  top_annotation = column_ha,

                                  # Add signal direction arrows.
                                  cell_fun = cell_fun.dev,

                                  # Specify some parameters row dendrogram aesthetics and row labels.
                                  row_title_side = "right",
                                  row_title_rot = 0,
                                  row_split = 5,
                                  row_dend_side = "right",
                                  row_names_side = "left",

                                  # Clustering parameters.
                                  cluster_columns = FALSE,
                                  cluster_rows = TRUE,
                                  clustering_distance_rows = "pearson",

                                  # Font sizes
                                  row_names_gp = grid::gpar(fontsize=8),
                                  column_names_gp = grid::gpar(fontsize=10),
                                  row_title_gp = grid::gpar(fontsize=10),
                                  column_title_gp = grid::gpar(fontsize=10)
)

# Save heatmap
png(filename = "Padilla_OP_Dev_BMC_heatmap.png", width = 22, height = 26, unit = "cm", res=300)
draw(htlist.dev, merge_legend = FALSE, annotation_legend_list = lgd_list,
     annotation_legend_side = "top", align_annotation_legend = "heatmap_center")
dev.off()


# Acute Heatmap ---------------------------------------------------


# Gather BMD data
ac <- c("avgS_L", "avgA_L", "avgJ_L","hbt_L",
        "strtlA","strtlAavg", "strtlF",
        "avgS_D","avgA_D", "avgJ_D", "hbt_D",
        "avgS_T","AUC_r")

## Gather BMC's
bmd.act <- lapply(tcplfits.act_n, function(chm) unlist(lapply(chm, function(fit) {
  hitcall <- fit[["summary"]]$hitcall
  bmd <- fit[["summary"]]$bmd
  if (!(hitcall>0.8 & !is.na(bmd) & bmd != 0)) {
    bmd <- 10000
  }
  bmd
})))
bmd.act <- do.call("rbind", bmd.act)
bmd.act <- log10(bmd.act[, ac])
# bmd <- bmd[,-c(1,9,14)]

# Identify active chemicals
actives.act <- rownames(bmd.act)[apply(bmd.act, 1, function(row) any(row < 4))]

# Isolate data of interest
to.fit <- bmd.act

# Layer matrix
# Create a matrix specifying if an up, down, up and down, or down and up arrow should be printed in cells.
to.fit.act <- lapply(tcplfits.act_n, function(chm) chm[-c(1,9,14)])
layer.mat.act <- lapply(to.fit.act, function(chm) unlist(lapply(chm, function(fit) {
  hitcall <- fit[["summary"]]$hitcall
  fit_method <- fit[["summary"]]$fit_method
  dir <- sign(fit[["summary"]]$top)
  if (hitcall>0.8 & fit_method=="gnls" & dir==1) {
    layer <- 0
  } else if (hitcall>0.8 & fit_method=="gnls" & dir==-1) {
    layer <- 1
  } else if (hitcall>0.8 & dir==1) {
    layer <- 2
  } else if (hitcall>0.8 & dir==-1) {
    layer <- 3
  } else layer <- 4
})))

layer.mat.act <- do.call("rbind", layer.mat.act)
layer.mat.act1 <- layer.mat.act[row.names(layer.mat.act) %in% actives.act,]

# Labels
# Column labels
col_labels <- c(expression("Average Speed in Light"^1), "Average Acceleration in Light", "Average Jerk in Light", expression("Habituation in Light"^3),
                "Startle Acceleration", "Startle Relative to Avg. Speed in Light", "Startle Fold-Change",
                expression("Average Speed in Dark"^1), "Average Acceleration in Dark", "Average Jerk in Dark", expression("Habituation in Dark"^3),
                "Average Speed in Both Phases", expression("AUC in Dark / AUC in Light Ratio"^2)) # Superscripts notate references in poster
# Legends

# Custom heat legend.
f2 = circlize::colorRamp2(seq(min(bmd.act), max(bmd.act), length = 8), rev(viridis(8)), space = "sRGB")
heat_lgd = ComplexHeatmap::Legend(col_fun = f2,
                                  title = paste0("BMC log(","\U03BC","M)"),
                                  title_position = "lefttop",
                                  legend_width = unit(4,"cm"),
                                  direction = "horizontal")
# Column annotation by phase legend.
ann_lgd = ComplexHeatmap::Legend(labels = c("Light","Transition","Dark","Light+Dark"),
                                 title = "Phase",
                                 title_position = "leftcenter",
                                 labels_gp = gpar(fontsize=8),
                                 title_gp = gpar(fontsize=8),
                                 legend_gp = grid::gpar(fill = c("white","grey","black","red")),
                                 border = TRUE,
                                 nrow = 1,
                                 column_gap = unit(5, 'mm'))
# Legend for cell signal arrows.
dir_lgd = ComplexHeatmap::Legend(labels = c(paste("\U2191","Gain"),
                                            paste("\U2193","Loss"),
                                            paste("\U21C5","GainLoss"),
                                            paste("\U21F5","LossGain")),
                                 title = "Signal Direction",
                                 labels_gp = gpar(fontsize=8),
                                 title_gp = gpar(fontsize=8),
                                 title_position = "leftcenter",
                                 nrow = 1,
                                 column_gap = unit(0, 'mm')) # Will produce warnings, don't worry.
lgd_list <- packLegend(ann_lgd, dir_lgd)

# Create column annotation indicating the phase of the LMR that is described.
column_ha <- ComplexHeatmap::HeatmapAnnotation(Phase = factor(c(rep("Light",4), rep("Transition",3),rep("Dark",4),rep("Light+Dark",2)),
                                                              levels=c("Light","Transition","Dark","Light+Dark")),
                                               border = TRUE,
                                               simple_anno_size = unit(0.25, 'cm'),
                                               col = list(Phase=c("Light"="white",
                                                                  "Transition"="grey",
                                                                  "Dark"="black",
                                                                  "Light+Dark"="red")),
                                               annotation_legend_param = list(nrow = 1),
                                               show_annotation_name = FALSE,
                                               show_legend = FALSE)

# Create function to add arrows indicating signal direction in cells.
cell_fun.act <- function(j, i, x, y, width, height, fill) {
  if (layer.mat.act1[i,j] == 0) {
    grid.text("\U21C5", x, y, gp=gpar(fontsize=14))
  } else if (layer.mat.act1[i,j] == 1) {
    grid.text("\U21F5", x, y, gp=gpar(fontsize=14))
  } else if (layer.mat.act1[i,j] == 2) {
    grid.text("\U2191", x, y, gp=gpar(fontsize=14))
  } else if (layer.mat.act1[i,j] == 3) {
    grid.text("\U2193", x, y, gp=gpar(fontsize=14))
  }
}

# Heatmap
# Create main heat map.
htlist.act <- ComplexHeatmap::Heatmap(bmd.act[row.names(bmd.act) %in% actives.act,],

                                      # Specify some parameters for heat legend.
                                      name = paste0("BMC log(","\U03BC","M)"),
                                      col = f2,
                                      border_gp = grid::gpar(col="black",lwd=1),
                                      rect_gp=grid::gpar(col="grey"),
                                      show_heatmap_legend = TRUE,
                                      heatmap_legend_param = list(legend_height = unit(3.5,"cm"),
                                                                  direction = "vertical",
                                                                  title_gp = grid::gpar(fontsize=10),
                                                                  labels_gp = grid::gpar(fontsize=10)),

                                      # Heatmap width and height
                                      # width = unit(16,"in"),
                                      # height = unit(18,"in"),

                                      # Column label parameters
                                      column_labels = col_labels,
                                      column_names_rot = 45,


                                      # Split columns by phase.
                                      column_split = factor(c(rep("Light",4), rep("Transition",3),rep("Dark",4),rep("Light+Dark",2)),
                                                            levels=c("Light","Transition","Dark","Light+Dark")),

                                      # Specify some column parameters.
                                      column_title = NULL,
                                      top_annotation = column_ha,

                                      # Add signal direction arrows.
                                      cell_fun = cell_fun.act,

                                      # Specify some parameters row dendrogram aesthetics and row labels.
                                      row_title_side = "right",
                                      row_title_rot = 0,
                                      row_split = 5,
                                      row_dend_side = "right",
                                      row_names_side = "left",

                                      # Clustering parameters.
                                      cluster_columns = FALSE,
                                      cluster_rows = TRUE,
                                      clustering_distance_rows = "pearson",

                                      # Font sizes
                                      row_names_gp = grid::gpar(fontsize=8),
                                      column_names_gp = grid::gpar(fontsize=10),
                                      row_title_gp = grid::gpar(fontsize=10),
                                      column_title_gp = grid::gpar(fontsize=10)
)

# Save heatmap
png(filename = "Padilla_OP_Acute_BMC_heatmap.png", width = 22, height = 26, unit = "cm", res=300)
draw(htlist.act, merge_legend = FALSE, annotation_legend_list = lgd_list,
     annotation_legend_side = "top", align_annotation_legend = "heatmap_center")
dev.off()

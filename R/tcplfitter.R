## ---------------------------
##
## Script Name: OP tcplfitter
##
## Purpose of Script: Fit tcplfit2 to OP endpoint data.
##
## Author: Zachary Rowson
##
## Date Created: 2022-02-08
##
## Email: Rowson.Zachary@epa.gov
##
## Working Directory: "/ccte/home2/zrowson/Desktop/Stephanie Projects/OP Analysis/pipelining"
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

library(data.table)
library(magrittr)
library(here)
library(gabi)
library(pheatmap)
library(viridis)

here::i_am("R/tcplfitter.R")


# Developmental Exposure --------------------------------------------------


# Load row data
load("~/Desktop/Stephanie Projects/OP Analysis/pipelining/pipelined data/Developmental Exposure/Padilla_OP_Dev_rows_n.rda")

# Get tcplfits
tcplfits.dev_n <- lapply(rows.dev_n, function (chm) {
  lapply(chm, function(row) {
    gabi::concRespCoreZR(row, do.plot = TRUE)})})


# Save fits
save(tcplfits.dev_n, file = "pipelined data/Developmental Exposure/Padilla_OP_Dev_tcplfits.rda")

# Acute Exposure ----------------------------------------------------------


# Load row data
load("~/Desktop/Stephanie Projects/OP Analysis/pipelining/pipelined data/Acute Exposure/Padilla_OP_Acute_rows_n.rda")

# Get tcplfits
tcplfits.act_n <- lapply(rows.act_n, function (chm) {
  lapply(chm, function(row) {
    gabi::concRespCoreZR(row, do.plot = TRUE)})})

# Save fits
save(tcplfits.act_n, file = "pipelined data/Acute Exposure/Padilla_OP_Acute_tcplfits.rda")

## ---------------------------
##
## Script Name: lmr0 Formatting of OP Data
##
## Purpose of Script: Format OP Data Prior to gabi Analysis
##
## Author: Zachary Rowson
##
## Date Created: 2022-02-08
##
## Email: Rowson.Zachary@epa.gov
##
## Working Directory:"/ccte/home2/zrowson/Desktop/Stephanie Projects/OP Analysis/pipelining"
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

library(data.table)
library(magrittr)
library(readxl)
library(here)

here::i_am("R/data upload.R")

# transform data to lmr0 format


# Developmental Exposure Data ---------------------------------------------


# upload data as data.table and copy

## upload raw data
raw.dev <- as.data.table(readxl::read_excel(here("raw data","Developmental Exposure","OP Dev Final Sept 22 2021.xlsx"), skip=1))
raw.dev1 <- data.table::copy(raw.dev)

# identify and change names

## new time period names
t.new <- paste0("vt", 11:50)

## change names
new.names <- names(raw.dev)
new.names[c(1,3,5:ncol(raw.dev))] <- c("wllq","apid","rowcol","cpid","conc",t.new)
data.table::setnames(raw.dev1, names(raw.dev), new.names)

# edit and/or create lmr0 columns

## edit wllq
raw.dev1[!is.na(`Removed due to tracking issues/errors`), wllq := 0][, `Removed due to tracking issues/errors` := NULL]
raw.dev1[!is.na(wllq), wllq := 0][is.na(wllq), wllq := 1]

# separate rowcol into rowi and coli

## isolate rowcol and then match letter row index with a value from 1:8
## save column index as an integer
rowcol <- raw.dev1[, rowcol]
row <- sapply(rowcol, function (pos) substring(pos, 1, 1)) %>%
  match(LETTERS)
col <- sapply(rowcol, function (pos) substring(pos, 2)) %>%
  as.integer()
raw.dev1[, `:=` (rowi = row, coli = col, rowcol = NULL)]

# create wllt column

## identify positive control, vehicle control, and test fish
raw.dev1[, wllt := "t"]
raw.dev1[cpid == "DMSO", wllt := "v"][cpid=="Chlorpyrifos+", wllt := "p"]
raw.dev1[wllt == "p", cpid := "Chlorpyrifos"]

## add source file and assay component information
srcf <- "OP Dev Final Sept 22 2021.xlsx"
acid <- "ZFlmrALD-20-40-40"
raw.dev1[, `:=` (srcf = srcf, acid = acid)]

# perform final touches

## remove rows that are blank
raw.dev2 <- raw.dev1[!is.na(cpid)]

## ensure concentration and behavior data is of numeric class
raw.dev2[, (c("conc", t.new)) := lapply(.SD, as.numeric), .SDcols = c("conc", t.new)]

# create lmr0 formatted table

## reorder columns and declare lmr0 object
lmr0.cols <- c("srcf", "acid", "cpid", "apid", "rowi",
               "coli", "wllt", "wllq", "conc", t.new)
lmr0.dev <- raw.dev2[, ..lmr0.cols]

# Find chemical|concentration groups and control groups with too many dead or malformed fish
lmr0.dev.perc.wllq.0.rmv.conc <- lmr0.dev[wllt=="t", .(perc.wllq.0 = length(which(wllq==0)) / .N), by = .(cpid,apid,conc)][perc.wllq.0 > 0.25, .(cpid,apid,conc)]
lmr0.dev.perc.wllq.0.rmv.apid <- lmr0.dev[wllt=="v", .(perc.wllq.0 = length(which(wllq==0)) / .N), by = .(cpid,apid,conc)][perc.wllq.0 > 0.20, apid]

# Save this information
save(lmr0.dev.perc.wllq.0.rmv.conc, file = here("raw data","Developmental Exposure","Padilla_OP_Dev_ChemConcs Excluded.Rdata"))
save(lmr0.dev.perc.wllq.0.rmv.apid, file = here("raw data","Developmental Exposure","Padilla_OP_Dev_Plates Excluded.Rdata"))

# Remove these chemical|concentration groups and plates
lmr0.dev1 <- lmr0.dev[!lmr0.dev.perc.wllq.0.rmv.conc, on = .(cpid,apid,conc)][!(apid %in% lmr0.dev.perc.wllq.0.rmv.apid)]

# save lmr0 formatted data
rm(lmr0.dev)
lmr0.dev <- data.table::copy(lmr0.dev1)
save(lmr0.dev, file = "raw data/Developmental Exposure/Padilla_OP_Dev_lmr0.Rdata")

# Clear environment prior to acute exposure data formatting
rm(list = ls())

# Acute Exposure Data -----------------------------------------------------


# NOTE Light Exposure is different
# 3-5-5 protocol. Or 6min-10min-10min protocol
path <- "raw data/Acute Exposure/OP Acute Final Sept 22 2021.xlsx"
raw.act <- as.data.table(readxl::read_excel(path))
raw.act1 <- data.table::copy(raw.act)

# identify and change names

## new time period names
t.new <- paste0("vt", 1:13)

## change names
new.names <- names(raw.act)
new.names[c(1,3,5:ncol(raw.act))] <- c("wllq","apid","rowcol","cpid","conc",t.new)
data.table::setnames(raw.act1, names(raw.act), new.names)

# edit and/or create lmr0 columns

## edit wllq
raw.act1[!is.na(`Removed due to tracking issues/errors`), wllq := 0][, `Removed due to tracking issues/errors` := NULL]
raw.act1[!is.na(wllq), wllq := 0][is.na(wllq), wllq := 1]

## remove sudo-egid column for now
raw.act1[, `Link plates tested together` := NULL]

## separate rowcol into rowi and coli

## isolate rowcol and then match letter row index with a value from 1:8
## save column index as an integer
rowcol <- raw.act1[, rowcol]
row <- sapply(rowcol, function (pos) substring(pos, 1, 1)) %>%
  match(LETTERS)
col <- sapply(rowcol, function (pos) substring(pos, 2)) %>%
  as.integer()
raw.act1[, `:=` (rowi = row, coli = col, rowcol = NULL)]

## create wllt column

## identify positive control, vehicle control, and test fish
raw.act1[, wllt := "t"]
raw.act1[cpid == "DMSO", wllt := "v"][cpid=="Chlropyrifos+", wllt := "p"]
raw.act1[wllt == "p", cpid := "Chlorpyrifos"]

## add source file and assay component information
srcf <- "OP Acute Final Sept 22 2021.xlsx"
acid <- "ZFlmrALD-6-10-10"
raw.act1[, `:=` (srcf = srcf, acid = acid)]

# perform final touches

## remove rows that are blank
raw.act2 <- raw.act1[!is.na(cpid)]

## ensure concentration and behavior data is of numeric class
raw.act2[, (c("conc", t.new)) := lapply(.SD, as.numeric), .SDcols = c("conc", t.new)]

# create lmr0 formatted table

## reorder columns and declare lmr0 object
lmr0.cols <- c("srcf", "acid", "cpid", "apid", "rowi",
               "coli", "wllt", "wllq", "conc", t.new)
lmr0.act <- raw.act2[, ..lmr0.cols]

# Find chemical|concentration groups and control groups with too many dead or malformed fish
lmr0.act.perc.wllq.0.rmv.conc <- lmr0.act[wllt=="t", .(perc.wllq.0 = length(which(wllq==0)) / .N), by = .(cpid,apid,conc)][perc.wllq.0 > 0.25, .(cpid,apid,conc)]
lmr0.act.perc.wllq.0.rmv.apid <- lmr0.act[wllt=="v", .(perc.wllq.0 = length(which(wllq==0)) / .N), by = .(cpid,apid,conc)][perc.wllq.0 > 0.20, apid]

# Save this information
save(lmr0.act.perc.wllq.0.rmv.conc, file = here("raw data","Acute Exposure","Padilla_OP_Acute_ChemConcs Excluded.Rdata"))
save(lmr0.act.perc.wllq.0.rmv.apid, file = here("raw data","Acute Exposure","Padilla_OP_Acute_Plates Excluded.Rdata"))

# Remove these chemical|concentration groups and plates
lmr0.act1 <- lmr0.act[!lmr0.act.perc.wllq.0.rmv.conc, on = .(cpid,apid,conc)][!(apid %in% lmr0.act.perc.wllq.0.rmv.apid)]

# save lmr0 formatted data
rm(lmr0.act)
lmr0.act <- data.table::copy(lmr0.act1)
# save lmr0 data
save(lmr0.act, file = "raw data/Acute Exposure/Padilla_OP_Acute_lmr0.Rdata")

rm(list = ls())

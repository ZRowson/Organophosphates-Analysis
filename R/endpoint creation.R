## ---------------------------
##
## Script Name: Create Endpoint Data for OP data
##
## Purpose of Script: Format Padilla OP data as mc0 prior to tcplfit2 analysis.
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

here::i_am("R/endpoint creation.R")


# Developmental Exposure --------------------------------------------------

load(here("raw data","Developmental Exposure","Padilla_OP_Dev_lmr0.Rdata"))

# remove positive control
data <- lmr0.dev[wllt!="p"]

## format as mc0 and transform to mc0_n with Box-Cox
ac <- c("AUC_L","avgS_L", "avgA_L", "avgJ_L","hbt_L",
        "strtlA","strtlAavg", "strtlF",
        "AUC_D", "avgS_D","avgA_D", "avgJ_D", "hbt_D",
        "AUC_T","avgS_T","AUC_r")
mc0.dev <- lapply(ac, function(endp) as.data.table(gabi::as_mc0(data, rval = endp, no.A = 0, no.L = 20, no.D = 20)))
names(mc0.dev) <- ac
mc0.dev_n <- lapply(mc0.dev, gabi::apply_bxcx)

# save mc0 and mc0_n
save(mc0.dev, file = here("pipelined data","Developmental Exposure","Padilla_OP_Dev_mc0.rda"))
save(mc0.dev_n, file = "pipelined data/Developmental Exposure/Padilla_OP_Dev_mc0_n.rda")

rm(list = ls())


# Acute Exposure ----------------------------------------------------------

load(here("raw data","Acute Exposure","Padilla_OP_Acute_lmr0.Rdata"))

# remove positive control
data <- lmr0.act[wllt!="p"]

## format as mc0 and transform to mc0_n with Box-Cox
ac <- c("AUC_L","avgS_L", "avgA_L", "avgJ_L","hbt_L",
        "strtlA","strtlAavg", "strtlF",
        "AUC_D", "avgS_D","avgA_D", "avgJ_D", "hbt_D",
        "AUC_T","avgS_T","AUC_r")
mc0.act <- lapply(ac, function(endp) as.data.table(gabi::as_mc0(data, rval = endp, no.A = 3, no.L = 5, no.D = 5)))
names(mc0.act) <- ac
mc0.act_n <- lapply(mc0.act, gabi::apply_bxcx)

# save mc0 and mc0_n
save(mc0.act, file = "pipelined data/Acute Exposure/Padilla_OP_Acute_mc0.rda")
save(mc0.act_n, file = "pipelined data/Acute Exposure/Padilla_OP_Acute_mc0_n.rda")


## ---------------------------
##
## Script Name: Create OP row Objects
##
## Purpose of Script: Create row objects from Padilla OP mc0 data prior to tcplfit2
##                    curve-fitting.
##
## Author: Zachary Rowson
##
## Date Created: 2022-02-14
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
library(here)
library(gabi)

here::i_am("R/row creation.R")


# Developmental Exposure --------------------------------------------------


# load endpoint data
load(here("pipelined data","Developmental Exposure","Padilla_OP_Dev_mc0_n.rda"))

# Create tcpl row objects

# Gather chemical names. Separate chemicals that have only one replicate or had a concentration group that was removed
chemicals <- mc0.dev_n[[1]][["data"]][wllt=="t",unique(cpid)]
names(chemicals) <- chemicals

# Create rows for chemicals with two replicates and no concentration groups removed
rows.dev_n <- lapply(chemicals, function(chm) {
  # Create row objects for each endpoint
  lapply(mc0.dev_n, function(list) {
    data <- list[["data"]]
    lam.hat <- list[["lam.hat"]]
    shift <- list[["shift"]]
    endp <- gsub("_ZFlmrALD-20-40-40","",data[,unique(acid)])
    gabi::as_row(data, endp=endp, chemical=chm, lam.hat = lam.hat, shift = shift)
  })
})

# Save row objects
save(rows.dev_n, file = "pipelined data/Developmental Exposure/Padilla_OP_Dev_rows_n.rda")


# Acute Exposure ----------------------------------------------------------


# load endpoint data
load(here("pipelined data","Acute Exposure","Padilla_OP_Acute_mc0_n.rda"))

# Create tcpl row objects

# Gather chemical names. Separate chemicals that have only one replicate or had a concentration group that was removed
chemicals <- mc0.act_n[[1]][["data"]][wllt=="t",unique(cpid)]
names(chemicals) <- chemicals

# Create rows for chemicals with two replicates and no concentration groups removed
rows.act_n <- lapply(chemicals, function(chm) {
  # Create row objects for each endpoint
  lapply(mc0.act_n, function(list) {
    data <- list[["data"]]
    lam.hat <- list[["lam.hat"]]
    shift <- list[["shift"]]
    endp <- gsub("_ZFlmrALD-6-10-10","",data[,unique(acid)])
    gabi::as_row(data, endp=endp, chemical=chm, lam.hat = lam.hat, shift = shift)
  })
})

# Save row objects
save(rows.act_n, file = "pipelined data/Acute Exposure/Padilla_OP_Acute_rows_n.rda")

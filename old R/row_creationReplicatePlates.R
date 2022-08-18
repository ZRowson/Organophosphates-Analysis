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
chms.loop <- mc0.dev_n[[1]][["data"]][wllt=="t",unique(cpid)]
names(chms.loop) <- chms.loop
not.fit <- c("Terbufos","Chlorpyrifos","Malathion","Ethoprop")
fit.later <- c("Trichlorfon","Phosmet")
to.fit <- chms.loop[!(names(chms.loop) %in% c(not.fit,fit.later))]

# Create rows for chemicals with two replicates and no concentration groups removed
rows.dev_n1 <- lapply(to.fit, function(chm) {
  # Create row objects for each endpoint
  lapply(mc0.dev_n, function(list) {
    data <- list[["data"]]
    lam.hat <- list[["lam.hat"]]
    shift <- list[["shift"]]
    endp <- gsub("_ZFpmrALD-20-40-40","",data[,unique(acid)])
    gabi::as_row(data, endp=endp, chemical=chm, lam.hat = lam.hat, shift = shift,
                 replicate.plate = TRUE)
  })
})

# Create rows for chemicals with one replicate and no concentration groups removed
rows.dev_n2 <- lapply(fit.later, function(chm) {
  # Create row objects for each endpoint
  lapply(mc0.dev_n, function(list) {
    data <- list[["data"]]
    lam.hat <- list[["lam.hat"]]
    shift <- list[["shift"]]
    endp <- gsub("_ZFpmrALD-20-40-40","",data[,unique(acid)])
    gabi::as_row(data, endp=endp, chemical=chm, lam.hat = lam.hat, shift = shift,
                 replicate.plate = FALSE)
  })
})
names(rows.dev_n2) <- fit.later

# Append lists together and save row objects
rows.dev_n <- append(rows.dev_n1, rows.dev_n2)
save(rows.dev_n, file = "pipelined data/Developmental Exposure/Padilla_OP_Dev_rows_n.rda")


# Acute Exposure ----------------------------------------------------------


# load endpoint data
load(here("pipelined data","Acute Exposure","Padilla_OP_Acute_mc0_n.rda"))

# Create tcpl row objects

# Gather chemical names. Separate chemicals that have only one replicate or had a concentration group that was removed
chms.loop <- mc0.act_n[[1]][["data"]][wllt=="t",unique(cpid)]
names(chms.loop) <- chms.loop
not.fit <- c("Acephate","Pirimiphos methyl","Tribufos")
fit.later <- c("Tebupirimiphos")
to.fit <- chms.loop[!(names(chms.loop) %in% c(not.fit,fit.later))]

# Create rows for chemicals with two replicates and no concentration groups removed
rows.act_n1 <- lapply(to.fit, function(chm) {
  # Create row objects for each endpoint
  lapply(mc0.act_n, function(list) {
    data <- list[["data"]]
    lam.hat <- list[["lam.hat"]]
    shift <- list[["shift"]]
    endp <- gsub("_ZFpmrALD-6-10-10","",data[,unique(acid)])
    gabi::as_row(data, endp=endp, chemical=chm, lam.hat = lam.hat, shift = shift,
                 replicate.plate = TRUE)
  })
})

# Create rows for chemicals with one replicate and no concentration groups removed
rows.act_n2 <- lapply(fit.later, function(chm) {
  # Create row objects for each endpoint
  lapply(mc0.act_n, function(list) {
    data <- list[["data"]]
    lam.hat <- list[["lam.hat"]]
    shift <- list[["shift"]]
    endp <- gsub("_ZFpmrALD-6-10-10","",data[,unique(acid)])
    gabi::as_row(data, endp=endp, chemical=chm, lam.hat = lam.hat, shift = shift,
                 replicate.plate = FALSE)
  })
})
names(rows.act_n2) <- fit.later

# Append lists together and save row objects
rows.act_n <- append(rows.act_n1, rows.act_n2)
save(rows.act_n, file = "pipelined data/Acute Exposure/Padilla_OP_Acute_rows_n.rda")

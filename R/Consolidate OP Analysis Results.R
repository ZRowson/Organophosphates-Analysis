## ---------------------------
##
## Script Name: Consolidate OP Anlayis Results
##
## Purpose of Script: Save a table full of summary data from tcplfit2 analysis and
##                    save graphs an .pngs in a graphics folder
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


# tcplfit2 output tables --------------------------------------------------


# Load developmental tcplfits
load("~/Desktop/Stephanie Projects/OP Analysis/pipelining/pipelined data/Developmental Exposure/Padilla_OP_Dev_tcplfits.rda")

# Create table of tcpl output data for developmental exposures
data <- lapply(tcplfits.dev_n, function(chm) {
  x <- lapply(chm, function(endp) {
        summary <- endp[["summary"]]
        summary[,!names(summary)%in%c("conc","resp")]})
  do.call('rbind', x)
})

OP_dev_tcpl_out <- do.call('rbind', data)

save(OP_dev_tcpl_out, file = "Results/Padilla_OP_Dev_tcpl_out.Rdata")


# Load acute tcplfits
load("~/Desktop/Stephanie Projects/OP Analysis/pipelining/pipelined data/Acute Exposure/Padilla_OP_Acute_tcplfits.rda")

# Create table of tcpl output data for Acute exposures
data <- lapply(tcplfits.act_n, function(chm) {
  x <- lapply(chm, function(endp) {
    summary <- endp[["summary"]]
    summary[,!names(summary)%in%c("conc","resp")]})
  do.call('rbind', x)
})

OP_act_tcpl_out <- do.call('rbind', data)

save(OP_act_tcpl_out, file = "Results/Padilla_OP_Acute_tcpl_out.Rdata")


# Save graphics -----------------------------------------------------------


setwd("../graphics")

# loop through chemicals for developmental chemicals
dir.create("Developmental Exposures")
setwd("Developmental Exposures")
chms <- names(tcplfits.dev_n)
lapply(chms, function (chm) {
  dir.create(chm)
  setwd(chm)
  to_fit <- tcplfits.dev_n[[chm]]
  lapply(to_fit, function(row) {
    summary <- row[["summary"]]
    png(paste0(paste("concResp", summary$name, summary$assay, sep= "_"),
               ".png"),
        width = 720,
        height = 720)
    print(row[["plot"]])
    dev.off()
  })
  setwd("..")
})

# loop through chemicals for Acute chemicals
dir.create("Acute Exposures")
setwd("Acute Exposures")
chms <- names(tcplfits.act_n)
lapply(chms, function (chm) {
              dir.create(chm)
              setwd(chm)
              to_fit <- tcplfits.act_n[[chm]]
              lapply(to_fit, function(row) {
                              summary <- row[["summary"]]
                              png(paste0(paste("concResp", summary$name, summary$assay, sep= "_"),
                                         ".png"),
                                  width = 720,
                                  height = 720)
                              print(row[["plot"]])
                              dev.off()
                              })
              setwd("..")
              })

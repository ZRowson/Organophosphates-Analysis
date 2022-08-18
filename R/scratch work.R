library(gabi)

pmr0.df <- as.data.frame(pmr0)
gabi::as_mc0(pmr0.df, rval = "AUC_D")

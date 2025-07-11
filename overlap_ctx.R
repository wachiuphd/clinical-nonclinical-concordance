library(ctxR)
library(readxl)
my_key <- ctx_key()
clinical.1 <- read_excel("Weitekamp et al. 2025-Supplemental Tables.xlsx",
                           sheet = "Table S3")
clinical.1.drugs <- unique(clinical.1$drug)

clinical.1.ctx.exact <- chemical_equal_batch(clinical.1.drugs)

found.1 <- str_replace_all(tolower(clinical.1.drugs),"-"," ") %in% tolower(clinical.1.ctx.exact$valid$searchValue)

clinical.1.ctx.exact.nameonly <- chemical_equal_batch(str_split_i(clinical.1.drugs[!found.1]," ",1))

clinical.1.ctx <- rbind(cbind(data.frame(drug=clinical.1.drugs[found.1]),clinical.1.ctx.exact$valid),
                        cbind(data.frame(drug=clinical.1.drugs[!found.1]),clinical.1.ctx.exact.nameonly$valid))

clinical.2 <- read_excel("Weitekamp et al. 2025-Supplemental Tables.xlsx",
                           sheet = "Table S4")
clinical.2.drugs <- unique(clinical.2$drug)

clinical.2.ctx.exact <- chemical_equal_batch(clinical.2.drugs)

found.2 <- str_replace_all(tolower(clinical.2.drugs),"-"," ") %in% tolower(clinical.2.ctx.exact$valid$searchValue)

# Sodium zirconium not found - use name from clinical.1
clinical.2.ctx.exact.NaZirconium <- chemical_equal_batch("sodium zirconium cyclosilicate")

clinical.2.ctx <- rbind(cbind(data.frame(drug=clinical.2.drugs[found.2]),clinical.2.ctx.exact$valid),
                        cbind(data.frame(drug=clinical.2.drugs[!found.2]),clinical.2.ctx.exact.NaZirconium$valid))


kvasnicka.dat <- as_tibble(read.csv(file.path("..","!Projects-Current","!USETOX","QSAR-PPRTV-Comparison","Two-Stage-ML-Results-Browser",
                                              "PODs800k.csv")))
kvasnicka.dat <- rename(kvasnicka.dat, casrn = CAS_RN_Dashboard)
kvasnicka.dat <- rename(kvasnicka.dat, dtxsid = DTXSID)
kvasnicka.dat$pod.min <- pmin(kvasnicka.dat$pod.gen,kvasnicka.dat$pod.repdev,na.rm=T)
kvasnicka.dat$pod.min.lb <- pmin(kvasnicka.dat$lb.gen,kvasnicka.dat$lb.repdev,na.rm=T)

clinical.1.ctx <- left_join(clinical.1.ctx,kvasnicka.dat,by="dtxsid")
clinical.2.ctx <- left_join(clinical.2.ctx,kvasnicka.dat,by="dtxsid")

write.csv(clinical.1.ctx,"Clinical Dataset 1 CTX.csv",row.names = FALSE)
write.csv(clinical.2.ctx,"Clinical Dataset 2 CTX.csv",row.names = FALSE)

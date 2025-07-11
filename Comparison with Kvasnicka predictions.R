library(ctxR)
library(tidyverse)
library(readxl)
my_key <- ctx_key()

endpoints <- read_excel("Weitekamp et al. 2025-Supplemental Tables.xlsx",
                        sheet = "Table S6")
endpoints.keep <- subset(endpoints,category=="NA")


clinical.1 <- read_excel("Weitekamp et al. 2025-Supplemental Tables.xlsx",
                         sheet = "Table S3")
clinical.1$effect <- trimws(clinical.1$effect) # one entry as extra ws at end
clinical.1$keep <- clinical.1$effect %in% endpoints.keep$effect
clinical.1$LOAEL <- as.numeric(clinical.1$mg_kg_d)
clinical.1.min <- aggregate(LOAEL ~ drug,data=subset(clinical.1,keep),FUN=min,na.rm=T) # min of "keep" endpoints
clinical.1.ctx <- read.csv("Clinical Dataset 1 CTX.csv")
clinical.1.min <- left_join(clinical.1.min,clinical.1.ctx,by="drug")
ggplot(clinical.1.min)+geom_point(aes(x=LOAEL,y=pod.min))+
  scale_x_log10(limits=c(1e-4,400))+scale_y_log10(limits=c(1e-4,400))+geom_abline(slope=1,intercept = 0)

clinical.2 <- read_excel("Weitekamp et al. 2025-Supplemental Tables.xlsx",
                         sheet = "Table S4")
clinical.2$keep <- clinical.2$effect %in% endpoints.keep$effect
clinical.2$LOAEL <- as.numeric(clinical.2$mg_kg_d)
clinical.2.min <- aggregate(LOAEL ~ drug,data=subset(clinical.2,keep),FUN=min,na.rm=T)
clinical.2.ctx <- read.csv("Clinical Dataset 2 CTX.csv")
clinical.2.min <- left_join(clinical.2.min,clinical.2.ctx,by="drug")
ggplot(clinical.2.min)+geom_point(aes(x=LOAEL,y=pod.min))+
  scale_x_log10(limits=c(1e-3,125))+scale_y_log10(limits=c(1e-3,125))+geom_abline(slope=1,intercept = 0)



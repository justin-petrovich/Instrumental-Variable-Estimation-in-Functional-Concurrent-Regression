### Load Data
load("LSE Data - cleaned.Rdata")

### Load packages
library(refund)
library(fdapace)
library(face)
library(tidyverse)

### Create grid
grid.out <- min(elastdf$year):max(elastdf$year)
K <- length(grid.out)

### Choose a subsample of the entire data set due to size
ids <- sort(unique(elastdf$id))

### randomly choose N people to subsample
set.seed(123)
N <- 50000

### female, married
fm_ids <- unique(filter(elastdf, sex==2, married==1)[['id']])
fm_samp <- sample(fm_ids,size = ifelse(length(fm_ids) < N,length(fm_ids),N),replace = F)
subdata_fm <- elastdf %>%
  filter(id %in% fm_samp)

### female, unmarried
fu_ids <- unique(filter(elastdf, sex==2, married==0)[['id']])
fu_samp <- sample(fu_ids,size = ifelse(length(fu_ids) < N,length(fu_ids),N),replace = F)
subdata_fu <- elastdf %>%
  filter(id %in% fu_samp)

### male, married
mm_ids <- unique(filter(elastdf, sex==1, married==1)[['id']])
mm_samp <- sample(mm_ids,size = ifelse(length(mm_ids) < N,length(mm_ids),N),replace = F)
subdata_mm <- elastdf %>%
  filter(id %in% mm_samp)

### male, unmarried
mu_ids <- unique(filter(elastdf, sex==1, married==0)[['id']])
mu_samp <- sample(mu_ids,size = ifelse(length(mu_ids) < N,length(mu_ids),N),replace = F)
subdata_mu <- elastdf %>%
  filter(id %in% mu_samp)

### Join groups into one subdata set
subdata <- full_join(subdata_fm,subdata_fu)
subdata <- full_join(subdata,subdata_mm)
subdata <- full_join(subdata,subdata_mu)
rm(subdata_fm,subdata_fu,subdata_mm,subdata_mu)
subdata <- subdata %>% rename(subj = id,argvals = year,y = lannhrs)

### Run separate models for male/female, married/unmarried
elast_est <- list("male_un" = NA,"male_mar" = NA,
                  "female_un" = NA,"female_mar" = NA)
source("Unmarried Males Model.R")
elast_est[["male_un"]] <- list(stage2 = stage2, lse = beta1hat_2sls)
source("Married Males Model.R")
elast_est[["male_mar"]] <- list(stage2 = stage2, lse = beta1hat_2sls)
source("Unmarried Females Model.R")
elast_est[["female_un"]] <- list(stage2 = stage2, lse = beta1hat_2sls)
source("Married Females Model.R")
elast_est[["female_mar"]] <- list(stage2 = stage2, lse = beta1hat_2sls)
# save(elast_est,file = "TSLS_models.RData")

### Plot estimates (Figure 4)
ylim = c(-2,1.5)
# png("LSE Estimates.png", width = 1200,height = 800)
par(mfrow = c(2,2),mar=c(5,5,4,1)+.1)
plot(elast_est$male_mar$stage2,select = 2,se = F,scale = 0,xlab = "Year",ylab = expression(hat(alpha)(t)),main = "Married Males",cex.lab = 1.25,cex.axis= 1.1,ylim = ylim,lwd = 2)
abline(h = 0,lty = 2)
plot(elast_est$male_un$stage2,select = 2,se = F,scale = 0,xlab = "Year",ylab = expression(hat(alpha)(t)),main = "Unmarried Males",cex.lab = 1.25,cex.axis= 1.1,ylim = ylim,lwd = 2)
abline(h = 0,lty = 2)
plot(elast_est$female_mar$stage2,select = 2,se = F,scale = 0,xlab = "Year",ylab = expression(hat(alpha)(t)),main = "Married Females",cex.lab = 1.25,cex.axis= 1.1,ylim = ylim,lwd = 2)
abline(h = 0,lty = 2)
plot(elast_est$female_un$stage2,select = 2,se = F,scale = 0,xlab = "Year",ylab = expression(hat(alpha)(t)),main = "Unmarried Females",cex.lab = 1.25,cex.axis= 1.1,ylim = ylim,lwd = 2)
abline(h = 0,lty = 2)
# dev.off()
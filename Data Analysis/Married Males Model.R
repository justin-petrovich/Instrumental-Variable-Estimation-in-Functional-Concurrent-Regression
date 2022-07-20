## Married Males
data <- subdata %>% filter(sex == 1, married == 1)

id <- sort(sample(unique(data$subj),replace = T))

lens <- tabulate(data$subj)[id]
idx <- sequence(lens,
                match(unique(data$subj),data$subj)[match(id,unique(data$subj))])
dfb <- data[idx,]
dfb <- dfb %>%
  mutate(subj_old = subj,
         subj = rep(seq_along(id), lens)) %>%
  relocate(subj_old,.after = subj)

#-------------------------------------------#
## Step 1: "Impute" sparse functional data ##
#-------------------------------------------#
Nstar <- ifelse(nrow(dfb)/2 < N, nrow(dfb)/2, N)
# Nstar <- dfb %>% group_by(subj) %>% n_groups()

## use FACE to estimate full curves
ydata <- dfb[,c("subj","argvals","y")]
latwage.data <- dfb[,c("subj","argvals","latwage")]; names(latwage.data)[3] <- "y"
nonwage.data <- dfb[,c("subj","argvals","pernonwage")]; names(nonwage.data)[3] <- "y"
age.data <- dfb[,c("subj","argvals","age")]; names(age.data)[3] <- "y"
kid.data <- dfb[,c("subj","argvals","kids")]; names(kid.data)[3] <- "y"

newdata <- data.frame(subj = rep(dfb$subj,each = K),
                      argvals = rep(grid.out,times = Nstar),
                      y = rep(NA,Nstar*K))
fit.y <- face.sparse(data = ydata,newdata = rbind(ydata,newdata),
                     argvals.new = grid.out,knots = 7)
fit.wage <- face.sparse(data = latwage.data,newdata = rbind(latwage.data,newdata),
                        argvals.new = grid.out,knots = 7)
fit.nonwage <- face.sparse(data = nonwage.data,newdata = rbind(nonwage.data,newdata),
                           argvals.new = grid.out,knots = 7)
fit.age <- face.sparse(data = age.data,newdata = rbind(age.data,newdata),
                       argvals.new = grid.out,knots = 7)

Ydense <- matrix(fit.y$y.pred[-c(1:nrow(ydata))],nrow = Nstar,ncol = K,byrow = T)
wagedense <- matrix(fit.wage$y.pred[-c(1:nrow(latwage.data))],nrow = Nstar,ncol = K,byrow = T)
nonwagedense <- matrix(fit.nonwage$y.pred[-c(1:nrow(nonwage.data))],nrow = Nstar,ncol = K,byrow = T)
agedense <- matrix(fit.age$y.pred[-c(1:nrow(age.data))],nrow = Nstar,ncol = K,byrow = T)

zcols <- c("agesq")
data1 <- vector(mode = "list",length = 3 + length(zcols))
data1[[1]] <- wagedense
data1[[2]] <- nonwagedense
data1[[3]] <- agedense
names(data1) <- c("latwage","pernonwage","age","agesq")

for(i in 1:length(zcols)){
  zdata <- dfb[,c("subj","argvals",zcols[i])]; names(zdata)[3] <- "y" 
  fit.z <- face.sparse(data = zdata,newdata = rbind(zdata,newdata),
                       argvals.new = grid.out,knots = 7)
  Zdense <- matrix(fit.z$y.pred[-c(1:nrow(zdata))],nrow = Nstar,ncol = K,byrow = T)
  data1[[i+3]] <- Zdense
}

data1$edu <- dfb %>% 
  group_by(subj) %>% 
  summarize(edu = first(edu)) %>% 
  pull(edu)

data1$occ <- dfb %>% 
  group_by(subj) %>% 
  summarize(occ = first(occ_short)) %>% 
  pull(occ)

data1$kids <- dfb %>% 
  group_by(subj) %>% 
  summarize(kids = first(kids)) %>% 
  pull(kids)

#------------------------------------------------------------------------------#
## Step 2.1: regress hwage on controls and instruments and obtain \hat{hwage} ##
#------------------------------------------------------------------------------#
# Regression stage 1
stage1 <- pffr(latwage ~ pernonwage + kids + age + c(occ) + c(edu) + agesq,
               data = data1,yind = grid.out)
latwagehat <- fitted(stage1)

#-------------------------------------------#
## Step 2.2: regress annhrs on \hat{hwage} ##
#-------------------------------------------#
# Regression stage 2
data2 <- list(lannhrs = Ydense,latwagehat = latwagehat, pernonwage = nonwagedense,
              kids = data1$kids, age = agedense,
              edu = data1$edu, occ = data1$occ)
stage2 <- pffr(lannhrs ~ latwagehat + pernonwage + kids + age + c(occ) + c(edu),
               data = data2,yind = grid.out,
               bs.yindex = list(bs = "ps", k = 10, m = c(2,1)))
beta1hat_2sls <- coef(stage2,n1 = K)$smterms[[2]][['value']]

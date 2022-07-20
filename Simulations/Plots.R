library(tidyverse)
library(ggthemes)
library(gridExtra)

# MISE Plot
results$msedf %>%
  pivot_longer(cols = 1:4,names_to = "Method",values_to = "MISE") %>%
  ggplot(.) +
  geom_boxplot(aes(y = MISE,x = Method,col = Method)) +
  theme_bw() +
  theme(panel.grid = element_blank())


# # Beta Plot
# betas1 <- as.data.frame(do.call("cbind",results$betas$s1))
# betas1$beta1 <- results$beta1
# betas1$ts <- results$ts
# betas1$sim <- "Simulation 1"
# 
# betas2 <- as.data.frame(do.call("cbind",results$betas$s2))
# betas2$beta1 <- results$beta1
# betas2$ts <- results$ts
# betas2$sim <- "Simulation 2"
# 
# betas3 <- as.data.frame(do.call("cbind",results$betas$s3))
# betas3$beta1 <- results$beta1
# betas3$ts <- results$ts
# betas3$sim <- "Simulation 3"
# 
# betas <- rbind(betas1,betas2,betas3)
# colnames(betas)[which(colnames(betas)=="beta1")] <- "truth"
# method.names <- c("naive.fcr","naive.fcreg","fc2sls.face","fc2sls.pace","truth")
# betas <- betas %>%
#   select(naive.fcr,naive.fcreg,fc2sls.face,fc2sls.pace,truth,ts,sim) %>%
#   pivot_longer(cols = 1:5, names_to = "Method",values_to = "Value") %>%
#   filter(Method %in% method.names) %>%
#   mutate(Method = factor(Method,levels = method.names))
# 
# method.colors <- c("red2","red2","blue3","blue3","black")
# method.type <- c("longdash","dotted","dashed","dotdash","solid")
# 
# 
# g <- ggplot(betas) +
#   geom_line(aes(x = ts,y = Value,color = Method,linetype = Method,size = Method)) +
#   facet_wrap(vars(sim)) +
#   scale_color_manual(values = method.colors) +
#   scale_size_manual(values = c(rep(0.75,4),1.5)) +
#   scale_linetype_manual(values = method.type) +
#   labs(y = "Beta 1", x = "t") +
#   # lims(y = c(-5,11)) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.title = element_text(size = 16),
#         axis.text = element_text(size = 14),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 16),
#         legend.title.align = 0.5,
#         strip.text = element_text(size = 16)
#         # strip.background = element_blank(),
#         )
# 
# png(filename = "betas_s3sc1_N800.png",width = 1000,height = 400)
# g
# dev.off()



# Beta Plot - all combined

load("./Scenario 1/Scheme 1/mse-N100-s1scheme1.RData")
s1design1_N100 <- results
load("./Scenario 1/Scheme 1/mse-N200-s1scheme1.RData")
s1design1_N200 <- results
load("./Scenario 1/Scheme 1/mse-N400-s1scheme1.RData")
s1design1_N400 <- results
load("./Scenario 1/Scheme 1/mse-N800-s1scheme1.RData")
s1design1_N800 <- results
load("./Scenario 2/Scheme 1/mse-N100-s2scheme1.RData")
s2design1_N100 <- results
load("./Scenario 2/Scheme 1/mse-N200-s2scheme1.RData")
s2design1_N200 <- results
load("./Scenario 2/Scheme 1/mse-N400-s2scheme1.RData")
s2design1_N400 <- results
load("./Scenario 2/Scheme 1/mse-N800-s2scheme1.RData")
s2design1_N800 <- results
load("./Scenario 3/Scheme 1/mse-N100-s3scheme1.RData")
s3design1_N100 <- results
load("./Scenario 3/Scheme 1/mse-N200-s3scheme1.RData")
s3design1_N200 <- results
load("./Scenario 3/Scheme 1/mse-N400-s3scheme1.RData")
s3design1_N400 <- results
load("./Scenario 3/Scheme 1/mse-N800-s3scheme1.RData")
s3design1_N800 <- results
rm(results)


betas1 <- as.data.frame(do.call("cbind",s1design1_N100$betas$s1))
betas1 <- rbind(betas1, as.data.frame(do.call("cbind",s1design1_N200$betas$s1)))
betas1 <- rbind(betas1, as.data.frame(do.call("cbind",s1design1_N400$betas$s1)))
betas1 <- rbind(betas1, as.data.frame(do.call("cbind",s1design1_N800$betas$s1)))
betas1 <- rbind(betas1, as.data.frame(do.call("cbind",s2design1_N100$betas$s1)))
betas1 <- rbind(betas1, as.data.frame(do.call("cbind",s2design1_N200$betas$s1)))
betas1 <- rbind(betas1, as.data.frame(do.call("cbind",s2design1_N400$betas$s1)))
betas1 <- rbind(betas1, as.data.frame(do.call("cbind",s2design1_N800$betas$s1)))
betas1 <- rbind(betas1, as.data.frame(do.call("cbind",s3design1_N100$betas$s1)))
betas1 <- rbind(betas1, as.data.frame(do.call("cbind",s3design1_N200$betas$s1)))
betas1 <- rbind(betas1, as.data.frame(do.call("cbind",s3design1_N400$betas$s1)))
betas1 <- rbind(betas1, as.data.frame(do.call("cbind",s3design1_N800$betas$s1)))
betas1$beta1 <- c(s1design1_N100$beta1,s1design1_N200$beta1,s1design1_N400$beta1,s1design1_N800$beta1,
                  s2design1_N100$beta1,s2design1_N200$beta1,s2design1_N400$beta1,s2design1_N800$beta1,
                  s3design1_N100$beta1,s3design1_N200$beta1,s3design1_N400$beta1,s3design1_N800$beta1)
betas1$ts <- c(s1design1_N100$ts,s1design1_N200$ts,s1design1_N400$ts,s1design1_N800$ts,
               s2design1_N100$ts,s2design1_N200$ts,s2design1_N400$ts,s2design1_N800$ts,
               s3design1_N100$ts,s3design1_N200$ts,s3design1_N400$ts,s3design1_N800$ts)
betas1$scenario <- as.factor(rep(1:3,each = 400))
betas1$N <- rep(rep(c(100,200,400,800),each = 100),times = 3)
betas1$sim <- "(a)"


betas2 <- as.data.frame(do.call("cbind",s1design1_N100$betas$s2))
betas2 <- rbind(betas2, as.data.frame(do.call("cbind",s1design1_N200$betas$s2)))
betas2 <- rbind(betas2, as.data.frame(do.call("cbind",s1design1_N400$betas$s2)))
betas2 <- rbind(betas2, as.data.frame(do.call("cbind",s1design1_N800$betas$s2)))
betas2 <- rbind(betas2, as.data.frame(do.call("cbind",s2design1_N100$betas$s2)))
betas2 <- rbind(betas2, as.data.frame(do.call("cbind",s2design1_N200$betas$s2)))
betas2 <- rbind(betas2, as.data.frame(do.call("cbind",s2design1_N400$betas$s2)))
betas2 <- rbind(betas2, as.data.frame(do.call("cbind",s2design1_N800$betas$s2)))
betas2 <- rbind(betas2, as.data.frame(do.call("cbind",s3design1_N100$betas$s2)))
betas2 <- rbind(betas2, as.data.frame(do.call("cbind",s3design1_N200$betas$s2)))
betas2 <- rbind(betas2, as.data.frame(do.call("cbind",s3design1_N400$betas$s2)))
betas2 <- rbind(betas2, as.data.frame(do.call("cbind",s3design1_N800$betas$s2)))
betas2$beta1 <- c(s1design1_N100$beta1,s1design1_N200$beta1,s1design1_N400$beta1,s1design1_N800$beta1,
                  s2design1_N100$beta1,s2design1_N200$beta1,s2design1_N400$beta1,s2design1_N800$beta1,
                  s3design1_N100$beta1,s3design1_N200$beta1,s3design1_N400$beta1,s3design1_N800$beta1)
betas2$ts <- c(s1design1_N100$ts,s1design1_N200$ts,s1design1_N400$ts,s1design1_N800$ts,
               s2design1_N100$ts,s2design1_N200$ts,s2design1_N400$ts,s2design1_N800$ts,
               s3design1_N100$ts,s3design1_N200$ts,s3design1_N400$ts,s3design1_N800$ts)
betas2$scenario <- as.factor(rep(1:3,each = 400))
betas2$N <- rep(rep(c(100,200,400,800),each = 100),times = 3)
betas2$sim <- "(b)"


betas3 <- as.data.frame(do.call("cbind",s1design1_N100$betas$s3))
betas3 <- rbind(betas3, as.data.frame(do.call("cbind",s1design1_N200$betas$s3)))
betas3 <- rbind(betas3, as.data.frame(do.call("cbind",s1design1_N400$betas$s3)))
betas3 <- rbind(betas3, as.data.frame(do.call("cbind",s1design1_N800$betas$s3)))
betas3 <- rbind(betas3, as.data.frame(do.call("cbind",s2design1_N100$betas$s3)))
betas3 <- rbind(betas3, as.data.frame(do.call("cbind",s2design1_N200$betas$s3)))
betas3 <- rbind(betas3, as.data.frame(do.call("cbind",s2design1_N400$betas$s3)))
betas3 <- rbind(betas3, as.data.frame(do.call("cbind",s2design1_N800$betas$s3)))
betas3 <- rbind(betas3, as.data.frame(do.call("cbind",s3design1_N100$betas$s3)))
betas3 <- rbind(betas3, as.data.frame(do.call("cbind",s3design1_N200$betas$s3)))
betas3 <- rbind(betas3, as.data.frame(do.call("cbind",s3design1_N400$betas$s3)))
betas3 <- rbind(betas3, as.data.frame(do.call("cbind",s3design1_N800$betas$s3)))
betas3$beta1 <- c(s1design1_N100$beta1,s1design1_N200$beta1,s1design1_N400$beta1,s1design1_N800$beta1,
                  s2design1_N100$beta1,s2design1_N200$beta1,s2design1_N400$beta1,s2design1_N800$beta1,
                  s3design1_N100$beta1,s3design1_N200$beta1,s3design1_N400$beta1,s3design1_N800$beta1)
betas3$ts <- c(s1design1_N100$ts,s1design1_N200$ts,s1design1_N400$ts,s1design1_N800$ts,
               s2design1_N100$ts,s2design1_N200$ts,s2design1_N400$ts,s2design1_N800$ts,
               s3design1_N100$ts,s3design1_N200$ts,s3design1_N400$ts,s3design1_N800$ts)
betas3$scenario <- as.factor(rep(1:3,each = 400))
betas3$N <- rep(rep(c(100,200,400,800),each = 100),times = 3)
betas3$sim <- "(c)"


betas <- rbind(betas1,betas2,betas3)
colnames(betas)[which(colnames(betas)=="beta1")] <- "Truth"
colnames(betas)[which(colnames(betas)=="naive.fcr")] <- "Naive-f"
colnames(betas)[which(colnames(betas)=="naive.fcreg")] <- "Naive-p"
colnames(betas)[which(colnames(betas)=="fc2sls.face")] <- "FC2SLS-f"
colnames(betas)[which(colnames(betas)=="fc2sls.pace")] <- "FC2SLS-p"
betas <- betas %>%
  select(`Naive-f`,`Naive-p`,`FC2SLS-f`,`FC2SLS-p`,Truth,ts,N,scenario,sim) %>%
  pivot_longer(cols = 1:5, names_to = "Method",values_to = "Value") %>%
  mutate(Method = as.factor(Method),
         sim = as.factor(sim),
         N = as.factor(N)) %>%
  arrange(N,sim,scenario,Method)

method.colors <- c("red2","red2","blue3","blue3","black")
method.type <- c("longdash","dotted","dashed","dotdash","solid")


g <- 
  betas %>%
  filter(scenario==2) %>%
  ggplot(.) +
  geom_line(aes(x = ts,y = Value,group = Method,color = Method,linetype = Method,size = Method)) +
  facet_grid(rows = vars(N),cols = vars(sim)) +
  scale_color_manual(values = method.colors) +
  scale_size_manual(values = c(rep(0.75,4),1.5)) +
  scale_linetype_manual(values = method.type) +
  labs(y = expression(alpha(t)), x = "t") +
  lims(y = c(-5,10)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.title.align = 0.5,
        strip.text = element_text(size = 16)
        # strip.background = element_blank(),
  )

g

g2 <- 
  betas %>%
  filter(scenario==3) %>%
  ggplot(.) +
  geom_line(aes(x = ts,y = Value,group = Method,color = Method,linetype = Method,size = Method)) +
  facet_grid(rows = vars(N),cols = vars(sim)) +
  scale_color_manual(values = method.colors) +
  scale_size_manual(values = c(rep(0.75,4),1.5)) +
  scale_linetype_manual(values = method.type) +
  labs(y = expression(alpha(t)), x = "t") +
  lims(y = c(-3,4)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.title.align = 0.5,
        strip.text = element_text(size = 16)
        # strip.background = element_blank(),
  )

g2


png(filename = "example_sims_s2sc1.png",width = 1200,height = 1000)
g
dev.off()

png(filename = "example_sims_s3sc1.png",width = 1200,height = 1000)
g2
dev.off()

load("LSE Data - 1988 2015.Rdata")

# Load packages
library(tidyverse)

summary(elastdf)
# notice annhrs has a negative value
elastdf %>%
  filter(annhrs1 < 0 | annhrs2 < 0) %>%
  select(annhrs1,annhrs2)
# people with negative values for annhrs. Drop these observations from the sample.

# atwage1 and atwage2 also have negative values
elastdf %>%
  filter(atwage1 < 0 | atwage2 < 0) %>%
  select(atwage1,atwage2)

# pernonwage1 and pernonwage2 also have negative values
elastdf %>%
  filter(pernonwage1 < 0 | pernonwage2 < 0) %>%
  select(pernonwage1,pernonwage2)

# many missing values for atwage1 and atwage2
elastdf %>%
  filter(is.na(atwage1) | is.na(atwage2)) %>%
  select(contains("wage"))

# missing wage data per year
elastdf %>%
  group_by(year) %>%
  summarise(n = n(),
            atwage1.mis = sum(is.na(atwage1))/n()*100,
            atwage2.mis = sum(is.na(atwage2))/n()*100,
            frate1.mis = sum(is.na(frate1))/n()*100,
            frate2.mis = sum(is.na(frate2))/n()*100,
            srate1.mis = sum(is.na(srate1))/n()*100,
            srate2.mis = sum(is.na(srate2))/n()*100,
            ficar1.mis = sum(is.na(ficar1))/n()*100,
            ficar2.mis = sum(is.na(ficar2))/n()*100,
            hwage1.mis = sum(is.na(hwage1))/n()*100,
            hwage2.mis = sum(is.na(hwage2))/n()*100,
            diff1 = atwage1.mis - hwage1.mis,
            diff2 = atwage2.mis - hwage2.mis) %>%
  print(n = 32)
### Note: atwage data is missing more than hwage data
### though both are missing a lot of data

# occupation recoding function
### These are based on groupings shown in the Occ groups.xlsx file
occ_rec <- function(occ){
  out <- rep(NA,length(occ))
  out[which(occ >= 3 & occ <= 22)] <- 1      # executive, admin, managerial
  out[which(occ >= 23 & occ <= 37)] <- 2     # management related occs
  out[which(occ >= 43 & occ <= 200)] <- 3    # professional specialty occs
  out[which(occ >= 203 & occ <= 235)] <- 4   # technicians & related support occs
  out[which(occ >= 243 & occ <= 258)] <- 5   # financial sales and related occs
  out[which(occ >= 274 & occ <= 283)] <- 6   # retail sales occs
  out[which(occ >= 303 & occ <= 389)] <- 7   # administrative support occs
  out[which(occ >= 405 & occ <= 472)] <- 8   # personal service occs
  out[which(occ >= 473 & occ <= 475)] <- 9   # farm operators and managers
  out[which(occ >= 479 & occ <= 498)] <- 10  # other agricultural and related occs
  out[which(occ >= 503 & occ <= 549)] <- 11  # mechanics and repairers
  out[which(occ >= 558 & occ <= 599)] <- 12  # construction trades
  out[which(occ >= 614 & occ <= 617)] <- 13  # extractive occs
  out[which(occ >= 628 & occ <= 699)] <- 14  # precision production occs
  out[which(occ >= 703 & occ <= 799)] <- 15  # machine operators, assemblers, inspectors
  out[which(occ >= 803 & occ <= 889)] <- 16  # transportation and material moving occs
  as.factor(out)
}

elastdf <-
  elastdf %>%
  select(id,annhrs1,annhrs2,hwage1,hwage2,age,sex,married,kids,year1,year2,
         pernonwage1,pernonwage2,atwage1,atwage2,occ1990dd,edu) %>%
  filter(age >= 25 & age <= 65,
         !is.na(annhrs1), !is.na(annhrs2),
         !is.na(hwage1), !is.na(hwage2),
         !is.na(occ1990dd),
         annhrs1 > 0, annhrs2 > 0,
         pernonwage1 >= 0, pernonwage2 >= 0,
         atwage1 > 0, atwage2 > 0) %>%
  pivot_longer(!c(id,sex,age,married,kids,occ1990dd,edu),
               names_to = c(".value","obs_num"),
               names_pattern = "(.*)(.)") %>%
  mutate(occ_short = occ_rec(occ1990dd),
         edu = as.factor(edu),
         agesq = age^2,
         lhwage = log(hwage),
         latwage = log(atwage),
         lannhrs = log(annhrs))
# save(elastdf,file = "LSE Data - cleaned.Rdata")


# Total number of respondents included in data set
length(unique(elastdf$id))
## 245531

# Number of respondents per year
elastdf %>%
  group_by(year) %>%
  count

# Figure 1
library(tidyr)
library(ggthemes)
# png(filename = "Obs_per_year.png",width = 1000,height = 800)
elastdf %>%
  group_by(year) %>%
  count %>%
  ggplot(.,aes(x = year,y = n)) +
  geom_bar(stat = "identity",color = "black",fill = "darkblue") +
  labs(x = "Year",y = "Respondents") +
  theme_clean() +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))
# dev.off()
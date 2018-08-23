rm(list=ls())

##READ THIS!###
#meta-data relevant to this file:
#species_origin = which spp did the dust originate from OR which species is the focal individual
#species_contacted = which spp received dust OR which species are contacting the focal individual
###mylu = M. lucifugus
###pesu = P. subflavus
###myse = M. septentrionalis
#sim_sum = the fraction of each species (from the species contacted column) at each site that an individual was connected to via dust, shared groups (clustering), or physical contact. This includes the 6 rearragenements (arousals) over the winter for physical contact and shared groups
#species_contacted_count = the count each species at a site (from the species_contacted column)
#dat.type = the parameter data type (e.g. dust, shared groups, physical contact, but coded differently)
### dust = UVF-dust
### cluster = physical contact
### neigh = shared groups


####### Contact data#########
library("brms")
library("rstan")
library(emmeans)
library(tidyverse)

#read in file - this includes the rearrangements for shared groups and physical contact
d1=read.csv("sim_comb_df6.csv")
#give number of "successes" (e.g. infections) and number of trials
d1$inf = d1$sim_sum*d1$species_contacted_count
d1$inf = round(d1$inf) #round these
d2=d1 #preserve d1 in case we make a mistake below
d2$inf [d2$inf=="NA"]=NA #if nothing in the inf column, NA
d2 = d2[!(is.na(d2$inf)),]

#make a data frame for each origin species
pesu.dat=subset(d2, species_origin=="pesu")
myse.dat=subset(d2, species_origin=="myse")
mylu.dat=subset(d2, species_origin=="mylu")

######PESU#######
pesu.dat$species_contacted=relevel(pesu.dat$species_contacted, ref="pesu") #relevel for easier contrast assessment
pesu.dat$dat.type=relevel(pesu.dat$dat.type, ref="dust")

get_prior(inf|trials(species_contacted_count) ~ dat.type * species_contacted + (1|site),
          data = pesu.dat, family = binomial()) #look at the default priors

prior=set_prior("normal(0,5)", class = "b") 
#set prior for all parameters of type "b" to be a normal distribution with sd of 5

#warning! this is the large bayesian analysos - it takes ~20 min to run on a 2017 macbook air
pesu2 <- brm(formula = inf|trials(species_contacted_count) ~ dat.type * species_contacted  + (1|site),
            data = pesu.dat, family = binomial(), prior = prior,chains = 4,
            control = list(adapt_delta = 0.97),cores=2,
            save_all_pars = TRUE);summary(pesu2, waic = TRUE)

marginal_effects(pesu2,"dat.type:species_contacted")# probs = c(0.4, 0.6)
#get the posteriors in response space (instead of contrasts)
emmeans(pesu2, ~dat.type:species_contacted)#
#use emmeans to get the pairwise comparisons (this makes table S7 in supplment and letters in fig 2)
pw_pesu = emmeans(pesu2, pairwise ~ dat.type:species_contacted)

pwc_pesu = pw_pesu$contrasts
pwc_pesu = as.data.frame(pwc_pesu)
names(pwc_pesu)

library(tidyverse)
pwc_pesu2 = pwc_pesu %>%
  arrange(estimate) 

#write.csv(pwc_pesu2,file ="stats output/pesu_bayes_comparisons.csv", row.names = F)

emmip(pesu2, ~species_contacted|dat.type,CIs = TRUE, lty = F)#
stanplot(pesu2, pars = "^b", type = "intervals")


######MYSE#######
myse.dat$species_contacted=relevel(myse.dat$species_contacted, ref="myse")
myse.dat$dat.type=relevel(myse.dat$dat.type, ref="dust")

get_prior(inf|trials(species_contacted_count) ~ dat.type * species_contacted + (1|site),
          data = myse.dat, family = binomial())

prior=set_prior("normal(0,5)", class = "b")

#warning - this Bayesian analysis takes ~1 hour to run on a 2017 macbook air
myse1 <- brm(formula = inf|trials(species_contacted_count) ~ dat.type * species_contacted  + (1|site),
             data = myse.dat, family = binomial(), prior = prior,chains = 4,
             control = list(adapt_delta = 0.97),cores=2,
             save_all_pars = TRUE);summary(myse1, waic = TRUE)

emmeans(myse1, ~dat.type:species_contacted)#
emmip(myse1, ~dat.type:species_contacted,CIs = TRUE)#

pw_myse = emmeans(myse1, pairwise ~ dat.type:species_contacted) # this makes table S6 in supplement and letters in figure 2
pwc_myse = as.data.frame(pw_myse$contrasts)
pwc_myse2 = pwc_myse %>%
  arrange(estimate)

#write.csv(pwc_myse2,file ="stats output/myse_bayes_comparisons.csv", row.names = F)

######MYLU#######
mylu.dat$species_contacted=relevel(mylu.dat$species_contacted, ref="mylu")
mylu.dat$dat.type=relevel(mylu.dat$dat.type, ref="dust")

get_prior(inf|trials(species_contacted_count) ~ dat.type * species_contacted + (1|site),
          data = mylu.dat, family = binomial())

prior=set_prior("normal(0,5)", class = "b")

#warning - this Bayesian analysis takes ~6 hours to run on a 2017 macbook air
mylu1 <- brm(formula = inf|trials(species_contacted_count) ~ dat.type * species_contacted  + (1|site),
             data = mylu.dat, family = binomial(), prior = prior,chains = 4,
             control = list(adapt_delta = 0.97),cores=2,
             save_all_pars = TRUE);summary(mylu1, waic = TRUE)

library(emmeans)
emmeans(mylu1, ~dat.type:species_contacted)#
emmip(mylu1, ~dat.type:species_contacted,CIs = TRUE)#
emmeans(mylu1, pairwise ~ species_contacted)
pw_mylu = emmeans(mylu1, pairwise ~ dat.type:species_contacted) # this makes table S5 in supplement and letters in figure 2

pwc_mylu = as.data.frame(pw_mylu$contrasts)

pwc_mylu2 = pwc_mylu %>%
  arrange(estimate)

#write.csv(pwc_mylu2,file ="stats output/mylu_bayes_comparisons.csv", row.names = F)



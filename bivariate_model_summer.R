#################################
# R script Author: Louis Bliard
# University of Zurich
# Contact 1: louis.bliard@evobio.eu 
# Contact 2: louis.bliard@uzh.ch 
# Last modified 27/03/2023
#################################

####################################################################
# Is flight initiation distance associated with longer-term survival
# in yellow-bellied marmots (Marmota flaviventer)?

# Dan T. Blumstein, McKenna Sanchez, Conner S. Philson, Louis Bliard

# Code for the multivariate model of FID and summer survival
####################################################################

#############################################
#### Clean environment and load packages ####
#############################################

rm(list = ls())  # clean the environment

library(mgcv)
library(car)
library(MASS)
library(psych)
library(data.table)
set.seed(105)   # for reproducibility, random number generator


library(ggplot2)
library(MCMCvis)
library(ggpolypath)
library(ks)
library(stringr)

library(brms)

# Load FID dataset (change path depending on where your file is stored)
fid <- read.csv("~/My Drive/phd/fieldwork/mckenna_project/for_submission/fid_data_osf.csv", row.names=1)
fid <- subset(fid, select=-c(jdate, AD)) # remove unnecessary columns

# Load survival dataset (change path depending on where your file is stored)
survival <- read.csv("~/My Drive/phd/fieldwork/mckenna_project/for_submission/revision/summer_data_osf_23March23.csv") 
survival$summer_survival[survival$summer_survival=="YES"] <- 1 # code survival as 0/1
survival$summer_survival[survival$summer_survival=="NO"] <- 0 # code survival as 0/1
survival$summer_survival <- as.numeric(survival$summer_survival)

colnames(fid)[11] <- "uid"
colnames(survival)[13] <- "uid"

#drop row without fid data and without distance burrow data
fid <- fid[!is.na(fid$FID),]
fid <- fid[!is.na(fid$dist_burrow),]

#drop row without mass data for survival
survival <- survival[!is.na(survival$massjun),]

#drop unused columns
survival <- survival[,-c(8,9)]

fid$SD <- scale(log(fid$SD))[,1]  # log Starting Distance, and scale it

#### get unique ID across the 2 datasets
all_id <- as.data.frame(cbind(unique(c(survival$uid, fid$uid)), 
                              c(1:length(unique(c(survival$uid, fid$uid))))))
all_id$V2 <- as.numeric(all_id$V2)
survival <- merge(survival, all_id, by.x = "uid", by.y = "V1", all.x = T)
fid <- merge(fid, all_id, by.x = "uid", by.y = "V1", all.x = T)
# Unique identifier is called V2



###############################
#### Prepare data for brms ####
###############################

#### prepare data for brms, because needs to be a single dataset (and not 2 as we have for the moment)

add_survival <- setdiff(colnames(fid), colnames(survival))
add_fid <- setdiff(colnames(survival), colnames(fid))

fid$summer_survival <- NA
fid$combined_colony <- NA
fid$predator_index <- NA
fid$massjun <- NA
survival$SD <- NA
survival$FID <- NA
survival$dist_burrow <- NA
survival$trial_number <- NA

#reorder columns
fid <- fid[,c(1,2,3,7,8,9,10,12,4,5,6,11,13,14,15,16)]
survival <- survival[,c(1,2,3,5,9,10,4,12,13,14,15,16,6,7,8,11)]

data_brms <- rbind(fid, survival)

data_brms$sub_fid <- 0
data_brms$sub_fid[!is.na(data_brms$FID)] <- 1
data_brms$sub_survival <- 0
data_brms$sub_survival[!is.na(data_brms$summer_survival)] <- 1

# data_brms is the final dataset
# rows with sub_fid=1 correspond to the fid part of the dataset
# rows with sub_survival=1 correspond to the survival part of the dataset 

data_brms$combined_colony[data_brms$col_area %in% c("picnic", "picnic_lower", "picnic_upper", "picnic_middle", "northpk")] <- "picnic"
data_brms$combined_colony[data_brms$col_area %in% c("mm", "mm_aspen", "mm_maintalus")] <- "mm"
data_brms$combined_colony[data_brms$col_area %in% c("river", "river_rivermound", "river_sagemound", "river_sprucemound", "river_southmound")] <- "river"
data_brms$combined_colony[data_brms$col_area %in% c("avalanche")] <- "avalanche"
data_brms$combined_colony[data_brms$col_area %in% c("bench")] <- "bench"
data_brms$combined_colony[data_brms$col_area %in% c("cliff")] <- "cliff"
data_brms$combined_colony[data_brms$col_area %in% c("stonefield")] <- "stonefield"
data_brms$combined_colony[data_brms$col_area %in% c("horsemound")] <- "horsemound"
data_brms$combined_colony[data_brms$col_area %in% c("boulder")] <- "boulder"
data_brms$combined_colony[data_brms$col_area %in% c("gothictown")] <- "gothictown"
data_brms$combined_colony[data_brms$col_area %in% c("rvannex")] <- "rvannex"


###########################
#### Check sample size ####
###########################

##FID

# number of observations
length(data_brms$FID[data_brms$sub_fid==1 & !is.na(data_brms$dist_burrow)])
# number of individuals
length(unique(data_brms$V2[data_brms$sub_fid==1 & !is.na(data_brms$dist_burrow)]))
# males
length(unique(data_brms$V2[data_brms$sub_fid==1 & data_brms$sex=="M" & !is.na(data_brms$dist_burrow)]))
# females
length(unique(data_brms$V2[data_brms$sub_fid==1 & data_brms$sex=="F" & !is.na(data_brms$dist_burrow)]))


##Summer survival

# number of observations
length(data_brms$summer_survival[data_brms$sub_survival==1 & !is.na(data_brms$massjun)])
# number of individuals
length(unique(data_brms$V2[data_brms$sub_survival==1 & !is.na(data_brms$massjun)]))
# males
length(unique(data_brms$V2[data_brms$sub_survival==1 & data_brms$sex=="M" & !is.na(data_brms$massjun)]))
# females
length(unique(data_brms$V2[data_brms$sub_survival==1 & data_brms$sex=="F" & !is.na(data_brms$massjun)]))

## Model in brms

####################
#### Full model ####
####################

# survival submodel
bf_survival <- bf(summer_survival | subset(sub_survival) ~ predator_index + valley_location + sex + scale(age) + scale(massjun) + (1|p|V2) + (1|year) + (1|combined_colony)) + bernoulli()
# fid submodel
bf_fid <- bf(log(FID) | subset(sub_fid) ~ SD + sex + scale(age) + scale(trial_number) + valley_location + scale(dist_burrow) + (1|p|V2) + (1|year) + (1|combined_colony))

# Run model
fit1 <- brm(bf_survival + bf_fid,
    data = data_brms, chains = 2, cores = 2, seed = 50,
    iter = 10000, thin = 1, warmup = 5000,
    control = list(adapt_delta = 0.9),
    save_pars = save_pars(group = F),
    prior = c(
      set_prior("normal(0,1.5)", class = "b", resp = c("summersurvival", "logFID")),
      set_prior("normal(0,3)", class = "Intercept", resp = c("summersurvival", "logFID")),
      set_prior("normal(0,1)", class = "sd", resp = c("summersurvival", "logFID")),
      set_prior("lkj(1)", class = "cor")),
    backend = "cmdstanr")
summary(fit1)

# posterior predictive check
pp_check(fit1, resp = "logFID")
pp_check(fit1, resp = "summersurvival")


# Posterior samples into a dataframe format
posterior <- as.matrix(fit1)
dat_plot <- as.data.frame(posterior)

# # calculate standardised selection gradient
# selection_gradient <- ((dat_plot$cor_V2__summersurvival_Intercept__logFID_Intercept * 
#                         dat_plot$sd_V2__summersurvival_Intercept *
#                         dat_plot$sd_V2__logFID_Intercept) /
#                         (dat_plot$sd_V2__logFID_Intercept^2)) *
#                      (dat_plot$sd_V2__logFID_Intercept /
#                         plogis(dat_plot$b_summersurvival_Intercept))
# mean(selection_gradient)
# quantile(selection_gradient, c(0.025, 0.975))



# extract posterior distribution for individual random intercepts from the model output
random_effect_surv <- dat_plot[,1067:1582]  # need to make sure the numbers are right if the model is modified
random_effect_fid <- dat_plot[,1583:2098]   # need to make sure the numbers are right if the model is modified

#keep only individuals present in both survival and fid datasets
id_to_keep <- intersect(unique(data_brms$V2[data_brms$sub_fid==1]), unique(data_brms$V2[data_brms$sub_survival==1]))

temporary_storage <- str_extract_all(names(random_effect_fid), "[0-9]+")
id_key <- sapply(temporary_storage,function(x) x[2])

colnames(random_effect_surv) <- id_key
colnames(random_effect_fid) <- id_key

random_effect_surv <- random_effect_surv[names(random_effect_surv) %in% as.character(id_to_keep)]
random_effect_fid <- random_effect_fid[names(random_effect_fid) %in% as.character(id_to_keep)]

# extract posterior distribution for the models intercepts
intercept_surv <- dat_plot[,1]
intercept_fid <- dat_plot[,2]

# predict individual position in the "FID-Survival space" (model intercept + random individual intercept)
re_fid_bt <- matrix(data=NA, nrow = 10000, ncol = dim(random_effect_fid)[2])
re_surv_bt <- matrix(data=NA, nrow = 10000, ncol = dim(random_effect_fid)[2])
for(i in 1:dim(random_effect_fid)[2]) {
  re_fid_bt[,i] <- intercept_fid + random_effect_fid[,i]
  re_surv_bt[,i] <- exp(intercept_surv+random_effect_surv[,i])/(1+exp(intercept_surv+random_effect_surv[,i]))
}
  

plot_data <- as.data.frame(cbind(c(re_fid_bt), c(re_surv_bt)))

# get probability contours for plot
kd <- ks::kde(plot_data, compute.cont=TRUE)
contour_90 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                    z=estimate, levels=cont["10%"])[[1]])
contour_90 <- data.frame(contour_90)

contour_80 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                    z=estimate, levels=cont["20%"])[[1]])
contour_80 <- data.frame(contour_80)

contour_70 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                    z=estimate, levels=cont["30%"])[[1]])
contour_70 <- data.frame(contour_70)

contour_60 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                    z=estimate, levels=cont["40%"])[[1]])
contour_60 <- data.frame(contour_60)

contour_50 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                    z=estimate, levels=cont["50%"])[[1]])
contour_50 <- data.frame(contour_50)

contour_40 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                    z=estimate, levels=cont["60%"])[[1]])
contour_40 <- data.frame(contour_40)

contour_30 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                    z=estimate, levels=cont["70%"])[[1]])
contour_30 <- data.frame(contour_30)

contour_20 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                    z=estimate, levels=cont["80%"])[[1]])
contour_20 <- data.frame(contour_20)

contour_10 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                    z=estimate, levels=cont["90%"])[[1]])
contour_10 <- data.frame(contour_10)


# Make Figure
p <- ggplot(plot_data, aes(x=V1, y=V2))+
  geom_polypath(aes(x, y), data=contour_90, fill="grey80") +
  geom_polypath(aes(x, y), data=contour_80, fill="grey70") +
  geom_polypath(aes(x, y), data=contour_70, fill="grey60") +
  geom_polypath(aes(x, y), data=contour_60, fill="grey50") +
  geom_polypath(aes(x, y), data=contour_50, fill="grey40") +
  geom_polypath(aes(x, y), data=contour_40, fill="grey30") +
  geom_polypath(aes(x, y), data=contour_30, fill="grey20") +
  geom_polypath(aes(x, y), data=contour_20, fill="grey10") +
  geom_polypath(aes(x, y), data=contour_10, fill="grey0") +
  xlab("Log(FID)")+
  ylab("Summer survival probability")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y=element_text(size=14),
        axis.text.x = element_text(size=14),
        plot.title = element_text(hjust = 0.5))



######################
#### Model female ####
######################

# survival submodel
bf_survival <- bf(summer_survival | subset(sub_survival) ~ predator_index + valley_location + scale(age) + scale(massjun) + (1|p|V2) + (1|year) + (1|combined_colony)) + bernoulli()
# fid submodel
bf_fid <- bf(log(FID) | subset(sub_fid) ~ SD + scale(age) + scale(trial_number) + valley_location + scale(dist_burrow) + (1|p|V2) + (1|year) + (1|combined_colony))

# Run model
fit2 <- brm(bf_survival + bf_fid,
            data = subset(data_brms, sex=="F"), chains = 2, cores = 2, seed = 5,
            iter = 10000, thin = 1, warmup = 5000,
            control = list(adapt_delta = 0.985),
            save_pars = save_pars(group = F),
            prior = c(
              set_prior("normal(0,1.5)", class = "b", resp = c("summersurvival", "logFID")),
              set_prior("normal(0,3)", class = "Intercept", resp = c("summersurvival", "logFID")),
              set_prior("normal(0,1)", class = "sd", resp = c("summersurvival", "logFID")),
              set_prior("lkj(1)", class = "cor")),
            backend = "cmdstanr")
summary(fit2)

# posterior predictive check
pp_check(fit2, resp = "logFID")
pp_check(fit2, resp = "summersurvival")


# Posterior samples into a dataframe format
posterior <- as.matrix(fit2)
dat_plot <- as.data.frame(posterior)

# # calculate standardised selection gradient
# selection_gradient <- ((dat_plot$cor_V2__summersurvival_Intercept__logFID_Intercept * 
#                           dat_plot$sd_V2__summersurvival_Intercept *
#                           dat_plot$sd_V2__logFID_Intercept) /
#                          (dat_plot$sd_V2__logFID_Intercept^2)) *
#   (dat_plot$sd_V2__logFID_Intercept /
#      plogis(dat_plot$b_summersurvival_Intercept))
# mean(selection_gradient)
# quantile(selection_gradient, c(0.025, 0.975))

####################
#### Model male ####
####################

# survival submodel
bf_survival <- bf(summer_survival | subset(sub_survival) ~ predator_index + valley_location + scale(age) + scale(massjun) + (1|p|V2) + (1|year) + (1|combined_colony)) + bernoulli()
# fid submodel
bf_fid <- bf(log(FID) | subset(sub_fid) ~ SD + scale(age) + scale(trial_number) + valley_location + scale(dist_burrow) + (1|p|V2) + (1|year) + (1|combined_colony))

# Run model
fit3 <- brm(bf_survival + bf_fid,
            data = subset(data_brms, sex=="M"), chains = 2, cores = 2, seed = 5,
            iter = 10000, thin = 1, warmup = 5000,
            control = list(adapt_delta = 0.9),
            save_pars = save_pars(group = F),
            prior = c(
              set_prior("normal(0,1.5)", class = "b", resp = c("summersurvival", "logFID")),
              set_prior("normal(0,3)", class = "Intercept", resp = c("summersurvival", "logFID")),
              set_prior("normal(0,1)", class = "sd", resp = c("summersurvival", "logFID")),
              set_prior("lkj(1)", class = "cor")),
            backend = "cmdstanr")
summary(fit3)

# posterior predictive check
pp_check(fit3, resp = "logFID")
pp_check(fit3, resp = "summersurvival")



# Posterior samples into a dataframe format
posterior <- as.matrix(fit3)
dat_plot <- as.data.frame(posterior)

# # calculate standardised selection gradient
# selection_gradient <- ((dat_plot$cor_V2__summersurvival_Intercept__logFID_Intercept * 
#                           dat_plot$sd_V2__summersurvival_Intercept *
#                           dat_plot$sd_V2__logFID_Intercept) /
#                          (dat_plot$sd_V2__logFID_Intercept^2)) *
#   (dat_plot$sd_V2__logFID_Intercept /
#      plogis(dat_plot$b_summersurvival_Intercept))
# mean(selection_gradient)
# quantile(selection_gradient, c(0.025, 0.975))


rm(list = ls())

library(knitr)
library(coda)
library(rjags)
library(kableExtra)
library(ggplot2)
library(ggthemes)
library(metafor)
library(meta)
library(tidyverse)

###########taking temp as a dataset with combination#############
# mortality data entered as  TAILOR-PCI, POPULAR GENETICS, TAILOR-PCI WHOLE, ADPT, IAC-PCI, PHARMCLO
Mc<-c(10, 19, 27, 9, 9, 0) #mortality event in non-expose (CONVENTIONAL)
Nc<-c(946, 1246, 2635, 255, 299, 440 ) #total in non-expose
Me<-c(6, 19, 25, 13, 1, 0 ) #mortality event in expose (GENOTYPE-GUIDED)
Ne<-c(903, 1242, 2641, 249, 301, 448 ) #total in expose
MACEc<-c(10, 47, 80, 6, 15, 88) #CV death, MI, stroke event in non-expose (CONVENTIONAL) -MACE
MACEe<-c(10, 36, 62, 9, 3, 54) #CV death, MI, stroke event in expose (GENOTYPE-GUIDED) -MACE
PrimaryOute <- c(35, 63, 113, 34, 81, 71) #Respective Primary Outcome event in expose (GENOTYPE-GUIDED) 
PrimaryOutc <- c(54, 73, 135, 26, 27, 114) #Respective Primary Outcome event in non-expose (CONVENTIONAL) 
CVDeathe <- c(4, 9, 20, 5, 1, 28) #CV Death Outcome event in expose (GENOTYPE-GUIDED) 
CVDeathc <- c(8, 10, 21, 2, 6, 34) #CV Death Primary Outcome event in non-expose (CONVENTIONAL) 
Strokee <- c(2, 8, 32, 2, 1, 5) #Stroke Outcome event in expose (GENOTYPE-GUIDED) 
Strokec <- c(4, 11, 47, 2, 2, 7) #Stroke Outcome event in non-expose (CONVENTIONAL) 
MIe <- c(11, 19, 10, 2, 1, 21) #MI Outcome event in expose (GENOTYPE-GUIDED) 
MIc <- c(14, 26, 12, 2, 7, 47) #MI Outcome event in non-expose (CONVENTIONAL) 
Stente <- c(2, 9, 14, 1, 2, 3) #Stent Thrombosis Outcome event in expose (GENOTYPE-GUIDED) 
Stentc <- c(8, 11, 19, 5, 9, 5) #Stent Thrombosis event in non-expose (CONVENTIONAL) 
Bleede <- c(26, 126, 61, 13, 4, 17) #Bleeding Thrombosis Outcome event in expose (GENOTYPE-GUIDED) 
Bleedc <- c(16, 163, 50, 13, 11, 26) #Bleeding Thrombosis event in non-expose (CONVENTIONAL) 





### Primary outcome difference of the TAILOR-PCI
temp <- data.frame(Mc=Mc,Me=Me,MACEc=MACEc, MACEe=MACEe, PrimaryOutc=PrimaryOutc, PrimaryOute=PrimaryOute, Nc=Nc, Ne=Ne, CVDeathc=CVDeathc, CVDeathe=CVDeathe,Strokec=Strokec, Strokee=Strokee, MIc=MIc, MIe=MIe, Stentc=Stentc, Stente=Stente, Bleedc=Bleedc, Bleede=Bleede)
temp$Study <- c("TAILOR-PCI", "POP-GENTEICS", "TAILOR-PCI Whole", "ADAPT", "IAC-PCI", "PHARMCLO")
temp <- temp[c(19,1:18)]

kable(temp, caption="Outcomes at 1 year") %>%
  kable_styling(bootstrap_options = "striped", full_width = F)
write.csv(temp, file = "data1.csv", row.names = FALSE)

###########################################################################################
# TAILOR-PCI WHOLE data-MACE
TPMACE_R_c <- temp[3,4]
TPMACE_R_e <- temp[3,5]
TPMACE_n_c <- temp[3,8]
TPMACE_n_e <- temp[3,9]

paste("Conirming that hot coded data ", (62/2641) - (80/2635), "is = to data read from file, ", (TPMACE_R_e/TPMACE_n_e) - (TPMACE_R_c/TPMACE_n_c))


set.seed(1234)
#Prior is beta(1,1)
# sampling 100,000 random variables from posterior
post_TPMACE_R_c <- rbeta(100000, TPMACE_R_c + 1, TPMACE_n_c - TPMACE_R_c + 1 )
post_TPMACE_R_e <- rbeta(100000, TPMACE_R_e + 1, TPMACE_n_e - TPMACE_R_e + 1 )

# calculting posterior of differences
post_TPMACE_R_diff <- post_TPMACE_R_e - post_TPMACE_R_c
paste("TAILOR-PCI data alone - Differences in MACE between Genotype-Guided vs Standard Therapy")
quantile(post_TPMACE_R_diff, probs = c(0.025, .5, 0.975))

post_TPMACE_R_RR <- post_TPMACE_R_e / post_TPMACE_R_c
quantile(post_TPMACE_R_RR, probs = c(0.025, .5, 0.975))

# probabilities >0 and >1
paste("TAILOR-PCI data alone - Probability Genotype-Guided vs Standard Therapy = ", sum(post_TPMACE_R_diff*100 <0)/100000)
paste("TAILOR-PCI data alone - Probability Genotype-Guided better than Standard Therapy by >-1% = ", sum(-1>post_TPMACE_R_diff*100)/100000)

# given large sample sizes, can verify answers with normal approximation
paste("With normal approximation, TAILOR-PCI data alone - Probability Genotype-Guided better than Standard Therapy =", round(pnorm(0, mean(post_TPMACE_R_diff*100), sd(post_TPMACE_R_diff*100)),3))
paste("With normal approximation, TAILOR-PCI data alone - Probability Genotype-Guided better than Standard Therapy by >-1% =", round(pnorm(-1, mean(post_TPMACE_R_diff*100), sd(post_TPMACE_R_diff*100)),3))

TPMACE_diff_df <- data.frame(post_TPMACE_R_diff)
TPMACE_diff_df$post_TPMACE_R_diff <- TPMACE_diff_df$post_TPMACE_R_diff*100

# plot showing good normal approximation to binomial histogram
ggplot(TPMACE_diff_df, aes(x= post_TPMACE_R_diff)) + 
  geom_histogram(aes(y=..density..), bins=100, colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") +
  scale_x_continuous(name = "MACE Outcome Difference (Genotype-Guided - Standard Therapy)") +
  scale_y_continuous(name = "Density") +
  ggtitle("TAILOR-PCI Probability density for risk difference in  MACE\n (with superimposed Gaussian kernal density estimate )") +
  theme_economist()

### TAILOR-PCI MACE outcome graph for the TREATMENT APPROACH
#### Figure 1c
ggplot(data.frame(x = c(-2.5, 2)), aes(x = x)) +
  stat_function(fun = dnorm, args = list(mean(post_TPMACE_R_diff*100), sd(post_TPMACE_R_diff*100)), colour = "deeppink") +
  scale_x_continuous(name = "MACE Outcome Difference (Genotype-Guided - Standard Therapy) -> Genotype-Guided Better") +
  scale_y_continuous(name = "Density") +
  ggtitle("TAILOR-PCI Probability Density Function \n (MACE Outcome Risk Difference )") +
  geom_vline(xintercept=mean(post_TPMACE_R_diff*100)) +
  annotate("text", label = "Black vertical line = mean outcome \n difference (0.8%) decreased with Genotype-Guided Therapy", x = 0.65, y = .60, color = "black") +
  annotate("text", label = "Grey AUC = probability (28%) \n Genotype-Guided > Standard Therapy outcome by > 1%", x =0.85, y = .35, color = "black") +
  annotate("text", label = "Grey + yellow AUC = probability (94%) \n Genotype-Guided > Standard Therapy  outcome", x = 1.25, y = 0.15, color = "black") +
  theme_economist() +
  stat_function(fun = dnorm, args = list(mean(post_TPMACE_R_diff*100), sd(post_TPMACE_R_diff*100)), xlim = c(0,-3), geom = "area", alpha = 0.2) +
  stat_function(fun = dnorm, args = list(mean(post_TPMACE_R_diff*100), sd(post_TPMACE_R_diff*100)), xlim = c(0,-1), geom = "area", alpha = 0.2, fill = "yellow") 

#################################################################################
#############################REPEAT SAME FOR ALL STUDY OUTCOMES##################

########################NOW UPDATING THE PRIORS FOR MACE OUTCOME#######################
#####THIS SECTION IS FOR GENERATING AND UPDATING PRIORS########
#################### MAIN STUDY OUTCOMES######################
#   MACE= CV DEATH + STROKE + MI                              #
#   CV DEATH                                                  #
#   MI                                                        #
#   STROKE                                                    #
#   STENT THROMBOSIS                                          #
#   BLEEDING (MAJOR OR MINOR BLEEDING)                        #
########################################################################################


post.normal.mean <- function(data.mean, data.var, prior.mean, prior.var)
{
  ####################################################################
  #  R function for Bayesian analysis of normal mean, variance known #
  #  Parameters included are:                                        #
  #                                                                  #
  #  Inputs:                                                         #
  #                                                                  #
  #   x = vector of data                                             #
  #   prior.mean = prior mean                                        #
  #   prior.var  = prior variance                                    #
  #   data.var   = assumed known variance of data                    #
  #                                                                  #
  #  Outputs:                                                        #
  #                                                                  #
  #   post.mean = posterior mean                                     #
  #   post.var  = posterior variance                                 #
  #                                                                  #
  ####################################################################
  post.mean.numerator <- prior.mean/prior.var + data.mean/data.var
  post.mean.denominator <- 1/prior.var + 1/data.var
  post.mean <-  post.mean.numerator/post.mean.denominator
  post.var <- (1/(1/prior.var + 1/data.var))
  a <- "Post mean = "
  b <- "Post Var = "
  c <- "Post SD = "
  cat(a, post.mean, ",", b, post.var, ",", c, sqrt(post.var), "\n" )
  newlist <- list(post.mean, post.var, sqrt(post.var))
  return(newlist)
}

#### Random effects model of previous studies #####
temp1<-temp[c(2,4,5,6,3,1),]
prior <- temp1[-c(5,6),]

############################MACE###########################################
RD_re <- rma(ai=MACEe, n1i= Ne, ci=MACEc, n2i=Nc, data=prior, measure="RD",
             slab=paste(Study), method="REML")

prior.diff.prim.mean <-RD_re$beta
prior.diff.prim.sd<- RD_re$se

prior.diff.prim.ARR <-rnorm(100000,mean=prior.diff.prim.mean , sd=prior.diff.prim.sd)
hist(prior.diff.prim.ARR)

quantile(prior.diff.prim.ARR, probs = c(0.025, .5, 0.975))

# same procedure for the TAILOR-PCI data (likelihood)
like.e.prim.mean <-temp[3,5]/temp[3,9] # like Genotype-guided (experimental) for MACE outcome
like.c.prim.mean <-temp[3,4]/temp[3,8]
like.diff.prim.mean <- like.e.prim.mean- like.c.prim.mean
like.e.prim.sd<-sqrt(like.e.prim.mean * (1-like.e.prim.mean) / temp[3,9])
like.c.prim.sd<-sqrt(like.c.prim.mean * (1-like.c.prim.mean) /temp[3,8])
like.diff.prim.sd<-sqrt(like.e.prim.sd^2 + like.c.prim.sd^2)
# TAILOR-PCI MACE outcome + random effects model of previous studies combined with normal conjugacy
data.mean <- like.diff.prim.mean # TAILOR-PCI data
data.var <- like.diff.prim.sd^2
prior.mean <- RD_re$beta  # rma estimates
prior.var <- RD_re$se^2 
post.primary <- post.normal.mean(data.mean, data.var, prior.mean, prior.var)
post.primary

df <- data.frame(100*rbind(c(prior.diff.prim.mean, prior.diff.prim.sd), c(like.diff.prim.mean, like.diff.prim.sd), c(post.primary[[1]], post.primary[[3]])))
df$type <- c("Prior", "TAILOR-PCI", "Posterior")
df <- df %>%
  rename(mean = X1, se = X2)

ggplot(data.frame(x = c(-6, 2)), aes(x = x)) +
  stat_function(fun = dnorm, args = list(post.primary[[1]]*100, post.primary[[3]]*100), colour = "deeppink") +
  stat_function(fun = dnorm, args = list(RD_re$beta*100, RD_re$se*100), colour = "blue") +
  stat_function(fun = dnorm, args = list(like.diff.prim.mean*100, like.diff.prim.sd*100), colour = "dark green") +
  ggtitle("Triplot of MACE Outcome \n(Prior, Likelihood (TAILOR-PCI) and Combined Posterior Distribution)") +
  xlab ("
        MACE Outcome Difference (Genotype-Guided - Standard Therapy)") +
  ylab ("Density
        
        ") +
  annotate("text", label = "Blue Line \n Prior Information \n (POPular Genetics, ADAPT, IAC-PCI \n& PHARMCLO) ", x = -4, y = .35, color = "black") +
  annotate("text", label = "Green Line \n Likelihood from TAILOR-PCI ", x = 1.1, y = .11, color = "black") + 
  annotate("text", label = "Red Line Posterior Distribution (Combined Data) \n Probability > 0 = 98% (yellow + grey) \n Probability > -1% = 43% (grey) ", x = -2.7, y = .65, color = "black") +
  theme_economist() +
  stat_function(fun = dnorm, args = list(post.primary[[1]]*100, post.primary[[3]]*100), 
                xlim = c(0,-3), geom = "area", alpha = 0.2) +
  stat_function(fun = dnorm, args = list(post.primary[[1]]*100, post.primary[[3]]*100), 
                xlim = c(0,-1), geom = "area", alpha = 0.2, fill = "yellow") 



paste(" Combined posterior mean of TAILOR-PCI data and prior studies while accounting for betweeen study variation = ", round(post.primary[[1]]*100,2), " (", round(post.primary[[3]]*100, 2), ") fewer outcomes with Genotype-Guided Therapy")
paste("95% CrI", round(post.primary[[1]]*100 - 1.96 * post.primary[[3]]*100,2), " - ", round(post.primary[[1]]*100 + 1.96 * post.primary[[3]]*100,2))  
paste(" Probability that MACE outcome Genotype-Guided > Standard Therapy")
pnorm(0, post.primary[[1]], post.primary[[3]])
paste(" Probability that MACE outcome Genotype-Guided > Standard Therapy risk by at least 1 / 100 treated")
pnorm(-.01, post.primary[[1]], post.primary[[3]])

############################MACE-Risk Ratio###########################################
##https://urldefense.proofpoint.com/v2/url?u=https-3A__www.rdocumentation.org_packages_metafor_versions_2.4-2D0_topics_rma.mh&d=DwIGAg&c=o3PTkfaYAd6-No7SurnLt5qpge1aKYwPQyBFS7c8AA0&r=M5XeGtAeas0lTps1_ivk4JpK173cg_DJZHyoH19U9kQ&m=yRSo6bjIMez-16btvUp1ADyC3btqW190oOM3j6Wwa6o&s=bHQZqa7_gV268u6DjhaqjuXeCDAKg8mUWqReLK4CSZY&e= 
##https://www.medcalc.org/manual/relativerisk_oddsratio.php

RR_re<-rma.mh(ai=MACEe, ci=MACEc, n1i=Ne, n2i=Nc, measure="RR", data=prior, slab=paste(Study))

prior.primary.logRR.mean <-RR_re$beta
prior.primary.logRR.sd<- RR_re$se
prior.primary.logRR.var <- prior.primary.logRR.sd^2


prior.primary.logRR <-rnorm(100000,mean=prior.primary.logRR.mean, sd=prior.primary.logRR.sd)
hist(prior.primary.logRR)

prior.primary.RR <- exp(prior.primary.logRR)
hist(prior.primary.RR)
prior.diff.primary.mean<-mean(prior.primary.RR)
prior.diff.primary.sd<-sd(prior.primary.RR)
quantile(prior.primary.RR, probs = c(0.025, .5, 0.975))

# same procedure for the TAILOR-PCI data (likelihood)
data.e.primary.mean <-temp[3,5]/temp[3,9] # like Genotype-guided (experimental) for MACE outcome
data.c.primary.mean <-temp[3,4]/temp[3,8]
data.RR.primary.mean <-data.e.primary.mean/data.c.primary.mean;

data.logRR.primary.mean <- log(data.RR.primary.mean)
data.logRR.primary.var <- 1/temp[3,5]+1/temp[3,4]-1/temp[3,9]-1/temp[3,8]
data.logRR.primary.sd<- sqrt(data.logRR.primary.var)
data.RR.primary.lowerCI<-exp(data.logRR.primary.mean-(1.96*data.logRR.primary.sd))
data.RR.primary.UpperCI<-exp(data.logRR.primary.mean+(1.96*data.logRR.primary.sd))
data.RR.primary.sd<-(data.RR.primary.mean-data.RR.primary.lowerCI)/1.96

data.RR.primary.mean
data.RR.primary.lowerCI
data.RR.primary.UpperCI

post.primary.logRR <- post.normal.mean(data.logRR.primary.mean, data.logRR.primary.var, prior.primary.logRR.mean, prior.primary.logRR.var)
post.primary.logRR

post.primary.logRR.mean <- post.primary.logRR[[1]]
post.primary.logRR.var <- post.primary.logRR[[2]]
post.primary.logRR.sd<- post.primary.logRR[[3]]

post.primary.logRR <-rnorm(100000,mean=post.primary.logRR.mean, sd=post.primary.logRR.sd)
hist(post.primary.logRR)

post.primary.RR <- exp(post.primary.logRR)
hist(post.primary.RR)
post.primary.RR.mean<-mean (post.primary.RR)
post.primary.RR.sd<-sd(post.primary.RR)
quantile(post.primary.RR, probs = c(0.025, .5, 0.975))


df <- data.frame(100*rbind(c(prior.diff.primary.mean, prior.diff.primary.sd), c(post.primary.RR.mean, post.primary.RR.sd), c(post.primary.RR.mean, post.primary.RR.sd)))
df$type <- c("Prior", "TAILOR-PCI", "Posterior")
df <- df %>%
rename(mean = X1, se = X2)
ggplot(data.frame(x = c(0.5, 1.1)), aes(x = x)) +
stat_function(fun = dnorm, args = list(prior.diff.primary.mean, prior.diff.primary.sd), colour = "blue") +
stat_function(fun = dnorm, args = list(data.RR.primary.mean, data.RR.primary.sd), colour = "dark green") +
  stat_function(fun = dnorm, args = list(post.primary.RR.mean, post.primary.RR.sd), colour = "deeppink")+
   theme_classic()+
  stat_function(fun = dnorm, args = list(post.primary.RR.mean, post.primary.RR.sd), 
                                 xlim = c(1.0,0.5), geom = "area", alpha = 0.2)

paste(" Probability that MACE outcome Genotype-Guided > Standard Therapy with RR<1")
pnorm(1, post.primary.RR.mean, post.primary.RR.sd)
paste(" Probability that MACE outcome Genotype-Guided > Standard Therapy with RR<0.9")
pnorm(0.9, post.primary.RR.mean, post.primary.RR.sd)
paste(" Probability that MACE outcome Genotype-Guided > Standard Therapy with RR<0.8")
pnorm(0.8, post.primary.RR.mean, post.primary.RR.sd)



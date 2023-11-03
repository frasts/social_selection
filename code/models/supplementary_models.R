###############################################

# R script for Supplementary Material analyses
# Associated to "Social interactions generate complex selection patterns in virtual worlds"
# Santostefano, Francesca; Fraser Franco, Maxime; Montiglio, Pierre-Olivier
# 01-Nov-2023
# See associated Readme file for metadata description

###############################################

# This R script is organized as follows:

# 1. Dependencies and setup
# 2. Mixed models points  
##  a. Model 1 (linear gradients)
##  b. Model 2
##  c. Model 3
##  d. Model 4 (quadratic and correlational)
# 3. Selection gradients 
# 4. Selection differentials
# 5. Contributions to selection
# 6. Model selection

#######################
# 1. Dependencies and setup #
######################

#"R version 4.2.3 (2023-03-15 ucrt)"

#Set working directory to you own path where the folder "social_selection" is
setwd("/path/to/folder/social_selection")

# List of packages needed
packages_list <- c("data.table","lme4","arm","dplyr","jtools" ,"sjPlot" , "performance", "MuMIn")

# Load packages - install them if not already present
lapply(packages_list, function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
})

#read in data file
df <- fread("data/df_social_selection.csv")


#######################
# 2.  Mixed models points  #
######################

##  a. Model 1 (linear gradients)

system.time(mod_bp_1.1<- glmer(cbind(bloodpoints, 32000-bloodpoints)  ~  Zgen+ Zcrouching + Zchase + Zunhook +
                                 Zsocial_gen+ Zsocial_crouching + Zsocial_chase + Zsocial_unhook+
                                 Zgame_duration + 
                                 (1|map_code), data = df, family = binomial (link = "logit"),
                               control = glmerControl(optimizer = "nloptwrap", nAGQ0initStep = TRUE) ))

ss <- getME(mod_bp_1.1,c("theta","fixef"))
mod_bp_1.1.upd <- update(mod_bp_1.1,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))

saveRDS(mod_bp_1.1.upd, file = "outputs/models_outputs/mod_bp_1.1.upd.rds")
#mod_bp_1.1.upd<- readRDS("outputs/models_outputs/mod_bp_1.1.upd.rds")

#models summary
summ(mod_bp_1.1.upd) #log odds
summ(mod_bp_1.1.upd, exp = TRUE, r.squared = FALSE) #odds ratio
coefs_bp_1.1.upd<-summ(mod_bp_1.1.upd)$coeftable

#function to extract probabilities
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

logit2prob(coefs_bp_1.1.upd[,1]) #probabilities


#sim  posterior simulations o beta from a glm object
nsim<-1000
bsim<-arm::sim(mod_bp_1.1.upd, n.sim=nsim)

#credibility interval for fixef
apply(bsim@fixef, 2, quantile, prob=c(0.025, 0.975))

#model coefficients
summary(mod_bp_1.1.upd)$coefficients

#function to copy model output 

table <- function(model) {
  a<- apply(bsim@fixef, 2, quantile, prob=c(0.025, 0.975))
  ta<-t(a)
  ta<-as.data.frame(ta)
  b<- summary(model)$coefficients
  b<-as.data.frame(b)
  b<- b[ -c(2,3) ]
  tab<-cbind(b,ta)  
  tab<-round(tab,2)
  
  tab
  
}

#table for fixed effects
fixefs_mod_bp_1.1.upd<-(table(mod_bp_1.1.upd)[,-2])
names(fixefs_mod_bp_1.1.upd)<-cbind("Selection gradient", "lower 2.5% CI",  "upper 97.5% CI")
fixefs_mod_bp_1.1.upd <- data.frame("Fixed effects" = c("Intercept", "Focal resource acquisition ", "Focal hiding", "Focal defense", "Focal rescue", "Social resource acquisition ", "Social hiding", "Social defense", "Social rescue",  "Game duration"), fixefs_mod_bp_1.1.upd )

#table for random effects
tab_map_code_mod_bp_1.1.upd<-round(t(as.data.frame(quantile(apply(bsim@ranef$map_code[ , , 1],1,var), prob=c(0.5, 0.025, 0.975)))),2)
ranefs_mod_bp_1.1.upd<-as.data.frame(tab_map_code_mod_bp_1.1.upd)
names(ranefs_mod_bp_1.1.upd)<-cbind("Variance", "lower 2.5% CI",  "upper 97.5% CI")
ranefs_mod_bp_1.1.upd <- data.frame("Random effects" = "Map", ranefs_mod_bp_1.1.upd)

#TABLE S4 (model 1 output)
        tab_dfs(list(fixefs_mod_bp_1.1.upd,ranefs_mod_bp_1.1.upd),
        titles = c("Model 1","Model 1"),
        col.header = col,
        file = "outputs/tables/table_mod_bp_1.1.upd.doc")


##  b. Model 2

system.time(mod_bp_2.1<- glmer (cbind(bloodpoints, 32000-bloodpoints) ~   Zgen+ Zcrouching + Zchase + Zunhook +  Zsocial_gen+ Zsocial_crouching + Zsocial_chase + Zsocial_unhook+
                                  Zgen*Zcrouching + Zgen*Zchase + Zgen*Zunhook +
                                  Zcrouching*Zchase + Zcrouching*Zunhook + Zchase*Zunhook +
                                  Zsocial_gen*Zsocial_crouching + Zsocial_gen*Zsocial_chase + Zsocial_gen*Zsocial_unhook +
                                  Zsocial_crouching*Zsocial_chase + Zsocial_crouching*Zsocial_unhook + Zsocial_chase*Zsocial_unhook +
                                  Zgen*Zsocial_gen + Zcrouching*Zsocial_crouching +  Zchase*Zsocial_chase +  Zunhook*Zsocial_unhook +
                                  Zgame_duration +
                                  (1|map_code),
                                data = df, family = binomial (link = "logit"),
                                control = glmerControl(optimizer = "nloptwrap", nAGQ0initStep = TRUE)) )

ss <- getME(mod_bp_2.1,c("theta","fixef"))
mod_bp_2.1.upd <- update(mod_bp_2.1,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))

saveRDS(mod_bp_2.1.upd, file = "outputs/models_outputs/mod_bp_2.1.upd.rds")
#mod_bp_2.1.upd<- readRDS("outputs/models_outputs/mod_bp_2.1.upd.rds")


##  c. Model 3

system.time(mod_bp_3.1<- glmer(cbind(bloodpoints, 32000-bloodpoints)  ~   Zgen + Zcrouching + Zchase+ Zunhook  + 
                                 Zsocial_gen+ Zsocial_crouching + Zsocial_chase + Zsocial_unhook+
                                 + I(Zgen^2) + I(Zcrouching^2) +  I(Zchase^2) +  I(Zunhook^2) +
                                 + I(Zsocial_gen^2) + I(Zsocial_crouching^2) +  I(Zsocial_chase^2) +  I(Zsocial_unhook^2) +
                                 Zgame_duration + 
                                 (1|map_code),
                               data = df, family = binomial (link = "logit"),
                               control = glmerControl(optimizer = "nloptwrap", nAGQ0initStep = TRUE)) )

ss <- getME(mod_bp_3.1,c("theta","fixef"))
mod_bp_3.1.upd <- update(mod_bp_3.1,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))

saveRDS(mod_bp_3.1.upd, file = "outputs/models_outputs/mod_bp_3.1.upd.rds")
#mod_bp_3.1.upd<- readRDS("outputs/models_outputs/mod_bp_3.1.upd.rds")


##  d. Model 4

system.time(mod_bp_4.1<- glmer(cbind(bloodpoints, 32000-bloodpoints) ~   Zgen + Zcrouching + Zchase+ Zunhook  + 
                                 Zsocial_gen+ Zsocial_crouching + Zsocial_chase + Zsocial_unhook+
                                 + I(Zgen^2) + I(Zcrouching^2) +  I(Zchase^2) +  I(Zunhook^2) +
                                 + I(Zsocial_gen^2) + I(Zsocial_crouching^2) +  I(Zsocial_chase^2) +  I(Zsocial_unhook^2) +
                                 Zgen*Zcrouching + Zgen*Zchase + Zgen*Zunhook +
                                 Zcrouching*Zchase + Zcrouching*Zunhook + Zchase*Zunhook +
                                 Zsocial_gen*Zsocial_crouching + Zsocial_gen*Zsocial_chase + Zsocial_gen*Zsocial_unhook +
                                 Zsocial_crouching*Zsocial_chase + Zsocial_crouching*Zsocial_unhook + Zsocial_chase*Zsocial_unhook +
                                 Zgen*Zsocial_gen +   Zchase*Zsocial_chase + Zcrouching*Zsocial_crouching + Zunhook*Zsocial_unhook +
                                 Zgame_duration + 
                                 (1|map_code),
                               data = df, family = binomial (link = "logit"),
                               control = glmerControl(optimizer = "nloptwrap", nAGQ0initStep = TRUE)) )

ss <- getME(mod_bp_4.1,c("theta","fixef"))
mod_bp_4.1.upd <- update(mod_bp_4.1,start=ss,control=glmerControl( optCtrl=list(maxfun=2e4)))

saveRDS(mod_bp_4.1.upd, file = "outputs/models_outputs/mod_bp_4.1.upd.rds")
#mod_bp_4.1.upd<- readRDS("outputs/models_outputs/mod_bp_4.1.upd.rds")

#model summaries
summ(mod_bp_4.1.upd) #log odds, to use for gradients
summ(mod_bp_4.1.upd, exp = TRUE) #odds ratio
coefs_bp_4.1.upd<-summ(mod_bp_4.1.upd)$coeftable
logit2prob(coefs_bp_4.1.upd[,1]) #probabilities

#sim  posterior simulations of beta from a glm object
nsim<-1000
bsim<-arm::sim(mod_bp_4.1.upd, n.sim=nsim)

#credibility interval for fixef
apply(bsim@fixef, 2, quantile, prob=c(0.025, 0.975))
#model coefficients
summary(mod_bp_4.1.upd)$coefficients

#table for fixed effects
fixefs_mod_bp_4.1.upd<-(table(mod_bp_4.1.upd)[,-2])
fixefs_mod_bp_4.1.upd <- data.frame("Fixed effects" = c("Intercept", "Focal resource acquisition", "Focal hiding", "Focal defense", "Focal rescue", "Social resource acquisition", "Social hiding", "Social defense", "Social rescue", 
                                                        "Quadratic Focal resource acquisition", "Quadratic Focal hiding", "Quadratic Focal defense", "Quadratic Focal rescue", "Quadratic Social resource acquisition ", "Quadratic Social hiding", "Quadratic Social defense", "Quadratic Social rescue",
                                                        "Game duration", "Focal resource acquisition*hiding", "Focal resource acquisition*defense", "Focal resource acquisition*rescue", "Focal hiding*defense", "Focal hiding*rescue", "Focal defense*rescue",
                                                        "Social resource acquisition*hiding", "Social resource acquisition*defense", "Social resource acquisition*rescue", "Social hiding*defense", "Social hiding*rescue", "Social defense*rescue",
                                                        "Focal*Social resource acquisition", "Focal*Social hiding", "Focal*Social defense", "Focal*Social rescue"
), fixefs_mod_bp_4.1.upd )


#quadratic terms need recalculation too and also the respective CI (2x multiply coefficients)

#function to get quadratic terms coefficients and CI for those
quad <- function(model) {
  aa<- as.data.frame(bsim@fixef)
  aa<-aa[c(10:17)] #this may change if we change fixed effects 
  aa2<-aa*2                 
  aaa2<-apply(aa2, 2, quantile, prob=c(0.025, 0.975))
  taaa2<-t(aaa2)
  taaa2<-as.data.frame(taaa2)
  
  bb<-summary(model)$coefficients
  bb<-as.data.frame(bb)
  bb<-bb$Estimate[c(10:17)]
  Estimate<-bb*2
  
  tabq<-cbind(Estimate,taaa2) 
  tabq<-round(tabq,2)
  tabq
}

quad_mod_bp_4.1.upd<-quad(mod_bp_4.1.upd)
quad_mod_bp_4.1.upd <- data.frame("Fixed effects" = c("Quadratic Focal resource acquisition", "Quadratic Focal hiding", "Quadratic Focal defense", "Quadratic Focal rescue", "Quadratic Social resource acquisition ", "Quadratic Social hiding", "Quadratic Social defense", "Quadratic Social rescue"), quad_mod_bp_4.1.upd )

#replacing with the quadratic coefficients x 2
fixefs_mod_bp_4.1.upd[match(quad_mod_bp_4.1.upd$Fixed.effects, fixefs_mod_bp_4.1.upd$Fixed.effects), ] <- quad_mod_bp_4.1.upd
names(fixefs_mod_bp_4.1.upd)<-cbind("Fixed Effects","Selection gradient", "lower 2.5% CI",  "upper 97.5% CI")

#table for random effects
tab_map_code_mod_bp_4.1.upd<-round(t(as.data.frame(quantile(apply(bsim@ranef$map_code[ , , 1],1,var), prob=c(0.5, 0.025, 0.975)))),2)
ranefs_mod_bp_4.1.upd<-as.data.frame(tab_map_code_mod_bp_4.1.upd)
names(ranefs_mod_bp_4.1.upd)<-cbind("Variance", "lower 2.5% CI",  "upper 97.5% CI")
ranefs_mod_bp_4.1.upd <- data.frame("Random effects" = ("Map"), ranefs_mod_bp_4.1.upd)

#TABLE S4 (model 4 output)
tab_dfs(list(fixefs_mod_bp_4.1.upd,ranefs_mod_bp_4.1.upd),
        titles = c("Model 1","Model 1"),
        col.header = col,
        file = "outputs/tables/table_mod_bp_4.1.upd.doc")


#######################
# 3. Selection gradients #
######################

coefs_bp_1.1.upd<-summ(mod_bp_1.1.upd)$coeftable
coefs_bp_2.1.upd<-summ(mod_bp_2.1.upd)$coeftable
coefs_bp_3.1.upd<-summ(mod_bp_3.1.upd)$coeftable
coefs_bp_4.1.upd<-summ(mod_bp_4.1.upd)$coeftable


#extracting log odds linear coefficients (betas - selection gradients)
coefs_bp_1.1.upd[,1]
B_bn_gen<-coefs_bp_1.1.upd[2,1]
B_bn_crouching<-coefs_bp_1.1.upd[3,1]
B_bn_chase<-coefs_bp_1.1.upd[4,1]
B_bn_unhook<-coefs_bp_1.1.upd[5,1]
B_bs_gen<-coefs_bp_1.1.upd[6,1]
B_bs_crouching<-coefs_bp_1.1.upd[7,1]
B_bs_chase<-coefs_bp_1.1.upd[8,1]
B_bs_unhook<-coefs_bp_1.1.upd[9,1]

#######################
# 4. Selection differentials #
######################

#EQ 11 from Wolf et al. 1999 paper: S=P*bN+CI*bS

#organizing vectors of gradients

B_bn = matrix(
  c( B_bn_gen,
     B_bn_crouching,
     B_bn_chase, 
     B_bn_unhook
  ),
  nrow=4,
  ncol=1,
  byrow = TRUE) 

B_bn

B_bs = matrix(
  c( B_bs_gen, 
     B_bs_crouching,
     B_bs_chase, 
     B_bs_unhook
  ),
  nrow=4,
  ncol=1,
  byrow = TRUE) 

B_bs


# (variance is 1 for all traits in the full dataset)
V=1  

#matrix of VCV for focal traits, P  
focaltraits <- select(df, Zgen, Zcrouching, Zchase, Zunhook)
P<-cov(focaltraits)
P

#CI is the covariance matrix of focal and social traits (interactant phenotypic covariances)
#The elements of the matrix CI are the covariance of the focal individual phenotype (the rows) with the mean social environment (the columns) experienced by that individual. 

#        Zsocial_gen    Zsocial_crouching        Zsocial_chase Zsocial_unhook
#Zgen     Cg                   Cgcr               Cgch           Cgu
#Zcrouching Ccrg               Ccr                 Ccrch         Ccru
#Zchase   Cchg                Cchcr               Cch           Cchu
#Zunhook  Cug                 Cucr                Cuch           Cu

Cg <-cov(df$Zgen, df$Zsocial_gen)
Ccr <-cov(df$Zcrouching, df$Zsocial_crouching)
Cch <-cov(df$Zchase, df$Zsocial_chase)
Cu <-cov(df$Zunhook, df$Zsocial_unhook)
Cgcr <-cov(df$Zgen, df$Zsocial_crouching)
Cgch <-cov(df$Zgen, df$Zsocial_chase)
Cgu <-cov(df$Zgen, df$Zsocial_unhook)
Ccrg <-cov(df$Zcrouching, df$Zsocial_gen)
Ccrch <-cov(df$Zcrouching, df$Zsocial_chase)
Ccru <-cov(df$Zcrouching, df$Zsocial_unhook)
Cchg <-cov(df$Zchase, df$Zsocial_gen)
Cchcr <-cov(df$Zchase, df$Zsocial_crouching)
Cchu <-cov(df$Zchase, df$Zsocial_unhook)
Cug <-cov(df$Zunhook, df$Zsocial_gen)
Cucr <-cov(df$Zunhook, df$Zsocial_crouching)
Cuch <-cov(df$Zunhook, df$Zsocial_chase)

#matrix for interacting phen traits, CI
CI  = matrix(
  c( Cg ,                Cgcr,               Cgch,           Cgu,
     Ccrg ,              Ccr ,               Ccrch ,        Ccru,
     Cchg,                Cchcr ,              Cch,           Cchu,
     Cug,                Cucr,                Cuch,           Cu
     
  ), # the data elements
  nrow=4,              # number of rows
  ncol=4,              # number of columns
  byrow = TRUE)        # fill matrix by rows

CI                     # print the matrix


#equation for the selection differentials S=P*bN+CI*bS

B_S=P%*%B_bn+CI%*%B_bs



# TABLE S2 with gradients and selection differentials

BP_gradients <- data.frame(Behaviour=as.character(c("Resource acquisition", "Hiding", "Defense", "Rescue")),
                           Interaction_type_predicted= as.character(c("Cooperation +/+", "Selfish +/-", "Altruistic -/+", "Altruistic -/+")),
                           ßN = B_bn,
                           ßS = B_bs,
                           Si = B_S)
BP_gradients  

tab_df(BP_gradients,
        titles =   "B) Selection gradients Points",
        col.header = col,
        file = "outputs/tables/table_gradients_points.doc")


#######################
# 5. Contributions to selection
######################

#natural selection contrib
P_B_bn<-P%*%B_bn
P_B_bn

#social sel contrib
CI_B_bs<-CI%*%B_bs
CI_B_bs

#Si
B_S

P
B_bn
CI
B_bs

#direct and indirect
#nat sel direct
D_B_bn_gen<-P[1,1]%*%B_bn[1,1]
D_B_bn_crouch<-P[2,2]%*%B_bn[2,1]
D_B_bn_chase<-P[3,3]%*%B_bn[3,1]
D_B_bn_unhook<-P[4,4]%*%B_bn[4,1]

#nat sel indirect
I_B_bn_gen<-P[1,2]%*%B_bn[2,1]+P[1,3]%*%B_bn[3,1]+P[1,4]%*%B_bn[4,1]
I_B_bn_crouch<-P[2,1]%*%B_bn[1,1]+P[2,3]%*%B_bn[3,1]+P[2,4]%*%B_bn[4,1]
I_B_bn_chase<-P[3,1]%*%B_bn[1,1]+P[3,2]%*%B_bn[2,1]+P[3,4]%*%B_bn[4,1]
I_B_bn_unhook<-P[4,1]%*%B_bn[1,1]+P[4,2]%*%B_bn[2,1]+P[4,3]%*%B_bn[3,1]


#soc sel direct
D_B_bs_gen<-CI[1,1]%*%B_bs[1,1]
D_B_bs_crouch<-CI[2,2]%*%B_bs[2,1]
D_B_bs_chase<-CI[3,3]%*%B_bs[3,1]
D_B_bs_unhook<-CI[4,4]%*%B_bs[4,1]

#soc sel indirect
I_B_bs_gen<-CI[1,2]%*%B_bs[2,1]+CI[1,3]%*%B_bs[3,1]+CI[1,4]%*%B_bs[4,1]
I_B_bs_crouch<-CI[2,1]%*%B_bs[1,1]+CI[2,3]%*%B_bs[3,1]+CI[2,4]%*%B_bs[4,1]
I_B_bs_chase<-CI[3,1]%*%B_bs[1,1]+CI[3,2]%*%B_bs[2,1]+CI[3,4]%*%B_bs[4,1]
I_B_bs_unhook<-CI[4,1]%*%B_bs[1,1]+CI[4,2]%*%B_bs[2,1]+CI[4,3]%*%B_bs[3,1]


#TABLE S5 - contributions to selection (also used for fig. S5 )

B_gen_contrib<-c(D_B_bn_gen,I_B_bn_gen,D_B_bs_gen,I_B_bs_gen )
B_crouch_contrib<-c(D_B_bn_crouch,I_B_bn_crouch,D_B_bs_crouch,I_B_bs_crouch )
B_chase_contrib<-c(D_B_bn_chase,I_B_bn_chase,D_B_bs_chase,I_B_bs_chase )
B_unhook_contrib<-c(D_B_bn_unhook,I_B_bn_unhook,D_B_bs_unhook,I_B_bs_unhook )


B_contributions <- data.frame(Behaviour=as.character(c("Resource acquisition", "Resource acquisition","Resource acquisition","Resource acquisition",
                                                       "Hiding", "Hiding","Hiding","Hiding",
                                                       "Defense","Defense","Defense","Defense",
                                                       "Rescue", "Rescue", "Rescue", "Rescue"
)),
Selection_component= as.character(c("Natural Direct", "Natural Indirect", "Social Direct", "Social Indirect",
                                    "Natural Direct", "Natural Indirect", "Social Direct", "Social Indirect",
                                    "Natural Direct", "Natural Indirect", "Social Direct", "Social Indirect",
                                    "Natural Direct", "Natural Indirect", "Social Direct", "Social Indirect")),
Value = round(as.numeric(c(B_gen_contrib, B_crouch_contrib, B_chase_contrib,B_unhook_contrib)),2),
Si = round(as.numeric(c(B_S[1,1],B_S[1,1],B_S[1,1],B_S[1,1],
                        B_S[2,1],B_S[2,1],B_S[2,1],B_S[2,1],
                        B_S[3,1],B_S[3,1],B_S[3,1],B_S[3,1],
                        B_S[4,1],B_S[4,1],B_S[4,1],B_S[4,1])),2)
)


write.csv(B_contributions, file = "outputs/tables/B_contributions.csv") 



#######################
# 6. Model selection #
######################

#reading in the model outputs from the main analyses on escape if they are not loaded already
mod_esc_1.1.upd <- readRDS("outputs/models_outputs/mod_esc_1.1.upd.rds")
mod_esc_2.1.upd <- readRDS("outputs/models_outputs/mod_esc_2.1.upd.rds")
mod_esc_3.1.upd <- readRDS("outputs/models_outputs/mod_esc_3.1.upd.rds")
mod_esc_4.1.upd <- readRDS("outputs/models_outputs/mod_esc_4.1.upd.rds")

#model performance checks
perf_mod_esc_1.1.upd<-performance::model_performance(mod_esc_1.1.upd) #gives R2 and AIC 
#check_mod_esc_1.1.upd<-performance::check_model(mod_esc_1.1.upd) 
perf_mod_esc_2.1.upd<-performance::model_performance(mod_esc_2.1.upd) 
#check_mod_esc_2.1.upd<-performance::check_model(mod_esc_2.1.upd) 
perf_mod_esc_3.1.upd<-performance::model_performance(mod_esc_3.1.upd) 
#check_mod_esc_3.1.upd<-performance::check_model(mod_esc_3.1.upd)
perf_mod_esc_4.1.upd<-performance::model_performance(mod_esc_4.1.upd) 
#check_mod_esc_4.1.upd<-performance::check_model(mod_esc_4.1.upd) 

#can do more detailed performance checks for the two models we use for linear gradients and for quadratic and correlational 
#checking model performance
#performance::check_overdispersion(mod_esc_1.1.upd) 
#performance::check_convergence(mod_esc_1.1.upd) 
#performance::r2(mod_esc_1.1.upd)
#performance::check_collinearity(mod_esc_1.1.upd)
#performance::binned_residuals(mod_esc_1.1.upd) 
#performance::check_autocorrelation(mod_esc_1.1.upd)
#performance::check_outliers(mod_esc_1.1.upd)

#checking model performance
#performance::check_overdispersion(mod_esc_4.1.upd)
#performance::check_convergence(mod_esc_4.1.upd)
#performance::r2(mod_esc_4.1.upd)
#performance::check_collinearity(mod_esc_4.1.upd)
#performance::binned_residuals(mod_esc_4.1.upd)
#performance::check_autocorrelation(mod_esc_4.1.upd)
#performance::check_outliers(mod_esc_4.1.upd)

perf_mod_bp_1.1.upd<-performance::model_performance(mod_bp_1.1.upd) 
#check_mod_bp_1.1.upd<-performance::check_model(mod_bp_1.1.upd) 
perf_mod_bp_2.1.upd<-performance::model_performance(mod_bp_2.1.upd) 
#check_mod_bp_2.1.upd<-performance::check_model(mod_bp_2.1.upd) 
perf_mod_bp_3.1.upd<-performance::model_performance(mod_bp_3.1.upd)
#check_mod_bp_3.1.upd<-performance::check_model(mod_bp_3.1.upd) 
perf_mod_bp_4.1.upd<-performance::model_performance(mod_bp_4.1.upd) 
#check_mod_bp_4.1.upd<-performance::check_model(mod_bp_4.1.upd) 

# model selection with AICc

#escape
options(na.action = "na.fail")
out.put.aic_esc<-model.sel(mod_esc_1.1.upd, mod_esc_2.1.upd, mod_esc_3.1.upd, mod_esc_4.1.upd)
out.put.aic_esc
out.put.bic_esc<-model.sel(mod_esc_1.1.upd, mod_esc_2.1.upd, mod_esc_3.1.upd, mod_esc_4.1.upd, rank = BIC)
out.put.bic_esc
#points
out.put.aic_bp<-model.sel(mod_bp_1.1.upd, mod_bp_2.1.upd, mod_bp_3.1.upd, mod_bp_4.1.upd)
out.put.aic_bp
out.put.bic_bp<-model.sel(mod_bp_1.1.upd, mod_bp_2.1.upd, mod_bp_3.1.upd, mod_bp_4.1.upd , rank = BIC)
out.put.bic_bp


#reorganize output to make tables

#survival
#bic table
out.put.bic_esc<-as.data.frame(out.put.bic_esc)
out.put.bic_esc2<- dplyr::select(out.put.bic_esc, df, logLik, BIC, delta, weight)
out.put.bic_esc2 <- cbind(rownames(out.put.bic_esc2), data.frame(out.put.bic_esc2, row.names=NULL))
attributes(out.put.bic_esc2) <- NULL
names(out.put.bic_esc2)<-cbind("model","df","logLik", "BIC","deltaBIC", "BICweight")
out.put.bic_esc2<-as.data.frame(out.put.bic_esc2)
#aic table
out.put.aic_esc<-as.data.frame(out.put.aic_esc)
out.put.aic_esc2<- dplyr::select(out.put.aic_esc, df, logLik, AICc, delta, weight)
out.put.aic_esc2 <- cbind(rownames(out.put.aic_esc2), data.frame(out.put.aic_esc2, row.names=NULL))
attributes(out.put.aic_esc2) <- NULL
names(out.put.aic_esc2)<-cbind("model","df","logLik", "AICc","deltaAIC", "AICcweight")
out.put.aic_esc2<-as.data.frame(out.put.aic_esc2)

#in order of best to worse
perf_mod_esc_all<-rbind(perf_mod_esc_4.1.upd, perf_mod_esc_2.1.upd, perf_mod_esc_3.1.upd, perf_mod_esc_1.1.upd)
mod_sel_esc<-cbind(out.put.aic_esc2, out.put.bic_esc2[4:6],perf_mod_esc_all[3:5] )

#points
#bic table
out.put.bic_bp<-as.data.frame(out.put.bic_bp)
out.put.bic_bp2<- dplyr::select(out.put.bic_bp, df, logLik, BIC, delta, weight)
out.put.bic_bp2 <- cbind(rownames(out.put.bic_bp2), data.frame(out.put.bic_bp2, row.names=NULL))
attributes(out.put.bic_bp2) <- NULL
names(out.put.bic_bp2)<-cbind("model","df","logLik", "BIC","deltaBIC", "BICweight")
out.put.bic_bp2<-as.data.frame(out.put.bic_bp2)
#aic table
out.put.aic_bp<-as.data.frame(out.put.aic_bp)
out.put.aic_bp2<- dplyr::select(out.put.aic_bp, df, logLik, AICc, delta, weight)
out.put.aic_bp2 <- cbind(rownames(out.put.aic_bp2), data.frame(out.put.aic_bp2, row.names=NULL))
attributes(out.put.aic_bp2) <- NULL
names(out.put.aic_bp2)<-cbind("model","df","logLik", "AICc","deltaAIC", "AICcweight")
out.put.aic_bp2<-as.data.frame(out.put.aic_bp2)

#in order of best to worse
perf_mod_bp_all<-rbind(perf_mod_bp_4.1.upd, perf_mod_bp_2.1.upd, perf_mod_bp_3.1.upd, perf_mod_bp_1.1.upd)
mod_sel_bp<-cbind(out.put.aic_bp2, out.put.bic_bp2[4:6],perf_mod_bp_all[3:5] )

# TABLE S3

tab_dfs(list(mod_sel_esc,mod_sel_bp),
        titles = c("Model selection survival","Model selection points"),
        col.header = col,
        file = "outputs/tables/table_mod_selection.doc")




###############################################

# R script for statistical analyses
# associated to "Social interactions generate complex selection patterns in virtual worlds"
# Santostefano, Francesca; Fraser Franco, Maxime; Montiglio, Pierre-Olivier
# 01-Nov-2023
# See associated README file for metadata description

#################################################

# This R script is organized as follows:

# 1. Dependencies and setup
# 2. Mixed models survival
##  a. Model 1 (linear gradients)
##  b. Model 2
##  c. Model 3
##  d. Model 4 (quadratic and correlational)
# 3. Selection gradients 
# 4. Selection differentials + Table 1
# 5. Contributions to selection



#################################################
# 1. Dependencies and setup
#################################################

#"R version 4.2.3 (2023-03-15 ucrt)"

# List of packages needed
packages_list <- c(
  "data.table", "lme4", "arm",
  "dplyr", "jtools", "sjPlot"
)

# Load packages - install them if not already present
lapply(
  packages_list, function(package) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package)
      library(package, character.only = TRUE)
    }
  }
)

# Import in data file
df <- fread("data/df_social_selection.csv")



#################################################
# 2. Mixed models survival
#################################################

## a. Model 1 (linear gradients)

f1.1 <- "escaped ~
    Zgen + Zcrouching + Zchase + Zunhook +
    Zsocial_gen+ Zsocial_crouching + Zsocial_chase +
    Zsocial_unhook + Zgame_duration +
    (1|map_code)"

system.time(
  mod_esc_1.1 <- glmer(
    formula = f1.1,
    data = df,
    family = binomial(link = "logit"),
    control = glmerControl(optimizer = "nloptwrap", nAGQ0initStep = TRUE)
  )
)

# Updating model with controls
ss <- getME(mod_esc_1.1, c("theta", "fixef"))
mod_esc_1.1.upd <- update(
  mod_esc_1.1,
  start = ss,
  control = glmerControl(optCtrl = list(maxfun = 2e4))
)

# Saving model object output
saveRDS(mod_esc_1.1.upd, file = "outputs/models_outputs/mod_esc_1.1.upd.rds")

# Model summary
summ(mod_esc_1.1.upd) # log odds
summ(mod_esc_1.1.upd, exp = TRUE, r.squared = FALSE) # odds ratio
coefs_esc_1.1.upd <- summ(mod_esc_1.1.upd)$coeftable

# Function to extract probabilities
logit2prob <- function(logit) {
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

logit2prob(coefs_esc_1.1.upd[, 1]) # probabilities

# Posterior simulations of beta from a glm object
nsim <- 1000
bsim <- arm::sim(mod_esc_1.1.upd, n.sim = nsim)

# Credibility interval for fixef
apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.975))

# Model coefficients
summary(mod_esc_1.1.upd)$coefficients

# Function to copy model output for exporting it later
table <- function(model) {
  a <- apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.975))
  ta <- t(a)
  ta <- as.data.frame(ta)
  b <- summary(model)$coefficients
  b <- as.data.frame(b)
  b <- b[-c(2, 3)]
  tab <- cbind(b, ta)
  tab <- round(tab, 2)
  tab
}

# Table for fixed effects
fixefs_mod_esc_1.1.upd <- (table(mod_esc_1.1.upd)[, -2])
names(fixefs_mod_esc_1.1.upd) <- cbind(
  "Selection gradient",
  "lower 2.5% CI",
  "upper 97.5% CI"
)
fixefs_mod_esc_1.1.upd <- data.frame(
  "Fixed effects" = c(
    "Intercept", "Focal resource acquisition ", "Focal hiding",
    "Focal defense", "Focal rescue", "Social resource acquisition ",
    "Social hiding", "Social defense", "Social rescue", "Game duration"
  ),
  fixefs_mod_esc_1.1.upd
)

# Table for random effects
tab_map_code_mod_esc_1.1.upd <- round(
  t(as.data.frame(quantile(apply(bsim@ranef$map_code[ , , 1], 1, var), prob = c(0.5, 0.025, 0.975)))),
  digits = 2
)
ranefs_mod_esc_1.1.upd <- tab_map_code_mod_esc_1.1.upd
names(ranefs_mod_esc_1.1.upd) <- cbind(
  "Variance", "lower 2.5% CI", "upper 97.5% CI"
)
ranefs_mod_esc_1.1.upd <- data.frame(
  "Random effects" = "Map",
  ranefs_mod_esc_1.1.upd
)

# TABLE S4 (model 1 output)
tab_dfs(
  list(fixefs_mod_esc_1.1.upd, ranefs_mod_esc_1.1.upd),
  titles = c("Model 1", "Model 1"),
  col.header = col,
  file = "outputs/tables/table_mod_esc_1.1.upd.doc"
)


## b. Model 2

f2.1 <- "escaped ~
    Zgen+ Zcrouching + Zchase + Zunhook +
    Zsocial_gen + Zsocial_crouching + Zsocial_chase + Zsocial_unhook +
    Zgen:Zcrouching + Zgen:Zchase + Zgen:Zunhook +
    Zcrouching:Zchase + Zcrouching:Zunhook + Zchase:Zunhook +
    Zsocial_gen:Zsocial_crouching + Zsocial_gen:Zsocial_chase +
    Zsocial_gen:Zsocial_unhook + Zsocial_crouching:Zsocial_chase +
    Zsocial_crouching:Zsocial_unhook + Zsocial_chase:Zsocial_unhook +
    Zgen:Zsocial_gen + Zcrouching:Zsocial_crouching + 
    Zchase:Zsocial_chase + Zunhook:Zsocial_unhook +
    Zgame_duration + (1|map_code)"

system.time(
  mod_esc_2.1 <- glmer(
    formula = f2.1,
    data = df,
    family = binomial(link = "logit"),
    control = glmerControl(optimizer = "nloptwrap", nAGQ0initStep = TRUE)
  )
)

ss <- getME(mod_esc_2.1, c("theta","fixef"))
mod_esc_2.1.upd <- update(
  mod_esc_2.1,
  start = ss,
  control = glmerControl(optCtrl = list(maxfun = 2e4))
)

saveRDS(mod_esc_2.1.upd, file = "outputs/models_outputs/mod_esc_2.1.upd.rds")
#mod_esc_2.1.upd <- readRDS("outputs/models_outputs/mod_esc_2.1.upd.rds")

summary(mod_esc_2.1.upd)


##  c. Model 3

system.time(mod_esc_3.1<- glmer(escaped ~   Zgen + Zcrouching + Zchase+ Zunhook  + 
                                Zsocial_gen+ Zsocial_crouching + Zsocial_chase + Zsocial_unhook+
                                + I(Zgen^2) + I(Zcrouching^2) +  I(Zchase^2) +  I(Zunhook^2) +
                                + I(Zsocial_gen^2) + I(Zsocial_crouching^2) +  I(Zsocial_chase^2) +  I(Zsocial_unhook^2) +
                                Zgame_duration + 
                                (1|map_code), data = df, family = binomial (link = "logit"),
                                control = glmerControl(optimizer = "nloptwrap", nAGQ0initStep = TRUE)) )

ss <- getME(mod_esc_3.1,c("theta","fixef"))
mod_esc_3.1.upd <- update(mod_esc_3.1,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))

saveRDS(mod_esc_3.1.upd, file = "outputs/models_outputs/mod_esc_3.1.upd.rds")
#mod_esc_3.1.upd <- readRDS("outputs/models_outputs/mod_esc_3.1.upd.rds")

summary(mod_esc_3.1.upd)


##  4. Model 4 (quadratic and correlational)

system.time(mod_esc_4.1<- glmer(escaped ~   Zgen + Zcrouching + Zchase+ Zunhook  + 
                                Zsocial_gen+ Zsocial_crouching + Zsocial_chase + Zsocial_unhook+
                                + I(Zgen^2) + I(Zcrouching^2) +  I(Zchase^2) +  I(Zunhook^2) +
                                + I(Zsocial_gen^2) + I(Zsocial_crouching^2) +  I(Zsocial_chase^2) +  I(Zsocial_unhook^2) +
                                Zgen*Zcrouching + Zgen*Zchase + Zgen*Zunhook +
                                Zcrouching*Zchase + Zcrouching*Zunhook + Zchase*Zunhook +
                                Zsocial_gen*Zsocial_crouching + Zsocial_gen*Zsocial_chase + Zsocial_gen*Zsocial_unhook +
                                Zsocial_crouching*Zsocial_chase + Zsocial_crouching*Zsocial_unhook + Zsocial_chase*Zsocial_unhook +
                                Zgen*Zsocial_gen +   Zchase*Zsocial_chase +  Zunhook*Zsocial_unhook +Zcrouching*Zsocial_crouching+
                                Zgame_duration + (1|map_code), data = df, family = binomial (link = "logit"),
                              control = glmerControl(optimizer = "nloptwrap", nAGQ0initStep = TRUE)) )

ss <- getME(mod_esc_4.1,c("theta","fixef"))
mod_esc_4.1.upd <- update(mod_esc_4.1,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))

saveRDS(mod_esc_4.1.upd, file = "outputs/models_outputs/mod_esc_4.1.upd.rds")
#mod_esc_4.1.upd <- readRDS("outputs/models_outputs/mod_esc_4.1.upd.rds")

#summary
summ(mod_esc_4.1.upd) #log odds
summ(mod_esc_4.1.upd, exp = TRUE) #odds ratio
coefs_esc_4.1.upd<-summ(mod_esc_4.1.upd)$coeftable
logit2prob(coefs_esc_4.1.upd[,1]) #probabilities

summary(mod_esc_4.1.upd)


#sim posterior simulations of beta from a glm object
nsim<-1000
bsim<-arm::sim(mod_esc_4.1.upd, n.sim=nsim)

#credibility interval for fixef
apply(bsim@fixef, 2, quantile, prob=c(0.025, 0.975))

#model coefficients
summary(mod_esc_4.1.upd)$coefficients

#table for fixed effects
fixefs_mod_esc_4.1.upd<-(table(mod_esc_4.1.upd)[,-2])

fixefs_mod_esc_4.1.upd <- data.frame("Fixed effects" = c("Intercept", "Focal resource acquisition", "Focal hiding", "Focal defense", "Focal rescue", "Social resource acquisition", "Social hiding", "Social defense", "Social rescue", 
                                                         "Quadratic Focal resource acquisition", "Quadratic Focal hiding", "Quadratic Focal defense", "Quadratic Focal rescue", "Quadratic Social resource acquisition ", "Quadratic Social hiding", "Quadratic Social defense", "Quadratic Social rescue",
                                                         "Game duration", "Focal resource acquisition*hiding", "Focal resource acquisition*defense", "Focal resource acquisition*rescue", "Focal hiding*defense", "Focal hiding*rescue", "Focal defense*rescue",
                                                         "Social resource acquisition*hiding", "Social resource acquisition*defense", "Social resource acquisition*rescue", "Social hiding*defense", "Social hiding*rescue", "Social defense*rescue",
                                                         "Focal*Social resource acquisition", "Focal*Social hiding", "Focal*Social defense", "Focal*Social rescue"
), fixefs_mod_esc_4.1.upd )


#function to get quadratic terms coefficients and CI for those, multiplies x2 the quadratic coefficients

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

quad_mod_esc_4.1.upd<-quad(mod_esc_4.1.upd)
quad_mod_esc_4.1.upd <- data.frame("Fixed effects" = c("Quadratic Focal resource acquisition", "Quadratic Focal hiding", "Quadratic Focal defense", "Quadratic Focal rescue", "Quadratic Social resource acquisition ", "Quadratic Social hiding", "Quadratic Social defense", "Quadratic Social rescue"), quad_mod_esc_4.1.upd )

#replacing with the quadratic coefficients x 2
fixefs_mod_esc_4.1.upd[match(quad_mod_esc_4.1.upd$Fixed.effects, fixefs_mod_esc_4.1.upd$Fixed.effects), ] <- quad_mod_esc_4.1.upd
names(fixefs_mod_esc_4.1.upd)<-cbind("Fixed Effects","Selection gradient", "lower 2.5% CI",  "upper 97.5% CI")

#table for random effects
tab_map_code_mod_esc_4.1.upd<-round(t(as.data.frame(quantile(apply(bsim@ranef$map_code[ , , 1],1,var), prob=c(0.5, 0.025, 0.975)))),2)
ranefs_mod_esc_4.1.upd<-tab_map_code_mod_esc_4.1.upd
names(ranefs_mod_esc_4.1.upd)<-cbind("Variance", "lower 2.5% CI",  "upper 97.5% CI")
ranefs_mod_esc_4.1.upd <- data.frame("Random effects" = "Map", ranefs_mod_esc_4.1.upd)

#TABLE S4 (model 4 output)
tab_dfs(list(fixefs_mod_esc_4.1.upd,ranefs_mod_esc_4.1.upd),
        titles = c("Model 4","Model 4"),
        col.header = col,
        file = "outputs/tables/table_mod_esc_4.1.upd.doc")



#################################################
# 3. Selection gradients
#################################################

#loading models and coefs if they are not loaded already
#mod_esc_1.1.upd <- readRDS("outputs/models_outputs/mod_esc_1.1.upd.rds")
#mod_esc_2.1.upd <- readRDS("outputs/models_outputs/mod_esc_2.1.upd.rds")
#mod_esc_3.1.upd <- readRDS("outputs/models_outputs/mod_esc_3.1.upd.rds")
#mod_esc_4.1.upd <- readRDS("outputs/models_outputs/mod_esc_4.1.upd.rds")
coefs_esc_1.1.upd<-summ(mod_esc_1.1.upd)$coeftable
coefs_esc_2.1.upd<-summ(mod_esc_2.1.upd)$coeftable
coefs_esc_3.1.upd<-summ(mod_esc_3.1.upd)$coeftable
coefs_esc_4.1.upd<-summ(mod_esc_4.1.upd)$coeftable

#extracting log odds linear coefficients (betas - selection gradients)
coefs_esc_1.1.upd[,1]
#need the linear gradients for focal and social
E_bn_gen<-coefs_esc_1.1.upd[2,1]
E_bn_crouching<-coefs_esc_1.1.upd[3,1]
E_bn_chase<-coefs_esc_1.1.upd[4,1]
E_bn_unhook<-coefs_esc_1.1.upd[5,1]
E_bs_gen<-coefs_esc_1.1.upd[6,1]
E_bs_crouching<-coefs_esc_1.1.upd[7,1]
E_bs_chase<-coefs_esc_1.1.upd[8,1]
E_bs_unhook<-coefs_esc_1.1.upd[9,1]



#################################################
# 4. Selection differentials
#################################################

#EQ 11 from Wolf et al. 1999 paper: S=P*bN+CI*bS

#organizing vectors of gradients
E_bn = matrix(
  c( E_bn_gen,
     E_bn_crouching,
     E_bn_chase, 
     E_bn_unhook
  ),
  nrow=4,
  ncol=1,
  byrow = TRUE) 

E_bn

E_bs = matrix(
  c( E_bs_gen, 
     E_bs_crouching,
     E_bs_chase, 
     E_bs_unhook
  ),
  nrow=4,
  ncol=1,
  byrow = TRUE) 

E_bs


# (variance is 1 for all traits in the full dataset)
V=1  

#matrix of VCV for focal traits, P  
focaltraits <- select(df, Zgen, Zcrouching, Zchase, Zunhook)
P<-cov(focaltraits)
P

write.csv(P, file = "outputs/tables/P.csv", col.names=FALSE)

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

write.csv(CI, file = "outputs/tables/CI.csv", col.names=FALSE) 

#equation for the selection differentials S=P*bN+CI*bS

E_S=P%*%E_bn+CI%*%E_bs
E_S

# TABLE 1 with gradients and selection differentials

Esc_rate_gradients <- data.frame(Behaviour=as.character(c("Resource acquisition", "Hiding", "Defense", "Rescue")),
                                Interaction_type_predicted= as.character(c("Cooperation +/+", "Selfish +/-", "Altruistic -/+", "Altruistic -/+")),
                                ßN = E_bn,
                                ßS = E_bs,
                                Si = E_S)
                                
Esc_rate_gradients 

tab_df(Esc_rate_gradients,
       titles =   "B) Selection gradients Survival",
       col.header = col,
       file = "outputs/tables/table_gradients_survival.doc")


#################################################
# 5. Contributions to selection
#################################################

E_S=P%*%E_bn+CI%*%E_bs

P
E_bn
CI
E_bs

#natural selection contrib
P_E_bn<-P%*%E_bn
P_E_bn

#social sel contrib
CI_E_bs<-CI%*%E_bs
CI_E_bs

#Si
E_S

#direct and indirect
#nat sel direct
D_E_bn_gen<-P[1,1]%*%E_bn[1,1]
D_E_bn_crouch<-P[2,2]%*%E_bn[2,1]
D_E_bn_chase<-P[3,3]%*%E_bn[3,1]
D_E_bn_unhook<-P[4,4]%*%E_bn[4,1]

#nat sel indirect
I_E_bn_gen<-P[1,2]%*%E_bn[2,1]+P[1,3]%*%E_bn[3,1]+P[1,4]%*%E_bn[4,1]
I_E_bn_crouch<-P[2,1]%*%E_bn[1,1]+P[2,3]%*%E_bn[3,1]+P[2,4]%*%E_bn[4,1]
I_E_bn_chase<-P[3,1]%*%E_bn[1,1]+P[3,2]%*%E_bn[2,1]+P[3,4]%*%E_bn[4,1]
I_E_bn_unhook<-P[4,1]%*%E_bn[1,1]+P[4,2]%*%E_bn[2,1]+P[4,3]%*%E_bn[3,1]

#soc sel direct
D_E_bs_gen<-CI[1,1]%*%E_bs[1,1]
D_E_bs_crouch<-CI[2,2]%*%E_bs[2,1]
D_E_bs_chase<-CI[3,3]%*%E_bs[3,1]
D_E_bs_unhook<-CI[4,4]%*%E_bs[4,1]

#soc sel indirect
I_E_bs_gen<-CI[1,2]%*%E_bs[2,1]+CI[1,3]%*%E_bs[3,1]+CI[1,4]%*%E_bs[4,1]
I_E_bs_crouch<-CI[2,1]%*%E_bs[1,1]+CI[2,3]%*%E_bs[3,1]+CI[2,4]%*%E_bs[4,1]
I_E_bs_chase<-CI[3,1]%*%E_bs[1,1]+CI[3,2]%*%E_bs[2,1]+CI[3,4]%*%E_bs[4,1]
I_E_bs_unhook<-CI[4,1]%*%E_bs[1,1]+CI[4,2]%*%E_bs[2,1]+CI[4,3]%*%E_bs[3,1]


#TABLE S5 - contributions to selection (also used for fig. 1 )

E_gen_contrib<-c(D_E_bn_gen,I_E_bn_gen,D_E_bs_gen,I_E_bs_gen )
E_crouch_contrib<-c(D_E_bn_crouch,I_E_bn_crouch,D_E_bs_crouch,I_E_bs_crouch )
E_chase_contrib<-c(D_E_bn_chase,I_E_bn_chase,D_E_bs_chase,I_E_bs_chase )
E_unhook_contrib<-c(D_E_bn_unhook,I_E_bn_unhook,D_E_bs_unhook,I_E_bs_unhook )

E_contributions <- data.frame(Behaviour=as.character(c("Resource acquisition", "Resource acquisition","Resource acquisition","Resource acquisition",
                                                       "Hiding", "Hiding","Hiding","Hiding",
                                                       "Defense","Defense","Defense","Defense",
                                                       "Rescue", "Rescue", "Rescue", "Rescue"
)),
Selection_component= as.character(c("Natural Direct", "Natural Indirect", "Social Direct", "Social Indirect",
                                    "Natural Direct", "Natural Indirect", "Social Direct", "Social Indirect",
                                    "Natural Direct", "Natural Indirect", "Social Direct", "Social Indirect",
                                    "Natural Direct", "Natural Indirect", "Social Direct", "Social Indirect")),
Value = round(as.numeric(c(E_gen_contrib, E_crouch_contrib, E_chase_contrib,E_unhook_contrib)),2),
Si = round(as.numeric(c(E_S[1,1],E_S[1,1],E_S[1,1],E_S[1,1],
                        E_S[2,1],E_S[2,1],E_S[2,1],E_S[2,1],
                        E_S[3,1],E_S[3,1],E_S[3,1],E_S[3,1],
                        E_S[4,1],E_S[4,1],E_S[4,1],E_S[4,1])),2)
)

E_contributions

write.csv(E_contributions, file = "outputs/tables/E_contributions.csv") 


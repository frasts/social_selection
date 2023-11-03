###############################################

# R script for plots
# Associated to "Social interactions generate complex selection patterns in virtual worlds"
# Santostefano, Francesca; Fraser Franco, Maxime; Montiglio, Pierre-Olivier
# 01-Nov-2023
# See associated Readme file for metadata description

###############################################

# This R script is organized as follows:

# 1. Dependencies and setup
# 2. Figure 1 - selection gradients survival 
##  a. Individual traits - natural selection
##  b. Individual traits - social selection
##  c. Paneling all traits 
# 3. Figure 2 - contributions to selection differential

#######################
# 1. Dependencies and setup #
######################

#"R version 4.2.3 (2023-03-15 ucrt)"

#Set working directory to you own path where the folder "social_selection" is
setwd("/path/to/folder/social_selection")

# Load packages - install them if not already present
packages_list <- c("sjPlot","data.table","cowplot", "ggplot2")

# Load packages
lapply(packages_list, function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
})


#read in data file
df <- fread("data/df_social_selection.csv")

#######################
# 2. Figure 1 - selection gradients survival 
######################

#plots of predicted values (marginal effects) for specific model terms
#Marginal effects measure the impact that an instantaneous unit change in one variable has on the outcome variable while all other variables are held constant.
#given the plots with very few extreme high Z values, cut the margins for better visualization and help interpretation

##  a. Individual traits - natural selection


#resource acquisition
mod_esc_4.1.upd <- readRDS("outputs/models_outputs/mod_esc_4.1.upd.rds")

mod_esc_4.1.upd_gen<-plot_model(mod_esc_4.1.upd, type = "pred", terms = "Zgen[all]") +
  theme_classic() + theme(plot.title = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y=element_text(size=10),
                          axis.text.x=element_text(size=10),
                          axis.title.x=element_text(size=10))+
  scale_x_continuous(breaks = (seq(0, 600, 100) - mean(df$gen_duration)) / sd(df$gen_duration), 
                     labels = seq(0, 600, 100))+
  scale_y_continuous(breaks = seq(0, 1, .20)) +
  xlab("\nResource acquisition (%)")

#saving full plot before cutting
#saveRDS(mod_esc_4.1.upd_gen, file = "outputs/models_plots/mod_esc_4.1.upd_gen.rds")
#mod_esc_4.1.upd_gen <- readRDS("outputs/models_plots/mod_esc_4.1.upd_gen.rds")

#cutting the plot, adding quantile
qZgen<-as.vector(quantile(df$Zgen, prob = 0.95))
plot_mod_esc_4.1.upd_gen<-mod_esc_4.1.upd_gen + 
  coord_cartesian(xlim =c(-1.5,3.5), ylim = c(0, 1))+
  geom_point(data=df, show.legend = F, alpha = 1/50, aes(x = Zgen, y = escaped))+
  geom_vline(xintercept=qZgen, colour="grey30", size= 0.5, linetype=2)

ggsave("plot_mod_esc_4.1.upd_gen.png", path = "outputs/models_plots")
#plot_mod_esc_4.1.upd_gen <- readPNG("outputs/models_plots/plot_mod_esc_4.1.upd_gen.png")


#active defense
mod_esc_4.1.upd_chase<-plot_model(mod_esc_4.1.upd, type = "pred", terms = "Zchase[all]")+ 
  theme_classic() + theme(plot.title = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y=element_text(size=10),
                          axis.text.x=element_text(size=10),
                          axis.title.x=element_text(size=10))+
  scale_x_continuous(breaks = (seq(0, 300, 50) - mean(df$chase_duration)) / sd(df$chase_duration), 
                     labels = seq(0, 300, 50))+
  scale_y_continuous(breaks = seq(0, 1, .20)) +
  xlab("\nChase duration (sec)")

#saving full plot before cutting
#saveRDS(mod_esc_4.1.upd_chase, file = "outputs/models_plots/mod_esc_4.1.upd_chase.rds")
#mod_esc_4.1.upd_chase <- readRDS("outputs/models_plots/mod_esc_4.1.upd_chase.rds")

#cutting the plot, adding quantile
qZchase<-as.vector(quantile(df$Zchase, prob = 0.95))
plot_mod_esc_4.1.upd_chase<-mod_esc_4.1.upd_chase + coord_cartesian(xlim =c(-1.5,5.5), ylim = c(0, 1))+
  geom_point(data=df, show.legend = F,  alpha = 1/50, aes(x = Zchase, y = escaped))+
  geom_vline(xintercept=qZchase, colour="grey30", size= 0.5, linetype=2)

ggsave("plot_mod_esc_4.1.upd_chase.png", path = "outputs/models_plots")
#plot_mod_esc_4.1.upd_chase <- readPNG("outputs/models_plots/plot_mod_esc_4.1.upd_chase.png")


#rescue events
mod_esc_4.1.upd_unhook<-plot_model(mod_esc_4.1.upd, type = "pred", terms = "Zunhook[all]")+ 
  theme_classic() + theme(plot.title = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y=element_text(size=10),
                          axis.text.x=element_text(size=10),
                          axis.title.x=element_text(size=10))+
  scale_x_continuous(breaks = (seq(0, 6, 1) - mean(df$unhook_count)) / sd(df$unhook_count), 
                     labels = seq(0, 6, 1))+
  scale_y_continuous(breaks = seq(0, 1, .20)) +
  xlab("\n# rescues") 

#saving full plot before cutting
#saveRDS(mod_esc_4.1.upd_unhook, file = "outputs/models_plots/mod_esc_4.1.upd_unhook.rds")
#mod_esc_4.1.upd_unhook <- readRDS("outputs/models_plots/mod_esc_4.1.upd_unhook.rds")

#cutting the plot, adding quantile
qZunhook<-as.vector(quantile(df$Zunhook, prob = 0.95))
plot_mod_esc_4.1.upd_unhook<-mod_esc_4.1.upd_unhook + coord_cartesian(xlim =c(-1,4.5), ylim = c(0, 1))+
  geom_point(data=df, show.legend = F, alpha = 1/50, aes(x = Zunhook, y = escaped))+
  geom_vline(xintercept=qZunhook, colour="grey30", size= 0.5, linetype=2)

ggsave("plot_mod_esc_4.1.upd_unhook.png", path = "outputs/models_plots")
#plot_mod_esc_4.1.upd_unhook <- readPNG("outputs/models_plots/plot_mod_esc_4.1.upd_unhook.png")


#hiding
mod_esc_4.1.upd_crouching<-plot_model(mod_esc_4.1.upd, type = "pred", terms = "Zcrouching[all]")+ 
  theme_classic() + theme(plot.title = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y=element_text(size=10),
                          axis.text.x=element_text(size=10),
                          axis.title.x=element_text(size=10))+
  scale_x_continuous(breaks = (seq(0, 1400, 200) - mean(df$crouching_duration)) / sd(df$crouching_duration), #change seq breaks
                     labels = seq(0, 1400, 200))+
  scale_y_continuous(breaks = seq(0, 1, .20)) +
  xlab("\n Hiding duration") 

#saving full plot before cutting
#saveRDS(mod_esc_4.1.upd_crouching, file = "outputs/models_plots/mod_esc_4.1.upd_crouching.rds")
#mod_esc_4.1.upd_crouching <- readRDS("outputs/models_plots/mod_esc_4.1.upd_crouching.rds")

#cutting the plot, adding quantile
qZcrouching<-as.vector(quantile(df$Zcrouching, prob = 0.95))
plot_mod_esc_4.1.upd_crouching<-mod_esc_4.1.upd_crouching + coord_cartesian(xlim =c(-1,6), ylim = c(0, 1))+
  geom_point(data=df, show.legend = F, alpha = 1/50, aes(x = Zcrouching, y = escaped))+
  geom_vline(xintercept=qZcrouching, colour="grey30", size= 0.5, linetype=2)

ggsave("plot_mod_esc_4.1.upd_crouching.png", path = "outputs/models_plots")
#plot_mod_esc_4.1.upd_crouching <- readPNG("outputs/models_plots/mod_esc_4.1.upd_crouching.png")


##  b. Individual traits - social selection

#resource acquisition
mod_esc_4.1.upd_gen_social<-plot_model(mod_esc_4.1.upd, type = "pred", terms = "Zsocial_gen[all]")+ 
   theme_classic() + theme(plot.title = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y=element_text(size=10),
                          axis.text.x=element_text(size=10),
                          axis.title.x=element_text(size=10))+
  scale_x_continuous(breaks = (seq(0, 400, 100) - mean(df$avg_social_gen)) / sd(df$avg_social_gen), #change seq breaks
                     labels = seq(0, 400, 100))+
  scale_y_continuous(breaks = seq(0, 1, .20)) +
   xlab("\nAvg resource acquisition by conspecifics (%)")

#saving full plot before cutting
#saveRDS(mod_esc_4.1.upd_gen_social, file = "outputs/models_plots/mod_esc_4.1.upd_gen_social.rds")
#mod_esc_4.1.upd_gen_social <- readRDS("outputs/models_plots/mod_esc_4.1.upd_gen_social.rds")

#cutting the plot, adding quantile
qZsocial_gen<-as.vector(quantile(df$Zsocial_gen, prob = 0.95))

plot_mod_esc_4.1.upd_gen_social<-mod_esc_4.1.upd_gen_social + coord_cartesian(xlim =c(-2.5,3), ylim = c(0, 1))+
  geom_point(data=df, show.legend = F,  alpha = 1/50, aes(x = Zsocial_gen, y = escaped))+
  geom_vline(xintercept=qZsocial_gen, colour="grey30", size= 0.5, linetype=2)

ggsave("plot_mod_esc_4.1.upd_gen_social.png", path = "outputs/models_plots")
#plot_mod_esc_4.1.upd_gen_social <- readPNG("outputs/models_plots/mod_esc_4.1.upd_gen_social.png")


#active defense
mod_esc_4.1.upd_chase_social<-plot_model(mod_esc_4.1.upd, type = "pred", terms = "Zsocial_chase[all]")+ 
  theme_classic() + theme(plot.title = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y=element_text(size=10),
                          axis.text.x=element_text(size=10),
                          axis.title.x=element_text(size=10))+
  scale_x_continuous(breaks = (seq(0, 200, 50) - mean(df$avg_social_chase)) / sd(df$avg_social_chase), #change seq breaks
                     labels = seq(0, 200, 50))+
  scale_y_continuous(breaks = seq(0, 1, .20)) +
  xlab("\nAvg chase duration by conspecifics (sec)")

#saving full plot before cutting
#saveRDS(mod_esc_4.1.upd_chase_social, file = "outputs/models_plots/mod_esc_4.1.upd_chase_social.rds")
#mod_esc_4.1.upd_chase_social <- readRDS("outputs/models_plots/mod_esc_4.1.upd_chase_social.rds")

#cutting the plot, adding quantile
qZsocial_chase<-as.vector(quantile(df$Zsocial_chase, prob = 0.95))
plot_mod_esc_4.1.upd_chase_social<-mod_esc_4.1.upd_chase_social + coord_cartesian(xlim =c(-2,4.5), ylim = c(0, 1))+
  geom_point(data=df, show.legend = F, alpha = 1/50, aes(x = Zsocial_chase, y = escaped))+
  geom_vline(xintercept=qZsocial_chase, colour="grey30", size= 0.5, linetype=2)

ggsave("plot_mod_esc_4.1.upd_chase_social.png", path = "outputs/models_plots")
#plot_mod_esc_4.1.upd_chase_social <- readPNG("outputs/models_plots/mod_esc_4.1.upd_chase_social.png")


#rescuing 
mod_esc_4.1.upd_unhook_social<-plot_model(mod_esc_4.1.upd, type = "pred", terms = "Zsocial_unhook[all]")+
  theme_classic() + theme(plot.title = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y=element_text(size=10),
                          axis.text.x=element_text(size=10),
                          axis.title.x=element_text(size=10))+
  scale_x_continuous(breaks = (seq(0, 2.7, 0.5) - mean(df$avg_social_unhook)) / sd(df$avg_social_unhook), 
                     labels = seq(0, 2.7, 0.5))+
  scale_y_continuous(breaks = seq(0, 1, .20)) +
  xlab("\nAvg # rescues by conspecifics")

#saving full plot before cutting
#saveRDS(mod_esc_4.1.upd_unhook_social, file = "outputs/models_plots/mod_esc_4.1.upd_unhook_social.rds")
#mod_esc_4.1.upd_unhook_social <- readRDS("outputs/models_plots/mod_esc_4.1.upd_unhook_social.rds")

#cutting the plot, adding quantile
qZsocial_unhook<-as.vector(quantile(df$Zsocial_unhook, prob = 0.95))
plot_mod_esc_4.1.upd_unhook_social<-mod_esc_4.1.upd_unhook_social +
  geom_point(data=df, show.legend = F,  alpha = 1/50, aes(x = Zsocial_unhook, y = escaped))+
  geom_vline(xintercept=qZsocial_unhook, colour="grey30", size= 0.5, linetype=2)

ggsave("plot_mod_esc_4.1.upd_unhook_social.png", path = "outputs/models_plots")
#plot_mod_esc_4.1.upd_unhook_social.png <- readPNG("C:/Users/ext_fsantostefano/Projects/Social selection/Social-selection/outputs_2/models_plots/mod_esc_4.1.upd_unhook_social.png")


#hiding
mod_esc_4.1.upd_crouching_social<-plot_model(mod_esc_4.1.upd, type = "pred", terms = "Zsocial_crouching[all]")+ 
  theme_classic() + theme(plot.title = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y=element_text(size=10),
                          axis.text.x=element_text(size=10),
                          axis.title.x=element_text(size=10))+
  scale_x_continuous(breaks = (seq(0, 1100, 200) - mean(df$avg_social_crouching)) / sd(df$avg_social_crouching),
                     labels = seq(0, 1100, 200))+
  scale_y_continuous(breaks = seq(0, 1, .20)) +
  xlab("\nAvg time spent hiding by conspecifics")

#saving full plot before cutting
#saveRDS(mod_esc_4.1.upd_crouching_social, file = "outputs/models_plots/mod_esc_4.1.upd_crouching_social.rds")
#mod_esc_4.1.upd_crouching_social <- readRDS("outputs/models_plots/mod_esc_4.1.upd_crouching_social.rds")

#cutting the plot, adding quantile
qZsocial_crouching<-as.vector(quantile(df$Zsocial_crouching, prob = 0.95))
plot_mod_esc_4.1.upd_crouching_social<-mod_esc_4.1.upd_crouching_social +
  geom_point(data=df, show.legend = F,  alpha = 1/50, aes(x = Zsocial_crouching, y = escaped))+
  geom_vline(xintercept=qZsocial_crouching, colour="grey30", size= 0.5, linetype=2)

ggsave("plot_mod_esc_4.1.upd_crouching_social.png", path = "outputs/models_plots")
#plot_mod_esc_4.1.upd_crouching_social.png <- readPNG("C:/Users/ext_fsantostefano/Projects/Social selection/Social-selection/outputs_2/models_plots/mod_esc_4.1.upd_crouching_social.png")


##  c. Paneling all traits 

panel_plot <- cowplot::plot_grid(plot_mod_esc_4.1.upd_gen,
                                 plot_mod_esc_4.1.upd_crouching,
                                 plot_mod_esc_4.1.upd_chase,
                                 plot_mod_esc_4.1.upd_unhook,
                                 plot_mod_esc_4.1.upd_gen_social,
                                 plot_mod_esc_4.1.upd_crouching_social,
                                 plot_mod_esc_4.1.upd_chase_social,
                                 plot_mod_esc_4.1.upd_unhook_social,
                                 ncol = 4, nrow = 2,
                                 labels = c("(A) Resource acquisition", "(B) Hiding", "(C) Defense", "(D) Rescue","","","",""),
                                 label_size = 10,
                                 label_x = 0,
                                 label_y = 1,
                                 hjust = -0.5,
                                 vjust = 0.5)

ggsave("panel_plot_mod_esc_4.1.upd.png", width = 10, height = 5, path = "outputs/models_plots")


#######################
# 3. Figure 2 - contributions to selection differentials
######################

E_contributions <- fread(file = "outputs/tables/E_contributions.csv")

plot_E_contributions<- ggplot(E_contributions, aes(x=Behaviour, y=Value,fill=Selection_component)) + 
  geom_bar (stat="identity",colour="gray45",width = 0.6) +
  scale_y_continuous(breaks = seq(-0.5,1,0.25),limits = c(-0.5,1))+
  geom_hline(yintercept=0, colour="black", linewidth= 1, linetype=2)+
  theme_bw()+
  scale_fill_manual(values=c("firebrick", "coral", "dodgerblue", "lightskyblue"), name = "Selection component")+
  geom_point(data = E_contributions, aes(x = Behaviour, y= Si), size=3)+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=14,face="bold"),
        legend.position="none")+#,
  labs(y = "Selection differential contribution")+
  coord_flip()

ggsave("plot_E_contributions.png", width = 10, height = 4, path = "outputs/models_plots")



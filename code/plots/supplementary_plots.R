###############################################

# R script for plots
# Associated to "Social interactions generate complex selection patterns in virtual worlds"
# Santostefano, Francesca; Fraser Franco, Maxime; Montiglio, Pierre-Olivier
# 01-Nov-2023
# See associated Readme file for metadata description

###############################################
###############################################

# This R script is organized as follows:

# 1. Dependencies and setup
# 2. Figure S4 - selection gradients points 
##  a. Individual traits - natural selection
##  b. Individual traits - social selection
##  c. Paneling all traits 
# 3. Figure S5 - contributions to selection differential points
# 4. Figure S3 - interactions survival
# 5. Figure S2 - P and CI matrix

#######################
# 1. Dependencies and setup #
######################

#"R version 4.2.3 (2023-03-15 ucrt)"

#Set working directory to you own path where the folder "social_selection" is
setwd("/path/to/folder/social_selection")

# Load packages - install them if not already present
packages_list <- c("sjPlot","data.table","cowplot","ggplot2","plotly", "htmlwidgets","lme4" ,"corrplot", "RColorBrewer")

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
# 2. Figure S4 - selection gradients points 
######################

##  a. Individual traits - natural selection

#resource acquisition
mod_bp_4.1.upd <- readRDS("outputs/models_outputs/mod_bp_4.1.upd.rds")

mod_bp_4.1.upd_gen<-plot_model(mod_bp_4.1.upd, type = "pred", terms = "Zgen[all]")+
  theme_classic() + theme(plot.title = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y=element_text(size=10),
                          axis.text.x=element_text(size=10),
                          axis.title.x=element_text(size=10))+
  scale_x_continuous(breaks = (seq(0, 615, 100) - mean(df$gen_duration)) / sd(df$gen_duration), 
                     labels = seq(0, 615, 100))+
  scale_y_continuous(breaks = seq(0, 1, .20),labels = seq(0, 32000, 6400), expand = c(0, 0),limits = seq(0, 32000))  +
  xlab("\nResource acquisition (%)")

#cutting the plot and adding raw data
plot_mod_bp_4.1.upd_gen<-mod_bp_4.1.upd_gen + coord_cartesian(xlim =c(-1.5,3.5), ylim = c(0, 1))+
  geom_point(data=df, show.legend = F, size = 1, alpha = 1/25, aes(x = Zgen, y = bloodpoints/32000))

ggsave("plot_mod_bp_4.1.upd_gen.png", path = "outputs/models_plots")


#active defense
mod_bp_4.1.upd_chase<-plot_model(mod_bp_4.1.upd, type = "pred", terms = "Zchase[all]", show.data = T) +
  theme_classic() + theme(plot.title = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y=element_text(size=10),
                          axis.text.x=element_text(size=10),
                          axis.title.x=element_text(size=10))+
  scale_x_continuous(breaks = (seq(0, 300, 50) - mean(df$chase_duration)) / sd(df$chase_duration), 
                     labels = seq(0, 300, 50))+
  scale_y_continuous(breaks = seq(0, 1, .20),labels = seq(0, 32000, 6400), expand = c(0, 0),limits = seq(0, 32000))  +
  xlab("\nChase duration (sec)")


#cutting the plot if needed
plot_mod_bp_4.1.upd_chase<-mod_bp_4.1.upd_chase+
  geom_point(data=df, show.legend = F, size = 1, alpha = 1/25, aes(x = Zchase, y = bloodpoints/32000))

ggsave("plot_mod_bp_4.1.upd_chase.png", path = "outputs/models_plots")


#rescue
mod_bp_4.1.upd_unhook<-plot_model(mod_bp_4.1.upd, type = "pred", terms = "Zunhook[all]", show.data = T )+ 
  theme_classic() + theme(plot.title = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y=element_text(size=10),
                          axis.text.x=element_text(size=10),
                          axis.title.x=element_text(size=10))+
  scale_x_continuous(breaks = (seq(0, 6, 1) - mean(df$unhook_count)) / sd(df$unhook_count), 
                     labels = seq(0, 6, 1))+
  scale_y_continuous(breaks = seq(0, 1, .20),labels = seq(0, 32000, 6400),limits = seq(0, 32000))  +
  xlab("\n# rescues") 

#cutting the plot if needed
plot_mod_bp_4.1.upd_unhook<-mod_bp_4.1.upd_unhook +
  geom_point(data=df, show.legend = F, size = 1, alpha = 1/25, aes(x = Zunhook, y = bloodpoints/32000))

ggsave("plot_mod_bp_4.1.upd_unhook.png", path = "outputs/models_plots")


#hiding
mod_bp_4.1.upd_crouching<-plot_model(mod_bp_4.1.upd, type = "pred", terms = "Zcrouching[all]", show.data = T )+
  theme_classic() + theme(plot.title = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y=element_text(size=10),
                          axis.text.x=element_text(size=10),
                          axis.title.x=element_text(size=10))+
  scale_x_continuous(breaks = (seq(0, 1500, 400) - mean(df$crouching_duration)) / sd(df$crouching_duration), #change seq breaks
                     labels = seq(0, 1500, 400))+
  scale_y_continuous(breaks = seq(0, 1, .20),labels = seq(0, 32000, 6400), expand = c(0, 0),limits = seq(0, 32000))  +
  xlab("\n Hiding duration (sec)") 

#cutting the plot if needed and adding raw data
plot_mod_bp_4.1.upd_crouching<-mod_bp_4.1.upd_crouching +
  geom_point(data=df, show.legend = F, size = 1, alpha = 1/25, aes(x = Zcrouching, y = bloodpoints/32000))

ggsave("plot_mod_bp_4.1.upd_crouching.png", path = "outputs/models_plots")


##  b. Individual traits - social selection

#resource acquisition
mod_bp_4.1.upd_gen_social<-plot_model(mod_bp_4.1.upd, type = "pred", terms = "Zsocial_gen[all]", show.data = T )+ 
  theme_classic() + theme(plot.title = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y=element_text(size=10),
                          axis.text.x=element_text(size=10),
                          axis.title.x=element_text(size=10))+
  scale_x_continuous(breaks = (seq(0, 450, 100) - mean(df$avg_social_gen)) / sd(df$avg_social_gen), #change seq breaks
                     labels = seq(0, 450, 100))+
  scale_y_continuous(breaks = seq(0, 1, .20),labels = seq(0, 32000, 6400), expand = c(0, 0),limits = seq(0, 32000))  +
  xlab("\nAvg resource acquisition by conspecifics (%)")

#cutting the plot if needed and adding raw data
plot_mod_bp_4.1.upd_gen_social<-mod_bp_4.1.upd_gen_social + coord_cartesian(xlim =c(-2.5,3), ylim = c(0, 1))+
  geom_point(data=df, show.legend = F, size = 1, alpha = 1/25, aes(x = Zsocial_gen, y = bloodpoints/32000))

ggsave("plot_mod_bp_4.1.upd_gen_social.png", path = "outputs/models_plots")


#active defense
mod_bp_4.1.upd_chase_social<-plot_model(mod_bp_4.1.upd, type = "pred", terms = "Zsocial_chase[all]", show.data = T )+
  theme_classic() + theme(plot.title = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y=element_text(size=10),
                          axis.text.x=element_text(size=10),
                          axis.title.x=element_text(size=10))+
  scale_x_continuous(breaks = (seq(0, 160, 50) - mean(df$avg_social_chase)) / sd(df$avg_social_chase), #change seq breaks
                     labels = seq(0, 160, 50))+
  scale_y_continuous(breaks = seq(0, 1, .20),labels = seq(0, 32000, 6400), expand = c(0, 0),limits = seq(0, 32000))  +
  xlab("\nAvg chase duration by conspecifics (sec)")

#cutting the plot if needed and adding raw data
plot_mod_bp_4.1.upd_chase_social<-mod_bp_4.1.upd_chase_social + coord_cartesian(xlim =c(-2,4.5), ylim = c(0, 1))+
  geom_point(data=df, show.legend = F, size = 1, alpha = 1/25, aes(x = Zsocial_chase, y = bloodpoints/32000))

ggsave("plot_mod_bp_4.1.upd_chase_social.png", path = "outputs/models_plots")


#rescue
mod_bp_4.1.upd_unhook_social<-plot_model(mod_bp_4.1.upd, type = "pred", terms = "Zsocial_unhook[all]", show.data = T )+
  theme_classic() + theme(plot.title = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y=element_text(size=10),
                          axis.text.x=element_text(size=10),
                          axis.title.x=element_text(size=10))+
  scale_x_continuous(breaks = (seq(0, 3, 0.5) - mean(df$avg_social_unhook)) / sd(df$avg_social_unhook), #change seq breaks
                     labels = seq(0, 3, 0.5))+
  scale_y_continuous(breaks = seq(0, 1, .20),labels = seq(0, 32000, 6400), expand = c(0, 0),limits = seq(0, 32000))  +
  xlab("\n Avg # rescues by conspecifics")

#cutting the plot if needed and adding raw data
plot_mod_bp_4.1.upd_unhook_social<-mod_bp_4.1.upd_unhook_social +
  geom_point(data=df, show.legend = F, size = 1, alpha = 1/25, aes(x = Zsocial_unhook, y = bloodpoints/32000))

ggsave("plot_mod_bp_4.1.upd_unhook_social.png", path = "outputs/models_plots")


#hiding
mod_bp_4.1.upd_crouching_social<-plot_model(mod_bp_4.1.upd, type = "pred", terms = "Zsocial_crouching[all]", show.data = T )+
  theme_classic() + theme(plot.title = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y=element_text(size=10),
                          axis.text.x=element_text(size=10),
                          axis.title.x=element_text(size=10))+
  scale_x_continuous(breaks = (seq(0, 1100, 400) - mean(df$avg_social_crouching)) / sd(df$avg_social_crouching), 
                     labels = seq(0, 1100, 400))+
  scale_y_continuous(breaks = seq(0, 1, .20),labels = seq(0, 32000, 6400), expand = c(0, 0),limits = seq(0, 32000))  +
  xlab("\n Avg time spent hiding by conspecifics")

plot_mod_bp_4.1.upd_crouching_social<-mod_bp_4.1.upd_crouching_social +
  geom_point(data=df, show.legend = F, size = 1, alpha = 1/25, aes(x = Zsocial_crouching, y = bloodpoints/32000))

ggsave("plot_mod_bp_4.1.upd_crouching_social.png", path = "outputs/models_plots")

##  c. Paneling all traits 

panel_plot <- plot_grid(plot_mod_bp_4.1.upd_gen,
                        plot_mod_bp_4.1.upd_crouching,
                        plot_mod_bp_4.1.upd_chase,
                        plot_mod_bp_4.1.upd_unhook,
                        plot_mod_bp_4.1.upd_gen_social,
                        plot_mod_bp_4.1.upd_crouching_social,
                        plot_mod_bp_4.1.upd_chase_social,
                        plot_mod_bp_4.1.upd_unhook_social,
                        ncol = 4, nrow = 2,
                        labels = c("(A) Resource acquisition", "(B) Hiding", "(C) Active defense", "(D) Rescue events","","","",""),
                        label_size = 6,
                        label_x = 0,
                        label_y = 1,
                        hjust = -0.5,
                        vjust = 0.5)

ggsave("panel_plot_mod_bp_4.1.upd.png",  width = 10, height = 5, path = "outputs/models_plots")


#######################
# 3. Figure S5 - contributions to selection differentials
######################

B_contributions<- fread(file = "outputs/tables/B_contributions.csv") 

plot_B_contributions<- ggplot(B_contributions, aes(x=Behaviour, y=Value,fill=Selection_component)) + 
  geom_bar (stat="identity",colour="gray45",width = 0.6) +
  scale_y_continuous(breaks = seq(-0.5,1,0.25),limits = c(-0.5,1))+
  geom_hline(yintercept=0, colour="black", linewidth= 1, linetype=2)+
  theme_bw()+
  scale_fill_manual(values=c("firebrick", "coral", "dodgerblue", "lightskyblue"), name = "Selection component")+
  geom_point(data = B_contributions, aes(x = Behaviour, y= Si), size=3)+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=14,face="bold"),
        legend.position="none")+#,
  labs(y = "Selection differential contribution")+
  coord_flip()

ggsave("plot_B_contributions.png", width = 10, height = 4, path = "outputs/models_plots")



#######################
# 4. Figure S3 - interactions survival
######################

#reading in the model outputs from the main analyses on escape if they are not loaded already
mod_esc_1.1.upd <- readRDS("outputs/models_outputs/mod_esc_1.1.upd.rds")
mod_esc_2.1.upd <- readRDS("outputs/models_outputs/mod_esc_2.1.upd.rds")
mod_esc_3.1.upd <- readRDS("outputs/models_outputs/mod_esc_3.1.upd.rds")
mod_esc_4.1.upd <- readRDS("outputs/models_outputs/mod_esc_4.1.upd.rds")

# Focal and social resource acquisition

# Select values for the surface
genf <- seq(min(df$Zgen), max(df$Zgen), length.out = 10)
gens <- seq(min(df$Zsocial_gen), max(df$Zsocial_gen), length.out = 10)

# Compute the z axis values
z1 <- outer(genf, 
            gens, 
            FUN = function(x, y) {plogis(fixef(mod_esc_4.1.upd)[1] + 
                                           (fixef(mod_esc_4.1.upd)[2] * x) +
                                           (fixef(mod_esc_4.1.upd)[6] * y) + 
                                           (fixef(mod_esc_4.1.upd)[10] * I(x^2)) +
                                           (fixef(mod_esc_4.1.upd)[14] * I(y^2)) + 
                                           (fixef(mod_esc_4.1.upd)[31] * x * y))
            }
)


# Plot
plot1 <- plot_ly(x = ~genf, 
                 y = ~gens, 
                 z = ~z1,
                 type = "contour",
                 coloraxis = "coloraxis",
                 autocontour = F,
                 contours = list(start = 0,
                                 end = 1,
                                 size = 0.1)) %>%
  
  layout(xaxis = list(title = "Focal resource acquisition"),
         yaxis = list (title = "Social resource acquisition"))


# Focal and social hiding

# Select values for the surface
crouchf <- seq(min(df$Zcrouching), max(df$Zcrouching), length.out = 10)
crouchs <- seq(min(df$Zsocial_crouching), max(df$Zsocial_crouching), length.out = 10)

# Compute the z axis values
z2 <- outer(crouchf, 
            crouchs, 
            FUN = function(x, y) {plogis(fixef(mod_esc_4.1.upd)[1] + 
                                           (fixef(mod_esc_4.1.upd)[3] * x) +
                                           (fixef(mod_esc_4.1.upd)[7] * y) + 
                                           (fixef(mod_esc_4.1.upd)[11] * I(x^2)) +
                                           (fixef(mod_esc_4.1.upd)[15] * I(y^2)) + 
                                           (fixef(mod_esc_4.1.upd)[34] * x * y))
            }
)

# Plot
plot2 <- plot_ly(x = ~crouchf, 
                 y = ~crouchs, 
                 z = ~z2,
                 type = "contour",
                 coloraxis = "coloraxis",
                 autocontour = F,
                 contours = list(start = 0,
                                 end = 1,
                                 size = 0.1)) %>%
  
  layout(xaxis = list(title = "Focal hiding"),
         yaxis = list(title = "Social hiding"))


# Focal and social hiding

# Select values for the surface
chasef <- seq(min(df$Zchase), max(df$Zchase), length.out = 10)
chases <- seq(min(df$Zsocial_chase), max(df$Zsocial_chase), length.out = 10)

# Compute the z axis values
z3 <- outer(chasef, 
            chases, 
            FUN = function(x, y) {plogis(fixef(mod_esc_4.1.upd)[1] + 
                                           (fixef(mod_esc_4.1.upd)[4] * x) +
                                           (fixef(mod_esc_4.1.upd)[8] * y) + 
                                           (fixef(mod_esc_4.1.upd)[12] * I(x^2)) +
                                           (fixef(mod_esc_4.1.upd)[16] * I(y^2)) + 
                                           (fixef(mod_esc_4.1.upd)[32] * x * y))
            }
)


# Plot
plot3 <- plot_ly(x = ~chasef, 
                 y = ~chases, 
                 z = ~z3,
                 type = "contour",
                 coloraxis = "coloraxis",
                 autocontour = F,
                 contours = list(start = 0,
                                 end = 1,
                                 size = 0.1)) %>%
  
  layout(xaxis = list(title = "Focal defense"),
         yaxis = list(title = "Social defense"))


# Focal and social rescue

# Select values for the surface
unhookf <- seq(min(df$Zunhook), max(df$Zunhook), length.out = 10)
unhooks <- seq(min(df$Zsocial_unhook), max(df$Zsocial_unhook), length.out = 10)

# Compute the z axis values
z4 <- outer(chasef, 
            chases, 
            FUN = function(x, y) {plogis(fixef(mod_esc_4.1.upd)[1] + 
                                           (fixef(mod_esc_4.1.upd)[5] * x) +
                                           (fixef(mod_esc_4.1.upd)[9] * y) + 
                                           (fixef(mod_esc_4.1.upd)[13] * I(x^2)) +
                                           (fixef(mod_esc_4.1.upd)[17] * I(y^2)) + 
                                           (fixef(mod_esc_4.1.upd)[33] * x * y))
            }
)

# Plot
plot4 <- plot_ly(x = ~unhookf, 
                 y = ~unhooks, 
                 z = ~z4,
                 type = "contour",
                 coloraxis = "coloraxis",
                 autocontour = F,
                 contours = list(start = 0,
                                 end = 1,
                                 size = 0.1)) %>%
  
  layout(xaxis = list(title = "Focal rescue"),
         yaxis = list(title = "Social rescue"))

#plots together in one figure
#https://plotly.com/r/subplots/

interactions_plot <- subplot(plot1, plot2, plot3, plot4,
                             nrows = 2,
                             titleY = TRUE,
                             titleX = TRUE,
                             margin = 0.07 ) %>%
  
  layout(title = 'Interactions between focal and social behaviours - Survival',coloraxis = list(colorscale = 'Viridis'))

interactions_plot 

#save it by hand as there is no package to export fixed image of plotly at the moment


#######################
# 5. Figure S2 - P and CI matrix
######################

#P 
#read in the P matrix created with the main models script
P<- fread("outputs/tables/P.csv")
P$V1<- NULL

P_mat<-as.matrix(P)

#replace the diagonal with 1
diag(P_mat)<-1

colnames(P_mat) = c('Resource acquisition', 'Hiding', 'Defense','Rescue')
rownames(P_mat) = c('Resource acquisition', 'Hiding', 'Defense','Rescue')

png(height=1800, width=1800, 
    file = "outputs/models_plots/PI_mat.png",
    type = "cairo")

corrplot(P_mat, type = "lower", method = "ellipse", 
         cl.pos = "r", cl.cex = 4,
         tl.pos = "lt", tl.col = "black", 
         tl.cex = 4, tl.srt = 45, number.digits = 3,
         addCoef.col = "black", 
         col = brewer.pal(n = 10, name = "RdBu"), 
         mar = c(0,0,0,0),
         number.cex = 4
)

dev.off()


#CI
#read in the CI matrix created with the main models script
CI <- fread("outputs/tables/CI.csv")
CI[,1]<- NULL

CI_mat<-as.matrix(CI)

colnames(CI_mat) = c('Social Resource acquisition', 'Social Hiding', 'Social Defense','Social Rescue')
rownames(CI_mat) = c('Focal Resource acquisition', 'Focal Hiding', 'Focal Defense','Focal Rescue')


png(height=1800, width=1800, 
    file = "outputs/models_plots/CI_mat.png",
    type = "cairo")

corrplot(CI_mat, is.corr = TRUE, type = "full", method = "ellipse", 
         cl.pos = "r", cl.cex = 4,
         tl.pos = "lt", tl.col = "black", 
         tl.cex = 4, tl.srt = 45, number.digits = 3,
         addCoef.col = "black", 
         col = brewer.pal(n = 10, name = "RdBu"), 
         mar = c(0,0,0,0),
         number.cex = 4)

dev.off()

rm(list = ls())
library("vars")
library("urca")
lrm(list = ls())
library("vars")
library("urca")
library("tseries")
library("timeSeries")
library("lmtest")
library(ggplot2)
library(gridExtra)
library(scales)
library(grid)

############################
## I.Préparation Des Données##
############################
mainpath <- "~/Google Drive/HEC/projet Macro/"
pathdata <- paste(mainpath,"Data/", sep="")
pathgraphs <- paste(mainpath,"Graphs/", sep="")
pathgraphs_analysis <- paste(mainpath,"redaction/graphs", sep="")

# Charger données #

setwd(pathdata)
getwd()
library(readxl)
data <- read_excel("data_suisse.xlsx", sheet = "Feuille 2")

summary(data)

#############################
## Préparation des Données ##
#############################

#######################
## Données en niveau ##
#######################
PIB_real <- data$`PIB réel
plot_PIB_real <- c(plot(PIB_real, main="PIB Réel"),lines(PIB_real))


log_PIB_real <- log(PIB_real)
plot_log_PIB_real <- c(plot(log_PIB_real, main = "Log PIB Réel"),lines(log_PIB_real))


r <- data$`Taux intérêt`
plot_r <- c(plot(r, main = "Taux d'intérêt"),lines(r))

################################
## Données en 1ère Différence ##
################################

D_PIB_real <- diff(log_PIB_real, lag = 1, trim=TRUE)
plot_D_PIB_real <- c(plot(D_PIB_real, main = "1ère Différence log PIB Réel"),lines(D_PIB_real))

##############################
## II.Test de Stationnarité ##
##############################

####################################
## Test de Dickey-Fuller Augmenté ##
####################################

# H0: le processus est non-stationnaire #
# Non-rejet de H0: t-stat >= t_tabulé

# test de stationnaire à seuil de 5%#
##########
# niveau #
##########
# modèle avec constante et dérive temporelle #
summary(ur.df(PIB_real, type="trend", selectlags = "AIC"))# Bon #
summary(ur.df(PIB_real, type="drift", selectlags = "AIC"))# Pas Bon #
T_log_PIB_t <- summary(ur.df(log_PIB_real, type="trend", selectlags = "AIC"))# Bon #
T_log_PIB_t@teststat
T_log_PIB_t@cval
T_log_PIB_d <- summary(ur.df(log_PIB_real, type="drift", selectlags = "AIC"))# Pas Bon #
T_log_PIB_d@teststat
T_log_PIB_d@cval

T_r_t <- summary(ur.df(r, type="trend", selectlags = "AIC")) #Bon + Stationnaire #
T_r_t@teststat
T_r_t@cval
T_r_d <- summary(ur.df(r, type="drift", selectlags = "AIC"))
T_r_d@teststat
T_r_d@cval
T_r_n <- summary(ur.df(r, type="none", selectlags = "AIC"))
T_r_n@teststat
T_r_n@cval

#######################
# Différence Première #
#######################
T_D_log_PIB_t <- summary(ur.df(D_PIB_real, type="trend", selectlags = "AIC"))# BON + stationnaire #
T_D_log_PIB_t@teststat
T_D_log_PIB_t@cval
T_D_log_PIB_d <- summary(ur.df(D_PIB_real, type="drift", selectlags = "AIC"))
T_D_log_PIB_d@teststat
T_D_log_PIB_d@cval
T_D_log_PIB_n <- summary(ur.df(D_PIB_real, type="none", selectlags = "AIC"))
T_D_log_PIB_n@teststat
T_D_log_PIB_n@cval

# tau3 = parameter of lags -> Stationarity si different de 0
# phi2 = paramètre de la constante  
# phi3 = paramètre de la tendance

########################################
## III. Analyse d'un VAR et d'un SVAR ##
########################################

#############################
## Préparation Des Données ##
#############################

attach(data)
data_var_n <- as.data.frame(cbind(PIB_real, r))
colnames(data_var) <- c("PIB réel","Taux intérêt")
varnames <- c("PIB réel","Taux intérêt")

data_var_log <- as.data.frame(cbind(log_PIB_real, r))
colnames(data_var_log) <- c("PIB réel","Taux intérêt")
varnames <- c("PIB réel","Taux intérêt")

data_var <- as.data.frame(cbind(D_PIB_real, r))# c'est celui qu'il faut choisir #
colnames(data_var) <- c("D_PIB_real", "r")
varnames <- c("D_PIB_real", "r")


####################
## Estimation VAR ##
####################

# Lag selection #

VARselect(data_var, lag.max = 5)

# Estimation #

VAR <- VAR(data_var, p = 1)
VAR


######################
## Fonctions utiles ##
######################

## fonction estimant un VAR ##

IRF_data_VAR <- function(data_var = NULL, p=4, ci=0.95) {
  VAR_est <- VAR(data_var,p=p) 
  IRF <- irf(VAR_est, ci = ci)
  return(IRF)
}

## fonction estimant un SVAR ##

IRF_data_SVAR <- function(data_var = NULL, p=1, ci=0.95, BM = matrix(c(1,0,1,NA), ncol = 2, byrow = T)) {
  VAR_est <- VAR(data_var,p=p) 
  SVAR_est <- SVAR(VAR_est, Bmat = BM, lrtest = FALSE)
  IRF <- irf(SVAR_est, ci = ci)
  return(IRF)
}

## fonction pour créer des plots IRF ##

plot_irf <- function(data_irf = NULL, var_shock = "", var_reac = "", col_ic = "blue", alpha = 0.5, ylim= c(NA,NA)){
  xtitle <- var_reac
  if (var_reac == 'D_PIB_real') {
    xtitle <- 'D_PIB_real'
  }
  if (var_reac == 'r') {
    xtitle <- 'r'
  }
  name <- paste(var_shock,".",var_reac, sep = "")
  col_IC <- adjustcolor(col_ic, alpha.f = alpha) 
  irf <- as.data.frame(data_irf$irf)[, name]
  lower <- as.data.frame(data_irf$Lower)[, name]
  upper <- as.data.frame(data_irf$Upper)[, name]
  plot_data <- as.data.frame(cbind(1:length(irf),irf,lower,upper))
  colnames(plot_data) <- c("lag", "irf", "lower", "upper")
  if (is.na(ylim[2])){
    ylim[2] <- max(upper)*1 + abs(min(upper)*0.5)
  }
  if (is.na(ylim[1])){
    ylim[1] <- min(lower)*1 - abs(min(lower)*0.5)
  }
  ymax <- ylim[2]
  ymin <- ylim[1]
  p1 <- ggplot(plot_data, aes(lag, irf)) +
    theme_bw() +
    theme(plot.title = element_text(size=9), 
          axis.line = element_line(colour = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.text.x = element_text(size = 6.5),
          axis.text.y = element_text(size = 6.5), 
          plot.margin = unit(c(0,0,-0.1,-0.1), "cm")) + 
    #geom_line(aes(x=lag,y=lower), linetype = "dashed") +
    #geom_line(aes(x=lag,y=upper), linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dotted") + 
    geom_ribbon(data=plot_data,aes(ymin=lower,ymax=upper),alpha=alpha, fill=col_ic) +
    geom_line() +
    xlab("") +
    ylab("") + 
    coord_cartesian(ylim = c(ymin, ymax)) +
    ggtitle(xtitle) 
  p1
}



################
## SVAR Model ##
################

## test mais pas utiliser dans le document rendu ##
attach(data_var_demid)
time <- 1:length(D_PIB_real)

BM = matrix(c(1,0,1,NA), ncol = 2, byrow = T)
IRF_SVAR <- IRF_data_SVAR(data_var_demid, p=1,ci=0.8)
IRF_SVAR
plot_irf(IRF_SVAR, var_shock = "r", var_reac = "D_PIB_real")

#Prévision faite par notre VAR#
predict(VAR)

# Matrice Variance-Covraiance#
vcov(VAR)

BQ(VAR)

# Decomposition de la variance #
fevd(VAR)

################
## Graphe IRF ##
################

plot((IRF_data_VAR(data_var_n)))
plot(IRF_data_VAR(data_var_log))
plot(IRF_data_VAR(data_var))# prendre celui la


p1 <- plot_irf(IRF_data_VAR(data_var), var_shock = "r", var_reac = "D_PIB_real")
p2 <- plot_irf(IRF_data_VAR(data_var), var_shock = "r", var_reac = "r")
pdf("irf_r.pdf")
irf_r <- grid.arrange(arrangeGrob(p1, p2, ncol=1), top = textGrob("Choc Taux d'intérêt",gp=gpar(fontsize=20,font=3)))
dev.off()

p3 <- plot_irf(IRF_data_VAR(data_var), var_shock = "D_PIB_real", var_reac = "D_PIB_real")
p4 <- plot_irf(IRF_data_VAR(data_var), var_shock = "r", var_reac = "r")
pdf("irf_D_PIB_real.pdf")
irf_PIB_real <- grid.arrange(arrangeGrob(p3, p4, ncol=1), top = textGrob("Choc PIB réel",gp=gpar(fontsize=20,font=3)))
dev.off()

setwd(pathgraphs)
pdf("irf_total.pdf")
irf_demid_total <- grid.arrange(irf_r, irf_PIB_real, ncol= 2, nrow= 1)
dev.off()

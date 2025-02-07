######################################################
# Wavelet Quantile-on-Quantile Regression (WQQR) With Pvalue
# Date: 1/01/2025
# Authors: Tomiwa Sunday Adebayo, Oktay Özkan, Dilber Uzun Ozsahin, Babatunde Sunday Eweade and Bright Akwasi Gyamfi
# Twaikline@gmail.com
# https://enveurope.springeropen.com/articles/10.1186/s12302-025-01059-z
# Journal Name: Environmental Sciences Europe
#######################################################
rm(list=ls(all=TRUE))

# set working directory for interactive run


library(readxl)
library(quantreg)
library(KernSmooth)
library(np)
library(waveslim)
library(ggplot2)
library(WaveletComp)
library("writexl")
library(plotly)
library("parallel")
options(warn=-1)
options(mc.cores=detectCores())

# Loading DATA, where it is located and name of file.

dataset <- read_excel("DATA.xlsx")

attach(dataset)

DEP <- CO2
IND <- ICT

d1=as.data.frame(mra(DEP, wf = "la8", J = 5, method = "modwt", boundary = "periodic"))
d2=as.data.frame(mra(IND, wf = "la8", J = 5, method = "modwt", boundary = "periodic"))

y <- d1[,4]+d1[,5]
x <- d2[,4]+d2[,5]

## for medium-run use d1[,3] and d2[,3]
## for long-run use d1[,4]+d1[,5] and d2[,4]+d2[,5]

# Set quantiles
quantiles <- 0.05 # Set the quantile interval
h <- 1 
num <- 1 / quantiles - 1


# Apply QR and get the coefficents
t <- nrow(dataset)
yt <- y[1:t]
xt <- x[1:t]
tau <- seq(from = quantiles,
            to = 1 - quantiles,
            by = quantiles)
QR.coef <- summary(rq(yt ~ xt, tau = tau, ci = F), se = "iid") # you can add control variables here
for (i in 1:num){
  QR.coef[[i]] <- QR.coef[[i]]$coefficients
}
QR.pval <- rep(NA, num)
for (i in 1:num){
  QR.pval[i] <- QR.coef[[i]][2,4]
}
QR.coef1 <- rep(NA, num)
for (i in 1:num){
  QR.coef1[i] <- QR.coef[[i]][2,1]
}
QR.coef <- QR.coef1
rm(QR.coef1)

QR_coef=data.frame(QR.coef)
write_xlsx(QR_coef," QR_coef.xlsx")

QR_pval=data.frame(QR.pval)
write_xlsx(QR_pval," QR_pval.xlsx")

# Add QQR Function
lprq <- function(x, y, m = num, tau = .5) {
  yt <- y[1:t]
  xt <- x[1:t]
  xx <- seq(min(xt), max(xt), length = m)
  pv <- xx
  dv <- xx
  Fn <- xt
  for (i in 1:length(xt)) {
    Fn[i] <- length(which(xt < xt[i])) / length(xt)
  }
  
  for (i in 1:length(xx)) {
    zt <- xt - xx[i]
    wx <-
      dnorm(zt / h) # solve the gaussian kernal as the weight
    r <-
      summary(rq(yt ~ zt, # you can add control variables here
                 weights = wx,
                 tau = tau,
                 ci = F), se = "iid")  # FALSE method = "br",
    pv[i] <- r$coefficients[2, 4]
    dv[i] <- r$coefficients[2, 1]
  }
  list(xx = xx, pv = pv, dv = dv)
}

#Create a matrix to save the QQR estimates
QQR.coef <- as.data.frame(matrix(0, ncol = num, nrow = num))
QQR.pval <- as.data.frame(matrix(0, ncol = num, nrow = num))

# Run QQR in a loop and save estimates in matrix "QQR.coef"
#Note: 9 (origanlly 0.05) in below loop is the bandwidth that can be adjusted
## which can be obtaned from the long-run variance estimates (in eviews) ...
for (i in 1:num) {
  res <- lprq(x = x,
              y = y,
              m = num,
              tau = tau[i])
  QQR.coef[, i] <- res$dv  ## "SJ","Silverman","YJ","CV"
  QQR.pval[, i] <- res$pv
  print(i)
}

#Save the QQR results

QQR_coef=data.frame(QQR.coef)
write_xlsx(QQR_coef, "QQR_coef.xlsx")

QQR_pval=data.frame(QQR.pval)
write_xlsx(QQR_pval, "QQR_pval.xlsx")

# HEATMAP

library(viridisLite)
library(ggplot2)
library(scales)
library(tidyverse)
library(dplyr)
library(tidyr)
library(reshape2)

QQR <- as.matrix(QQR_coef)
transpoze_QQR <- t(QQR)

tau_formatted <- sprintf("%.2f", tau)  

p <- plot_ly(z=~transpoze_QQR, x=~tau, y=~tau)%>% 
add_surface(colorscale='Rainbow', contours = list(
          x = list(show = TRUE, color = "black", width = 1, start = 0,
                                end = 1, size = 0.05),                               
          y = list(show = TRUE, color = "black", width = 1, start = 0,
                                end = 1, size = 0.05 
                                ))) %>%
  layout(
    showlegend = F,
    scene = list(
      xaxis = list(nticks=10, range=c(0,1), ticks="outside", showline= T, autorange = "reversed",title = 'lnFD'), 
      yaxis = list(nticks=10, range=c(0,1), ticks="outside", showline= T, autorange = "reversed",title = 'lnEP'),
      zaxis = list(showline= T, title = ''),
    aspectratio = list(x = 0.70, y = 0.70, z = 0.70)  # E?it boyutlar i?in aspectratio ekleyin
  )
) %>% 

 colorbar(len=0.35, x=0.73, thickness=8, title = "", y=0.5)      
      
p



veri <- as.matrix(QQR_pval)
veri[veri > 0.05] <- 1
veri[veri <= 0.05] <- 0

veri <- as.data.frame(t(veri))

veris <- veri %>% 
  pivot_longer(cols = everything(), names_to = "Sutunlar", values_to = "Degerler")


data <- expand.grid(X = tau, Y = tau)
data$Z <- veris$Degerler
data <- data[complete.cases(data), ]

pallete <- c("firebrick1","white")

p <- ggplot(data, aes(x = as.factor(X), y = as.factor(Y), fill = Z)) +
  geom_tile(color = "steelblue3", size = 0.125) +
  scale_x_discrete(breaks = as.character(tau), labels = sprintf("%.2f", tau)) +
  scale_y_discrete(breaks = as.character(tau), labels = sprintf("%.2f", tau)) +  
  scale_fill_gradientn(colours = pallete, limits=c(0,1)) +
  labs(x = "CO2", y = "ICT", fill = "") +
  theme(axis.text.x = element_text(angle = 90))

p <- p + theme(
    legend.position = "none"
  )

print(p)


# WQR vs AWQQR 

df <- data.frame(tau = as.factor(tau), results = QR.coef)
ggplot(df, aes(x = tau, y = results, group = 1.25)) +
  geom_line(color = "firebrick1", size=1.0) +
  geom_line(aes(y = colMeans(QQR.coef)), linetype = "dashed", color = "steelblue3", size=1.25) +
  labs(x = "Quantiles", y = "Slope Coefficients")+
  theme(axis.text.x = element_text(angle = 45))



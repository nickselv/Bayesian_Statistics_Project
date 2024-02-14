
##### MODELING c_i

#install spatio-tmemporal model packages
install.packages("CARBayes")
install.packages("CARBayesST")

#load needed packages
library(CARBayes)
library(sf)
library(spdep)
library(sfheaders)
library(dplyr)
library(CARBayesST)
library(coda)
library(ggplot2)

#load dataset
data_CAR22 <- read.csv("~/Bayesian statistics/Project/data_CAR22.csv")

#DATA EXPLORATION: 

#transform the dataset into a spatial dataset (long and lat as geographical coordinates)
data_sf <- st_as_sf(data_CAR22, coords = c("longitude", "latitude"))
#setting crs (coordinate reference system)
data_sf <- st_set_crs(data_sf, 4326) #for latitude and longitude we use WGS84
#manipulation
data_sf$id_varco <- as.factor(data_sf$id_varco)
data_CAR22$id_varco <- as.factor(data_CAR22$id_varco)
data_CAR22$data <- rep(1:365, each=40)

#average entries in each gate plot
data4plot <- data_sf[1:40,]
ingressi.avg <- summarise(group_by(data_stan, id_varco), 
                          ingressi.mean=mean(numero_transiti))
data4plot$numero_transiti <- ingressi.avg$ingressi.mean

library(leaflet)
colours <- colorNumeric(palette = "YlOrRd", domain = data4plot$numero_transiti)
leaflet(data = data4plot) %>%
  addTiles() %>%
  addCircleMarkers(fillColor = ~colours(data4plot$numero_transiti), 
                   radius = data4plot$numero_transiti/500,
                   fillOpacity = 0.8, color = "gray") %>%
  addLegend(pal = colours, values = data4plot$numero_transiti, opacity = 1, 
            title = "passages") %>%
  addScaleBar(position = "bottomleft")

#SPATIO-TERMPORAL MODEL for c_i:

#proximity matrix W
#define 3 different proximity matrices

#1) we assign w=1 if the gates are adjacent, 0 otherwise
W1 <- matrix(0, 40, 40)
W1[1,c(2,40)] <- c(1,1) #first row
W1[40,c(1,39)] <- c(1,1) #last row
for (i in 2:39) {
  W1[i,c(i-1,i+1)] <- c(1,1)
}

#2) we assign decreasing values starting from 1 for adjacent gates, 0.5, 0.25 etc
W2 <- matrix(0, 40, 40)
w <- c()
w[1] <- 0
for (i in 1:20) {
  w[i+1] <- 1/2^(abs(i-1)) 
}
w <- c(w, rev(w[2:20]))
W2[1,] <- w #first row
for (i in 2:40) {
  W2[i,] <- c(w[(41-i+1):40],w[1:(41-i)])
}

#3) we assign values according to the distance between gates 
dgates <- function(long1, lat1, long2, lat2) {
  d <- sqrt((lat1-lat2)^2 + (long1-long2)^2)
  return(d)
}
data4dist <- data_CAR22[1:40,c(1,6,7)]
dist <- matrix(0, nrow = 40, ncol = 40)
for (i in 1:39) {
  for (j in (i+1):40) {
    long1 <- data4dist[i,2]
    lat1 <- data4dist[i,3]
    long2 <- data4dist[j,2]
    lat2 <- data4dist[j,3]
    dist[i,j] <- dgates(long1, lat1, long2, lat2)
    dist[j,i] <- dist[i,j]
  }
}
W3 <- matrix(nrow = 40, ncol = 40)
for (i in 1:40) {
  for(j in 1:40) {
    if(i!=j)
      W3[i,j] <- 1/(1000*dist[i,j])
    else
      W3[i,j] <- 0
  }
}

#dataset further preparation
ingressi.avg <- summarise(group_by(data_CAR22, id_varco), 
                          ingressi.avg=mean(numero_transiti))
ingressi.avg <- rep(ingressi.avg$ingressi.avg, each=304)
data.reg1 <- cbind(data_CAR22, ingressi.avg)
colnames(data.reg1)[4] <- "avg_transiti"

#Poisson models (one for each W)
#1) W1
st.model1.1 <- ST.CARar(numero_transiti ~ offset(log(avg_transiti)), 
                        family = "poisson", data = data.reg1, W = W1, 
                        burnin = 1000, n.sample = 5000, AR = 1)
print(st.model1.1)
st.model1.1$modelfit

data4ci <- cbind(data_CAR22[,1:2], st.model1.1$fitted.values)
colnames(data4ci)[3] <- "c_i.poi1"
c_i.df1.poi <- summarise(group_by(data4ci, id_varco), c_i.mean=mean(c_i.poi1))
colnames(c_i.df1.poi)[2] <- "c_i.poi1"
hist(c_i.df1.poi$c_i.poi1, probability = TRUE, col = "lightblue", 
     main = "c_i estimates distribution", xlab = "c_i estimates")

#2) W2
st.model2.1 <- ST.CARar(numero_transiti ~ offset(log(avg_transiti)), 
                        family = "poisson", data = data.reg1, W = W2, 
                        burnin = 1000, n.sample = 5000, AR = 1)
print(st.model2.1)
st.model2.1$modelfit

data4ci <- cbind(data_CAR22[,1:2], st.model2.1$fitted.values)
colnames(data4ci)[3] <- "c_i.poi2"
c_i.df2.poi <- summarise(group_by(data4ci, id_varco), c_i.mean=mean(c_i.poi2))
colnames(c_i.df2.poi)[2] <- "c_i.poi2"
hist(c_i.df2.poi$c_i.poi2, probability = TRUE, col = "lightblue", 
     main = "c_i estimates distribution", xlab = "c_i estimates")

#3) W3
st.model3.1 <- ST.CARar(numero_transiti ~ offset(log(avg_transiti)), 
                        family = "poisson", data = data.reg1, W = W3, 
                        burnin = 1000, n.sample = 5000, AR = 1)
print(st.model3.1)
st.model3.1$modelfit

data4ci <- cbind(data_CAR22[,1:2], st.model3.1$fitted.values)
colnames(data4ci)[3] <- "c_i.poi3"
c_i.df3.poi <- summarise(group_by(data4ci, id_varco), c_i.mean=mean(c_i.poi3))
colnames(c_i.df3.poi)[2] <- "c_i.poi3"
hist(c_i.df3.poi$c_i.poi3, probability = TRUE, col = "lightblue", 
     main = "c_i estimates distribution", xlab = "c_i estimates")

#Gaussian models with standardization (one for each W)
#instead of using the Poisson we use the Gaussian distribution but first the 
#number of entries are standardized
ingressi.sd <- summarise(group_by(data_CAR22, id_varco), ingressi.sd=sd(numero_transiti))
ingressi.sd <- rep(ingressi.sd$ingressi.sd, each=304)
transiti.std <- (data_CAR22$numero_transiti - data.reg1$avg_transiti)/ingressi.sd
hist(transiti.std, probability = TRUE, col = "lightblue", xlab = " ",
     main = "Standardized entries", xlim = c(-3, 3))
data.reg1 <- cbind(data.reg1, transiti.std)

#1) W1
st.model1.3 <- ST.CARar(transiti.std ~ 1, family = "gaussian", data = data.reg1, 
                        W = W1, burnin = 1000, n.sample = 2000, AR = 1)
print(st.model1.3)
st.model1.3$modelfit

data4ci <- cbind(data_CAR22[,1:2], st.model1.3$fitted.values)
colnames(data4ci)[3] <- "c_i.gauss1"
c_i.df1 <- summarise(group_by(data4ci, id_varco), c_i.mean=mean(c_i.gauss1))
colnames(c_i.df1)[2] <- "c_i.gauss1"
hist(c_i.df1$c_i.gauss1, probability = TRUE, col = "lightblue", 
     main = "c_i estimates distribution", xlab = "c_i estimates")

#2) W2
st.model2.3 <- ST.CARar(transiti.std ~ 1, family = "gaussian", data = data.reg1, 
                        W = W2, burnin = 1000, n.sample = 2000, AR = 1)
print(st.model2.3)
st.model2.3$modelfit

data4ci <- cbind(data_CAR22[,1:2], st.model2.3$fitted.values)
colnames(data4ci)[3] <- "c_i.gauss2"
c_i.df2 <- summarise(group_by(data4ci, id_varco), c_i.mean=mean(c_i.gauss2))
colnames(c_i.df2)[2] <- "c_i.gauss2"
hist(c_i.df2$c_i.gauss2, probability = TRUE, col = "lightblue", 
     main = "c_i estimates distribution", xlab = "c_i estimates")

#3) W3
st.model3.3 <- ST.CARar(transiti.std ~ 1, family = "gaussian", data = data.reg1, 
                        W = W3, burnin = 1000, n.sample = 2000, AR = 1)
print(st.model3.3)
st.model3.3$modelfit

data4ci <- cbind(data_CAR22[,1:2], st.model3.3$fitted.values)
colnames(data4ci)[3] <- "c_i.gauss3"
c_i.df3 <- summarise(group_by(data4ci, id_varco), c_i.mean=mean(c_i.gauss3))
colnames(c_i.df3)[2] <- "c_i.gauss3"
hist(c_i.df3$c_i.gauss3, probability = TRUE, col = "lightblue", 
     main = "c_i estimates distribution", xlab = "c_i estimates")

#create c_i estimates dataset to export in Python
c_i.df <- cbind(c_i.df1[,1], c_i.df1[,2], c_i.df1.poi[,2], c_i.df2[,2], 
                c_i.df2.poi[,2], c_i.df3[,2], c_i.df3.poi[,2])
write.csv(c_i.df, file = "C:/Users/HP//Bayesian statistics/Project/c_i22.csv", 
          row.names = FALSE)

#NOTE: only Poisson c_i will be used 

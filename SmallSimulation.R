#Small Sim 

SmallSim <- function(scale, a, b, c){
  coordx <- rep(1:20,20)
  coordy <- c()
  for(i in 1:20){
    coordy <- c(coordy, rep(i,20))
  }
  coordinates = cbind(coordx, coordy)
  library("CompRandFld")
  library(RandomFields)
  #and scale 0.1 and sill = 1
  A = RFsim(coordx = coordx, coordy = coordy, corrmodel = "exponential", model = "BinaryGauss", threshold = 0, param = list(nugget = 0, mean = 0, scale = scale, sill = 0.6))
  #A = matrix(A$data, nrow = 50, ncol = 50, byrow = T)
  #MEdiation simulation
  
  #correlation isn't so high 0.2, maybe collinearity shouldn't be a huge problem in that scenario
  model <- "spherical"
  parameters <- list(mean = 0.5, sill = 1, nugget = 0.5, scale = scale)
  Con = RFsim(coordx = coordx, coordy = coordy, corrmodel = model, model = "Gaussian", param = parameters)
  #Simulate the A dependent part 
  M_A <- a * A$data + 1.5*Con$data
  #Simulation the covariance function
  #correlation options, cauchy, exponential, gauss, gencauchy, spherical, stable, matern, wave
  model <- "spherical"
  
  CorrelationParam("spherical")
  parameters <- list(mean = 0, sill = 1, nugget = 0, scale = scale)
  
  M_cov = RFsim(coordx = coordx, coordy = coordy, corrmodel = model, model = "Gaussian", param = parameters)
  #M_cov = matrix(M_cov$data, nrow = 50,ncol = 50, byrow = T)
  
  #Final M 
  
  M <- M_A + M_cov$data
  
  
  #Simulate part of Y dependent on M and A 
  
  Y_AM = b * A$data + c * M + 1.0*Con$data + ((b+c)/2)*A$data*M
  
  #Covariance part of Y 
  model = "exponential"
  CorrelationParam(model)
  parameters <- list(mean = 0, sill = 1, nugget = 0, scale = scale)
  
  Y_cov = RFsim(coordx = coordx, coordy = coordy, corrmodel = model, model = "Gaussian", param = parameters)
  #Y_cov = matrix(Y_cov$data, nrow = 50,ncol = 50, byrow = T)
  
  Y = Y_AM + Y_cov$data
  
  Ydf <- as.data.frame(cbind(coordinates, Y))
  names(Ydf) <- c("X", "Y", "data_resp")
  library(raster)
  # create spatial points data frame
  spg <- Ydf
  coordinates(spg) <- ~X + Y
  # coerce to SpatialPixelsDataFrame
  gridded(spg) <- TRUE
  # coerce to raster
  Yraster <- raster(spg)
  
  Mdf <- as.data.frame(cbind(coordinates, M))
  names(Mdf) <- c("X", "Y", "data_M")
  library(raster)
  # create spatial points data frame
  spg <- Mdf
  coordinates(spg) <- ~X + Y
  # coerce to SpatialPixelsDataFrame
  gridded(spg) <- TRUE
  # coerce to raster
  Mraster <- raster(spg)
  
  Adf <- as.data.frame(cbind(coordinates, A$data))
  names(Adf) <- c("X", "Y", "data_A")
  library(raster)
  # create spatial points data frame
  spg <- Adf
  coordinates(spg) <- ~X + Y
  # coerce to SpatialPixelsDataFrame
  gridded(spg) <- TRUE
  # coerce to raster
  Araster <- raster(spg)
  
  Cdf <- as.data.frame(cbind(coordinates, Con$data))
  names(Cdf) <- c("X", "Y", "data_C")
  library(raster)
  # create spatial points data frame
  spg <- Cdf
  coordinates(spg) <- ~X + Y
  # coerce to SpatialPixelsDataFrame
  gridded(spg) <- TRUE
  # coerce to raster
  Craster <- raster(spg)
  
  raster_stack <- stack(Yraster, Araster, Mraster, Craster)
  return(raster_stack)
}

S1 = SmallSim(1, 0.1012, .6, .252)
S2 = SmallSim(3, 0.1012, .6, .252)
S3 = SmallSim(10, 0.1012, .6, .252) 
S4 = SmallSim(1, 0.3, 0.4, 0.75)
S5 = SmallSim(3, 0.3, 0.4, 0.75)
S6 = SmallSim(10, 0.3, 0.4, 0.75)
S7 = SmallSim(1, .32, 6.0, 0.8)
S8 = SmallSim(3, .32, 6.0, 0.8)
S9 = SmallSim(10, .32, 6.0, 0.8)
S10 = SmallSim(1, 0.94, 4.0, 2.35)
S11 = SmallSim(3, 0.94, 4.0, 2.35)
S12 = SmallSim(10, 0.94, 4.0, 2.35)

SimulationAnalysis <- function(raster_stack){
  
  df_data <- data.frame(rasterToPoints(raster_stack))
  coords <- df_data[,c(1,2)]
  library("spmoran")
  start = Sys.time()
  #Moran SVC
  y_mod <- besf_vc(df_data$data_resp, df_data[,c("data_A", "data_M", "data_C")], coords = coords, x_sel = FALSE)
  rmse_svc = sqrt(mean((y_mod$pred - df_data$data_resp) ^ 2, na.rm = T))
  
  Avar <- df_data[,c(4,6)]
  Mvar <- df_data$data_M
  m_mod = besf_vc(df_data$data_M, Avar, coords = coords,x_sel = FALSE)
  
  rmse_m_svc = sqrt(mean((m_mod$pred - df_data$data_M) ^ 2, na.rm = T))
  
  df_data$NDE = y_mod$b_vc[,2]
  
  df_data$NIE = y_mod$b_vc[,3] * m_mod$b_vc[,2]
  
  end_svc <- Sys.time() - start # 1.
  
  data_svc=rasterFromXYZ(df_data)
  library("spgwr")
  start <- Sys.time()
  mediation_spatialpoints <- rasterToPoints(raster_stack, spatial = T)
  med_spdf <- as(mediation_spatialpoints, "SpatialPointsDataFrame")
  col.bw <- gwr.sel(data_resp ~ data_A + data_M + data_C, data = med_spdf)
  
  gwr_outcome <- gwr(data_resp ~ data_A + data_M + data_C + data_A*data_M, data = med_spdf, bandwidth = col.bw)
  #predicted_surface = spplot(gwr_outcome$SDF, "pred")
  #true_sim_outcome = spplot(med_spdf, "data_resp")
  #grid.arrange(predicted_surface,true_sim_outcome, nrow = 1, ncol = 2)
  
  
  #RMSE
  rmse = sqrt(mean((gwr_outcome$SDF$pred - med_spdf$data_resp) ^ 2, na.rm = T))
  
  med.bw <- gwr.sel(data_M ~ data_A + data_C, data = med_spdf)
  #med.bw <- 5
  gwr_mediator <- gwr(data_M ~ data_A + data_C, data = med_spdf, bandwidth = med.bw)
  rmse_M = sqrt(mean((gwr_mediator$SDF$pred - med_spdf$data_M) ^ 2, na.rm = T))
  #0.262
  #rmse_Med[2] = rmse_M
  
  #tot.bw <- gwr.sel(data_resp ~ data_A , data = med_spdf)
  #tot.bw <- 5
  #gwr_tot <- gwr(data_resp ~ data_A , data = med_spdf, bandwidth = tot.bw)
  
  
  #spplot(gwr_outcome$SDF, "data_A") #NDE is data_A from this model
  gwr_outcome$SDF$NIE = gwr_outcome$SDF$data_M * gwr_mediator$SDF$data_A
  #spplot(gwr_outcome$SDF, "NIE")
  
  library("raster")
  
  library("akima")
  r <- rasterFromXYZ(as.data.frame(gwr_outcome$SDF)[, c("x", "y", "pred", "NIE", "data_A")])
  end_gwr <- Sys.time() - start #5.78
  
  return(list("time_svc" = end_svc, "time_gwr" = end_gwr, "gwr_results" = gwr_outcome, "svc_results" = data_svc, "rmse_svc" = c(rmse_svc, rmse_m_svc), "rmse_gwr" = c(rmse, rmse_M)))
}
data_list = list(S1, S2, S3,S4,S5,S6,S7,S8,S9,S10,S11,S12)
j = 1
results = list()
for(i in data_list){
  results[[j]] = SimulationAnalysis(i)
  j = j + 1
}

#Plots

#Pros and cons, Bayesian Analysis Take Manual Priors estimation

mediation_spatialpoints <- rasterToPoints(raster_stack, spatial = T)
med_spdf <- as(mediation_spatialpoints, "SpatialPointsDataFrame")
med_spdf
v1 = variogram(data_resp ~ data_A + data_M, data = med_spdf, cutoff = 8)
fit.variogram(v1)
plot(v1)

#Distance 2, nugget = 0.8, Sill 1.2
#Tau sq is nugget
# Sigma is sill 
#phi is range
r = 2 
#starting are the values coming from the variogram we fit
starting = list("phi"=rep(2, r), "sigma.sq"=rep(1.2, r), "tau.sq"=0.8)
#Tuning are relatively small (I think it is how much it moves)
#The value portion of each tag defines the variance of the Metropolis sampler Normal proposal distribution.
tuning = list("phi"=rep(0.1, r), "sigma.sq"=rep(0.1, r), "tau.sq"=0.1)
priors <- list("phi.Unif" = list(rep(0.1, r), rep(10, r)),
               "sigma.sq.IG" = list(rep(2,r), rep(1, r)),
               "tau.sq.IG" = c(3, 0.5))
corr = "spherical"
#Even with 400 data points and 10000 MCMC samples takes 10 minutes (problem)
bayesian_SVC <- function(corr = "spherical", parameters, starting, tuning, coords, data){
  library("rgeos")
  library("spBayes")
  n.samples = 10000
  #Mediator Model
  coords = coordinates(S1$data_resp)
  m.1 <- spSVC(data_resp ~ data_A + data_M + data_C, coords=coords, starting=starting, svc.cols=c(1,2),
               tuning=tuning, priors=priors, cov.model=corr,
               n.samples=n.samples, n.omp.threads=8, data = med_spdf)
  
  #recover posterior samples
  m.1 <- spRecover(m.1, start=floor(0.75*n.samples), thin=2, n.omp.threads=4)
  
  ##check fitted values
  quant <- function(x){quantile(x, prob=c(0.025, 0.5, 0.975))}
  
  beta.samples <- m.1$p.beta.recover.samples
  w.samples = m.1$p.w.recover.samples
  
  #Random effects so what are the total coefficients
 data_A_re = apply(w.samples, 1, mean)[1:400]
 data_M_re =  apply(w.samples, 1, mean)[401:800]
  y.hat <- apply(m.1$p.y.samples, 1, quant)
  
  #Y model 
  m.2 <- spSVC(data_M ~ data_A + data_C, coords=coords, starting=starting, svc.cols=c(1,2),
               tuning=tuning, priors=priors, cov.model=corr,
               n.samples=n.samples, n.omp.threads=4, data = med_spdf)
  
  #recover posterior samples
  m.2 <- spRecover(m.2, start=floor(0.75*n.samples), thin=2, n.omp.threads=4)
  
  m.hat <- apply(m.2$p.y.samples, 1, quant)
  
  
  
}

#gFormula quick test

Avar <- df_data[,c(4,6)]
Mvar <- df_data$data_M
m_mod = besf_vc(df_data$data_M, Avar, coords = coords,x_sel = FALSE)


#NEW A data
A_0 = raster(matrix(data = rep(0, 400), nrow = 20, ncol = 20))
extent(A_0) = c(0,20,0,20)
C = df_data$data_C
newdata = data.frame(cbind(rep(0,400), C))
test = lm(data_M ~ data_A + data_C, data = df_data)
names(newdata) = c("data_A", "data_C")
newdata$data_M = predict(test, newdata)

test2 = lm(data_resp ~ data_A + data_M + data_C, data = df_data)
data_resp = predict(test2, newdata)
newdata$data_resp = data_resp

newdata_1 = data.frame(cbind(rep(1,400), C))
names(newdata_1) = c("data_A", "data_C")
newdata_1$data_M = predict(test, newdata_1)
newdata_1$data_resp = predict(test2, newdata_1)

newdata$x = newdata_1$x = df_data$x
newdata$y = newdata_1$y = df_data$y
df_data$A_0 = NULL
df_data = rbind(df_data, newdata, newdata_1)

#Error with GFormula 
#Error in optim(fn = lik_resf_vc, par00, par0 = par0, ev = ev, M = MM,  : 
      #           L-BFGS-B needs finite values of 'fn'


time_svc = c(results[[1]]$time_svc, results[[2]]$time_svc, results[[3]]$time_svc, results[[4]]$time_svc, results[[5]]$time_svc, results[[6]]$time_svc, results[[7]]$time_svc, results[[8]]$time_svc, results[[9]]$time_svc, results[[10]]$time_svc, results[[11]]$time_svc, results[[12]]$time_svc)

time_gwr = c(results[[1]]$time_gwr, results[[2]]$time_gwr, results[[3]]$time_gwr, results[[4]]$time_gwr, results[[5]]$time_gwr, results[[6]]$time_gwr, results[[7]]$time_gwr, results[[8]]$time_gwr, results[[9]]$time_gwr, results[[10]]$time_gwr, results[[11]]$time_gwr, results[[12]]$time_gwr)

rmse_svc = c(results[[1]]$rmse_svc, results[[2]]$rmse_svc, results[[3]]$rmse_svc, results[[4]]$rmse_svc, results[[5]]$rmse_svc, results[[6]]$rmse_svc, results[[7]]$rmse_svc, results[[8]]$rmse_svc, results[[9]]$rmse_svc, results[[10]]$rmse_svc, results[[11]]$rmse_svc, results[[12]]$rmse_svc)
ys = abs(c(summary(S1$data_resp)[3],summary(S1$data_M)[3], summary(S2$data_resp)[3],summary(S2$data_M)[3], summary(S3$data_resp)[3],summary(S3$data_M)[3], summary(S4$data_resp)[3],summary(S4$data_M)[3],summary(S5$data_resp)[3],summary(S5$data_M)[3], summary(S6$data_resp)[3],summary(S6$data_M)[3], 
           summary(S7$data_resp)[3],summary(S7$data_M)[3], summary(S8$data_resp)[3],summary(S8$data_M)[3], summary(S9$data_resp)[3],summary(S9$data_M)[3], summary(S10$data_resp)[3],summary(S10$data_M)[3], summary(S11$data_resp)[3],summary(S11$data_M)[3],summary(S12$data_resp)[3],summary(S12$data_M)[3]))
rmse_gwr = c(results[[1]]$rmse_gwr, results[[2]]$rmse_gwr, results[[3]]$rmse_gwr, results[[4]]$rmse_gwr, results[[5]]$rmse_gwr, results[[6]]$rmse_gwr, results[[7]]$rmse_gwr, results[[8]]$rmse_gwr, results[[9]]$rmse_gwr, results[[10]]$rmse_gwr, results[[11]]$rmse_gwr, results[[12]]$rmse_gwr)

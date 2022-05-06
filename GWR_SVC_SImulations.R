
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Mediation Analysis")



#Low total Effect high proportion mediated Scale 0.5#######
FullSim <- function(scale,a,b,c){
  library("spmoran")
  library("spgwr")
  source("CalSim.R")
  
  raster_stack = CalSim(scale,a,b,c)
  start <- Sys.time()
  df_data <-data.frame(rasterToPoints(raster_stack))
  Xvar <- df_data[,c(4,5)]
  coords <- df_data[,c(1,2)]
  y_mod <- besf_vc(df_data$data_resp, df_data[,c("data_A", "data_M", "data_C")], coords = coords, x_sel = FALSE)
  rmse_svc = sqrt(mean((y_mod$pred - df_data$data_resp) ^ 2, na.rm = T))
  
  #y_mod <- besf_vc(df_data$data_resp, df_data$data_A, coords = coords)
  
  Avar <- df_data[,c(4,6)]
  Mvar <- df_data$data_M
  m_mod = besf_vc(df_data$data_M, Avar, coords = coords)
  
  rmse_m_svc = sqrt(mean((m_mod$pred - df_data$data_M) ^ 2, na.rm = T))
  
  df_data$NDE = y_mod$b_vc[,2]
  
  df_data$NIE = y_mod$b_vc[,3] * m_mod$b_vc[,2]
  
  df_data$Y_pred = y_mod$pred
  
  end2 <- Sys.time() - start # 1.79 min
  #r <- rasterFromXYZ(df_data)
  
  start <- Sys.time()
  mediation_spatialpoints <- rasterToPoints(raster_stack, spatial = T)
  med_spdf <- as(mediation_spatialpoints, "SpatialPointsDataFrame")
  col.bw <- gwr.sel(data_resp ~ data_A + data_M + data_C, data = med_spdf)
  
  gwr_outcome <- gwr(data_resp ~ data_A + data_M + data_C, data = med_spdf, bandwidth = col.bw)
  #predicted_surface = spplot(gwr_outcome$SDF, "pred")
  #true_sim_outcome = spplot(med_spdf, "data_resp")
  #grid.arrange(predicted_surface,true_sim_outcome, nrow = 1, ncol = 2)
  
  #RMSE
  rmse = sqrt(mean((gwr_outcome$SDF$pred - med_spdf$data_resp) ^ 2, na.rm = T))
  #rmse/mean(med_spdf$data_resp)
  #0.210
  hist(rnorm(100))
  med.bw <- gwr.sel(data_M ~ data_A + data_C, data = med_spdf)
  #med.bw <- 5
  gwr_mediator <- gwr(data_M ~ data_A + data_C, data = med_spdf, bandwidth = 0.1297171)
  rmse_M = sqrt(mean((gwr_mediator$SDF$pred - med_spdf$data_M) ^ 2, na.rm = T))
  #rmse_M/abs(mean(med_spdf$data_M))
  #0.262
  #rmse_Med = rmse_M
  
  #tot.bw <- gwr.sel(data_resp ~ data_A , data = med_spdf)
  #tot.bw <- 5
  #gwr_tot <- gwr(data_resp ~ data_A , data = med_spdf, bandwidth = tot.bw)
  
  
  #spplot(gwr_outcome$SDF, "data_A") #NDE is data_A from this model
  gwr_outcome$SDF$NIE = gwr_outcome$SDF$data_M * gwr_mediator$SDF$data_A
  gwr_outcome$SDF$Y = med_spdf$data_resp
  #spplot(gwr_outcome$SDF, "NIE")
  
  library("raster")
  
  library("akima")
  r <- rasterFromXYZ(as.data.frame(gwr_outcome$SDF)[, c("x", "y", "pred", "NIE", "data_A")])
  end_gwr <- Sys.time() - start #5.55
  return(list("data_svc" = df_data, "data_gwr" = gwr_outcome$SDF, "rmse_gwr" = rmse, "rmse_m_gwr" = rmse_M, "rmse_svc" = rmse_svc, "rmse_m_svc" = rmse_m_svc, "time_svc" = end2, "time_gwr" = end_gwr))
}


SS1 = FullSim(0.5, 0.1012, .6, .252)
SS2 = FullSim(2, 0.1012, .6, .252)
SS3 = FullSim(5, 0.1012, .6, .252) #This setting is very difficult
SS4 = FullSim(0.5, 0.3, 0.4, 0.75)
SS5 = FullSim(2, 0.3, 0.4, 0.75)
SS6 = FullSim(5, 0.3, 0.4, 0.75)
SS7 = FullSim(0.5, .32, 6.0, 0.8)
SS8 = FullSim(2, .32, 6.0, 0.8)
SS9 = FullSim(5, .32, 6.0, 0.8)
SS10 = FullSim(0.5, 0.94, 4.0, 2.35)
SS11 = FullSim(2, 0.94, 4.0, 2.35)
SS12 = FullSim(5, 0.94, 4.0, 2.35)


sd(c(S1$time_gwr, S2$time_gwr, S3$time_gwr, S4$time_gwr, S5$time_gwr, S6$time_gwr, S7$time_gwr, S8$time_gwr, S9$time_gwr, S10$time_gwr, S11$time_gwr, S12$time_gwr))
#5.79 minutes
#sd 1.00 minutes
sd(c(S1$time_svc, S2$time_svc, S3$time_svc, S4$time_svc, S5$time_svc, S6$time_svc, S7$time_svc, S8$time_svc, S9$time_svc, S10$time_svc, S11$time_svc, S12$time_svc))
#2.44 minutes
#0.9 minutes


# S1_GWR = rasterFromXYZ(as.data.frame(S1$data_gwr)[, c("x", "y", "pred", "NIE", "data_A", "Y")])
# S1_SVC = rasterFromXYZ(S1$data_svc)
# 
# 
# 
# S2_GWR = rasterFromXYZ(as.data.frame(S2$data_gwr)[, c("x", "y", "pred", "NIE", "data_A", "Y")])
# S2_SVC = rasterFromXYZ(S2$data_svc)
# 
# S3_GWR = rasterFromXYZ(as.data.frame(S3$data_gwr)[, c("x", "y", "pred", "NIE", "data_A", "Y")])
# S3_SVC = rasterFromXYZ(S3$data_svc)
# 
# S4_GWR = rasterFromXYZ(as.data.frame(S4$data_gwr)[, c("x", "y", "pred", "NIE", "data_A", "Y")])
# S4_SVC = rasterFromXYZ(S4$data_svc)
# 
# S5_GWR = rasterFromXYZ(as.data.frame(S5$data_gwr)[, c("x", "y", "pred", "NIE", "data_A", "Y")])
# S5_SVC = rasterFromXYZ(S5$data_svc)
# 
# 
# S6_GWR = rasterFromXYZ(as.data.frame(S6$data_gwr)[, c("x", "y", "pred", "NIE", "data_A", "Y")])
# S6_SVC = rasterFromXYZ(S6$data_svc)
# 
# 
# S7_GWR = rasterFromXYZ(as.data.frame(S7$data_gwr)[, c("x", "y", "pred", "NIE", "data_A", "Y")])
# S7_SVC = rasterFromXYZ(S7$data_svc)
# 
# 
# S8_GWR = rasterFromXYZ(as.data.frame(S8$data_gwr)[, c("x", "y", "pred", "NIE", "data_A", "Y")])
# S8_SVC = rasterFromXYZ(S8$data_svc)
# 
# S9_GWR = rasterFromXYZ(as.data.frame(S9$data_gwr)[, c("x", "y", "pred", "NIE", "data_A", "Y")])
# S9_SVC = rasterFromXYZ(S9$data_svc)
# 
S10_GWR = rasterFromXYZ(as.data.frame(SS10$data_gwr)[, c("x", "y", "pred", "NIE", "data_A", "Y")])
S10_SVC = rasterFromXYZ(SS10$data_svc)
# 
# S11_GWR = rasterFromXYZ(as.data.frame(S11$data_gwr)[, c("x", "y", "pred", "NIE", "data_A", "Y")])
# S11_SVC = rasterFromXYZ(S11$data_svc)
# 
# S12_GWR = rasterFromXYZ(as.data.frame(S12$data_gwr)[, c("x", "y", "pred", "NIE", "data_A", "Y")])
# S12_SVC = rasterFromXYZ(S12$data_svc)


plotting <- function(S, Sname){
  S_GWR = rasterFromXYZ(as.data.frame(S$data_gwr)[, c("x", "y", "pred", "NIE", "data_A", "Y")])
  # plot 1 - GWE NDE
  main = paste(Sname, "_GWR_NDE.png", sep = "")
  png(main, width = 6.5, height = 6, units = 'in', res = 300)
  plot(S_GWR$data_A, main = paste0("GWR NDE Estimate"))
  text(-118, 40, "E(NDE) = 4.00")
  dev.off()
  
  # plot 2 - GWR NIE
  main = paste(Sname, "_GWR_NIE.png", sep = "")
  png(main, width = 6.5, height = 6, units = 'in', res = 300)
  plot(S_GWR$NIE, main = paste0("GWR NIE Estimate"))
  text(-118, 40, "E(NIE) = 2.21")
  dev.off()
  
  # plot 3 - GWR Truth - PRed
  main = paste(Sname, "_GWR_TruthPred.png", sep = "")
  png(main, width = 6.5, height = 6, units = 'in', res = 300)
  plot(S_GWR$Y - S_GWR$pred, main = paste0("GWR Truth - Pred"), zlim = c(-2,2))
  dev.off()
  S_SVC = rasterFromXYZ(S$data_svc)
  # plot 4 - SVC Truth - PRed
  main = paste(Sname, "_SVC_TruthPred.png", sep = "")
  png(main, width = 6.5, height = 6, units = 'in', res = 300)
  plot(S_SVC$data_resp - S_SVC$Y_pred, main = paste0("SVC Truth - Pred"), zlim = c(-2,2))
  dev.off()
  
  # plot 5 - SVC NDE
  main = paste(Sname, "-SVC_NDE.png", sep = "")
  png(main, width = 6.5, height = 6, units = 'in', res = 300)
  plot(S_SVC$NDE, main =paste0("SVC NDE Estimate"))
  text(-118, 40, "E(NDE) = 4.00")
  dev.off()
  
  # plot 6 - SVC NIE
  main = paste(Sname, "_SVC_NIE.png", sep = "")
  png(main, width = 6.5, height = 6, units = 'in', res = 300)
  plot(S_SVC$NIE, main = paste0("SVC NIE Estimate"))
  text(-118, 40, "E(NIE) = 2.21")
  dev.off()
}


plotting(S1, "S1")
plotting(S2, "S2")
plotting(S3, "S3")
plotting(S4, "S4")
plotting(S5, "S5")
plotting(S6, "S6")
plotting(S7, "S7")
plotting(S8, "S8")
plotting(S9, "S9")
plotting(SS10, "S10")
plotting(S11, "S11")
plotting(S12, "S12")



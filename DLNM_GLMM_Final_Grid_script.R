################################################################################
# Title: Modelling the Effect of Climatic Drivers and Cross-Border Travel on 
#        Sporadic Outbreaks of Dengue in Mauritius

# Aim: To quantify the association between climatic variables 
#      (temperature, rainfall and vapor pressure deficit) 
#      and international resident arrivals and dengue incidence in Mauritius


# Author: Mohabeer Teeluck, Oyelola Adegboye*, Benn Sartorius, Eloise Skinner, 
#         Stephan Karl, Alec Henderson, Dr. M.F. Khodaboccus, Diana P. Iyaloo, 
#         Emma McBryde*
# *Lead supervisors: Oyelola Adegboye and Emma McBryde
# Institution: University of Queensland and James Cook University
# Date: 27/04/2026
# R version: 4.4.1
################################################################################

# DLNM grid search - Temperature, rainfall, VPD and resident arrivals
# Temperature: lag 4,8,12  | df_var 2,3,4 | df_lag 2,3
# Rainfall:    lag 8,10,12 | df_var 2,3,4 | df_lag 2,3
# VPD:         lag 8,12    | df_var 2,3,4 | df_lag 2,3
# Arrivals:    lag 4,8,12  | df_var 2,3   | df_lag 2,3
################################################################################
# Replace tmin throughout with tmean or tmax to run other temperature models.

rm(list = ls())
gc()

library(glmmTMB)
library(dlnm)
library(splines)
library(dplyr)
library(tidyr)

set.seed(202)
N_WORKERS <- 10

# 1) Data preprocessing

data <- read.csv(file.choose(),check.names=TRUE)

data <- data %>%
  filter(!is.na(population) & population>0,!is.na(local_cases),!is.na(tmin),
         !is.na(rain),!is.na(VPD),!is.na(imported_cases),!is.na(res_arr),
         !is.na(state_name),year>=2014,year<=2024) %>%
  mutate(week=as.integer(week),year=as.integer(year),
         local_cases=as.integer(local_cases),imported_cases=as.numeric(imported_cases),
         res_arr_scaled=as.numeric(res_arr)/1000,tmin=as.numeric(tmin),
         rain=as.numeric(rain),VPD=as.numeric(VPD)) %>%
  filter(year<2020 | (year==2020 & week<14) | (year==2021 & week>40) | year>2021) %>%
  arrange(state_name,year,week) %>%
  group_by(state_name) %>%
  mutate(time_index=row_number()) %>%
  ungroup() %>%
  mutate(analysis_segment=if_else(year<=2020,"pre_gap","post_gap")) %>%
  group_by(state_name,analysis_segment) %>%
  mutate(lag_imported=lag(imported_cases,1),
         lagged_local_cases=lag(local_cases,1),
         outbreak_level=case_when(
           year %in% c(2020,2023,2024) ~ "high",
           year %in% c(2014,2015,2019) ~ "moderate",
           TRUE ~ "none")) %>%
  ungroup() %>%
  mutate(dlnm_group=interaction(state_name,analysis_segment,drop=TRUE),
         state_name=factor(state_name),year=factor(year),
         outbreak_level=factor(outbreak_level,levels=c("none","moderate","high")))


# 2) Common analysis sample
# Maximum lag anywhere in grid = 12 weeks.
# Common sample defined WITHIN each DLNM segment.

MAX_LAG <- 12

data <- data %>%
  group_by(dlnm_group) %>%
  mutate(grid_position=row_number()) %>%
  ungroup()

common_rows <- data$grid_position>MAX_LAG
data_grid <- droplevels(data[common_rows,,drop=FALSE])

stopifnot(all(complete.cases(
  data_grid[,c("local_cases","population","tmin","rain","VPD","res_arr_scaled",
               "lagged_local_cases","lag_imported","time_index","outbreak_level",
               "state_name","year")]
)))

cat("\nCommon-sample N =",nrow(data_grid),"\n")


# 3) Candidate specifications

tmin_specs <- expand.grid(lag=c(4,8,12),df_var=c(2,3,4),df_lag=c(2,3))
rain_specs <- expand.grid(lag=c(8,10,12),df_var=c(2,3,4),df_lag=c(2,3))
vpd_specs <- expand.grid(lag=c(8,12),df_var=c(2,3,4),df_lag=c(2,3))
arrival_specs <- expand.grid(lag=c(4,8,12),df_var=c(2,3),df_lag=c(2,3))

tmin_specs$id <- seq_len(nrow(tmin_specs))
rain_specs$id <- seq_len(nrow(rain_specs))
vpd_specs$id <- seq_len(nrow(vpd_specs))
arrival_specs$id <- seq_len(nrow(arrival_specs))


# 4) Cross-basis setup

make_cb <- function(variable,lag,df_var,df_lag){
  crossbasis(data[[variable]],lag=lag,
             argvar=list(fun="ns",df=df_var),
             arglag=list(fun="ns",df=df_lag),
             group=data$dlnm_group)
}

build_basis_list <- function(specs,variable){
  lapply(seq_len(nrow(specs)),function(i)
    make_cb(variable,specs$lag[i],specs$df_var[i],specs$df_lag[i]))
}

cat("\nBuilding cross-bases...\n")

cb_tmin_all <- build_basis_list(tmin_specs,"tmin")
cb_rain_all <- build_basis_list(rain_specs,"rain")
cb_vpd_all <- build_basis_list(vpd_specs,"VPD")
cb_res_all <- build_basis_list(arrival_specs,"res_arr_scaled")

# Safety check: every basis must be complete on the common sample
basis_ok <- function(x) all(complete.cases(x[common_rows,,drop=FALSE]))

stopifnot(
  all(vapply(cb_tmin_all,basis_ok,logical(1))),
  all(vapply(cb_rain_all,basis_ok,logical(1))),
  all(vapply(cb_vpd_all,basis_ok,logical(1))),
  all(vapply(cb_res_all,basis_ok,logical(1)))
)


# 5) Model grid

grid <- expand.grid(
  tmin_id=seq_len(nrow(tmin_specs)),
  rain_id=seq_len(nrow(rain_specs)),
  vpd_id=seq_len(nrow(vpd_specs)),
  res_id=seq_len(nrow(arrival_specs))
)

grid$model_id <- seq_len(nrow(grid))

cat("Total candidate models:",nrow(grid),"\n")

# 6) Model control

ctrl <- glmmTMBControl(
  optimizer=optim,
  optArgs=list(method="BFGS"),
  optCtrl=list(maxit=2000,reltol=1e-8)
)

# 7) Output files and resuming

progress_file   <- "DLNM_tmin_GRID_46656_PROGRESS.csv"
final_file      <- "DLNM_tmin_GRID_46656_ALL.csv"
ranked_file     <- "DLNM_tmin_GRID_46656_AIC_RANKED.csv"
shortlist2_file <- "DLNM_tmin_GRID_DELTA_AIC2.csv"
shortlist4_file <- "DLNM_tmin_GRID_DELTA_AIC4.csv"
shortlist10_file <- "DLNM_tmin_GRID_DELTA_AIC10.csv"

if(file.exists(progress_file)){
  results <- read.csv(progress_file,check.names=FALSE)
  completed_ids <- results$model_id
  cat("\nRESUMING:",length(completed_ids),"models already completed.\n")
} else {
  results <- data.frame()
  completed_ids <- integer(0)
}



# 8) Model fitting

fit_one <- function(i){
  
  g <- grid[i, ]
  
  ts <- tmin_specs[g$tmin_id, ]
  rs <- rain_specs[g$rain_id, ]
  vs <- vpd_specs[g$vpd_id, ]
  as <- arrival_specs[g$res_id, ]
  
  cb_tmin <- cb_tmin_all[[g$tmin_id]][common_rows,,drop=FALSE]
  cb_rain  <- cb_rain_all[[g$rain_id]][common_rows,,drop=FALSE]
  cb_vpd   <- cb_vpd_all[[g$vpd_id]][common_rows,,drop=FALSE]
  cb_res   <- cb_res_all[[g$res_id]][common_rows,,drop=FALSE]
  
  fit <- tryCatch(
    suppressWarnings(
      glmmTMB(
        local_cases ~ cb_tmin + cb_rain + cb_vpd + cb_res +
          ns(time_index,df=11) + outbreak_level +
          lagged_local_cases + lag_imported +
          (1|state_name) + (1|year),
        family=nbinom1(),
        offset=log(population),
        data=data_grid,
        control=ctrl
      )
    ),
    error=function(e) NULL
  )
  
  if(is.null(fit)){
    return(data.frame(
      model_id=g$model_id,
      lag_tmin=ts$lag, dfvar_tmin=ts$df_var, dflag_tmin=ts$df_lag,
      lag_rain=rs$lag, dfvar_rain=rs$df_var, dflag_rain=rs$df_lag,
      lag_vpd=vs$lag, dfvar_vpd=vs$df_var, dflag_vpd=vs$df_lag,
      lag_res=as$lag, dfvar_res=as$df_var, dflag_res=as$df_lag,
      N=NA, convergence=NA, pdHess=FALSE,
      AIC=NA, BIC=NA, Pearson=NA
    ))
  }
  
  conv <- fit$fit$convergence
  pd <- isTRUE(fit$sdr$pdHess)
  valid <- conv==0 && pd
  
  if(valid){
    rp <- residuals(fit,type="pearson")
    aic <- AIC(fit)
    bic <- BIC(fit)
    pearson <- sum(rp^2,na.rm=TRUE)/df.residual(fit)
  } else {
    aic <- NA_real_
    bic <- NA_real_
    pearson <- NA_real_
  }
  
  data.frame(
    model_id=g$model_id,
    lag_tmin=ts$lag, dfvar_tmin=ts$df_var, dflag_tmin=ts$df_lag,
    lag_rain=rs$lag, dfvar_rain=rs$df_var, dflag_rain=rs$df_lag,
    lag_vpd=vs$lag, dfvar_vpd=vs$df_var, dflag_vpd=vs$df_lag,
    lag_res=as$lag, dfvar_res=as$df_var, dflag_res=as$df_lag,
    N=nobs(fit), convergence=conv, pdHess=pd,
    AIC=aic, BIC=bic, Pearson=pearson
  )
}


# 9) Parallel grid search

start_time <- Sys.time()

cl <- parallel::makeCluster(N_WORKERS)

parallel::clusterEvalQ(cl,{
  library(glmmTMB)
  library(dlnm)
  library(splines)
})

parallel::clusterExport(
  cl,
  c("fit_one","grid","tmin_specs","rain_specs",
    "vpd_specs","arrival_specs","cb_tmin_all",
    "cb_rain_all","cb_vpd_all","cb_res_all",
    "common_rows","data_grid","ctrl"),
  envir=.GlobalEnv
)

pending_i <- which(!grid$model_id %in% completed_ids)

cat("\nUsing",N_WORKERS,"parallel workers\n")
cat("Models remaining:",length(pending_i),"of",nrow(grid),"\n\n")

# Run in batches of 50, saving progress after each batch
BATCH_SIZE <- 50
batches <- split(pending_i,ceiling(seq_along(pending_i)/BATCH_SIZE))

for(b in seq_along(batches)){
  
  batch_i <- batches[[b]]
  batch_results <- parallel::parLapply(cl,batch_i,fit_one)
  
  results <- bind_rows(results,bind_rows(batch_results)) %>%
    distinct(model_id,.keep_all=TRUE)
  
  write.csv(results,progress_file,row.names=FALSE)
  
  elapsed <- as.numeric(difftime(Sys.time(),start_time,units="mins"))
  n_completed <- nrow(results)
  pct <- 100*n_completed/nrow(grid)
  
  cat(sprintf(
    "Completed %d / %d models | %.2f%% | %.1f min elapsed\n",
    n_completed,nrow(grid),pct,elapsed
  ))
}

parallel::stopCluster(cl)
cat("\nParallel fitting complete.\n")



# 10) Save results

results <- results %>%
  distinct(model_id,.keep_all=TRUE) %>%
  arrange(model_id)

write.csv(results,final_file,row.names=FALSE)


# 11) Model ranking

results <- results %>%
  mutate(valid=convergence==0 & pdHess &
           N==nrow(data_grid) &
           is.finite(AIC) & is.finite(BIC))

valid_results <- results %>%
  filter(valid) %>%
  arrange(AIC,BIC) %>%
  mutate(
    delta_AIC=AIC-min(AIC),
    delta_BIC=BIC-min(BIC),
    akaike_weight=exp(-0.5*delta_AIC)/sum(exp(-0.5*delta_AIC))
  )


# 12) AIC shortlists

shortlist2 <- valid_results %>% filter(delta_AIC<=2)
shortlist4 <- valid_results %>% filter(delta_AIC<=4)
shortlist10 <- valid_results %>% filter(delta_AIC<=10)

write.csv(valid_results,ranked_file,row.names=FALSE)
write.csv(shortlist2,shortlist2_file,row.names=FALSE)
write.csv(shortlist4,shortlist4_file,row.names=FALSE)
write.csv(shortlist10,shortlist10_file,row.names=FALSE)


# 13) Summary

cat("\nTemperature grid search complete\n")
cat("Total models:",nrow(grid),"\n")
cat("Common-sample N:",nrow(data_grid),"\n")
cat("Valid models:",nrow(valid_results),"\n")
cat("Invalid / failed:",nrow(grid)-nrow(valid_results),"\n")
cat("Delta-AIC <= 2:",nrow(shortlist2),"\n")
cat("Delta-AIC <= 4:",nrow(shortlist4),"\n")
cat("Delta-AIC <= 10:",nrow(shortlist10),"\n\n")

cat("BEST AIC MODEL:\n")
print(head(valid_results,1))

cat("\nTOP 20 BY AIC:\n")
print(head(valid_results,20))

cat("\nFILES SAVED:\n",
    final_file,"\n",ranked_file,"\n",
    shortlist2_file,"\n",shortlist4_file,"\n",
    shortlist10_file,"\n")

cat("\nDONE\n")

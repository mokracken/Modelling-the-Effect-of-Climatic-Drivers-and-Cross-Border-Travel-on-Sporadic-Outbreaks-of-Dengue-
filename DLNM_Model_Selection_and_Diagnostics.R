############################################################
# SCRIPT 2 — DLNM model selection and diagnostics
# Grid CSV first; weekly dengue data second. Temperature is detected automatically.
# Refits competitive models and runs the original inferential and residual checks.
# PURPOSE
# 1. Validate the targeted model grid file.
# 2. Verify the known-good positive-control model.
# 3. Refit ALL potentially competitive grid models.
# 4. Require exact grid/refit reproducibility.
# 5. Check VCOV, fixed-effect rank, RR/CI and full RR curves.
# 6. Run DHARMa + temporal residual diagnostics on the
#    competitive inferentially usable models.
#
# DOES NOT modify or rerun the model grid.
#
# PARTIAL GRID:
#   Results are explicitly PROVISIONAL.
#
# COMPLETE GRID:
#   Results become the final candidate screen.
############################################################

rm(list=ls()); gc()
library(glmmTMB)
library(dlnm)
library(DHARMa)
library(splines)
library(dplyr)
set.seed(202)
grid_results <- read.csv(file.choose(), check.names=FALSE)
options <- c("tmin", "tmean", "tmax")
found <- options[sapply(options, function(x)
  all(paste0(c("lag_", "dfvar_", "dflag_"), x) %in% names(grid_results)) )]
if(length(found) != 1)
  stop("Cannot identify temperature variable in grid-search CSV")
temp_var <- found[1]
cat("Detected temperature:", temp_var, "\n")
temp_settings <- list( tmin=list( label="Tmin", expected_grid=31104, lag=c(8,12), df_var=c(2,3,4),
                                  df_lag=c(2,3) ), tmax=list( label="Tmax", expected_grid=46656, lag=c(4,8,12), df_var=c(2,3,4),
                                                              df_lag=c(2,3) ), tmean=list( label="Tmean", expected_grid=46656, lag=c(4,8,12), df_var=c(2,3,4),
                                                                                           df_lag=c(2,3) ) )
temp_set <- temp_settings[[temp_var]]
if(anyNA(c( temp_set$expected_grid, temp_set$lag, temp_set$df_var, temp_set$df_lag ))){
  stop("Grid settings have not yet been entered for ",temp_var)
}

MAX_LAG          <- 12
EXPECTED_GRID    <- temp_set$expected_grid
AIC_SCREEN_DELTA <- 10
MATCH_TOL        <- 1e-4
DHARMA_SIMS      <- 2500
output_dir <- file.path("outputs",temp_var)
canonical_dir <- file.path(output_dir,"canonical_models")
dir.create( canonical_dir, showWarnings=FALSE, recursive=TRUE )
out <- function(x) file.path(output_dir,x)

# 1) LOAD + PREPROCESS ORIGINAL WEEKLY DATA
data_full <- read.csv(file.choose(),check.names=TRUE)
stopifnot(temp_var %in% names(data_full))
data_full <- data_full %>% mutate( temp=as.numeric(.data[[temp_var]]), week=as.integer(week),
                                   year=as.integer(year), local_cases=as.integer(local_cases), imported_cases=as.numeric(imported_cases),
                                   res_arr_scaled=as.numeric(res_arr)/1000, rain=as.numeric(rain), VPD=as.numeric(VPD) ) %>% filter(
                                     !is.na(population) & population>0, !is.na(local_cases), !is.na(temp), !is.na(rain), !is.na(VPD),
                                     !is.na(imported_cases), !is.na(res_arr), !is.na(state_name), year>=2014, year<=2024 ) %>% filter(
                                       year<2020 | (year==2020 & week<14) | (year==2021 & week>40) | year>2021 ) %>%
  arrange(state_name,year,week) %>% group_by(state_name) %>% mutate(time_index=row_number()) %>% ungroup() %>%
  mutate( analysis_segment=if_else( year<=2020,"pre_gap","post_gap" ) ) %>%
  group_by(state_name,analysis_segment) %>% mutate( lag_imported=lag(imported_cases,1),
                                                    lagged_local_cases=lag(local_cases,1), outbreak_level=case_when( year %in% c(2020,2023,2024) ~ "high",
                                                                                                                     year %in% c(2014,2015,2019) ~ "moderate", TRUE ~ "none" ) ) %>% ungroup() %>% mutate(
                                                                                                                       dlnm_group=interaction( state_name,analysis_segment,drop=TRUE ), state_name=factor(state_name),
                                                                                                                       year=factor(year), outbreak_level=factor( outbreak_level, levels=c("none","moderate","high") ) ) %>%
  group_by(dlnm_group) %>% mutate(grid_position=row_number()) %>% ungroup()
common_rows <- data_full$grid_position>MAX_LAG
data_grid <- droplevels( data_full[common_rows,,drop=FALSE] )
data <- data_grid
stopifnot( nrow(data_full)==4329, nrow(data_grid)==4113, sum(data_grid$local_cases)==7261,
           nlevels(data_full$dlnm_group)==18 )
stopifnot( all(complete.cases( data_grid[,c( "local_cases","population","temp","rain","VPD",
                                             "res_arr_scaled","lagged_local_cases","lag_imported", "time_index","outbreak_level","state_name","year"
)] )) )
cat( "\n============================================\n", "DATA CHECK\n",
     "============================================\n", "Full rows =",nrow(data_full),"\n",
     "Common N =",nrow(data_grid),"\n", "Cases =",sum(data_grid$local_cases),"\n",
     "Districts =",nlevels(data_grid$state_name),"\n", "DLNM segments =",nlevels(data_full$dlnm_group),"\n" )

# 2) RECONSTRUCT TARGET GRID
temp_specs <- expand.grid( lag=temp_set$lag, df_var=temp_set$df_var, df_lag=temp_set$df_lag,
                           KEEP.OUT.ATTRS=FALSE )
rain_specs <- expand.grid( lag=c(8,10,12), df_var=c(2,3,4), df_lag=c(2,3), KEEP.OUT.ATTRS=FALSE )
vpd_specs <- expand.grid( lag=c(8,12), df_var=c(2,3,4), df_lag=c(2,3), KEEP.OUT.ATTRS=FALSE )
arrival_specs <- expand.grid( lag=c(4,8,12), df_var=c(2,3), df_lag=c(2,3), KEEP.OUT.ATTRS=FALSE )
target_grid <- expand.grid( temp_id=seq_len(nrow(temp_specs)), rain_id=seq_len(nrow(rain_specs)),
                            vpd_id=seq_len(nrow(vpd_specs)), res_id=seq_len(nrow(arrival_specs)), KEEP.OUT.ATTRS=FALSE )
target_grid$model_id <- seq_len(nrow(target_grid))
target_grid <- target_grid %>% mutate( lag_temp=temp_specs$lag[temp_id],
                                       dfvar_temp=temp_specs$df_var[temp_id], dflag_temp=temp_specs$df_lag[temp_id],
                                       lag_rain=rain_specs$lag[rain_id], dfvar_rain=rain_specs$df_var[rain_id],
                                       dflag_rain=rain_specs$df_lag[rain_id], lag_vpd=vpd_specs$lag[vpd_id], dfvar_vpd=vpd_specs$df_var[vpd_id],
                                       dflag_vpd=vpd_specs$df_lag[vpd_id], lag_res=arrival_specs$lag[res_id],
                                       dfvar_res=arrival_specs$df_var[res_id], dflag_res=arrival_specs$df_lag[res_id] )
stopifnot(nrow(target_grid)==EXPECTED_GRID)
cat( "\nTarget grid reconstructed correctly:", nrow(target_grid),"models\n",
     "Temperature =",temp_set$label,"\n" )

# 3) VALIDATE GRID RESULTS
temp_cols <- paste0(c("lag_","dfvar_","dflag_"),temp_var)
if(!all(temp_cols %in% names(grid_results)))
  stop("Selected grid file does not match temp_var = ",temp_var)
names(grid_results)[match(temp_cols,names(grid_results))] <- c("lag_temp","dfvar_temp","dflag_temp")
required_cols <- c( "model_id", "lag_temp","dfvar_temp","dflag_temp", "lag_rain","dfvar_rain","dflag_rain",
                    "lag_vpd","dfvar_vpd","dflag_vpd", "lag_res","dfvar_res","dflag_res", "N","convergence","pdHess","AIC","BIC"
)
missing_cols <- setdiff( required_cols, names(grid_results) )
if(length(missing_cols)){
  stop( "WRONG GRID FILE: missing columns: ", paste(missing_cols,collapse=", ") )
}
if(anyDuplicated(grid_results$model_id)){
  stop("GRID FILE CONTAINS DUPLICATE model_id VALUES.")
}
if(any( !is.finite(grid_results$model_id) | grid_results$model_id<1 | grid_results$model_id>EXPECTED_GRID |
        grid_results$model_id!=as.integer(grid_results$model_id) )){
  stop("GRID FILE CONTAINS INVALID model_id VALUES.")
}
if(any( !is.na(grid_results$N) & grid_results$N!=4113 )){
  stop( "WRONG GRID FILE / WRONG SAMPLE: ", "non-missing N values are not all 4113." )
}

expected_rows <- target_grid[ match( grid_results$model_id, target_grid$model_id ), , drop=FALSE ]
spec_cols <- c( "lag_temp","dfvar_temp","dflag_temp", "lag_rain","dfvar_rain","dflag_rain",
                "lag_vpd","dfvar_vpd","dflag_vpd", "lag_res","dfvar_res","dflag_res" )
bad_mapping <- rep(FALSE,nrow(grid_results))
for(v in spec_cols){
  bad_mapping <- bad_mapping | is.na(grid_results[[v]]) | grid_results[[v]]!=expected_rows[[v]]
}
if(any(bad_mapping)){
  bad_ids <- grid_results$model_id[ which(bad_mapping)[1:min(10,sum(bad_mapping))] ]
  stop( "GRID FILE DOES NOT MATCH THE TARGET GRID. ", "First mismatching model IDs: ",
        paste(bad_ids,collapse=", ") )
}

n_completed <- nrow(grid_results)
complete_grid <- n_completed==EXPECTED_GRID && setequal( grid_results$model_id, seq_len(EXPECTED_GRID) )
TAG <- if(complete_grid) "FINAL" else "PROVISIONAL"
cat( "\n============================================\n", "GRID AUDIT\n",
     "============================================\n", "Rows completed =",n_completed,"/",EXPECTED_GRID,"\n",
     "Grid status =",TAG,"\n", "Exact model-ID/specification mapping = PASS\n",
     "Non-missing sample sizes = 4113\n" )
if(!complete_grid){
  cat( "\nIMPORTANT: GRID IS PARTIAL.\n", "All rankings produced below are PROVISIONAL.\n" )
}

# 4) HELPERS
ctrl <- glmmTMBControl( optimizer=optim, optArgs=list(method="BFGS"), optCtrl=list( maxit=2000, reltol=1e-8 )
)
make_cb <- function(x,lag,dfv,dfl){
  crossbasis( x, lag=lag, argvar=list(fun="ns",df=dfv), arglag=list(fun="ns",df=dfl),
              group=data_full$dlnm_group )
}

get_sd <- function(fit,g){
  tryCatch( as.numeric( attr( VarCorr(fit)$cond[[g]], "stddev" )[1] ), error=function(e) NA_real_ )
}

get_grad <- function(fit){
  x <- tryCatch( fit$sdr$gradient.fixed, error=function(e) NULL )
  if(is.null(x) || !length(x))
    return(NA_real_)
  x <- x[is.finite(x)]
  if(length(x))
    max(abs(x))
  else
    NA_real_
}

safe_cor <- function(x,y){
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if( length(x)<3 || sd(x)==0 || sd(y)==0 ){
    return(NA_real_)
  }
  cor(x,y)
}

safe_max <- function(x){
  x <- x[is.finite(x)]
  if(length(x))
    max(x)
  else
    NA_real_
}

safe_min <- function(x){
  x <- x[is.finite(x)]
  if(length(x))
    min(x)
  else
    NA_real_
}

safe_p <- function(expr){
  tryCatch( expr$p.value, error=function(e) NA_real_ )
}

# RR at requested percentiles
rr_one <- function( cb,prefix,vec,probs,b,V,label ){
  q <- as.numeric( quantile( vec,probs, na.rm=TRUE, names=FALSE ) )
  cen <- median( vec, na.rm=TRUE )
  ind <- grep( paste0("^",prefix), names(b) )
  labs <- sprintf( "P%02d", round(probs*100) )
  out <- setNames( rep( NA_real_, 1+4*length(probs) ), c( paste0(label,"_Reference"), unlist( lapply( labs,
                                                                                                      function(z)c( paste0(label,"_",z), paste0(label,"_RR_",z), paste0(label,"_LCL_",z),
                                                                                                                    paste0(label,"_UCL_",z) ) ) ) ) )
  out[1] <- cen
  for(k in seq_along(q)){
    out[ paste0(label,"_",labs[k]) ] <- q[k]
  }
  if( length(ind)!=ncol(cb) || any(!is.finite(b[ind])) || any(!is.finite( V[ind,ind,drop=FALSE] )) ){
    return(out)
  }
  
  cp <- tryCatch( crosspred( cb, coef=b[ind], vcov=V[ind,ind,drop=FALSE], model.link="log", cen=cen, at=q ),
                  error=function(e) NULL )
  if(is.null(cp))
    return(out)
  for(k in seq_along(q)){
    j <- which.min( abs( as.numeric(cp$predvar)-q[k] ) )
    out[ paste0(label,"_RR_",labs[k]) ] <- cp$allRRfit[j]
    out[ paste0(label,"_LCL_",labs[k]) ] <- cp$allRRlow[j]
    out[ paste0(label,"_UCL_",labs[k]) ] <- cp$allRRhigh[j]
  }
  out
}

# Whole cumulative curve stability

# NO arbitrary RR cutoff.

# Only checks whether curve and CI are finite/positive.
curve_check <- function( cb,prefix,vec,b,V,label, probs=c(.05,.95), n=201 ){
  q <- as.numeric( quantile( vec,probs, na.rm=TRUE, names=FALSE ) )
  cen <- median( vec, na.rm=TRUE )
  at <- seq( q[1], q[2], length.out=n )
  ind <- grep( paste0("^",prefix), names(b) )
  empty <- c( setNames(q[1],paste0(label,"_Curve_Low")), setNames(q[2],paste0(label,"_Curve_High")),
              setNames(NA_real_,paste0(label,"_Curve_Min_RR")), setNames(NA_real_,paste0(label,"_Curve_Max_RR")),
              setNames(NA_real_,paste0(label,"_Curve_Min_LCL")), setNames(NA_real_,paste0(label,"_Curve_Max_UCL")),
              setNames(0,paste0(label,"_Curve_Finite")) )
  if( length(ind)!=ncol(cb) || any(!is.finite(b[ind])) || any(!is.finite( V[ind,ind,drop=FALSE] )) ){
    return(empty)
  }
  
  cp <- tryCatch( crosspred( cb, coef=b[ind], vcov=V[ind,ind,drop=FALSE], model.link="log", cen=cen, at=at ),
                  error=function(e) NULL )
  if(is.null(cp))
    return(empty)
  rr <- as.numeric(cp$allRRfit)
  lo <- as.numeric(cp$allRRlow)
  hi <- as.numeric(cp$allRRhigh)
  finite_ok <- all(is.finite(c(rr,lo,hi))) && all(c(rr,lo,hi)>0)
  c( setNames(q[1],paste0(label,"_Curve_Low")), setNames(q[2],paste0(label,"_Curve_High")),
     setNames(safe_min(rr),paste0(label,"_Curve_Min_RR")),
     setNames(safe_max(rr),paste0(label,"_Curve_Max_RR")),
     setNames(safe_min(lo),paste0(label,"_Curve_Min_LCL")),
     setNames(safe_max(hi),paste0(label,"_Curve_Max_UCL")), setNames( as.integer(finite_ok),
                                                                      paste0(label,"_Curve_Finite") ) )
}

# 5) FIND POSITIVE CONTROL AUTOMATICALLY
control_candidates <- grid_results %>% filter( N==4113, convergence==0, pdHess==TRUE, is.finite(AIC),
                                               is.finite(BIC) ) %>% arrange(AIC,BIC) %>% mutate(grid_rank=row_number())
if(!nrow(control_candidates)){
  stop("NO CORE-VALID GRID MODEL AVAILABLE FOR CONTROL SEARCH.")
}
cat( "\n============================================\n", "AUTOMATIC POSITIVE-CONTROL SEARCH\n",
     "============================================\n", "Core-valid models available =",
     nrow(control_candidates),"\n", "Models will be tested in ascending AIC order.\n" )
control_fit <- NULL
control_grid_raw <- NULL
control_attempts <- data.frame()
for(i in seq_len(nrow(control_candidates))){
  g <- control_candidates[i,]
  cat( "\n----------------------------------------------------\n", "Control attempt",i,
       "| grid rank",g$grid_rank, "| model",g$model_id, "| grid AIC",g$AIC,"\n", temp_set$label,
       g$lag_temp,g$dfvar_temp,g$dflag_temp, "| Rain", g$lag_rain,g$dfvar_rain,g$dflag_rain, "| VPD",
       g$lag_vpd,g$dfvar_vpd,g$dflag_vpd, "| Res", g$lag_res,g$dfvar_res,g$dflag_res, "\n" )
  
  # Full cross-bases
  cb_temp_full <- make_cb( data_full$temp, g$lag_temp, g$dfvar_temp, g$dflag_temp )
  cb_rain_full <- make_cb( data_full$rain, g$lag_rain, g$dfvar_rain, g$dflag_rain )
  cb_vpd_full <- make_cb( data_full$VPD, g$lag_vpd, g$dfvar_vpd, g$dflag_vpd )
  cb_res_arr_full <- make_cb( data_full$res_arr_scaled, g$lag_res, g$dfvar_res, g$dflag_res )
  
  # Common sample
  cb_temp <- cb_temp_full[ common_rows,,drop=FALSE ]
  cb_rain <- cb_rain_full[ common_rows,,drop=FALSE ]
  cb_vpd <- cb_vpd_full[ common_rows,,drop=FALSE ]
  cb_res_arr <- cb_res_arr_full[ common_rows,,drop=FALSE ]
  fit_warnings <- character()
  fit <- withCallingHandlers( tryCatch( glmmTMB( local_cases ~ cb_temp + cb_rain + cb_vpd + cb_res_arr +
                                                   ns(time_index,df=11) + outbreak_level + lagged_local_cases + lag_imported + (1|state_name) +
                                                   (1|year), family=nbinom1(), offset=log(population), data=data_grid, control=ctrl ),
                                        error=function(e)e ), warning=function(w){
                                          fit_warnings <<- c( fit_warnings, conditionMessage(w) )
                                          invokeRestart("muffleWarning")
                                        } )
  if(inherits(fit,"error")){
    control_attempts <- bind_rows( control_attempts, data.frame( grid_rank=g$grid_rank, model_id=g$model_id,
                                                                 Grid_AIC=g$AIC, Refit_AIC=NA_real_, status="FIT_ERROR", Warning_count=length(fit_warnings),
                                                                 detail=conditionMessage(fit) ) )
    write.csv( control_attempts, out("CONTROL_DISCOVERY_ATTEMPTS.csv"), row.names=FALSE )
    next
  }
  
  b <- fixef(fit)$cond
  V <- suppressWarnings( tryCatch( vcov(fit)$cond, error=function(e) NULL ) )
  X <- tryCatch( model.matrix(fit), error=function(e) NULL )
  n_fit <- nobs(fit)
  conv <- fit$fit$convergence
  pd <- isTRUE(fit$sdr$pdHess)
  aic <- suppressWarnings(AIC(fit))
  bic <- suppressWarnings(BIC(fit))
  aic_diff <- aic-g$AIC
  bic_diff <- bic-g$BIC
  n_ok <- n_fit==4113
  core_ok <- conv==0 && pd && is.finite(aic) && is.finite(bic)
  pipeline_match <- is.finite(aic_diff) && is.finite(bic_diff) && abs(aic_diff)<=MATCH_TOL &&
    abs(bic_diff)<=MATCH_TOL
  coef_ok <- length(b)>0 && all(is.finite(b))
  v_ok <- !is.null(V) && length(V)>0 && all(is.finite(V))
  rank_ok <- !is.null(X) && qr(X)$rank==ncol(X)
  status <- if(!n_ok){
    "N_MISMATCH"
  } else if(!core_ok){
    "CORE_FIT_FAIL"
  } else if(!pipeline_match){
    "PIPELINE_MISMATCH"
  } else if(!coef_ok){
    "NONFINITE_COEFFICIENT"
  } else if(!v_ok){
    "NONFINITE_VCOV"
  } else if(!rank_ok){
    "RANK_DEFICIENT"
  } else {
    "OK"
  }
  
  rr_all <- numeric(0)
  curves_all <- numeric(0)
  
  # RR / CI / whole-curve checks
  if(status=="OK"){
    rr_all <- c( rr_one( cb_temp_full,"cb_temp", data$temp,c(.05,.95), b,V,temp_set$label ), rr_one(
      cb_rain_full,"cb_rain", data$rain,c(.05,.95), b,V,"Rain" ), rr_one( cb_vpd_full,"cb_vpd",
                                                                          data$VPD,c(.05,.95), b,V,"VPD" ), rr_one( cb_res_arr_full,"cb_res_arr", data$res_arr_scaled,
                                                                                                                    c(.05,.10,.90,.95), b,V,"Res" ) )
    curves_all <- c( curve_check( cb_temp_full,"cb_temp", data$temp,b,V,temp_set$label ), curve_check(
      cb_rain_full,"cb_rain", data$rain,b,V,"Rain" ), curve_check( cb_vpd_full,"cb_vpd", data$VPD,b,V,"VPD"
      ), curve_check( cb_res_arr_full,"cb_res_arr", data$res_arr_scaled,b,V,"Res", probs=c(.10,.90) ) )
    rr_ci_names <- grep( "_(RR|LCL|UCL)_P[0-9]+$", names(rr_all), value=TRUE )
    curve_flags <- grep( "_Curve_Finite$", names(curves_all), value=TRUE )
    rr_ok <- length(rr_ci_names)>0 && all(is.finite(rr_all[rr_ci_names])) && all(rr_all[rr_ci_names]>0) &&
      length(curve_flags)==4 && all(curves_all[curve_flags]==1)
    if(!rr_ok)
      status <- "RR_OR_CURVE_FAIL"
  }
  
  control_attempts <- bind_rows( control_attempts, data.frame( grid_rank=g$grid_rank, model_id=g$model_id,
                                                               Grid_AIC=g$AIC, Refit_AIC=aic, status=status, Warning_count=length(fit_warnings), detail="" ) )
  write.csv( control_attempts, out("CONTROL_DISCOVERY_ATTEMPTS.csv"), row.names=FALSE )
  cat( "Refit AIC =",aic, "| Δ vs grid =",aic_diff, "| pdHess =",pd, "| VCOV =",v_ok, "| rank =",
       if(is.null(X)) NA else qr(X)$rank, "/", if(is.null(X)) NA else ncol(X), "| STATUS =",status, "\n" )
  
  # First fully usable model becomes positive control
  if(status=="OK"){
    control_fit <- fit
    control_grid_raw <- g
    control_b <- b
    control_V <- V
    control_X <- X
    control_rr <- rr_all
    control_curves <- curves_all
    control_warnings <- fit_warnings
    control_AIC <- aic
    control_BIC <- bic
    break
  }
  gc()
}
if(is.null(control_fit)){
  stop( "NO FULLY USABLE POSITIVE-CONTROL MODEL WAS FOUND." )
}

KNOWN_GOOD_ID <- as.integer( control_grid_raw$model_id )
ct <- c( control_grid_raw$lag_temp, control_grid_raw$dfvar_temp, control_grid_raw$dflag_temp )
cr <- c( control_grid_raw$lag_rain, control_grid_raw$dfvar_rain, control_grid_raw$dflag_rain )
cv <- c( control_grid_raw$lag_vpd, control_grid_raw$dfvar_vpd, control_grid_raw$dflag_vpd )
ca <- c( control_grid_raw$lag_res, control_grid_raw$dfvar_res, control_grid_raw$dflag_res )
control_path <- file.path( canonical_dir, paste0( "model_", KNOWN_GOOD_ID, "_positive_control.rds" ) )
saveRDS( control_fit, control_path )
cat( "\n============================================\n", "POSITIVE CONTROL FOUND\n",
     "============================================\n", "Model =",KNOWN_GOOD_ID,"\n",
     "Grid rank =",control_grid_raw$grid_rank,"\n", "AIC =",control_AIC,"\n", "BIC =",control_BIC,"\n",
     temp_set$label, ct[1],"/",ct[2],"/",ct[3],"\n", "Rain", cr[1],"/",cr[2],"/",cr[3],"\n", "VPD",
     cv[1],"/",cv[2],"/",cv[3],"\n", "Arrivals", ca[1],"/",ca[2],"/",ca[3],"\n" )

# 7) IDENTIFY CORE-VALID GRID MODELS
grid_core <- grid_results %>% filter( N==4113, convergence==0, pdHess==TRUE, is.finite(AIC), is.finite(BIC)
) %>% arrange(AIC,BIC) %>% mutate( grid_rank=row_number(), grid_delta_AIC=AIC-min(AIC),
                                   grid_delta_BIC=BIC-min(BIC) )
cat( "\nCore-valid grid models available =", nrow(grid_core),"\n" )

# 8) SAFE CANDIDATE WINDOW

# Control is already a valid model.

# Therefore:

# cannot possibly fall within Delta-AIC <= 10 of the

# best usable model.

# This replaces arbitrary TOP-50 / stop-after-10 rules.
screen_limit <- control_AIC + AIC_SCREEN_DELTA
candidate_pool <- grid_core %>% filter( AIC<=screen_limit+MATCH_TOL ) %>% arrange(AIC,BIC)
write.csv( candidate_pool, out( paste0( TAG, "_CANDIDATE_POOL.csv" ) ), row.names=FALSE )
cat( "\n============================================\n", "CANDIDATE SCREEN\n",
     "============================================\n", "Known-good control AIC =",control_AIC,"\n",
     "Screen boundary = control AIC +", AIC_SCREEN_DELTA, "=",screen_limit,"\n",
     "Grid candidates inside boundary =", nrow(candidate_pool),"\n", "NO top-N limit.\n",
     "NO stop-after-failures rule.\n", "NO stop-after-successes rule.\n" )

# 9) CREATE POSITIVE-CONTROL RESULT ROW
control_grid_rank <- NA_integer_
control_grid_AIC <- NA_real_
control_grid_BIC <- NA_real_
if(nrow(control_grid_raw)==1){
  control_grid_AIC <- control_grid_raw$AIC
  control_grid_BIC <- control_grid_raw$BIC
  tmp_rank <- grid_core %>% filter(model_id==KNOWN_GOOD_ID) %>% pull(grid_rank)
  if(length(tmp_rank)==1)
    control_grid_rank <- tmp_rank
}

control_pearson <- sum( residuals( control_fit, type="pearson" )^2 ) / df.residual(control_fit)
control_row <- data.frame( Grid_status=TAG, Candidate_source="POSITIVE_CONTROL", grid_rank=control_grid_rank,
                           model_id=KNOWN_GOOD_ID, lag_temp=ct[1], dfvar_temp=ct[2], dflag_temp=ct[3], lag_rain=cr[1],
                           dfvar_rain=cr[2], dflag_rain=cr[3], lag_vpd=cv[1], dfvar_vpd=cv[2], dflag_vpd=cv[3], lag_res=ca[1],
                           dfvar_res=ca[2], dflag_res=ca[3], Grid_AIC=control_grid_AIC, Grid_BIC=control_grid_BIC,
                           Refit_AIC=control_AIC, Refit_BIC=control_BIC, AIC_difference= ifelse( is.finite(control_grid_AIC),
                                                                                                 control_AIC-control_grid_AIC, NA_real_ ), BIC_difference= ifelse( is.finite(control_grid_BIC),
                                                                                                                                                                   control_BIC-control_grid_BIC, NA_real_ ), N=nobs(control_fit), convergence=control_fit$fit$convergence,
                           pdHess=isTRUE(control_fit$sdr$pdHess), Coefficients_finite=TRUE, VCOV_finite=TRUE,
                           X_rank=qr(control_X)$rank, X_columns=ncol(control_X), Full_rank=TRUE, NB1_alpha=sigma(control_fit),
                           District_SD=get_sd(control_fit,"state_name"), Year_SD=get_sd(control_fit,"year"),
                           Max_abs_gradient=get_grad(control_fit), Pearson=control_pearson, Warning_count=length(control_warnings),
                           Warning_text=paste( unique(control_warnings), collapse=" || " ), status="OK", Fit_error="",
                           Fit_path=control_path, check.names=FALSE )
control_row <- cbind( control_row, as.data.frame( as.list( c( control_rr, control_curves ) ),
                                                  check.names=FALSE ) )
write.csv( control_row, out("POSITIVE_CONTROL_RESULTS.csv"), row.names=FALSE )

# 10) REFIT EVERY COMPETITIVE GRID CANDIDATE
screen_results <- data.frame()
fit_store <- list()
fit_store[[as.character(KNOWN_GOOD_ID)]] <- control_fit
pool_refit <- candidate_pool %>% filter(model_id != KNOWN_GOOD_ID)
cat( "\nModels requiring canonical refit =", nrow(pool_refit), "\n" )
for(i in seq_len(nrow(pool_refit))){
  g <- pool_refit[i,]
  cat( "\n----------------------------------------------------\n", "Canonical refit",i,"/",nrow(pool_refit),
       "| grid rank",g$grid_rank, "| model",g$model_id, "| grid AIC",g$AIC,"\n",
       temp_set$label,g$lag_temp,g$dfvar_temp,g$dflag_temp, "| Rain",g$lag_rain,g$dfvar_rain,g$dflag_rain,
       "| VPD",g$lag_vpd,g$dfvar_vpd,g$dflag_vpd, "| Res",g$lag_res,g$dfvar_res,g$dflag_res, "\n" )
  cb_temp_full <- make_cb( data_full$temp, g$lag_temp, g$dfvar_temp, g$dflag_temp )
  cb_rain_full <- make_cb( data_full$rain, g$lag_rain, g$dfvar_rain, g$dflag_rain )
  cb_vpd_full <- make_cb( data_full$VPD, g$lag_vpd, g$dfvar_vpd, g$dflag_vpd )
  cb_res_arr_full <- make_cb( data_full$res_arr_scaled, g$lag_res, g$dfvar_res, g$dflag_res )
  
  # Exact common sample for model fitting
  cb_temp <- cb_temp_full[ common_rows,,drop=FALSE ]
  cb_rain <- cb_rain_full[ common_rows,,drop=FALSE ]
  cb_vpd <- cb_vpd_full[ common_rows,,drop=FALSE ]
  cb_res_arr <- cb_res_arr_full[ common_rows,,drop=FALSE ]
  fit_warnings <- character()
  fit <- withCallingHandlers( tryCatch( glmmTMB( local_cases ~ cb_temp + cb_rain + cb_vpd + cb_res_arr +
                                                   ns(time_index,df=11) + outbreak_level + lagged_local_cases + lag_imported + (1|state_name) +
                                                   (1|year), family=nbinom1(), offset=log(population), data=data_grid, control=ctrl ),
                                        error=function(e)e ), warning=function(w){
                                          fit_warnings <<- c( fit_warnings, conditionMessage(w) )
                                          invokeRestart("muffleWarning")
                                        } )
  base_row <- data.frame( Grid_status=TAG, Candidate_source="GRID_REFIT", grid_rank=g$grid_rank,
                          model_id=g$model_id, lag_temp=g$lag_temp, dfvar_temp=g$dfvar_temp, dflag_temp=g$dflag_temp,
                          lag_rain=g$lag_rain, dfvar_rain=g$dfvar_rain, dflag_rain=g$dflag_rain, lag_vpd=g$lag_vpd,
                          dfvar_vpd=g$dfvar_vpd, dflag_vpd=g$dflag_vpd, lag_res=g$lag_res, dfvar_res=g$dfvar_res,
                          dflag_res=g$dflag_res, Grid_AIC=g$AIC, Grid_BIC=g$BIC, stringsAsFactors=FALSE )
  if(inherits(fit,"error")){
    row <- base_row
    row$Refit_AIC <- NA_real_
    row$Refit_BIC <- NA_real_
    row$AIC_difference <- NA_real_
    row$BIC_difference <- NA_real_
    row$N <- NA_integer_
    row$convergence <- NA_integer_
    row$pdHess <- FALSE
    row$Coefficients_finite <- FALSE
    row$VCOV_finite <- FALSE
    row$X_rank <- NA_integer_
    row$X_columns <- NA_integer_
    row$Full_rank <- FALSE
    row$NB1_alpha <- NA_real_
    row$District_SD <- NA_real_
    row$Year_SD <- NA_real_
    row$Max_abs_gradient <- NA_real_
    row$Pearson <- NA_real_
    row$Warning_count <- length(fit_warnings)
    row$Warning_text <- paste( unique(fit_warnings), collapse=" || " )
    row$status <- "FIT_ERROR"
    row$Fit_error <- conditionMessage(fit)
    row$Fit_path <- NA_character_
    screen_results <- bind_rows( screen_results, row )
    write.csv( screen_results, out("CANONICAL_SCREEN_PROGRESS.csv"), row.names=FALSE )
    next
  }
  
  b <- fixef(fit)$cond
  V <- suppressWarnings( tryCatch( vcov(fit)$cond, error=function(e) NULL ) )
  X <- tryCatch( model.matrix(fit), error=function(e) NULL )
  n_fit <- nobs(fit)
  conv <- fit$fit$convergence
  pd <- isTRUE(fit$sdr$pdHess)
  aic <- suppressWarnings( AIC(fit) )
  bic <- suppressWarnings( BIC(fit) )
  aic_diff <- aic-g$AIC
  bic_diff <- bic-g$BIC
  n_ok <- n_fit==4113
  core_ok <- conv==0 && pd && is.finite(aic) && is.finite(bic)
  aic_match <- is.finite(aic_diff) && abs(aic_diff)<=MATCH_TOL
  bic_match <- is.finite(bic_diff) && abs(bic_diff)<=MATCH_TOL
  pipeline_match <- aic_match && bic_match
  coef_ok <- length(b)>0 && all(is.finite(b))
  v_ok <- !is.null(V) && length(V)>0 && all(is.finite(V))
  x_rank <- if(is.null(X))
    NA_integer_
  else
    qr(X)$rank
  x_cols <- if(is.null(X))
    NA_integer_
  else
    ncol(X)
  rank_ok <- !is.null(X) && x_rank==x_cols
  pearson <- tryCatch( sum( residuals( fit, type="pearson" )^2, na.rm=TRUE ) / df.residual(fit),
                       error=function(e) NA_real_ )
  
  # RR / CI / curve checks only once core inference is usable
  rr_all <- numeric(0)
  curves_all <- numeric(0)
  rr_ok <- FALSE
  if( n_ok && core_ok && pipeline_match && coef_ok && v_ok && rank_ok ){
    rr_temp <- rr_one( cb_temp_full, "cb_temp", data$temp, c(.05,.95), b,V, temp_set$label )
    rr_rain <- rr_one( cb_rain_full, "cb_rain", data$rain, c(.05,.95), b,V, "Rain" )
    rr_vpd <- rr_one( cb_vpd_full, "cb_vpd", data$VPD, c(.05,.95), b,V, "VPD" )
    rr_res <- rr_one( cb_res_arr_full, "cb_res_arr", data$res_arr_scaled, c(.05,.10,.90,.95), b,V, "Res" )
    rr_all <- c( rr_temp, rr_rain, rr_vpd, rr_res )
    curves_all <- c( curve_check( cb_temp_full, "cb_temp", data$temp, b,V, temp_set$label ), curve_check(
      cb_rain_full, "cb_rain", data$rain, b,V, "Rain" ), curve_check( cb_vpd_full, "cb_vpd", data$VPD, b,V,
                                                                      "VPD" ), curve_check( cb_res_arr_full, "cb_res_arr", data$res_arr_scaled, b,V, "Res", probs=c(.10,.90)
                                                                      ) )
    rr_ci_names <- grep( "_(RR|LCL|UCL)_P[0-9]+$", names(rr_all), value=TRUE )
    curve_flags <- grep( "_Curve_Finite$", names(curves_all), value=TRUE )
    rr_ok <- length(rr_ci_names)>0 && all(is.finite( rr_all[rr_ci_names] )) && all( rr_all[rr_ci_names]>0 ) &&
      length(curve_flags)==4 && all( curves_all[curve_flags]==1 )
  }
  
  status <- if(!n_ok){
    "N_MISMATCH"
  } else if(!core_ok){
    "CORE_FIT_FAIL"
  } else if(!pipeline_match){
    "PIPELINE_MISMATCH"
  } else if(!coef_ok){
    "NONFINITE_COEFFICIENT"
  } else if(!v_ok){
    "NONFINITE_VCOV"
  } else if(!rank_ok){
    "RANK_DEFICIENT"
  } else if(!rr_ok){
    "RR_OR_CURVE_FAIL"
  } else {
    "OK"
  }
  
  row <- base_row
  row$Refit_AIC <- aic
  row$Refit_BIC <- bic
  row$AIC_difference <- aic_diff
  row$BIC_difference <- bic_diff
  row$N <- n_fit
  row$convergence <- conv
  row$pdHess <- pd
  row$Coefficients_finite <- coef_ok
  row$VCOV_finite <- v_ok
  row$X_rank <- x_rank
  row$X_columns <- x_cols
  row$Full_rank <- rank_ok
  row$NB1_alpha <- tryCatch( sigma(fit), error=function(e) NA_real_ )
  row$District_SD <- get_sd( fit,"state_name" )
  row$Year_SD <- get_sd( fit,"year" )
  row$Max_abs_gradient <- get_grad(fit)
  row$Pearson <- pearson
  row$Warning_count <- length(fit_warnings)
  row$Warning_text <- paste( unique(fit_warnings), collapse=" || " )
  row$status <- status
  row$Fit_error <- ""
  fit_path <- NA_character_
  if(status=="OK"){
    fit_path <- file.path( canonical_dir, sprintf( "gridrank_%05d_model_%05d.rds", g$grid_rank, g$model_id ) )
    saveRDS( fit, fit_path )
    fit_store[[as.character(g$model_id)]] <- fit
  }
  row$Fit_path <- fit_path
  if(length(rr_all)){
    row <- cbind( row, as.data.frame( as.list( c( rr_all, curves_all ) ), check.names=FALSE ) )
  }
  
  screen_results <- bind_rows( screen_results, row )
  write.csv( screen_results, out("CANONICAL_SCREEN_PROGRESS.csv"), row.names=FALSE )
  cat( "Refit AIC =",aic, "| Δ vs grid =",aic_diff, "| pdHess =",pd, "| VCOV =",v_ok,
       "| rank =",x_rank,"/",x_cols, "| year SD =",get_sd(fit,"year"), "| STATUS =",status, "\n" )
  gc()
}

# 11) COMBINE GRID REFITS + POSITIVE CONTROL
all_screen <- bind_rows( control_row, screen_results )
all_screen <- all_screen %>% arrange( Refit_AIC, Refit_BIC )
write.csv( all_screen, out( paste0( TAG, "_CANONICAL_SCREEN_RESULTS.csv" ) ), row.names=FALSE )
cat("\nSCREEN STATUS SUMMARY:\n")
print( all_screen %>% count(status) )

# 12) USABLE MODELS + TRUE COMPETITIVE SET
usable <- all_screen %>% filter( status=="OK", N==4113, convergence==0, pdHess, Coefficients_finite,
                                 VCOV_finite, Full_rank, is.finite(Refit_AIC), is.finite(Refit_BIC) ) %>% arrange( Refit_AIC, Refit_BIC )
if(!nrow(usable)){
  stop( "NO USABLE MODEL EXISTS EVEN THOUGH ", "THE POSITIVE CONTROL PASSED. ", "THIS SHOULD NOT OCCUR." )
}

best_usable_AIC <- min( usable$Refit_AIC )
competitive <- usable %>% mutate( delta_AIC= Refit_AIC-best_usable_AIC ) %>% filter(
  delta_AIC<=AIC_SCREEN_DELTA+MATCH_TOL ) %>% mutate( delta_BIC_within_set= Refit_BIC-min(Refit_BIC),
                                                      AIC_within_2= delta_AIC<=2, AIC_within_4= delta_AIC<=4, AIC_within_10= delta_AIC<=10 ) %>% arrange(
                                                        Refit_AIC, Refit_BIC )
write.csv( competitive, out( paste0( TAG, "_COMPETITIVE_USABLE_MODELS.csv" ) ), row.names=FALSE )
cat( "\n============================================\n", "USABLE MODEL SUMMARY\n",
     "============================================\n", "Usable screened models =", nrow(usable),"\n",
     "Best usable AIC =", best_usable_AIC,"\n", "Competitive ΔAIC <= 2 =", sum(competitive$delta_AIC<=2),"\n",
     "Competitive ΔAIC <= 4 =", sum(competitive$delta_AIC<=4),"\n", "Competitive ΔAIC <= 10 =",
     nrow(competitive),"\n" )

# 13) FULL DIAGNOSTICS FOR COMPETITIVE USABLE MODELS
diagnose_fit <- function( fit, model_id ){
  rp <- residuals( fit, type="pearson" )
  pearson_ratio <- sum(rp^2,na.rm=TRUE) / df.residual(fit)
  pred <- fitted( fit, type="response" )
  obs <- model.response( model.frame(fit) )
  rmse <- sqrt( mean( (obs-pred)^2, na.rm=TRUE ) )
  mae <- mean( abs(obs-pred), na.rm=TRUE )
  r2 <- safe_cor( obs,pred )^2
  set.seed(123)
  sim <- tryCatch( DHARMa::simulateResiduals( fit, n=DHARMA_SIMS, refit=FALSE ), error=function(e)e )
  if(inherits(sim,"error")){
    summary_row <- data.frame( model_id=model_id, Pearson_dispersion_ratio=pearson_ratio,
                               DHARMa_p_dispersion=NA_real_, DHARMa_p_uniformity=NA_real_, DHARMa_p_zero_inflation=NA_real_,
                               DHARMa_p_outliers=NA_real_, Residual_Lag1_pooled=NA_real_, Residual_Lag1_mean=NA_real_,
                               Residual_Lag1_median=NA_real_, Residual_Lag1_max_abs=NA_real_, RMSE=rmse, MAE=mae, R2_pred=r2,
                               Diagnostic_error=conditionMessage(sim) )
    return( list( summary=summary_row, segment=data.frame() ) )
  }
  
  p_disp <- safe_p( DHARMa::testDispersion(sim) )
  p_unif <- safe_p( DHARMa::testUniformity(sim) )
  p_zero <- safe_p( DHARMa::testZeroInflation(sim) )
  p_out <- safe_p( DHARMa::testOutliers( sim, type="bootstrap" ) )
  resid_df <- data_grid %>% transmute( state_name, analysis_segment, time_index, DHARMa_resid=
                                         sim$scaledResiduals ) %>% arrange( state_name, analysis_segment, time_index ) %>% group_by(
                                           state_name, analysis_segment ) %>% mutate( previous_resid= lag(DHARMa_resid), previous_time=
                                                                                        lag(time_index), consecutive= time_index-previous_time==1 ) %>% ungroup()
  pairs <- resid_df %>% filter( consecutive, is.finite(DHARMa_resid), is.finite(previous_resid) )
  segment_acf <- pairs %>% group_by( state_name, analysis_segment ) %>% summarise( N_pairs=n(), Lag1_corr=
                                                                                     safe_cor( DHARMa_resid, previous_resid ), .groups="drop" ) %>% mutate( model_id=model_id, .before=1 )
  acf_vals <- segment_acf$Lag1_corr[ is.finite( segment_acf$Lag1_corr ) ]
  pooled_acf <- safe_cor( pairs$DHARMa_resid, pairs$previous_resid )
  summary_row <- data.frame( model_id=model_id, Pearson_dispersion_ratio= pearson_ratio, DHARMa_p_dispersion=
                               p_disp, DHARMa_p_uniformity= p_unif, DHARMa_p_zero_inflation= p_zero, DHARMa_p_outliers= p_out,
                             Residual_Lag1_pooled= pooled_acf, Residual_Lag1_mean= if(length(acf_vals)) mean(acf_vals) else NA_real_,
                             Residual_Lag1_median= if(length(acf_vals)) median(acf_vals) else NA_real_, Residual_Lag1_max_abs=
                               if(length(acf_vals)) max(abs(acf_vals)) else NA_real_, RMSE=rmse, MAE=mae, R2_pred=r2,
                             Diagnostic_error="" )
  list( summary=summary_row, segment=segment_acf )
}

diag_summary <- data.frame()
diag_segment <- data.frame()
cat( "\n============================================\n", "FULL DIAGNOSTICS\n",
     "============================================\n", "Competitive usable models to diagnose =",
     nrow(competitive),"\n" )
for(i in seq_len(nrow(competitive))){
  mid <- competitive$model_id[i]
  cat( "\nDiagnostics",i,"/", nrow(competitive), "| model",mid, "| ΔAIC", competitive$delta_AIC[i], "\n" )
  fit <- fit_store[[as.character(mid)]]
  if(is.null(fit)){
    fit_path <- competitive$Fit_path[ competitive$model_id==mid ][1]
    if( !is.na(fit_path) && file.exists(fit_path) ){
      fit <- readRDS(fit_path)
    }
  }
  if(is.null(fit)){
    d <- list( summary=data.frame( model_id=mid, Pearson_dispersion_ratio=NA_real_,
                                   DHARMa_p_dispersion=NA_real_, DHARMa_p_uniformity=NA_real_, DHARMa_p_zero_inflation=NA_real_,
                                   DHARMa_p_outliers=NA_real_, Residual_Lag1_pooled=NA_real_, Residual_Lag1_mean=NA_real_,
                                   Residual_Lag1_median=NA_real_, Residual_Lag1_max_abs=NA_real_, RMSE=NA_real_, MAE=NA_real_,
                                   R2_pred=NA_real_, Diagnostic_error="Fit object unavailable" ), segment=data.frame() )
  } else {
    d <- diagnose_fit( fit, mid )
  }
  
  diag_summary <- bind_rows( diag_summary, d$summary )
  diag_segment <- bind_rows( diag_segment, d$segment )
  write.csv( diag_summary, out("DIAGNOSTIC_PROGRESS.csv"), row.names=FALSE )
}
write.csv( diag_summary, out( paste0( TAG, "_COMPETITIVE_MODEL_DIAGNOSTICS.csv" ) ), row.names=FALSE )
write.csv( diag_segment, out( paste0( TAG, "_DISTRICT_SEGMENT_RESIDUAL_LAG1.csv" ) ), row.names=FALSE )

# 14) FINAL MODEL-SELECTION TABLE
selection_table <- competitive %>% left_join( diag_summary, by="model_id" ) %>% arrange( Refit_AIC, Refit_BIC
) %>% mutate( Selection_status= ifelse( complete_grid, "FINAL_GRID", "PROVISIONAL_PARTIAL_GRID" ) )
write.csv( selection_table, out( paste0( TAG, "_MODEL_SELECTION_TABLE.csv" ) ), row.names=FALSE )

# 15) CONSOLE SUMMARY
cat( "\n============================================\n", "MODEL-SELECTION SCREEN COMPLETE\n",
     "============================================\n", "Grid status:",TAG,"\n", "Grid rows available:",
     n_completed,"/",EXPECTED_GRID,"\n", "Control model AIC:", control_AIC,"\n", "Candidate screen boundary:",
     screen_limit,"\n", "Candidates canonical-refitted:", nrow(pool_refit),"\n", "Usable models:",
     nrow(usable),"\n", "Competitive usable models (ΔAIC <= 10):", nrow(competitive),"\n\n" )
cat("STATUS COUNTS:\n")
print( all_screen %>% count(status) )
cat("\nCOMPETITIVE MODEL TABLE:\n")
print( selection_table %>% select( Selection_status, grid_rank, model_id, lag_temp,dfvar_temp,dflag_temp,
                                   lag_rain,dfvar_rain,dflag_rain, lag_vpd,dfvar_vpd,dflag_vpd, lag_res,dfvar_res,dflag_res, Refit_AIC,
                                   delta_AIC, Refit_BIC, delta_BIC_within_set, VCOV_finite, Full_rank, Year_SD, Pearson_dispersion_ratio,
                                   DHARMa_p_dispersion, DHARMa_p_uniformity, DHARMa_p_zero_inflation, DHARMa_p_outliers,
                                   Residual_Lag1_pooled, Residual_Lag1_max_abs, RMSE, MAE, R2_pred ) )
if(!complete_grid){
  cat( "\n============================================\n", "IMPORTANT\n",
       "============================================\n", "This was a PARTIAL grid.\n",
       "The screen is useful for checking currently available ",
       "models, but DO NOT select the final model yet.\n", "Run this same Script 2 again when all ",
       EXPECTED_GRID," grid models are complete.\n" )
} else {
  cat( "\n============================================\n", "FULL GRID VERIFIED AND SCREENED\n",
       "============================================\n", "Use FINAL_MODEL_SELECTION_TABLE.csv ",
       "for the final candidate comparison.\n" )
}
cat( "\nFILES WRITTEN TO ",output_dir,"\n", "DONE\n", sep="" )

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
# 1) Load packages

suppressPackageStartupMessages({
  library(glmmTMB)
  library(dlnm)
  library(DHARMa)
  library(splines)
  library(dplyr)
  library(ggplot2)
  library(car)         
  library(corrplot)    
  library(mvtnorm)
})

set.seed(202)

# Temperature model
temp_var <- "tmin"   # change to "tmean" or "tmax"

temp_settings <- list(
  tmin  = list(label="Minimum temperature", lag=8, df_var=2, df_lag=3),
  tmean = list(label="Mean temperature",    lag=8, df_var=2, df_lag=3),
  tmax  = list(label="Maximum temperature", lag=12, df_var=4, df_lag=3)
)

temp_set <- temp_settings[[temp_var]]

# Separate outputs for each temperature model
output_dir <- file.path("outputs_synthetic", temp_var)
if(!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)
out <- function(x) file.path(output_dir, x)

# 2) Functions

compute_vif <- function(data, vars) {
  formula_str <- paste("local_cases ~", paste(vars, collapse = " + "))
  lm_model <- lm(as.formula(formula_str), data = data)
  print(round(car::vif(lm_model), 2))
}

evaluate_model <- function(model){
  
  mf <- model$frame
  id <- match(rownames(mf),rownames(data_grid))
  
  set.seed(123)
  sim <- DHARMa::simulateResiduals(model,n=1000,refit=FALSE)
  
  a <- data.frame(
    state=data_grid$state_name[id],
    segment=data_grid$analysis_segment[id],
    time=data_grid$time_index[id],
    resid=sim$scaledResiduals
  ) %>%
    arrange(state,segment,time) %>%
    group_by(state,segment) %>%
    mutate(prev=lag(resid),pt=lag(time)) %>%
    ungroup() %>%
    filter(time-pt==1,is.finite(resid),is.finite(prev))
  
  obs <- mf$local_cases
  pred <- fitted(model,type="response")
  rp <- residuals(model,type="pearson")
  
  list(
    AIC=round(AIC(model),3),
    BIC=round(BIC(model),3),
    Pearson_dispersion_ratio=round(sum(rp^2)/df.residual(model),4),
    DHARMa_p_dispersion=DHARMa::testDispersion(sim,plot=FALSE)$p.value,
    DHARMa_p_uniformity=DHARMa::testUniformity(sim,plot=FALSE)$p.value,
    DHARMa_p_zero_inflation=DHARMa::testZeroInflation(sim,plot=FALSE)$p.value,
    DHARMa_p_outliers=DHARMa::testOutliers(
      sim,type="bootstrap",nBoot=1000,
      plot=FALSE,plotBoostrap=FALSE)$p.value,
    Lag1_ACF_pooled=round(cor(a$resid,a$prev),4),
    R2_pred=cor(obs,pred)^2,
    RMSE=sqrt(mean((obs-pred)^2)),
    MAE=mean(abs(obs-pred))
  )
}

# MRV = interior local minimum (decreasing -> increasing); boundaries excluded

find_mrv <- function(x,rr,target=NULL){
  ok <- is.finite(x)&is.finite(rr)&rr>0; x <- x[ok]; y <- log(rr[ok])
  if(length(x)<5) return(NA_real_)
  d <- diff(y); k <- which(d[-length(d)]<0 & d[-1]>0)+1L
  if(!length(k)) return(NA_real_)
  
  if(is.null(target)) return(x[k[which.min(y[k])]])
  
  x[k[which.min(abs(x[k]-target))]]
}

mrv_point <- function(cb,b,V,x,within=c(.05,.95),ngrid=2001){
  lo <- quantile(x,within[1],na.rm=TRUE); hi <- quantile(x,within[2],na.rm=TRUE)
  at <- seq(lo,hi,length.out=ngrid); cen <- median(x[x>=lo & x<=hi],na.rm=TRUE)
  cp <- tryCatch(crosspred(cb,coef=b,vcov=V,model.link="log",at=at,cen=cen),
                 error=function(e) NULL)
  if(is.null(cp)) return(NA_real_)
  find_mrv(at,as.numeric(cp$allRRfit))
}

mrv_ci_se <- function(cb,b,V,x,within=c(.05,.95),ngrid=501,nsim=1000){
  hat <- mrv_point(cb,b,V,x,within,ngrid)
  lo <- quantile(x,within[1],na.rm=TRUE); hi <- quantile(x,within[2],na.rm=TRUE)
  at <- seq(lo,hi,length.out=ngrid); cen <- median(x[x>=lo & x<=hi],na.rm=TRUE)
  sim <- mvtnorm::rmvnorm(nsim,as.numeric(b),as.matrix(V))
  
  ms <- vapply(seq_len(nsim),function(i){
    cp <- tryCatch(crosspred(cb,coef=sim[i,],vcov=V,model.link="log",at=at,cen=cen),
                   error=function(e) NULL)
    if(is.null(cp)) return(NA_real_)
    z <- find_mrv(at,as.numeric(cp$allRRfit),
                  target=if(is.finite(hat)) hat else NULL)
    if(length(z)==1 && is.finite(z)) as.numeric(z) else NA_real_
  },numeric(1))
  
  v <- ms[is.finite(ms)]
  
  list(
    mrv=hat,
    ci=if(is.finite(hat)&&length(v)>10) unname(quantile(v,c(.025,.975))) else c(NA,NA),
    se=if(is.finite(hat)&&length(v)>1) sd(v) else NA_real_,
    valid_prop=mean(is.finite(ms)),
    draws=ms
  )
}

summarize_rr <- function(var_label, vec, cb, coef_all, vcov_all,
                         center_fn,
                         prob_pairs = list(c(0.05, 0.95), c(0.04, 0.96), c(0.06, 0.94))) {
  ind <- grep(paste0("^", var_label), names(coef_all)); stopifnot(length(ind) > 0)
  cen <- if (is.function(center_fn)) center_fn(vec, na.rm = TRUE) else as.numeric(center_fn)
  is_bad_row <- function(rr, low, high) {
    any(!is.finite(c(rr, low, high))) || (abs(rr - 1) < 1e-12 && abs(low - 1) < 1e-12 && abs(high - 1) < 1e-12)
  }
  used <- NULL; pred_ok <- NULL
  for (pair in prob_pairs) {
    p_lo <- unname(as.numeric(quantile(vec, pair[1], na.rm = TRUE)))
    p_hi <- unname(as.numeric(quantile(vec, pair[2], na.rm = TRUE)))
    if (is.finite(cen) && (abs(p_lo - cen) < 1e-12 || abs(p_hi - cen) < 1e-12)) next
    pred_try <- crosspred(cb, coef = coef_all[ind], vcov = vcov_all[ind, ind], model.link = "log", cen = cen, at = c(p_lo, p_hi))
    rr_lo <- unname(pred_try$allRRfit[as.character(p_lo)]); lcl_lo <- unname(pred_try$allRRlow[as.character(p_lo)]); ucl_lo <- unname(pred_try$allRRhigh[as.character(p_lo)])
    rr_hi <- unname(pred_try$allRRfit[as.character(p_hi)]); lcl_hi <- unname(pred_try$allRRlow[as.character(p_hi)]); ucl_hi <- unname(pred_try$allRRhigh[as.character(p_hi)])
    if (!is_bad_row(rr_lo, lcl_lo, ucl_lo) && !is_bad_row(rr_hi, lcl_hi, ucl_hi)) { used <- pair; pred_ok <- pred_try; break }
  }
  if (is.null(used)) return(tibble::tibble(variable = var_label, p_lo_prob = NA, p_hi_prob = NA))
  pull_row <- function(x, prob_used, side_lab) {
    rr <- unname(pred_ok$allRRfit[as.character(x)]); low <- unname(pred_ok$allRRlow[as.character(x)]); high <- unname(pred_ok$allRRhigh[as.character(x)])
    logRR <- log(rr); SE <- (log(high) - log(low)) / (2 * 1.96)
    z <- if (is.finite(SE) && SE > 0) logRR / SE else NA_real_
    p <- if (is.finite(z)) 2 * pnorm(abs(z), lower.tail = FALSE) else NA_real_
    tibble::tibble(
      !!paste0("p_", side_lab, "_prob") := prob_used, !!paste0("p_", side_lab, "_value") := x,
      !!paste0("RR_", side_lab) := rr, !!paste0("LCL_", side_lab) := low, !!paste0("UCL_", side_lab) := high, !!paste0("p_", side_lab) := p
    )
  }
  p_lo <- unname(as.numeric(quantile(vec, used[1], na.rm = TRUE)))
  p_hi <- unname(as.numeric(quantile(vec, used[2], na.rm = TRUE)))
  dplyr::bind_cols(
    tibble::tibble(variable = var_label, center_value = cen),
    pull_row(p_lo, used[1], "lo"),
    pull_row(p_hi, used[2], "hi")
  )
}

lagwise_rr <- function(var_label, vec, cb, coef_all, vcov_all,
                       center_fn, probs_used = NULL,
                       prob_pairs = list(c(0.05, 0.95), c(0.04, 0.96), c(0.06, 0.94)),
                       bylag = 0.25) {
  ind <- grep(paste0("^", var_label), names(coef_all)); stopifnot(length(ind) > 0)
  cen <- if (is.function(center_fn)) center_fn(vec, na.rm = TRUE) else as.numeric(center_fn)
  if (is.null(probs_used)) {
    for (pair in prob_pairs) {
      p_lo <- unname(as.numeric(quantile(vec, pair[1], na.rm = TRUE)))
      p_hi <- unname(as.numeric(quantile(vec, pair[2], na.rm = TRUE)))
      if (abs(p_lo - cen) > 1e-12 && abs(p_hi - cen) > 1e-12) { probs_used <- pair; break }
    }
  }
  if (is.null(probs_used)) stop("Could not find valid percentiles for lagwise_rr().")
  p_lo <- unname(as.numeric(quantile(vec, probs_used[1], na.rm = TRUE)))
  p_hi <- unname(as.numeric(quantile(vec, probs_used[2], na.rm = TRUE)))
  pred <- crosspred(cb, coef = coef_all[ind], vcov = vcov_all[ind, ind], model.link = "log", cen = cen, at = c(p_lo, p_hi), bylag = bylag)
  z_cols <- ncol(pred$matRRfit); raw_lag <- pred$lag
  lag_grid <- if (length(raw_lag) == z_cols) raw_lag else seq(min(raw_lag, na.rm = TRUE), max(raw_lag, na.rm = TRUE), length.out = z_cols)
  nearest_row <- function(x, label) {
    i <- which.min(abs(pred$predvar - x))
    tibble::tibble(
      variable = var_label, exposure_quantile = label,
      at_value = pred$predvar[i], lag = lag_grid,
      RR = as.numeric(pred$matRRfit[i, ]), LCL = as.numeric(pred$matRRlow[i, ]), UCL = as.numeric(pred$matRRhigh[i, ])
    ) |>
      dplyr::mutate(
        logRR = log(RR), SE = (log(UCL) - log(LCL)) / (2 * 1.96),
        z = dplyr::if_else(is.finite(SE) & SE > 0, logRR / SE, NA_real_),
        p = dplyr::if_else(is.finite(z), 2 * pnorm(abs(z), lower.tail = FALSE), NA_real_),
        p_used = ifelse(x == p_lo, probs_used[1], probs_used[2])
      )
  }
  dplyr::bind_rows(
    nearest_row(p_lo, paste0("p", round(probs_used[1] * 100))),
    nearest_row(p_hi, paste0("p", round(probs_used[2] * 100)))
  )
}

# Joint Wald test for full DLNM cross-basis
# Gasparrini & Armstrong: Wald test for multivariate association
wald_dlnm <- function(var_label, coef_all, vcov_all) {
  ind <- grep(paste0("^", var_label), names(coef_all))
  beta <- coef_all[ind]; V <- vcov_all[ind, ind, drop=FALSE]
  W <- as.numeric(t(beta) %*% solve(V, beta)); df <- length(beta)
  
  data.frame(
    variable=var_label, Wald_chisq=W, df=df,
    Wald_p=pchisq(W, df=df, lower.tail=FALSE)
  )
}

# Joint Wald test for lag-response curve at a specified exposure vs reference
wald_lagcurve <- function(var_label, cb, coef_all, vcov_all,
                          center_value, exposure_value, exposure_quantile) {
  ind <- grep(paste0("^", var_label), names(coef_all))
  stopifnot(length(ind)>0)
  
  red <- crossreduce(
    basis=cb, type="var", value=exposure_value,
    coef=coef_all[ind], vcov=vcov_all[ind,ind,drop=FALSE],
    model.link="log", cen=center_value
  )
  
  beta <- as.numeric(red$coefficients); V <- as.matrix(red$vcov)
  stopifnot(all(is.finite(beta)), all(is.finite(V)))
  W <- as.numeric(t(beta) %*% solve(V, beta)); df <- length(beta)
  
  data.frame(
    variable=var_label, exposure_quantile=exposure_quantile,
    exposure_value=exposure_value, reference_value=center_value,
    Joint_lagcurve_Wald_chisq=W, Joint_lagcurve_df=df,
    Joint_lagcurve_p=pchisq(W, df=df, lower.tail=FALSE)
  )
}

# Joint Wald test for cumulative exposure-response curve
wald_cumulative_curve <- function(var_label, cb, coef_all, vcov_all,
                                  center_value) {
  
  ind <- grep(paste0("^", var_label), names(coef_all))
  stopifnot(length(ind) > 0)
  
  red <- crossreduce(
    basis = cb,
    type = "overall",
    coef = coef_all[ind],
    vcov = vcov_all[ind, ind, drop = FALSE],
    model.link = "log",
    cen = center_value
  )
  
  beta <- as.numeric(red$coefficients)
  V <- as.matrix(red$vcov)
  
  stopifnot(all(is.finite(beta)), all(is.finite(V)))
  
  W <- as.numeric(t(beta) %*% solve(V, beta))
  df <- length(beta)
  
  data.frame(
    variable = var_label,
    Cumulative_curve_Wald_chisq = W,
    Cumulative_curve_df = df,
    Cumulative_curve_Wald_p = pchisq(W, df = df, lower.tail = FALSE)
  )
}


xlab_for <- function(var)
  switch(var,
         temp           = paste0(temp_set$label, " (°C)"),
         rain           = "Rainfall (mm)",
         VPD            = "VPD (kPa)",
         res_arr_scaled = "International resident arrivals (×1000)",
         var
  )

overall_title_for <- function(var) paste0("Cumulative RR vs ", xlab_for(var))
lag_title_for <- function(var) paste0("Lag-response: ", xlab_for(var))

# Cumulative exposure-response plot using exposure-specific reference values
plot_cumulative_rr_publication <- function(d, rr_row, var_label, save_path=NULL) {
  stopifnot(nrow(rr_row)==1)
  
  ref <- as.numeric(rr_row$center_value[1])
  low_lab <- sprintf("P%02d", round(100*rr_row$p_lo_prob[1]))
  high_lab <- sprintf("P%02d", round(100*rr_row$p_hi_prob[1]))
  ref_lab <- sprintf("Reference = %.2f", ref)
  
  pts <- tibble::tibble(
    Contrast=factor(c(low_lab,high_lab), levels=c(low_lab,high_lab)),
    value=c(rr_row$p_lo_value[1],rr_row$p_hi_value[1]),
    RR=c(rr_row$RR_lo[1],rr_row$RR_hi[1]),
    LCL=c(rr_row$LCL_lo[1],rr_row$LCL_hi[1]),
    UCL=c(rr_row$UCL_lo[1],rr_row$UCL_hi[1])
  )
  
  cols <- c("Cumulative RR"="black")
  cols[low_lab] <- "blue"; cols[high_lab] <- "red"
  ltypes <- c("RR = 1"="dashed"); ltypes[ref_lab] <- "dotted"
  
  p <- ggplot(d,aes(value,RR)) +
    geom_ribbon(aes(ymin=LCL,ymax=UCL,fill="95% CI"),alpha=.25,colour=NA) +
    geom_line(aes(colour="Cumulative RR"),linewidth=.9) +
    geom_hline(data=data.frame(y=1,lt="RR = 1"),aes(yintercept=y,linetype=lt),
               inherit.aes=FALSE,colour="black",linewidth=.6) +
    geom_vline(data=data.frame(x=ref,lt=ref_lab), aes(xintercept=x,linetype=lt),
               inherit.aes=FALSE,colour="black",linewidth=.8) +
    geom_vline(xintercept=pts$value,linetype="dotted",colour="grey45",
               linewidth=.45,show.legend=FALSE) +
    geom_errorbar(data=pts,aes(x=value,ymin=LCL,ymax=UCL,colour=Contrast),
                  inherit.aes=FALSE,width=0,linewidth=.6) +
    geom_point(data=pts,aes(x=value,y=RR,colour=Contrast),
               inherit.aes=FALSE,size=2.8) +
    scale_colour_manual(
      values=cols, breaks=c("Cumulative RR",low_lab,high_lab),
      labels=c("Cumulative RR",paste0(low_lab," contrast"),
               paste0(high_lab," contrast")), name=NULL
    ) +
    scale_fill_manual(values=c("95% CI"="grey80"),name=NULL) +
    scale_linetype_manual(values=ltypes,breaks=c("RR = 1",ref_lab),name=NULL) +
    labs(title=overall_title_for(var_label),x=xlab_for(var_label),
         y="Relative Risk (RR)") +
    theme_minimal(base_size=13) +
    theme(legend.position="bottom",legend.box="vertical",
          panel.grid.minor=element_blank()) +
    guides(colour=guide_legend(order=1),
           fill=guide_legend(order=2),
           linetype=guide_legend(order=3))
  
  print(p)
  if(!is.null(save_path))
    ggsave(save_path,p,width=7,height=5.5,dpi=600,bg="white")
  invisible(p)
}

# Lag-response plot at selected low/high exposures vs reference
plot_lagcurve_publication <- function(d, var_label, reference_value, save_path=NULL) {
  d <- d %>% mutate(Contrast=sprintf("P%02d",round(100*p_used)))
  lev <- sprintf("P%02d",round(100*sort(unique(d$p_used))))
  d$Contrast <- factor(d$Contrast,levels=lev)
  
  cols <- setNames(c("blue","red")[seq_along(lev)],lev)
  labs_leg <- setNames(paste0(lev," vs reference (line + 95% CI)"),lev)
  ref_lab <- sprintf("RR = 1 (reference = %.2f)",as.numeric(reference_value))
  
  p <- ggplot(d,aes(lag,RR,colour=Contrast,fill=Contrast)) +
    geom_ribbon(aes(ymin=LCL,ymax=UCL),alpha=.20,colour=NA) +
    geom_line(linewidth=1) +
    geom_hline(data=data.frame(y=1,lt=ref_lab),aes(yintercept=y,linetype=lt),
               inherit.aes=FALSE,colour="black",linewidth=.6) +
    scale_colour_manual(values=cols,breaks=lev,labels=labs_leg,
                        name="Exposure contrast (line + shaded 95% CI)") +
    scale_fill_manual(values=cols,breaks=lev,labels=labs_leg,
                      name="Exposure contrast (line + shaded 95% CI)") +
    scale_linetype_manual(values=setNames("dashed",ref_lab),name=NULL) +
    labs(title=lag_title_for(var_label),x="Lag (weeks)",y="Relative Risk (RR)") +
    theme_minimal(base_size=13) +
    theme(legend.position="bottom",legend.box="vertical",
          panel.grid.minor=element_blank()) +
    guides(colour=guide_legend(order=1),
           fill=guide_legend(order=1),
           linetype=guide_legend(order=2))
  
  print(p)
  if(!is.null(save_path))
    ggsave(save_path,p,width=7,height=5.5,dpi=600,bg="white")
  invisible(p)
}


# 3) Data preprocessing

csv_path <- Sys.getenv("DLNM_DATA",unset=NA)
if(is.na(csv_path)||!nzchar(csv_path)){if(interactive()) csv_path<-file.choose() else stop("No CSV specified.")}

data_full <- read.csv(csv_path,header=TRUE,check.names=TRUE)

stopifnot(temp_var %in% names(data_full))

data_full <- data_full %>%
  mutate(temp=as.numeric(.data[[temp_var]]),
         week=as.integer(week),year=as.integer(year),local_cases=as.integer(local_cases),
         imported_cases=as.numeric(imported_cases),res_arr_scaled=as.numeric(res_arr)/1000,
         rain=as.numeric(rain),VPD=as.numeric(VPD)) %>%
  filter(!is.na(population)&population>0,!is.na(local_cases),!is.na(temp),!is.na(rain),!is.na(VPD),
         !is.na(imported_cases),!is.na(res_arr),!is.na(state_name),year>=2014,year<=2024) %>%
  
  filter(year<2020|(year==2020&week<14)|(year==2021&week>40)|year>2021) %>%
  arrange(state_name,year,week) %>%
  group_by(state_name) %>% mutate(time_index=row_number()) %>% ungroup() %>%
  mutate(analysis_segment=if_else(year<=2020,"pre_gap","post_gap")) %>%
  group_by(state_name,analysis_segment) %>%
  mutate(lag_imported=lag(imported_cases,1),lagged_local_cases=lag(local_cases,1),
         outbreak_level=case_when(year%in%c(2020,2023,2024)~"high",
                                  year%in%c(2014,2015,2019)~"moderate",TRUE~"none")) %>%
  ungroup() %>%
  mutate(dlnm_group=interaction(state_name,analysis_segment,drop=TRUE),
         state_name=factor(state_name),year=factor(year),
         outbreak_level=factor(outbreak_level,levels=c("none","moderate","high")))

MAX_LAG <- 12
data_full <- data_full %>% group_by(dlnm_group) %>% mutate(grid_position=row_number()) %>% ungroup()
common_rows <- data_full$grid_position>MAX_LAG
data_grid <- droplevels(data_full[common_rows,,drop=FALSE])

stopifnot(nrow(data_full)==4329,nrow(data_grid)==4113)
stopifnot(all(complete.cases(data_grid[,c("local_cases","population","temp","rain","VPD","res_arr_scaled",
                                          "lagged_local_cases","lag_imported","time_index","outbreak_level",
                                          "state_name","year")])))
cat("\nFull rows =",nrow(data_full),"\nCommon-sample N =",nrow(data_grid),
    "\nCases =",sum(data_grid$local_cases),"\nDLNM segments =",nlevels(data_full$dlnm_group),"\n")


data <- data_grid


# Tables 3A and 3B
mean_sd <- function(x) sprintf("%.2f (%.2f)", mean(x), sd(x))
med_range <- function(x) sprintf("%.2f (%.2f-%.2f)", median(x), min(x), max(x))

by_district <- function(f) bind_rows(
  f(data_grid) %>% mutate(Admin="MAURITIUS", .before=1),
  f(group_by(data_grid, state_name)) %>% rename(Admin=state_name)
)

table3A <- by_district(function(d) summarise(d,
                                             across(c(local_cases, imported_cases),
                                                    list(Total=sum, Mean_SD=mean_sd, Median_Range=med_range,
                                                         Zero_pct=~round(mean(.x==0)*100,1))), .groups="drop"))

table3B <- by_district(function(d) summarise(d,
                                             across(c(temp, rain, VPD, res_arr_scaled),
                                                    list(Mean_SD=mean_sd, Median_Range=med_range)), .groups="drop"))

write.csv(table3A, out("Table_3A_Dengue_cases.csv"), row.names=FALSE)
write.csv(table3B, out("Table_3B_Climate_arrivals.csv"), row.names=FALSE)


# 4) Covariate checks

vars <- data %>% select(temp,rain,VPD,res_arr_scaled,lagged_local_cases,lag_imported,time_index)
cor_matrix <- cor(vars,use="complete.obs"); print(round(cor_matrix,2))
corrplot(cor_matrix,method="color",type="lower",addCoef.col="black",
         tl.col="black",number.cex=.8,diag=FALSE)

compute_vif(data,vars=c("temp","rain","VPD","res_arr_scaled",
                        "lagged_local_cases","lag_imported","time_index"))


# 5) Cross-basis setup

cb_specs <- list(
  temp=list(var="temp",lab="cb_temp",lag=temp_set$lag,
            df_var=temp_set$df_var,df_lag=temp_set$df_lag,within=c(.05,.95)),
  rain=list(var="rain",lab="cb_rain",lag=10,df_var=4,df_lag=3,within=c(.05,.95)),
  VPD=list(var="VPD",lab="cb_vpd",lag=12,df_var=2,df_lag=3,within=c(.05,.95)),
  res_arr=list(var="res_arr_scaled",lab="cb_res_arr",lag=12,df_var=2,df_lag=3,
               within=c(.10,.90))
)

make_cb <- function(x,lag,df_var,df_lag,group)
  crossbasis(x,lag=lag,argvar=list(fun="ns",df=df_var),
             arglag=list(fun="ns",df=df_lag),group=group)

for(nm in names(cb_specs)){
  s <- cb_specs[[nm]]
  cb_full <- make_cb(data_full[[s$var]],s$lag,s$df_var,s$df_lag,data_full$dlnm_group)
  assign(paste0(s$lab,"_full"),cb_full,envir=.GlobalEnv)
  assign(s$lab,cb_full[common_rows,,drop=FALSE],envir=.GlobalEnv)
}

stopifnot(
  all(complete.cases(cb_temp)), all(complete.cases(cb_rain)),
  all(complete.cases(cb_vpd)), all(complete.cases(cb_res_arr))
)

# 6) Fit DLNM-GLMM model

glmm_model_final <- glmmTMB(
  local_cases~cb_temp+cb_rain+cb_vpd+cb_res_arr+ns(time_index,df=11)+
    outbreak_level+lagged_local_cases+lag_imported+(1|state_name)+(1|year),
  family=nbinom1(),offset=log(population),data=data_grid,
  control=glmmTMBControl(
    optimizer=optim,optArgs=list(method="BFGS"),
    optCtrl=list(maxit=2000,reltol=1e-8))
)


cat("\n================ FIT CHECK ================\n",
    "N =",nobs(glmm_model_final),"\n",
    "AIC =",AIC(glmm_model_final),"\n",
    "BIC =",BIC(glmm_model_final),"\n",
    "Convergence =",glmm_model_final$fit$convergence,"\n",
    "pdHess =",isTRUE(glmm_model_final$sdr$pdHess),"\n",
    "NB1 alpha =",sigma(glmm_model_final),"\n",
    "===========================================\n")

V <- vcov(glmm_model_final)$cond
stopifnot(
  nobs(glmm_model_final)==4113,
  glmm_model_final$fit$convergence==0,
  isTRUE(glmm_model_final$sdr$pdHess),
  all(is.finite(V))
)

X <- model.matrix(glmm_model_final)
cat("VCOV =",nrow(V),"x",ncol(V),
    "| Non-finite =",sum(!is.finite(V)),
    "| Rank =",qr(X)$rank,"/",ncol(X),"\n")

print(summary(glmm_model_final))
print(VarCorr(glmm_model_final))


# 6A) Model diagnostics

final_eval <- evaluate_model(glmm_model_final)

print(final_eval)

set.seed(123)
sim_final <- DHARMa::simulateResiduals(glmm_model_final,n=2500,refit=FALSE)

png(out("DHARMa_diagnostics_600dpi.png"),
    width=7,height=7,units="in",res=600)
plot(sim_final)
dev.off()

model_no_crossbasis <- glmmTMB(
  local_cases ~  temp + rain + VPD + res_arr_scaled + ns(time_index,df=11) +
    outbreak_level + lagged_local_cases + lag_imported +
    (1|state_name) + (1|year),
  family=nbinom1(),
  offset=log(population),
  data=data_grid,
  control=glmmTMBControl(
    optimizer=optim,
    optArgs=list(method="BFGS"),
    optCtrl=list(maxit=2000,reltol=1e-8)
  )
)

AIC(glmm_model_final,model_no_crossbasis)
BIC(glmm_model_final,model_no_crossbasis)

evaluate_model(model_no_crossbasis)

# 7) Wald tests

coef_all <- fixef(glmm_model_final)$cond
vcov_all <- vcov(glmm_model_final)$cond

wald_results <- bind_rows(
  wald_dlnm("cb_temp",coef_all,vcov_all),
  wald_dlnm("cb_rain",coef_all,vcov_all),
  wald_dlnm("cb_vpd",coef_all,vcov_all),
  wald_dlnm("cb_res_arr",coef_all,vcov_all)
)

print(wald_results)


full_cb <- function(s) get(paste0(s$lab,"_full"))
median_ref <- setNames(sapply(cb_specs,\(s) median(data[[s$var]],na.rm=TRUE)),
                       sapply(cb_specs,`[[`,"var"))

# Cumulative exposure-response curve Wald tests

cumulative_wald_results <- bind_rows(
  lapply(cb_specs, function(s) {
    
    wald_cumulative_curve(
      var_label = s$lab,
      cb = full_cb(s),
      coef_all = coef_all,
      vcov_all = vcov_all,
      center_value = median_ref[s$var]
    )
    
  })
)

print(cumulative_wald_results)


# 8) MRV estimation + uncertainty check for ALL exposures

MRV <- lapply(cb_specs,function(s){
  cb <- full_cb(s); j <- grep(paste0("^",s$lab),names(coef_all))
  z <- mrv_ci_se(cb,coef_all[j],vcov_all[j,j,drop=FALSE],
                 data[[s$var]],within=s$within,ngrid=501,nsim=1000)
  
  if(is.finite(z$mrv)){
    p <- crosspred(cb,coef=coef_all[j],vcov=vcov_all[j,j,drop=FALSE],
                   model.link="log",cen=median_ref[s$var],at=z$mrv)
    rr <- c(p$allRRfit[1],p$allRRlow[1],p$allRRhigh[1])
  } else rr <- rep(NA_real_,3)
  
  list(
    summary=data.frame(
      variable=s$var,MRV=z$mrv,LCL=z$ci[1],UCL=z$ci[2],SE=z$se,
      Valid_prop=z$valid_prop,median_ref=median_ref[s$var],
      RR_MRV_vs_median=rr[1],RR_LCL=rr[2],RR_UCL=rr[3]),
    draws=z$draws
  )
})

mrv_table <- do.call(rbind,lapply(MRV,`[[`,"summary"))
print(mrv_table)


# PRIMARY REFERENCES
# Rainfall = MRV
# Temperature, VPD, arrivals = observed median


primary_ref <- median_ref

rain_mrv <- mrv_table$MRV[mrv_table$variable == "rain"]

stopifnot(length(rain_mrv) == 1, is.finite(rain_mrv))

primary_ref["rain"] <- rain_mrv

print(primary_ref)

# 8A) DLNM contour plots

plot_contour <- function(s,fn){
  cb <- full_cb(s); j <- grep(paste0("^",s$lab),names(coef_all))
  cen <- if(s$var=="rain") mrv_table$MRV[mrv_table$variable=="rain"] else median(data[[s$var]],na.rm=TRUE)
  
  x <- data[[s$var]]
  q <- quantile(x,s$within,na.rm=TRUE)
  at <- seq(q[1],q[2],length.out=100)
  
  p <- crosspred(cb,coef=coef_all[j],vcov=vcov_all[j,j,drop=FALSE],
                 model.link="log",cen=cen,at=at,bylag=.25)
  
  z <- p$matRRfit; lev <- pretty(z,20)
  cols <- c(colorRampPalette(c("#006400","gold"))(sum(lev<1)),
            colorRampPalette(c("gold","red"))(sum(lev>1)))
  
  png(out(paste0(fn,"_contour_600dpi.png")),7,5.5,units="in",res=600)
  filled.contour(p$predvar,seq(min(p$lag),max(p$lag),length.out=ncol(z)),z,
                 levels=lev,col=cols,
                 plot.title=title(main=paste("Relative Risk (RR)\n",xlab_for(s$var)),
                                  xlab=xlab_for(s$var),ylab="Lag (weeks)",font.main=2),
                 key.title=title(main="RR",font.main=2))
  dev.off()
}

for(nm in names(cb_specs)){
  s <- cb_specs[[nm]]
  plot_contour(s,if(nm=="temp") temp_var else nm)
}


# MRV simulation distributions
png(out("INTERIOR_MRV_simulation_distributions.png"),
    width=9,height=7,units="in",res=300)

par(mfrow=c(2,2))
for(i in seq_along(MRV)){
  x <- MRV[[i]]$draws
  x <- x[is.finite(x)]
  hist(x,breaks=30,
       main=names(cb_specs)[i],
       xlab="Simulated MRV")
  abline(v=mrv_table$MRV[i],lwd=2,lty=2)
}
par(mfrow=c(1,1))
dev.off()


# 9) Exposure-specific RR estimates

summ_one <- function(s){
  probs <- if(s$var=="res_arr_scaled") list(c(.10,.90)) else list(c(.05,.95))
  summarize_rr(
    s$lab,data[[s$var]],full_cb(s),coef_all,vcov_all,
    center_fn=unname(primary_ref[s$var]),prob_pairs=probs
  )
}

rr_list <- lapply(cb_specs,summ_one)
rr_overall <- bind_rows(rr_list)

print(rr_overall)

# 9A) Lag-curve wald tests

wald_pair <- function(s,low_lab,high_lab){
  r <- rr_overall %>% filter(variable==s$lab)
  
  bind_rows(
    wald_lagcurve(
      s$lab,full_cb(s),coef_all,vcov_all,
      r$center_value,r$p_lo_value,low_lab
    ),
    wald_lagcurve(
      s$lab,full_cb(s),coef_all,vcov_all,
      r$center_value,r$p_hi_value,high_lab
    )
  )
}

lagcurve_wald_results <- bind_rows(
  wald_pair(cb_specs$temp,"P05","P95"),
  wald_pair(cb_specs$rain,"P05","P95"),
  wald_pair(cb_specs$VPD,"P05","P95"),
  wald_pair(cb_specs$res_arr,"P10","P90")
)

print(lagcurve_wald_results)


# 9B) Table 4 — Relative risk estimates

fmt_p <- function(p) ifelse(is.na(p), "p=NA",
                            ifelse(p<.001, "p<0.001", sprintf("p=%.3f",p)))

report_table <- rr_overall %>%
  left_join(cumulative_wald_results, by="variable") %>%
  transmute(
    Variable=dplyr::recode(variable, cb_temp=temp_set$label, cb_rain="Rain",
                    cb_vpd="VPD", cb_res_arr="ResArr (x1000)"),
    Reference=if_else(variable=="cb_rain",
                      sprintf("%.2f (median); %.2f (local minimum)",median_ref["rain"],rain_mrv),
                      sprintf("%.2f (median)",center_value)),
    p05=sprintf("P%02d: %.2f",round(p_lo_prob*100),p_lo_value),
    RR_p05=sprintf("%.2f (%.2f-%.2f), %s",RR_lo,LCL_lo,UCL_lo,fmt_p(p_lo)),
    p95=sprintf("P%02d: %.2f",round(p_hi_prob*100),p_hi_value),
    RR_p95=sprintf("%.2f (%.2f-%.2f), %s",RR_hi,LCL_hi,UCL_hi,fmt_p(p_hi)),
    Overall_cumulative_Wald=sprintf("%.2f (%d df), %s",
                                    Cumulative_curve_Wald_chisq,Cumulative_curve_df,fmt_p(Cumulative_curve_Wald_p))
  )

print(report_table)
write.csv(report_table, out("Table_4_DLNMs_RR.csv"), row.names=FALSE)


# 10) Cumulative RR curves
curve_one <- function(s){
  x <- data[[s$var]]
  at <- seq(quantile(x,s$within[1],na.rm=TRUE),
            quantile(x,s$within[2],na.rm=TRUE),length.out=201)
  cen <- unname(primary_ref[s$var])
  
  cb <- full_cb(s); j <- grep(paste0("^",s$lab),names(coef_all))
  p <- crosspred(cb,coef=coef_all[j],vcov=vcov_all[j,j,drop=FALSE],
                 model.link="log",cen=cen,at=at)
  
  data.frame(variable=s$lab,value=at,Reference=cen,
             RR=p$allRRfit,LCL=p$allRRlow,UCL=p$allRRhigh)
}

dlnm_curves <- bind_rows(lapply(cb_specs,curve_one))
write.csv(dlnm_curves,out("overall_RR_curve_primary_reference.csv"),row.names=FALSE)

rr_curve <- bind_rows(lapply(cb_specs,function(s){
  probs <- if(s$var=="res_arr_scaled") list(c(.10,.90)) else list(c(.05,.95))
  summarize_rr(s$lab,data[[s$var]],full_cb(s),coef_all,vcov_all,
               center_fn=unname(primary_ref[s$var]),prob_pairs=probs)
}))

for(nm in names(cb_specs)){
  s <- cb_specs[[nm]]
  fn <- if(nm=="temp") temp_var else nm
  
  plot_cumulative_rr_publication(
    filter(dlnm_curves,variable==s$lab),
    filter(rr_curve,variable==s$lab),
    s$var,
    out(paste0(fn,"_cumulative_RR_600dpi.png"))
  )
}

# 11) Lag-response curves

probs_used <- lapply(rr_list,\(x) c(x$p_lo_prob[1],x$p_hi_prob[1]))

lag_one <- function(s,p){
  cen <- if(s$var=="rain") mrv_table$MRV[mrv_table$variable=="rain"] else median_ref[s$var]
  lagwise_rr(s$lab,data[[s$var]],full_cb(s),coef_all,vcov_all,
             center_fn=cen,probs_used=p)
}

lag_rr_all <- bind_rows(Map(lag_one,cb_specs,probs_used))
write.csv(lag_rr_all,out("lagwise_RR_primary_reference.csv"),row.names=FALSE)

for(nm in names(cb_specs)){
  s <- cb_specs[[nm]]
  fn <- if(nm=="temp") temp_var else nm
  
  plot_lagcurve_publication(
    filter(lag_rr_all,variable==s$lab),
    s$var,
    if(s$var=="rain") mrv_table$MRV[mrv_table$variable=="rain"] else median_ref[s$var],
    out(paste0(fn,"_lag_RR_600dpi.png"))
  )
}


# 12) Observed vs Predicted plots

nat_pred <- data.frame(
  year=data_grid$year,
  week=data_grid$week,
  observed=data_grid$local_cases,
  predicted=fitted(glmm_model_final,type="response")
) %>%
  group_by(year,week) %>%
  summarise(
    observed=sum(observed),
    predicted=sum(predicted),
    .groups="drop"
  ) %>%
  mutate(time=as.numeric(as.character(year))+(week-1)/52.18)

# 1) Observed vs predicted through time
p_time <- ggplot(nat_pred,aes(time)) +
  geom_line(aes(y=observed,colour="Observed"),linewidth=1.2) +
  geom_line(aes(y=predicted,colour="Predicted"),linewidth=1.2,linetype="dashed") +
  annotate(
    "label",x=2014.3,y=Inf,vjust=1.2,hjust=0,
    label=sprintf(
      "District-week model fit:\nR² = %.2f\nRMSE = %.2f",
      final_eval$R2_pred,
      final_eval$RMSE
    )
  ) +
  scale_x_continuous(breaks=2014:2024) +
  labs(
    title="Observed vs Predicted Weekly Dengue Cases",
    x="Time (weeks, 2014–2024)",
    y="Number of cases",
    colour=NULL
  ) +
  theme_minimal(base_size=13) +
  theme(legend.position="top")

print(p_time)

ggsave(
  out("Observed_vs_Predicted_weekly.png"),
  p_time,width=9,height=5.5,dpi=600,bg="white"
)

# 2) Observed vs predicted scatter

d <- data.frame(
  observed=glmm_model_final$frame$local_cases,
  predicted=fitted(glmm_model_final,type="response")
)

p_scatter <- ggplot(d,aes(predicted,observed)) +
  geom_point(alpha=.35) +
  geom_abline(slope=1,intercept=0,linetype="dashed",colour="red") +
  geom_smooth(method="lm",se=FALSE) +
  annotate(
    "label",x=Inf,y=Inf,hjust=1.05,vjust=1.2,
    label=sprintf("R² = %.2f\nRMSE = %.2f",
                  final_eval$R2_pred,
                  final_eval$RMSE)
  ) +
  labs(
    title="Observed vs Predicted (district-week)",
    x="Predicted weekly cases",
    y="Observed weekly cases"
  ) +
  theme_minimal(base_size=13)

print(p_scatter)

ggsave(
  out("Observed_vs_Predicted_scatter.png"),
  p_scatter,width=6,height=5.4,dpi=600,bg="white"
)
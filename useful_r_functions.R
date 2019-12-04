## Useful R functions
library(RandomFields)
library(FactoMineR)
library(spatstat)


# Estimate optimal range parameter in a gp smooth
optimal_range <- function(min_dist, max_dist, length.out = 100, model_data, k=-1){
  
REML <- r <- seq(min_dist, max_dist, length.out = length.out)
for (i in seq_along(r)) {
  m <- gam(cbind(n_pos, n_neg) ~ s(x, y, k = k, bs = "gp", m = c(3, r[i])), 
           family = "binomial", data = model_data, method = "REML")
  REML[i] <- m$gcv.ubre
}
return(list(REML = REML,
            best_m = r[which.min(REML)]))
}

## Get exceedance probabilities
exceedance_prob <- function(gam_model, prediction_data, n_sims, threshold){
  Cg <- predict(gam_model, prediction_data, type = "lpmatrix")
  sims <- rmvn(n_sims, mu = coef(gam_model), V = vcov(gam_model, unconditional = TRUE))
  fits <- Cg %*% t(sims)
  
  # Add residual error
  #error_samp <- sample(resid(gam_model), nrow(fits)*ncol(fits), replace = TRUE)
  fits <- fits + error_samp
  fits_prev <- exp(fits) / (1 + exp(fits))
  apply(fits_prev, 1, function(x) {sum(x>threshold)/n_sims})
}


posterior_mean <- function(gam_model, prediction_data, n_sims){
  
  Cg <- predict(gam_model, prediction_data, type = "lpmatrix")
  sims <- rmvn(n_sims, mu = coef(gam_model), V = vcov(gam_model, unconditional = TRUE))
  fits <- Cg %*% t(sims)
  fits_prev <- exp(fits) / (1 + exp(fits))
  apply(fits_prev, 1, mean)
}

validate_posterior <- function(gam_model, prediction_data, n_sims, prob_threshold, prob_width){
  
  Cg <- predict(gam_model, prediction_data, type = "lpmatrix")
  pred_mean <- predict(gam_model, prediction_data, type = "response")
  sims <- rmvn(n_sims, mu = coef(gam_model), V = vcov(gam_model, unconditional = TRUE))
  fits <- Cg %*% t(sims)
  
  # Add residual error
  #error_samp <- sample(gam_model$residuals, nrow(fits)*ncol(fits), replace = TRUE)
  #error_samp <- rnorm(nrow(fits)*ncol(fits), sd(gam_model$residuals))
  fits <- fits + error_samp
  fits_prev <- exp(fits) / (1 + exp(fits))
  
  # Calc prevalence threshold from probability threshold
  prev_threshold <- apply(fits_prev, 1, function(x){quantile(x, prob = prob_threshold)})
  exceedance_perf <- mean(prediction_data$prev >= prev_threshold)
  
  # Calc proportion of true prevalence values fall within the middle 
  # prob_width
  quant_probs <- c((1 - prob_width)/2, (1 - (1 - prob_width)/2))
  prev_quantiles <- apply(fits_prev, 1, function(x){quantile(x, prob = quant_probs)})
  posterior_perf <- mean(prediction_data$prev >= prev_quantiles[1,] & prediction_data$prev <= prev_quantiles[2,])
  return(list(exceedance_perf = exceedance_perf,
              posterior_perf = posterior_perf,
              prev_quantiles = t(prev_quantiles)))
}


validate_posterior_spaMM <- function(spaMM_model, prediction_data, n_sims, prob_threshold, prob_width){
  
  fits_prev <- simulate(spaMM_model, 
                   type = "(ranef|response)", 
                   nsim = n_sims,
                   newdata =prediction_data)
  
  # Convert to prevalence
  fits_prev <- fits_prev/100
  
  # Calc prevalence threshold from probability threshold
  prev_threshold <- apply(fits_prev, 1, function(x){quantile(x, prob = prob_threshold)})
  exceedance_perf <- mean(prediction_data$prev >= prev_threshold)
  
  # Calc proportion of true prevalence values fall within the middle 
  # prob_width
  quant_probs <- c((1 - prob_width)/2, (1 - (1 - prob_width)/2))
  prev_quantiles <- apply(fits_prev, 1, function(x){quantile(x, prob = quant_probs)})
  posterior_perf <- mean(prediction_data$prev >= prev_quantiles[1,] & prediction_data$prev <= prev_quantiles[2,])
  return(list(exceedance_perf = exceedance_perf,
              posterior_perf = posterior_perf,
              prev_quantiles = t(prev_quantiles)))
}

# Simulate risk
simulate_risk <- function(seed, var, scale, mean, nrow=256, ncol=256){ 
  
  # Simluate a risk surface based on elevation + GP
  set.seed(seed)
  model <- RMexp(var=var, scale=scale)
  
  set.seed(seed)
  simu <- RandomFields::RFsimulate(model, x=1:ncol, 
                     y=1:nrow, RFoptions(spConform=FALSE))

  # Convert to raster
  simu_raster <- raster(nrows = nrow, ncol = ncol, xmn=0, xmx=1, ymn=0, ymx=1)
  simu_raster[] <- as.vector(simu)
  
  # Add mean
  log_odds_raster <- mean + simu_raster 

  # Convert to probability
  prev_raster <- exp(log_odds_raster) / (1 + exp(log_odds_raster))
  names(prev_raster) <- "prev_raster"
  return(prev_raster)
}




## Simulate some villages
simulate_villages <- function(seed, kappa, scale, mu, ref_raster){

  set.seed(seed)
  owin_area <- owin(extent(ref_raster)[1:2], extent(ref_raster)[3:4])
  set.seed(seed)
  villages_pp <- rMatClust(kappa, scale, mu, win=owin_area)
  
  # Convert to SPDF
  villages_spdf <- SpatialPointsDataFrame(cbind(villages_pp$x, villages_pp$y),
                                          data=data.frame(id=1:villages_pp$n,
                                                          x=villages_pp$x,
                                                          y=villages_pp$y))
  villages_spdf <- cbind(villages_spdf, raster::extract(ref_raster, villages_spdf))
  names(villages_spdf) <- c("id", "x", "y", names(ref_raster))
  return(villages_spdf)
}


# Take initial sample
initial_survey <- function(seed, villages_SP, n, n_ind){

  set.seed(seed)
  selected_villages <- sample(1:nrow(villages_SP), n)
  survey_villages <- villages_SP[selected_villages,]
  survey_villages$n_pos <- rbinom(length(selected_villages), n_ind, survey_villages$prev)
  survey_villages$n_neg <- n_ind - survey_villages$n_pos
  return(survey_villages)
}


# Generate stack of raster with lat and lng
gen_pred_stack <- function(ref_raster){
  
  x <- y <- ref_raster
  x[] <- coordinates(x)[,1]
  y[] <- coordinates(y)[,2]
  stacked <- stack(x, y)
  names(stacked) <- c("x", "y")
  return(stacked)
}




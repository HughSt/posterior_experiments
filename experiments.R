library(mgcv)
library(geostatsp)
RFsimulate <- RandomFields::RFsimulate # REquired as geostatsp masks this function

# Experimenting with different spatial prediction methods for 
# binomial data. Source useful functions
source("https://raw.githubusercontent.com/HughSt/posterior_experiments/master/useful_r_functions.R")

# Simluate some data
risk_raster <- simulate_risk(seed=1981, 
              var=0.5, 
              scale=50, 
              mean= -2)
names(risk_raster) <- "prev"


# Simulate some villages
villages <- simulate_villages(seed=1981, 
                  kappa=10, 
                  scale=0.2, 
                  mu=200, 
                  ref_raster = risk_raster)
plot(risk_raster)
points(villages)

# Take a sample
model_data <- initial_survey(1981, villages, n=100, n_ind=100)

# Fit GAM model. First calculate optimal 'range' parameter of a
# gp model (using matern covariance model)
REML_estimates <- optimal_range(min_dist = 0.01, 
              max_dist = max(dist(model_data@data[,c("x", "y")])),
              model_data = model_data)

# fit models
gam_mod_gp <- mgcv::gam(cbind(n_pos, n_neg) ~ s(x, y, k=-1, bs="gp", m = c(3,REML_estimates$best_m)),
                        data = model_data, family="binomial", method = "REML")



# Predict
pred_raster <- gen_pred_stack(risk_raster)
predictions <- predict(pred_raster, gam_mod_gp, type="response")
plot(predictions)

# Get exceedance probabilities
exceedance_probs <- exceedance_prob(gam_mod_gp, as.data.frame(coordinates(pred_raster)), 500, 0.1)
exceedance_probs_raster <- predictions
exceedance_probs_raster[] <- exceedance_probs
plot(exceedance_probs_raster)

# Validate the posterior in 2 ways. 1) Calculate the proportion of times 
# the 'true' prevalence value was above threshold and compare to 
# exceedance probabilities and 2) calc the proportion of times
# the true prevalence value lies within defined quantile range of the posterior 
validation_results <- list(exceedance_perf=NULL,
                           posterior_perf=NULL)
for(i in 1:99){
  validation <- validate_posterior(gam_mod_gp, 
                     villages@data, 
                     n_sims = 500, 
                     prob_threshold = i/100,
                     prob_width = i/100)

  validation_results$exceedance_perf <- c(validation_results$exceedance_perf,
                                          validation$exceedance_perf)
  validation_results$posterior_perf <- c(validation_results$posterior_perf,
                                         validation$posterior_perf)
}

# Plot exceedance probabilities
plot(1 - (1:99/100), validation_results$exceedance_perf ,
     xlab = "Exceedance probability",
     ylab = "Proportion correct")

# Plot coverage
plot(1:99/100, validation_results$posterior_perf,
     xlab = "Centred prediction quantile",
     ylab = "Proportion correct")

# plot 95% prediction intervals
prediction_interval <- validate_posterior(gam_mod_gp, 
                   villages@data, 
                   n_sims = 500, 
                   prob_threshold = NULL,
                   prob_width = 0.95)

# First get an order for plotting
plot_order <- order(villages@data$prev)

# Plot the first n points with confidence intervals
plot(villages@data$prev[plot_order],cex=0.2,pch=16, xlim=c(1,100))
for(i in 1:nrow(prediction_interval$prev_quantiles)){
  lines(rep(i, 2), c(prediction_interval$prev_quantiles[plot_order[i],1],
                     prediction_interval$prev_quantiles[plot_order[i],2]), col="gray80")
}; points(villages@data$prev[plot_order],cex=0.2,pch=16)

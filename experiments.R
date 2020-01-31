library(mgcv)
library(RANN)
library(RandomFields)
library(FactoMineR)
library(spatstat)
library(raster)

# Experimenting with different spatial prediction methods for 
# binomial data. Source useful functions
source("https://raw.githubusercontent.com/HughSt/posterior_experiments/master/useful_r_functions.R")

# Simluate some data
risk_raster <- simulate_risk(seed=1, 
              var=0.5, 
              scale=50, 
              mean= -2)
names(risk_raster) <- "prev"


# Simulate some villages
villages <- simulate_villages(seed=1, 
                  kappa=10, 
                  scale=0.2, 
                  mu=200, 
                  ref_raster = risk_raster)
plot(risk_raster)
points(villages)

# Take a sample
model_data <- initial_survey(1, villages, n=100, n_ind=100)

# Fit GAM model. First calculate optimal 'range' parameter of a
# gp model (using matern covariance model)
REML_estimates <- optimal_range(min_dist = 0.01, 
                                k = 99,
              max_dist = max(dist(model_data@data[,c("x", "y")])),
              model_data = model_data)

# fit model with optimal range
gam_mod_gp_1 <- mgcv::gam(cbind(n_pos, n_neg) ~ 
                          s(x, y, bs="gp", 
                            k=99,
                            m = c(3,REML_estimates$best_m)),
                          #method = "REML",
                          #method = "GCV"
                        data = model_data, family="binomial")

# Plot obs v exp
plot(model_data$prev, predict(gam_mod_gp_1, type="response"))

# Plot across a single lng
pred_sing_lng <- data.frame(x = rep(mean(model_data$x), 500), 
                            y = seq(min(model_data$y), max(model_data$y), length.out = 500))
plot(pred_sing_lng$y, predict(gam_mod_gp_1, pred_sing_lng, type="response"),
     type="l", lwd=2)


# Predict on grid
pred_raster <- gen_pred_stack(risk_raster)
predictions_1 <- predict(pred_raster, gam_mod_gp_1, type="response")
plot(predictions_1)

# Get exceedance probabilities
exceedance_probs <- exceedance_prob(gam_mod_gp_1, 
                                    as.data.frame(coordinates(pred_raster)), 
                                    n_sims = 500, 
                                    threshold = 0.5)
exceedance_probs_raster <- predictions_1
exceedance_probs_raster[] <- exceedance_probs
plot(exceedance_probs_raster)

# Validate the posterior in 2 ways. 1) Calculate the proportion of times 
# the 'true' prevalence value was above threshold and compare to 
# exceedance probabilities and 2) calc the proportion of times
# the true prevalence value lies within defined quantile range of the posterior 
# First give the prediction_points a cluster
villages@data$cluster <- kmeans(coordinates(villages), 50)$cluster
validation_results <- list(exceedance_perf=NULL,
                           posterior_perf=NULL,
                           cluster_exceedance_perf=NULL,
                           cluster_posterior_perf=NULL)
for(i in 1:99){
  validation <- validate_posterior(gam_mod_gp_1, 
                     villages@data, 
                     n_sims = 1000, 
                     prob_threshold = i/100,
                     prob_width = i/100)

  validation_results$exceedance_perf <- c(validation_results$exceedance_perf,
                                          validation$exceedance_perf)
  validation_results$posterior_perf <- c(validation_results$posterior_perf,
                                         validation$posterior_perf)
  validation_results$cluster_exceedance_perf <- c(validation_results$cluster_exceedance_perf,
                                          validation$cluster_exceedance_perf)
  validation_results$cluster_posterior_perf <- c(validation_results$cluster_posterior_perf,
                                         validation$cluster_posterior_perf)
}

# Plot exceedance probabilities
plot(1:99/100, validation_results$exceedance_perf ,
     xlab = "Exceedance probability",
     ylab = "Proportion correct")
abline(1,-1, col="blue", lwd=2)

# Plot coverage
plot(1:99/100, validation_results$posterior_perf,
     xlab = "Centred prediction quantile",
     ylab = "Proportion correct")
abline(0,1, col="blue", lwd=2)

# Plot exceedance probabilities for clusters
plot(1:99/100, validation_results$cluster_exceedance_perf,
     xlab = "Exceedance probability",
     ylab = "Proportion correct")
abline(1,-1, col="blue", lwd=2)

# Plot coverage for clusters
plot(1:99/100, validation_results$cluster_posterior_perf,
     xlab = "Centred prediction quantile",
     ylab = "Proportion correct")
abline(0,1, col="blue", lwd=2)

# plot 95% prediction intervals
prediction_interval <- validate_posterior(gam_mod_gp_1, 
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




############################### ############################### ############################### 
############################### try spaMM  ############################### 
############################### ############################### ############################### 
library(spaMM)
lfit <- fitme(cbind(n_pos, n_neg) ~
                   Matern(1|x+y),
                  fixed = list(nu = 3/2),
                  data=model_data@data,
                  family=binomial())

pred_raster <- gen_pred_stack(risk_raster)
predictions_lfit <- predict(pred_raster, lfit, type="response")
plot(predictions_lfit)

# Plot for single lat
plot(pred_sing_lng$y, predict(lfit, pred_sing_lng, type="response"),
     type="l", lwd=2)

sims <- simulate(lfit, 
                 type = "(ranef|response)", 
                 nsim = 200,
                 newdata = villages@data)
hist(sims[1,])
pal <- colorNumeric(tim.colors(), c(0,60))
plot(villages$x, villages$y, col = pal(sims[,1]), pch= 16)

res <- validate_posterior_spaMM(lfit, villages@data, 100, 0.2, 0.2)
res$posterior_perf

# Validate the Posterior of the spaMM model. Get a cup of tea...
validation_results_spaMM <- list(exceedance_perf=NULL,
                           posterior_perf=NULL)

prev_cutoffs <- seq(1,99,5)
for(i in prev_cutoffs){
  # validation_spaMM <- validate_posterior_spaMM(lfit, 
  #                                  villages@data, 
  #                                  n_sims = 500, 
  #                                  prob_threshold = i/100,
  #                                  prob_width = i/100)
  
                        validation_spaMM <- validate_posterior_spaMM(lfit, 
                                               villages@data, 
                                               n_sims = 500, 
                                               prob_threshold = i/100,
                                               prob_width = i/100)
  
  validation_results_spaMM$exceedance_perf <- c(validation_results_spaMM$exceedance_perf,
                                                validation_spaMM$exceedance_perf)
  validation_results_spaMM$posterior_perf <- c(validation_results_spaMM$posterior_perf,
                                               validation_spaMM$posterior_perf)
}

# Plot exceedance probabilities
plot(prev_cutoffs/100, validation_results_spaMM$exceedance_perf ,
     xlab = "Exceedance probability",
     ylab = "Proportion correct")

# Plot coverage
plot(prev_cutoffs/100, validation_results_spaMM$posterior_perf,
     xlab = "Centred prediction quantile",
     ylab = "Proportion correct"); abline(0,1)

# plot 95% prediction intervals
prediction_interval_spaMM <-  validate_posterior_spaMM(lfit, 
                                                 villages@data, 
                                                 n_sims = 500, 
                                                 prob_threshold = NULL,
                                                 prob_width = 0.95)


# Plot the first n points with confidence intervals
plot(villages@data$prev[plot_order],cex=0.2,pch=16)
for(i in 1:nrow(prediction_interval_spaMM$prev_quantiles)){
  lines(rep(i, 2), c(prediction_interval_spaMM$prev_quantiles[plot_order[i],1],
                     prediction_interval_spaMM$prev_quantiles[plot_order[i],2]), col="gray80")
}; points(villages@data$prev[plot_order],cex=0.2,pch=16)



### Try RF
# For observation irst calc distance to nearest n
# points and their respecitve values
nn <- nn2(model_data@data[,c("x", "y")],
    model_data@data[,c("x", "y")],
    k=50)

model_data@data <- cbind(model_data@data,
                         nn$nn.dists[,-1],
                         matrix((model_data@data$n_pos/100)[nn$nn.idx[,-1]],
                                nrow = nrow(model_data)))
                         
names(model_data@data)[7:55] <- paste0("dist_near", 1:49)
names(model_data@data)[56:104] <- paste0("prev_near", 1:49)

# Fit model
Y <- factor(c(rep(0, nrow(model_data@data)),
              rep(1, nrow(model_data@data))))

rf_formula <- as.formula(paste("Y ~", paste0("dist_near", 1:49, collapse = "+"),
                                   "+", 
                                   paste0("prev_near", 1:49, collapse = "+")))

rf_mod <- ranger(rf_formula,
                 seed = 1981,
                 data = model_data@data,
                 probability = TRUE,
                 case.weights = c(model_data@data$n_neg,
                                  model_data@data$n_pos))
preds <- predict(rf_mod, model_data@data, type = "response")
pred_df_rf <- coordinates(pred_raster[[1]])

nn_pred <- nn2(model_data@data[,c("x", "y")],
               pred_df_rf,
          k=50)

pred_df_rf <- as.data.frame(cbind(nn_pred$nn.dists[,-1],
                         matrix((model_data@data$n_pos/100)[nn_pred$nn.idx[,-1]],
                                nrow = nrow(pred_df_rf))))

names(pred_df_rf)[1:49] <- paste0("dist_near", 1:49)
names(pred_df_rf)[50:98] <- paste0("prev_near", 1:49)

pred_for_raster_rf <- predict(rf_mod, pred_df_rf, type = "response")
pred_raster_rf <- pred_raster[[1]]
pred_raster_rf[] <- pred_for_raster_rf$predictions[,2]
plot(pred_raster_rf)




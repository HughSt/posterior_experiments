library(mgcv)
library(geostatsp)
library(RANN)
RFsimulate <- RandomFields::RFsimulate # REquired as geostatsp masks this function

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
              max_dist = max(dist(model_data@data[,c("x", "y")])),
              model_data = model_data)

# fit models
gam_mod_gp_1 <- mgcv::gam(cbind(n_pos, n_neg) ~ 
                          s(x, y, bs="gp", 
                            #m = c(3, 0.03)),
                            m = c(3,REML_estimates$best_m)),
                          #s(as.factor(model_data$id), bs="re"),
                        data = model_data, family="binomial")

gam_mod_gp_2 <- mgcv::gamm(cbind(n_pos, n_neg) ~ 
                           s(x, y, bs="gp", 
                             m = c(3, 0.03)),
                         #m = c(3,REML_estimates$best_m)),
                         #s(id, bs="re"),
                         correlation=corGaus(.1,form=~x+y),
                         data = model_data@data, family="binomial")
plot(model_data$prev, predict(gam_mod_gp_1, type="response"))
points(model_data$prev, predict(gam_mod_gp_2$gam, type="response"), pch=16)

# Fit GP model with INLA using the geostatsp package
# glgm_mod <- glgm(formula = n_pos ~ 1, 
#      data = model_data, grid = 150, family = "binomial",
#      Ntrials = (model_data$n_neg + model_data$n_neg), shape = 1, buffer = 0.1,
#      priorCI = list(sd = c(0.2, 4), range = c(0.5, 5e+05)))

# Predict
pred_raster <- gen_pred_stack(risk_raster)
predictions <- predict(pred_raster, gam_mod_gp_1, type="response")
plot(predictions)

# Get exceedance probabilities
exceedance_probs <- exceedance_prob(gam_mod_gp_1, as.data.frame(coordinates(pred_raster)), 500, 0.1)
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
  validation <- validate_posterior(gam_mod_gp_1, 
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
plot(1:99/100, validation_results$exceedance_perf ,
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




############################### ############################### ############################### 
############################### try spaMM  ############################### 
############################### ############################### ############################### 
library(spaMM)
lfit <- fitme(cbind(n_pos, n_neg) ~
                   Matern(1|x+y),
                  data=model_data@data,
                  family=binomial())

pred_raster <- gen_pred_stack(risk_raster)
predictions_lfit <- predict(pred_raster, lfit, type="response")
plot(predictions_lfit)

sims <- simulate(lfit, 
                 type = "(ranef|response)", 
                 nsim = 200,
                 newdata = villages@data)
hist(sims[1,])
pal <- colorNumeric(tim.colors(), c(0,60))
plot(villages$x, villages$y, col = pal(sims[,1]), pch= 16)

res <- validate_posterior_spaMM(lfit, villages@data, 100, 0.2, 0.2)


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
  
                        validation_spaMM <- vp(lfit, 
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
     ylab = "Proportion correct")

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




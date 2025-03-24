##################################
## Predicting pollutant levels ###
##  from methylation profiles  ###
##################################

# Based on Dr Eamonn Mallon's epigenetic clock script
## HOUSEKEEPING ##
# Remove unused objects and free up memory
rm(list = ls())
gc()  # Force garbage collection
options(max.print=1000)  # Limit console print output
options(stringsAsFactors = FALSE)  # Set default for character data to avoid memory issues with factors

set.seed(123)  # Set a seed for reproducibility

predict_pollutant <- function(pollutant) {
  
  message("Welcome! This ML machine aims to predict levels of various pollutants using methylation profiles.")
  
  valid_pollutants <- c("pesticide", "herbicide", "insecticide", "fungicide")
  if (!pollutant %in% valid_pollutants) {
    stop(paste0("Invalid pollutant. Must be one of: ", paste(valid_pollutants, collapse = ", ")))
  }
  
  ## PACKAGE LOADING ##
  library(dplyr)
  library(purrr)
  library(tidyverse)
  library(ggfortify)
  library(glmnet)
  library(ggplot2)
  library(ggpubr)
  library(matrixStats)
  library(vegan)
  library(modeest)
  library(Metrics)
  library(vip)
  library(viridis)
  library(broom)
  library(lme4)
  library(emmeans)
  library(data.table)
  library(glmmTMB)
  library(RColorBrewer)
  library(scales)
  library(doParallel)
  library(foreach)
  library(caret)
  library(report)
  library(stringr)
  library(rlang)
  
  setwd(paste0("~/final_project/github/ML/", pollutant))
  pollutant_title <- str_to_title(pollutant)
  
  # =============================== ML TRAINING  =========================================
  # This script will:
  # (1) Pre-process methylation data
  # (2) Trains 'Elastic Net' regression model
  # (3) Evaluates ability to predict chemical pollutant (herbicide, pesticide, fungicide, insecticide) levels.
  
  # ========== Step 1: Load & Prepare Data ==========
  # Loading pollutant training data:
  data <- read.csv(paste0(pollutant, "_lvl_ml_training_data.csv"))
  
  # Separating predictor and response variable
  methylation_lvl <- as.matrix(data[,-1]) #just methylation profiles; predictor
  pollutant_lvl <- data[, 1] #just pollutant levels; response variable

  # ========== Step 2: Define Training Control ==========
  train_control <- trainControl(
    method = "repeatedcv",      # repeated K-fold cross validation
    number = 10,                # 10-fold
    repeats = 10,               # Repeat CV 10 times for stability
    savePredictions = "final",  # Save cross-validated predictions
    returnResamp = "final"      # Return re-sampling results
  )
  
  # ========== Step 3: Set Up Parallel Processing ==========
  num_cores <- detectCores() -1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # ========== Step 4: Standardize Features (Center & Scale) ==========
  # Compute pre-processing parameters using training data
  preProc_control <- preProcess(methylation_lvl, method = c("center", "scale"))
  
  # Apply standardisation to the data
  meth_scaled <- predict(preProc_control, methylation_lvl)
  
  # ========== Step 5: Train Elastic Net Model ==========
  # Define a rough lambda search using cross-validation
  cv_fit <- cv.glmnet(as.matrix(meth_scaled), pollutant_lvl, alpha = 0.5, nfolds = 5)
  
  # Extract best lambda estimates
  best_lambda <- cv_fit$lambda.min
  lambda_1se <- cv_fit$lambda.1se
  
  # Fine-tune lambda grid
  # This is the final, started much wider and refined
  refined_lambda_grid <- exp(seq(log(0.01), log(1), length.out = 50))
  tuneGrid <- expand.grid(alpha = 0.5, lambda = refined_lambda_grid)
  print("First few combinations of the lambda grid:")
  print(head(tuneGrid))  # Check first few combinations
  #     alpha lambda
  # 1   0.5 0.01000000
  # 2   0.5 0.01098541
  # 3   0.5 0.01206793
  # 4   0.5 0.01325711
  # 5   0.5 0.01456348
  # 6   0.5 0.01599859
  
  # Train cross validation model with refined hyper-parameters
  cv_model <- train(
    x = meth_scaled,  # Use PCA-transformed data
    y = pollutant_lvl,
    method = "glmnet",
    trControl = train_control,
    tuneGrid = tuneGrid,
    metric = "Rsquared"  # Use RMSE for regressionRsquared
  )
  
  # Evaluate cross-validation results
  cv_results <- cv_model$results
  print("Cross validation model results:")
  print(cv_results)
  write.csv(cv_results, file = paste0(pollutant,"_cross_validation_results.csv"))
  
  plot(cv_model)
  pdf(file = paste0(pollutant, "_cross_validation_model.pdf"), paper = "a4r")
  print(plot(cv_model))
  dev.off()
  
  # Extract best hyper-parameters
  best_lambda <- cv_model$bestTune$lambda
  best_alpha <- cv_model$bestTune$alpha
  print("Best hyper-parameters:")
  print(cv_model$bestTune)
  
  # Final model trained on full dataset with best parameters
  final_model <- train(
    x = meth_scaled,  # Use PCA-transformed data
    y = pollutant_lvl,
    method = "glmnet",
    tuneGrid = expand.grid(alpha = best_alpha, lambda = best_lambda),
    trControl = trainControl(method = "none")
  )
  
  # Stop parallel processing
  stopCluster(cl)
  registerDoSEQ()
  
  # ========== Step 6: Make Predictions ==========
  y_pred <- predict(final_model, newdata = meth_scaled) 
  
  # Calculate performance metrics
  MAE <- mean(abs(y_pred - pollutant_lvl))
  MSE <- mean((y_pred - pollutant_lvl)^2)
  RMSE <- sqrt(MSE)
  R2 <- cor(y_pred, pollutant_lvl)^2
  
  # Report metrics
  print(paste("MAE:", MAE))
  print(paste("RMSE:", RMSE))
  print(paste("R²:", R2))
  
  # Merge predictions back into full dataset
  graph_data <- data %>%
    mutate(!!sym(paste0("pred_", pollutant)) := y_pred)
  
  # ========== Step 7: Visualize Predictions ==========
  actual_pollutant_col <- pollutant
  predicted_pollutant_col <- paste0("pred_", pollutant)
  
  actual_predicted_plot <- ggplot(graph_data, aes(x = .data[[actual_pollutant_col]], y = .data[[predicted_pollutant_col]])) +
  geom_point(color = "blue") +        # Scatter plot of actual vs. predicted
   geom_smooth(method = "lm", color = "red", se = FALSE) +  # Regression line
   labs(
     title = paste0(pollutant_title, " Actual vs Predicted Values"),
    x = "Actual Values",
    y = "Predicted Values"
   ) +
   theme_minimal()
  
  ggsave(paste0(pollutant, "_actual_vs_predicted.pdf"), 
         plot =  actual_predicted_plot, 
         width = 297, height = 210, units = "mm") #A4 landscape
  
  prediction_bypop_plot <- ggplot(graph_data, aes(x = .data[[actual_pollutant_col]], y = .data[[predicted_pollutant_col]], color = population, shape = population)) +
    geom_jitter(width = 0.4, height = 0.4, alpha = 0.7, size = 2.5) +  # Increase jitter width/height, reduce size
    geom_smooth(method = "lm", se = FALSE, linetype = "solid", linewidth = 1.2) +  # Use linewidth for line thickness
    labs(
      title = paste0(pollutant_title, " Levels Predictions by Population"),
      x = paste0("Actual ",pollutant," Levels"),
      y = paste0("Predicted",pollutant,"Levels"),
      color = "Lake Ring Population",  # Change legend title to be more descriptive
      shape = "Population"   # Add shape legend title
    ) +
    scale_color_brewer(palette = "Set1") +  # Use a colorblind-friendly palette
    scale_shape_manual(values = c(16, 17, 18, 15)) +  # Specify 4 distinct shapes for each treatment
    theme_classic(base_size = 14) +  # Use a classic theme with base font size adjustment
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),  # Center and bold title
      axis.title = element_text(size = 14),  # Increase axis title font size
      axis.text = element_text(size = 12),   # Increase axis text font size
      legend.title = element_text(size = 14),  # Adjust legend title font size
      legend.text = element_text(size = 12),   # Adjust legend text font size
      legend.position = "top"  # Move legend to the top
    )
  
  ggsave(paste0(pollutant, "_predictions_bypop.pdf"), 
         plot =  prediction_bypop_plot, 
         width = 297, height = 210, units = "mm")
  
  # ========== Step 8: Save & Export Results ==========
  saveRDS(final_model, file = paste0(pollutant, "_final_model.rds"))
  
  # Extract coefficients from the final glmnet model within the caret object
  final_glmnet <- final_model$finalModel
  
  # Extract non-zero coefficients for the selected lambda value
  non_zero_coefs <- as.data.frame(as.matrix(coef(final_glmnet, s = best_lambda)))
  non_zero_coefs <- non_zero_coefs[non_zero_coefs != 0, , drop = FALSE]  # Only keep non-zero coefficients
  
  # Exporting to CSV with row names
  write.csv(non_zero_coefs, paste0(pollutant, "_loci.csv"), row.names = TRUE)
  
  # Exporting to CSV all the data from model output
  write_csv(graph_data, paste0(pollutant, "_graph_data.csv"))
  
  # ========== Step 9: Feature Importance Analysis ==========
  # Extract feature importance (which CpGs)
  feature_importance <- varImp(final_model, scale = FALSE)  # Get importance
  feature_df <- data.frame(
    Feature = rownames(feature_importance$importance),
    Importance = feature_importance$importance[, 1]
  )
  
  # Sort by importance (highest first)
  feature_df <- feature_df %>%
    arrange(desc(Importance)) %>%
    slice(1:30)  # Top 30 features (adjust as needed)
  
  # Plot with ggplot2
  features_importance_plot <- ggplot(feature_df, aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_col(fill = "steelblue") +
    coord_flip() +  # Horizontal bars
    theme_minimal(base_size = 14) +
    labs(
      title = paste0(pollutant_title, " Lvl Feature Importance (Top 30)"),
      x = "Feature",
      y = "Importance"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 10)
    )
  
  ggsave(paste0(pollutant, "_feature_importance_top30.pdf"), 
         plot =  features_importance_plot, 
         width = 297, height = 210, units = "mm")
  
  # ========== Step 10: Analysing the results of the model ==========
  # Dynamically create the formula based on the pollutant input
  formula_anova <- as.formula(paste0("pred_", pollutant, " ~ ", pollutant, " * population"))
  
  analysis<-lm(formula_anova, data = graph_data)
  anova(analysis)
  
  # Filter the graph_data for the Pesticide group
  pesticide_data <- subset(graph_data, population == "PP")
  
  # Fit a linear model using only the Pesticide data
  formula_lm <- as.formula(paste0("pred_", pollutant, " ~ ", pollutant))
  
  pesticide_model <- lm(formula_lm, data = pesticide_data)
  
  # Summarize the model to see the significance of the pollutant factor
  summary(pesticide_model)
  
  # Average performance metrics across all resamples
  print("Mean performance metrics from bootstrap resampling:")
  mean_results <- cv_model$results
  print(mean_results)
  
  # Check if Rsquared column is valid
  if (all(is.na(mean_results$Rsquared))) {
    print("All Rsquared values are NA. No best model metrics available.")
  } else {
    # Use which.max to find the row index with the best Rsquared
    best_model_index <- which.max(mean_results$Rsquared)
    
    # Extract the best model's metrics
    best_model_metrics <- mean_results[best_model_index, ]
    
    # Print the best model metrics
    print("Best model metrics:")
    print(best_model_metrics)
  }
  
  message("Complete! Thank you for using this script. -NM")
}

#Select from herbicide, pesticide, insecticide or fungicide
pollutant <- readline(prompt = "Enter the pollutant name (pesticide, herbicide, insecticide, fungicide): ")

predict_pollutant(pollutant)

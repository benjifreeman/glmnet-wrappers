
tuned_glmnet_cont <- function(data.train,
                              data.test = NULL,
                              n.folds = 10,
                              n.repeats = 5,
                              n.levels = c(100,25),
                              cv.metric = "mae",
                              seed = 18,
                              outcome.var){
    
    # Predicts a continous outcome (specified by `outcome.var`)
    # w/ cross-validated linear regression on training data (`data.train`)
    # If testing data (`data.test`) is provided, assesses models w/ standard metrics
    
    require(tidymodels)
    
    # set seed for reproducibility
    set.seed(seed)
    
    if (is.null(data.test)){
        message("No testing data provided - was this intentional?")
    } else {
        # Check that training and testing data have the same variables
        if (!setequal(colnames(data.train),colnames(data.test))){
            stop("training and testing data do not have the same variables")
        }
        
        # Reconstruct test/train split (if testing data was provided)
        data <- bind_rows(data.train, data.test)
        test.train.indices <- list(analysis = seq(nrow(data.train)),
                                   assessment = nrow(data.train) + seq(nrow(data.test)))
        
        data.split <- make_splits(test.train.indices, data)
        
    }

    # Define preprocessing steps:
        # assign roles to variables
        # normalize all numeric variables (convert to z-score)
        # convert nominal predictors to dummy variables
    data.rec <- recipe(data.train) %>% 
        update_role(eval(outcome.var),new_role = "outcome") %>% 
        update_role(-all_outcomes(),new_role = "predictor") %>% 
        update_role(ID,new_role = "id variable") %>% 
        step_normalize(all_numeric()) %>% 
        step_dummy(all_nominal_predictors())
        
    # Define model specification:
        # linear regression (continuous outcome)
        # penalty and mixture will be cross-validated over
        # engine is elastic net regression (glmnet)
    tune.spec <- linear_reg(penalty = tune(), mixture = tune()) %>% 
        set_engine("glmnet")
    
    # Create grid for penalty and mixture values:
    lambda.grid <- grid_regular(extract_parameter_set_dials(tune.spec),levels=n.levels)
    
    # Set up folds for cross-validation:
    cv.folds <- vfold_cv(data.train, v = n.folds, repeats = n.repeats)
    
    # Specify possibly cross-validation metrics:
    model.metrics <- metric_set(mae,yardstick::rmse,mape)
    
    # Define workflow (specify preprocessor and mode type):
    wf <- workflow() %>% 
        add_recipe(data.rec) %>% 
        add_model(tune.spec)
    
    # Cross-validate parameters over the grid:
    model.grid <- tune_grid(
        wf,
        resamples = cv.folds,
        grid = lambda.grid,
        metrics = model.metrics)
    
    # Collect cross-validation metrics
    cv.metrics <- model.grid %>% 
        collect_metrics()
    
    # Choose best parameters
    best.regs <- model.grid %>% 
        select_best(metric = cv.metric)
    
    # Finalize model with best penalty and mixture
    final.model <- finalize_workflow(
        wf,
        best.regs)
    
    # Fit model to the training data
    trained.fit <- final.model %>% 
        fit(data.train)
    
    # Extract model coefficients
    model.coeffs <- trained.fit %>% 
        extract_fit_parsnip() %>%
        tidy() %>% 
        arrange(desc(abs(estimate)))
    
    # If there is no testing data, stop here.
    if (is.null(data.test)){
        
        return(list(fit = model.coeffs,
                    best.params = best.regs,
                    final.model = trained.fit,
                    cv.scores = cv.metrics))
        
        # Otherwise, proceed to testing
    } else if (!is.null(data.test)){
        
        
        # Assess predictions on the testing data
        fit.metrics <- last_fit(final.model,
                                data.split,
                                metrics = metric_set(mae,rsq)) %>%
            collect_metrics()
        
        return(list(cv.scores = cv.metrics,
                    best.params=best.regs,
                    final.model=trained.fit,
                    fit=model.coeffs,
                    metrics=fit.metrics))
    }
}






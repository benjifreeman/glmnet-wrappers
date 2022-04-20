tuned_glmnet_bin <- function(data.train,
                             data.test = NULL,
                             n.folds = 10,
                             n.repeats = 5,
                             n.levels = c(100,25),
                             cv.metric = "accuracy",
                             outcome.var){
    
    # Binarizes a continous outcome (specified by `outcome.var`)
    # and performs cross-validated logistic regression on training data (`data.train`)
    # if testing data (`data.test`) is provided, provides ROC curve and assesses models
    
    require(tidymodels)
    
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
        # discretize (binarize) the outcome (immune signature)
    rec <- recipe(data.train) %>% 
        update_role(eval(outcome.var),new_role = "outcome") %>% 
        update_role(-all_outcomes(),new_role = "predictor") %>% 
        update_role(ID, new_role = "id variable") %>% 
        step_normalize(all_numeric()) %>% 
        step_dummy(all_nominal_predictors()) %>% 
        step_discretize(all_outcomes(),
                        num_breaks = 2,
                        # minimum number of samples per bin
                        min_unique = 1,
                        # `option` is necessary so that no bin_missing (empty) is created
                        options = list(keep_na = F, 
                                       # smallest to largest
                                       labels = c("outcome_lo","outcome_hi"),
                                       prefix = ""))
    
    
    # Define model specification:
        # logistic regression (classifier) for outcome low and high
        # penalty and mixture will be cross-validated over
        # engine is elastic net regression (glmnet)
    tune.spec <- logistic_reg(penalty = tune(), mixture = tune()) %>% 
        set_engine("glmnet")
    
    # Create grid for penalty and mixture values:
    lambda.grid <- grid_regular(penalty(),mixture(),levels=n.levels)
    
    # Set up folds for cross-validation:
    cv.folds <- vfold_cv(data.train, v = n.folds, repeats = n.repeats)
    
    # Define workflow (specify preprocessor and mode type):
    wf <- workflow() %>% 
        add_recipe(rec) %>% 
        add_model(tune.spec)
    
    # Cross-validate parameters over the grid
    model.grid <- tune_grid(wf,
                            resamples = cv.folds,
                            grid = lambda.grid)
    
    # Collect cross-validation metrics
    cv.metrics <- model.grid %>% 
        collect_metrics()
    
    # Choose best parameters
    # (select_best() not working as desired
    # due to poor loss resolution with accuracy metric)
    
    if (cv.metric == "accuracy"){
        best.regs <- cv.metrics %>% 
            filter(.metric == "accuracy") %>% 
            slice_max(mean) %>% 
            slice_min(mixture) %>% 
            slice_max(penalty)
    } else if (cv.metric == "roc_auc"){
        best.regs <- model.grid %>% 
            select_best(metric = "roc_auc")
    }

    # Finalize model with best penalty and mixture
    final.model <- finalize_workflow(wf, best.regs)
    
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
        # Generate/extract ROC curve:
            # Use model fit to training data
            # Predict outcome from testing data
            # Compare to truth (obtained from preprocessed data)
            # Generate ROC curve
        
        roc.curve <- predict(final.model %>%
                                 fit(data.train) %>%
                                 extract_fit_parsnip(),
                             bake(prep(rec),new_data = data.test) %>%
                                 select(-eval(outcome.var), -ID),
                             type="prob") %>%
            bind_cols(bake(prep(rec),new_data = data.test) %>% select(eval(outcome.var))) %>%
            rename(truth = eval(outcome.var)) %>%
            roc_curve(truth,.pred_outcome_lo)
        
        # Assess predictions on the testing data
        fit.metrics <- last_fit(final.model,
                                data.split,
                                metrics = metric_set(accuracy,roc_auc,spec,sens)) %>%
            collect_metrics()
        
        return(list(cv.scores = cv.metrics,
                    best.params=best.regs,
                    final.model=trained.fit,
                    fit=model.coeffs,
                    metrics=fit.metrics,
                    roc.curve = roc.curve))
        
    }

}


plot_cv_accuracy <- function(cv.scores){
    
    # plots countours of cross-validation accuracy
    # red dot flags the maximum accuracy
    
    require(ggplot2)
    cv.scores %>% 
        filter(.metric == "accuracy") %>%
        ggplot(aes(x=log10(penalty),
                   y=mixture,
                   z=mean))+
        geom_contour_filled()+
        geom_point(data=. %>% slice_max(mean), color="red")+
        theme_classic()
}



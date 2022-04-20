
tuned_glmnet_cox <- function(data.train,
                             data.test = NULL,
                             cv.metric = "deviance",
                             n.repeats = 5,
                             n.levels = c(100,25),
                             n.folds = 10){
    
    # Predicts a time-to-event outcome (specified by `data.train$time and data.train$event`)
    # w/ cross-validated logistic regression on training data (`data.train`)
    # If testing data (`data.test`) is provided, assesses models w/ C-index and deviance
    
    require(glmnet)
    require(tidymodels)
    
    # Define preprocessing steps:
        # assign roles to variables
        # normalize all numeric predictors (convert to z-score)
        # convert nominal predictors to dummy variables
    preprocess_data <- function(raw.data){
        recipe(raw.data) %>% 
            update_role(c(event,time),new_role = "outcome") %>%
            update_role(-all_outcomes(),new_role = "predictor") %>%
            update_role(ID,new_role = "id variable") %>%
            step_normalize(all_numeric_predictors()) %>%
            step_dummy(all_nominal_predictors()) %>%
            # run the above recipe on the data
            prep() %>% 
            juice() %>%
            # remove ID column
            select(-ID)
    }
    
    # Preprocess training data
    data.train <- data.train %>% preprocess_data()
    
    # Separate training outcomes and predictors
    # (creating survival object in the process)
    out.train <- Surv(event = data.train$event,
                      time = data.train$time)
    
    pred.train <- data.train %>% 
        select(-event,-time) %>% 
        as.matrix()
    
    
    if (is.null(data.test)){
        message("No testing data provided - was this intentional?")
    } else {
        # Check that training and testing data have the same variables
        if (!setequal(colnames(data.train),colnames(data.test))){
            stop("training and testing data do not have the same variables")
        }
        
        # Preprocess testing data
        data.test <- data.test %>% preprocess_data()
        
        # Reconstruct test/train split (if testing data was provided)
        data <- bind_rows(data.train, data.test)
        test.train.indices <- list(analysis = seq(nrow(data.train)),
                                   assessment = nrow(data.train) + seq(nrow(data.test)))
        
        data.split <- make_splits(test.train.indices, data)
        
        # Separate testing outcomes and predictors
        out.test <- Surv(event = data.test$event,
                         time = data.test$time)
        
        pred.test <- data.test %>% 
            select(-event,-time) %>% 
            as.matrix()
        
    }
    
    
    # Create grid for penalty and mixture values:
    cv.grid <- grid_regular(penalty(),mixture(),levels = n.levels)

    # Define fold indices: 
        # get a sequence of one to n.folds,
        # repeat enough times to cover all samples, 
        # truncate length to number of training samples,
        # randomly shuffle
    fold.id <- rep(1:n.folds,
                   times = ceiling(length(out.train)/n.folds))[1:length(out.train)] %>% 
        sample(size=length(out.train),
               replace = F)
    
    # Loop through mixture (alpha) values because cv.glmnet will not do so automatically
    # note: using cv.glmnet here because tidymodels::tune() does not yet support cox regressions
    
    cv.scores <- map_dfr(unique(cv.grid$mixture),function(alpha.cv){
        
        cv.res.tmp <- cv.glmnet(y=out.train,
                                x=pred.train,
                                family = "cox",
                                type.measure = cv.metric,
                                lambda = unique(cv.grid$penalty),
                                foldid = fold.id,
                                alpha=alpha.cv)
        
        # Store lambda, mixture, and deviance values
        tibble(penalty = cv.res.tmp$lambda,
               mixture = alpha.cv,
               .metric = cv.metric,
               mean = cv.res.tmp$cvm,
               std_err = cv.res.tmp$cvsd
        )
    })
    
    # Select best penalty and mixture
    if (cv.metric == "deviance"){
        best.lambda <- cv.res %>% slice_min(order_by = mean) %>% slice_head(n=1) %>% pull(penalty)
        best.mixture <- cv.res %>% slice_min(order_by = mean) %>% slice_head(n=1) %>% pull(mixture)
    } else if (cv.metric == "C"){
        best.lambda <- cv.res %>% slice_max(order_by = mean) %>% slice_head(n=1) %>% pull(penalty)
        best.mixture <- cv.res %>% slice_max(order_by = mean) %>% slice_head(n=1) %>% pull(mixture)
    }
    
    # Refit with best lambda and mixture
    trained.fit <- glmnet(y=out.train,
                       x=pred.train,
                       family="cox",
                       lambda = best.lambda,
                       alpha = best.mixture)
    
    # Extract model coefficients
    model.coeffs <- as.matrix(coef(trained.fit)) %>%
        as.data.frame() %>% 
        rownames_to_column(var="term") %>% 
        as_tibble() %>% 
        rename(estimate = s0) %>% 
        arrange(desc(abs(estimate)))
    
    # If there is no testing data, stop here.
    if (is.null(data.test)){
        
        return(list(fit = model.coeffs,
                    best.params = best.regs,
                    final.model = trained.fit,
                    cv.scores = cv.metrics))
        
    # Otherwise, proceed to testing
    } else if (!is.null(data.test)){
        
        # Assessing model on testing data
        fit.metrics.list <- assess.glmnet(trained.fit,
                                          newx = pred.test,
                                          newy = out.test,
                                          family = "cox")
        
        fit.metrics <- tibble(.metric = names(fit.metrics.list),
                              .estimate = sapply(fit.metrics.list,simplify = T,
                                                 function(x){x}))
        
        return(list(cv.scores = cv.metrics,
                    best.params=best.regs,
                    final.model=trained.fit,
                    fit=model.coeffs,
                    metrics=fit.metrics))
    }
}

plot_cv_deviance <- function(cv.scores){
    
    # plots contours of cross-validation MAE
    # red dot flags minimum deviance
    cv.scores %>% 
        filter(.metric == "deviance") %>%
        ggplot(aes(x=log10(penalty),
                   y=mixture,
                   z=mean))+
        geom_contour_filled()+
        geom_point(data= . %>% slice_min(mean), color="red")+
        theme_classic()
}
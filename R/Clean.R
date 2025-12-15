#######################################################################################################################################################################
#######################################################################################################################################################################
###################################################################################### Utility functions ##############################################################
#######################################################################################################################################################################
#######################################################################################################################################################################
## Dark Causality
# data The input dataset, rows are samples and columns are genes.
# cause	The number of integer positions of the cause variables in the dataset.
# effect The number of integer positions of the target variables in the dataset.
# E The embedding dimension, which influences the number of dimensions in which the time series is reconstructed for analysis.
# tau The time delay used in reconstructing the time series in the embedded space.
# metric A character string indicating the distance metric to be used (e.g., 'euclidean', 'maximum').
# h The prediction horizon, representing the number of steps ahead for which predictions are needed.
# weighted A logical indicating whether to use a weighted approach in the causality strength calculations.
# tpb A bool parameter to show progress bar.
# num.cores The numbers of CPU cores to run.
darkcausality_parallel <- function(data, cause, effect, E = 3, tau = 1, metric = "euclidean", h = 1, weighted = TRUE, tpb = FALSE, num.cores){
    
    allnames <- colnames(data)
    causenames <- allnames[cause]
    effectnames <- allnames[effect]
    index <- expand.grid(cause, effect)

    # get number of cores to run
    cl <- makeCluster(num.cores)
    registerDoParallel(cl)

    res.tmp <- foreach(i = 1:nrow(index), .packages = c("patterncausality")) %dopar% {
        patternCausality_all <- pcLightweight(data[, index[i, 1]], data[, index[i, 2]], E = E, tau = tau, metric = metric, h = h, weighted = weighted, tpb = tpb)	
    }

    # shut down the workers
    stopCluster(cl)
    stopImplicitCluster()
    
    res <- do.call(rbind, res.tmp)
    res_total <- matrix(res[, 1], ncol = length(cause), byrow = T)
    res_positive <- matrix(res[, 2], ncol = length(cause), byrow = T)
    res_negative <- matrix(res[, 3], ncol = length(cause), byrow = T)
    res_dark <- matrix(res[, 4], ncol = length(cause), byrow = T)
    colnames(res_total) <- colnames(res_positive) <- colnames(res_negative) <- colnames(res_dark) <- causenames
    rownames(res_total) <- rownames(res_positive) <- rownames(res_negative) <- rownames(res_dark) <- effectnames
    
    return(list(res_total, res_positive, res_negative, res_dark))
}

## Obtain significant positive, negative or dark causal patterns
# P_matrix Postive causality matrix.
# N_matrix Negative causality matrix.
# D_matrix Dark causality matrix.
# strength_cutoff The cutoff of causal strength.
sigPC <- function(P_matrix, N_matrix, D_matrix, strength_cutoff = 0.60){

  m <- nrow(P_matrix)
  n <- ncol(P_matrix)

  res_P <- res_N <- res_D <- matrix(0, m, n)
  rownames(res_P) <- rownames(res_N) <- rownames(res_D) <- rownames(P_matrix)
  colnames(res_P) <- colnames(res_N) <- colnames(res_D) <- colnames(P_matrix)

  for(i in seq(m)){
    for(j in seq(n)){
      if(P_matrix[i, j] > N_matrix[i, j] & P_matrix[i, j] > D_matrix[i, j] & P_matrix[i, j] > strength_cutoff){
        res_P[i, j] <- 1
      } else if(N_matrix[i, j] > P_matrix[i, j] & N_matrix[i, j] > D_matrix[i, j] & N_matrix[i, j] > strength_cutoff){
        res_N[i, j] <- 1
      } else if(D_matrix[i, j] > P_matrix[i, j] & D_matrix[i, j] > N_matrix[i, j] & D_matrix[i, j] > strength_cutoff){
        res_D[i, j] <- 1
      }
    }
  }

  return(list(res_P, res_N, res_D))
}

## Granger Causality
# data The input dataset, rows are samples and columns are genes.
# cause	The number of integer positions of the cause variables in the dataset.
# effect The number of integer positions of the target variables in the dataset.
# pvalue_cutoff p-value cutoff for rejecting Null Hypothesis (Time series X does not cause time series Y to Granger-cause itself.).
# num.cores The numbers of CPU cores to run.
granger_parallel <- function(data, cause, effect, pvalue_cutoff, num.cores){
    
    allnames <- colnames(data)
    causenames <- allnames[cause]
    effectnames <- allnames[effect]
    index <- expand.grid(cause, effect)

    # get number of cores to run
    cl <- makeCluster(num.cores)
    registerDoParallel(cl)

    res.tmp <- foreach(i = 1:nrow(index), .packages = c("lmtest", "vars", "dynlm")) %dopar% {
        x <- data[, index[i, 1]]
	y <- data[, index[i, 2]]  
	x <- ts(x)
        y <- ts(y)
	m1 <- dynlm(x ~ L(x, 1:1) + L(y, 1:1))
        m2 <- dynlm(x ~ L(x, 1:1))
	# Run on aliased variables when granger test does not
	# Functionally the same as running the granger test (uses the waldtest by default)
        granger_pvalue  <- anova(m1, m2, test = "F")$Pr[2]        
    }

    # shut down the workers
    stopCluster(cl)
    stopImplicitCluster()

    res <- matrix(unlist(res.tmp), ncol = length(cause), byrow = T)
    res <- (res < pvalue_cutoff)
    res[res] <- 1
    colnames(res) <- causenames
    rownames(res) <- effectnames
    
    return(res)
}

## Sequential Invariant Causal Prediction
# data The input dataset, rows are samples and columns are genes.
# cause	The number of integer positions of the cause variables in the dataset.
# effect The number of integer positions of the target variables in the dataset.
# pvalue_cutoff p-value cutoff for rejecting Null Hypothesis (Time series X does not cause time series Y).
# method Sequential Invariant Causal Prediction (seqICP) or Non-linear Invariant Causal Prediction (seqICPnl).
# num.cores The numbers of CPU cores to run.
seqICP_parallel <- function(data, cause, effect, pvalue_cutoff, method = c("seqICP", "seqICPnl"), num.cores){    
    
    allnames <- colnames(data)
    causenames <- allnames[cause]
    effectnames <- allnames[effect]
    index <- expand.grid(cause, effect)

    # get number of cores to run
    cl <- makeCluster(num.cores)
    registerDoParallel(cl)

    res.tmp <- foreach(i = 1:nrow(index), .packages = c("seqICP")) %dopar% {
        if(method == "seqICP"){
	    set.seed(123)
            seqICP_pvalue <- seqICP(data[, index[i, 1]], data[, index[i, 2]])$p.values
	}else{
	    set.seed(123)
	    seqICPnl_pvalue <- seqICPnl(data[, index[i, 1]], data[, index[i, 2]])$p.values
	}
    }

    # shut down the workers
    stopCluster(cl)
    stopImplicitCluster()

    res <- matrix(unlist(res.tmp), ncol = length(cause), byrow = TRUE)
    res <- (res < pvalue_cutoff)
    res[res] <- 1
    colnames(res) <- causenames
    rownames(res) <- effectnames
    
    return(res)
}

## Calculating heterogeneity matrix between two list of networks
# net1 List object, the first list of network.
# net2 List object, the second list of network.
Het.network <- function(net1, net2){

    if(class(net1)!="list" | class(net2)!="list") {
    stop("Please check your input network! The input network should be list object! \n")
    }

    m <- length(net1)
    n <- length(net2)
    Het <- matrix(NA, m, n)
    for (i in seq(m)){
        for (j in seq(n)){
	    overlap_interin <- nrow(as_data_frame(net1[[i]] %s% net2[[j]]))
	    Het[i, j] <- 1 - overlap_interin/min(nrow(as_data_frame(net1[[i]])), nrow(as_data_frame(net2[[j]])))
	}
    }
        
    return(Het)
}

## Calculating similarity matrix between two list of networks
# net1: List object, the first list of network.
# net2: List object, the second list of network.
Sim.network <- function(net1, net2){

    if(class(net1)!="list" | class(net2)!="list") {
    stop("Please check your input network! The input network should be list object! \n")
    }

    m <- length(net1)
    n <- length(net2)
    Sim <- matrix(NA, m, n)
    for (i in seq(m)){
        for (j in seq(n)){
	    overlap_interin <- ecount(net1[[i]] %s% net2[[j]])
	    Sim[i, j] <- overlap_interin/min(ecount(net1[[i]]), ecount(net2[[j]]))
	}
    }
        
    return(Sim)
}

## Identifying the overlap between multiple hubs
# hub List object, the list of hub.
# overlap.num The minimum number of each hub existing in multiple hubs.
Overlap.hub <- function(hub, overlap.num = 1, type = c("equal", "least")){

    if(class(hub)!="list") {
    stop("Please check your input hub! The input hub should be list object! \n")
    }

    hub.transform <- unlist(hub)
    hub.table <- table(hub.transform)
    
    if(type == "least"){
        hub.overlapped <- which(hub.table >= overlap.num)
    } else if(type == "equal"){
        hub.overlapped <- which(hub.table == overlap.num)
    }

    return(names(hub.overlapped))
}

# Function to identify gained and lost edges for each lncRNA
# g_asd An igraph object of ASD network.
# g_normal An igraph object of normal network.
calculate_lncrna_rewiring <- function(g_asd, g_normal) {
  
  # Get edge lists
  edges_asd <- igraph::as_data_frame(g_asd)
  edges_normal <- igraph::as_data_frame(g_normal)
  
  # Add edge identifiers
  if (nrow(edges_asd) > 0) {
    edges_asd$edge_id <- paste(edges_asd$from, edges_asd$to, sep = "_")
  } else {
    edges_asd$edge_id <- character(0)
  }
  
  if (nrow(edges_normal) > 0) {
    edges_normal$edge_id <- paste(edges_normal$from, edges_normal$to, sep = "_")
  } else {
    edges_normal$edge_id <- character(0)
  }
  
  # Identify gained edges (in ASD but not in normal)
  if (nrow(edges_asd) > 0) {
    gained_edges <- edges_asd[!edges_asd$edge_id %in% edges_normal$edge_id, ]
  } else {
    gained_edges <- data.frame()
  }
  
  # Identify lost edges (in normal but not in ASD)
  if (nrow(edges_normal) > 0) {
    lost_edges <- edges_normal[!edges_normal$edge_id %in% edges_asd$edge_id, ]
  } else {
    lost_edges <- data.frame()
  }
  
  # Get all lncRNAs (source nodes from both networks)
  # Assuming lncRNAs are the source nodes in the edges
  all_lncrnas <- unique(c(edges_asd$from, edges_normal$from))
  
  # Count gained and lost edges per lncRNA (source node only)
  if (nrow(gained_edges) > 0) {
    gained_count <- as.data.frame(table(gained_edges$from))
    colnames(gained_count) <- c("node", "gained_edges")
  } else {
    gained_count <- data.frame(node = character(), gained_edges = numeric())
  }
  
  if (nrow(lost_edges) > 0) {
    lost_count <- as.data.frame(table(lost_edges$from))
    colnames(lost_count) <- c("node", "lost_edges")
  } else {
    lost_count <- data.frame(node = character(), lost_edges = numeric())
  }
  
  # Create base dataframe with lncRNAs only
  rewiring_df <- data.frame(node = all_lncrnas)
  
  # Merge with gained and lost counts
  rewiring_df <- merge(rewiring_df, gained_count, by = "node", all.x = TRUE)
  rewiring_df <- merge(rewiring_df, lost_count, by = "node", all.x = TRUE)
  
  # Replace NA with 0
  rewiring_df[is.na(rewiring_df$gained_edges), "gained_edges"] <- 0
  rewiring_df[is.na(rewiring_df$lost_edges), "lost_edges"] <- 0
  
  # Calculate total rewiring
  rewiring_df$total_rewiring <- rewiring_df$gained_edges + rewiring_df$lost_edges  
  
  return(list(
    rewiring_df = rewiring_df,
    gained_edges = gained_edges,
    lost_edges = lost_edges
  ))
}

# Function to calculate risk scores based on z-score of rewiring metrics
# rewiring_df Rewiring edges including gained and lost edges.
calculate_risk_scores_zscore <- function(rewiring_df) {
  
  risk_scores <- rewiring_df
  
  # Calculate z-scores for each rewiring metric
  if (sd(risk_scores$gained_edges) > 0) {
    risk_scores$z_gained <- scale(risk_scores$gained_edges)[,1]
  } else {
    risk_scores$z_gained <- 0
  }
  
  if (sd(risk_scores$lost_edges) > 0) {
    risk_scores$z_lost <- scale(risk_scores$lost_edges)[,1]
  } else {
    risk_scores$z_lost <- 0
  }
  
  if (sd(risk_scores$total_rewiring) > 0) {
    risk_scores$z_total <- scale(risk_scores$total_rewiring)[,1]
  } else {
    risk_scores$z_total <- 0
  }  
  
  risk_scores$risk_score <- risk_scores$z_total  
  
  # Handle NA values
  risk_scores$risk_score[is.na(risk_scores$risk_score)] <- 0
  
  # Rank by risk score (higher z-score = higher risk)
  risk_scores <- risk_scores[order(-risk_scores$risk_score), ]
  
  return(risk_scores)
}

# Identify significant risk lncRNAs using statistical testing
# risk_scores Risk score of each lncRNA.
# p_value_threshold The threshold of p value for identifying significant risk lncRNAs.
identify_significant_risk_lncrnas <- function(risk_scores, 
                                             p_value_threshold = 0.05) {
  
  # Calculate empirical p-values based on z-score distribution
  # Assuming normal distribution of z-scores under null hypothesis
  risk_scores$p_value <- pnorm(risk_scores$risk_score, lower.tail = FALSE)  
  
  # Identify significant risk lncRNAs
  risk_scores$significant <- risk_scores$p_value < p_value_threshold
  
  return(risk_scores)
}

## SVM training and evaluation with e1071 
# data Input expression data.
# target_col The column of class.
# n_folds Number of folds for corss validation.
# parallel Logical value, TRUE for parallel, FALSE for non-parallel.
# n_cores Number of cpu cores to set.
svm_cv <- function(data, target_col = "Class", n_folds = 10, 
                           parallel = TRUE, n_cores = NULL) {  
  
  # Clean up any existing parallel sessions first
  try(stopImplicitCluster(), silent = TRUE)
  try(registerDoSEQ(), silent = TRUE)
  Sys.sleep(1)  # Brief pause to ensure cleanup
  
  # Prepare data as matrix
  features <- as.matrix(data[, -which(names(data) == target_col)])
  labels <- data[[target_col]]
  
  # Convert to numeric matrix for better performance
  storage.mode(features) <- "numeric"
  
  # Create folds
  folds <- createFolds(labels, k = n_folds)
  
  # Setup parallel processing with safer defaults
  if (parallel) {
    if (is.null(n_cores)) {
      n_cores <- min(parallel::detectCores() - 1, 48)  # Limit cores to avoid connection issues
    }
    
    # Clean any existing clusters
    try(stopCluster(cl), silent = TRUE)
    
    # Create new cluster with explicit cleanup
    cl <- makePSOCKcluster(n_cores)
    registerDoParallel(cl)
    
    # Ensure cleanup on exit
    on.exit({
      try(stopCluster(cl), silent = TRUE)
      try(registerDoSEQ(), silent = TRUE)
    })
  }
  
  cat("Performing", n_folds, "fold cross-validation...\n")
  
  # Pre-allocate results storage
  all_predictions <- character(nrow(data))
  all_true_labels <- character(nrow(data))
  all_probabilities <- matrix(0, nrow = nrow(data), ncol = length(levels(labels)))
  colnames(all_probabilities) <- levels(labels)
  
  # Process folds with better error handling
  fold_results <- foreach(fold = 1:n_folds, 
                         .packages = c("e1071", "pROC"),
                         .export = c("features", "labels", "folds"),
                         .errorhandling = "pass") %dopar% {
    
    tryCatch({
      # Get fold indices
      test_indices <- folds[[fold]]
      train_indices <- setdiff(1:nrow(data), test_indices)
      
      # Extract data subsets
      train_features <- features[train_indices, , drop = FALSE]
      train_labels <- labels[train_indices]
      test_features <- features[test_indices, , drop = FALSE]
      test_labels <- labels[test_indices]
      
      # Train SVM
      svm_model <- svm(
        x = train_features,
        y = train_labels,
        probability = TRUE,
        kernel = "radial",
        cache_size = 100,  # Reduced cache size
        scale = FALSE,     # Disable scaling for speed
        tolerance = 0.001  # Slightly higher tolerance
      )
      
      # Make predictions
      predictions <- predict(svm_model, test_features, probability = TRUE)
      prob_matrix <- attr(predictions, "probabilities")
      
      # Return compact results
      list(
        test_indices = test_indices,
        true_labels = as.character(test_labels),
        predictions = as.character(predictions),
        probabilities = prob_matrix,
        success = TRUE
      )
    }, error = function(e) {
      # Return error information but don't break
      warning("Fold ", fold, " failed: ", e$message)
      list(
        test_indices = folds[[fold]],
        true_labels = as.character(labels[folds[[fold]]]),
        predictions = rep(NA, length(folds[[fold]])),
        probabilities = matrix(NA, nrow = length(folds[[fold]]), ncol = length(levels(labels))),
        success = FALSE,
        error_message = e$message
      )
    })
  }
  
  # Clean up parallel backend immediately
  if (parallel) {
    try(stopCluster(cl), silent = TRUE)
    try(registerDoSEQ(), silent = TRUE)
  }
  
  # Assemble all results from fold_results
  successful_folds <- 0
  for (fold in 1:n_folds) {
    result <- fold_results[[fold]]
    
    if (is.list(result) && !is.null(result$success) && result$success) {
      indices <- result$test_indices
      all_predictions[indices] <- result$predictions
      all_true_labels[indices] <- result$true_labels
      
      # Ensure probability matrix alignment
      prob_cols <- colnames(result$probabilities)
      if (all(prob_cols %in% colnames(all_probabilities))) {
        all_probabilities[indices, prob_cols] <- result$probabilities[, prob_cols]
      }
      successful_folds <- successful_folds + 1
    } else {
      warning("Fold ", fold, " failed during processing")
    }
  }
  
  cat("Successfully processed", successful_folds, "out of", n_folds, "folds\n")
  
  # Convert to factors
  all_predictions <- factor(all_predictions, levels = levels(labels))
  all_true_labels <- factor(all_true_labels, levels = levels(labels))
  
  return(list(
    fold_results = fold_results,
    all_predictions = all_predictions,
    all_true_labels = all_true_labels,
    all_probabilities = all_probabilities,
    n_folds = n_folds,
    successful_folds = successful_folds
  ))
}

## Performance Evaluation
# Binary classification evaluation
# true_labels True labels of each sample.
# predictions The prediction results.
# probabilities The probabilities.
evaluate_binary <- function(true_labels, predictions, probabilities) {
  
  # Confusion matrix calculation
  cm <- table(predictions, true_labels)
  TN <- cm[1,1]
  FP <- cm[1,2]
  FN <- cm[2,1]
  TP <- cm[2,2]
  
  # Calculate metrics directly
  accuracy <- (TP + TN) / sum(cm)
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  precision <- TP / (TP + FP)
  f1_score <- 2 * (precision * sensitivity) / (precision + sensitivity)
  kappa <- (accuracy - sum(rowSums(cm) * colSums(cm)) / sum(cm)^2) / 
           (1 - sum(rowSums(cm) * colSums(cm)) / sum(cm)^2)
  
  # AUC calculation
  auc_value <- NA
  roc_obj <- NULL
  
  tryCatch({
    positive_class <- levels(true_labels)[2]
    if (positive_class %in% colnames(probabilities)) {
      prob_positive <- as.numeric(probabilities[, positive_class])
      numeric_response <- as.numeric(true_labels == positive_class)
      roc_obj <- roc(response = numeric_response, predictor = prob_positive, 
                     quiet = TRUE)  # Suppress output for speed
      auc_value <- auc(roc_obj)
    }
  }, error = function(e) {
    # Silent error handling
  })
  
  return(list(
    confusion_matrix = cm,
    accuracy = accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    f1_score = f1_score,
    kappa = kappa,
    auc = auc_value,
    roc_object = roc_obj
  ))
}

# Multiclass evaluation
# true_labels True labels of each sample.
# predictions The prediction results.
# probabilities The probabilities.
evaluate_multiclass <- function(true_labels, predictions, probabilities) {

  cm <- table(predictions, true_labels)
  
  # Calculate metrics
  n_classes <- nrow(cm)
  class_metrics <- matrix(0, nrow = n_classes, ncol = 3)
  colnames(class_metrics) <- c("Precision", "Recall", "F1")
  rownames(class_metrics) <- rownames(cm)
  
  for (i in 1:n_classes) {
    TP <- cm[i,i]
    FP <- sum(cm[i,]) - TP
    FN <- sum(cm[,i]) - TP
    TN <- sum(cm) - TP - FP - FN
    
    precision <- ifelse(TP + FP > 0, TP / (TP + FP), 0)
    recall <- ifelse(TP + FN > 0, TP / (TP + FN), 0)
    f1 <- ifelse(precision + recall > 0, 2 * precision * recall / (precision + recall), 0)
    
    class_metrics[i,] <- c(precision, recall, f1)
  }
  
  macro_precision <- mean(class_metrics[, "Precision"])
  macro_recall <- mean(class_metrics[, "Recall"])
  macro_f1 <- mean(class_metrics[, "F1"])
  accuracy <- sum(diag(cm)) / sum(cm)
  
  # Calculate kappa
  observed_accuracy <- accuracy
  expected_accuracy <- sum((rowSums(cm)/sum(cm)) * (colSums(cm)/sum(cm)))
  kappa <- ifelse(expected_accuracy < 1, 
                 (observed_accuracy - expected_accuracy) / (1 - expected_accuracy), 
                 1)
  
  # Multiclass AUC
  auc_values <- c()
  class_names <- levels(true_labels)  
  
  for (class in class_names) {
    tryCatch({
      if (class %in% colnames(probabilities)) {
        binary_labels <- as.numeric(true_labels == class)
        class_probs <- as.numeric(probabilities[, class])
        
        if (length(unique(binary_labels)) > 1) {
          roc_obj <- roc(binary_labels, class_probs, quiet = TRUE)
          auc_values <- c(auc_values, auc(roc_obj))
        }
      }
    }, error = function(e) NULL)
  }
  
  macro_auc <- if (length(auc_values) > 0) mean(auc_values) else NA
  
  return(list(
    confusion_matrix = cm,
    accuracy = accuracy,
    precision = macro_precision,
    recall = macro_recall,
    f1_score = macro_f1,
    kappa = kappa,
    auc = macro_auc,
    detailed_metrics = class_metrics
  ))
}

## Cross-Validation Analysis
# cv_results Cross-validation results.
# problem_type The type of classification, binary or multiclass.
analyze_cv_results <- function(cv_results, problem_type = "binary") {
  
  # Evaluation functions
  if (problem_type == "binary") {
    performance <- evaluate_binary(cv_results$all_true_labels, 
                                      cv_results$all_predictions, 
                                      cv_results$all_probabilities)
  } else {
    performance <- evaluate_multiclass(cv_results$all_true_labels, 
                                          cv_results$all_predictions, 
                                          cv_results$all_probabilities)
  }
  
  # Calculate fold metrics
  fold_accuracies <- sapply(cv_results$fold_results, function(fold) {
    if (is.list(fold) && !is.null(fold$success) && fold$success) {
      correct <- sum(fold$predictions == fold$true_labels)
      total <- length(fold$true_labels)
      correct / total
    } else {
      NA
    }
  })
  
  # Remove NA values from failed folds
  fold_accuracies <- fold_accuracies[!is.na(fold_accuracies)]
  
  # Calculate fold AUCs only for binary classification
  fold_aucs <- NULL
  if (problem_type == "binary") {
    fold_aucs <- sapply(cv_results$fold_results, function(fold) {
      if (is.list(fold) && !is.null(fold$success) && fold$success) {
        tryCatch({
          positive_class <- levels(cv_results$all_true_labels)[2]
          if (positive_class %in% colnames(fold$probabilities)) {
            prob_positive <- as.numeric(fold$probabilities[, positive_class])
            numeric_response <- as.numeric(fold$true_labels == positive_class)
            if (length(unique(numeric_response)) > 1) {
              roc_obj <- roc(numeric_response, prob_positive, quiet = TRUE)
              return(auc(roc_obj))
            }
          }
          return(NA)
        }, error = function(e) NA)
      } else {
        NA
      }
    })
    fold_aucs <- fold_aucs[!is.na(fold_aucs)]
  }
  
  return(list(
    overall_performance = performance,
    fold_accuracies = fold_accuracies,
    fold_aucs = fold_aucs,
    mean_accuracy = if (length(fold_accuracies) > 0) mean(fold_accuracies) else NA,
    sd_accuracy = if (length(fold_accuracies) > 0) sd(fold_accuracies) else NA,
    mean_auc = if (!is.null(fold_aucs) && length(fold_aucs) > 0) mean(fold_aucs, na.rm = TRUE) else NA,
    sd_auc = if (!is.null(fold_aucs) && length(fold_aucs) > 0) sd(fold_aucs, na.rm = TRUE) else NA
  ))
}


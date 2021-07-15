#' Pseudo-Likelihood estimate of the Ising model
#' 
#' Fit an Ising model to binary data without covariates. Computes the threshold
#' and interaction parameter using the pseudo-likelihood method. In addition, it
#' computes the corresponding score for each individual on each estimated parameter.
#' 
#' Used internally in when in type="Ising" and method="mob"
#' 
#' @param y p x n matrix or dataframe including binary observations of all n individuals on all k parameters.
#' @param x Not used yet
#' @param start Not used yet
#' @param weights Not used yet
#' @param offset Not used yet
#' @param model specifies which parameters will be returned to networktree and tested
#' for measurement invariance. "correlation" indicates all choose(p, 2) interactions
#' between nodes; "main" indicates the p threshold parameters.  
#' @param estfun Logical. Should the estimated score functions of the Ising model
#' be returned?
#' @param object Not used yet
#' @param ... Not used yet
#'
#' @export
#'
isingfit <- function(y, start = NULL, weights = NULL,
                       offset = NULL, model = c("correlation", "mean"),
                       estfun = TRUE, object = FALSE, ...) {
  
  nodevars <- y 
  
  k <- ncol(nodevars) # number of nodes
  n <- nrow(nodevars) # number of observations
  varnames <- colnames(nodevars)
  nodevars <- as.matrix(nodevars)
  
  ### put dots in a list
  dotlist <- list(...)
  
  # --- STEP 1: 
  
  fit <- tryCatch(fit_pseudoposterior(x = as.matrix(nodevars), prior_var = Inf),
                  error = function(e) e, warning = function(war) war)
  val <- colSums(nodevars) # Checking for low variablity leading to errors
  if(any(val < 2) ){
    if(inherits(fit, c("error", "warning"))) {
      
      loglik <- -Inf # Ising model can't be estimated, thus Log-Likelihood set to infinity;
      
      # Just open some empty vectors, otherwise R produces errors
      scores <- c()
      coef <- c()
      vc <- NULL
      
      warning('Ising model cannot be estimated in (sub-)model, presumably because of low variability')
    }
    if(!inherits(fit, c("error", "warning"))) {
      
      warning('Structural change test cannot be computed for (sub-)model, because of low variability')
      
      # Setup objects to avoid errors 
      loglik <- -Inf
      scores <- c()
      vc <- NULL
      
      # Get parameter estimates for submodel
      coef <- c(diag(fit$sigma), as.vector(fit$sigma)[as.vector(lower.tri(fit$sigma))])
      ynam <- if(is.null(varnames)) 1L:p else varnames
      objnames <- c(paste0("main_", ynam), combn(p, 2, function(x) paste0("cor_", x[1], "_", x[2])))
      
      id <- NULL
      if(any("main"        == model)) id <- c(id, 1:k)
      if(any("correlation" == model)) id <- c(id, 1:(k*(k-1)/2) + k)
      
      coef <- coef[id]
      objnames <- objnames[id]
      names(coef) <- objnames
    }
    
  } else {
    # Maartens Function 
    # main effect of variables on diagonal
    main <- diag(fit$sigma)
    
    # interactions on off-diagonal
    R <- fit$sigma
    inter <- as.vector(fit$sigma)[as.vector(lower.tri(fit$sigma))]
    diag(R) <- 0
    
    # --- STEP 2: Compute matrixes to make actual computation easier
    nodevars <- t(nodevars)
    
    M <- R%*%nodevars # compute the sums 
    A <- exp(main + M)
    D <- 1 + A
    E <- A/ D
    
    # --- STEP 3: Compute Log Likelihood of Ising function
    S <- main + M
    loglik <- sum(nodevars*S) - sum(log(D))
    
    # --- STEP 4: Score Functions
    
    # Score Function for main effect
    score_main <- nodevars - E # 0 up to the 3rd decimal
    
    
    # Score Function for interaction 
    score_rho <- combn(k, 2,
                       function(x) (2*(nodevars[x[1] ,] * nodevars[x[2], ]) - (E[x[1], ] * nodevars[x[2], ]) - (E[x[2], ] * nodevars[x[1], ])) )
    
    # --- STEP 5: Store objects depending which measures to consider
    
    # Set-Up
    coef <- c(main, inter)
    scores <- cbind(t(score_main), score_rho)
    ynam <- if(is.null(varnames)) 1L:p else varnames
    objnames <- c(paste0("main_", ynam), combn(k, 2, function(x) paste0("cor_", x[1], "_", x[2])))
    
    id <- NULL
    if(any("main"        == model)) id <- c(id, 1:k)
    if(any("correlation" == model)) id <- c(id, 1:(k*(k-1)/2) + k)
    
    coef <- coef[id]
    scores   <- scores[, id]
    objnames <- objnames[id]
    
    # Naming the coefficients and scores
    names(coef) <- objnames
    colnames(scores) <- objnames
    
    # --- STEP 6: Compute inverse of Hessian
    index <- t(rbind(1: (k * (k - 1)/2), combn(k, 2)))
    inverse.hessian <- invert_hessian(sigma = fit$sigma, index = index,
                                      x = t(nodevars), prior_var = Inf)
    vc <- - inverse.hessian
  }
  
  if(estfun == FALSE) {
    scores <- NULL
  }
  
  
  res <- list(coefficients = coef,
              objfun = -loglik,
              estfun = scores,
              object = vc)
  
  return(res)
}

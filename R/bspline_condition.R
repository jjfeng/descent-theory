# B-spline - lets understand the condition number

library(splines)
library(fda)
library(mgcv)

# Fits a single B-spline
Fit_bspline <- function(train_y, train_bspline_matrix, val_y, val_bspline_matrix, penmat, lam) {
    n_train <- length(train_y)
    
    fitted_gamma <- lsfit(
        1/n_train * t(train_bspline_matrix) %*% train_bspline_matrix + lam * penmat,
        t(train_bspline_matrix) %*% train_y,
        intercept = F
    )        
    coef <- 1/n_train * fitted_gamma$coefficients
    
    list(
        fitted_coef=coef,
        train=mean((train_y - train_bspline_matrix %*% coef)^2),
        validation=mean((val_y - val_bspline_matrix %*% coef)^2)
    )
}

Make_data <- function(n, ord=4, snr=2) {
    n <- floor(n)
    ## Generate X data
    x <- runif(n)
    epsilon <- rnorm(n=length(x), sd=1)
    y <- sin(x) + epsilon * sqrt(sum(sin(x)^2)) / sqrt(sum(epsilon^2)) / snr
    
    # plot(x,y)
    
    ## Split train and validation    
    train_idx <- sample(seq(1, length(x)), length(x)/2, replace=F)
    train_x <- x[train_idx]
    train_y <- y[train_idx]
    val_x <- x[-train_idx]
    val_y <- y[-train_idx]
    
    ## Create basis function matrices
    spline_xs <- sort(sample(train_x, sqrt(length(train_x)), replace = F))
    basisobj <- create.bspline.basis(
        c(-0.1,1.1),
        breaks=c(rep(0, ord), spline_xs, rep(1, ord)),
        norder=ord
    )
    dropind <- c(seq(1, ord - 1), basisobj$nbasis - seq(ord - 2, 0))
    train_bspline_matrix <- getbasismatrix(train_x, basisobj)[,-dropind]
    val_bspline_matrix <- getbasismatrix(val_x, basisobj)[,-dropind]
    penmat <- bsplinepen(basisobj, Lfdobj=2)[-dropind,-dropind]
    list(
        train_y=train_y,
        train_bspline_matrix=train_bspline_matrix,
        val_y=val_y,
        val_bspline_matrix=val_bspline_matrix,
        penmat=penmat
    )
}

Do_bspline_CV <- function(dataset, lambdas=10^seq(-13, 0)) {
    # ## Generate X data
    # x <- runif(n)
    # epsilon <- rnorm(n=length(x), sd=1)
    # y <- sin(x) + epsilon * sqrt(sum(sin(x)^2)) / sqrt(sum(epsilon^2)) / snr
    # 
    # # plot(x,y)
    # 
    # ## Split train and validation    
    # train_idx <- sample(seq(1, length(x)), length(x)/2, replace=F)
    # train_x <- x[train_idx]
    # train_y <- y[train_idx]
    # val_x <- x[-train_idx]
    # val_y <- y[-train_idx]
    # 
    # ## Create basis function matrices
    # spline_xs <- sort(sample(train_x, sqrt(length(train_x)), replace = F))
    # basisobj <- create.bspline.basis(
    #     c(-0.1,1.1),
    #     breaks=c(rep(0, ord), spline_xs, rep(1, ord)),
    #     norder=ord
    # )
    # dropind <- c(seq(1, ord - 1), basisobj$nbasis - seq(ord - 2, 0))
    # train_bspline_matrix <- getbasismatrix(train_x, basisobj)[,-dropind]
    # val_bspline_matrix <- getbasismatrix(val_x, basisobj)[,-dropind]
    # penmat <- bsplinepen(basisobj, Lfdobj=2)[-dropind,-dropind]
    
    ## Fit the model!
    if (length(lambdas) == 1 & lambdas[1] == 0) {
        res <- Fit_bspline(
            dataset$train_y,
            dataset$train_bspline_matrix,
            dataset$val_y,
            dataset$val_bspline_matrix,
            dataset$penmat,
            lam=0
        )
        fitted_coef <- res$fitted_coef
        train_loss <- res$train
        validation_loss <- res$validation
        best_lam <- 0
    } else {
        res <- lapply(lambdas, function(lam) {
            Fit_bspline(
                dataset$train_y,
                dataset$train_bspline_matrix,
                dataset$val_y,
                dataset$val_bspline_matrix,
                dataset$penmat,
                lam
            )
        })
        validation_errs <- unlist(lapply(res, "[[", "validation"))
        best_lam_idx <- which.min(validation_errs)
        # plot(lambdas, validation_errs, log="x")
        fitted_coef <- res[[best_lam_idx]]$fitted_coef
        best_lam <- lambdas[best_lam_idx]
        train_loss <- res[[best_lam_idx]]$train
        validation_loss <- res[[best_lam_idx]]$validation
    }
    
    data.frame(
        n=length(dataset$train_y) + length(dataset$val_y),
        train_loss=train_loss,
        validation_loss=validation_loss,
        best_lam=best_lam,
        fitted_values_sq=mean((dataset$val_bspline_matrix %*% fitted_coef)^2),
        penalty_sqrt=t(fitted_coef) %*% dataset$penmat %*% fitted_coef,
        condition_nbr=kappa(dataset$penmat)
    )
}

Do_bspline_CV_repl <- function(reps, n, lambdas, snr=2) {
    res <- replicate(reps, {
        dataset <- Make_data(n, snr=snr)
        Do_bspline_CV(dataset, lambdas = lambdas)
    }, simplify = T)
    res <- data.frame(t(res))
    data.frame(sapply(res, as.numeric))
}

Do_bspline_CV_repl(3, 100, lambdas=c(1, 2), snr=2)

# Think about the tail lambda behavior
cv_tail_shape <- lapply(c(1, 2, 4, 8), function(snr) {
    val_loss_diff <- replicate(reps, {
        dataset <- Make_data(n, snr=snr)
        res_small <- Do_bspline_CV(dataset, lambdas = c(1e-15))
        res_regular <- Do_bspline_CV(dataset, lambdas = c(1e-7, 1e-6, 1e-5))
        val_loss_diff <- res_small$validation_loss - res_regular$validation_loss
    }, simplify = T)
    data.frame(
        snr=snr,
        val_loss_diff=val_loss_diff
    )
})
aggregate(val_loss_diff ~ snr, do.call("rbind", cv_tail_shape), FUN = median)

# Think about the tail lambda behavior - changing SNR
cv_tail_shape <- lapply(c(1, 2, 4, 8), function(snr) {
    val_loss_diff <- replicate(20, {
        dataset <- Make_data(n = 100, snr=snr)
        res_small <- Do_bspline_CV(dataset, lambdas = c(1e-15))
        res_regular <- Do_bspline_CV(dataset, lambdas = c(1e-7, 1e-6, 1e-5))
        val_loss_diff <- res_small$validation_loss - res_regular$validation_loss
    }, simplify = T)
    data.frame(
        snr=snr,
        val_loss_diff=val_loss_diff
    )
})
aggregate(val_loss_diff ~ snr, do.call("rbind", cv_tail_shape), FUN = median)

# Think about the tail lambda behavior - changing sample size
cv_tail_shape_v_n <- lapply(2 * 10^seq(1, 4, 0.5), function(n) {
    val_loss_diff <- replicate(10, {
        dataset <- Make_data(n, snr=2)
        res_small <- Do_bspline_CV(dataset, lambdas = c(0))
        res_regular <- Do_bspline_CV(dataset, lambdas = c(1e-8, 1e-7, 1e-6))
        val_loss_diff <- res_small$validation_loss - res_regular$validation_loss
    }, simplify = T)
    data.frame(
        n=n,
        val_loss_diff=val_loss_diff
    )
})
aggregate(val_loss_diff ~ n, do.call("rbind", cv_tail_shape_v_n), FUN = median)

########## Lambda is always zero
reps <- 50
bspline_lam0 <- lapply(3 * 10^seq(1, 3, 0.25), function(n) {
    n <- floor(n)
    res <- Do_bspline_CV_repl(reps, n, lambdas=c(0))
    data.frame(
        n=n,
        best_lam=mean(res$best_lam),
        pen_to_val=mean(res$penalty_sqrt/res$fitted_values_sq)
    )
})
bspline_lam0 <- do.call("rbind", bspline_lam0)

# Condition number is growing in n for sure.
plot(
    bspline_lam0$n,
    bspline_lam0$pen_to_val,
    log="xy",
    xlab="Number of samples",
    ylab="penalty:norm",
    main="Penalty to norm ratio for sobolev, choose lambda = 0"
)

############ Lambda selected by CV
reps <- 50
cv_lambdas <- 10^seq(-15, -3)
bspline_CVs <- lapply(3 * 10^seq(1, 3, 0.25), function(n) {
    n <- floor(n)
    res <- Do_bspline_CV_repl(reps, n, lambdas=cv_lambdas)
    print(n)
    data.frame(
        n=n,
        best_lam=10^mean(log10(res$best_lam)),
        lam_std=10^sqrt(var(log10(res$best_lam))/resp),
        num_small_lam=sum(res$best_lam <= cv_lambdas[3]),
        pen_to_val=mean(res$penalty_sqrt/res$fitted_values_sq)
    )
})
bspline_CVs <- do.call("rbind", bspline_CVs)

# Ratio of penalty with fitted values does not tend to grow with n
plot(
    bspline_CVs$n,
    bspline_CVs$pen_to_val,
    log="xy",
    xlab="Number of samples",
    ylab="penalty:norm",
    main="Penalty to norm ratio for sobolev, choose lambda via Train/Val split"
)
plot(
    bspline_CVs$n,
    bspline_CVs$best_lam,
    log="y",
    xlab="Number of samples",
    ylab="Best lambda",
    main="Best Lambda for sobolev chosen via Train/Val split"
)
# small lambdas are chosen a lot with smaller sample sizes
# small lambdas are rarely chosen with big sample sizes
plot(
    bspline_CVs$n,
    bspline_CVs$num_small_lam,
    xlab="Number of samples",
    ylab="Num small lambdas",
    main="Num Small Lambdas for sobolev chosen via Train/Val split"
)























# B-spline

library(splines)
library(mgcv)
library(fda)

# Fits a single B-spline
Fit_bspline <- function(train_y, train_bspline_matrix, penmat, lam) {
    n_train <- length(train_y)
    fitted_gamma <- solve(
        1/n_train * t(train_bspline_matrix) %*% train_bspline_matrix + lam * penmat,
        t(train_bspline_matrix) %*% train_y
    )
    coef <- 1/n_train * fitted_gamma
}

Get_bspline_validation_loss <- function(dataset, lam) {
    coef <- Fit_bspline(dataset$train_y, dataset$train_bspline_matrix, dataset$penmat, lam)
    
    list(
        fitted_coef=coef,
        train=mean((dataset$train_y - dataset$train_bspline_matrix %*% coef)^2),
        validation=mean((dataset$val_y - dataset$val_bspline_matrix %*% coef)^2),
        true_validation_loss=mean((dataset$val_true_y - dataset$val_bspline_matrix %*% coef)^2)
    )
}

Get_bspline_oracle_loss <- function(dataset, lam) {
    coef <- Fit_bspline(dataset$train_y, dataset$train_bspline_matrix, dataset$penmat, lam)
    
    list(
        fitted_coef=coef,
        train=mean((dataset$train_y - dataset$train_bspline_matrix %*% coef)^2),
        test=mean((dataset$test_y - dataset$test_bspline_matrix %*% coef)^2),
        true_validation_loss=mean((dataset$val_true_y - dataset$val_bspline_matrix %*% coef)^2)
    )
}

# currently truth is the sin function
Make_data <- function(n_train, n_validate, n_test, ord=4, snr=2, xmin=0, xmax=6) {
    ## Generate X data
    x <- runif(n_train + n_validate + n_test, min=xmin, max=xmax)
    epsilon <- rnorm(n=length(x), sd=1)
    true_y <- sin(x)
    y <- true_y + epsilon * sqrt(sum(sin(x)^2)) / sqrt(sum(epsilon^2)) / snr
    
    ## Split train and validation
    shuffled_idx <- sample(seq(1, length(x)), length(x), replace=F)
    train_idx <- shuffled_idx[seq(n_train)]
    val_idx <- shuffled_idx[seq(n_train, n_train + n_validate)]
    test_idx <- shuffled_idx[seq(n_train + n_validate, n_train + n_validate + n_test)]
    train_x <- x[train_idx]
    train_y <- y[train_idx]
    val_x <- x[val_idx]
    val_y <- y[val_idx]
    val_true_y <- true_y[val_idx]
    test_x <- x[test_idx]
    test_y <- y[test_idx]
    
    ## Create basis function matrices
    spline_xs <- sort(sample(train_x, sqrt(length(train_x)), replace = F))
    basisobj <- create.bspline.basis(
        c(-0.1 + xmin, 0.1 + xmax),
        breaks=c(rep(xmin, ord), spline_xs, rep(xmax, ord)),
        norder=ord
    )
    dropind <- c(seq(1, ord - 1), basisobj$nbasis - seq(ord - 2, 0))
    train_bspline_matrix <- getbasismatrix(train_x, basisobj)[,-dropind]
    val_bspline_matrix <- getbasismatrix(val_x, basisobj)[,-dropind]
    test_bspline_matrix <- getbasismatrix(test_x, basisobj)[,-dropind]
    penmat <- bsplinepen(basisobj, Lfdobj=2)[-dropind,-dropind]
    list(
        train_y=train_y,
        train_bspline_matrix=train_bspline_matrix,
        val_y=val_y,
        val_true_y=val_true_y,
        val_bspline_matrix=val_bspline_matrix,
        test_y=test_y,
        test_bspline_matrix=test_bspline_matrix,
        penmat=penmat
    )
}

Do_bspline_CV <- function(dataset, lambdas=10^seq(-13, 0)) {
    ## Fit the model!
    res <- lapply(lambdas, function(lam) {
        Get_bspline_validation_loss(dataset, lam)
    })
    validation_errs <- unlist(lapply(res, "[[", "validation"))
    best_lam_idx <- which.min(validation_errs)
    # plot(lambdas, validation_errs, log="x")
    fitted_coef <- res[[best_lam_idx]]$fitted_coef
    best_lam <- lambdas[best_lam_idx]
    train_loss <- res[[best_lam_idx]]$train
    validation_loss <- res[[best_lam_idx]]$validation
    true_validation_loss <- res[[best_lam_idx]]$true_validation_loss
    
    data.frame(
        # n=length(dataset$train_y) + length(dataset$val_y),
        # cv_train_loss=train_loss,
        # cv_validation_loss=validation_loss,
        cv_true_validation_loss=true_validation_loss,
        cv_lam=best_lam
        # fitted_values_sq=mean((dataset$val_bspline_matrix %*% fitted_coef)^2),
        # penalty_sqrt=t(fitted_coef) %*% dataset$penmat %*% fitted_coef
    )
}

Do_bspline_oracle <- function(dataset, lambdas=10^seq(-13, 0)) {
    ## Fit the model!
    res <- lapply(lambdas, function(lam) {
        Get_bspline_oracle_loss(dataset, lam)
    })
    test_errs <- unlist(lapply(res, "[[", "test"))
    best_lam_idx <- which.min(test_errs)
    # plot(lambdas, test_errs, log="x")
    fitted_coef <- res[[best_lam_idx]]$fitted_coef
    best_lam <- lambdas[best_lam_idx]
    train_loss <- res[[best_lam_idx]]$train
    test_loss <- res[[best_lam_idx]]$test
    true_validation_loss <- res[[best_lam_idx]]$true_validation_loss
    
    data.frame(
        # n=length(dataset$train_y) + length(dataset$val_y),
        # oracle_train_loss=train_loss,
        # oracle_test_loss=test_loss,
        oracle_true_validation_loss=true_validation_loss,
        oracle_lam=best_lam
        # fitted_values_sq=mean((dataset$val_bspline_matrix %*% fitted_coef)^2),
        # penalty_sqrt=t(fitted_coef) %*% dataset$penmat %*% fitted_coef
    )
}

Do_bspline_cv_oracle_repl <- function(reps, n_train, n_validate, n_test, lambdas, snr=2) {
    res <- replicate(reps, {
        dataset <- Make_data(n_train, n_validate, n_test, snr=snr)
        cv_res <- Do_bspline_CV(dataset, lambdas = lambdas)
        oracle_res <- Do_bspline_oracle(dataset, lambdas = lambdas)
        res <- data.frame(cv_res, oracle_res)
        res$loss_diff <- res$cv_true_validation_loss - res$oracle_true_validation_loss
        res
    }, simplify = T)
    res <- data.frame(t(res))
    data.frame(sapply(res, as.numeric))
}

set.seed(10)

lambdas <- 10^seq(from=-7, to=-2, by=0.1)
n_sizes <- seq(from=10, to=100, by=10)
n_reps <- 50

cv_to_oracle_compare <- lapply(n_sizes, function(n) {
    print(n)
    cv_oracle <- Do_bspline_cv_oracle_repl(
        reps=n_reps, 
        n_train=n,
        n_validate=n,
        n_test=n * 20,
        lambdas=lambdas
    )
    data.frame(
        n=n,
        true_cv_loss=mean(cv_oracle$cv_true_validation_loss),
        true_oracle_loss=mean(cv_oracle$oracle_true_validation_loss),
        mean_diff=mean(cv_oracle$loss_diff),
        oracle_lambda=mean(cv_oracle$oracle_lam),
        cv_lambda=mean(cv_oracle$cv_lam)
    )
})
cv_to_oracle_compare <- do.call("rbind", cv_to_oracle_compare)

plot(
    cv_to_oracle_compare$n, cv_to_oracle_compare$true_cv_loss, type = "l", col="red",
    ylim = c(min(cv_to_oracle_compare$true_cv_loss, cv_to_oracle_compare$true_oracle_loss), max(cv_to_oracle_compare$true_cv_loss, cv_to_oracle_compare$true_oracle_loss))
)
lines(cv_to_oracle_compare$n, cv_to_oracle_compare$true_oracle_loss, col="green")

plot(cv_to_oracle_compare$n, 1/(cv_to_oracle_compare$n), ylim=c(0, 0.1), type = "l")
lines(cv_to_oracle_compare$n, cv_to_oracle_compare$mean_diff)

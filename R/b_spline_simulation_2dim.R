# B-spline, 2dim
library(doParallel)
library(foreach)
library(splines)
library(mgcv)
library(fda)

# n_cores <- 3
# cluster <- makeCluster(n_cores)
# registerDoParallel(cluster)

w <- 0.0001

Get_norm2 <- function(h) { sqrt(mean(h^2)) }
# Fits a single B-spline
Fit_bspline <- function(train_y, train_bspline_matrix, penmats, lams) {
    n_train <- length(train_y)
    gamma_size <- dim(train_bspline_matrix)[2]/length(lams)
    lam_penmats <- lapply(seq(length(penmats)), function(i) {
        lams[i] * penmats[[i]]
    })
    BB <- t(train_bspline_matrix) %*% train_bspline_matrix
    
    # also calculate the additional ridge perturbation
    lam_BBs <- lapply(seq(length(lams)), function(i) {
        Bi <- train_bspline_matrix[, seq((i-1) * gamma_size + 1, i * gamma_size)]
        lams[i] * t(Bi) %*% Bi
    })
    BB_diag <- bdiag(lam_BBs)
    
    fitted_gamma <- solve(
        1/n_train * BB + bdiag(lam_penmats) + w * BB_diag,
        t(train_bspline_matrix) %*% train_y
    )
    coef <- 1/n_train * fitted_gamma
}

Get_bspline_losses <- function(dataset, lams) {
    coef <- Fit_bspline(dataset$train_y, dataset$train_bspline_matrix, dataset$penmats, lams)
    # plot(dataset$X[[1]][dataset$val_idx], dataset$val_bspline_matrix %*% coef, col="green")  
    # points(dataset$X[[1]][dataset$val_idx], dataset$val_true_y)
    list(
        fitted_coef=coef,
        train=Get_norm2(dataset$train_y - dataset$train_bspline_matrix %*% coef),
        validation=Get_norm2(dataset$val_y - dataset$val_bspline_matrix %*% coef),
        test=Get_norm2(dataset$test_y - dataset$test_bspline_matrix %*% coef),
        true_validation_loss=Get_norm2(dataset$val_true_y - dataset$val_bspline_matrix %*% coef)
    )
}

Make_bspline_matrix <- function(X, idx, basis_objs, ord) {
    bspline_matrices <- mapply(function(x, basisobj) {
        dropind <- c(seq(1, ord - 1), basisobj$nbasis - seq(ord - 2, 0))
        getbasismatrix(x[idx], basisobj)[,-dropind]
    }, X, basis_objs, SIMPLIFY = F)
    do.call("cbind", bspline_matrices)
}

Make_bspline_objects <- function(X, train_idx, ord, xmin, xmax) {
    bspline_objs <- lapply(X, function(x) {
        spline_xs <- sort(sample(x[train_idx], sqrt(length(train_idx)), replace = F))
        create.bspline.basis(
            c(-0.1 + xmin, 0.1 + xmax),
            breaks=c(rep(xmin, ord), spline_xs, rep(xmax, ord)),
            norder=ord
        )    
    })
    bspline_objs
}

Make_bspline_pen <- function(bspline_objs, ord) {
    penmats <- lapply(bspline_objs, function(basisobj) {
        dropind <- c(seq(1, ord - 1), basisobj$nbasis - seq(ord - 2, 0))
        bsplinepen(basisobj, Lfdobj=2)[-dropind,-dropind]
    })
    penmats
}


# f1 = sin function, f2 = sin(4x + 1) function
# sqrt n_train knots
Make_data <- function(n_train, n_validate, n_test, ord=4, snr=2, xmin=0, xmax=6) {
    ## Generate X data
    n <- n_train + n_validate + n_test
    x1 <- runif(n, min=xmin, max=xmax)
    x2 <- runif(n, min=xmin, max=xmax)
    X <- list(x1, x2)
    # epsilon <- rnorm(n=n, sd=1)
    epsilon <- runif(n=n, min=-1, max=1)
    true_y <- sin(x1) + sin(x2 * 4 + 1)
    y <- true_y + epsilon * Get_norm2(true_y) / Get_norm2(epsilon) / snr
    
    ## Split train and validation
    shuffled_idx <- sample(seq(1, n), n, replace=F)
    train_idx <- shuffled_idx[seq(n_train)]
    val_idx <- shuffled_idx[seq(n_train, n_train + n_validate)]
    test_idx <- shuffled_idx[seq(n_train + n_validate, n_train + n_validate + n_test)]
    train_y <- y[train_idx]
    val_y <- y[val_idx]
    val_true_y <- true_y[val_idx]
    test_y <- y[test_idx]
    
    ## Create basis function matrices
    bspline_objs <- Make_bspline_objects(X, train_idx, ord, xmin, xmax)
    train_bspline_matrix <- Make_bspline_matrix(X, train_idx, bspline_objs, ord)
    val_bspline_matrix <- Make_bspline_matrix(X, val_idx, bspline_objs, ord)
    test_bspline_matrix <- Make_bspline_matrix(X, test_idx, bspline_objs, ord)
    penmats <- Make_bspline_pen(bspline_objs, ord)
    list(
        X=X,
        train_idx=train_idx,
        val_idx=val_idx,
        test_idx=test_idx,
        train_y=train_y,
        train_bspline_matrix=train_bspline_matrix,
        val_y=val_y,
        val_true_y=val_true_y,
        val_bspline_matrix=val_bspline_matrix,
        test_y=test_y,
        test_bspline_matrix=test_bspline_matrix,
        penmats=penmats
    )
}

Do_bspline_CV_oracle <- function(dataset, lambdas1, lambdas2) {
    ## Fit the model!
    lambda_grid <- expand.grid(lambdas1, lambdas2)
    res <- apply(lambda_grid, 1, function(lams) {
        Get_bspline_losses(dataset, lams)
    })
    validation_errs <- unlist(lapply(res, "[[", "validation"))
    cv_lam_idx <- which.min(validation_errs)
    cv_best_lams <- lambda_grid[cv_lam_idx,]
    cv_true_validation_loss <- res[[cv_lam_idx]]$true_validation_loss
    cv_validation_loss <- res[[cv_lam_idx]]$validation
    
    test_errs <- unlist(lapply(res, "[[", "test"))
    oracle_lam_idx <- which.min(test_errs)
    oracle_best_lams <- lambda_grid[oracle_lam_idx,]
    oracle_true_validation_loss <- res[[oracle_lam_idx]]$true_validation_loss
    oracle_validation_loss <- res[[oracle_lam_idx]]$validation
    
    print(paste("oracle_best_lams", oracle_best_lams))
    print(paste("cv_best_lams", cv_best_lams))
    if (cv_lam_idx == 1 || cv_lam_idx == length(lambda_grid)) {
        print(paste("cv_lam_idx", cv_lam_idx))
    }
    if (oracle_lam_idx == 1 || oracle_lam_idx == length(lambda_grid)) {
        print(paste("oracle_lam_idx", oracle_lam_idx))
    }
    
    
    data.frame(
        eps=mean((dataset$val_y - dataset$val_true_y)^2),
        cv_true_validation_loss=cv_true_validation_loss,
        cv_lams=cv_best_lams,
        cv_validation_loss=cv_validation_loss,
        oracle_true_validation_loss=oracle_true_validation_loss,
        oracle_validation_loss=oracle_validation_loss,
        oracle_lams=oracle_best_lams,
        loss_diff=cv_true_validation_loss - oracle_true_validation_loss
    )
}

Do_bspline_cv_oracle_repl <- function(reps, n_train, n_validate, n_test, lambdas1, lambdas2, snr=2) {
    res <- replicate(reps, {
        dataset <- Make_data(n_train, n_validate, n_test, snr=snr)
        Do_bspline_CV_oracle(dataset, lambdas1, lambdas2)
    }, simplify = T)
    # res <- foreach(
    #     rep=seq(reps),
    #     .combine= rbind,
    #     .export=ls(envir=globalenv())
    # ) %dopar% {
    #     dataset <- Make_data(n_train, n_validate, n_test, snr=snr)
    #     Do_bspline_CV_oracle(dataset, lambdas = lambdas)
    # }
    res <- data.frame(t(res))
    data.frame(sapply(res, as.numeric))
}

set.seed(10)


lambdas1 <- 10^seq(from=-5, to=-0.5, by=0.2)
lambdas2 <- 10^seq(from=-9, to=-3, by=0.2)
n_sizes <- seq(from=5, to=30, by=5)
n_reps <- 25


## Important: To see useful trends in the empirical process term,
# we should vary only the number of validation samples
# and keep all other terms constant!
cv_to_oracle_compare_w <- lapply(n_sizes, function(n) {
    print(n)
    cv_oracle <- Do_bspline_cv_oracle_repl(
        reps=n_reps, 
        n_train=30,
        n_validate=n,
        n_test=400,
        lambdas1=lambdas1,
        lambdas2=lambdas2
    )
    print(colMeans(cv_oracle))
    data.frame(
        n=n,
        t(colMeans(cv_oracle))
    )
})
cv_to_oracle_compare_w <- do.call("rbind", cv_to_oracle_compare_w)

ylim <- c(min(
        cv_to_oracle_compare_w$oracle_true_validation_loss,
        cv_to_oracle_compare_w$cv_true_validation_loss
    ), max(
        cv_to_oracle_compare_w$oracle_true_validation_loss,
        cv_to_oracle_compare_w$cv_true_validation_loss
))

plot(
    cv_to_oracle_compare_w$n, cv_to_oracle_compare_w$cv_true_validation_loss, type = "l", col="red",
    ylim = c(0.4, 1.0), #ylim
    ylab="Validation Loss",
    xlab="Validation Set Size"
)
lines(cv_to_oracle_compare_w$n, cv_to_oracle_compare_w$oracle_true_validation_loss, col="green")
legend(20,1,c("Training/Validation Split", "Oracle"),lty=c(1,1), lwd=c(2.5,2.5), col=c("red","green"))

plot(
    (log(cv_to_oracle_compare_w$n)/cv_to_oracle_compare_w$n)^0.5,
    abs(cv_to_oracle_compare_w$loss_diff),
    xlab="Expected Convergence Rate",
    ylab="Empirical Validation Loss Difference"
)
regl <- lm(I(abs(loss_diff)) ~ I((log(n)/n)^.5), cv_to_oracle_compare_w)
summary(regl)
abline(regl)

# B-spline, 2dim
library(doParallel)
library(foreach)
library(splines)
library(mgcv)
library(fda)

n_cores <- 3

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
        1/n_train * BB + bdiag(lam_penmats),
        t(train_bspline_matrix) %*% train_y
    )
    coef <- 1/n_train * fitted_gamma
}

Get_bspline_losses <- function(dataset, lams) {
    coef <- Fit_bspline(dataset$train_y, dataset$train_bspline_matrix, dataset$penmats, lams)
    list(
        coef=coef,
        train=Get_norm2(dataset$train_y - dataset$train_bspline_matrix %*% coef),
        validation=Get_norm2(dataset$val_y - dataset$val_bspline_matrix %*% coef),
        true_validation=Get_norm2(dataset$val_true_y - dataset$val_bspline_matrix %*% coef)
    )
}

Make_bspline_matrix <- function(X, idx, basis_objs, ord) {
    bspline_matrices <- mapply(function(x, basisobj) {
        dropind <- c(seq(1, ord - 1), basisobj$nbasis - seq(ord - 2, 0))
        getbasismatrix(x[idx], basisobj)[,-dropind]
    }, X, basis_objs, SIMPLIFY = F)
    cols <- lapply(bspline_matrices, function(m){ncol(m)})
    bspline_matrices <- do.call("cbind", bspline_matrices)
}

# choose sqrt n_train for spline knots
Make_bspline_objects <- function(X, train_idx, ord, xmin, xmax) {
    bspline_objs <- lapply(X, function(x) {
        sorted_X <- sort(x)
        sorted_trainX <- sort(x[train_idx])
        n <- length(train_idx)
        num_knots <- floor(sqrt(n))
        spline_xs <- sorted_trainX[seq(1, n, by=floor(n/num_knots))]
        spl <- create.bspline.basis(
            c(-0.1 + xmin, 0.1 + xmax),
            breaks=c(rep(sorted_X[1], ord), spline_xs, rep(sorted_X[length(x)], ord)),
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


# sqrt n_train knots
Make_data <- function(n_train, n_validate, ord=4, snr=2, xmin=-1, xmax=1) {
    ## Generate X data
    n <- n_train + n_validate
    x1 <- runif(n, min=xmin, max=xmax)
    x2 <- runif(n, min=xmin, max=xmax)
    
    X <- list(x1, x2)
    
    epsilon <- rnorm(n=n, sd=1)
    true_y <- exp(x1) + (x2)^2
    y <- true_y + epsilon * Get_norm2(true_y) / Get_norm2(epsilon) / snr
    
    ## Split train and validation
    train_idx <- seq(n_train)
    val_idx <- seq(n_train, n_train + n_validate)
    train_y <- y[train_idx]
    val_y <- y[val_idx]
    val_true_y <- true_y[val_idx]
    
    ## Create basis function matrices
    bspline_objs <- Make_bspline_objects(X, train_idx, ord, xmin, xmax)
    train_bspline_matrix <- Make_bspline_matrix(X, train_idx, bspline_objs, ord)
    val_bspline_matrix <- Make_bspline_matrix(X, val_idx, bspline_objs, ord)
    penmats <- Make_bspline_pen(bspline_objs, ord)
    list(
        X=X,
        train_idx=train_idx,
        val_idx=val_idx,
        train_y=train_y,
        train_bspline_matrix=train_bspline_matrix,
        val_y=val_y,
        val_true_y=val_true_y,
        val_bspline_matrix=val_bspline_matrix,
        penmats=penmats
    )
}

Get_best_lambdas <- function(res, lambda_grid, col) {
    errs <- unlist(lapply(res, "[[", col))
    lam_idx <- which.min(errs)
    list(
        lambdas=lambda_grid[lam_idx,],
        lam_idx=lam_idx,
        true_validation=res[[lam_idx]]$true_validation,
        validation=res[[lam_idx]]$validation
    )
    
}

Make_lambda_grid <- function(range1, range2, interval) {
    lambdas1 <- 10^seq(range1[1], range1[2], by=interval)
    lambdas2 <- 10^seq(range2[1], range2[2], by=interval)
    list(
        lambdas1=lambdas1,
        lambdas2=lambdas2,
        grid=expand.grid(lambdas1, lambdas2)
    )
}

Eval_losses <- function(dataset, grid) {
    apply(grid, 1, function(lams) {
        Get_bspline_losses(dataset, lams)
    })
}

Do_bspline_CV_oracle <- function(dataset, lambda1_range, lambda2_range, grid_int) {
    ## Fit the model!
    lam_grid <- Make_lambda_grid(
        lambda1_range,
        lambda2_range,
        grid_int
    )
    grid_res <- Eval_losses(dataset, lam_grid$grid)
    
    cv_best_lams <- Get_best_lambdas(grid_res, lam_grid$grid, "validation")
    print(paste("cv_best_lams", cv_best_lams$lambdas))
    
    oracle_best_lams <- Get_best_lambdas(grid_res, lam_grid$grid, "true_validation")
    print(paste("oracle_best_lams", oracle_best_lams$lambdas))
    
    # Plot_bpline_fits(dataset, dataset$train_idx, dataset$train_bspline_matrix, res[[oracle_lam_idx]], res[[cv_lam_idx]])
    # Plot_bpline_fits(dataset, dataset$val_idx, dataset$val_bspline_matrix, grid_res[[oracle_best_lams$lam_idx]], grid_res[[cv_best_lams$lam_idx]])
    
    data.frame(
        eps=mean((dataset$val_y - dataset$val_true_y)^2),
        cv_true_validation=cv_best_lams$true_validation,
        cv_lams=cv_best_lams$lambdas,
        cv_validation_loss=cv_best_lams$validation,
        oracle_true_validation_loss=oracle_best_lams$true_validation,
        oracle_validation_loss=oracle_best_lams$validation,
        oracle_lams=oracle_best_lams$lambdas,
        loss_diff_sq=cv_best_lams$true_validation^2 - oracle_best_lams$true_validation^2
    )
}


Plot_bpline_fits <- function(dataset, observations_idx, dataset_bspline_matrix, oracle_res, cv_res) {
    num_knots <- ncol(dataset_bspline_matrix)/2
    par(mfrow=c(1,2))
    X1 <- dataset$X[[1]][observations_idx]
    plot(X1, exp(X1), ylim = c(-5, 5)) # truth
    points(X1, dataset_bspline_matrix[,seq(num_knots)] %*% cv_res$coef[seq(num_knots)], col="red")
    points(X1, dataset_bspline_matrix[,seq(num_knots)] %*% oracle_res$coef[seq(num_knots)], col="green")
    
    X2 <- dataset$X[[2]][observations_idx]
    plot(X2, (X2)^2, ylim = c(-4, 5)) # truth
    points(X2, dataset_bspline_matrix[,seq(num_knots + 1, num_knots*2)] %*% cv_res$coef[seq(num_knots + 1, num_knots*2)], col="red")
    points(X2, dataset_bspline_matrix[,seq(num_knots + 1, num_knots*2)] %*% oracle_res$coef[seq(num_knots + 1, num_knots*2)], col="green")
}

Do_bspline_cv_oracle_repl <- function(reps, n_train, n_validate, lambda1_range, lambda2_range, grid_int, snr=2) {
    datasets <- lapply(seq(reps), function(i){
        Make_data(n_train, n_validate, snr=snr)
    })
    
    cl <- makeCluster(n_cores)
    clusterExport(cl, 
                  c("lambda1_range", "lambda2_range", "grid_int", "datasets", "bdiag",
                      "Do_bspline_CV_oracle", "Get_bspline_losses", "Make_lambda_grid", "Eval_losses", "Get_best_lambdas", "Fit_bspline", "Get_norm2")
    )
    res <- parLapply(cl, seq(reps), function(i){
        dataset <- datasets[[i]]
        Do_bspline_CV_oracle(dataset, lambda1_range, lambda2_range, grid_int)
    })
    res <- do.call("rbind", res)
    stopCluster(cl)
    
    data.frame(sapply(res, as.numeric))
}

set.seed(1)

n_trains <- 100
n_test <- 10
lambda1_range <- c(-9, -2)
lambda2_range <- c(-9, -2)
grid_int <- 0.05
n_sizes <- floor(20 * 2^seq(7, 0, by=-0.2))
n_reps <- 20
snr <- 2

cv_to_oracle_all <- lapply(n_sizes, function(n_val) {
    cv_to_oracle_ntrains <- lapply(n_trains, function(n_train){
        print(paste("train", n_train, "val", n_val))
        cv_oracle <- Do_bspline_cv_oracle_repl(
            reps=n_reps, 
            n_train=n_train,
            n_validate=n_val,
            lambda1_range=lambda1_range, 
            lambda2_range=lambda2_range, 
            grid_int=grid_int, 
            snr=snr
        )
        print(colMeans(cv_oracle))
        data.frame(
            n=n_val,
            n_train=n_train,
            cv_oracle
        )
    })
    do.call("rbind", cv_to_oracle_ntrains)
})
cv_to_oracle_all <- do.call("rbind", cv_to_oracle_all)
save(cv_to_oracle_all, file = "cv_to_oracle_all_final.RData")
cv_to_oracle_not_tiny <- cv_to_oracle_all[cv_to_oracle_all$loss_diff_sq > 1e-10,]
cv_to_oracle_compare_w <- aggregate(loss_diff_sq ~ n, cv_to_oracle_all, FUN = mean)

pdf('figures/validation_size_loss_diff_final.pdf', width=10, height=6)
par(mar=c(5,5,1,1), mfrow=c(1,1))
plot(
    cv_to_oracle_not_tiny$n,
    cv_to_oracle_not_tiny$loss_diff_sq,
    ylab="Validation Loss Difference",
    xlab="Validation Set Size",
    cex=0.3,
    cex.axis=1.25,
    cex.lab=1.25,
    log="xy"
)
lm_data <- cv_to_oracle_not_tiny
fit <- lm(log(loss_diff_sq) ~ log(n), lm_data)
fit_summ <- summary(fit)
lines(
    lm_data$n,
    exp(fit$fitted.values)
)
legend(
    "bottomleft", 
    legend=paste(
        "Exponent:", 
        format(fit$coefficients[2], digits=4), 
        "+/-", 
        format(1.96 * fit_summ$coefficients[2,2], digits=4)
    )
)
dev.off()


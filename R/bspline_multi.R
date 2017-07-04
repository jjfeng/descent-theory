# B-spline, 2dim
library(doParallel)
library(foreach)
library(splines)
library(mgcv)
library(fda)
library(ggplot2)

jitter_log <- function(vals, scaler=0.05) {
    noise <- rnorm(length(vals), mean=0, sd=vals*scaler)
    vals + noise
}


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

Get_bspline_val_losses <- function(dataset, lams) {
    coef <- Fit_bspline(dataset$train_y, dataset$train_bspline_matrix, dataset$penmats, lams)
    Get_norm2(dataset$val_y - dataset$val_bspline_matrix %*% coef)
}

Get_bspline_true_val_losses <- function(dataset, lams) {
    coef <- Fit_bspline(dataset$train_y, dataset$train_bspline_matrix, dataset$penmats, lams)
    Get_norm2(dataset$val_true_y - dataset$val_bspline_matrix %*% coef)
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
Make_data <- function(n_train, n_validate, num_p, true_func, ord=4, snr=2, xmin=-2, xmax=2) {
    ## Generate X data
    n <- n_train + n_validate
    
    X <- lapply(seq(num_p), function(x){runif(n, min=xmin, max=xmax)})
    
    epsilon <- rnorm(n=n, sd=1)
    true_y <- true_func(X)
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

Eval_losses <- function(dataset, grid) {
    apply(grid, 1, function(lams) {
        Get_bspline_losses(dataset, lams)
    })
}

Create_full_lams <- function(lambda_minis, n_lams, num_p) {
    lam_logs <- rep(lambda_minis, each=num_p/n_lams)
    10^lam_logs
}

Search_lam <- function(dataset, n_lams, loss_fun, additional_init=NA) {
    optim_fn <- function(lam_logs_mini) {
        lams <- Create_full_lams(lam_logs_mini, n_lams, num_p)
        loss_val <- loss_fun(dataset, lams)
        if (is.na(loss_val)) {
            1000
        } else {
            loss_val
        }
    }
    init_lam_vals = c(0, -1, -2)
    optim_cv_lam <- rep(init_lam_vals[1], n_lams)
    best_cv_res <- 1000
    for (init_val in init_lam_vals) {
        init_lams = rep(init_val, n_lams)
        cv_res <- nlm(
            optim_fn,
            init_lams
        )
        if (cv_res$minimum < best_cv_res) {
            best_cv_res <- cv_res$minimum
            optim_cv_lam <- cv_res$estimate
            print(init_val)
            print(optim_cv_lam)
        }
    }
    
    if (!is.na(additional_init)) {
        cv_res <- nlm(
            optim_fn,
            additional_init
        )
        if (cv_res$minimum < best_cv_res) {
            best_cv_res <- cv_res$minimum
            optim_cv_lam <- cv_res$estimate
        }
    }
    
    print(optim_cv_lam)
    optim_cv_lam
}

Do_bspline_CV_oracle <- function(dataset, n_lams, plt_func) {
    ## Fit the model!
    num_p = length(dataset$X)
    
    optim_cv_lam <- Search_lam(dataset, n_lams, Get_bspline_val_losses)
    print(optim_cv_lam)
    optim_oracle_lam <- Search_lam(dataset, n_lams, Get_bspline_true_val_losses, additional_init = optim_cv_lam)
    print(optim_oracle_lam)
    
    cv_lams <- Create_full_lams(optim_cv_lam, n_lams, num_p)
    oracle_lams <- Create_full_lams(optim_oracle_lam, n_lams, num_p)
    cv_fit_res <- Get_bspline_losses(dataset, cv_lams)
    oracle_fit_res <- Get_bspline_losses(dataset, oracle_lams)
    
    plt_func(dataset, dataset$train_idx, dataset$train_bspline_matrix, cv_fit_res, oracle_fit_res)
    # Plot_bpline_fits(dataset, dataset$val_idx, dataset$val_bspline_matrix, grid_res[[oracle_best_lams$lam_idx]], grid_res[[cv_best_lams$lam_idx]])
    
    data.frame(
        eps=mean((dataset$val_y - dataset$val_true_y)^2),
        cv_true_validation=cv_fit_res$true_validation,
        cv_validation_loss=cv_fit_res$validation,
        oracle_true_validation_loss=oracle_fit_res$true_validation,
        oracle_validation_loss=oracle_fit_res$validation,
        loss_diff_sq=cv_fit_res$true_validation^2 - oracle_fit_res$true_validation^2
    )
}

Do_bspline_cv_oracle_repl <- function(reps, n_train, n_validate, num_p, n_lams, true_func, plot_fnc, snr=2) {
    datasets <- lapply(seq(reps), function(i){
        Make_data(n_train, n_validate, num_p, true_func, snr=snr)
    })
    
    res <- lapply(datasets, function(dataset){
        Do_bspline_CV_oracle(dataset, n_lams, plot_fnc)
    })
    res <- do.call("rbind", res)
    
    data.frame(sapply(res, as.numeric))
}


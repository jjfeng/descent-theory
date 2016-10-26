# B-spline, 2dim
library(doParallel)
library(foreach)
library(splines)
library(mgcv)
library(fda)

w <- 0

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
    list(
        coef=coef,
        train=Get_norm2(dataset$train_y - dataset$train_bspline_matrix %*% coef),
        validation=Get_norm2(dataset$val_y - dataset$val_bspline_matrix %*% coef),
        true_test_loss=Get_norm2(dataset$test_true_y - dataset$test_bspline_matrix %*% coef),
        true_validation_loss=Get_norm2(dataset$val_true_y - dataset$val_bspline_matrix %*% coef)
    )
}

Make_bspline_matrix <- function(X, idx, basis_objs, ord) {
    bspline_matrices <- mapply(function(x, basisobj) {
        dropind <- c(seq(1, ord - 1), basisobj$nbasis - seq(ord - 2, 0))
        getbasismatrix(x[idx], basisobj)[,-dropind]
        
        #eval.basis(x[idx], basisobj)
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
        ) # plot(spl)
    })
    bspline_objs
}

Make_bspline_pen <- function(bspline_objs, ord) {
    penmats <- lapply(bspline_objs, function(basisobj) {
        dropind <- c(seq(1, ord - 1), basisobj$nbasis - seq(ord - 2, 0))
        bsplinepen(basisobj, Lfdobj=2)[-dropind,-dropind]
        
        #bsplinepen(basisobj, Lfdobj=2)
    })
    penmats
}


# sqrt n_train knots
Make_data <- function(n_train, n_validate, n_test, ord=4, snr=2, xmin=0, xmax=1) {
    ## Generate X data
    n <- n_train + n_validate + n_test
    x1 <- runif(n, min=xmin, max=xmax)
    x1_order <- order(x1)
    x2 <- runif(n, min=xmin, max=xmax)
    
    x1 <- x1[x1_order]
    x2 <- x2[x1_order]
    X <- list(x1, x2)
    
    epsilon <- rnorm(n=n, sd=1)
    # true_y <- sin(x1) + 0.5 * sin(x2 * 2 + 1)
    true_y <- exp(x1) #+ sin(x2 * 5)
    y <- true_y + epsilon * Get_norm2(true_y) / Get_norm2(epsilon) / snr
    
    ## Split train and validation
    shuffled_idx <- sample(seq(1, n), n, replace=F)
    train_idx <- sort(shuffled_idx[seq(n_train)])
    val_idx <- sort(shuffled_idx[seq(n_train, n_train + n_validate)])
    test_idx <- sort(shuffled_idx[seq(n_train + n_validate, n_train + n_validate + n_test)])
    train_y <- y[train_idx]
    val_y <- y[val_idx]
    val_true_y <- true_y[val_idx]
    test_y <- y[test_idx]
    test_true_y <- true_y[test_idx]
    
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
        test_true_y=test_true_y,
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
    cv_true_test_loss <- res[[cv_lam_idx]]$true_test_loss
    
    # test_errs <- unlist(lapply(res, "[[", "true_validation_loss"))
    test_errs <- unlist(lapply(res, "[[", "true_test_loss"))
    oracle_lam_idx <- which.min(test_errs)
    oracle_best_lams <- lambda_grid[oracle_lam_idx,]
    oracle_true_validation_loss <- res[[oracle_lam_idx]]$true_validation_loss
    oracle_validation_loss <- res[[oracle_lam_idx]]$validation
    oracle_true_test_loss <- res[[oracle_lam_idx]]$true_test_loss
    
    print(paste("oracle_best_lams", oracle_best_lams))
    print(paste("cv_best_lams", cv_best_lams))
    if (cv_lam_idx == 1 || cv_lam_idx == length(lambda_grid)) {
        print(paste("cv_lam_idx", cv_lam_idx))
    }
    if (oracle_lam_idx == 1 || oracle_lam_idx == length(lambda_grid)) {
        print(paste("oracle_lam_idx", oracle_lam_idx))
    }
    print(paste("loss diff", cv_true_validation_loss - oracle_true_validation_loss))
    
    # Plot_bpline_fits(dataset, dataset$train_idx, dataset$train_bspline_matrix, res[[oracle_lam_idx]], res[[cv_lam_idx]])
    Plot_bpline_fits(dataset, dataset$val_idx, dataset$val_bspline_matrix, res[[oracle_lam_idx]], res[[cv_lam_idx]])

    data.frame(
        eps=mean((dataset$val_y - dataset$val_true_y)^2),
        cv_true_validation_loss=cv_true_validation_loss,
        cv_lams=cv_best_lams,
        cv_true_test_loss=cv_true_test_loss,
        cv_validation_loss=cv_validation_loss,
        oracle_true_validation_loss=oracle_true_validation_loss,
        oracle_validation_loss=oracle_validation_loss,
        oracle_true_test_loss=oracle_true_test_loss,
        oracle_lams=oracle_best_lams,
        loss_diff=cv_true_validation_loss - oracle_true_validation_loss
    )
}


Plot_bpline_fits <- function(dataset, observations_idx, dataset_bspline_matrix, oracle_res, cv_res) {
    num_knots <- ncol(dataset_bspline_matrix)/2
    
    X1 <- dataset$X[[1]][observations_idx]
    plot(X1, exp(X1), ylim = c(0, 3)) # truth
    points(X1, dataset_bspline_matrix[,seq(num_knots)] %*% cv_res$coef[seq(num_knots)], col="red")
    points(X1, dataset_bspline_matrix[,seq(num_knots)] %*% oracle_res$coef[seq(num_knots)], col="green")
    
    X2 <- dataset$X[[2]][observations_idx]
    plot(X2, sin(X2 * 5), ylim = c(-3, 2)) # truth
    points(X2, dataset_bspline_matrix[,seq(num_knots + 1, num_knots*2)] %*% cv_res$coef[seq(num_knots + 1, num_knots*2)], col="red")
    points(X2, dataset_bspline_matrix[,seq(num_knots + 1, num_knots*2)] %*% oracle_res$coef[seq(num_knots + 1, num_knots*2)], col="green")
}

Do_bspline_cv_oracle_repl <- function(reps, n_train, n_validate, n_test, lambdas1, lambdas2, snr=2) {
    res <- replicate(reps, {
        dataset <- Make_data(n_train, n_validate, n_test, snr=snr)
        Do_bspline_CV_oracle(dataset, lambdas1, lambdas2)
    }, simplify = T)
    res <- data.frame(t(res))
    data.frame(sapply(res, as.numeric))
}

set.seed(10)

# very good settings!
n_trains <- c(500) #25, 50, 75, 100) #800?
n_test <- 10000 #3000
lambdas1 <- c(10^-11, 10^-9, 10^-6, 10^-3) #10^seq(from=-7, to=-3, by=0.2)
lambdas2 <- c(10^-11, 10^-9, 10^-6, 10^-3) # 10^seq(from=-7, to=-3, by=0.2)
n_sizes <- c(100, 500) # floor(seq(0.15, 0.45, by=0.02)^-4)
n_reps <- 2
snr <- 8

## Important: To see useful trends in the empirical process term,
# we should vary only the number of validation samples
# and keep all other terms constant!
cv_to_oracle_all <- lapply(n_sizes, function(n_val) {
    cv_to_oracle_ntrains <- lapply(n_trains, function(n_train){
        print(paste("train", n_train, "val", n_val))
        cv_oracle <- Do_bspline_cv_oracle_repl(
            reps=n_reps, 
            n_train=n_train,
            n_validate=n_val,
            n_test=n_test,
            lambdas1=lambdas1,
            lambdas2=lambdas2,
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
# save(cv_to_oracle_all_pos, file = "cv_to_oracle_all_pos.RData")
cv_to_oracle_compare_w <- aggregate(loss_diff ~ n, cv_to_oracle_all, FUN = mean)

pdf('figures/validation_size_loss_diff.pdf', width=10, height=6)
par(mar=c(5,5,1,1))
plot(
    cv_to_oracle_compare_w$n,
    cv_to_oracle_compare_w$loss_diff,
    type = "b",
    # ylim = c(0.1, 0.2),
    ylab="Validation Loss Diff",
    xlab="Validation Set Size",
    # xaxt="n",
    # log="xy",
    cex.axis=1.25,
    cex.lab=1.25
)
# axis(1, at = cv_to_oracle_compare_w$n, las=2)
dev.off()

# cv_to_oracle_all_pos <- rbind(cv_to_oracle_all, do.call("rbind", cv_to_oracle_all0))
# cv_to_oracle_all_pos <- cv_to_oracle_all_pos[cv_to_oracle_all_pos$loss_diff > 0,]

cv_to_oracle_some <- cv_to_oracle_all
cv_to_oracle_compare_w_pos <- aggregate(loss_diff ~ n, cv_to_oracle_some, FUN = mean)
cv_to_oracle_some$other <- cv_to_oracle_some$loss_diff * (cv_to_oracle_some$n)^(1/4)
cv_to_oracle_some <- merge(cv_to_oracle_compare_w_pos, cv_to_oracle_some, by="n")
cv_to_oracle_some <- cv_to_oracle_some[cv_to_oracle_some$loss_diff.y > 0,]

negligible_fit <- lm(log(loss_diff.y) ~ log(n), cv_to_oracle_some)
summary(negligible_fit)

negligible_fit_big_errs <- lm(log(loss_diff.y) ~ log(n), cv_to_oracle_some[cv_to_oracle_some$loss_diff.x < cv_to_oracle_some$loss_diff.y,])
summary(negligible_fit_big_errs)

nv_negligible_fit <- lm(other ~ I(n^(-1/4)) + oracle_true_validation_loss + 0, cv_to_oracle_some)
summary(nv_negligible_fit)

nv_negligible_fit_big_errs <- lm(other ~ I(n^(-1/4))+ oracle_true_validation_loss + 0, cv_to_oracle_some[cv_to_oracle_some$loss_diff.x < cv_to_oracle_some$loss_diff.y,])
summary(nv_negligible_fit_big_errs)

cv_to_oracle_compare_max <- aggregate(loss_diff ~ n, cv_to_oracle_all, FUN = max)
cv_to_oracle_compare_max <- merge(cv_to_oracle_all, cv_to_oracle_compare_max, by=c("loss_diff", "n"))
cv_to_oracle_compare_max$other <- cv_to_oracle_compare_max$n^(0.25) * cv_to_oracle_compare_max$loss_diff

nv_negligible_fit_max <- lm(other ~ I(n^(-1/4)) + oracle_true_validation_loss + 0, cv_to_oracle_compare_max)
summary(nv_negligible_fit_max)

summary(lm(log(loss_diff) ~ log(n), cv_to_oracle_compare_max))

plot((cv_to_oracle_some$n), cv_to_oracle_some$loss_diff.y)

plot((cv_to_oracle_some$n)^(-1/4), cv_to_oracle_some$other)

plot((cv_to_oracle_compare_max$n)^(-1/4), cv_to_oracle_compare_max$other)
lines(cv_to_oracle_compare_max$n^(-0.25), nv_negligible_fit_max$fitted.values)

cv_to_oracle_all

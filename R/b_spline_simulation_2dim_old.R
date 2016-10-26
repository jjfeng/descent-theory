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
    }, X, basis_objs, SIMPLIFY = F)
    cols <- lapply(bspline_matrices, function(m){ncol(m)})
    bspline_matrices <- do.call("cbind", bspline_matrices)
}

Make_bspline_objects <- function(X, train_idx, ord, xmin, xmax) {
    # knots = sqrt(train size)
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
    true_y <- exp(x1) + sin(x2 * 5)
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
    
    test_errs <- unlist(lapply(res, "[[", "true_validation_loss"))
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
    
    train_X1 <- dataset$X[[1]][observations_idx]
    plot(train_X1, exp(train_X1), ylim = c(0, 3)) # truth
    points(train_X1, dataset_bspline_matrix[,seq(num_knots)] %*% cv_res$coef[seq(num_knots)], col="red")
    points(train_X1, dataset_bspline_matrix[,seq(num_knots)] %*% oracle_res$coef[seq(num_knots)], col="green")
    
    train_X2 <- dataset$X[[2]][observations_idx]
    plot(train_X2, sin(train_X2 * 5), ylim = c(-3, 2)) # truth
    points(train_X2, dataset_bspline_matrix[,seq(num_knots + 1, num_knots*2)] %*% cv_res$coef[seq(num_knots + 1, num_knots*2)], col="red")
    points(train_X2, dataset_bspline_matrix[,seq(num_knots + 1, num_knots*2)] %*% oracle_res$coef[seq(num_knots + 1, num_knots*2)], col="green")
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
n_trains <- c(25, 50, 100)
n_test <- 800
lambdas1 <- 10^seq(from=-7, to=-3, by=0.1)
lambdas2 <- 10^seq(from=-7, to=-3, by=0.1)
n_sizes <- c(20, 40, 80, 160)
n_reps <- 15
snr <- 2

# trying more reps and bigger sizes
# n_train <- 100
# n_test <- 800
# lambdas1 <- 10^seq(from=-10, to=-1, by=0.25)
# lambdas2 <- 10^seq(from=-10, to=-1, by=0.25)
# n_sizes <- floor(10 * 1.25^seq(0, 18))
# n_reps <- 40
# snr <- 2

# jean test
# n_train <- 100
# n_test <- 800
# lambdas1 <- 10^seq(from=-9, to=-1, by=0.2)
# lambdas2 <- 10^seq(from=-9, to=-1, by=0.2)
# n_sizes <- seq(20, 400, by=20)
# n_reps <- 60
# snr <- 4

## Important: To see useful trends in the empirical process term,
# we should vary only the number of validation samples
# and keep all other terms constant!
cv_to_oracle_all <- lapply(n_sizes, function(n) {
    print(n)
    cv_oracle <- Do_bspline_cv_oracle_repl(
        reps=n_reps, 
        n_train=n_train,
        n_validate=n,
        n_test=n_test,
        lambdas1=lambdas1,
        lambdas2=lambdas2,
        snr=snr
    )
    print(colMeans(cv_oracle))
    data.frame(
        n=n,
        cv_oracle
    )
})
cv_to_oracle_all <- do.call("rbind", cv_to_oracle_all)
# save(cv_to_oracle_all, file = "cv_to_oracle_all.RData")
cv_to_oracle_compare_w <- aggregate(loss_diff ~ n, cv_to_oracle_all, FUN = mean)

# ylim <- c(min(
#         cv_to_oracle_compare_w$oracle_true_validation_loss,
#         cv_to_oracle_compare_w$cv_true_validation_loss
#     ), max(
#         cv_to_oracle_compare_w$oracle_true_validation_loss,
#         cv_to_oracle_compare_w$cv_true_validation_loss
# ))
# 
# # pdf('figures/validation_size_loss.pdf', width=7, height=5)
# plot(
#     cv_to_oracle_compare_w$n, cv_to_oracle_compare_w$cv_true_validation_loss, type = "b", col="red",
#     ylim = c(0.1, 0.2),
#     ylab="Validation Loss",
#     xlab="Validation Set Size"
# )
# lines(
#     cv_to_oracle_compare_w$n,
#     # rep(mean(cv_to_oracle_compare_w$oracle_true_validation_loss), length(n_sizes)),
#     cv_to_oracle_compare_w$oracle_true_validation_loss,
#     col="green"
# )
# legend(50,0.2,c("Training/Validation Split", "Oracle"),lty=c(1,1), lwd=c(2.5,2.5), col=c("red","green"))
# dev.off()

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

pdf('figures/validation_size_loss_diff_poster.pdf', width=18, height=7)
par(mar=c(6,6,1,1), oma=c(0,0,0,0))
plot(
    cv_to_oracle_compare_w$n,
    cv_to_oracle_compare_w$loss_diff,
    type = "l",
    lwd=3,
    # ylim = c(0.1, 0.2),
    ylab="Validation Loss Diff",
    xlab="Validation Set Size",
    # log="xy",
    xaxt="n",
    cex.axis=2,
    cex.lab=3
)
axis(1, at = cv_to_oracle_compare_w$n, las=0, cex.axis=2)
dev.off()

# pdf('figures/qqplot.pdf', width=5, height=5)
# plot(
#     expected_loss_diff,
#     abs(cv_to_oracle_compare_w$loss_diff),
#     xlab="Expected Convergence Rate",
#     ylab="Empirical Validation Loss Difference"
# )
# dev.off()

cv_to_oracle_all_pos <- cv_to_oracle_all[cv_to_oracle_all$loss_diff > 0,]
cv_to_oracle_compare_w_pos <- aggregate(cbind(loss_diff, oracle_true_validation_loss) ~ n, cv_to_oracle_all_pos, FUN = mean)

# cv_to_oracle_all_pos$expected_loss_diff1 <- v_rate
# cv_to_oracle_all_pos$expected_loss_diff2 <- sqrt(v_rate) * sqrt(cv_to_oracle_all_pos$oracle_true_validation_loss)

# cv_to_oracle_some <- cv_to_oracle_all_pos[cv_to_oracle_all_pos$loss_diff > 0 & cv_to_oracle_all_pos$n > 100,]
# cv_to_oracle_some <- cv_to_oracle_all_pos[cv_to_oracle_all_pos$loss_diff > 0,]
cv_to_oracle_some <- cv_to_oracle_all_pos
cv_to_oracle_some$other <- cv_to_oracle_some$loss_diff * (cv_to_oracle_some$n)^(1/4)

cv_to_oracle_some <- merge(cv_to_oracle_compare_w_pos, cv_to_oracle_some, by="n")
# should i compare only the big loss diffs?
cv_to_oracle_some <- cv_to_oracle_some[cv_to_oracle_some$loss_diff.x < cv_to_oracle_some$loss_diff.y,]

negligible_fit <- lm(log(loss_diff.y) ~ log(n), cv_to_oracle_some)
summary(negligible_fit)

# negligible_fit <- lm(log(loss_diff) ~ log(n), cv_to_oracle_some)
# summary(negligible_fit)

nv_negligible_fit <- lm(log(other) ~ log(n), cv_to_oracle_some)
summary(nv_negligible_fit)

cv_to_oracle_compare_w_pos$other <- cv_to_oracle_compare_w_pos$loss_diff * cv_to_oracle_compare_w_pos$n^(1/4)
nv_negligible_fit <- lm(log(other) ~ log(n), cv_to_oracle_compare_w_pos)
summary(nv_negligible_fit)

# b <- lm((other) ~ (expected_diff1) + 0, cv_to_oracle_some)
# summary(b)
# b1 <- lm((other) ~ (expected_diff1), cv_to_oracle_some)
# summary(b1)
# anova(b, b1)
# a <- lm((other) ~ (expected_diff2), cv_to_oracle_some)
# summary(a)
# a2 <- lm((other) ~ (expected_diff2) + (expected_diff1) + 0, cv_to_oracle_some)
# summary(a2)
# anova(b, a2)
# 
# b <- lm(log(other) ~ log(expected_diff1) + 0, cv_to_oracle_some)
# summary(b)
# b1 <- lm(log(other) ~ log(expected_diff1), cv_to_oracle_some)
# summary(b1)
# anova(b, b1)
# a <- lm(log(other) ~ log(expected_diff2) , cv_to_oracle_some)
# summary(a)
# a2 <- lm(log(other) ~ log(expected_diff2) + log(expected_diff1) + 0, cv_to_oracle_some)
# summary(a2)
# a2 <- lm(log(other) ~ log(expected_diff2) + log(expected_diff1), cv_to_oracle_some)
# summary(a2)
# anova(b, a2)
# 
# cor(cv_to_oracle_all_pos$expected_diff1, cv_to_oracle_all_pos$expected_diff2)
# plot(cv_to_oracle_all_pos$n, cv_to_oracle_all_pos$expected_diff2)
# 
# mean(log(cv_to_oracle_some$loss_diff, base = cv_to_oracle_some$n))
# 
# # regl1 <- lm(loss_diff ~ I((n)^-0.25) + I((n)^-0.5), cv_to_oracle_compare_w)
# # summary(regl1)
# # 
# # regl2 <- lm(loss_diff ~ I((n)^-0.5), cv_to_oracle_compare_w)
# # summary(regl2)
# # 
# # anova(regl2, regl1)
# # 
# # plot(
# #     (1/cv_to_oracle_compare_w$n)^0.5,
# #     abs(cv_to_oracle_compare_w$loss_diff),
# #     xlab="Expected Convergence Rate (1/n-1/2)",
# #     ylab="Empirical Validation Loss Difference"
# # )

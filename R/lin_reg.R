Make_linreg_dataset <- function(n, beta, eps_factor) {
    p <- length(beta)
    
    X <- runif(n * p)
    X <- matrix(X, nrow = n)
    epsilon <- rnorm(n)
    true_y <- X %*% beta
    y <- true_y + epsilon * eps_factor
    
    ## Split train and validation
    list(
        X=X,
        true_y=true_y,
        y=y
    )
}


Get_beta_hat <- function(dataset, lambda) {
    train_X <- dataset$X
    train_y <- dataset$y
    p <- ncol(train_X)
    n <- length(train_y)
    1/n * solve(1/n * t(train_X) %*% train_X + lambda * diag(p), t(train_X) %*% train_y)
}

Get_losses <- function(beta_fitted, dataset) {
    list(
        train_loss = Get_norm2(dataset$train_y - dataset$train_X %*% beta_fitted),
        train_true_loss = Get_norm2(dataset$train_true_y - dataset$train_X %*% beta_fitted),
        val_loss = Get_norm2(dataset$val_y - dataset$val_X %*% beta_fitted),
        val_true_loss = Get_norm2(dataset$val_true_y - dataset$val_X %*% beta_fitted),
        test_loss = Get_norm2(dataset$test_y - dataset$test_X %*% beta_fitted),
        test_true_loss = Get_norm2(dataset$test_true_y - dataset$test_X %*% beta_fitted),
        beta_diff = Get_norm2(beta_fitted - dataset$beta)
    )
}

Get_losses_for_lam_idx <- function(res, lam_idx) {
    list(
        lambda=res$lambda[lam_idx],
        train_loss=res$train_loss[lam_idx],
        train_true_loss=res$train_true_loss[lam_idx],
        val_true_loss=res$val_true_loss[lam_idx],
        val_loss=res$val_loss[lam_idx],
        test_true_loss=res$test_true_loss[lam_idx]
    )
}

Do_linreg_CV <- function(dataset, lambdas, beta_hats) {
    res <- lapply(seq(length(lambdas)), function(i) {
        lambda <- lambdas[i]
        beta_hat <- beta_hats[[i]]
        loss_res <- Get_losses(beta_hat, dataset)
        data.frame(
            lambda=as.numeric(lambda),
            loss_res
        )
    })
    res <- do.call("rbind", res)
    
    cv_lam_idx <- which.min(res$val_loss)
    cv_results <- Get_losses_for_lam_idx(res, cv_lam_idx)
    
    oracle_lam_idx <- which.min(res$val_true_loss)
    oracle_results <- Get_losses_for_lam_idx(res, oracle_lam_idx)
    
    train_oracle_lam_idx <- which.min(res$train_true_loss)
    train_oracle_results <- Get_losses_for_lam_idx(res, train_oracle_lam_idx)
    
    data.frame(
        n_train=dataset$n_train,
        n_val=dataset$n_val,
        cv=cv_results,
        oracle=oracle_results,
        train_oracle=train_oracle_results,
        loss_diff=cv_results$val_true_loss - oracle_results$val_true_loss,
        loss_sq_diff=cv_results$val_true_loss^2 - oracle_results$val_true_loss^2,
        lambda_diff=abs(cv_results$lambda - oracle_results$lambda)
    )
}

Do_val_sims <- function(n_val, train_dataset, test_dataset, beta, lambdas, beta_hats, eps_factor, num_reps) {
    val_sim_res <- lapply(seq(num_reps), function(i){
        val_dataset <- Make_linreg_dataset(n_val, beta, eps_factor = eps_factor)
        dataset <- list(
            beta=beta,
            n_train=length(train_dataset$y),
            n_val=length(val_dataset$y),
            train_X=train_dataset$X,
            val_X=val_dataset$X,
            test_X=test_dataset$X,
            train_y=train_dataset$y,
            train_true_y=train_dataset$true_y,
            val_y=val_dataset$y,
            val_true_y=val_dataset$true_y,
            test_y=test_dataset$y,
            test_true_y=test_dataset$true_y
        )
        
        r <- Do_linreg_CV(dataset, lambdas, beta_hats)
        print(r)
    })
    val_sim_res <- do.call("rbind", val_sim_res)
}

Do_train_val_res <- function(n_train, n_vals, beta, eps_factor, num_reps) {
    # do simulations given training set size and validation set sizes
    # draw training set and validation set new every single time
    
    # currently test dataset is unused.
    test_dataset <- Make_linreg_dataset(n_test, beta, eps_factor = eps_factor)
    
    n_vals_rep <- rep(n_vals, num_reps)
    vals_res <- lapply(n_vals_rep, function(n_val) {
        train_dataset <- Make_linreg_dataset(n_train, beta, eps_factor = eps_factor)
        print(paste("n_val", n_val))
        lambdas <- (n_train)^-1 * 10^seq(-2, 1, by=1/n_val)
        print(paste("range(lambdas)", range(lambdas)))
        beta_hats <- lapply(lambdas, function(lambda) {
            Get_beta_hat(train_dataset, lambda)
        })
        Do_val_sims(n_val, train_dataset, test_dataset, beta, lambdas, beta_hats, eps_factor, num_reps = 1)    
    })
    vals_res <- do.call("rbind", vals_res)
}

set.seed(10)

# train_vals_res1, p 100
# lambdas <- (n_train)^-1 * 10^seq(-3.5, 1.5, by=2/n_val)

n_trains <- floor(20 * 2^seq(1, 4, by=0.5))
n_vals <- floor(20 * 2^seq(0, 8, by=0.5))
n_test <- 4
snr <- 2
p <- 100
non_zero_beta <- 30
eps_factor <- 0.25 #non_zero_beta/snr/2
num_val_reps <- 5 #20

# try bigger p?
n_trains <- 160 #floor(20 * 2^seq(3, 4, by=0.5))
n_vals <- floor(20 * 2^seq(4, 7, by=0.5))
n_test <- 4
snr <- 2
p <- 200
non_zero_beta <- 30
eps_factor <- 0.25 #non_zero_beta/snr/2
num_val_reps <- 5 #20

beta <- c(rep(1, non_zero_beta), rep(0, p - non_zero_beta))

train_vals_res2 <- lapply(n_trains, function(n_train) {
    print(paste("n_train", n_train))
    res <- Do_train_val_res(n_train, n_vals, beta, eps_factor, num_val_reps)
})
train_vals_res2 <- do.call("rbind", train_vals_res2)
# train_vals_res <- rbind(train_vals_res1, train_vals_res2)
summary(lm(
    log(loss_sq_diff) ~ log(n_val) + log(oracle.val_true_loss),
    train_vals_res2[
        train_vals_res2$loss_sq_diff > 0,
    ]
))

agg_train_vals_res <- aggregate(loss_sq_diff ~ n_val + n_train, train_vals_res1, FUN=mean)
mm <- merge(train_vals_res1, agg_train_vals_res, by=c("n_val", "n_train"))

mean(train_vals_res1$loss_diff == 0)

lm_fit <- lm(log(loss_diff) ~ log(n_val) + log(n_train), train_vals_res1[train_vals_res1$loss_diff > 0, ])
summary(lm_fit)

summary(lm(log(loss_sq_diff) ~ log(n_val) + log(oracle.val_true_loss), train_vals_res1[train_vals_res1$loss_sq_diff > 0 &  train_vals_res1$n_train > 300, ]))


lm_sq_fit <- lm(log(loss_sq_diff) ~ log(n_val) + log(n_train), train_vals_res1[train_vals_res1$loss_sq_diff > 0, ])
summary(lm_sq_fit)

lm_sq_fit_big <- lm(log(loss_sq_diff.x) ~ log(n_val) + log(n_train), mm[mm$loss_sq_diff.x > mm$loss_sq_diff.y, ])
summary(lm_sq_fit_big)

summary(lm(log(lambda_diff) ~ log(n_val) + log(n_train), train_vals_res1[train_vals_res1$lambda_diff > 0, ]))

# train_lambda oracle rate should be -1
summary(lm(log(train_oracle.lambda) ~ log(n_train), train_vals_res1))
plot(train_vals_res1$n_train, train_vals_res1$oracle.lambda)

# train_oracle.train_true_loss rate should be -0.5 - good!
summary(lm(log(train_oracle.train_true_loss) ~ log(n_train), train_vals_res1))
plot(train_vals_res1$n_train, train_vals_res1$train_oracle.train_true_loss)

###

# n_vals <- floor(400 * 2^seq(1, 5.5, by=0.5))
# n_vals <- floor(2 * 2^seq(1, 5.5, by=0.5))
# 
# vals_res <- lapply(n_vals, function(n_val) {
#     print(paste("n_val", n_val))
#     lambdas <- 10^seq(0, 3, by=4/n_val)
#     beta_hats <- lapply(lambdas, function(lambda) {
#         Get_beta_hat(train_dataset, lambda)
#     })
#     Do_val_sims(n_val, train_dataset, test_dataset, beta, lambdas, beta_hats, eps_factor)    
# })
# vals_res <- do.call("rbind", vals_res)
# # vals_res_new <- vals_res
# # vals_res_all <- rbind(vals_res_old, vals_res)
# # vals_res <- vals_res_all
# 
# plot(vals_res$n_val, vals_res$loss_diff)
# print(vals_res)
# 
# summary(lm((loss_diff) ~ I((n_val)^-0.5), vals_res))
# summary(lm((loss_diff) ~ I((n_val)^-0.25), vals_res))
# summary(lm((loss_diff * I((n_val)^0.25)) ~ I((n_val)^-0.25), vals_res))
# 
# # i think our approx that loss_diff <= sq(loss_sq_diff) is very very loose
# # so we're getting wonky results
# lm_fit <- lm(log(loss_diff) ~ log(n_val), vals_res[vals_res$loss_diff > 0, ])
# summary(lm_fit)
# 
# # as expected, it is about n^-1
# lm_sq_fit <- lm(log(loss_sq_diff) ~ log(n_val), vals_res[vals_res$loss_diff > 0, ])
# summary(lm_sq_fit)
# 
# vals_res_max <- aggregate(loss_diff ~ n_val, vals_res, FUN = max )
# lm_fit_max <- lm(log(loss_diff) ~ log(n_val), vals_res_max)
# summary(lm_fit_max)
# 
# plot(log(vals_res$n_val[vals_res$loss_diff > 0]), log(vals_res$loss_diff[vals_res$loss_diff > 0]))
# lines(log(vals_res$n_val[vals_res$loss_diff > 0]), lm_fit$fitted.values)
# #plot(lm(log(loss_diff) ~ log(n_val), vals_res[vals_res$loss_diff > 0, ]))
# 
# # see if lambdas are convering super super fast - yes they are.
# summary(l <- lm(
#     log(lambda_diff) ~ log(n_val), 
#     data = vals_res[vals_res$lambda_diff > 0,]
# ))
# plot(log(vals_res$n_val), log(abs(vals_res$cv.lambda - vals_res$oracle.lambda)))
# lines(log(vals_res$n_val[abs(vals_res$cv.lambda - vals_res$oracle.lambda) > 0]), l$fitted.values)
# 
# # is the val true loss actually locally convex in lambda?
# # It looks very convex locally!
# set.seed(10)
# small_train_dataset <- Make_linreg_dataset(5, beta, eps_factor = eps_factor)
# n_val <- 3
# val_dataset <- Make_linreg_dataset(n_val, beta, eps_factor = eps_factor)
# lambdas <- seq(0.01, 6, 0.5)
# beta_hats <- lapply(lambdas, function(lambda) {
#     Get_beta_hat(small_train_dataset, lambda)
# })
# dataset <- list(
#     beta=beta,
#     n_train=length(small_train_dataset$y),
#     n_val=length(val_dataset$y),
#     train_X=small_train_dataset$X,
#     val_X=val_dataset$X,
#     test_X=test_dataset$X,
#     train_y=small_train_dataset$y,
#     val_y=val_dataset$y,
#     val_true_y=val_dataset$true_y,
#     test_y=test_dataset$y,
#     test_true_y=test_dataset$true_y
# )
# res <- lapply(seq(length(lambdas)), function(i) {
#     lambda <- lambdas[i]
#     beta_hat <- beta_hats[[i]]
#     loss_res <- Get_losses(beta_hat, dataset)
#     data.frame(
#         lambda=as.numeric(lambda),
#         loss_res
#     )
# })
# ress <- do.call("rbind", res)
# plot(ress$lambda, ress$val_true_loss)
# poly_fit <- lm(val_true_loss ~ lambda + I(lambda^2), ress)
# lines(ress$lambda, poly_fit$fitted.values, col="red")
# 
# # oracle_lambdas_fixed_res <- lapply(n_vals, function(n_val) {
# #     print(paste("n_val", n_val))
# #     lambdas <- 10^seq(1.8, 2.25, by=2/n_val)
# #     beta_diffs <- lapply(lambdas, function(lambda) {
# #         beta_hat <- Get_beta_hat(train_dataset, lambda)
# #         Get_norm2(beta_hat - beta)
# #     })
# #     oracle_lam <- lambdas[which.min(beta_diffs)]
# #     best_beta_hat <- Get_beta_hat(train_dataset, oracle_lam)
# #     data.frame(
# #         n_val=n_val,
# #         oracle_lam_beta_diff=oracle_lam,
# #         oracle_loss_beta_diff=Get_norm2(test_dataset$true_y - test_dataset$X %*% best_beta_hat)  
# #     )
# # })
# # oracle_lambdas_fixed_res <- do.call("rbind", oracle_lambdas_fixed_res)
# # m <- merge(vals_res, oracle_lambdas_fixed_res, by="n_val")
# # m$loss_diff_beta_diff <- m$cv.val_true_loss - m$oracle_loss_beta_diff
# #     
# # summary(lm(log(loss_diff_beta_diff) ~ log(n_val), m[m$loss_diff_beta_diff > 0, ]))
# 

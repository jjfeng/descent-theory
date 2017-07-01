source("bspline_multi.R")

####### A HOMOGENEOUS FUNCTION HERE
do_homogeneous_func <- function(X) {
    y = 0
    for (i in seq(length(X))) {
        y = y +  sin(X[[i]])
    }
    y
}

Plot_bpline_true_fit <- function(dataset, observations_idx, dataset_bspline_matrix, oracle_res, cv_res) {
    num_knots <- ncol(dataset_bspline_matrix)
    par(mfrow=c(1,2))
    X1 <- dataset$X[[1]][observations_idx]
    true_y <- sin(X1)
    plot(X1, true_y, ylim = c(-5, 5)) # truth
    points(X1, dataset_bspline_matrix[,seq(num_knots)] %*% cv_res$coef[seq(num_knots)], col="red")
    plot(X1, true_y, ylim = c(-5, 5)) # truth
    points(X1, dataset_bspline_matrix[,seq(num_knots)] %*% oracle_res$coef[seq(num_knots)], col="green")
    
    # X2 <- dataset$X[[2]][observations_idx]
    # plot(X2, dataset_bspline_matrix[,seq(num_knots)] %*% cv_res$coef[seq(num_knots)], col="red")
    # plot(X2, dataset_bspline_matrix[,seq(num_knots)] %*% oracle_res$coef[seq(num_knots)], col="green")
}

set.seed(1)

num_p <- 8
n_trains <- 200
true_func <- do_homogeneous_func
num_lams <- c(1,2,4,8)
n_reps <- 10
snr <- 2

n_cores <- 3

cl <- makeCluster(n_cores)
clusterExport(cl,
              c("bdiag",
                "Do_bspline_CV_oracle", 
                "Get_bspline_losses", 
                "Make_lambda_grid", 
                "Eval_losses", 
                "Get_best_lambdas", 
                "Fit_bspline", 
                "Get_norm2")
)


cv_to_oracle_hom <- lapply(num_lams, function(n_lams) {
    cv_to_oracle_ntrains <- lapply(n_trains, function(n_train){
        n_val = n_train
        print(paste("train", n_train, "val", n_val))
        cv_oracle <- Do_bspline_cv_oracle_repl(
            reps=n_reps, 
            n_train=n_train,
            n_validate=n_val,
            num_p=num_p, 
            n_lams=n_lams,
            true_func=true_func, 
            plot_fnc = Plot_bpline_true_fit,
            snr=snr
        )
        print(colMeans(cv_oracle))
        data.frame(
            n=n_lams,
            cv_oracle
        )
    })
    do.call("rbind", cv_to_oracle_ntrains)
})
cv_to_oracle_hom <- do.call("rbind", cv_to_oracle_hom)
save(cv_to_oracle_hom, file = "cv_to_oracle_hom_homogeneous.RData")
print(cv_to_oracle_hom)

stopCluster(cl)

pdf('figures/validation_size_loss_diff_homogeneous.pdf', width=8, height=5)
par(mar=c(5,5,1,1), mfrow=c(1,1))
plot(
    jitter(cv_to_oracle_hom$n, 0.15),
    cv_to_oracle_hom$loss_diff_sq,
    ylab="Validation Loss Difference",
    xlab="Number of penalty parameters",
    cex=0.7,
    cex.axis=1.2,
    cex.lab=1.4,
    log="xy",
    xaxt = "n"
)
axis(1, at=num_lams, labels=num_lams, cex.axis=1.2)
lmdata <- cv_to_oracle_hom[cv_to_oracle_hom$loss_diff_sq > 1e-5,]
fit <- lm(log(loss_diff_sq) ~ log(n), lmdata)
fit_summ <- summary(fit)
lines(
    lmdata$n,
    exp(fit$fitted.values)
)
legend(
    "bottomright",
    legend=paste(
        "Slope:",
        format(fit_summ$coefficients[2, 1], digits=4),
        "+/-",
        format(1.96 * fit_summ$coefficients[2,2], digits=4)
    ),
    cex = 1.2
)
dev.off()

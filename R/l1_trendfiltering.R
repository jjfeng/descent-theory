MIN_X <- 0
MAX_X <- 1
BY_X <- 0.01
EPSILON_SD <- 0.01
SUBDIVISIONS <- 2000
TOL <- 0.001

# Code copied from http://www.r-bloggers.com/filtering-data-with-l1-regularisation/
l1filter.optim <- function (observed_y, lambda = 0.0) {
    n <- length(observed_y)
    
    Get_penalty <- function(fitted_y) {
        P2 = 0
        for (i in 2:(n-1)) {
            P2 = P2 + abs(fitted_y[i-1] - 2 * fitted_y[i] + fitted_y[i+1])
        }
        P2
    }
    Get_loss <- function(fitted_y) {
        0.5/n * sum((fitted_y - observed_y)**2)
    }
    Get_objective <- function(fitted_y, lambda) {
        P1 = Get_loss(fitted_y)
        P2 = Get_penalty(fitted_y)
        P1 + lambda * P2
    }
    
    fit = optim(observed_y, Get_objective, lambda = lambda, method = "BFGS", control = list(maxit = 100000, type = 3))
    
    if (fit$convergence != 0) {
        warning(sprintf("Optimisation failed to converge! (lambda = %f)", lambda))
        print(fit)
    }
    if (fit$value != Get_loss(fit$par) + lambda * Get_penalty(fit$par)) {
        warning(sprintf("Optimisation values dont match! (lambda = %f)", lambda))
        print(fit$value)
        print(Get_loss(fit$par) + lambda * Get_penalty(fit$par))
    }
    
    data.frame(
        lambda=lambda,
        loss=Get_loss(fit$par),
        penalty=Get_penalty(fit$par)
    )
}

set.seed(5)
# Generate sawtooth data
y <- rep(0:10, 5)
y <- y + rnorm(length(y))

# Fit functions for different lambda values
LAMBDAS <- c(0,10**(-6:2))

l1_trendfilter_res <- lapply(LAMBDAS, function(lam) {l1filter.optim(y, lam)})
l1_trendfilter_res <- do.call(rbind, l1_trendfilter_res)

# why the fuck is this so bumpy. penalty is supposed to be monotonic. the optimization function must not be working
plot(
    l1_trendfilter_res$lambda, 
    l1_trendfilter_res$penalty, 
    "b", 
    log="y", 
    ylab="Penalty",
    xlab="Lambda",
    main="L1 Trend Filtering"
)
plot(
    l1_trendfilter_res$lambda, 
    l1_trendfilter_res$loss, 
    "b", 
    ylab="Loss", 
    xlab="Lambda",
    main="L1 Trend Filtering"
)



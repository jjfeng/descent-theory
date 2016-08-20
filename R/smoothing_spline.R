MIN_X <- 0
MAX_X <- 10
MAX_N <- 1000
EPSILON_SD <- 0.1
SUBDIVISIONS <- 2000
TOL <- 0.001

set.seed(10)
# Generate data
x <- seq(MIN_X, MAX_X, (MAX_X - MIN_X)/MAX_N)

epsilon <- rnorm(n=length(x), sd=EPSILON_SD)
y <- sin(x) + epsilon

shuffled_indices <- sample(seq(1, MAX_N), MAX_N)

# Single Example
sub_ind <- shuffled_indices[1:20]
test_f <- smooth.spline(x[sub_ind], y[sub_ind], spar=0, cv=NA, all.knots = T)
plot(x[sub_ind], y[sub_ind])
lines(smoothed_f)
test_penalty <- integrate(function(u){
    predict(test_f, u, deriv=2)$y^2
}, min_x, max_x, subdivisions=SUBDIVISIONS, rel.tol=TOL)

Fit_model <- function(spar, indices) {
    y_t <- y[indices]
    x_t <- x[indices]
    smoothed_f <- smooth.spline(x_t, y_t, spar=spar, cv=NA, all.knots = T)
    
    Get_Penalty_Integrand <- function(u){
        predict(smoothed_f, u, deriv=2)$y^2
    }
    penalty <- integrate(Get_Penalty_Integrand, min_x, max_x, subdivisions=SUBDIVISIONS, rel.tol=TOL)
    
    fitted_y_t <- predict(smoothed_f, x = x_t)
    train_loss <- mean((y_t - fitted_y_t$y)^2)
    
    y_v <- y[-indices]
    x_v <- x[-indices]
    fitted_y_z <- predict(smoothed_f, x = x_v)
    v_loss <- mean((y_v - fitted_y_z$y)^2)
    
    data.frame(
        train_loss=train_loss,
        v_loss=v_loss,
        n=length(indices),
        penalty=penalty$value,
        lam=smoothed_f$lambda
    )
}

# Now run multiple lambdas with fixed sample size to see what happens
SPAR_LIST <- seq(-0.5, 1.5, 0.1)
NUM_TRAIN <- 100
smooth_res <- do.call("rbind", lapply(SPAR_LIST, Fit_model,  indices=shuffled_indices[1:NUM_TRAIN]))

plot(
    smooth_res$lam, 
    smooth_res$penalty, 
    type="l",
    log="y", 
    ylab="Penalty", 
    xlab="Lambda",
    main="Smoothing Spline (Penalty=2nd-order Sobolev norm)"
)
lines(
    smooth_res$lam, 
    smooth_res$train_loss, 
    col="red"
)
lines(
    smooth_res$lam, 
    smooth_res$v_loss, 
    col="green"
)


# Now run multiple sample sizes with fixed lambda to see what happens
TRAIN_SIZE_LIST <- seq(100, MAX_N, 100)
SINGLE_SPAR <- -0.1
sample_size_res <- do.call("rbind", lapply(TRAIN_SIZE_LIST, function(train_size) {
    Fit_model(SINGLE_SPAR, shuffled_indices[1:train_size])
}))

plot(
    sample_size_res$n, 
    sample_size_res$penalty, 
    "b", 
    log="y", 
    ylab="Penalty", 
    xlab="Sample Size",
    main="Smoothing Spline (Penalty=2nd-order Sobolev norm)"
)

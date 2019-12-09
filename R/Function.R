##' Functions used in ISTA algrithm
##'
##' @title Functions used in ISTA algrithm
##' @author Yuan Zhou, Lingsong Meng
##' @export

# Simulation function
func.sim <- function(index, x, noise) {
    set.seed(index)
    p <- ncol(x)
    
    ## Polynomial
    e <- c(1/3, 2, 3)
    c <- list()
    y <- 0
    
    for (i in 1:3) {
        c[[i]] <- as.matrix(runif(p, -2, 2))
        y <- cbind(y, (x^e[i]) %*% c[[i]])
    }
    
    vars <- matrix(rnorm(length(y), 0, 1), nrow = nrow(y))
    y <- as.matrix(y) + vars * noise
    
    return(y)
}

# Split datasets
subs <- function(objects, rows) {
    if (is.matrix(objects)) 
        return(objects[rows, , drop = FALSE])
    if (is.array(objects)) 
        return(objects[rows, , ])
}

split <- function(objects, ratio, name1, name2) {
    n <- nrow(objects[[1]])
    idx <- sample(1:n)
    assign(name1, idx[1:ceiling(ratio * n)])
    assign(name2, idx[(ceiling(ratio * n) + 1):n])
    names.in.list <- mget(c(name1, name2))
    for (name in names(objects)) {
        list2env(objects, envir = environment())
        assign(paste0(name, ".", name1), subs(get(name), 
            get(name1)))
        assign(paste0(name, ".", name2), subs(get(name), 
            get(name2)))
        names.in.list.new <- mget(c(paste0(name, 
            ".", name1), paste0(name, ".", name2)))
        names.in.list <- c(names.in.list, names.in.list.new)
    }
    return(names.in.list)
}

# Objective function
Obj.func <- function(Y, X, beta, step_size) {
    g.value <- tcrossprod(t(Y - X %*% beta))/(2 * 
        step_size)
    # g.null <- tcrossprod(t(Y-
    # mean(Y)))/(2*step_size)
    return(g.value)
}


# Soft-thresholding operator
S <- function(lambda, X) {
    s <- function(lambda0, x) {
        if (x > lambda0) {
            result0 <- x - lambda0
        } else if (x < (-1) * lambda0) {
            result0 <- x + lambda0
        } else result0 <- 0
        
        return(result0)
    }
    result <- apply(X, 1, s, lambda0 = lambda)
    return(result)
}

# Gradient
Grad <- function(Y, X, parameters, lambda, step_size) {
    grad <- (-1) * t(X) %*% (Y - X %*% parameters)
    grad.g <- (parameters - S(lambda, X = parameters - 
        step_size * grad))/step_size
    Grads <- mget(c("grad", "grad.g"))
    return(Grads)
}

# Backtrcking line search
backtracking <- function(Y, X, parameters, grad, 
    grad.g, lambda, step_size, constant) {
    # list2env(Grads, envir = environment())
    # list2env(grad, envir = environment())
    # list2env(grad.g, envir = environment())
    
    while (Obj.func(Y, X, beta = parameters - 
        step_size * grad.g, step_size = step_size) * 
        step_size > Obj.func(Y, X, parameters, 
        step_size = step_size) * step_size - 
        step_size * t(grad) %*% grad.g + (step_size/2) * 
        tcrossprod(t(grad.g))) {
        step_size <- constant * step_size
        grad.g <- (parameters - S(lambda, (parameters - 
            step_size * grad)))/step_size
    }
    return(step_size)
}

# Update parameters
update.parameters <- function(lambda, grad, parameters, 
    step_size) {
    parameters_new = S(lambda, parameters - step_size * 
        grad)
    return(parameters_new)
}

# Predict
pred <- function(Y.test, X.test, parameters) {
    Y.test.hat <- X.test %*% parameters
    test.error <- Obj.func(Y.test, X.test, parameters, 
        step_size = 1)
    return(mget(c("Y.test.hat", "test.error")))
}

# Main algrithm
ISTA0 <- function(X.subtrain, Y.subtrain, X.valid, 
    Y.valid, epoch = 10000, patience, lambda) {
    # --------- Initializa parameter ---------
    p <- ncol(X.subtrain)
    parameters <- runif(p, -0.5, 0.5)
    step_size <- 0.001
    valid.error.best <- Obj.func(Y.valid, X.valid, 
        parameters, step_size)
    parameters.best <- parameters
    k = 0
    
    # --------- Calculate gradient and update
    # parameters -----------
    for (i in 1:epoch) {
        Grads <- Grad(Y.subtrain, X.subtrain, 
            parameters, lambda = lambda, step_size)
        step_size_update <- backtracking(Y.subtrain, 
            X.subtrain, Grads$grad, Grads$grad.g, 
            parameters, lambda = lambda, step_size, 
            constant = 1/2)
        parameters <- update.parameters(lambda = lambda, 
            Grads$grad, parameters, step_size_update)
        valid.error <- Obj.func(Y.valid, X.valid, 
            parameters, step_size)
        
        if (valid.error < valid.error.best) {
            valid.error.best <- valid.error
            parameters.best <- parameters
        } else {
            k = k + 1
        }
        # cat(i, valid.error, valid.error.best, k,
        # '\n')
        if (k == patience) {
            return(mget(c("i", "valid.error.best", 
                "parameters.best")))
        }
    }
    return(mget(c("i", "valid.error.best", "parameters.best")))
}


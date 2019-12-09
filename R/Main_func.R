##' Iterative Soft_thresholding Algorithm (ISTA)
##'
##' ISTA algrithm combined with Backtracking line search, using Cross-validation/AIC/BIC to select optimal parameter
##' @title Main function
##' @param data.X Predictor matrix
##' @param data.Y Response matrix
##' @param scale Whether the matrix should be scaled when applying algrithm, default is TRUE
##' @param lambda Parameter of LASSO regression, can be a vector with length > 1, default is 1e-2
##' @param method Method used to select optimal parameters, can only chosen from 'CV', 'AIC' or 'BIC', default is NULL
##' @param epoch The maximum number of steps allowed when training the model, default is 1e4
##' @param patience Used to decide when end training after a minima occured, default is 50
##' @return MSE between prediction values and true values of response'
##' @author Yuan Zhou, Lingsong Meng
##' @export

# library(MASS)
# library(Metrics)
# source("Function.R")

ISTA.main <- function(data.X, data.Y, scale = T,
    lambda = 0.01, method = NULL, epoch = 10000,
    patience = 50) {

    stopifnot(!missing(data.X), !missing(data.Y),
        is.matrix(data.X), is.matrix(data.Y),
        ncol(data.Y) == 1)

    if (scale == T) {
        data.X <- scale(data.X)
        data.Y <- scale(data.Y)
    }

    list2env(split(mget(c("data.X", "data.Y")),
        0.8, "train", "test"), envir = environment())

    if (is.null(method)) {
        for (i in 1:length(lambda)) {
            parameters.best <- ISTA0(data.X.train,
                data.Y.train, data.X.test, data.Y.test,
                epoch = epoch, patience = patience,
                lambda = lambda[i])$parameters.best
            Y.test.hat <- pred(data.Y.test, data.X.test,
                parameters.best)$Y.test.hat
            test.mse <- mse(data.Y.test, Y.test.hat)
            cat("Test MSE in test set with lambda",
                lambda[i], "is", test.mse, "\n")
        }
        return(test.mse)
    } else {
        stopifnot(method %in% c("CV", "AIC",
            "BIC"))
        p <- ncol(data.X)
        result.list <- list()
        # --------- Split datasets and find best
        # labmda ------------
        if ("CV" %in% method) {
            fold = 10
            result.cv <- data.frame(lambda = seq(0,
                0, len = length(lambda)), valid.error.mean = seq(0,
                0, len = length(lambda)), epoch = seq(0,
                0, len = length(lambda)))
            valid.error.cv <- numeric(fold)
            training.epoch <- numeric(fold)
            for (i in 1:length(lambda)) {
                data.train <- data.frame(data.X.train,
                  data.Y.train)
                data.train.s <- split.data.frame(data.train,
                  factor(1:fold))
                for (j in 1:fold) {
                  data.X.valid <- as.matrix(data.train.s[[j]][,
                    1:p])
                  data.Y.valid <- as.matrix(data.train.s[[j]][,
                    "data.Y.train"])
                  data.X.subtrain <- as.matrix(data.train[-as.numeric(row.names(data.X.valid)),
                    1:p])
                  data.Y.subtrain <- as.matrix(data.train[-as.numeric(row.names(data.X.valid)),
                    "data.Y.train"])

                  training.result <- ISTA0(data.X.subtrain,
                    data.Y.subtrain, data.X.valid,
                    data.Y.valid, epoch = epoch,
                    patience = patience, lambda = lambda[i])
                  valid.error.cv[j] <- training.result$valid.error.best
                  training.epoch[j] <- training.result$i
                }
                result.cv[i, "lambda"] <- lambda[i]
                result.cv[i, "valid.error.mean"] <- mean(valid.error.cv)
                result.cv[i, "epoch"] <- mean(training.epoch)
            }

            cat("Best lambda under cross validation criteria is",
                result.cv[which.min(result.cv$valid.error.mean),
                  "lambda"], "\n")
            result.list[["CV"]] <- result.cv
        }
        if ("AIC" %in% method || "BIC" %in% method) {
            list2env(split(list(data.X = data.X.train,
                data.Y = data.Y.train), 0.75,
                "subtrain", "valid"), envir = environment())
            result.aic <- data.frame(lambda = seq(0,
                0, len = length(lambda)), AIC = seq(0,
                0, len = length(lambda)), epoch = seq(0,
                0, len = length(lambda)))
            result.bic <- data.frame(lambda = seq(0,
                0, len = length(lambda)), BIC = seq(0,
                0, len = length(lambda)), epoch = seq(0,
                0, len = length(lambda)))

            for (i in 1:length(lambda)) {
                training.result <- ISTA0(data.X.subtrain,
                  data.Y.subtrain, data.X.valid,
                  data.Y.valid, epoch = epoch,
                  patience = patience, lambda = lambda[i])
                parameters.best <- training.result$parameters.best

                if ("AIC" %in% method) {
                  aic <- tcrossprod(t(data.Y.valid -
                    data.X.valid %*% parameters.best))/(nrow(data.X.valid) *
                    var(data.Y.valid)) + 2 *
                    (length(which(parameters.best ==
                      0)))/nrow(data.X.valid)
                  result.aic[i, "lambda"] <- lambda[i]
                  result.aic[i, "AIC"] <- aic
                  result.aic[i, "epoch"] <- training.result$i

                }
                if ("BIC" %in% method) {
                  bic <- tcrossprod(t(data.Y.valid -
                    data.X.valid %*% parameters.best))/(nrow(data.X.valid) *
                    var(data.Y.valid)) + log(nrow(data.X.valid)) *
                    (length(which(parameters.best ==
                      0)))/nrow(data.X.valid)
                  result.bic[i, "lambda"] <- lambda[i]
                  result.bic[i, "BIC"] <- bic
                  result.bic[i, "epoch"] <- training.result$i

                }
            }
            if ("AIC" %in% method) {
                cat("Best lambda under AIC criteria is",
                  result.aic[which.min(result.aic$AIC),
                    "lambda"], "\n")
                result.list[["AIC"]] <- result.aic
            }
            if ("BIC" %in% method) {
                cat("Best lambda under BIC criteria is",
                  result.bic[which.min(result.bic$BIC),
                    "lambda"], "\n")
                result.list[["BIC"]] <- result.bic
            }
        }

        # ---------- Using model in test data
        # --------------
        result.mse <- c("MSE")
        for (med in method) {
            result <- result.list[[med]]
            lambda.best <- result[which.min(result[,
                2]), "lambda"]
            parameters.best <- ISTA0(data.X.subtrain,
                data.Y.subtrain, data.X.valid,
                data.Y.valid, epoch = epoch,
                patience = patience, lambda = lambda.best)$parameters.best
            Y.test.hat <- pred(data.Y.test, data.X.test,
                parameters.best)$Y.test.hat
            test.mse <- mse(data.Y.test, Y.test.hat)
            cat("Test MSE in test set with best lambda under",
                med, "is", test.mse, "\n")
            result.mse <- cbind(result.mse, test.mse)
            colnames(result.mse)[ncol(result.mse)] <- med
        }
    }
    return(result.mse)
}



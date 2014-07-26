library("Hmisc")

# obs and pred need to be numeric for all functions!

# function to calculate the concordance index
cIDX <- function(pred, obs) {
  # calculate the concordance index 
  rcorr.cens(pred, obs)[1]
}

# function to calculate the root mean square error
# observed (independent variable, x)
# predicted (dependent variable, y)
RMSE <- function(obs, pred) {
  sqrt(mean((obs - pred)^2, na.rm = TRUE)) 
}

# coefficient of determination / R squared
R2 <- function(obs, pred) {
  # total variation in the observed data
  SSyy = sum((obs-mean(obs))^2)
  # sqared error of prediction
  SSE = sum((pred-obs)^2)
  1 - SSE / SSyy
}
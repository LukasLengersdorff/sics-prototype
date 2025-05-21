### Helpers ###

logit = function(y) log(y/(1-y))

expit = function(x) 1/(1+exp(-x))

log_expminexp = function(x,y) log(1-exp(y-x)) + x

pinv = function(X) solve(t(X) %*% X) %*% t(X)

logdiff = function(x,y) { # x > y
  x + log(1-exp(y-x))
}

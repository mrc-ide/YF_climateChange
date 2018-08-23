##################################################################################
### fun_glm_lnL ### FROM fun.glm.lnL ###
##################################################################################

# Calculating lnL of glm regression:
fun_glm_lnL = function(beta, x, y) {
  # beta are the model coefficients,
  # x the independent covariates (needs one column for the Intercept),
  # and y the binary outcome_
  # model predictions pi = p, 1-pi = q
  eta = as.numeric(x %*% beta)
  logq = -exp(eta) # -exp(X*beta) (= log(1-q) )
  logp = log(1-exp(logq)) # 
  #q = exp(logq)
  #p = 1-q
  logl = sum(logp[y==1]) + sum(logq[y==0])
  return (logl)
}
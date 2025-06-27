# =========================================================================
# make a link function
# =========================================================================
makeLinkStructure <- function(link = link){
  
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
    if (linktemp == "link") {
      linktemp <- eval(link)
    }
  }
  
  if (any(linktemp == c("logit", "probit", "cloglog"))) {
    
    stats <- make.link(linktemp)
    
  } else if (linktemp == "loglog") {
    
    stats <- list()
    
    stats$linkfun <- function(mu) -log(-log(mu))
    stats$linkinv <- function(eta) exp(-exp(-eta))
    stats$mu.eta <- function(eta) exp(-exp(-eta)) * exp(-eta)
    
  } else {
    stop(paste(
      linktemp, "link not available, available links are \"logit\", ", 
      "\"probit\" and \"cloglog\""
    ))
  }
  
  link1 <- structure(list(
    link = linktemp,
    linkfun = stats$linkfun,
    linkinv = stats$linkinv,
    mu.eta = stats$mu.eta
  ))
  
  linkfun <- link1$linkfun
  linkinv <- link1$linkinv
  mu.eta <- link1$mu.eta
  
  return(list(linkfun = linkfun, 
              linkinv = linkinv, 
              mu.eta = mu.eta))
  
}
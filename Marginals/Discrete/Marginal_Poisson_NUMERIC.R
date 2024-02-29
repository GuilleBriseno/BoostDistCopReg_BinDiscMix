################################################################################
### DESCRIPTION
### This file defines the Log-Normal-Marginal function. It contains the 
### expressions of the pdf, cdf, its derivatives wrt the two parameters, names,
### response functions, offset functions, check_y functions required to create an mboost Family object.
### This function is handed over to a copula function and the expressions and functions
### are merged with the respective copula via the rlang package, in order to create 
### appropriate Families objects from gamboostLSS.


### libraries 
library(mboost)
library(gamboostLSS)
library(rlang)



# 
# # helper function for discrete margins CDFs:
# 
# # This code creates a matrix that is very useful for computing the derivatives of the CDF of a discrete margin w.r.t. the distribution parameters
# y1 <- rpois(100, lambda = 30)
# 
# aux_matrix_disc <- function(response){
#   
#   ly1 <- length(response)
#   y1m <- list()
#   my1 <- max(response)
#   
#   for(i in 1:ly1){ 
#     
#     y1m[[i]] <- seq(0, response[i]); 
#     
#     length(y1m[[i]]) <- my1 + 1
#     
#   }
#   
#   
#   y1m <- do.call(rbind, y1m) 
#   
#   
#   return(y1m)
# }
# 
# 
# aux_matrix_disc(y1)
# 
# 
# ly1 <- length(y1)
# y1m <- list()
# my1 <- max(y1)
# 
# for(i in 1:ly1){ 
#   
#   y1m[[i]] <- seq(0, y1[i]); # Create a sequence of integers from 0 to y_i
#   
#   # Fill up this row up to the maximum observed in the discrete response
#   length(y1m[[i]]) <- my1 + 1 # Add 1 to accommodate P( Y = 0) (since it is not zero.truncated)
#   
#   } 
# 
# 
# 
# # zero truncated poisson
# # for(i in 1:ly1){ y1m[[i]] <- seq(1, y1[i]); length(y1m[[i]]) <- my1} 
# 
# y1m <- do.call(rbind, y1m)   
# 
# 
# 
# derpdf2.dermu2FUNCpo <- function(y2, mu2) exp(-mu2) * (mu2^(y2 - 1) * y2 - mu2^y2)/factorial(y2) 
# 
# 
# 
# prec <- pmax(53, Rmpfr::getPrec(mu), Rmpfr::getPrec(y1))
# 
# mu <- Rmpfr::mpfr(mu, prec)
# y1  <- Rmpfr::mpfr( y1, prec)        
# 
# 
# # This apply statement generates these derivatives, in fact, it is IDENTICAL to the output of the function!
# sapply(y1, function(i) 
#   sum( sapply(0:i, function(j) 
#   exp(-mu) * (mu^(j - 1) * j - mu^j)/factorial(j) )
#   ))
# 
# 
# # Evaluate gradient of PDF
# derpdf2.dermu2       <- as.numeric( derpdf2.dermu2FUNCpo(y1, mu) )
# 
# # Gradient of CDF
# derp2.dermu2           <- rowSums( matrix(as.numeric(derpdf2.dermu2FUNCpo(y1m, mu)), dim(y1m)[1],dim(y1m)[2]), na.rm = TRUE ) 



Poisson_Mar <- function(loc = NULL, offset_mu = NULL, aux_mat = NULL, dim_aux_mat = NULL){
  
  
  # check for appropriate offset values for parameter sigma
  if ((!is.null(offset_mu) && offset_mu <= 0))
    stop(sQuote("mu"), paste(" must be greater than zero in Marginal", loc))
  
  
  
  # creating the location specific expressions
  y <- parse_expr(paste("y[,", loc, "]", sep = ""))                      
  mu <- parse_expr(paste("mu", loc, sep = ""))                           
  
  
  
  ### generic functions
  
  # pdf
  
  pdf_gen <- expr( dpois(x = !!y, lambda = !!mu) )
  # cdf
  cdf_gen <- expr( ppois(q = !!y, lambda = !!mu) )  
  
  
  # NOT USED
  cdfmpdf_gen <- expr( pdffz( ppois(q = (!!y), lambda = !!mu) - dpois(x = (!!y), lambda = !!mu) ) )
  
  
  l_generic <- list(pdf = pdf_gen,
                    cdf = cdf_gen,
                    cdfmpdf = cdfmpdf_gen)
  
  
  
  ### mu functions
  
  # pdf
  pdf_mu  <- expr( dpois(x = !!y, lambda = exp(f) ) )
  
  # der pdf 
  derpdf1.deretamu   <- expr({ 
    
    para <- exp(f) 
    
    grd <- grad(func  = function(para) dpois(!!y, lambda = para), x = para)
    
    grd_full <- grd * exp(f)
    
    return(grd_full)
    
  })
  
  
  
  # cdf
  cdf_mu  <- expr( ppois(q = !!y, lambda = exp(f) ) )  
  
  # derivative cdf 
  #dercdf.deretamu <-  expr( rowSums(t(sapply(1:(!!dim_aux_mat)[1], function(i) derPOISSON_PDF.dermu2(!!aux_mat[i,], exp(f)[i]))), na.rm = TRUE) * exp(f) )
  dercdf.deretamu <- expr({ 
    
    para <- exp(f) 
    
    grd <- grad(func  = function(para) ppois(!!y, lambda = para), x = para)
    
    grd_full <- grd * exp(f)
    
    return(grd_full)
    
  })
  
  
  # NOT USED
  cdfmpdf_mu <- expr(  pdffz( ppois(q = !!y, lambda = exp(f) ) - dpois(x = !!y, lambda = exp(f) ) ) )
  
  # NOT USED
  dercdfmpdf.deretamu <- expr( sapply((!!y), function(i) sum( sapply(0:i, function(j) 
    ( ( exp( -exp(f) ) * (exp(f)^(j - 1) * j - exp(f)^j ) ) / factorial(j) ) * exp(f) )  ) 
  ) - 
    ( ( exp( -exp(f) ) * (exp(f)^(!!y - 1 ) * (!!y) - exp(f)^(!!y) ) ) / factorial(!!y) ) * exp(f)
  )
  
  l_mu <- list(pdf = pdf_mu,
               #derpdf = derpdf.deretamu,
               derpdf = derpdf1.deretamu,
               cdf = cdf_mu,
               dercdf = dercdf.deretamu,
               cdfmpdf = cdfmpdf_mu, 
               dercdfmpdf = dercdfmpdf.deretamu)
  
  
  
  ### response functions 
  
  response_mu <- function(f) exp(f)
  
  l_response <- list(mu = response_mu)
  
  ### offset functions     
  
  offset_mu <- expr({
    if (!is.null(!!mu)) {
      RET <- log(!!mu)
    }
    else {
      
      #RET <- weighted.mean((log(!!y) + weighted.mean(log(!!y), w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
      #temp <-  weighted.mean((log(!!y) + weighted.mean(log(!!y), w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
      #RET <- temp
      
      temp <- weighted.mean((!!y + weighted.mean(!!y, w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE) # weighted version
      
      RET <- log(temp)     
      
    }
    return(RET)
  })
  
  
  
  l_offset <- list(mu = offset_mu)
  
  ### names 
  
  name_mu <- "Poisson distribution: mu(log link)"
  
  
  l_names <- list(mu = name_mu)
  
  ### check y function   
  
  #check_y <- function(y) y                                                      
  check_y_mu <- expr({
    if (!is.numeric(!!y))
      stop("response is not but ", sQuote("Poisson()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("Poisson()"))
    y
  })
  
  
  l_check_y <- list(mu = check_y_mu)
  
  
  # return list
  
  l = list(parameters = c(mu),
           parameter_names = c("mu"),
           
           generic = l_generic,
           mu = l_mu,
           
           response = l_response,
           offset = l_offset,
           name = l_names,
           check_y = l_check_y,
           marg_name = "PoissonMarg")
  
  
  return(l)
  
}

# marg <- LogNormal_Mar(loc = 1)

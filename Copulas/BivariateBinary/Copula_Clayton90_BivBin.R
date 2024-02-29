################################################################################
### DESCRIPTION
### Gauss_Cop() defines the function to create a Families object from gamboostLSS
### for the Gaussian copula with arbitrary marginal distributions.
### Marginal functions (like the ones implemented in the Marginals folder) 
### are handed over to the arguments marg1 and marg2. 
### Over the course of the function run, first generic bodies and arguments of the functions 
### for the loss, risk and negative gradient are defined for each parameter submodel, 
### i.e. the marginal distribution parameters and 
### the copula dependence parameter (also via the rlang package). 
### Thereafter, the parameter specific functions and Family objects (see mboost) 
### are created via the family_gen() function 
### (see Marginals folder) for each parameter sub-model. 
### Finally, a gamboostLSS Families object is created that suits the gamboostLSS 
### boosting algorithm implemented in the gamboostLSS package.


### load libraries 
library(mboost)
library(gamboostLSS)
library(VGAM)
library(rlang)


# load Marginal functions

source("Marginals/family_gen.R")
# source("Marginals/Marginal_LogNormal.R")
# source("Marginals/Marginal_LogLogistic.R")
source("Marginals/Binary/Marginal_probit.R")
source("Marginals/Binary/Marginal_logit_newoffset.R")
source("Marginals/Binary/Marginal_cloglog.R")
#source("BoosCopContMarg/Marginals/family_gen.R")
# source("BoosCopContMarg/Marginals/Marginal_LogNormal.R")
# source("BoosCopContMarg/Marginals/Marginal_LogLogistic.R")
# source("BoostCopBinMarg/Marginals/Marginal_probit.R")


### Gauss-Copula function (for bivariate binary margins)

Clayton90_Cop_BivBinary <- function(marg1 = NULL, marg2 = NULL,
                                  mu1 = NULL, sigma1 = NULL, nu1 = NULL, tau1 = NULL,
                                  mu2 = NULL, sigma2 = NULL, nu2 = NULL, tau2 = NULL,
                                  rho = NULL,
                                  stabilization = c("none", "MAD", "L2")){
  
  ################################################################################
  ########################## marginal checks #####################################
  
  if(is.null(marg1)) stop("First marginal distribution not defined.")
  if(is.null(marg2)) stop("Second marginal distribution not defined.")
  
  if(!(marg1 %in% c("LOGIT", "PROBIT", "CLOGLOG", "GEVLINK"))) stop("First Marginal distribution not available.")
  if(!(marg2 %in% c("LOGIT", "PROBIT", "CLOGLOG", "GEVLINK"))) stop("Second Marginal distribution not available.")
  
  
  ################################################################################
  #################### calling the marginal distributions ########################
  
  if(marg1 == "NO"){
    marg1 <- Gaussian_Mar(loc = 1, offset_sigma = sigma1)
  }else if(marg1 == "SM") {
    marg1 <- SM_Mar(loc = 1, offset_mu = mu1, offset_sigma = sigma1, offset_tau = tau1)
  } else if(marg1 == "LOGNO") {
    marg1 <- LogNormal_Mar(loc = 1, offset_sigma = sigma1)
  } else if(marg1 == "GA") {
    marg1 <- Gamma_Mar(loc = 1, offset_mu = mu1, offset_sigma = sigma1)
  } else if(marg1 == "DAGUM") {
    marg1 <- Dagum_Mar(loc = 1, offset_mu = mu1, offset_sigma = sigma1, offset_nu = nu1)
  } else if(marg1 == "LOGLOG") {
    marg1 <- LogLogistic_Mar(loc = 1, offset_mu = mu1, offset_sigma = sigma1)
  } else if(marg1 == "LOGIT") {                                                 # BINARY RESPONSES
    marg1 <- BernoulliLogit_Mar(loc = 1, offset_mu = mu1)
  } else if(marg1 == "PROBIT") {
    marg1 <- BernoulliProbit_Mar(loc = 1, offset_mu = mu1)
  } else if(marg1 == "CLOGLOG") {
    marg1 <- BernoulliCloglog_Mar(loc = 1, offset_mu = mu1)
  }
  
  
  if(marg2 == "NO"){
    marg2 <- Gaussian_Mar(loc = 2, offset_sigma = sigma2)
  } else if(marg2 == "SM") {
    marg2 <- SM_Mar(loc = 2, offset_mu = mu2, offset_sigma = sigma2, offset_tau = tau2)
  }else if(marg2 == "LOGNO") {
    marg2 <- LogNormal_Mar(loc = 2, offset_sigma = sigma2)
  }else if(marg2 == "GA") {
    marg2 <- Gamma_Mar(loc = 2, offset_mu = mu2, offset_sigma = sigma2)
  }else if(marg2 == "DAGUM") {
    marg2 <- Dagum_Mar(loc = 2, offset_mu = mu2, offset_sigma = sigma2, offset_nu = nu2)
  }else if(marg2 == "LOGLOG") {
    marg2 <- LogLogistic_Mar(loc = 2, offset_mu = mu2, offset_sigma = sigma2)
  } else if(marg2 == "LOGIT") {                                                 # BINARY RESPONSES
    marg2 <- BernoulliLogit_Mar(loc = 2, offset_mu = mu2)
  } else if(marg2 == "PROBIT") {
    marg2 <- BernoulliProbit_Mar(loc = 2, offset_mu = mu2)
  } else if(marg2 == "CLOGLOG") {
    marg2 <- BernoulliCloglog_Mar(loc = 2, offset_mu = mu2) }
  
  
  ################################################################################
  ##################### check offsets of copula parameter ########################
  
  
  if ((!is.null(rho) && rho <= 0))
    stop(sQuote("rho"), " must be greater than 0.")
  
  
  ################################################################################
  ############################## Helpers #########################################
  ################### To be deleted when integration into package ################
  
  check_stabilization <- function(stabilization = c("none", "MAD", "L2")) {
    stabilization <- match.arg(stabilization)
    ## check if old stabilization interface is used and issue a warning
    if (getOption("gamboostLSS_stab_ngrad")) {
      warning("Usage of ", sQuote("options(gamboostLSS_stab_ngrad = TRUE)"),
              " is deprecated.\n", "Use argument ", sQuote("stabilization"),
              " in the fitting family. See ?Families for details.")
      if (stabilization == "none")
        warning(sQuote("stabilization"), " is set to ", dQuote("MAD"))
    }
    stabilization
  }
  
  
  ################################################################################
  ########################## Check stabilization #################################
  
  stabilization <- check_stabilization(stabilization)
  
  
  
  
  ################################################################################
  ############### Explanation of the automatic function creation: ################
  
  ### Aim: Automatically create loss, risk, gradient, offset and check_y function,
  ###      which are necessary for the family object construction (see mboost).
  
  ### Remark: Each parameter (marginals and copulas) has its individual functions.
  ###         The loss, risk and gradients consist partly of copula and partly of marginal elements,
  ###         which need to be combined appropriately.
  
  ### Procedure: All functions are created automatically. 
  ###            To do so, it is necessary to define appropriate arguments and appropriate bodies.
  ###            In the following, firstly the functions arguments are created, 
  ###            secondly the functions bodies (at least the generic functions) are created 
  ###            and finally for each parameter the appropriate functions and family object are constructed (via the family_gen function).
  
  ### for code implementation and the idea see body() and formals() function.
  
  
  
  ################################################################################
  ####################### Arguments for function creation ########################
  
  
  # Note: I decided to put the arguments creation in the copula function (and not in the family_gen function)
  #       because of potential problems with parameter specific loss creation. 
  #       Also conceptually it makes sense to have the whole function creation procedure in one place.
  
  # code explanation: Create the arguments for the loss, risk, gradient, offset and check_y function.
  # The risk, gradient, offset and check_y arguments are equal for all parameters.
  # For the loss function parameters specific arguments need to be created via the args_loss_creator. 
  # for further information see also https://stackoverflow.com/questions/17751862/create-a-variable-length-alist
  
  
  # arguments for gradient and risk functions
  args_grad_risk <- c("y", "f", "w")
  l_args_grad_risk <- rep(list(expr()), length(args_grad_risk)) 
  names(l_args_grad_risk) <- args_grad_risk
  l_args_grad_risk[["w"]] <- 1
  
  # arguments for offset functions
  args_offset <- c("y", "w")
  l_args_offset <- rep(list(expr()), length(args_offset)) 
  names(l_args_offset) <- args_offset
  
  # arguments for check_y functions
  args_ckeck_y <- c("y")
  l_args_check_y <- rep(list(expr()), length(args_ckeck_y)) 
  names(l_args_check_y) <- args_ckeck_y
  
  # generic arguments for loss functions
  args_loss_gen <- c("y", 
                     paste(marg1$parameter_names, "1", sep = ""),
                     paste(marg2$parameter_names, "2", sep = ""), 
                     "rho")
  
  # counter required to create parameter specific loss arguments via args_loss_creator()
  ind_arg <- 2
  
  # creates the parameter specific loss arguments
  args_loss_creator <- function(arg_names, ind_arg){    
    arg_names[ind_arg] <- "f"
    l_args <- rep(list(expr()), length(arg_names))
    names(l_args) <- arg_names
    ind_arg <<- ind_arg + 1     
    return(list(l_args,
                arg_names))
  }
  
  
  
  
  ################################################################################
  ############################# generic functions ################################
  
  # rho parameter definition: 
  # 1. rho_simp means simple rho for marginal parameters  
  # 2. rho_gen means rho generic for the f-expression 
  rho_simp <- expr(rho)
  #rho_gen <- expr(tanh(f))
  rho_gen <- expr({  exp( f ) })
  
  # generic functions for the parameter specific creation of loss, risk and gradients 
  
  # Gradient of likelihood w.r.t. first margin
  gradient_marg_gen <- function(cdf1, cdf2, dercdf1){
    
    expr({
      
      F1 <- 1 - !!cdf1
      
      F2 <- !!cdf2
      
      thet <- rho
      
      P11 <- pdffz( F2 - ( ( F1^-thet + F2^-thet - 1 )^( -1/thet ) ) )
      
      P00 <- pdffz(1 - F1 - F2 + P11)
      
      P10 <- pdffz( F1 - P11 )
      
      P01 <- pdffz( F2 - P11 )
      
      
      derp11_dercdf1 <- dvffz( (F1^(-thet) + F2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (F1^((-thet) - 1) * (-thet))) )
      
      ngr <- (
        y[,1]*y[,2]         * ( ( P11 )^(-1)                          *       derp11_dercdf1     ) +
          y[,1]*(1-y[,2])     * ( (P10 )^(-1)              * (1 -  derp11_dercdf1 )   ) +
          (1-y[,1])*y[,2]     * ( (P01 )^(-1)              * -     derp11_dercdf1     ) +
          (1-y[,1])*(1-y[,2]) * ( ( P00 )^(-1) * (-1 + derp11_dercdf1 )   )
      ) * !!dercdf1
      
      return(ngr)
    }
    )
  }
  
  # Gradient of likelihood w.r.t. second margin
  gradient_margin2_gen <- function(cdf1, cdf2, dercdf2){
    
    expr({
      
      F1 <- 1 - !!cdf1
      
      F2 <- !!cdf2
      
      thet <- rho
      
      P11 <- pdffz( F2 - ( ( F1^-thet + F2^-thet - 1 )^( -1/thet ) ) )
      
      P00 <- pdffz(1 - F1 - F2 + P11)
      
      P10 <- pdffz( F1 - P11 )
      
      P01 <- pdffz( F2 - P11 )
      
      
      derp11_dercdf2 <- dvffz( (F1^(-thet) + F2^(-thet) - 1)^((-1/thet) - 1) * ((-1/thet) * (F2^((-thet) - 1) * (-thet))) )
      
      derp11_dercdf2 <- 1 - derp11_dercdf2
      
      ngr <- (
        y[,1]*y[,2]         * ( P11^(-1)                           *       derp11_dercdf2  ) +
          y[,1]*(1-y[,2])     * ( (P10)^(-1)              * -     derp11_dercdf2  ) +
          (1-y[,1])*y[,2]     * ( (P01)^(-1)              * (1 -  derp11_dercdf2 )) +
          (1-y[,1])*(1-y[,2]) * ( (P00)^(-1) * (-1 + derp11_dercdf2 ))
      ) * !!dercdf2
      
      
      return(ngr)
    }
    )
  }
  
  gradient_cop_gen <- function(cdf1, cdf2, rho){
    
    expr({
      
      F1 <- 1 - !!cdf1
      
      F2 <- !!cdf2
      
      thet <- !!rho
      
      
      P11 <- pdffz( F2 - ( ( F1^-thet + F2^-thet - 1 )^( -1/thet ) ) )
      
      
      P00 <- pdffz(1 - F1 - F2 + P11)
      
      P10 <- pdffz( F1 - P11 )
      
      P01 <- pdffz(F2 - P11)
      
      
      derp11_derrho <- dvffz( ((-1 + F1^-thet + F2^-thet)^(-1/thet) *((thet *(F2^thet *log(F1) + F1^thet* log(F2)))/(F2^thet - F1^thet* (-1 + F2^thet)) + log(-1 + F1^-thet + F2^-thet)))/thet^2  )
      
      derp11_derrho <- - derp11_derrho
      
      ngr <-  (
        y[,1]*y[,2]    * ( P11^(-1)                  *   derp11_derrho ) +
          y[,1]*(1-y[,2])     * ( (P10)^(-1) * - derp11_derrho ) +
          (1-y[,1])*y[,2]     * ( (P01)^(-1) * - derp11_derrho ) +
          (1-y[,1])*(1-y[,2]) * ( (P00)^(-1) * derp11_derrho )
      ) * thet
      
      return(ngr)
    }
    )
  }
  
  
  loss_gen <- function(cdf1, cdf2, rho){
    
    expr( {
      
      F1 <- 1 - !!cdf1
      
      F2 <- !!cdf2
      
      thet <- !!rho
      
      P11 <- pdffz( F2 - ( ( F1^-thet + F2^-thet - 1 )^( -1/thet ) ) )
      
      P00 <- pdffz(1 - F1 - F2 + P11)
      
      P10 <- pdffz( F1 - P11 )
      
      P01 <- pdffz(F2 - P11)

      
      return(- ( y[,1]*y[,2] * log( P11 ) +
                   y[,1]*(1-y[,2]) * log( P10 ) +
                   (1-y[,1])*y[,2] * log( P01 ) +
                   (1-y[,1])*(1-y[,2]) * log( P00 ) 
      ))
    } 
    )
  }
  
  
  risk_gen <- function(param){
    
    a <- param
    b <- paste("=", param)
    param <- paste(a, b, collapse = ", ")
    
    loss <- parse_expr(paste0("loss(", param, ")"))  
    
    expr(sum(w*(!!loss)))
  }
  
  
  ################################################################################
  ############### initializing the gamboostLSS Families object ###################
  
  name_fam <- paste("Rotated Clayton (90°) copula", marg1$marg_name, marg2$marg_name)
  # quantile_function <- ...                                                      ### What are suitable quantile functions for copulas? 
  ### not yet required, rather for prediction, I think.
  GaussCopFam <- Families(name = name_fam)
  
  ################################################################################
  ############# creating the Family-objects for the marginals ####################
  
  ### first marginal
  for(i in marg1$parameter_names){
    
    #i <- marg1$parameter_names[1]
    
    # body and arguments for loss
    
    loss_body <- loss_gen(cdf1 = marg1[[i]]$cdf,
                          cdf2 = marg2[["generic"]]$cdf,
                          rho = rho_simp)
    
    loss_args <- args_loss_creator(arg_names = args_loss_gen,
                                   ind_arg = ind_arg)
    
    
    # body and arguments for risk
    
    risk_body <- risk_gen(loss_args[[2]]) 
    risk_args <- l_args_grad_risk
    
    
    # body and arguments for gradient
    
    grad_body <- gradient_marg_gen(cdf1 = marg1[[i]]$cdf, 
                                   cdf2 = marg2[["generic"]]$cdf, 
                                   dercdf1 = marg1[[i]]$dercdf) 
    
    grad_body <- exprs(!!grad_body,
                       ngr <- stabilize_ngradient(ngr, w = w, stabilization),
                       return(ngr))
    
    grad_args <- l_args_grad_risk
    
    
    
    # creating the Family object
    fam_obj <- family_gen(mu1 = mu1, sigma1 = sigma1, nu1 = nu1, tau1 = tau1,
                          mu2 = mu2, sigma2 = sigma2, nu2 = nu2, tau2 = tau2, 
                          rho = rho,
                          stabilization = stabilization,
                          loss_body = loss_body, loss_args = loss_args[[1]], 
                          risk_body = risk_body, risk_args = risk_args,
                          grad_body = grad_body, grad_args = grad_args,
                          offset_body = marg1$offset[[i]], offset_args = l_args_offset,
                          response = marg1$response[[i]],
                          name = marg1$name[[i]],
                          check_y_body = marg1$check_y[[i]], check_y_args = l_args_check_y)
    
    
    # saving the parameter specific family object in overall families object
    GaussCopFam[[paste(i, "1", sep = "")]] <- fam_obj
    
    # removing family object
    rm(fam_obj)
    
  }
  
  
  ### second marginal
  for(i in marg2$parameter_names){
    
    # body and arguments for loss
    loss_body <- loss_gen(cdf1 = marg1[["generic"]]$cdf,
                          cdf2 = marg2[[i]]$cdf,
                          rho = rho_simp)
    
    loss_args <- args_loss_creator(arg_names = args_loss_gen,
                                   ind_arg = ind_arg)
    
    
    # body and arguments for risk
    risk_body <- risk_gen(loss_args[[2]])
    
    risk_args <- l_args_grad_risk
    
    
    # creating the gradient-function
    grad_body <- gradient_margin2_gen(cdf1 = marg1[["generic"]]$cdf,
                                      cdf2 = marg2[[i]]$cdf,
                                      dercdf2 = marg2[[i]]$dercdf)
    
    grad_body <- exprs(!!grad_body,
                       ngr <- stabilize_ngradient(ngr, w = w, stabilization),
                       return(ngr))
    
    grad_args <- l_args_grad_risk
    
    
    
    # creating the family object
    fam_obj <- family_gen(mu1 = mu1, sigma1 = sigma1, nu1 = nu1, tau1 = tau1,
                          mu2 = mu2, sigma2 = sigma2, nu2 = nu2, tau2 = tau2, 
                          rho = rho,
                          stabilization = stabilization,
                          loss_body = loss_body, loss_args = loss_args[[1]], 
                          risk_body = risk_body, risk_args = risk_args,
                          grad_body = grad_body, grad_args = grad_args,
                          offset_body = marg2$offset[[i]], offset_args = l_args_offset,
                          response = marg2$response[[i]],
                          name = marg2$name[[i]],
                          check_y_body = marg2$check_y[[i]], check_y_args = l_args_check_y)
    
    
    # saving the parameter specific family object in overall families object
    GaussCopFam[[paste(i, "2", sep = "")]] <- fam_obj
    
    
    # removing family object
    rm(fam_obj)
    
  }
  
  
  ################################################################################
  ############# creating the Family-objects for the Copula Para ##################
  
  # body and arguments for loss
  loss_body <- loss_gen(cdf1 = marg1[["generic"]]$cdf,
                        cdf2 = marg2[["generic"]]$cdf,
                        rho = rho_gen)
  
  loss_args <- args_loss_creator(arg_names = args_loss_gen,
                                 ind_arg = ind_arg)
  
  
  # body and arguments for risk
  risk_body <- risk_gen(loss_args[[2]])
  
  risk_args <- l_args_grad_risk
  
  
  # body and arguments for gradient
  grad_body <- gradient_cop_gen(cdf1 = marg1[["generic"]]$cdf, 
                                cdf2 = marg2[["generic"]]$cdf, 
                                rho = rho_gen)
  grad_body <- exprs(!!grad_body,
                     ngr <- stabilize_ngradient(ngr, w = w, stabilization),
                     return(ngr))
  
  grad_args <- l_args_grad_risk
  
  
  # definition of offset, response, check_y and name for copula parameter
  offset_cop <- expr({
    if (!is.null(rho)){
      RET <- log(rho) #rho/sqrt(1 - rho^2)
    } else {
      RET <- 0 # rho = 0/(sqrt(1 + 0)) = 0
      # pear_cor <- wdm(x = y[,1], y = y[,2], method = "pearson", weights = w, remove_missing = F) 
      # RET <- pear_cor/sqrt(1-pear_cor^2)
    }
    return(RET)
  })
  
  # offset_cop <- expr({
  #   if (!is.null(rho)){
  #     RET <- rho/sqrt(1 - rho^2)
  #   } else {
  #     if (is.null(rho))
  #     RET <- 0
  #   }
  #   return(RET)
  # })
  
  response_cop <- function(f) exp(f)
  
  check_y_cop <- expr({
    y
  })
  
  name_cop <- "Rotated Clayton (90°) copula for bivariate binary responses: rho(... link)"                                       
  
  
  # creating the family object
  fam_obj <- family_gen(mu1 = mu1, sigma1 = sigma1, nu1 = nu1, tau1 = tau1,
                        mu2 = mu2, sigma2 = sigma2, nu2 = nu2, tau2 = tau2, 
                        rho = rho,
                        stabilization = stabilization,
                        loss_body = loss_body, loss_args = loss_args[[1]], 
                        risk_body = risk_body, risk_args = risk_args,
                        grad_body = grad_body, grad_args = grad_args,
                        offset_body = offset_cop, offset_args = l_args_offset,
                        response = response_cop,
                        name = name_cop,
                        check_y_body = check_y_cop, check_y_args = l_args_check_y)
  
  
  # saving the parameter specific family object in overall families object
  GaussCopFam[["rho"]] <- fam_obj
  
  # removing family object
  rm(fam_obj)
  
  
  ### return final Families object to boost with
  return(GaussCopFam)
  
}

# Gauss_Cop <- Gauss_Cop(marg1 = "LOGLOG", marg2 = "LOGLOG")
# Gauss_Cop$rho@offset


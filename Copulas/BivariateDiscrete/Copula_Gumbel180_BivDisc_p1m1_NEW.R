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
library(gamlss.dist)


# load Marginal functions

source("Marginals/family_gen.R")
#source("Marginals/Discrete/Marginal_Poisson_DMat.R")
# source("Marginals/Discrete/Marginal_Poisson_DMat.R")
# source("Marginals/Discrete/Marginal_Geometric_DMat.R")
# source("Marginals/Discrete/Marginal_ZeroInflatedPoisson_DMat.R")
# source("Marginals/Discrete/Marginal_NegativeBinomial_DMat.R")
# source("Marginals/Discrete/Marginal_ZeroAlteredLogarithmic_DMat.R")
# source("Marginals/Discrete/Marginal_PoissonInverseGaussian_Num.R")
# source("Marginals/Discrete/Marginal_ZeroAlteredNegativeBinomial_DMat.R")
# source("Marginals/Discrete/Marginal_ZeroInflatedNegativeBinomial_DMat.R")

#source("Marginals/Discrete/Marginal_Poisson_DMat.R")
source("Marginals/Discrete/Marginal_Poisson_DMat.R")
source("Marginals/Discrete/Marginal_Geometric_DMat.R")
source("Marginals/Discrete/Marginal_ZeroInflatedPoisson_DMat.R")

#source("Marginals/Discrete/Marginal_NegativeBinomial_DMatTest.R")
source("Marginals/Discrete/Marginal_NegativeBinomial_DMat.R")

#source("Marginals/Discrete/Marginal_ZeroAlteredLogarithmic_DMatTest.R")
source("Marginals/Discrete/Marginal_ZeroAlteredLogarithmic_DMat.R")

source("Marginals/Discrete/Marginal_PoissonInverseGaussian_Num.R")

source("Marginals/Discrete/Marginal_ZeroAlteredNegativeBinomial_DMat.R")
#source("Marginals/Discrete/Marginal_ZeroAlteredNegativeBinomial_DMatTest.R")

source("Marginals/Discrete/Marginal_ZeroInflatedNegativeBinomial_DMat.R")
#source("Marginals/Discrete/Marginal_ZeroInflatedNegativeBinomial_DMatTest.R")

source("Marginals/Discrete/Marginal_ZeroInflatedPoissonInverseGaussian_Num.R")

### Gauss-Copula function (for bivariate binary margins)

Gumbel180_Cop_BivDiscrete <- function(marg1 = NULL, marg2 = NULL,
                                       mu1 = NULL, sigma1 = NULL, nu1 = NULL, tau1 = NULL,
                                       mu2 = NULL, sigma2 = NULL, nu2 = NULL, tau2 = NULL,
                                       rho = NULL,
                                       auxiliary_material = NULL,
                                       stabilization = c("none", "MAD", "L2")){
  
  ################################################################################
  ########################## marginal checks #####################################
  
  if(is.null(marg1)) stop("First marginal distribution not defined.")
  if(is.null(marg2)) stop("Second marginal distribution not defined.")
  
  if(!(marg1 %in% c("POIS", "ZIP", "NBI", "NBII", "GEOM", "ZALG", "ZINBI", "ZANBI", "PIG"))) stop("First Marginal distribution not available.")
  if(!(marg2 %in% c("POIS", "ZIP", "NBI", "NBII", "GEOM", "ZALG", "ZINBI", "ZANBI", "PIG"))) stop("Second Marginal distribution not available.")
  
  ################################################################################
  #################### calling the marginal distributions ########################
  
  if(marg1 == "POIS") {
    marg1 <- Poisson_Mar(loc = 1, offset_mu = mu1, aux_mat = auxiliary_material[[1]], dim_aux_mat = dim(auxiliary_material[[1]]))
  } else if(marg1 == "GEOM") {
    marg1 <- Geometric_Mar(loc = 1, offset_mu = mu1, aux_mat = auxiliary_material[[1]], dim_aux_mat = dim(auxiliary_material[[1]]))
  } else if(marg1 == "ZIP") {
    #marg1 <- ZIPoisson_Mar(loc = 1, offset_mu = mu1, offset_sigma = sigma1, aux_mat = auxiliary_material[[1]]$zeroesIndx, dim_aux_mat = dim(auxiliary_material[[1]][[1]]), aux_matrix = auxiliary_material[[1]][[1]])
    marg1 <- ZIPoisson_Mar(loc = 1, offset_mu = mu1, offset_sigma = sigma1)
  }  else if(marg1 == "NBI") {
    #marg1 <- NBI_Mar(loc = 1, offset_mu = mu1, offset_sigma = sigma1, aux_mat = auxiliary_material[[1]], dim_aux_mat = dim(auxiliary_material[[1]]))
    marg1 <- NBI_Mar(loc = 1, offset_mu = mu1, offset_sigma = sigma1, aux_mat = NULL, dim_aux_mat = NULL)
  } else if(marg1 == "ZALG") {
    marg1 <- ZALG_Mar(loc = 1, offset_mu = mu1, offset_sigma = sigma1, aux_mat = NULL, dim_aux_mat = NULL)
  }  else if(marg1 == "ZINBI") {
    marg1 <- ZINBI_Mar(loc = 1, offset_mu = mu1, offset_sigma = sigma1, offset_nu = nu1, aux_mat = NULL, dim_aux_mat = NULL)
  } else if(marg1 == "ZANBI") {
    marg1 <- ZANBI_Mar(loc = 1, offset_mu = mu1, offset_sigma = sigma1, offset_nu = nu1, aux_mat = NULL, dim_aux_mat = NULL)
  } else if(marg1 == "PIG") {
    marg1 <- PIG_Mar(loc = 1, offset_mu = mu1, offset_sigma = sigma1, aux_mat = NULL, dim_aux_mat = NULL)
  }
  
  
  
  
  if(marg2 == "POIS") {
    marg2 <- Poisson_Mar(loc = 2, offset_mu = mu2, aux_mat = auxiliary_material[[2]], dim_aux_mat = dim(auxiliary_material[[2]])) 
  }else if(marg2 == "GEOM") {
    marg2 <- Geometric_Mar(loc = 2, offset_mu = mu2, aux_mat = auxiliary_material[[2]], dim_aux_mat = dim(auxiliary_material[[2]]))
  }else if(marg2 == "ZIP") {
    #marg2 <- ZIPoisson_Mar(loc = 2, offset_mu = mu2, offset_sigma = sigma2, aux_mat = auxiliary_material[[2]]$zeroesIndx, dim_aux_mat = dim(auxiliary_material[[2]][[1]]), aux_matrix = auxiliary_material[[2]][[1]])
    marg2 <- ZIPoisson_Mar(loc = 2, offset_mu = mu2, offset_sigma = sigma2)
  }else if(marg2 == "NBI") {
    #marg2 <- NBI_Mar(loc = 2, offset_mu = mu2, offset_sigma = sigma2, aux_mat = auxiliary_material[[2]], dim_aux_mat = dim(auxiliary_material[[2]]))
    marg2 <- NBI_Mar(loc = 2, offset_mu = mu2, offset_sigma = sigma2, aux_mat = NULL, dim_aux_mat = NULL)
  } else if(marg2 == "ZALG") {
    marg2 <- ZALG_Mar(loc = 2, offset_mu = mu2, offset_sigma = sigma2, aux_mat = NULL, dim_aux_mat = NULL)
  }  else if(marg2 == "ZINBI") {
    marg2 <- ZINBI_Mar(loc = 2, offset_mu = mu2, offset_sigma = sigma2, offset_nu = nu2, aux_mat = NULL, dim_aux_mat = NULL)
  } else if(marg2 == "ZANBI") {
    marg2 <- ZANBI_Mar(loc = 2, offset_mu = mu2, offset_sigma = sigma2, offset_nu = nu2, aux_mat = NULL, dim_aux_mat = NULL)
  } else if(marg2 == "PIG") {
    marg2 <- PIG_Mar(loc = 2, offset_mu = mu2, offset_sigma = sigma2, aux_mat = NULL, dim_aux_mat = NULL)
  }
  
  
  
  ################################################################################
  ##################### check offsets of copula parameter ########################
  
  
  if ((!is.null(rho) && rho <= 1))
    stop(sQuote("rho"), " must be greater than 1.")
  
  
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
  rho_gen <- expr({  ( exp( f ) + 1 ) })
  
  # generic functions for the parameter specific creation of loss, risk and gradients 
  
  # Gradient of likelihood w.r.t. first margin
  gradient_marg_gen <- function(cdf1, cdf2, pdf1, pdf2, dercdf1, derpdf1){
    
    expr({
      
      F1 <- pdffz(!!cdf1)
      
      F2 <- pdffz(!!cdf2)
      
      CDF1m1 <- pdffz(F1 - (!!pdf1))
      
      CDF2m1 <- pdffz(F2 - (!!pdf2))
      
      thet <- rho
      
      ###############################################################################################################################################
      FTilde_1 <- pdffz(1 - F1)
      
      FTilde_2 <- pdffz(1 - F2)
      
      CDF1_Tildem1 <- pdffz( 1 - (F1 - (!!pdf1)) )
      
      CDF2_Tildem1 <- pdffz( 1 - (F2 - (!!pdf2)) )
      
      ###############################################################################################################################################
      T1 <- pdffz( F1 + F2 - 1 + ( exp( - ( (-log(FTilde_1))^thet + (-log(FTilde_2))^thet )^(1/thet) ) ) )
      
      T2 <- pdffz( CDF1m1 + F2 - 1 + ( exp( - ( (-log(CDF1_Tildem1))^thet + (-log(FTilde_2))^thet )^(1/thet) ) )  ) 
      
      T3 <- pdffz( F1 + CDF2m1 - 1 + ( exp( - ( (-log(FTilde_1))^thet + (-log(CDF2_Tildem1))^thet )^(1/thet) ) )  ) 
      
      T4 <- pdffz( CDF1m1 + CDF2m1 - 1 + ( exp( - ( (-log(CDF1_Tildem1))^thet + (-log(CDF2_Tildem1))^thet )^(1/thet) ) )  ) 
      ###############################################################################################################################################
      
      L_Total <- pdffz( T1 - T2 - T3 + T4 )
      
      # Partial derivatives:
      ######################## ########################
      dF1_deta <- dvffz(!!dercdf1)
      
      dpdf1_deta <- dvffz(!!derpdf1)
      
      
      dT1_dF1 <- dvffz( 1 - ( (exp(-((-log(FTilde_1))^thet + (-log(FTilde_2))^thet)^((1/thet)))* (-log(FTilde_1))^(-1 + thet)* ((-log(FTilde_1))^thet + (-log(FTilde_2))^thet)^(-1 + 1/thet))/FTilde_1 ) )
      
      dT2_dF1 <- dvffz( 1 - ( (exp(-((-log(CDF1_Tildem1))^thet + (-log(FTilde_2))^thet)^((1/thet)))* (-log(CDF1_Tildem1))^(-1 + thet)* ((-log(CDF1_Tildem1))^thet + (-log(FTilde_2))^thet)^(-1 + 1/thet))/CDF1_Tildem1 ) )
      
      dT3_dF1 <- dvffz( 1 - ( (exp(-((-log(FTilde_1))^thet + (-log(CDF2_Tildem1))^thet)^((1/thet)))* (-log(FTilde_1))^(-1 + thet)* ((-log(FTilde_1))^thet + (-log(CDF2_Tildem1))^thet)^(-1 + 1/thet))/FTilde_1 ) )
      
      dT4_dF1 <- dvffz( 1 - ( (exp(-((-log(CDF1_Tildem1))^thet + (-log(CDF2_Tildem1))^thet)^((1/thet)))* (-log(CDF1_Tildem1))^(-1 + thet)* ((-log(CDF1_Tildem1))^thet + (-log(CDF2_Tildem1))^thet)^(-1 + 1/thet))/CDF1_Tildem1 ) )
      
      ngr <- ( ( L_Total )^(-1) * ( ( dT1_dF1 - dT3_dF1 ) * (dF1_deta) + ( - dT2_dF1 + dT4_dF1  ) * (dF1_deta - dpdf1_deta) )  )
      
      return(ngr)
      
    })
  }
  
  # Gradient of likelihood w.r.t. second margin
  gradient_margin2_gen <- function(cdf1, cdf2, pdf1, pdf2, dercdf2, derpdf2){
    
    expr({
      
      F1 <- pdffz(!!cdf1)
      
      F2 <- pdffz(!!cdf2)
      
      CDF1m1 <- pdffz(F1 - (!!pdf1))
      
      CDF2m1 <- pdffz(F2 - (!!pdf2))
      
      thet <- rho
      
      ###############################################################################################################################################
      FTilde_1 <- pdffz(1 - F1)
      
      FTilde_2 <- pdffz(1 - F2)
      
      CDF1_Tildem1 <- pdffz( 1 - (F1 - (!!pdf1)) )
      
      CDF2_Tildem1 <- pdffz( 1 - (F2 - (!!pdf2)) )
      
      ###############################################################################################################################################
      T1 <- pdffz( F1 + F2 - 1 + ( exp( - ( (-log(FTilde_1))^thet + (-log(FTilde_2))^thet )^(1/thet) ) ) )
      
      T2 <- pdffz( CDF1m1 + F2 - 1 + ( exp( - ( (-log(CDF1_Tildem1))^thet + (-log(FTilde_2))^thet )^(1/thet) ) )  ) 
      
      T3 <- pdffz( F1 + CDF2m1 - 1 + ( exp( - ( (-log(FTilde_1))^thet + (-log(CDF2_Tildem1))^thet )^(1/thet) ) )  ) 
      
      T4 <- pdffz( CDF1m1 + CDF2m1 - 1 + ( exp( - ( (-log(CDF1_Tildem1))^thet + (-log(CDF2_Tildem1))^thet )^(1/thet) ) )  ) 
      ###############################################################################################################################################
      
      L_Total <- pdffz( T1 - T2 - T3 + T4 )
      
      
      # Partial derivatives:
      ######################## ########################
      dF2_deta <- dvffz(!!dercdf2)
      
      dpdf2_deta <- dvffz(!!derpdf2)
      
      dT1_dF2 <- dvffz( 1 - ( (exp(-((-log(FTilde_2))^thet + (-log(FTilde_1))^thet)^((1/thet)))* (-log(FTilde_2))^(-1 + thet)* ((-log(FTilde_2))^thet + (-log(FTilde_1))^thet)^(-1 + 1/thet))/FTilde_2 ) )
      
      dT2_dF2 <- dvffz( 1 - ( (exp(-((-log(FTilde_2))^thet + (-log(CDF1_Tildem1))^thet)^((1/thet)))* (-log(FTilde_2))^(-1 + thet)* ((-log(FTilde_2))^thet + (-log(CDF1_Tildem1))^thet)^(-1 + 1/thet))/FTilde_2 ) )
      
      dT3_dF2 <- dvffz( 1 - ( (exp(-((-log(CDF2_Tildem1))^thet + (-log(FTilde_1))^thet)^((1/thet)))* (-log(CDF2_Tildem1))^(-1 + thet)* ((-log(CDF2_Tildem1))^thet + (-log(FTilde_1))^thet)^(-1 + 1/thet))/CDF2_Tildem1 ) )
      
      dT4_dF2 <- dvffz( 1 - ( (exp(-((-log(CDF2_Tildem1))^thet + (-log(CDF1_Tildem1))^thet)^((1/thet)))* (-log(CDF2_Tildem1))^(-1 + thet)* ((-log(CDF2_Tildem1))^thet + (-log(CDF1_Tildem1))^thet)^(-1 + 1/thet))/CDF2_Tildem1 ) )
      
      ######################## ########################
      ngr <- ( ( L_Total )^(-1) * ( ( dT1_dF2 - dT2_dF2 ) * (dF2_deta) + ( -dT3_dF2 + dT4_dF2 ) * (dF2_deta - dpdf2_deta) ) )
      
      return(ngr)
      
    })
  }
  
  gradient_cop_gen <- function(cdf1, cdf2, pdf1, pdf2, rho){
    
    expr({
      
      F1 <- pdffz(!!cdf1)
      
      F2 <- pdffz(!!cdf2)
      
      CDF1m1 <- pdffz(F1 - (!!pdf1))
      
      CDF2m1 <- pdffz(F2 - (!!pdf2))
      
      thet <- !!rho
      
      
      ###############################################################################################################################################
      FTilde_1 <- pdffz(1 - F1)
      
      FTilde_2 <- pdffz(1 - F2)
      
      CDF1_Tildem1 <- pdffz( 1 - (F1 - (!!pdf1)) )
      
      CDF2_Tildem1 <- pdffz( 1 - (F2 - (!!pdf2)) )
      
      ###############################################################################################################################################
      T1 <- pdffz( F1 + F2 - 1 + ( exp( - ( (-log(FTilde_1))^thet + (-log(FTilde_2))^thet )^(1/thet) ) ) )
      
      T2 <- pdffz( CDF1m1 + F2 - 1 + ( exp( - ( (-log(CDF1_Tildem1))^thet + (-log(FTilde_2))^thet )^(1/thet) ) )  ) 
      
      T3 <- pdffz( F1 + CDF2m1 - 1 + ( exp( - ( (-log(FTilde_1))^thet + (-log(CDF2_Tildem1))^thet )^(1/thet) ) )  ) 
      
      T4 <- pdffz( CDF1m1 + CDF2m1 - 1 + ( exp( - ( (-log(CDF1_Tildem1))^thet + (-log(CDF2_Tildem1))^thet )^(1/thet) ) )  ) 
      ###############################################################################################################################################
      
      L_Total <- pdffz( T1 - T2 - T3 + T4 )
      
      ########## Partial derivatives of copula terms. 
      dT1_dtheta <- dvffz( (1/(thet^2))*exp(-((-log(FTilde_1))^thet + (-log(FTilde_2))^thet)^((1/thet)))* ((-log(FTilde_1))^thet + (-log(FTilde_2))^thet)^(-1 + 1/thet)* (-thet* (-log(FTilde_1))^thet* log(-log(FTilde_1)) + ((-log(FTilde_1))^thet + (-log(FTilde_2))^thet)* log((-log(FTilde_1))^thet + (-log(FTilde_2))^thet) - thet *(-log(FTilde_2))^thet* log(-log(FTilde_2))) )
      
      dT2_dtheta <- dvffz( (1/(thet^2))*exp(-((-log(CDF1_Tildem1))^thet + (-log(FTilde_2))^thet)^((1/thet)))* ((-log(CDF1_Tildem1))^thet + (-log(FTilde_2))^thet)^(-1 + 1/thet)* (-thet* (-log(CDF1_Tildem1))^thet* log(-log(CDF1_Tildem1)) + ((-log(CDF1_Tildem1))^thet + (-log(FTilde_2))^thet)* log((-log(CDF1_Tildem1))^thet + (-log(FTilde_2))^thet) - thet *(-log(FTilde_2))^thet* log(-log(FTilde_2))) )
      
      dT3_dtheta <- dvffz( (1/(thet^2))*exp(-((-log(FTilde_1))^thet + (-log(CDF2_Tildem1))^thet)^((1/thet)))* ((-log(FTilde_1))^thet + (-log(CDF2_Tildem1))^thet)^(-1 + 1/thet)* (-thet* (-log(FTilde_1))^thet* log(-log(FTilde_1)) + ((-log(FTilde_1))^thet + (-log(CDF2_Tildem1))^thet)* log((-log(FTilde_1))^thet + (-log(CDF2_Tildem1))^thet) - thet *(-log(CDF2_Tildem1))^thet* log(-log(CDF2_Tildem1))) )
      
      dT4_dtheta <- dvffz( (1/(thet^2))*exp(-((-log(CDF1_Tildem1))^thet + (-log(CDF2_Tildem1))^thet)^((1/thet)))* ((-log(CDF1_Tildem1))^thet + (-log(CDF2_Tildem1))^thet)^(-1 + 1/thet)* (-thet* (-log(CDF1_Tildem1))^thet* log(-log(CDF1_Tildem1)) + ((-log(CDF1_Tildem1))^thet + (-log(CDF2_Tildem1))^thet)* log((-log(CDF1_Tildem1))^thet + (-log(CDF2_Tildem1))^thet) - thet *(-log(CDF2_Tildem1))^thet* log(-log(CDF2_Tildem1))) )
      
      
      ########## FINAL GRADIENT
      ngr <-  ( ( L_Total )^(-1) * ( dT1_dtheta - dT2_dtheta - dT3_dtheta + dT4_dtheta ) * (exp(f))
      )
      
      return(ngr)
    }
    )
  }
  
  loss_gen <- function(cdf1, cdf2, pdf1, pdf2, rho){
    
    expr({
      
      F1 <- pdffz(!!cdf1)
      
      F2 <- pdffz(!!cdf2)
      
      CDF1m1 <- pdffz(F1 - (!!pdf1))
      
      CDF2m1 <- pdffz(F2 - (!!pdf2))
      
      ###############################################################################################################################################
      FTilde_1 <- pdffz(1 - F1)
      
      FTilde_2 <- pdffz(1 - F2)
      
      CDF1_Tildem1 <- pdffz( 1 - (F1 - (!!pdf1)) )
      
      CDF2_Tildem1 <- pdffz( 1 - (F2 - (!!pdf2)) )
      
      ###############################################################################################################################################
      thet <- !!rho
      
      T1 <- pdffz( F1 + F2 - 1 + ( exp( - ( (-log(FTilde_1))^thet + (-log(FTilde_2))^thet )^(1/thet) ) ) )
      
      T2 <- pdffz( CDF1m1 + F2 - 1 + ( exp( - ( (-log(CDF1_Tildem1))^thet + (-log(FTilde_2))^thet )^(1/thet) ) )  ) 
      
      T3 <- pdffz( F1 + CDF2m1 - 1 + ( exp( - ( (-log(FTilde_1))^thet + (-log(CDF2_Tildem1))^thet )^(1/thet) ) )  ) 
      
      T4 <- pdffz( CDF1m1 + CDF2m1 - 1 + ( exp( - ( (-log(CDF1_Tildem1))^thet + (-log(CDF2_Tildem1))^thet )^(1/thet) ) )  ) 
      ###############################################################################################################################################
      
      L_Total <- pdffz( T1 - T2 - T3 + T4 )
      
      return(- (  log( L_Total ) ) ) 
      
    })
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
  
  name_fam <- paste("RotatedGumbelCop180", marg1$marg_name, marg2$marg_name)
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
                          pdf1 = marg1[[i]]$pdf,
                          pdf2 = marg2[["generic"]]$pdf,
                          rho = rho_simp)
    
    loss_args <- args_loss_creator(arg_names = args_loss_gen,
                                   ind_arg = ind_arg)
    
    
    # body and arguments for risk
    
    risk_body <- risk_gen(loss_args[[2]]) 
    risk_args <- l_args_grad_risk
    
    
    # body and arguments for gradient
    
    grad_body <- gradient_marg_gen(cdf1 = marg1[[i]]$cdf, 
                                   cdf2 = marg2[["generic"]]$cdf, 
                                   pdf1 = marg1[[i]]$pdf,
                                   pdf2 = marg2[["generic"]]$pdf,
                                   dercdf1 = marg1[[i]]$dercdf, 
                                   derpdf1 = marg1[[i]]$derpdf) 
    
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
                          pdf1 = marg1[["generic"]]$pdf,
                          pdf2 = marg2[[i]]$pdf,
                          rho = rho_simp)
    
    loss_args <- args_loss_creator(arg_names = args_loss_gen,
                                   ind_arg = ind_arg)
    
    
    # body and arguments for risk
    risk_body <- risk_gen(loss_args[[2]])
    
    risk_args <- l_args_grad_risk
    
    
    # creating the gradient-function
    grad_body <- gradient_margin2_gen(cdf1 = marg1[["generic"]]$cdf,
                                      cdf2 = marg2[[i]]$cdf,
                                      pdf1 = marg1[["generic"]]$pdf,
                                      pdf2 = marg2[[i]]$pdf,
                                      dercdf2 = marg2[[i]]$dercdf,
                                      derpdf2 = marg2[[i]]$derpdf)
    
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
                        pdf1 = marg1[["generic"]]$pdf,
                        pdf2 = marg2[["generic"]]$pdf,
                        rho = rho_gen)
  
  loss_args <- args_loss_creator(arg_names = args_loss_gen,
                                 ind_arg = ind_arg)
  
  
  # body and arguments for risk
  risk_body <- risk_gen(loss_args[[2]])
  
  risk_args <- l_args_grad_risk
  
  
  # body and arguments for gradient
  grad_body <- gradient_cop_gen(cdf1 = marg1[["generic"]]$cdf, 
                                cdf2 = marg2[["generic"]]$cdf, 
                                pdf1 = marg1[["generic"]]$pdf,
                                pdf2 = marg2[["generic"]]$pdf,
                                rho = rho_gen)
  grad_body <- exprs(!!grad_body,
                     ngr <- stabilize_ngradient(ngr, w = w, stabilization),
                     return(ngr))
  
  grad_args <- l_args_grad_risk
  
  
  # definition of offset, response, check_y and name for copula parameter
  offset_cop <- expr({
    if (!is.null(rho)){
      RET <- log(rho - 1)
    } else {
      RET <- 0 # rho = 0/(sqrt(1 + 0)) = 0
      #RET <- rep(0, nrow(y))
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
  
  response_cop <- function(f) (exp(f) + 1)
  
  check_y_cop <- expr({
    y
  })
  
  name_cop <- "Gumbel copula (180°) for bivariate discrete / count responses: rho(... link)"                                       
  
  
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



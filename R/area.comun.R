area.comun <-
function(x,y,model="Normal-HalfNormal",method="Kernel"){
    
    if(model=="Normal-Exponencial"){
      
      fit<-lm(y~x)
      
      error<-fit$residuals
      
      densi.exp<-function(error,par){(1/par[1])*exp((error/par[1])+((par[2]^2)/2*par[1]^2))*pnorm((-error/par[2])-(par[2]/par[1]))}
      
      if(method=="Adaptative-Kernel"){
        area_comun_exp <- function (error) { 
          F <- function (x,par) integrate (function(t) densi.exp(t,par), -Inf, x) $ value
          u.trans <- function(par) pmin (sapply (error, function(x) F(x,par)),
                                         0.9999999) 
          AC <- function(par){
            alfa<-par[3]
            resid_trans    <- qnorm(u.trans(par))
            densidad_resid <- adkernel(resid_trans,alpha=alfa) 
            normal         <- dnorm(densidad_resid$x)
            fmin           <- pmin(normal,densidad_resid$y)
            long    <- length(densidad_resid$y)
            bases   <- diff(densidad_resid$x)
            alturas <- (fmin[-1]+fmin[-long])/2
            -sum(bases*alturas)
          }
          area <- optim (c(1,1,0.5), AC,method='L-BFGS-B',lower=0,upper = 1)
          sigma_u<-area$par[1]
          sigma_v<-area$par[2]
          landa_ac<-sigma_u/sigma_v
          list(sigma_u=sigma_u,sigma_v=sigma_v,landa_ac=landa_ac)
        }
        res<-area_comun_exp(error)
        return(res)
      }
      if(method=="Kernel"){
        
        area_comun_exp <- function (error) { 
          F <- function (x,par) integrate (function(t) densi.exp(t,par), -Inf, x) $ value
          u.trans <- function(par) pmin (sapply (error, function(x) F(x,par)),
                                         0.9999999) 
          AC <- function(par){
            resid_trans    <- qnorm(u.trans(par))
            densidad_resid <- density(resid_trans) 
            normal         <- dnorm(densidad_resid$x)
            fmin           <- pmin(normal,densidad_resid$y)
            long    <- length(densidad_resid$y)
            bases   <- diff(densidad_resid$x)
            alturas <- (fmin[-1]+fmin[-long])/2
            -sum(bases*alturas)
          }
          area <- optim (c(1,1), AC)
          sigma_u<-area$par[1]
          sigma_v<-area$par[2]
          landa_ac<-sigma_u/sigma_v
          resid_trans<-qnorm(u.trans(area$par))
          list(landa_ac=landa_ac,resid_trans=resid_trans)
        }
        res<-area_comun_exp(error)
        return(res)
      }
      
    }
    if(model=="Normal-HalfNormal"){
      
      modelo_mv <- suppressWarnings(frontier::sfa(y ~ x))
      e.mv <- modelo_mv$resid
      asimetria_mv <- skewness(e.mv)
      param_mv <- modelo_mv$mleParam
      sigma2_mv <- param_mv[["sigmaSq"]]
      gamma_mv <- param_mv[["gamma"]]
      sigma2u_mv <- gamma_mv * sigma2_mv
      sigmau_mv <- sqrt(sigma2u_mv)
      sigma2v_mv <- (1 - gamma_mv) * sigma2_mv
      sigmav_mv <- sqrt(sigma2v_mv)
      sigma2epsi_mv <- sigma2v_mv + sigma2u_mv * ((pi - 2)/pi)
      landa_mv <- sigmau_mv/sigmav_mv
      densidad.epsi<-function (x, L, s2) 
      {
        2/sqrt(s2 * (1 + (2/(pi/L^2 + (pi - 2))))) * dnorm(x/sqrt(s2 * 
                                                                    (1 + (2/(pi/L^2 + (pi - 2)))))) * pnorm(-x * L/sqrt(s2 * 
                                                                                                                          (1 + (2/(pi/L^2 + (pi - 2))))))
      }
      if(method=="Adaptative-Kernel"){
        area.comun_norm <- function (e.mv, s2) { 
          F <- function (x,L,s2) integrate (function(t) densidad.epsi(t,L,s2), -Inf, x) $ value
          u.trans <- function(L) pmin (sapply (e.mv, function(x) F(x,L,s2)),
                                       0.9999999) 
          AC <- function(par){
            L<-par[1]
            alfa<-par[2]
            
            resid_trans    <- qnorm(u.trans(L))
            densidad_resid <- adkernel(resid_trans,alpha=alfa)
            normal         <- dnorm(densidad_resid$x)
            fmin           <- pmin(normal,densidad_resid$y) 
            long    <- length(densidad_resid$y) 
            bases   <- diff(densidad_resid$x)
            alturas <- (fmin[-1]+fmin[-long])/2
            -sum(bases*alturas)
          }
          area <- optim (c(1,0.5),AC,method='L-BFGS-B',lower=0,upper = 1)
          landa_ac <- area$par[1]
          asimetria_ac <- (sqrt(2)*(pi-4)*landa_ac^3)/((pi+(pi-2)*landa_ac^2)^(3/2))  #Asimetria teorica dado un landa
          resid_trans    <- qnorm(u.trans(landa_ac))
          list(landa_ac=landa_ac,asimetria_ac=asimetria_ac,resid_trans=resid_trans)
        }
      }
      if(method=="Kernel"){
        area.comun_norm <- function (e.mv, s2) { 
          F <- function (x,L,s2) integrate (function(t) densidad.epsi(t,L,s2), -Inf, x) $ value
          u.trans <- function(L) pmin (sapply (e.mv, function(x) F(x,L,s2)),
                                       0.9999999) 
          AC <- function(L){
            resid_trans    <- qnorm(u.trans(L))
            densidad_resid <- density(resid_trans) 
            normal         <- dnorm(densidad_resid$x)
            fmin           <- pmin(normal,densidad_resid$y) 
            long    <- length(densidad_resid$y) 
            bases   <- diff(densidad_resid$x)
            alturas <- (fmin[-1]+fmin[-long])/2
            sum(bases*alturas)
          }
          area <- optimize (AC, c(0,20), tol=0.0001, maximum=TRUE)
          landa_ac <- area[[1]]
          asimetria_ac <- (sqrt(2)*(pi-4)*landa_ac^3)/((pi+(pi-2)*landa_ac^2)^(3/2))  #Asimetria teorica dado un landa
          resid_trans    <- qnorm(u.trans(landa_ac))
          list(landa_ac=landa_ac,asimetria_ac=asimetria_ac,resid_trans=resid_trans)
        }
      }
      
      res<-area.comun_norm(e.mv,sigma2epsi_mv)
      return(res)
    }
  }

ng.area.comun <-
function(x,y){
    
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
    negentropy<-function(x){(1/12)*(mean(x^3)^2)+(1/48)*kurtosis(x)^2}
    
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
    
    landas<-NULL
    ng<-NULL
    resid<-e.mv
    asim1<-e1071::skewness(resid)
    asim<-e1071::skewness(resid)
    for(i in 1:(ceiling(length(resid)*0.05)+1)){
      if(asim>0){
        area<-area.comun_norm(resid,sigma2epsi_mv)
        landas[i]<-area$landa_ac
        ng[i]<-negentropy(area$resid_trans)
        ng_all<-sapply(1:length(area$resid_trans),function(j)negentropy(area$resid_trans[-j]))
        quitar<-which.min(ng_all)
        resid<-resid[-quitar]
        asim<-e1071::skewness(resid)
      }
    }
    
    if(asim1<0){
      area<-area.comun_norm(resid,sigma2epsi_mv)
      return(area)}else{
        return(c(landa_ng=landas[which.min(ng)],landa_ac=landas[1]))
      }
  }

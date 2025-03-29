SUBROUTINE COMPUTE_INDX(lam,spec,indice)

  !The absorption-line index calculations, Worthey (1994) and Trager et al. (1998)
  USE params; USE sub_func                       
  IMPLICIT NONE
  
  INTEGER :: i, lmin, lmax
  REAL(KIND(1.d0)) :: s_low, s_high, flux_b, flux_r, spec_b, spec_r
  REAL(KIND(1.d0)) :: flux_count
  REAL(KIND(1.d0)) :: lambda_b, lambda_r
  REAL(KIND(1.d0)), INTENT(IN), DIMENSION(nspec) :: lam, spec
  REAL(KIND(1.d0)), INTENT(inout), DIMENSION(n_index) :: indice
  
  
  !---------------------------------
  DO i=1, n_index


     
     !Blue continuum bandpass
     
     lmin = locate(lam,index_info(i,3))
     lmax = locate(lam,index_info(i,4))
     
     s_low = spec(lmin)+(spec(lmin+1)-spec(lmin))*&
                (index_info(i,3)-lam(lmin))/(lam(lmin+1)-lam(lmin))
                
     s_high = spec(lmax)+(spec(lmax+1)-spec(lmax))*&
                (index_info(i,4)-lam(lmax))/(lam(lmax+1)-lam(lmax))
            
     IF(lmin.EQ.lmax)THEN
      
         flux_b = (s_low+s_high)/2.*(index_info(i,4)-index_info(i,3))
     ELSE 
           
        flux_b = SUM( ABS(lam(lmin+2:lmax)-lam(lmin+1:lmax-1))* &
                           (spec(lmin+2:lmax) + spec(lmin+1:lmax-1))/2.0)   
                                                
                
        flux_b = flux_b + ((index_info(i,4)-lam(lmax))*(s_high+spec(lmax))/2.0)       
        flux_b = flux_b + ((lam(lmin+1)-index_info(i,3))*(s_low+spec(lmin+1))/2.0)
         
     ENDIF
      spec_b = flux_b/(index_info(i,4)-index_info(i,3))
      lambda_b = (index_info(i,4)+index_info(i,3))/2.0
      
   
      !----------------------------------------------------
      !Red continuum bandpass
      
      lmin = locate(lam,index_info(i,5))
      lmax = locate(lam,index_info(i,6))
           
      s_low = spec(lmin)+(spec(lmin+1)-spec(lmin))*&
                (index_info(i,5)-lam(lmin))/(lam(lmin+1)-lam(lmin))
      s_high = spec(lmax)+(spec(lmax+1)-spec(lmax))*&
                (index_info(i,6)-lam(lmax))/(lam(lmax+1)-lam(lmax))          
        
     IF(lmin.EQ.lmax)THEN
      
         flux_r = (s_low+s_high)/2.*(index_info(i,6)-index_info(i,5))
     ELSE          
                
        flux_r = SUM( ABS(lam(lmin+2:lmax)-lam(lmin+1:lmax-1))* &
                (spec(lmin+2:lmax) + spec(lmin+1:lmax-1))/2.0)     
                
                           
         
       flux_r = flux_r + ((index_info(i,6)-lam(lmax))*(s_high+spec(lmax))/2.0) 
       flux_r = flux_r + ((lam(lmin+1)-index_info(i,5))*(s_low+spec(lmin+1))/2.0)  
      
     ENDIF 
      spec_r = flux_r/(index_info(i,6)-index_info(i,5))
      lambda_r = (index_info(i,6)+index_info(i,5))/2.0
      
      !--------------------------------------------
      
      lmin = locate(lam,index_info(i,1))
      lmax = locate(lam,index_info(i,2))
      
      s_low = spec(lmin)/(spec_b+(spec_r-spec_b)*(lam(lmin)-lambda_b)/(lambda_r-lambda_b))&
                  +(spec(lmin+1)/(spec_b+(spec_r-spec_b)*(lam(lmin+1)-lambda_b)/(lambda_r-lambda_b))&
                  -spec(lmin)/(spec_b+(spec_r-spec_b)*(lam(lmin)-lambda_b)/(lambda_r-lambda_b)))*&
                  (index_info(i,1)-lam(lmin))/(lam(lmin+1)-lam(lmin))
                  
      s_high = spec(lmax)/(spec_b+(spec_r-spec_b)*(lam(lmax)-lambda_b)/(lambda_r-lambda_b))&
                   +(spec(lmax+1)/(spec_b+(spec_r-spec_b)*(lam(lmax+1)-lambda_b)/(lambda_r-lambda_b))&
                   -spec(lmax)/(spec_b+(spec_r-spec_b)*(lam(lmax)-lambda_b)/(lambda_r-lambda_b)))*&
                   (index_info(i,2)-lam(lmax))/(lam(lmax+1)-lam(lmax)) 
                   
                   
        
    IF(lmin.EQ.lmax)THEN
      
       indice(i) = (s_low+s_high)/2.*(index_info(i,2)-index_info(i,1))
    ELSE           
              
      indice(i) = SUM( ABS(lam(lmin+2:lmax)-lambda(lmin+1:lmax-1)) * &
                     ((spec(lmin+2:lmax)/(spec_b+(spec_r-spec_b)*(lam(lmin+2:lmax)-lambda_b)/(lambda_r-lambda_b))&
                      + spec(lmin+1:lmax-1)/(spec_b+(spec_r-spec_b)*(lam(lmin+1:lmax-1)-lambda_b)/(lambda_r-lambda_b)))/2.0))
                    
                                    
                
      
                          
      indice(i) = indice(i) +  ((index_info(i,2)-lam(lmax))*&
                      (s_high+spec(lmax)/(spec_b+(spec_r-spec_b)*(lam(lmax)-lambda_b)/(lambda_r-lambda_b)))/2.0)
                           
                           
      indice(i) = indice(i) + ((lam(lmin+1)-index_info(i,1))*&
                      (s_low+spec(lmin+1)/(spec_b+(spec_r-spec_b)*(lam(lmin+1)-lambda_b)/(lambda_r-lambda_b)))/2.0)               
    
    ENDIF   
      
       IF (index_info(i,7).EQ.1.) THEN
        !for magnitude units
        indice(i) = -2.5*LOG10(indice(i)/(index_info(i,2)-index_info(i,1)))
      ELSE IF (index_info(i,7).EQ.2.) THEN
        !for Angstrom units
        indice(i) = (index_info(i,2)-index_info(i,1)) - indice(i)
        
      ENDIF 
      
     

  ENDDO




END SUBROUTINE COMPUTE_INDX

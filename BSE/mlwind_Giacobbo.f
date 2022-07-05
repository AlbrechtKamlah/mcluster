***
      real*8 FUNCTION mlwind(kw,lum,r,mt,mc,rl,z)
      implicit none
      integer kw! Albrecht 12.08.2020 - MOCCA-BSE (LevelD)
      integer mdflag ! mdflag Nbody6++GPU parameter - Albrecht 23.07.2020
      common /FLAGS2/ mdflag  ! include mdflag - Albrecht 23.07.2020
      real*8 lum,r,mt,mc,rl,z
      real*8 fz ! include fz for proper 
      real*8 dml,dms,dmt,p0,x,mew,neta,bwind,hewind,mxns
      real*8 lum0,kap,V1,V2,Zsun,eddington,gsun,flbv 
      parameter(lum0=7.0d+04,kap=-0.5d0,V1=1.3d0,V2=2.6d0,Zsun=0.02,
     &          gsun=1.0d0,flbv = 1.5d0) ! Albrecht 12.08.2020 - MOCCA-BSE (LevelD)
      common /value1/ neta,bwind,hewind,mxns ! Albrecht 12.08.2020 - MOCCA-BSE
      real*8 teff,teff1,t40,vw ! include teff,teff1,t40,vw,flbv - Albrecht 23.07.2020
      real*8 ind1,ind2,ind3,ind4,ind5 ! Albrecht 12.08.2020 - MOCCA-BSE (LevelD)
      real*8 ind,Xh,edfac,teffR,geffR ! Albrecht 12.08.2020 - MOCCA-BSE (LevelD)
*
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                                                          
* mdflag - Calculate stellar wind mass loss          y                          *                                
* Level A        = 1 - SSE basic rates [Hurley et al. (2000)]                  *
* Level A        = 2 - SSE + LBV added [Humphreys & Davidson (1994)]           *
* Level B        = 3 - [Belczynski et al. (2010)].                             *
* Level C        = 4 - [Belczynski et al. (2010)] no bi-stability jump.        *
* Level D        = 5 - [Giacobbo et al. (2018)]                                *                                                                                                                                     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC          
CCCCCCC ADAPTED FROM MOCCA-BSE BY ALBRECHT - 12.08.2020 - SEE BELOW 
      teff1 = 0.0d0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate stellar wind mass loss.                                            *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(mdflag.le.4)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Use fixed ALPHAS [Giacobbo et al. (2018), Eq. 7] in the prescriptions        *  
* by [Belczynski et al. (2010)]                                                *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        ind1 = 0.85d0
        ind2 = 0.85d0
        ind3 = 0.0d0
        ind4 = 0.86d0
        ind5 = 0.5d0
      else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate Hydrogen fraction                                                  *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        if(kw.ge.7 .and.kw.lt.10)then
         Xh = 0.d0
        else
         if(kw.eq.0 .or.kw.eq.1)then
            Xh = 0.7d0
         elseif(kw.eq.2)then
            Xh = 0.6d0
         elseif(kw.eq.3)then
            Xh = 0.5d0
         elseif(kw.eq.4)then
            Xh = 0.4d0
         elseif(kw.eq.5.or.kw.eq.6)then
            Xh = 0.2d0
         endif
       endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate Eddington factor based on Eq. 8 in [Graefener et al. (2011)]       *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
        edfac = 10.d0**(-4.813d0 + Log10(1.d0 + Xh) + Log10(lum) -
     &               Log10(mt))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate exponent of the dependence on metallicity: (Z/Zsun)^ind,           *
* based on the Eddington factor [Chen et al. (2015)]                           *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        ind = 0.85d0
        if(edfac.gt.2.d0/3.d0)then
          ind = (2.45d0 - 2.4d0*edfac)
        endif
*
        ind1 = ind
        ind2 = ind
        ind3 = ind
        ind4 = ind
        ind5 = ind
      endif
*
CCCCCCC ADAPTED FROM MOCCA-BSE ALBRECHT - 12.08.2020 - SEE ABOVE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
*
      dms = 0.d0
*
      teff1 = 0.0d0
      teff = 3.762d0 + 0.25d0*log10(lum) - 0.5d0*log10(r)
      teff = 10.d0**teff
      teff = MIN(teff,50000.d0)
      x = 1.0d-5*r*SQRT(lum) ! ADDED BY ALBRECHT - 12.08.2020
*
      if(lum.gt.4000.d0.and.mdflag.le.2)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Apply mass loss of Nieuwenhuijzen & de Jager (1990, A&A, 231, 134)         *
* for massive stars over the entire HRD with a metallicity factor            *
* from Kudritzki et al. (1989, A&A, 219, 205).                               *
* mdflag =< 2 - Level A                                                      *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(kw.eq.1)then
            x = 1.d0
         else
            x = MIN(1.d0,(lum-4000.d0)/500.d0)
         endif
         dms = 9.6d-15*x*(r**0.81d0)*(lum**1.24d0)*(mt**0.16d0)
         dms = dms*(z/0.02d0)**ind5
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate standard 'Reimers' mass loss giants [Kudritzki & Reimers (1978)] *
C Level B  if neta > 0, then follows [Reimers (1975)] = 0.477                *
C Level D  if neta < 0, then follows [Schroeder & Cuntz (2005)] = -0.172     *
* mdflag =< 2 - Level A                                                      *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*
      if(kw.ge.2)then
*
        if(neta.gt.0.d0)then
          dml = neta*4.0d-13*r*lum/mt
        else
          teffR = 1000.d0*((1130.d0*lum/
     &                    (r**2.d0))**(1.d0/4.d0))
          geffR = mt/(r*r)
          dml = abs(neta)*4.0d-13*(r*lum/mt)
          dml=dml*( teffR/4000.d0 )**3.5d0     ! Dependence with temperature
          dml=dml*( 1 + gsun/(4300.d0*geffR) ) ! Dependence with surface gravity
        endif
      endif
*
      if(mdflag.le.2)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Check for any tidally enhanced mass loss in binary systems (optional):     *
* see [Tout & Eggleton (1988)]                                               *
* mdflag =< 2 - Level A                                                      *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(rl.gt.0.d0)then
            dml = dml*(1.d0 + bwind*(MIN(0.5d0,(r/rl)))**6)
         endif
         dms = MAX(dms,dml)
      endif
*
      if(kw.le.6.and.mdflag.le.2)then
*
         if(kw.eq.5.or.kw.eq.6)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Enhanced Mass Loss during AGB phase                                        *
* Apply mass loss of Vassiliadis & Wood (1993, ApJ, 413, 641)                *
* for high pulsation periods on AGB.                                         *
* mdflag =< 2 - Level A                                                      *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            p0 = -2.07d0 - 0.9d0*log10(mt) + 1.94d0*log10(r)
            p0 = 10.d0**p0
            p0 = MIN(p0,2000.d0)
            dmt = -11.4d0+0.0125d0*(p0-100.d0*MAX(mt-2.5d0,0.d0))
            dmt = 10.d0**dmt
            dmt = 1.d0*MIN(dmt,1.36d-09*lum)
            dms = MAX(dms,dmt)
         endif
*
         mew = ((mt-mc)/mt)*MIN(5.d0,MAX(1.2d0,(lum/lum0)**kap))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Apply the reduced WR-like mass loss for small H-envelope mass              *
* as described in the [Hurley et al. (2000)] SSE paper.                      *
* mdflag =< 2 - Level A                                                      *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(mew.lt.1.d0)then
            dml = 1.0d-13*lum**(3.d0/2.d0)*(1.d0 - mew)
            dms = MAX(dms,dml)
         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* LBV-like mass loss beyond the Humphreys-Davidson limit                     *
* (see Humphreys & Davidson 1994 and Belczynski et al. 2010).                *
* mdflag =< 2 - Level A                                                      *     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         x = 1.0d-5*r*SQRT(lum)
         if(mdflag.eq.2.and.lum.gt.6.0d+05.and.x.gt.1.d0)then
            dml = 0.1d0*(x-1.d0)**3*(lum/6.0d+05-1.d0)
            dms = dms + dml
         endif
*
      elseif(kw.le.6.and.mdflag.gt.2)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
* LBV-like mass loss beyond the Humphreys-Davidson limit                     *
* (see Humphreys & Davidson 1994 and Belczynski et al. 2010).                *
* mdflag > 2  - Level B, C                                                   *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
         x = 1.0d-5*r*SQRT(lum)
         if(lum.gt.6.0d+05.and.x.gt.1.d0)then
            dms = 1.0d-04*flbv
         else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Apply mass loss for hot, massive H-rich O/B stars following                *
* Vink et al. (2001, A&A, 369,574). mdflag = 3 - Level B, C                  *
* mdflag = 4  - Level C (taking into account bi-stability jump at around     *
*                         25000K.                                            *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            if(teff.ge.12500.0.and.teff.le.25000.0.and.mdflag.eq.3)then 
               teff1 = MIN(teff,22500.d0)
               vw = 1.3d0/2.d0
               dml = -6.688d0 + 2.21d0*log10(lum/1.0d+05) -  
     &             1.339d0*log10(mt/30.d0) - 1.601d0*log10(vw) + 
     &             0.85d0*log10(z/0.02d0) + 1.07d0*log10(teff1/20000.d0)
               dms = 10.d0**dml
            elseif(teff.gt.12500.0.and.teff.le.50000.1)then
               teff1 = MAX(teff,27500.d0)
               vw = 2.6d0/2.d0
               t40 = log10(teff1/40000.d0)
               dml = -6.697d0 + 2.194d0*log10(lum/1.0d+05) -  
     &             1.313d0*log10(mt/30.d0) - 1.226d0*log10(vw) + 
     &             0.85d0*log10(z/0.02d0) + 
     &             0.933d0*t40*(1.d0 - 11.704d0*t40)
               dms = 10.d0**dml
            else
               if(lum.gt.4000.d0)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Apply mass loss of Nieuwenhuijzen & de Jager (1990, A&A, 231, 134)         *
* for massive stars over the entire HRD with a metallicity factor            *
* from Kudritzki et al. (1989, A&A, 219, 205). mdflag >= 3 - Level B, C      *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                x = MIN(1.d0,(lum-4000.d0)/500.d0)
                dms = 9.6d-15*x*(r**0.81d0)*(lum**1.24d0)*(mt**0.16d0)
                dms = dms*(z/0.02d0)**(1.d0/2.d0)
               endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Check for any tidally enhanced mass loss in binary systems (optional):     *
* see Tout & Eggleton (1988, MNRAS, 231, 823). mdflag >= 3 - Level B, C      *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               if(rl.gt.0.d0)then
                 dml = dml*(1.d0 + bwind*(MIN(0.5d0,(r/rl)))**6)
               endif
               dms = MAX(dms,dml)
*
               if(kw.eq.5.or.kw.eq.6)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC			   
* Apply mass loss of Vassiliadis & Wood (1993, ApJ, 413, 641)                *
* for high pulsation periods on AGB. mdflag >= 3 - Level B, C                *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                  p0 = -2.07d0 - 0.9d0*log10(mt) + 1.94d0*log10(r)
                  p0 = 10.d0**p0
                  p0 = MIN(p0,2000.d0)
                  dmt = -11.4d0+0.0125d0*(p0-100.d0*MAX(mt-2.5d0,0.d0))
                  dmt = 10.d0**dmt
                  dmt = 1.d0*MIN(dmt,1.36d-09*lum)
                  dms = MAX(dms,dmt)
               endif
*
               mew = ((mt-mc)/mt)*MIN(5.d0,MAX(1.2d0,(lum/lum0)**kap))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC			   
* Apply the reduced WR-like mass loss for small H-envelope mass              *
* as described in the Hurley, Pols & Tout (2000) SSE paper.                  *
* mdflag >= 3  - Level B, C                                                  *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               if(mew.lt.1.d0)then
                  dml = 1.0d-13*lum**(3.d0/2.d0)*(1.d0 - mew)
                  dms = MAX(dms,dml)
               endif
            endif
         endif
      endif
*
      if(kw.gt.6)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Add metallicity factor from Vink & de Koter (2005, A&A, 442, 587).         *
* mdflag >= 3  - Level B, C                                                  *
* MISSING IN NBODY6++GPU  - 23.07.2020                                       *
* PROPERLY IMPLEMENTED NOW!                                                  *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if((z/0.02d0).ge.10.d0)then 
           fz = 10.d0**ind4
         elseif((z/0.02d0).ge.0.001d0.and.(z/0.02d0).lt.10.d0)then
           fz = (z/0.02d0)**ind4                            
         else 
           fz = 0.001d0**ind4
         endif
*
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	  
* Apply mass loss of Hamann & Koesterke (1998, A&A, 335, 1003)               *
* for WR (naked helium) stars.                                               *
* mdflag >= 3  - Level B, C                                                  *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	
         dml = 1.0d-13*lum**(3.d0/2.d0)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	
* Add metallicity factor from Vink & de Koter (2005, A&A, 442, 587).         *
* mdflag = 4  - Level C (taking into account bi-stability jump).             *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	
         if(mdflag.ge.3)then 
           dml = dml*fz   ! This is fz in MOCCA BSE: fz=(z/0.02d0)**0.86d0
           dms = MAX(dms,dml)
        endif 
      endif
*
      mlwind = dms
*
      return
      end
***
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* SOURCES:                                                                   *               
* Belczynski K., Bulik T., Fryer C. L. et al. (2010, ApJ, 714, 1217-)        * 
* Chen Y., Bressan A., Girardi L. et al. (2015, MNRAS, 452, 1068-)           *               
* Giacobbo et al. (2018, MNRAS, 474, 2959–)                                  * 
* Graefener G., Vink J. S., de Koter A. et al. (2011, A&A, 535, 14-)         * 
* Hamann W. R. & Koesterke L. (1998, A&A, 335, 1003-)                        *              
* Humphreys R. M. & Davidson K. (1979, ApJ, 232, 409–)                       *
* Humphreys R. M. & Davidson K. (1994, PASP, 106, 1025-)                     *           
* Hurley J. R., Pols, O. R., & Tout, C. A. (2000, MNRAS, 315, 543-)  	     *
* Hurley J. R., Tout C. A., Pols O. R. et al. (2002, MNRAS, 329, 897-)       *                    
* Kudritzki R. P. & Reimers D. (1978, A&A, 70, 227-)                         *
* Kudritzki et al. (1989, A&A, 219, 205-)                                    *         
* McDonald I. & Zijlstra A. A. (2015, MNRAS, 448, 502-)                      *
* Nieuwenhuijzen H. & de Jager C. (1990, A&A, 231, 134-)                     *                    
* Vassiliadis A. & Wood P. R. (1993, ApJ, 413, 641-)                         *                    
* Vink J. S., de Koter A., Lamers H. J. G. L. M. (2001, A&A, 369, 574-)      *                    
* Vink J. S. & de Koter A. (2005, A&A, 442, 587-)                            *               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



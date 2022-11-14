***
SUBROUTINE evolv1(kw,mass,mt,r,lum,mc,rc,menv,renv,ospin,
   &                  epoch,tm,tphys,tphysf,dtp,z,zpars,vkick,
   &                vs)
c-------------------------------------------------------------c
c
c     Evolves a single star.
c     Mass loss is an option.
c     The timestep is not constant but determined by certain criteria.
c     Plots the HRD and variables as a function of time.
c
c     Written by Jarrod Hurley 26/08/97 at the Institute of
c     Astronomy, Cambridge.
c
c-------------------------------------------------------------c
c
c     STELLAR TYPES - KW
c
c        0 - deeply or fully convective low mass MS star
c        1 - Main Sequence star
c        2 - Hertzsprung Gap
c        3 - First Giant Branch
c        4 - Core Helium Burning
c        5 - First Asymptotic Giant Branch
c        6 - Second Asymptotic Giant Branch
c        7 - Main Sequence Naked Helium star
c        8 - Hertzsprung Gap Naked Helium star
c        9 - Giant Branch Naked Helium star
c       10 - Helium White Dwarf
c       11 - Carbon/Oxygen White Dwarf
c       12 - Oxygen/Neon White Dwarf
c       13 - Neutron Star
c       14 - Black Hole
c       15 - Massless Supernova
c
c-------------------------------------------------------------c
    implicit none
*
    integer kw,it,ip,jp,j,kwold,rflag
    integer nv
    parameter(nv=50000)
*
    real*8 mass,z,aj
    real*8 epoch,tphys,tphys2,tmold,tbgold
    real*8 mt,tm,tn,tphysf,dtp,tsave
    real*8 tscls(20),lums(10),GB(10),zpars(20)
    real*8 r,lum,mc,teff,rc,menv,renv,vs(3), vkick
    real*8 ospin,jspin,djt,djmb,k2,k3
    parameter(k3=0.21d0)
    real*8 m0,r1,lum1,mc1,rc1,menv1,renv1,k21
    real*8 dt,dtm,dtr,dr,dtdr,dms,dml,mt2,rl
    real*8 tol,tiny
    parameter(tol=1.0d-10,tiny=1.0d-14)
    real*8 ajhold,rm0,eps,alpha2
    parameter(eps=1.0d-06,alpha2=0.09d0)
    real*8 mlwind,vrotf
    external mlwind,vrotf
    logical iplot,isave
    REAL*8 neta,bwind,hewind,mxns
    COMMON /VALUE1/ neta,bwind,hewind,mxns
    REAL*8 pts1,pts2,pts3
    COMMON /POINTS/ pts1,pts2,pts3
    REAL scm(50000,14),spp(20,3)
    COMMON /SINGLE/ scm,spp
    INTEGER BHSPIN
    COMMON /VALUE6/ BHSPIN
    INTEGER KMECH
    REAL*8 FBFAC,FBTOT,MCO
    INTEGER ECS
    COMMON /FBACK/ FBFAC,FBTOT,MCO,ECS
    COMMON /FLAGS3/ kmech
    integer ceflag,tflag,ifflag,nsflag,wdflag,ecflag,
   &  psflag  
    REAL*8 aspin, aconst, bconst, alow, mone, mtwo, mbh ! Albrecht - 14.08.2020 - BH natal spins
    REAL*8 aone, bone, atwo, btwo ! Albrecht - 14.08.2020 - BH natal spins
    REAL*8 G,M_sun,rsunkm,parsec,Km,Kmps,cspeed,yearsc
*
*******************************************************************************
    PARAMETER (G=6.6743D-08, M_sun=1.9884D+33)
    PARAMETER (rsunkm=6.955D+10, parsec=3.0856776D+18)
    PARAMETER (Km=1.0D+05, Kmps=1.0D+05)
    PARAMETER (cspeed=3.0D+10, yearsc=3.154D+07)
*******************************************************************************
*
**************************************************************************
*
    dtm = 0.d0
    r = 0.d0
    lum = 0.d0
    mc = 0.d0
    mc1 = 0.d0
    rc = 0.d0
    rl = 0.d0
    if(ospin.le.0.d0)then
       ospin = 1.0d-10
       jspin = 1.0d-10
    endif
    k2 = 0.15d0
    rflag = 0
    vs(1) = 0.d0
    vs(2) = 0.d0
    vs(3) = 0.d0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Setup variables which control the output (if it is required).                *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    ip = 0
    jp = 0
    tsave = tphys
    isave = .true.
    iplot = .false.
    if(dtp.le.0.d0)then
       iplot = .true.
       isave = .false.
       tsave = tphysf
    elseif(dtp.gt.tphysf)then
       isave = .false.
       tsave = tphysf
    endif
* 
    do 10 , j = 1,nv
*
       if(neta.gt.tiny.and.j.gt.1)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate mass loss from the previous timestep.                              *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          dt = 1.0d+06*dtm
          dms = mlwind(kw,lum,r,mt,mc,rl,z)*dt
          if(kw.lt.10)then
             dml = mt - mc
             if(dml.lt.dms)then
                dtm = (dml/dms)*dtm
                dms = dml
             endif
          endif
       else
          dms = 0.d0
       endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Limit to 1% mass loss.                                                       *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       if(dms.gt.0.01d0*mt)then
          dtm = 0.01d0*mt*dtm/dms
          dms = 0.01d0*mt
       endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate the rate of angular momentum loss due to magnetic braking          *
* and/or mass loss.                                                            *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       if(j.gt.1)then
          djt = (2.d0/3.d0)*(dms/(1.0d+06*dtm))*r*r*ospin
          if(mt.gt.0.35d0.and.kw.lt.10)then
             djmb = 5.83d-16*menv*(r*ospin)**3/mt
             djt = djt + djmb
          endif
       endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Update mass and time and reset epoch for a MS (and possibly a HG) star.      *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       if(dms.gt.0.d0)then
          mt = mt - dms
          if(kw.le.2.or.kw.eq.7)then
             m0 = mass
             mc1 = mc
             mass = mt
             tmold = tm
             tbgold = tscls(1)
             CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
             if(kw.eq.2)then
                if(GB(9).lt.mc1.or.m0.gt.zpars(3))then
                   mass = m0
                else
                   epoch = tm + (tscls(1) - tm)*(ajhold-tmold)/
   &                            (tbgold - tmold)
                   epoch = tphys - epoch
                endif
             else
                epoch = tphys - ajhold*tm/tmold
             endif
          endif
       endif
       tphys2 = tphys
       tphys = tphys + dtm
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Find the landmark luminosities and timescales as well as setting             *
* the GB parameters.                                                           *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*         write(*,*) aj,tphys,epoch
       aj = tphys - epoch
       CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Find the current radius, luminosity, core mass and stellar type              *
* given the initial mass, current mass, metallicity and age                    *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       kwold = kw
       CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
   &               r,lum,kw,mc,rc,menv,renv,k2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* If mass loss has occurred and no type change then check that we              * 
* have indeed limited the radius change to 10%.                                *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       if(kw.eq.kwold.and.dms.gt.0.d0.and.rflag.ne.0)then
          mt2 = mt + dms
          dml = dms/dtm
          it = 0
20         dr = r - rm0
          if(ABS(dr).gt.0.1d0*rm0)then
             it = it + 1
             if(it.eq.20.and.kw.eq.4) goto 30
             if(it.gt.30)then
                WRITE(99,*)' DANGER1! ',it,kw,mass,dr,rm0
                WRITE(*,*)' STOP: EVOLV1 FATAL ERROR 1'
*                  CALL exit(0)
*                  STOP 
             endif
             dtdr = dtm/ABS(dr)
             dtm = alpha2*MAX(rm0,r)*dtdr
             if(it.ge.20) dtm = 0.5d0*dtm
             if(dtm.lt.1.0d-07*aj) goto 30
             dms = dtm*dml
             mt = mt2 - dms
             if(kw.le.2.or.kw.eq.7)then
                mass = mt
                CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
                if(kw.eq.2)then
                   if(GB(9).lt.mc1.or.m0.gt.zpars(3))then
                      mass = m0
                   else
                      epoch = tm + (tscls(1) - tm)*(ajhold-tmold)/
   &                               (tbgold - tmold)
                      epoch = tphys2 - epoch
                   endif
                else
                   epoch = tphys2 - ajhold*tm/tmold
                endif
             endif
             tphys = tphys2 + dtm
             aj = tphys - epoch
             mc = mc1
             CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
             CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
   &                     r,lum,kw,mc,rc,menv,renv,k2)
             goto 20
          endif
30         continue
       endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Initialize or adjust the spin of the star.                                   *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       if(j.eq.1)then
          if(tphys.lt.tiny.and.ospin.lt.0.001d0)then
*               write (*,*) vrotf(mt,0),mt,r
             ospin = 45.35d0*vrotf(mt,0)/r
          endif
          jspin = ospin*(k2*r*r*(mt-mc)+k3*rc*rc*mc)
       else
          jspin = MAX(1.0d-10,jspin - djt*1.0d+06*dtm)
          ospin = jspin/(k2*r*r*(mt-mc)+k3*rc*rc*mc)
       endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Test for changes in evolution type.                                          *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       if(j.eq.1.or.kw.ne.kwold)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Force new NS or BH to have a one second period.                              *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          if(kw.eq.14)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
** ADDITIONS in accordance with Nbody6++GPU - Albrecht - 14.08.2020 
******* BH Kerr Metric spin parameter *****
             if(BHSPIN.EQ.1) then
             WRITE(*,*) 'Using bhspin = 1 - Geneva BH spins'	   				
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	   
* BH natal spin from Geneva models (experimental) [Belczynski et al. (2020)]   *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                 if (z.lt.0.001D0) then
                 alow = 0.25D0
                 mtwo = 38.8D0
                 mone = 32.0D0
                 aconst = -0.088D0
                 bconst = 3.666D0
                 elseif (z.ge.0.001D0.and.z.lt.0.004D0) then
                 alow = 0.0D0
                 mtwo = 27.7D0
                 mone = 18.0D0
                 aconst = -0.088D0
                 bconst = 2.434D0
                 elseif (z.ge.0.004D0.and.z.lt.0.01D0) then
                 alow = 0.25D0
                 mtwo = 37.8D0
                 mone = 31.0D0
                 aconst = -0.088D0
                 bconst = 3.578D0
                 else
                 alow = 0.13D0
                 mtwo = 24.2D0
                 mone = 16.0D0
                 aconst = -0.088D0
                 bconst = 2.258D0
                 endif
                 if (MCO.le.mone) then
                    aspin = 0.85D0
                 elseif (MCO.gt.mone.and.MCO.lt.mtwo) then
                    aspin = (aconst*MCO) + bconst
                 else
                    aspin = alow
                 endif
               if (aspin.lt.0.0D0) aspin = 0.0D0
*********
              elseif(BHSPIN.EQ.2) then
              WRITE(*,*) 'Using bhspin = 2 - MESA BH spins'	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	   
* BH natal spin from MESA models (experimental) [Belczynski et al. (2020)]     *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                 if (z.lt.0.001D0) then
                 aone = -0.0010D0
                 bone = 0.125D0
                 atwo = 0.0D0
                 btwo = 0.0D0
                 mone = 1.0E+10
              elseif (z.ge.0.001D0.and.z.lt.0.004D0) then
                 aone = 0.0076D0
                 bone = 0.050D0
                 atwo = -0.0019D0
                 btwo = 0.165D0
                 mone = 12.09D0
              elseif (z.ge.0.004D0.and.z.lt.0.01D0) then
                 aone = -0.0006D0
                 bone = 0.105D0
                 atwo = 0.0D0
                 btwo = 0.0D0
                 mone = 1.0E+10
              else
                 aone = -0.0016D0
                 bone = 0.115D0
                 atwo = 0.0D0
                 btwo = 0.0D0
                 mone = 1.0D+10
              endif
              if (MCO.le.mone) then
                 aspin = (aone*MCO) + bone
              else
                 aspin = (atwo*MCO) + btwo
              endif
              if (aspin.lt.0.0D0) then
                   aspin = 0.0D0
              endif 
*********
            elseif(BHSPIN.EQ.0) then
               WRITE(*,*) 'Using bhspin = 0 - Fuller BH spins'	
*		   
* ------------ Zero BH natal spins --------------------------------------
*
              aspin = 0.0D0
            endif
*		   
 
*    
*		   
* ------------ gm -------------------------------------------------------
*
             mbh = mt
             WRITE(*,*) 'MBH', mbh
             mbh = mt*M_sun   ! Doublecheck mt!
*   WRITE(*,*) 'MBH in g', mbh
*		   
* ------------ gm cm^2 / s ----------------------------------------------
*		   
             jspin = (G*aspin*mbh**2)/cspeed
*    WRITE(*,*) 'jspin: (gm cm^2 / s) = ', jspin			   
*		   
* ------------ Msun Rsun^2 / year ---------------------------------------
*
             jspin = jspin*(yearsc/(M_sun*(rsunkm**2)))
*    WRITE(*,*) 'jspin: (Msun Rsun^2 / year) =', jspin		
        ospin = jspin/(k3*rc*rc*mc)
*    WRITE(*,*) 'Spin after kick (year^-1)', ospin
*		   
* ------------ Scaled ---------------------------------------------------
*
        WRITE(*,*) 'aspin', aspin
             CALL kick(kw,mass,mt,0.d0,0.d0,-1.d0,0.d0,vs)  ! Natal Spins BHs Kamlah 19.02.2021
          vkick = dsqrt(vs(1)*vs(1)+vs(2)*vs(2)+vs(3)*vs(3))
             WRITE(*,*) 'Spin (1/year)', ospin  ! This does not work properly
          WRITE(*,*) 'Using BHSPIN =', bhspin
            endif
       if(kw.eq.13)then 
         ospin = 2.0d+08
              jspin = k3*rc*rc*mc*ospin
          CALL kick(kw,mass,mt,0.d0,0.d0,-1.d0,0.d0,vs)  ! NS seperate
              vkick = dsqrt(vs(1)*vs(1)+vs(2)*vs(2)+vs(3)*vs(3))
            endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* BELOW WD KICK CALL - ALbrecht Kamlah - 10.02.2021                            *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
             if(kw.eq.10.or.kw.eq.11.or.kw.eq.12)then
*              ospin = 2.0d+08  
             jspin = k3*rc*rc*mc*ospin
             CALL kick(kw,mass,mt,0.d0,0.d0,-1.d0,0.d0,vs)
             vkick = dsqrt(vs(1)*vs(1)+vs(2)*vs(2)+vs(3)*vs(3))
          endif
*
          jp = jp + 1
          spp(jp,1) = tphys
          spp(jp,2) = float(kw)
          if(kw.eq.15)then
             spp(jp,3) = mass 
             goto 90
          else
             spp(jp,3) = mt
          endif
*
       endif

       if(kw.eq.15)then
         goto 90
       endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* record values for plotting and reset epoch.                                  *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       epoch = tphys - aj
       if((isave.and.tphys.ge.tsave).or.iplot)then
*
          ip = ip + 1
          scm(ip,1) = tphys
          scm(ip,2) = float(kw)
          scm(ip,3) = mass
          scm(ip,4) = mt
          scm(ip,5) = log10(lum)
          scm(ip,6) = log10(r)
          teff = 1000.d0*((1130.d0*lum/(r**2.d0))**(1.d0/4.d0))
          scm(ip,7) = log10(teff)
          scm(ip,8) = mc
          scm(ip,9) = rc
          scm(ip,10) = menv
          scm(ip,11) = renv
          scm(ip,12) = epoch
          scm(ip,13) = ospin
          if(isave) tsave = tsave + dtp
          if(tphysf.lt.tiny)then
             ip = ip + 1
             do 35 , it = 1,13
                scm(ip,it) = scm(ip-1,it)
35            continue
          endif
*
       endif
*
       if(tphys.ge.tphysf)then
*
          jp = jp + 1
          spp(jp,1) = tphys
          spp(jp,2) = float(kw)
          spp(jp,3) = mt
*
          goto 90
       endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Record radius and current age.                                               * 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       rm0 = r
       ajhold = aj
       if(kw.ne.kwold) kwold = kw
       CALL deltat(kw,aj,tm,tn,tscls,dtm,dtr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Check for type change.                                                       *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       it = 0
       m0 = mass
       if((dtr-dtm).le.tol.and.kw.le.9)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Check final radius for too large a jump.                                     *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          aj = MAX(aj,aj*(1.d0-eps)+dtr)
          mc1 = mc 
          CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
   &                  r1,lum1,kw,mc1,rc1,menv1,renv1,k21)
          dr = r1 - rm0
          if(ABS(dr).gt.0.1d0*rm0)then
             dtm = dtr - ajhold*eps
             dtdr = dtm/ABS(dr)
             dtm = alpha2*MAX(r1,rm0)*dtdr
             goto 40
          else
             dtm = dtr
             goto 50
          endif
       endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Limit to a 10% increase in radius assuming no further mass loss              *
* and thus that the pertubation functions due to small envelope mass           *
* will not change the radius.                                                  *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
40      aj = ajhold + dtm
       mc1 = mc 
       CALL hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
   &               r1,lum1,kw,mc1,rc1,menv1,renv1,k21)
       dr = r1 - rm0
       it = it + 1
       if(it.eq.20.and.kw.eq.4) goto 50
       if(it.gt.1000)then
          WRITE(99,*)' DANGER2! ',it,kw,mass,dr,rm0
          goto 50
*            WRITE(*,*)' STOP: EVOLV1 FATAL ERROR '
*            CALL exit(0)
*            STOP 
       endif
       if(ABS(dr).gt.0.1d0*rm0)then
          dtdr = dtm/ABS(dr)
          dtm = alpha2*MAX(rm0,r1)*dtdr
          if(it.ge.20) dtm = 0.5d0*dtm
          goto 40
       endif
*
50      continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Ensure that change of type has not occurred during radius check.             *
* This is rare but may occur for HG stars of ZAMS mass > 50 Msun.              *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       if(kw.ne.kwold)then
          kw = kwold
          mass = m0
          CALL star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
       endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Choose minimum of time-scale and remaining interval (> 100 yrs).             *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       dtm = MAX(dtm,1.0d-07*aj)
       dtm = MIN(dtm,tsave-tphys)
*
10   continue
*
90   continue
*
    tphysf = tphys
    scm(ip+1,1) = -1.0
    spp(jp+1,1) = -1.0
    if(ip.ge.nv)then
       WRITE(99,*)' EVOLV1 ARRAY ERROR ',mass
       WRITE(*,*)' STOP: EVOLV1 ARRAY ERROR 3'
*         CALL exit(0)
*         STOP
    endif
*
    RETURN
    END
***

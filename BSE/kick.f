***
      subroutine kick(kw,m1,m1n,m2,ecc,sep,jorb,vs)
      implicit none
      INTEGER kw,k
      INTEGER idum
      COMMON /VALUE3/ idum
      REAL*8 m1,m2,m1n,ecc,sep,jorb,ecc2,z
      REAL*8 pi,twopi,gmrkm
      REAL*8 mm,em,dif,der,del,r
      REAL*8 u1,u2,vk,v(4),s,theta,phi
      REAL*8 sphi,cphi,stheta,ctheta,salpha,calpha
      REAL*8 vr,vr2,vk2,vn2,hn2
      REAL*8 mu,cmu,vs(3),v1,v2,mx1,mx2
      REAL ran3,xx
      EXTERNAL ran3
      REAL*8 DISP,DISP0
      INTEGER BHFLAG 
      INTEGER KMECH
      REAL*8 FBFAC,FBTOT,MCO
      INTEGER ECS
      COMMON /FBACK/ FBFAC,FBTOT,MCO,ECS
      COMMON /FLAGS3/ KMECH
      COMMON /VALUE4/ DISP,BHFLAG
      integer ceflag,tflag,ifflag,nsflag,wdflag,ecflag,
     &  psflag
      COMMON /FLAGS/ ceflag,tflag,ifflag,nsflag,wdflag
      COMMON /FLAGS1/ psflag,ecflag
      REAL*8 aspin, aconst, bconst, alow, mone, mtwo, JSPIN, MBH ! Albrecht - 14.08.2020 - BH natal spins
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
      real*8 CONVF,MNUEFF,MNSTYP
      real*8 ECSIG,WDSIG1,WDSIG2,WDKMAX
      CONVF = 5.0D0
      MNUEFF = 7.0D0
      MNSTYP = 1.4D0
      ECSIG = 3.0D0
      WDSIG1 = 2.0D0
      WDSIG2 = 2.0D0
      WDKMAX = 6.0D0

*
*******************************************************************************
*      ------------------------------------------------------   
*      | Velocity kick for WD, neutron stars or black holes |
*      ------------------------------------------------------
*
* There are various choices that can be made for kicks. 
* Make sure you are comfortable with the choices set (see below) 
* as this will critically affect retention statistics. 
*
* For regular NSs the kick choice is controlled by the value of 
* DISP (sigma in BSE). Choices are: 
*
*    DISP < 0.0 - kick drawn randomly between 0 - ABS(DISP) km/s
*    DISP = 0.0 - no kick
*    DISP > 0.0 - kick drawn from a Maxwellian with dispersion DISP. 
*
* You may also choose to have the kick velocity distribution scaled 
* by VSTAR (i.e. scaled by the initial escape velocity). 
* To do this set VFAC to a non-zero value and VFAC*VSTAR will be 
* either the maximum of the flat distribution (DISP < 0) or 
* the dispersion of the Maxwellian (DISP > 0). 
* DO NOT USE USE VSTAR VFAC FOR STANDALONE BSE!!!
*
* Then for an electron capture supernova or an accretion-induced 
* collapse the choice is determined by the value of ECSIG set 
* internally here. Choices are: 
*
*    ECSIG = 0.0 - no kick
*    ECSIG > 0.0 - kick drawn from a Maxwellian, dispersion ECSIG 
*    ECSIG < 0.0 - same as for the regular NSs but scaled by ABS(ECSIG).
*
* These ECSNe are identified by their mass of 1.26 Msun. 
*
* For BHs the kick choice is controlled by the value of BHFLAG.  
* Choices are: 
*
*    BHFLAG = 0 - no kick
*    BHFLAG = 1 - same as for the regular NSs
*    BHFLAG = 2 - same as for the regular NSs but scaled by fallback. 
*
* Small kicks for WDs can also be set if KZ(25) > 0 in the input file. 
* In this case you can distinguish: 
*
*    WDSIG1 - He and COWDs
*    WDSIG2 - ONeWDs 
*
* as the dispersion in a Maxwellian for the different WD types. 
* A limit of WDKMAX is set. 
* See Fellhauer et al. (2003, ApJ, 595, L53) for more on WD kicks.
*
*    KMECH = 1 - standard momentum-conserving,
*    KMECH = 2 - convection-asymmetry-driven,
*    KMECH = 3 - collapse-asymmetry-driven,
*    KMECH = 4 - neutrino-driven
*
* It is assumed that one of these four mechanisms is the primary driver
* of SN natal kicks that we observe
*
* CONVF: convective boost factor larger CO core masses 
*        in the case of convection-asymmerty-driven
*        kick mechanism (typical range: 2.0-10.0)
*
* MNUEFF: in case of neutrino-driven kick mechanism, the 
*         effective remnant mass beyond which the neutrino emission does not
*         enhance significantly as the remnant (baryonic) mass is increased
*         (typical range: 5.0-10.0 Msun)
*
* MNSTYP: typical mass of a neutron star with the input dispersion velocity 'DISP'
*
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* bhflag - BH kicks and natal spins, same as for CC NSs but reduced for momentum cons. if fallback
* Level A           = 0: no BH kick.
* Level B			= 1: same as for the regular NSs but scaled by fallback. 
* Level D			= 2: Fallback as bhflag = 1 + BH Kerr Metric spin parameter 
*                                                 MESA models [Belczynski et al. (2020)]
* Level D			= 3: Fallback as bhflag = 1 + BH Kerr Metric spin parameter
*                                                 Geneva models [Belczynski et al. (2020)]
* kmech - Kick mechanism for NSs and BHs [Fryer suggestions at around 2018]
* Level A           = 1: standard momentum-conserving 
* Level B			= 2: convection-asymmetry-driven
* Level B			= 3: collapse-asymmetry-driven
* Level C			= 4: neutrino-driven
* disp - Dispersion in a Maxwellian velocity distribution for kick distributions (WDs,NSs,BHs)
* Level A           = 190 km/s following [Hobbs et al. (2005)]
* Level B           = 265 km/s following [Belczynski et al. (2008)]			
* WDSIG1 - kicks for He and COWDs
* Level B           = 2.0 km/s [Fellhauer et al. (2003)]
* WDSIG2 - kicks for ONeWDs 
* Level B           = 2.0 km/s [Fellhauer et al. (2003)]
* WDKMAX - maximum kick velocity for all WDs
* Level B           = 6.0 km/s [Fellhauer et al. (2003)]
* ECSIG - Electron capture supernova (ECSNe) or an accretion-induced supernove (AISNe)
* Level A           = 190 km/s following [Hobbs et al. (2005)]
* Level B           = 3.0 km/s [Gessner & Janka (2018)]
* CONVF - Convective boost factor larger CO core masses (Convection-asymmetry driven kick) 2 =< CONVF =< 10)
* Level B           = 5.0 
* MNUEFF -  - effective remnant mass beyond which the neutrino emission does not 
*             enhance significantly as the remnant (baryonic) mass is increased 
*             (Neutrino-driven driven kick) 5 M* =< MNUEFF =< 10 M* 
* Level B           = 7.0 M*
* MNSTYP - typical mass of a neutron star with the input dispersion velocity 'disp'
* Level B           = 1.4 M*
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*
*
      do k = 1,3
         vs(k) = 0.d0
      enddo
*     if(kw.eq.14.and.bhflag.eq.0) goto 95
*
      pi = ACOS(-1.d0)
      twopi = 2.d0*pi
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Conversion factor to ensure velocities are in km/s using mass and
* radius in solar units.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      gmrkm = 1.906125d+05
*
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Find the initial separation by randomly choosing a mean anomaly.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(sep.gt.0.d0.and.ecc.ge.0.d0)then
         xx = RAN3(idum)
         mm = xx*twopi
         em = mm
 2       dif = em - ecc*SIN(em) - mm
         if(ABS(dif/mm).le.1.0d-04) goto 3
         der = 1.d0 - ecc*COS(em)
         del = dif/der
         em = em - del
         goto 2
 3       continue
         r = sep*(1.d0 - ecc*COS(em))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Find the initial relative velocity vector.                                   *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         salpha = SQRT((sep*sep*(1.d0-ecc*ecc))/(r*(2.d0*sep-r)))
         calpha = (-1.d0*ecc*SIN(em))/SQRT(1.d0-ecc*ecc*COS(em)*COS(em))
         vr2 = gmrkm*(m1+m2)*(2.d0/r - 1.d0/sep)
         vr = SQRT(vr2)
      else
         vr = 0.d0
         vr2 = 0.d0
         salpha = 0.d0
         calpha = 0.d0
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* ADDITIONS in accordance with Nbody6++GPU - Albrecht - 31.07.2020
      DISP0 = MAX(DISP,0.D0)
      IF(KW.EQ.10.OR.KW.EQ.11) DISP0 = MAX(WDSIG1,0.D0)
      IF(KW.EQ.12) DISP0 = MAX(WDSIG2,0.D0)
      IF(ECS.EQ.1)THEN
        IF(ECSIG.LT.-0.01)THEN
           DISP0 = DISP0*ABS(ECSIG)
        ELSE
           DISP0 = MAX(ECSIG,0.D0)
        ENDIF
      ENDIF
      IF(KW.EQ.14.AND.BHFLAG.EQ.0) DISP0 = 0.D0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Generate Kick Velocity using Maxwellian Distribution (Phinney 1992).         * 
* Use Henon's method for pairwise components (Douglas Heggie 22/5/97).         *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do 20 k = 1,2
         u1 = RAN3(idum)
         u2 = RAN3(idum)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Generate two velocities from polar coordinates S & THETA.                    *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         s = DISP0*SQRT(-2.d0*LOG(1.d0 - u1))
         theta = twopi*u2
         v(2*k-1) = s*COS(theta)
         v(2*k) = s*SIN(theta)
 20   continue
*      
	 IF(DISP0.GT.0.001D0)THEN ! Albrecht additions
         vk2 = v(1)**2 + v(2)**2 + v(3)**2
         vk = SQRT(vk2)		 
	  ELSE 
         vk2 = 0.d0
         vk = 0.d0
         if(kw.lt.0) kw = 13
      endif
	  v(4) = vk
	  sphi = -1.d0 + 2.d0*u1
      phi = ASIN(sphi)
      cphi = COS(phi)
      stheta = SIN(theta)
      ctheta = COS(theta)
*     WRITE(66,*)' KICK VK PHI THETA ',vk,phi,theta
      if(sep.le.0.d0.or.ecc.lt.0.d0) goto 90
*
* Impose the maximum WD kick velocity. 


      IF(KW.GE.10.AND.KW.LE.12.AND.VK.GT.WDKMAX)THEN
         VK = WDKMAX
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
** ADDITIONS in accordance with Nbody6++GPU - Albrecht - 31.07.2020 
* Restrict the BH / NS kick velocity by fallback. 
* This could be done better but in the N-body code we only have 
* limited information.
* See arXiv: 1902.07718 (Banerjee, S., Belczynski, C., Fryer, C., et al.)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       IF (BHFLAG.GT.1) THEN
         IF (KW.EQ.13.OR.KW.EQ.14) THEN
*
* ------------ Skip ECS-NS -------------------------------------------
*
            IF(ECS.EQ.0)THEN
*
* ------------ Standard momentum-conserving kick ---------------------
*
               VK = v(4)*(1.0D0 - FBFAC)
*
* ------------ Convection-asymmetry-driven kick ----------------------
*     
               IF ((KMECH.EQ.2).AND.(MCO.LE.3.5D0))
     &              VK = SQRT(vk2)*(MNSTYP/m1n)
               IF ((KMECH.EQ.2).AND.(MCO.GT.3.5D0))
     &              VK = SQRT(vk2)*(MNSTYP/m1n)*CONVF
*
* ------------ Collapse-asymmetry-driven kick ------------------------
*
               IF ((KMECH.EQ.3).AND.(MCO.LE.3.0D0))
     &              VK = SQRT(vk2)*(MNSTYP/m1n)
               IF ((KMECH.EQ.3).AND.(MCO.GT.3.0D0))
     &              VK = SQRT(vk2)*(MNSTYP/m1n)*0.1D0
*
* ------------ Neutrino-driven kick ----------------------------------
*
               IF (KMECH.EQ.4)
     &              VK = v(4)*(MIN(m1n,MNUEFF)/m1n)
	        ENDIF
               WRITE(*,*) "Using fallback scaled kicks"
	    ENDIF
	 ENDIF
*
** ADDITIONS in accordance with Nbody6++GPU - Albrecht - 14.08.2020 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(sep.le.0.d0.or.ecc.lt.0.d0) goto 90
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Determine the magnitude of the new relative velocity.                        *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      vn2 = vk2+vr2-2.d0*vk*vr*(ctheta*cphi*salpha-stheta*cphi*calpha)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate the new semi-major axis.                                           *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      sep = 2.d0/r - vn2/(gmrkm*(m1n+m2))
      sep = 1.d0/sep
*     if(sep.le.0.d0)then
*        ecc = 1.1d0
*        goto 90
*     endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Determine the magnitude of the cross product of the separation vector        *
* and the new relative velocity.                                               *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      v1 = vk2*sphi*sphi
      v2 = (vk*ctheta*cphi-vr*salpha)**2
      hn2 = r*r*(v1 + v2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate the new eccentricity.                                              *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ecc2 = 1.d0 - hn2/(gmrkm*sep*(m1n+m2))
      ecc2 = MAX(ecc2,0.d0)
      ecc = SQRT(ecc2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate the new orbital angular momentum taking care to convert            *
* hn to units of Rsun^2/yr.                                                    *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      jorb = (m1n*m2/(m1n+m2))*SQRT(hn2)*(yearsc/rsunkm)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Determine the angle between the new and old orbital angular                  *
* momentum vectors.                                                            *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      cmu = (vr*salpha-vk*ctheta*cphi)/SQRT(v1 + v2)
      mu = ACOS(cmu)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate the components of the velocity of the new centre-of-mass.          *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 90   continue
      if(ecc.le.1.0)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate the components of the velocity of the new centre-of-mass.          *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         mx1 = vk*m1n/(m1n+m2)
         mx2 = vr*(m1-m1n)*m2/((m1n+m2)*(m1+m2))
         vs(1) = mx1*ctheta*cphi + mx2*salpha
         vs(2) = mx1*stheta*cphi + mx2*calpha
         vs(3) = mx1*sphi
      else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate the relative hyperbolic velocity at infinity (simple method).      *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         sep = r/(ecc-1.d0)
*        cmu = SQRT(ecc-1.d0)
*        mu = ATAN(cmu)
         mu = ACOS(1.d0/ecc)
         vr2 = gmrkm*(m1n+m2)/sep
         vr = SQRT(vr2)
         vs(1) = vr*SIN(mu)
         vs(2) = vr*COS(mu)
         vs(3) = 0.d0
         ecc = MIN(ecc,99.99d0)
      endif
*
 95   continue
*
      RETURN
      END
***

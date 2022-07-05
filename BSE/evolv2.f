***
      SUBROUTINE evolv2(kstar,mass0,mass,rad,lumin,massc,radc,
     &                  menv,renv,ospin,epoch,tms,
     &                  tphys,tphysf,dtp,z,zpars,tb,ecc,vkick,
     &                  vs1,vs2)
      implicit none
*      INCLUDE 'const_bse.h'
***
*
*           B I N A R Y
*           ***********
*
*       Roche lobe overflow.
*       --------------------
*
*       Developed by Jarrod Hurley, IOA, Cambridge.
*       .........................................................
*
*       Advice by Christopher Tout, Onno Pols & Sverre Aarseth.
*       ++++++++++++++++++++++++++++++++++++++++++++++++++
*
* Adapted from Aarseth's code 21st September 1996.
* Fully revised on 27th November 1996 to remove vestiges of N-body code and
* incorporate corrections.
* Fully revised on 1st April 1998 to include new stellar evolution formulae
* and associated binary evolution changes.
* Fully revised on 4th July 1998 to include eccentricity, tidal 
* circularization, wind accretion, velocity kicks for supernovae and all
* associated orbital momentum changes.
* Revised on 31st October 2000 to upgrade restrictions imposed on the 
* timestep owing to magnetic braking and orbital angular momentum changes. 
*
***
*
* See Tout et al., 1997, MNRAS, 291, 732 for a description of many of the
* processes in this code as well as the relevant references mentioned
* within the code. 
*
* Reference for the stellar evolution formulae is Hurley, Pols & Tout, 
* 2000, MNRAS, 315, 543 (SSE paper).  
* Reference for the binary evolution algorithm is Hurley, Tout & Pols, 
* 2002, MNRAS, 329, 897 (BSE paper). 
*
***
*
* March 2001 *
* Changes since version 3, i.e. since production of Paper3:  
*
* 1) The Eddington limit flag (on/off) has been replaced by an 
*    Eddington limit multiplicative factor (eddfac). So if you 
*    want to neglect the Eddington limit you would set eddfac 
*    to a large value. 
*
* 2) To determine whether material transferred during RLOF forms 
*    an accretion disk around the secondary or hits the secondary 
*    in a direct stream we calculate a minimum radial distance, rmin, 
*    of the mass stream from the secondary. This is taken from eq.(1) 
*    of Ulrich & Burger (1976, ApJ, 206, 509) which they fitted to 
*    the calculations of Lubow & Shu (1974, ApJ, 198, 383).
*    If rmin is less than the radius of the secondary then an 
*    accretion disk is not formed. 
*    Note that the formula for rmin given by Ulrich & Burger is valid 
*    for all q whereas that given by Nelemans et al. (2001, A&A, 
*    submitted) in their eq.(6) is only valid for q < 1 where 
*    they define q = Mdonor/Maccretor, i.e. DD systems. 
*
* 3) The changes to orbital and spin angular momentum owing to 
*    RLOF mass transfer have been improved, and an new input option 
*    now exists. 
*    When mass is lost from the system during RLOF there are now 
*    three choices as to how the orbital angular momentum is 
*    affected: a) the lost material carries with it a fraction 
*    gamma of the orbital angular momentum, i.e. 
*    dJorb = gamma*dm*a^2*omega_orb; b) the material carries with it 
*    the specific angular momentum of the primary, i.e.
*    dJorb = dm*a_1^2*omega_orb; or c) assume the material is lost 
*    from the system as if a wind from the secondary, i.e.
*    dJorb = dm*a_2^2*omega_orb. 
*    The parameter gamma is an input option. 
*    Choice c) is used if the mass transfer is super-Eddington 
*    or the system is experiencing novae eruptions. 
*    In all other cases choice a) is used if gamma > 0.0, b) if 
*    gamma = -1.0 and c) is used if gamma = -2.0. 
*    The primary spin angular momentum is reduced by an amount 
*    dm1*r_1^2*omega_1 when an amount of mass dm1 is transferred 
*    from the primary. 
*    If the secondary accretes through a disk then its spin 
*    angular momentum is altered by assuming that the material 
*    falls onto the star from the inner edge of a Keplerian 
*    disk and that the system is in a steady state, i.e. 
*    an amount dm2*SQRT(G*m_2*r_2). 
*    If there is no accretion disk then we calculate the angular 
*    momentum of the transferred material by using the radius at 
*    at which the disk would have formed (rdisk = 1.7*rmin, see 
*    Ulrich & Burger 1976) if allowed, i.e. the angular momentum 
*    of the inner Lagrangian point, and add this directly to 
*    the secondary, i.e. an amount dm2*SQRT(G*m_2*rdisk). 
*    Total angular momentum is conserved in this model. 
*
* 4) Now using q_crit = 3.0 for MS-MS Roche systems (previously we 
*    had nothing). This corresponds roughly to R proportional to M^5 
*    which should be true for the majority of the MS (varies from 
*    (M^17 -> M^2). If q > q_crit then contact occurs. 
*    For CHeB primaries we also take q_crit = 3.0 and allow 
*    common-envelope to occur if this is exceeded. 
*
* 5) The value of lambda used in calculations of the envelope binding 
*    energy for giants in common-envelope is now variable (see function 
*    in zfuncs). The lambda function has been fitted by Onno to detailed 
*    models ... he will write about this soon!
*
* 6) Note that eq.42 in the paper is missing a SQRT around the
*    MR^2/a^5 part. This needs to be corrected in any code update
*    paper with a thanks to Jeremy Sepinsky (student at NorthWestern).
*    It is ok in the code.
*
* March 2003 *
* New input options added:  
*
*    ifflag - for the mass of a WD you can choose to use the mass that 
*             results from the evolution algorithm (basically a competition
*             between core-mass growth and envelope mass-loss) or use the IFMR 
*             proposed by Han, Podsiadlowski & Eggleton, 1995, MNRAS, 272, 800
*             [>0 activates HPE IFMR]. 
*
*    wdflag - for the cooling of WDs you can choose to use either the standard 
*             Mestel cooling law (see SSE paper) or a modified-Mestel law that 
*             is better matched to detailed models (provided by Brad Hansen 
*             ... see Hurley & Shara, 2003, ApJ, May 20, in press) 
*             [>0 activates modified-Mestel]. 
*
*    bhflag - choose whether or not black holes should get velocity kicks 
*             at formation 
*             [0= no kick; >0 kick]. 
*
*    nsflag - for the mass of neutron stars and black holes you can use either 
*             the SSE prescription or the prescription presented by 
*             Belczynski et al. 2002, ApJ, 572, 407 who found that SSE was 
*             underestimating the masses of these stars. In either case you also 
*             need to set the maximum NS mass (mxns) for the prescription  
*             [0= SSE, mxns=1.8; >0 Belczynski, mxns=3.0].
*
* Sept 2004 *
* Input options added/changed:
*
*    ceflag - set to 3 this uses de Kool (or Podsiadlowski) CE prescription, 
*             other options, such as Yungelson, could be added as well. 
*
*    hewind - factor to control the amount of He star mass-loss, i.e.
*             1.0e-13*hewind*L^(2/3) gives He star mass-loss.
*
* NOTE: some versions may have contained a bug in the calculation of the 
*       f factor for convective tides. The incorrect line was: 
*             f = MIN(1.d0,(ttid/(2.d0*tc)**2))
*
*
*       ++++++++++++++++++++++++++++++++++++++++++++++++++
***
*
      INTEGER loop,iter,intpol,k,ip,jp,j1,j2
      PARAMETER(loop=20000000)
      INTEGER kstar(2),kw,kst,kw1,kw2,kmin,kmax
      INTEGER ktype(0:14,0:14)
      COMMON /TYPES/ ktype
      INTEGER ceflag,tflag,ifflag,nsflag,wdflag
      COMMON /FLAGS/ ceflag,tflag,ifflag,nsflag,wdflag
*
      REAL*8 km,km0,tphys,tphys0,dtm0,tphys00
      REAL*8 tphysf,dtp,tsave
      REAL*8 aj(2),aj0(2),epoch(2),tms(2),tbgb(2),tkh(2),dtmi(2)
      REAL*8 mass0(2),mass(2),massc(2),menv(2),mass00(2),mcxx(2)
      REAL*8 rad(2),rol(2),rol0(2),rdot(2),radc(2),renv(2),radx(2)
      REAL*8 lumin(2),k2str(2),q(2),dms(2),dmr(2),dmt(2),vkick(2)
      REAL*8 vs1(3),vs2(3)
      REAL*8 dml,vorb2,vwind2,omv2,ivsqm,lacc,vs(3)
      REAL*8 sep,dr,tb,dme,tdyn,taum,dm1,dm2,dmchk,qc,dt,pd,rlperi
      REAL*8 m1ce,m2ce,mch,tmsnew,dm22,mew
      PARAMETER(mch=1.44d0)
      REAL*8 yeardy,aursun,yearsc
      PARAMETER(yeardy=365.24d0,aursun=214.95d0,yearsc=3.1557d+07)
      REAL*8 acc1,tiny
      PARAMETER(acc1=3.920659d+08,tiny=1.0d-14)
      REAL*8 ecc,ecc1,tc,tcirc,ttid,ecc2,omecc2,sqome2,sqome3,sqome5
      REAL*8 f1,f2,f3,f4,f5,f,raa2,raa6,eqspin,rg2,tcqr
      REAL*8 k3,mr23yr,twopi
      PARAMETER(k3=0.21d0,mr23yr=0.4311d0)
      REAL*8 jspin(2),ospin(2),jorb,oorb,jspbru,ospbru
      REAL*8 delet,delet1,dspint(2),djspint(2),djtx(2)
      REAL*8 dtj,djorb,djgr,djmb,djt,djtt,rmin,rdisk
      REAL*8 neta,bwind,hewind,mxns
      COMMON /VALUE1/ neta,bwind,hewind,mxns
      REAL*8 beta,xi,acc2,epsnov,eddfac,gamma
      COMMON /VALUE5/ beta,xi,acc2,epsnov,eddfac,gamma
*
      REAL*8 z,tm,tn,m0,mt,rm,lum,mc,rc,me,re,k2,age,dtm,dtr
      REAL*8 tscls(20),lums(10),GB(10),zpars(20)
      REAL*8 zero,ngtv,ngtv2,mt2,rrl1,rrl2,mcx,teff1,teff2
      REAL*8 mass1i,mass2i,tbi,ecci
      LOGICAL coel,com,prec,inttry,change,snova,sgl,bsymb,esymb,bss
      LOGICAL supedd,novae,disk
      LOGICAL isave,iplot
      REAL*8 rl,mlwind,vrotf,corerd
      EXTERNAL rl,mlwind,vrotf,corerd
      REAL bcm(5000000,45),bpp(8000000,10)
      COMMON /BINARY/ bcm,bpp
      INTEGER BHSPIN
      COMMON /VALUE6/ BHSPIN
      REAL*8 aspin, aconst, b_const, alow, mone, mtwo, mbh ! Albrecht - 14.08.2020 - BH natal spins
      REAL*8 aone, bone, atwo, btwo ! Albrecht - 14.08.2020 - BH natal spins
      REAL*8 G_grav,M_sun,rsunkm,parsec,K_m,Kmps,cspeed
      REAL*8 FBFAC,FBTOT,MCO
      INTEGER ECS
      COMMON /FBACK/ FBFAC,FBTOT,MCO,ECS
*
*******************************************************************************
      PARAMETER (G_grav=6.6743D-08, M_sun=1.9884D+33)
      PARAMETER (rsunkm=6.955D+10, parsec=3.0856776D+18)
       PARAMETER (K_m=1.0D+05)
*      PARAMETER (Kmps=1.0D+05)
      PARAMETER (cspeed=3.0D+10)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Save the initial state.                                                      *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      mass1i = mass0(1)
      mass2i = mass0(2)
      tbi = tb
      ecci = ecc
*
      vs1(1) = 0.d0
      vs1(2) = 0.d0
      vs1(3) = 0.d0
      vs2(1) = 0.d0
      vs2(2) = 0.d0
      vs2(3) = 0.d0
      vkick(1)=0.d0
      vkick(2)=0.d0
      zero = 0.d0
      ngtv = -1.d0
      ngtv2 = -2.d0
      twopi = 2.d0*ACOS(-1.d0)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Initialize the parameters.                                                   *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      kmin = 1
      kmax = 2
      sgl = .false.
      mt2 = MIN(mass(1),mass(2))
      kst = 0
*
      if(mt2.lt.tiny.or.tb.le.0.d0)then
         sgl = .true.
         if(mt2.lt.tiny)then
            mt2 = 0.d0
            if(mass(1).lt.tiny)then
               if(tphys.lt.tiny)then
                  mass0(1) = 0.01d0
                  mass(1) = mass0(1)
                  kst = 1
               else
                  kmin = 2
                  lumin(1) = 1.0d-10
                  rad(1) = 1.0d-10
                  massc(1) = 0.d0
                  dmt(1) = 0.d0
                  dmr(1) = 0.d0
               endif
               ospin(1) = 1.0d-10
               jspin(1) = 1.0d-10
            else
               if(tphys.lt.tiny)then
                  mass0(2) = 0.01d0
                  mass(2) = mass0(2)
                  kst = 2
               else
                  kmax = 1
                  lumin(2) = 1.0d-10
                  rad(2) = 1.0d-10
                  massc(2) = 0.d0
                  dmt(2) = 0.d0
                  dmr(2) = 0.d0
               endif
               ospin(2) = 1.0d-10
               jspin(2) = 1.0d-10
            endif
         endif
         ecc = -1.d0
         tb = 0.d0
         sep = 1.0d+10
         oorb = 0.d0
         jorb = 0.d0
         if(ospin(1).lt.0.0) ospin(1) = 1.0d-10
         if(ospin(2).lt.0.0) ospin(2) = 1.0d-10
         q(1) = 1.0d+10
         q(2) = 1.0d+10
         rol(1) = 1.0d+10
         rol(2) = 1.0d+10
      else
         tb = tb/yeardy
         sep = aursun*(tb*tb*(mass(1) + mass(2)))**(1.d0/3.d0)
         oorb = twopi/tb
         jorb = mass(1)*mass(2)/(mass(1)+mass(2))
     &          *SQRT(1.d0-ecc*ecc)*sep*sep*oorb
         if(ospin(1).lt.0.d0) ospin(1) = oorb
         if(ospin(2).lt.0.d0) ospin(2) = oorb
      endif
*
      do 500 , k = kmin,kmax
         age = tphys - epoch(k)
         mc = massc(k)
         rc = radc(k)
         CALL star(kstar(k),mass0(k),mass(k),tm,tn,tscls,lums,GB,zpars)
         CALL hrdiag(mass0(k),age,mass(k),tm,tn,tscls,lums,GB,zpars,
     &               rm,lum,kstar(k),mc,rc,me,re,k2)
         aj(k) = age
         epoch(k) = tphys - age
         rad(k) = rm
         lumin(k) = lum
*         teff(k) = 1000.d0*((1130.d0*lumin(k)/
*     &                    (rad(k)**2.d0))**(1.d0/4.d0))
         massc(k) = mc
         radc(k) = rc
         menv(k) = me
         renv(k) = re
         k2str(k) = k2
         tms(k) = tm
         tbgb(k) = tscls(1)
*
         if(tphys.lt.tiny.and.ospin(k).le.0.001d0)then
            ospin(k) = 45.35d0*vrotf(mass(k))/rm
         endif
         jspin(k) = ospin(k)*(k2*rm*rm*(mass(k)-mc)+k3*rc*rc*mc)
         if(.not.sgl)then
            q(k) = mass(k)/mass(3-k)
            rol(k) = rl(q(k))*sep
         endif
         rol0(k) = rol(k)
         dmr(k) = 0.d0
         dmt(k) = 0.d0
         djspint(k) = 0.d0
         dtmi(k) = 1.0d+06
*
 500  continue
*
      if(mt2.lt.tiny)then
         sep = 0.d0
         if(kst.gt.0)then
            mass0(kst) = 0.d0
            mass(kst) = 0.d0
            kmin = 3 - kst
            kmax = kmin
         endif
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* On the first entry the previous timestep is zero to prevent mass loss.       *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      dtm = 0.d0   !changed from 0.0d0   Kamlah - 17.02.2021
      delet = 0.d0
      djorb = 0.d0
      bss = .false.
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
      if(tphys.ge.tphysf) goto 140
*
 4    iter = 0
      intpol = 0
      inttry = .false.
      change = .false.
      prec = .false.
      snova = .false.
      coel = .false.
      com = .false.
      bsymb = .false.
      esymb = .false.
      tphys0 = tphys
      ecc1 = ecc
      j1 = 1
      j2 = 2
      if(kstar(1).ge.10.and.kstar(1).le.14) dtmi(1) = 0.01d0
      if(kstar(2).ge.10.and.kstar(2).le.14) dtmi(2) = 0.01d0
      dm1 = 0.d0
      dm2 = 0.d0
*
 5    kw1 = kstar(1)
      kw2 = kstar(2)
*
      dt = 1.0d+06*dtm
      eqspin = 0.d0
      djtt = 0.d0
*
      if(intpol.eq.0.and.ABS(dtm).gt.tiny.and..not.sgl)then
         vorb2 = acc1*(mass(1)+mass(2))/sep
         ivsqm = 1.d0/SQRT(1.d0-ecc*ecc)
         do 501 , k = 1,2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate wind mass loss from the previous timestep.                         *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            if(neta.gt.tiny)then
               rlperi = rol(k)*(1.d0-ecc)
               dmr(k) = mlwind(kstar(k),lumin(k),rad(k),mass(k),
     &                         massc(k),rlperi,z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate how much of wind mass loss from companion will be                  *
* accreted (Boffin & Jorissen, A&A 1988, 205, 155).                            *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               if(beta.lt.0.d0)then !PK. following startrack
                  beta = 0.125
                  if(kstar(k).le.1)then
                     if(mass(k).gt.120.d0)then
                        beta = 7.d0
                     elseif(mass(k).le.1.4d0)then
                        beta = 0.5
                     else
                        beta = 7.d0*((mass(k)-1.4d0)/(120.d0-1.4d0))
     &                         + 0.5d0
                     endif
                  elseif(kstar(k).ge.7.and.kstar(k).le.9)then
                     if(mass(k).gt.120.d0)then
                        beta = 7.d0
                     elseif(mass(k).le.10.d0)then
                        beta = 0.125
                     else
                        beta = 7.d0*((mass(k)-10.d0)/(120.d0-10.d0))
     &                               + 0.125d0
                     endif
                  endif
               endif
               vwind2 = 2.d0*beta*acc1*mass(k)/rad(k)
               omv2 = (1.d0 + vorb2/vwind2)**(3.d0/2.d0)
               dmt(3-k) = ivsqm*acc2*dmr(k)*((acc1*mass(3-k)/vwind2)**2)
     &                    /(2.d0*sep*sep*omv2)
               dmt(3-k) = MIN(dmt(3-k),0.8d0*dmr(k))
            else
               dmr(k) = 0.d0
               dmt(3-k) = 0.d0
            endif
 501     continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Diagnostic for Symbiotic-type stars.                                         *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(neta.gt.tiny.and..not.esymb)then
            lacc = 3.14d+07*mass(j2)*dmt(j2)/rad(j2)
            lacc = lacc/lumin(j1)
            if((lacc.gt.0.01d0.and..not.bsymb).or.
     &         (lacc.lt.0.01d0.and.bsymb))then
               jp = MIN(80,jp + 1)
               bpp(jp,1) = tphys
               bpp(jp,2) = mass(1)
               bpp(jp,3) = mass(2)
               bpp(jp,4) = float(kstar(1))
               bpp(jp,5) = float(kstar(2))
               bpp(jp,6) = sep
               bpp(jp,7) = ecc
               bpp(jp,8) = rad(1)/rol(1)
               bpp(jp,9) = rad(2)/rol(2)
               if(bsymb)then
                  bpp(jp,10) = 13.0
                  esymb = .true.
               else
                  bpp(jp,10) = 12.0
                  bsymb = .true.
               endif
            endif
         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate orbital angular momentum change due to wind mass loss.             *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         ecc2 = ecc*ecc
         omecc2 = 1.d0 - ecc2
         sqome2 = SQRT(omecc2)
*
         djorb = ((dmr(1)+q(1)*dmt(1))*mass(2)*mass(2) + 
     &            (dmr(2)+q(2)*dmt(2))*mass(1)*mass(1))*
     &           sep*sep*sqome2*oorb/(mass(1)+mass(2))**2
         delet = ecc*(dmt(1)*(0.5d0/mass(1) + 1.d0/(mass(1)+mass(2))) +
     &                dmt(2)*(0.5d0/mass(2) + 1.d0/(mass(1)+mass(2))))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* For very close systems include angular momentum loss owing to                *
* gravitational radiation.                                                     *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(sep.le.10.d0)then
            djgr = 8.315d-10*mass(1)*mass(2)*(mass(1)+mass(2))/
     &             (sep*sep*sep*sep)
            f1 = (19.d0/6.d0) + (121.d0/96.d0)*ecc2
            sqome5 = sqome2**5
            delet1 = djgr*ecc*f1/sqome5
            djgr = djgr*jorb*(1.d0+0.875d0*ecc2)/sqome5
            djorb = djorb + djgr
            delet = delet + delet1
         endif
*
         do 502 , k = 1,2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate change in the intrinsic spin of the star.                          *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            djtx(k) = (2.d0/3.d0)*xi*dmt(k)*rad(3-k)*rad(3-k)*ospin(3-k)
            djspint(k) = (2.d0/3.d0)*(dmr(k)*rad(k)*rad(k)*ospin(k)) -
     &                   djtx(k)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Include magnetic braking for stars that have appreciable convective          *
* envelopes. This includes MS stars with M < 1.25, HG stars near the GB        *
* and giants. MB is not allowed for fully convective MS stars.                 *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            if(mass(k).gt.0.35d0.and.kstar(k).lt.10)then
               djmb = 5.83d-16*menv(k)*(rad(k)*ospin(k))**3/mass(k)
               djspint(k) = djspint(k) + djmb
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Limit to a 3% angular momentum change for the star owing to MB.              *
* This is found to work best with the maximum iteration of 20000,              *
* i.e. does not create an excessive number of iterations, while not            *
* affecting the evolution outcome when compared with a 2% restriction.         *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               if(djmb.gt.tiny)then
                  dtj = 0.03d0*jspin(k)/ABS(djmb)
                  dt = MIN(dt,dtj)
               endif
            endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate circularization, orbital shrinkage and spin up.                    *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            dspint(k) = 0.d0
            if(((kstar(k).le.9.and.rad(k).ge.0.01d0*rol(k)).or.
     &         (kstar(k).ge.10.and.k.eq.j1)).and.tflag.gt.0)then
*
               raa2 = (rad(k)/sep)**2
               raa6 = raa2**3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Hut's polynomials.                                                           *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               f5 = 1.d0+ecc2*(3.d0+ecc2*0.375d0)
               f4 = 1.d0+ecc2*(1.5d0+ecc2*0.125d0)
               f3 = 1.d0+ecc2*(3.75d0+ecc2*(1.875d0+ecc2*7.8125d-02))
               f2 = 1.d0+ecc2*(7.5d0+ecc2*(5.625d0+ecc2*0.3125d0))
               f1 = 1.d0+ecc2*(15.5d0+ecc2*(31.875d0+ecc2*(11.5625d0
     &                  +ecc2*0.390625d0)))
*
               if((kstar(k).eq.1.and.mass(k).ge.1.25d0).or.
     &            kstar(k).eq.4.or.kstar(k).eq.7)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Radiative damping (Zahn, 1977, A&A, 57, 383 and 1975, A&A, 41, 329).         *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                  tc = 1.592d-09*(mass(k)**2.84d0)
                  f = 1.9782d+04*SQRT((mass(k)*rad(k)*rad(k))/sep**5)*
     &                tc*(1.d0+q(3-k))**(5.d0/6.d0)
                  tcqr = f*q(3-k)*raa6
                  rg2 = k2str(k)
               elseif(kstar(k).le.9)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Convective damping (Hut, 1981, A&A, 99, 126).                                *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                  tc = mr23yr*(menv(k)*renv(k)*(rad(k)-0.5d0*renv(k))/
     &                 (3.d0*lumin(k)))**(1.d0/3.d0)
                  ttid = twopi/(1.0d-10 + ABS(oorb - ospin(k)))
*                  f = MIN(1.d0,(ttid/(2.d0*tc)**2))
                  f = MIN(1.d0,(ttid/(2.d0*tc))**2)
                  tcqr = 2.d0*f*q(3-k)*raa6*menv(k)/(21.d0*tc*mass(k))
                  rg2 = (k2str(k)*(mass(k)-massc(k)))/mass(k)
               else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Degenerate damping (Campbell, 1984, MNRAS, 207, 433)                         *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                  f = 7.33d-09*(lumin(k)/mass(k))**(5.d0/7.d0)
                  tcqr = f*q(3-k)*q(3-k)*raa2*raa2/(1.d0+q(3-k))
                  rg2 = k3
               endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Circularization.                                                             *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               sqome3 = sqome2**3
               delet1 = 27.d0*tcqr*(1.d0+q(3-k))*raa2*(ecc/sqome2**13)*
     &                  (f3 - (11.d0/18.d0)*sqome3*f4*ospin(k)/oorb)
               tcirc = ecc/(ABS(delet1) + 1.0d-20)
               delet = delet + delet1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Spin up of star.                                                             *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               dspint(k) = (3.d0*q(3-k)*tcqr/(rg2*omecc2**6))*
     &                     (f2*oorb - sqome3*f5*ospin(k))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate the equilibrium spin at which no angular momentum                  *
* can be transferred.                                                          *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               eqspin = oorb*f2/(sqome3*f5)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate angular momentum change for the star owing to tides.               *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               djt = (k2str(k)*(mass(k)-massc(k))*rad(k)*rad(k) +
     &                k3*massc(k)*radc(k)*radc(k))*dspint(k)
               if(kstar(k).le.6.or.ABS(djt)/jspin(k).gt.0.1d0)then
                  djtt = djtt + djt
               endif
            endif
 502     continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Limit to 2% orbital angular momentum change.                                 *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         djtt = djtt + djorb 
         if(ABS(djtt).gt.tiny)then
            dtj = 0.02d0*jorb/ABS(djtt)
            dt = MIN(dt,dtj)
         endif
         dtm = dt/1.0d+06
*
      elseif(ABS(dtm).gt.tiny.and.sgl)then
         do 503 , k = kmin,kmax
            if(neta.gt.tiny)then
               rlperi = 0.d0
               dmr(k) = mlwind(kstar(k),lumin(k),rad(k),mass(k),
     &                         massc(k),rlperi,z)
            else
               dmr(k) = 0.d0
            endif
            dmt(k) = 0.d0
            djspint(k) = (2.d0/3.d0)*dmr(k)*rad(k)*rad(k)*ospin(k)
            if(mass(k).gt.0.35d0.and.kstar(k).lt.10)then
               djmb = 5.83d-16*menv(k)*(rad(k)*ospin(k))**3/mass(k)
               djspint(k) = djspint(k) + djmb
               if(djmb.gt.tiny)then
                  dtj = 0.03d0*jspin(k)/ABS(djmb)
                  dt = MIN(dt,dtj)
               endif
            endif
 503     continue
         dtm = dt/1.0d+06
      endif
*
      do 504 , k = kmin,kmax
*
         dms(k) = (dmr(k) - dmt(k))*dt
         if(kstar(k).lt.10)then
            dml = mass(k) - massc(k)
            if(dml.lt.dms(k))then
               dml = MAX(dml,2.d0*tiny)
               dtm = (dml/dms(k))*dtm
               if(k.eq.2) dms(1) = dms(1)*dml/dms(2)
               dms(k) = dml
               dt = 1.0d+06*dtm
            endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Limit to 1% mass loss.                                                       *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            if(dms(k).gt.0.01d0*mass(k))then
               dtm = 0.01d0*mass(k)*dtm/dms(k)
               if(k.eq.2) dms(1) = dms(1)*0.01d0*mass(2)/dms(2)
               dms(k) = 0.01d0*mass(k)
               dt = 1.0d+06*dtm
            endif
         endif
*
 504  continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Update mass and intrinsic spin (checking that the star is not spun           *
* past the equilibrium) and reset epoch for a MS (and possibly a HG) star.     *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do 505 , k = kmin,kmax
*
         if(eqspin.gt.0.d0.and.ABS(dspint(k)).gt.tiny)then
            if(intpol.eq.0)then
               if(dspint(k).ge.0.d0)then
                  dspint(k) = MIN(dspint(k),(eqspin-ospin(k))/dt)
               else
                  dspint(k) = MAX(dspint(k),(eqspin-ospin(k))/dt)
               endif
               djt = (k2str(k)*(mass(k)-massc(k))*rad(k)*rad(k) +
     &                k3*massc(k)*radc(k)*radc(k))*dspint(k)
               djorb = djorb + djt
               djspint(k) = djspint(k) - djt
            endif
         endif
*
         jspin(k) = MAX(1.0d-10,jspin(k) - djspint(k)*dt)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Ensure that the star does not spin up beyond break-up.                       *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         ospbru = twopi*SQRT(mass(k)*aursun**3/rad(k)**3)
         jspbru = (k2str(k)*(mass(k)-massc(k))*rad(k)*rad(k) +
     &             k3*massc(k)*radc(k)*radc(k))*ospbru
         if(jspin(k).gt.jspbru.and.ABS(dtm).gt.tiny)then
            mew = 1.d0
            if(djtx(k).gt.0.d0)then
               mew = MIN(mew,(jspin(k) - jspbru)/djtx(k))
            endif
            jspin(k) = jspbru
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* If excess material should not be accreted, activate next line.               *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*           dms(k) = dms(k) + (1.d0 - mew)*dmt(k)*dt
         endif
*
         if(ABS(dms(k)).gt.tiny)then
            mass(k) = mass(k) - dms(k)
            if(kstar(k).le.2.or.kstar(k).eq.7)then
               m0 = mass0(k)
               mass0(k) = mass(k)
               CALL star(kstar(k),mass0(k),mass(k),tm,tn,tscls,
     &                   lums,GB,zpars)
               if(kstar(k).eq.2)then
                  if(GB(9).lt.massc(k).or.m0.gt.zpars(3))then
                     mass0(k) = m0
                  else
                     epoch(k) = tm + (tscls(1) - tm)*(aj(k)-tms(k))/
     &                               (tbgb(k) - tms(k))
                     epoch(k) = tphys - epoch(k)
                  endif
               else
                  epoch(k) = tphys - aj(k)*tm/tms(k)
               endif
            endif
         endif
*
 505  continue
*
      if(.not.sgl)then
*
         ecc1 = ecc1 - delet*dt
         ecc = MAX(ecc1,0.d0)
         if(ecc.lt.1.0d-10) ecc = 0.d0
*
         if(ecc.ge.1.d0) goto 135
*
         jorb = jorb - djorb*dt
         sep = (mass(1) + mass(2))*jorb*jorb/
     &         ((mass(1)*mass(2)*twopi)**2*aursun**3*(1.d0-ecc*ecc))
         tb = (sep/aursun)*SQRT(sep/(aursun*(mass(1)+mass(2))))
         oorb = twopi/tb
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Advance the time.                                                            *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(intpol.eq.0)then
         tphys0 = tphys
         dtm0 = dtm
      endif
      tphys = tphys + dtm
*
      do 6 , k = kmin,kmax
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Acquire stellar parameters (M, R, L, Mc & K*) at apparent evolution age.     *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         age = tphys - epoch(k)
         aj0(k) = age
         kw = kstar(k)
         m0 = mass0(k)
         mt = mass(k)
         mc = massc(k)
         if(intpol.eq.0) mcxx(k) = mc
         if(intpol.gt.0) mc = mcxx(k)
         mass00(k) = m0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Masses over 100Msun should probably not be trusted in the                    *
* evolution formulae. ! Commented out - Albrecht 12.08.2020                    *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*         if(mt.gt.100.d0)then
*            WRITE(99,*)' MASS EXCEEDED ',mass1i,mass2i,tbi,ecci,mt
*            goto 140
*         endif
*
         CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars)
         CALL hrdiag(m0,age,mt,tm,tn,tscls,lums,GB,zpars,
     &               rm,lum,kw,mc,rc,me,re,k2)
*
         if(kw.ne.15)then
            ospin(k) = jspin(k)/(k2*(mt-mc)*rm*rm+k3*mc*rc*rc)
         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* At this point there may have been a supernova.                               *
* MODIFIED - ALBRECHT KAMLAH - 10.20.2021 - include WD kicks                   *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*         if(kw.ne.kstar(k).and.kstar(k).le.12.and.
          if(kw.ne.kstar(k).and.kstar(k).le.9.and.
     &      (kw.gt.9.and.kw.lt.15))then
            if(sgl)then
               CALL kick(kw,mass(k),mt,0.d0,0.d0,-1.d0,0.d0,vs)
               vkick(k) = dsqrt(vs(1)*vs(1)+vs(2)*vs(2)+vs(3)*vs(3))
               if(k.eq.1) then
                   vs1(1)=vs(1)
                   vs1(2)=vs(2)
                   vs1(3)=vs(3)
               elseif(k.eq.2.) then
                   vs2(1)=vs(1)
                   vs2(2)=vs(2)
                   vs2(3)=vs(3)
               endif
            else
               CALL kick(kw,mass(k),mt,mass(3-k),ecc,sep,jorb,vs)
               vkick(k) = dsqrt(vs(1)*vs(1)+vs(2)*vs(2)+vs(3)*vs(3))
               if(k.eq.1) then
                   vs1(1)=vs(1)
                   vs1(2)=vs(2)
                   vs1(3)=vs(3)
               elseif(k.eq.2.) then
                   vs2(1)=vs(1)
                   vs2(2)=vs(2)
                   vs2(3)=vs(3)
               endif
               if(ecc.gt.1.d0)then
                  kstar(k) = kw
                  mass(k) = mt
                  epoch(k) = tphys - age
                  goto 135
               endif
               tb = (sep/aursun)*SQRT(sep/(aursun*(mt+mass(3-k))))
               oorb = twopi/tb
            endif
            snova = .true.
         endif
*
         if(kw.ne.kstar(k))then
            change = .true.
            mass(k) = mt
            dtmi(k) = 0.01d0
            if(kw.eq.15)then
               kstar(k) = kw
               goto 135
            endif
            mass0(k) = m0
            epoch(k) = tphys - age
            if(kw.gt.6.and.kstar(k).le.6)then
               bsymb = .false.
               esymb = .false.
            endif
         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Force new NS or BH to have a second period.                                  *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(kstar(k).eq.13)then
            if(tphys-epoch(k).lt.tiny)then
               ospin(k) = 2.0d+08
               jspin(k) = k3*rc*rc*mc*ospin(k)
            endif
          endif
          if(kstar(k).eq.14)then
            if(tphys-epoch(k).lt.tiny)then
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
                   b_const = 3.666D0
                   elseif (z.ge.0.001D0.and.z.lt.0.004D0) then
                   alow = 0.0D0
                   mtwo = 27.7D0
                   mone = 18.0D0
                   aconst = -0.088D0
                   b_const = 2.434D0
                   elseif (z.ge.0.004D0.and.z.lt.0.01D0) then
                   alow = 0.25D0
                   mtwo = 37.8D0
                   mone = 31.0D0
                  aconst = -0.088D0
                   b_const = 3.578D0
                   else
                   alow = 0.13D0
                   mtwo = 24.2D0
                   mone = 16.0D0
                   aconst = -0.088D0
                   b_const = 2.258D0
                   endif
                   if (MCO.le.mone) then
                      aspin = 0.85D0
                   elseif (MCO.gt.mone.and.MCO.lt.mtwo) then
                      aspin = (aconst*MCO) + b_const
                   else
                     aspin = alow
                  endif
 		          if (aspin.lt.0.0D0) then
                   aspin = 0.0D0
				  endif
			 endif
*********
              if(BHSPIN.EQ.2) then
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
			endif 
            if(BHSPIN.EQ.0) then
            WRITE(*,*) 'Using bhspin = 0 - Fuller BH spins'	
*		   
* ------------ Zero BH natal spins --------------------------------------
*
                aspin = 0.0D0
             endif    
*		   
* ------------ gm -------------------------------------------------------
*
               mbh = mt
			   WRITE(*,*) 'MBH', mbh
               mbh = mt*M_sun   ! Doublecheck mt!
			   WRITE(*,*) 'MBH in g', mbh
*		   
* ------------ gm cm^2 / s ----------------------------------------------
*		   
               jspin(k) = (G_grav*aspin*mbh**2)/cspeed
			   WRITE(*,*) 'jspin: (gm cm^2 / s) = ', jspin			   
*		   
* ------------ Msun Rsun^2 / year ---------------------------------------
*
               jspin(k) = jspin(k)*(yearsc/(M_sun*(rsunkm**2)))
			   WRITE(*,*) 'jspin: (Msun Rsun^2 / year) =', jspin	
		       ospin(k) = jspin(k)/(k3*rc*rc*mc)
			   WRITE(*,*) 'Spin after kick (year^-1)', ospin
*		   
* ------------ Scaled ---------------------------------------------------
*
       WRITE(*,*) 'aspin', aspin
*       WRITE(*,*) 'Spin after kick', ospin, jspin   ! This does not work properly
*      WRITE(*,*) 'Using BHSPIN =', bhspin	
            endif
          endif
*       endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Set radius derivative for later interpolation.                               *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(ABS(dtm).gt.tiny)then
            rdot(k) = ABS(rm - rad(k))/dtm
         else
            rdot(k) = 0.d0
         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Base new time scale for changes in radius & mass on stellar type.            *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         dt = dtmi(k)
         CALL deltat(kw,age,tm,tn,tscls,dt,dtr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Choose minimum of time-scale and remaining interval.                         *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         dtmi(k) = MIN(dt,dtr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Save relevent solar quantities.                                              *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         aj(k) = age
         kstar(k) = kw
         rad(k) = rm
         lumin(k) = lum
         massc(k) = mc
         radc(k) = rc
         menv(k) = me
         renv(k) = re
         k2str(k) = k2
         tms(k) = tm
         tbgb(k) = tscls(1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Check for blue straggler formation.                                          *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(kw.le.1.and.tm.lt.tphys.and..not.bss)then
            bss = .true.
            jp = MIN(8000000,jp + 1)
            bpp(jp,1) = tphys
            bpp(jp,2) = mass(1)
            bpp(jp,3) = mass(2)
            bpp(jp,4) = float(kstar(1))
            bpp(jp,5) = float(kstar(2))
            bpp(jp,6) = sep
            bpp(jp,7) = ecc
            bpp(jp,8) = rad(1)/rol(1)
            bpp(jp,9) = rad(2)/rol(2)
            bpp(jp,10) = 14.0
         endif
*
 6    continue
*
      if(.not.sgl)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Determine the mass ratios.                                                   *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         do 506 , k = 1,2
            q(k) = mass(k)/mass(3-k)
 506     continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Determine the Roche lobe radii and adjust the radius derivative.             *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         do 507 , k = 1,2
            rol(k) = rl(q(k))*sep
            if(ABS(dtm).gt.tiny)then
               rdot(k) = rdot(k) + (rol(k) - rol0(k))/dtm
               rol0(k) = rol(k)
            endif
 507     continue
      else
         do 508 , k = kmin,kmax
            rol(k) = 10000.d0*rad(k)
 508     continue
      endif
*
      if((tphys.lt.tiny.and.ABS(dtm).lt.tiny.and.
     &    (mass2i.lt.0.1d0.or..not.sgl)).or.snova)then
         jp = MIN(80,jp + 1)
         bpp(jp,1) = tphys
         bpp(jp,2) = mass(1)
         bpp(jp,3) = mass(2)
         bpp(jp,4) = float(kstar(1))
         bpp(jp,5) = float(kstar(2))
         bpp(jp,6) = sep
         bpp(jp,7) = ecc
         bpp(jp,8) = rad(1)/rol(1)
         bpp(jp,9) = rad(2)/rol(2)
         bpp(jp,10) = 1.0
         if(snova)then
            bpp(jp,10) = 2.0
            dtm = 0.d0
            goto 4
         endif
      endif
*
      if((isave.and.tphys.ge.tsave).or.iplot)then
         if(sgl.or.(rad(1).lt.rol(1).and.rad(2).lt.rol(2)).
     &      or.tphys.lt.tiny)then
            ip = ip + 1
            bcm(ip,1) = tphys
            bcm(ip,2) = float(kstar(1))
            bcm(ip,3) = mass0(1)
            bcm(ip,4) = mass(1)
            bcm(ip,5) = log10(lumin(1))
            bcm(ip,6) = log10(rad(1))
            teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                       (rad(1)**2.d0))**(1.d0/4.d0))
            bcm(ip,7) = log10(teff1)
            bcm(ip,8) = massc(1)
            bcm(ip,9) = radc(1)
            bcm(ip,10) = menv(1)
            bcm(ip,11) = renv(1)
            bcm(ip,12) = epoch(1)
            bcm(ip,13) = ospin(1)
            bcm(ip,14) = dmt(1) - dmr(1)
            bcm(ip,15) = rad(1)/rol(1)
            bcm(ip,16) = float(kstar(2))
            bcm(ip,17) = mass0(2)
            bcm(ip,18) = mass(2)
            bcm(ip,19) = log10(lumin(2))
            bcm(ip,20) = log10(rad(2))
            teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                       (rad(2)**2.d0))**(1.d0/4.d0))
            bcm(ip,21) = log10(teff2)
            bcm(ip,22) = massc(2)
            bcm(ip,23) = radc(2)
            bcm(ip,24) = menv(2)
            bcm(ip,25) = renv(2)
            bcm(ip,26) = epoch(2)
            bcm(ip,27) = ospin(2)
            bcm(ip,28) = dmt(2) - dmr(2)
            bcm(ip,29) = rad(2)/rol(2)
            bcm(ip,30) = tb
            bcm(ip,31) = sep
            bcm(ip,32) = ecc
            if(isave) tsave = tsave + dtp
         endif
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* If not interpolating set the next timestep.                                  *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(intpol.eq.0)then
         dtm = MAX(1.0d-07*tphys,MIN(dtmi(1),dtmi(2)))
         dtm = MIN(dtm,tsave-tphys)
         if(iter.eq.0) dtm0 = dtm
      endif
      if(sgl) goto 98
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Set j1 to the donor - the primary                                            *
* and j2 to the accretor - the secondary.                                      *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(intpol.eq.0)then
         if(rad(1)/rol(1).ge.rad(2)/rol(2))then
            j1 = 1
            j2 = 2
         else
            j1 = 2
            j2 = 1
         endif
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Test whether Roche lobe overflow has begun.                                  *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(rad(j1).gt.rol(j1))then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Interpolate back until the primary is just filling its Roche lobe.           *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(rad(j1).ge.1.002d0*rol(j1))then
            if(intpol.eq.0) tphys00 = tphys
            intpol = intpol + 1
            if(iter.eq.0) goto 7
            if(inttry) goto 7
            if(intpol.ge.100)then
               WRITE(99,*)' INTPOL EXCEEDED ',mass1i,mass2i,tbi,ecci
               goto 140 
            endif
            dr = rad(j1) - 1.001d0*rol(j1)
            if(ABS(rdot(j1)).lt.tiny.or.prec)then
               goto 7
            endif
            dtm = -dr/ABS(rdot(j1))
            if(ABS(tphys0-tphys).gt.tiny) dtm = MAX(dtm,tphys0-tphys)
            if(kstar(1).ne.kw1)then
               kstar(1) = kw1
               mass0(1) = mass00(1)
               epoch(1) = tphys - aj0(1)
            endif
            if(kstar(2).ne.kw2)then
               kstar(2) = kw2
               mass0(2) = mass00(2)
               epoch(2) = tphys - aj0(2)
            endif
            change = .false.
         else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Enter Roche lobe overflow                                                    *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            if(tphys.ge.tphysf) goto 140
            goto 7
         endif
      else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Check if already interpolating.                                              *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(intpol.gt.0)then
            intpol = intpol + 1
            if(intpol.ge.80)then
               inttry = .true.
            endif
            if(ABS(rdot(j1)).lt.tiny)then
               prec = .true.
               dtm = 1.0d-07*tphys
            else
               dr = rad(j1) - 1.001d0*rol(j1)
               dtm = -dr/ABS(rdot(j1))
            endif
            if((tphys+dtm).ge.tphys00)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* If this occurs then most likely the star is a high mass type 4               *
* where the radius can change very sharply or possibly there is a              *
* discontinuity in the radius as a function of time and HRDIAG                 *
* needs to be checked!                                                         *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               dtm = 0.5d0*(tphys00 - tphys0)
               dtm = MAX(dtm,1.0d-10)
               prec = .true.
            endif
            tphys0 = tphys
         endif
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Check for collision at periastron.                                           *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      pd = sep*(1.d0 - ecc)
      if(pd.lt.(rad(1)+rad(2)).and.intpol.eq.0) goto 130
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Go back for the next step or interpolation.                                  *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 98   continue
      if(tphys.ge.tphysf.and.intpol.eq.0) goto 140
      if(change)then
         change = .false.
         jp = MIN(8000000,jp + 1)
         bpp(jp,1) = tphys
         bpp(jp,2) = mass(1)
         bpp(jp,3) = mass(2)
         bpp(jp,4) = float(kstar(1))
         bpp(jp,5) = float(kstar(2))
         bpp(jp,6) = sep
         bpp(jp,7) = ecc
         bpp(jp,8) = rad(1)/rol(1)
         bpp(jp,9) = rad(2)/rol(2)
         bpp(jp,10) = 2.0
      endif
*
      iter = iter + 1
*
      if(iter.ge.loop)then
         WRITE(99,*)' MAXIMUM ITER EXCEEDED ',mass1i,mass2i,tbi,ecci
         goto 140
      endif
      goto 5
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Set the nuclear timescale in years and slow-down factor.                     *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 7    km0 = dtm0*1.0d+03/tb
      if(km0.lt.tiny) km0 = 0.5d0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Force co-rotation of primary and orbit to ensure that the tides do not       *
* lead to unstable Roche (not currently used).                                 *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*     if(ospin(j1).gt.1.05d0*oorb)then
*        ospin(j1) = oorb
*        jspin(j1) = (k2str(j1)*rad(j1)*rad(j1)*(mass(j1)-massc(j1))+
*    &                k3*radc(j1)*radc(j1)*massc(j1))*ospin(j1)
*     endif
*
      iter = 0
      coel = .false.
      change = .false.
      radx(j1) = MAX(radc(j1),rol(j1))
      radx(j2) = rad(j2)
*
      jp = MIN(8000000,jp + 1)
      bpp(jp,1) = tphys
      bpp(jp,2) = mass(1)
      bpp(jp,3) = mass(2)
      bpp(jp,4) = float(kstar(1))
      bpp(jp,5) = float(kstar(2))
      bpp(jp,6) = sep
      bpp(jp,7) = ecc
      bpp(jp,8) = rad(1)/rol(1)
      bpp(jp,9) = rad(2)/rol(2)
      bpp(jp,10) = 3.0
*
      if(iplot.and.tphys.gt.tiny)then
         ip = ip + 1
         bcm(ip,1) = tphys
         bcm(ip,2) = float(kstar(1))
         bcm(ip,3) = mass0(1)
         bcm(ip,4) = mass(1)
         bcm(ip,5) = log10(lumin(1))
         bcm(ip,6) = log10(rad(1))
         teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                    (rad(1)**2.d0))**(1.d0/4.d0))
         bcm(ip,7) = log10(teff1)
         bcm(ip,8) = massc(1)
         bcm(ip,9) = radc(1)
         bcm(ip,10) = menv(1)
         bcm(ip,11) = renv(1)
         bcm(ip,12) = epoch(1)
         bcm(ip,13) = ospin(1)
         bcm(ip,14) = 0.0
         bcm(ip,15) = rad(1)/rol(1)
         bcm(ip,16) = float(kstar(2))
         bcm(ip,17) = mass0(2)
         bcm(ip,18) = mass(2)
         bcm(ip,19) = log10(lumin(2))
         bcm(ip,20) = log10(rad(2))
         teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                    (rad(2)**2.d0))**(1.d0/4.d0))
         bcm(ip,21) = log10(teff2)
         bcm(ip,22) = massc(2)
         bcm(ip,23) = radc(2)
         bcm(ip,24) = menv(2)
         bcm(ip,25) = renv(2)
         bcm(ip,26) = epoch(2)
         bcm(ip,27) = ospin(2)
         bcm(ip,28) = 0.0
         bcm(ip,29) = rad(2)/rol(2)
         bcm(ip,30) = tb
         bcm(ip,31) = sep
         bcm(ip,32) = ecc
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Eddington limit for accretion on to the secondary in one orbit.              *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 8    dme = 2.08d-03*eddfac*(1.d0/(1.d0 + zpars(11)))*rad(j2)*tb
      supedd = .false.
      novae = .false.
      disk = .false.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Determine whether the transferred material forms an accretion                *
* disk around the secondary or hits the secondary in a direct                  *
* stream, by using eq.(1) of Ulrich & Burger (1976, ApJ, 206, 509)             *
* fitted to the calculations of Lubow & Shu (1974, ApJ, 198, 383).             *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*     if(kstar(j2).ge.10) disk = .true.
      rmin = 0.0425d0*sep*(q(j2)*(1.d0+q(j2)))**(1.d0/4.d0)
      if(rmin.gt.rad(j2)) disk = .true.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Kelvin-Helmholtz time from the modified classical expression.                *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do 13 , k = 1,2
         tkh(k) = 1.0d+07*mass(k)/(rad(k)*lumin(k))
         if(kstar(k).le.1.or.kstar(k).eq.7.or.kstar(k).ge.10)then
            tkh(k) = tkh(k)*mass(k)
         else
            tkh(k) = tkh(k)*(mass(k) - massc(k))
         endif
 13   continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Dynamical timescale for the primary.                                         *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      tdyn = 5.05d-05*SQRT(rad(j1)**3/mass(j1))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Identify special cases.                                                      *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(kstar(j1).eq.2)then
         qc = 4.d0
      elseif(kstar(j1).eq.3.or.kstar(j1).eq.5.or.kstar(j1).eq.6)then
         qc = (1.67d0-zpars(7)+2.d0*(massc(j1)/mass(j1))**5)/2.13d0 !Kamlah 10.02.2021
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Alternatively use condition of Hjellming & Webbink, 1987, ApJ, 318, 794.     *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*         qc = 0.362 + 1.0/(3.0*(1.0 - massc(j1)/mass(j1)))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Or allow all cases to avoid common-envelope.                                 *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*        qc = 100.d0
      elseif(kstar(j1).eq.8.or.kstar(j1).eq.9)then
         qc = 0.784d0
      else
         qc = 3.d0
      endif
*
      if(kstar(j1).eq.0.and.q(j1).gt.0.695d0)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* This will be dynamical mass transfer of a similar nature to                  *
* common-envelope evolution.  The result is always a single                    *
* star placed in *2.                                                           *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         taum = SQRT(tkh(j1)*tdyn)
         dm1 = mass(j1)
         if(kstar(j2).le.1)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Restrict accretion to thermal timescale of secondary.                        *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            dm2 = taum/tkh(j2)*dm1
            mass(j2) = mass(j2) + dm2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Rejuvenate if the star is still on the main sequence.                        *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            mass0(j2) = mass(j2)
            CALL star(kstar(j2),mass0(j2),mass(j2),tmsnew,tn,
     &                tscls,lums,GB,zpars)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* If the star has no convective core then the effective age decreases,         *
* otherwise it will become younger still.                                      *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            if(mass(j2).lt.0.35d0.or.mass(j2).gt.1.25d0)then
               aj(j2) = tmsnew/tms(j2)*aj(j2)*(mass(j2) - dm2)/mass(j2)
            else
               aj(j2) = tmsnew/tms(j2)*aj(j2)
            endif
            epoch(j2) = tphys - aj(j2)
         elseif(kstar(j2).le.6)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Add all the material to the giant's envelope.                                *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            dm2 = dm1
            mass(j2) = mass(j2) + dm2
            if(kstar(j2).eq.2)then
               mass0(j2) = mass(j2)
               CALL star(kstar(j2),mass0(j2),mass(j2),tmsnew,tn,tscls,
     &                   lums,GB,zpars)
               aj(j2) = tmsnew + tscls(1)*(aj(j2)-tms(j2))/tbgb(j2)
               epoch(j2) = tphys - aj(j2)
            endif
         elseif(kstar(j2).le.12)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Form a new giant envelope.                                                   *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            dm2 = dm1
            kst = ktype(kstar(j1),kstar(j2))
            if(kst.gt.100) kst = kst - 100
            if(kst.eq.4)then
               aj(j2) = aj(j2)/tms(j2)
               massc(j2) = mass(j2)
            endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Check for planets or low-mass WDs.                                           *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            if((kstar(j2).eq.10.and.mass(j2).lt.0.05d0).or.
     &         (kstar(j2).ge.11.and.mass(j2).lt.0.5d0))then
               kst = kstar(j1)
               mass(j1) = mass(j2) + dm2
               mass(j2) = 0.d0
            else
               mass(j2) = mass(j2) + dm2
               CALL gntage(massc(j2),mass(j2),kst,zpars,
     &                     mass0(j2),aj(j2))
               epoch(j2) = tphys - aj(j2)
            endif
            kstar(j2) = kst
         else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* The neutron star or black hole simply accretes at the Eddington rate.        *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            dm2 = MIN(dme*taum/tb,dm1)
            if(dm2.lt.dm1) supedd = .true. 
            mass(j2) = mass(j2) + dm2
         endif
         coel = .true.
         if(mass(j2).gt.0.d0)then
            mass(j1) = 0.d0
            kstar(j1) = 15
         else
            kstar(j1) = kstar(j2)
            kstar(j2) = 15
            mass(j2) = 0.d0   !Kamlah addition 10.02.2021
         endif
         goto 135
      elseif(((ABS(ABS(2*kstar(j1)-11)-3).eq.2.or.kstar(j1).eq.9).
     &        and.(q(j1).gt.qc.or.radx(j1).le.radc(j1))).or.
     &        (kstar(j1).eq.2.and.q(j1).gt.qc).or.
     &        (kstar(j1).eq.4.and.q(j1).gt.qc))then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Common-envelope evolution.                                                   *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         m1ce = mass(j1)
         m2ce = mass(j2)
         CALL comenv(mass0(j1),mass(j1),massc(j1),aj(j1),jspin(j1),
     &               kstar(j1),mass0(j2),mass(j2),massc(j2),aj(j2),
     &               jspin(j2),kstar(j2),zpars,ecc,sep,jorb,coel)
*
         jp = MIN(80,jp + 1)
         bpp(jp,1) = tphys
         bpp(jp,2) = mass(1)
         if(kstar(1).eq.15) bpp(jp,2) = mass0(1)
         bpp(jp,3) = mass(2)
         if(kstar(2).eq.15) bpp(jp,3) = mass0(2)
         bpp(jp,4) = float(kstar(1))
         bpp(jp,5) = float(kstar(2))
         bpp(jp,6) = sep
         bpp(jp,7) = ecc
         bpp(jp,8) = rad(1)/rol(1)
         bpp(jp,9) = rad(2)/rol(2)
         bpp(jp,10) = 7.0
*
         epoch(j1) = tphys - aj(j1)
         if(coel)then
            com = .true.
            goto 135
         endif
         epoch(j2) = tphys - aj(j2)
         if(ecc.gt.1.d0)then
            if(kstar(1).ge.13)then
               rc = corerd(kstar(1),mass(1),mass(1),zpars(2))
               ospin(1) = jspin(1)/(k3*rc*rc*mass(1))
            endif
            if(kstar(2).ge.13)then
               rc = corerd(kstar(2),mass(2),mass(2),zpars(2))
               ospin(2) = jspin(2)/(k3*rc*rc*mass(2))
            endif
            goto 135
         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Next step should be made without changing the time.                          *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         dm1 = m1ce - mass(j1)
         dm2 = mass(j2) - m2ce
         dm22 = dm2
         dtm = 0.d0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Reset orbital parameters as separation may have changed.                     *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         tb = (sep/aursun)*SQRT(sep/(aursun*(mass(1)+mass(2))))
         oorb = twopi/tb
      elseif(kstar(j1).ge.10.and.kstar(j1).le.12.and.
     &       q(j1).gt.0.628d0)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Dynamic transfer from a white dwarf.  Secondary will have KW > 9.            *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         taum = SQRT(tkh(j1)*tdyn)
         dm1 = mass(j1)
         if(eddfac.lt.10.d0)then
            dm2 = MIN(dme*taum/tb,dm1)
            if(dm2.lt.dm1) supedd = .true. 
         else
            dm2 = dm1
         endif
         mass(j2) = mass(j2) + dm2
*
         if(kstar(j1).eq.10.and.kstar(j2).eq.10)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Assume the energy released by ignition of the triple-alpha reaction          *
* is enough to destroy the star.                                               *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            kstar(j2) = 15
            mass(j2) = 0.d0
         elseif(kstar(j1).eq.10.or.kstar(j2).eq.10)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Should be helium overflowing onto a CO or ONe core in which case the         *
* helium swells up to form a giant envelope so a HeGB star is formed.          *
* Allowance for the rare case of CO or ONe flowing onto He is made.            *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            kst = 9
            if(kstar(j2).eq.10) massc(j2) = dm2
            CALL gntage(massc(j2),mass(j2),kst,zpars,mass0(j2),aj(j2))
            kstar(j2) = kst
            epoch(j2) = tphys - aj(j2)
         elseif(kstar(j2).le.12)then
            mass0(j2) = mass(j2)
            if(kstar(j1).eq.12.and.kstar(j2).eq.11)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Mixture of ONe and CO will result in an ONe product.                         *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               kstar(j2) = 12
            endif
         endif
         kstar(j1) = 15
         mass(j1) = 0.d0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Might be a supernova that destroys the system.                               *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(kstar(j2).le.11.and.mass(j2).gt.mch)then
            kstar(j2) = 15
            mass(j2) = 0.d0
         endif
         coel = .true.
         goto 135
      elseif(kstar(j1).eq.13)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Gamma ray burster?                                                           *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         dm1 = mass(j1)
         mass(j1) = 0.d0
         kstar(j1) = 15
         dm2 = dm1
         mass(j2) = mass(j2) + dm2
*
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Gamma ray burster? ---> Kamlah 21.02.2021 addition                           *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(mass(j2).gt.mxns)then
           kstar(j2) = 14
         else
	   kstar(j2) = 13
         endif
*
         kstar(j2) = 14
         coel = .true.
         goto 135
      elseif(kstar(j1).eq.14)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Both stars are black holes.  Let them merge quietly.                         *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         dm1 = mass(j1)
         mass(j1) = 0.d0
         kstar(j1) = 15
         dm2 = dm1
         mass(j2) = mass(j2) + dm2
         coel = .true.
         goto 135
      else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Mass transfer in one Kepler orbit.                                           *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         dm1 = 3.0d-06*tb*(LOG(rad(j1)/rol(j1))**3)*
     &         MIN(mass(j1),5.d0)**2
         if(kstar(j1).eq.2)then
            mew = (mass(j1) - massc(j1))/mass(j1)
            dm1 = MAX(mew,0.01d0)*dm1
         elseif(kstar(j1).ge.10)then
*           dm1 = dm1*1.0d+03/MAX(rad(j1),1.0d-04)
            dm1 = dm1*1.0d+03*mass(j1)/MAX(rad(j1),1.0d-04)
         endif
         kst = kstar(j2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Possibly mass transfer needs to be reduced if primary is rotating            *
* faster than the orbit (not currently implemented).                           *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*        spnfac = MIN(3.d0,MAX(ospin(j1)/oorb,1.d0))
*        dm1 = dm1/spnfac**2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Limit mass transfer to the thermal rate for remaining giant-like stars       *
* and to the dynamical rate for all others.                                    *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(kstar(j1).ge.2.and.kstar(j1).le.9.and.kstar(j1).ne.7)then
***
* JH_temp ... this may be good for HG RLOF??
*           if(kstar(j1).eq.2)then
*              mew = rad(j1)/rol(j1) - 1.d0
*              mew = 2.d0*mew
*              dm1 = dm1*10.d0**mew
*           endif
***
            dm1 = MIN(dm1,mass(j1)*tb/tkh(j1))
         elseif(rad(j1).gt.10.d0*rol(j1).or.(kstar(j1).le.1.and.
     &          kstar(j2).le.1.and.q(j1).gt.qc))then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Allow the stars to merge with the product in *1.                             *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            m1ce = mass(j1)
            m2ce = mass(j2)
            CALL mix(mass0,mass,aj,kstar,zpars)
            dm1 = m1ce - mass(j1)
            dm2 = mass(j2) - m2ce
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Next step should be made without changing the time.                          *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            dtm = 0.d0
            epoch(1) = tphys - aj(1)
            coel = .true.
            goto 135
         else
            dm1 = MIN(dm1,mass(j1)*tb/tdyn)
         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate wind mass loss from the stars during one orbit.                    *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         vorb2 = acc1*(mass(1)+mass(2))/sep
         ivsqm = 1.d0/SQRT(1.d0-ecc*ecc)
         do 14 , k = 1,2
            if(neta.gt.tiny)then
               rlperi = rol(k)*(1.d0-ecc)
               dmr(k) = mlwind(kstar(k),lumin(k),radx(k),
     &                         mass(k),massc(k),rlperi,z)
               vwind2 = 2.d0*beta*acc1*mass(k)/radx(k)
               omv2 = (1.d0 + vorb2/vwind2)**(3.d0/2.d0)
               dmt(3-k) = ivsqm*acc2*dmr(k)*((acc1*mass(3-k)/vwind2)**2)
     &                    /(2.d0*sep*sep*omv2)
               dmt(3-k) = MIN(dmt(3-k),dmr(k))
            else
               dmr(k) = 0.d0
               dmt(3-k) = 0.d0
            endif
 14      continue
*
         do 15 , k = 1,2
            dms(k) = (dmr(k)-dmt(k))*tb
 15      continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Increase time-scale to relative mass loss of 0.5% but not more than twice.   *
* KM is the number of orbits for the timestep.                                 *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         km = MIN(2.d0*km0,5.0d-03/
     &            MAX(ABS(dm1+dms(j1))/mass(j1),dms(j2)/mass(j2)))
         km0 = km
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Modify time-step & mass loss terms by speed-up factor.                       *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         dt = km*tb
         dtm = dt/1.0d+06
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Take the stellar evolution timestep into account but don't let it            *
* be overly restrictive for long lived phases.                                 *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(iter.le.1000) dtm = MIN(dtm,dtmi(1),dtmi(2)) 
         dtm = MIN(dtm,tsave-tphys)
         dt = dtm*1.0d+06
         km = dt/tb
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Decide between accreted mass by secondary and/or system mass loss.           *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         taum = mass(j2)/dm1*tb
         if(kstar(j2).le.2.or.kstar(j2).eq.4)then 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Limit according to the thermal timescale of the secondary.                   *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            dm2 = MIN(1.d0,10.d0*taum/tkh(j2))*dm1
         elseif(kstar(j2).ge.7.and.kstar(j2).le.9)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Naked helium star secondary swells up to a core helium burning star          *
* or SAGB star unless the primary is also a helium star.                       *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            if(kstar(j1).ge.7)then
               dm2 = MIN(1.d0,10.d0*taum/tkh(j2))*dm1
            else
               dm2 = dm1
               dmchk = dm2 - 1.05d0*dms(j2)
               if(dmchk.gt.0.d0.and.dm2/mass(j2).gt.1.0d-04)then
                  kst = MIN(6,2*kstar(j2)-10)
                  if(kst.eq.4)then
                     aj(j2) = aj(j2)/tms(j2)
                     mcx = mass(j2)
                  else
                     mcx = massc(j2)
                  endif
                  mt2 = mass(j2) + km*(dm2 - dms(j2))
                  CALL gntage(mcx,mt2,kst,zpars,mass0(j2),aj(j2))
                  epoch(j2) = tphys + dtm - aj(j2)
*
                  jp = MIN(8000000,jp + 1)
                  bpp(jp,1) = tphys
                  bpp(jp,2) = mass(j1)
                  bpp(jp,3) = mt2
                  bpp(jp,4) = float(kstar(j1))
                  bpp(jp,5) = float(kst)
                  bpp(jp,6) = sep
                  bpp(jp,7) = ecc
                  bpp(jp,8) = rad(1)/rol(1)
                  bpp(jp,9) = rad(2)/rol(2)
                  bpp(jp,10) = 8.0
                  if(j1.eq.2)then
                     bpp(jp,2) = mt2
                     bpp(jp,3) = mass(j1)
                     bpp(jp,4) = float(kst)
                     bpp(jp,5) = float(kstar(j1))
                  endif
*
               endif
            endif            
         elseif(kstar(j1).le.6.and. 
     &           (kstar(j2).ge.10.and.kstar(j2).le.12))then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* White dwarf secondary.                                                       *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            if(dm1/tb.lt.2.71d-07)then
               if(dm1/tb.lt.1.03d-07)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Accrete until a nova explosion blows away most of the accreted material.     *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                  novae = .true. 
                  dm2 = MIN(dm1,dme)
                  if(dm2.lt.dm1) supedd = .true. 
                  dm22 = epsnov*dm2
               else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Steady burning at the surface                                                *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                  dm2 = dm1
               endif
            else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Make a new giant envelope.                                                   *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               dm2 = dm1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Check for planets or low-mass WDs.                                           *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               if((kstar(j2).eq.10.and.mass(j2).lt.0.05d0).or.
     &            (kstar(j2).ge.11.and.mass(j2).lt.0.5d0))then
                  kst = kstar(j2)
               else
                  kst = MIN(6,3*kstar(j2)-27)
                  mt2 = mass(j2) + km*(dm2 - dms(j2))
                  CALL gntage(massc(j2),mt2,kst,zpars,mass0(j2),aj(j2))
                  epoch(j2) = tphys + dtm - aj(j2)
*
                  jp = MIN(8000000,jp + 1)
                  bpp(jp,1) = tphys
                  bpp(jp,2) = mass(j1)
                  bpp(jp,3) = mt2
                  bpp(jp,4) = float(kstar(j1))
                  bpp(jp,5) = float(kst)
                  bpp(jp,6) = sep
                  bpp(jp,7) = ecc
                  bpp(jp,8) = rad(1)/rol(1)
                  bpp(jp,9) = rad(2)/rol(2)
                  bpp(jp,10) = 8.0
                  if(j1.eq.2)then
                     bpp(jp,2) = mt2
                     bpp(jp,3) = mass(j1)
                     bpp(jp,4) = float(kst)
                     bpp(jp,5) = float(kstar(j1))
                  endif
*
               endif
*
            endif
         elseif(kstar(j2).ge.10)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Impose the Eddington limit.                                                  *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            dm2 = MIN(dm1,dme)
            if(dm2.lt.dm1) supedd = .true. 
*
         else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* We have a giant whose envelope can absorb any transferred material.          *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            dm2 = dm1
         endif
         if(.not.novae) dm22 = dm2
*
         if(kst.ge.10.and.kst.le.12)then
            mt2 = mass(j2) + km*(dm22 - dms(j2))
            if(kstar(j1).le.10.and.kst.eq.10.and.mt2.ge.0.7d0)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* HeWD can only accrete helium-rich material up to a mass of 0.7 when          *
* it is destroyed in a possible Type 1a SN.                                    *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               mass(j1) = mass(j1) - km*(dm1 + dms(j1))
               if (mass(j1).lt.0.0d0) then     !Kamlah addition - 10.02.2021
                 mass(j1) = 0.0d0
                 kstar(j1) = 15
               endif
               mass(j2) = 0.d0
               kstar(j2) = 15
               goto 135
            elseif(kstar(j1).le.10.and.kst.ge.11)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* CO and ONeWDs accrete helium-rich material until the accumulated             *
* material exceeds a mass of 0.15 when it ignites. For a COWD with             *
* mass less than 0.95 the system will be destroyed as an ELD in a              *
* possible Type 1a SN. COWDs with mass greater than 0.95 and ONeWDs            *
* will survive with all the material converted to ONe (JH 30/09/99).           *
** Now changed to an ELD for all COWDs when 0.15 accreted (JH 11/01/00).       *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               if((mt2-mass0(j2)).ge.0.15d0)then
                  if(kst.eq.11)then
                     mass(j1) = mass(j1) - km*(dm1 + dms(j1))
                     if (mass(j1).lt.0.0d0) then     !Kamlah addition - 10.02.2021
                       mass(j1) = 0.0d0
                       kstar(j1) = 15
                     endif
                     mass(j2) = 0.d0
                     kstar(j2) = 15
                     goto 135
                  endif
                  mass0(j2) = mt2
               endif
            else
               mass0(j2) = mt2
            endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* If the Chandrasekhar limit is exceeded for a white dwarf then destroy        *
* the white dwarf in a supernova. If the WD is ONe then a neutron star         *
* will survive the supernova and we let HRDIAG take care of this when          *
* the stars are next updated.                                                  *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            if(kst.eq.10.or.kst.eq.11)then
               if(mt2.ge.mch)then
                  dm1 = mch - mass(j2) + km*dms(j2)
                  mass(j1) = mass(j1) - dm1 - km*dms(j1)
                  if (mass(j1).lt.0.0d0) then         !Kamlah addition - 10.02.2021
                    mass(j1) = 0.0d0
                    kstar(j1) = 15
                  endif
                  mass(j2) = 0.d0
                  kstar(j2) = 15
                  goto 135
               endif
            endif
         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Modify mass loss terms by speed-up factor.                                   *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         dm1 = km*dm1
         dm2 = km*dm2
         dm22 = km*dm22
         dme = km*dme
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate orbital angular momentum change due to system mass loss.           *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         djorb = ((dmr(1)+q(1)*dmt(1))*mass(2)*mass(2) +
     &            (dmr(2)+q(2)*dmt(2))*mass(1)*mass(1))/
     &           (mass(1)+mass(2))**2
         djorb = djorb*dt
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* For super-Eddington mass transfer rates, for gamma = -2.0,                   *
* and for novae systems, assume that material is lost from                     *
* the system as if a wind from the secondary.                                  *
* If gamma = -1.0 then assume the lost material carries with it                *
* the specific angular momentum of the primary and for all                     * 
* gamma > 0.0 assume that it takes away a fraction gamma of                    *
* the orbital angular momentum.                                                *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(supedd.or.novae.or.gamma.lt.-1.5d0)then
            djorb = djorb + (dm1 - dm22)*mass(j1)*mass(j1)/
     &              (mass(1)+mass(2))**2
         elseif(gamma.ge.0.d0)then
            djorb = djorb + gamma*(dm1 - dm2)
         else
            djorb = djorb + (dm1 - dm2)*mass(j2)*mass(j2)/
     &              (mass(1)+mass(2))**2
         endif
*
         ecc2 = ecc*ecc
         omecc2 = 1.d0 - ecc2
         sqome2 = SQRT(omecc2)
*
         djorb = djorb*sep*sep*sqome2*oorb
         delet = 0.d0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* For very close systems include angular momentum loss mechanisms.             *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(sep.le.10.d0)then
            djgr = 8.315d-10*mass(1)*mass(2)*(mass(1)+mass(2))/
     &             (sep*sep*sep*sep)
            f1 = (19.d0/6.d0) + (121.d0/96.d0)*ecc2
            sqome5 = sqome2**5
            delet1 = djgr*ecc*f1/sqome5
            djgr = djgr*jorb*(1.d0+0.875d0*ecc2)/sqome5
            djorb = djorb + djgr*dt
            delet = delet + delet1*dt
         endif
*
         do 602 , k = 1,2
*
            dms(k) = km*dms(k)
            if(kstar(k).lt.10) dms(k) = MIN(dms(k),mass(k) - massc(k))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate change in the intrinsic spin of the star.                          *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            djspint(k) = (2.d0/3.d0)*(dmr(k)*radx(k)*radx(k)*ospin(k) -
     &                   xi*dmt(k)*radx(3-k)*radx(3-k)*ospin(3-k))
            djspint(k) = djspint(k)*dt
*
            if(mass(k).gt.0.35d0.and.kstar(k).lt.10)then
               djmb = 5.83d-16*menv(k)*(rad(k)*ospin(k))**3/mass(k)
               djspint(k) = djspint(k) + djmb*dt
            endif
*
 602     continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Adjust the spin angular momentum of each star owing to mass transfer         *
* and conserve total angular momentum.                                         *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         djt = dm1*radx(j1)*radx(j1)*ospin(j1)
         djspint(j1) = djspint(j1) + djt
         djorb = djorb - djt
         if(disk)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Alter spin of the degenerate secondary by assuming that material             *
* falls onto the star from the inner edge of a Keplerian accretion             *
* disk and that the system is in a steady state.                               *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            djt = dm2*twopi*aursun*SQRT(aursun*mass(j2)*radx(j2)) 
            djspint(j2) = djspint(j2) - djt 
            djorb = djorb + djt
*
         else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* No accretion disk.                                                           *
* Calculate the angular momentum of the transferred material by                *
* using the radius of the disk (see Ulrich & Burger) that would                *
* have formed if allowed.                                                      *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            rdisk = 1.7d0*rmin
            djt = dm2*twopi*aursun*SQRT(aursun*mass(j2)*rdisk) 
            djspint(j2) = djspint(j2) - djt
            djorb = djorb + djt
*
         endif
         djtx(2) = djt
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Adjust the secondary spin if a nova eruption has occurred.                   *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(novae)then
            djt = (dm2 - dm22)*radx(j2)*radx(j2)*ospin(j2) 
            djspint(j2) = djspint(j2) + djt 
            djtx(2) = djtx(2) - djt
         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Calculate circularization, orbital shrinkage and spin up.                    *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         do 603 , k = 1,2
*
            dspint(k) = 0.d0
            if(((kstar(k).le.9.and.rad(k).ge.0.01d0*rol(k)).or.
     &         (kstar(k).ge.10.and.k.eq.j1)).and.tflag.gt.0)then
*
               raa2 = (radx(k)/sep)**2
               raa6 = raa2**3
*
               f5 = 1.d0+ecc2*(3.d0+ecc2*0.375d0)
               f4 = 1.d0+ecc2*(1.5d0+ecc2*0.125d0)
               f3 = 1.d0+ecc2*(3.75d0+ecc2*(1.875d0+ecc2*7.8125d-02))
               f2 = 1.d0+ecc2*(7.5d0+ecc2*(5.625d0+ecc2*0.3125d0))
               f1 = 1.d0+ecc2*(15.5d0+ecc2*(31.875d0+ecc2*(11.5625d0
     &                  +ecc2*0.390625d0)))
*
               if((kstar(k).eq.1.and.mass(k).ge.1.25d0).or.
     &            kstar(k).eq.4.or.kstar(k).eq.7)then
                  tc = 1.592d-09*(mass(k)**2.84d0)
                  f = 1.9782d+04*SQRT((mass(k)*radx(k)*radx(k))/sep**5)*
     &                tc*(1.d0+q(3-k))**(5.d0/6.d0)
                  tcqr = f*q(3-k)*raa6
                  rg2 = k2str(k)
               elseif(kstar(k).le.9)then
                  renv(k) = MIN(renv(k),radx(k)-radc(k))
                  renv(k) = MAX(renv(k),1.0d-10)
                  tc = mr23yr*(menv(k)*renv(k)*(radx(k)-0.5d0*renv(k))/
     &                 (3.d0*lumin(k)))**(1.d0/3.d0)
                  ttid = twopi/(1.0d-10 + ABS(oorb - ospin(k)))
                  f = MIN(1.d0,(ttid/(2.d0*tc)**2))
                  tcqr = 2.d0*f*q(3-k)*raa6*menv(k)/(21.d0*tc*mass(k))
                  rg2 = (k2str(k)*(mass(k)-massc(k)))/mass(k)
               else
                  f = 7.33d-09*(lumin(k)/mass(k))**(5.d0/7.d0)
                  tcqr = f*q(3-k)*q(3-k)*raa2*raa2/(1.d0+q(3-k))
                  rg2 = k3
               endif
               sqome3 = sqome2**3
               delet1 = 27.d0*tcqr*(1.d0+q(3-k))*raa2*(ecc/sqome2**13)*
     &                  (f3 - (11.d0/18.d0)*sqome3*f4*ospin(k)/oorb)
               tcirc = ecc/(ABS(delet1) + 1.0d-20)
               delet = delet + delet1*dt
               dspint(k) = (3.d0*q(3-k)*tcqr/(rg2*omecc2**6))*
     &                     (f2*oorb - sqome3*f5*ospin(k))
               eqspin = oorb*f2/(sqome3*f5)
               if(dt.gt.0.d0)then
                  if(dspint(k).ge.0.d0)then
                     dspint(k) = MIN(dt*dspint(k),eqspin-ospin(k))/dt
                  else
                     dspint(k) = MAX(dt*dspint(k),eqspin-ospin(k))/dt
                  endif
               endif
               djt = (k2str(k)*(mass(k)-massc(k))*radx(k)*radx(k) +
     &                k3*massc(k)*radc(k)*radc(k))*dspint(k)
               djorb = djorb + djt*dt
               djspint(k) = djspint(k) - djt*dt
*
            endif
*
            jspin(k) = MAX(1.0d-10,jspin(k) - djspint(k))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Ensure that the star does not spin up beyond break-up, and transfer          *
* the excess angular momentum back to the orbit.                               *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            ospbru = twopi*SQRT(mass(k)*aursun**3/radx(k)**3)
            jspbru = (k2str(k)*(mass(k)-massc(k))*radx(k)*radx(k) +
     &                k3*massc(k)*radc(k)*radc(k))*ospbru
            if(jspin(k).gt.jspbru)then
               mew = 1.d0
               if(djtx(2).gt.0.d0)then
                  mew = MIN(mew,(jspin(k) - jspbru)/djtx(2))
               endif
               djorb = djorb - (jspin(k) - jspbru)
               jspin(k) = jspbru
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* If excess material should not be accreted, activate next line.               *
*              dm22 = (1.d0 - mew)*dm22                                        *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            endif
*
 603     continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Update the masses.                                                           *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         kstar(j2) = kst
         mass(j1) = mass(j1) - dm1 - dms(j1)
         if(kstar(j1).le.1.or.kstar(j1).eq.7) mass0(j1) = mass(j1)
         mass(j2) = mass(j2) + dm22 - dms(j2)
         if(kstar(j2).le.1.or.kstar(j2).eq.7) mass0(j2) = mass(j2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* For a HG star check if the initial mass can be reduced.                      *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(kstar(j1).eq.2.and.mass0(j1).le.zpars(3))then
            m0 = mass0(j1)
            mass0(j1) = mass(j1)
            CALL star(kstar(j1),mass0(j1),mass(j1),tmsnew,tn,tscls,
     &                lums,GB,zpars)
            if(GB(9).lt.massc(j1))then
               mass0(j1) = m0
            endif
         endif
         if(kstar(j2).eq.2.and.mass0(j2).le.zpars(3))then
            m0 = mass0(j2)
            mass0(j2) = mass(j2)
            CALL star(kstar(j2),mass0(j2),mass(j2),tmsnew,tn,tscls,
     &                lums,GB,zpars)
            if(GB(9).lt.massc(j2))then
               mass0(j2) = m0
            endif
         endif
*
         ecc = ecc - delet
         ecc = MAX(ecc,0.d0)
         if(ecc.lt.1.0d-10) ecc = 0.d0
*
         if(ecc.ge.1.d0) goto 135
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Ensure that Jorb does not become negative which could happen if the          *
* primary overfills its Roche lobe initially. In this case we simply           *
* allow contact to occur.                                                      *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         jorb = MAX(1.d0,jorb - djorb)
         sep = (mass(1) + mass(2))*jorb*jorb/
     &         ((mass(1)*mass(2)*twopi)**2*aursun**3*(1.d0-ecc*ecc))
         tb = (sep/aursun)*SQRT(sep/(aursun*(mass(1)+mass(2))))
         oorb = twopi/tb
*
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Always rejuvenate the secondary and age the primary if they are on           *
* the main sequence.                                                           *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(kstar(j1).le.2.or.kstar(j1).eq.7)then
         CALL star(kstar(j1),mass0(j1),mass(j1),tmsnew,tn,tscls,
     &             lums,GB,zpars)
         if(kstar(j1).eq.2)then
            aj(j1) = tmsnew + (tscls(1) - tmsnew)*(aj(j1)-tms(j1))/
     &                        (tbgb(j1) - tms(j1))
         else
            aj(j1) = tmsnew/tms(j1)*aj(j1)
         endif
         epoch(j1) = tphys - aj(j1)
      endif
*
      if(kstar(j2).le.2.or.kstar(j2).eq.7)then
         CALL star(kstar(j2),mass0(j2),mass(j2),tmsnew,tn,tscls,
     &             lums,GB,zpars)
         if(kstar(j2).eq.2)then
            aj(j2) = tmsnew + (tscls(1) - tmsnew)*(aj(j2)-tms(j2))/
     &                        (tbgb(j2) - tms(j2))
         elseif((mass(j2).lt.0.35d0.or.mass(j2).gt.1.25d0).
     &           and.kstar(j2).ne.7)then
            aj(j2) = tmsnew/tms(j2)*aj(j2)*(mass(j2) - dm22)/mass(j2)
         else
            aj(j2) = tmsnew/tms(j2)*aj(j2)
         endif
         epoch(j2) = tphys - aj(j2)
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Obtain the stellar parameters for the next step.                             *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      tphys = tphys + dtm
      do 90 , k = 1,2
         age = tphys - epoch(k)
         m0 = mass0(k)
         mt = mass(k)
         mc = massc(k)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Masses over 100Msun should probably not be trusted in the                    *
* evolution formulae. ! Commented out - Albrecht 12.08.2020                    *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*         if(mt.gt.100.d0)then
*            WRITE(99,*)' MASS EXCEEDED ',mass1i,mass2i,tbi,ecci,mt
*            goto 140
*         endif
         kw = kstar(k)
         CALL star(kw,m0,mt,tm,tn,tscls,lums,GB,zpars)
         CALL hrdiag(m0,age,mt,tm,tn,tscls,lums,GB,zpars,
     &               rm,lum,kw,mc,rc,me,re,k2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Check for a supernova and correct the semi-major axis if so.                 *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*************************************************************************
*   
*   SN kicks not allowed for accreted WD and NS COMMENTED OUT!!! 
*   KAMLAH - 10.02.2021
*         if(kw.ne.kstar(k).and.kstar(k).le.12.and.
*     &      (kw.eq.13.or.kw.eq.14))then
*
         if(kw.ne.kstar(k).and.kstar(k).le.9.and.
     &      (kw.gt.9.and.kw.lt.15))then
            dms(k) = mass(k) - mt
            CALL kick(kw,mass(k),mt,mass(3-k),ecc,sep,jorb,vs)
            vkick(k) = dsqrt(vs(1)*vs(1)+vs(2)*vs(2)+vs(3)*vs(3))
            if(k.eq.1) then
                vs1(1)=vs(1)
                vs1(2)=vs(2)
                vs1(3)=vs(3)
            elseif(k.eq.2.) then
                vs2(1)=vs(1)
                vs2(2)=vs(2)
                vs2(3)=vs(3)
            endif	  
            if(ecc.gt.1.d0)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
** ADDITIONS in accordance with Nbody6++GPU - Albrecht - 14.08.2020 
******* BH Kerr Metric spin parameter *****
              if(KW.eq.14) then
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
                   b_const = 3.666D0
                   elseif (z.ge.0.001D0.and.z.lt.0.004D0) then
                   alow = 0.0D0
                   mtwo = 27.7D0
                   mone = 18.0D0
                   aconst = -0.088D0
                   b_const = 2.434D0
                   elseif (z.ge.0.004D0.and.z.lt.0.01D0) then
                   alow = 0.25D0
                   mtwo = 37.8D0
                   mone = 31.0D0
                  aconst = -0.088D0
                   b_const = 3.578D0
                   else
                   alow = 0.13D0
                   mtwo = 24.2D0
                   mone = 16.0D0
                   aconst = -0.088D0
                   b_const = 2.258D0
                   endif
                   if (MCO.le.mone) then
                      aspin = 0.85D0
                   elseif (MCO.gt.mone.and.MCO.lt.mtwo) then
                      aspin = (aconst*MCO) + b_const
                   else
                     aspin = alow
                  endif
 		          if (aspin.lt.0.0D0) then
                   aspin = 0.0D0
				  endif
			 endif
*********
              if(BHSPIN.EQ.2) then
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
			endif 
                if(BHSPIN.EQ.0) then
                WRITE(*,*) 'Using bhspin = 0 - Fuller BH spins'	
*		   
* ------------ Zero BH natal spins --------------------------------------
*
                aspin = 0.0D0
                endif    
*		   
* ------------ gm -------------------------------------------------------
*
               mbh = mt
			   WRITE(*,*) 'MBH', mbh
               mbh = mt*M_sun   ! Doublecheck mt!
			   WRITE(*,*) 'MBH in g', mbh
*		   
* ------------ gm cm^2 / s ----------------------------------------------
*		   
               jspin(k) = (G_grav*aspin*mbh**2)/cspeed
*			   WRITE(*,*) 'jspin: (G*aspin*mbh**2)/cspeed (gm cm^2 / s) = ', jspin			   
*		   
* ------------ Msun Rsun^2 / year ---------------------------------------
*
               jspin(k) = jspin(k)*(yearsc/(M_sun*(rsunkm**2)))
*			   WRITE(*,*) 'jspin: (yearsc/(M_sun*(rsunkm**2))) (Msun Rsun^2 / year) =', jspin		
		       ospin(k) = jspin(k)/(k3*rc*rc*mc)
*			   WRITE(*,*) 'Spin after kick (year^-1)', ospin(k)
*		   
* ------------ Scaled ---------------------------------------------------
*
		      WRITE(*,*) 'aspin', aspin
              WRITE(*,*) 'Spin after kick', ospin, jspin   ! This does not work properly
			  WRITE(*,*) 'Using BHSPIN =', bhspin	
               kstar(k) = kw
               mass(k) = mt
               epoch(k) = tphys - age
               goto 135
              endif	
              elseif(kw.gt.9.and.kw.lt.14)then		
               ospin(k) = 2.0d+08                !KAMLAH 10.02.2021
               jspin(k) = k3*rc*rc*mc*ospin(k)   !Kamlah 10.02.2021
               kstar(k) = kw
               mass(k) = mt
               epoch(k) = tphys - age
               goto 135
            endif
            tb = (sep/aursun)*SQRT(sep/(aursun*(mt+mass(3-k))))
            oorb = twopi/tb
         endif
         if(kw.ne.kstar(k))then
            change = .true.
            if((kw.eq.13.or.kw.eq.14).and.kstar(k).le.12)then  !MIGHT NEED TO CHANGE
               snova = .true.
            endif
            mass(k) = mt
            if(kw.eq.15)then
               kstar(k) = kw
               goto 135
            endif
            mass0(k) = m0
            epoch(k) = tphys - age
            if(kw.gt.9.and.kw.lt.14)then         !KAMLAH 23.02.2021
               ospin(k) = 2.0d+08
               jspin(k) = k3*rc*rc*mc*ospin(k)
            endif
         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Determine stellar evolution timescale for nuclear burning types.             *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(kw.le.9)then
            CALL deltat(kw,age,tm,tn,tscls,dt,dtr)
            dtmi(k) = MIN(dt,dtr)
*           dtmi(k) = dtr
            dtmi(k) = MAX(1.0d-07,dtmi(k))
         else
            dtmi(k) = 1.0d+10
         endif
*        dtmi(k) = MAX((tn-age),1.0d-07)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Save relevent solar quantities.                                              *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         aj(k) = age
         kstar(k) = kw
         rad(k) = rm
         radx(k) = rm
         lumin(k) = lum
         massc(k) = mc
         radc(k) = rc
         menv(k) = me
         renv(k) = re
         k2str(k) = k2
         tms(k) = tm
         tbgb(k) = tscls(1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Check for blue straggler formation.                                          *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(kw.le.1.and.tm.lt.tphys.and..not.bss)then
            bss = .true.
            jp = MIN(8000000,jp + 1)
            bpp(jp,1) = tphys
            bpp(jp,2) = mass(1)
            bpp(jp,3) = mass(2)
            bpp(jp,4) = float(kstar(1))
            bpp(jp,5) = float(kstar(2))
            bpp(jp,6) = sep
            bpp(jp,7) = ecc
            bpp(jp,8) = rad(1)/rol(1)
            bpp(jp,9) = rad(2)/rol(2)
            bpp(jp,10) = 14.0
         endif
*
 90   continue
*
      do 100 , k = 1,2
         q(k) = mass(k)/mass(3-k)
         rol(k) = rl(q(k))*sep
 100  continue
      if(rad(j1).gt.rol(j1)) radx(j1) = MAX(radc(j1),rol(j1))
      do 110 , k = 1,2
         ospin(k) = jspin(k)/(k2str(k)*(mass(k)-massc(k))*radx(k)*
     &              radx(k) + k3*massc(k)*radc(k)*radc(k))
 110  continue
*
      if((isave.and.tphys.ge.tsave).or.iplot)then
         ip = ip + 1
         bcm(ip,1) = tphys
         bcm(ip,2) = float(kstar(1))
         bcm(ip,3) = mass0(1)
         bcm(ip,4) = mass(1)
         bcm(ip,5) = log10(lumin(1))
         bcm(ip,6) = log10(rad(1))
         teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                    (rad(1)**2.d0))**(1.d0/4.d0))
         bcm(ip,7) = log10(teff1)
         bcm(ip,8) = massc(1)
         bcm(ip,9) = radc(1)
         bcm(ip,10) = menv(1)
         bcm(ip,11) = renv(1)
         bcm(ip,12) = epoch(1)
         bcm(ip,13) = ospin(1)
         bcm(ip,15) = rad(1)/rol(1)
         bcm(ip,16) = float(kstar(2))
         bcm(ip,17) = mass0(2)
         bcm(ip,18) = mass(2)
         bcm(ip,19) = log10(lumin(2))
         bcm(ip,20) = log10(rad(2))
         teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                    (rad(2)**2.d0))**(1.d0/4.d0))
         bcm(ip,21) = log10(teff2)
         bcm(ip,22) = massc(2)
         bcm(ip,23) = radc(2)
         bcm(ip,24) = menv(2)
         bcm(ip,25) = renv(2)
         bcm(ip,26) = epoch(2)
         bcm(ip,27) = ospin(2)
         bcm(ip,29) = rad(2)/rol(2)
         bcm(ip,30) = tb
         bcm(ip,31) = sep
         bcm(ip,32) = ecc
         dt = MAX(dtm,1.0d-12)*1.0d+06
         if(j1.eq.1)then
            bcm(ip,14) = (-1.0*dm1 - dms(1))/dt
            bcm(ip,28) = (dm2 - dms(2))/dt
         else
            bcm(ip,14) = (dm2 - dms(1))/dt
            bcm(ip,28) = (-1.0*dm1 - dms(2))/dt
         endif
         if(isave) tsave = tsave + dtp
      endif
*
      if(tphys.ge.tphysf) goto 140
*
      if(change)then
         change = .false.
         jp = MIN(8000000,jp + 1)
         bpp(jp,1) = tphys
         bpp(jp,2) = mass(1)
         bpp(jp,3) = mass(2)
         bpp(jp,4) = float(kstar(1))
         bpp(jp,5) = float(kstar(2))
         bpp(jp,6) = sep
         bpp(jp,7) = ecc
         bpp(jp,8) = rad(1)/rol(1)
         bpp(jp,9) = rad(2)/rol(2)
         bpp(jp,10) = 2.0
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Test whether the primary still fills its Roche lobe.                         *   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(rad(j1).gt.rol(j1).and..not.snova)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Test for a contact system                                                    *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(rad(j2).gt.rol(j2)) goto 130
         iter = iter + 1
         goto 8
      else
         jp = MIN(8000000,jp + 1)
         bpp(jp,1) = tphys
         bpp(jp,2) = mass(1)
         bpp(jp,3) = mass(2)
         bpp(jp,4) = float(kstar(1))
         bpp(jp,5) = float(kstar(2))
         bpp(jp,6) = sep
         bpp(jp,7) = ecc
         bpp(jp,8) = rad(1)/rol(1)
         bpp(jp,9) = rad(2)/rol(2)
         bpp(jp,10) = 4.0
         dtm = 0.d0
         goto 4
      endif
*
 130  continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Contact system.                                                              *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      coel = .true.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* If *1 or *2 is giant-like this will be common-envelope evolution.            *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      m1ce = mass(j1)
      m2ce = mass(j2)
      rrl1 = MIN(999.999d0,rad(1)/rol(1))
      rrl2 = MIN(999.999d0,rad(2)/rol(2))
*
      jp = MIN(8000000,jp + 1)
      bpp(jp,1) = tphys
      bpp(jp,2) = mass(1)
      bpp(jp,3) = mass(2)
      bpp(jp,4) = float(kstar(1))
      bpp(jp,5) = float(kstar(2))
      bpp(jp,6) = sep
      bpp(jp,7) = ecc
      bpp(jp,8) = rrl1
      bpp(jp,9) = rrl2
      bpp(jp,10) = 5.0
*
      if(kstar(j1).ge.2.and.kstar(j1).le.9.and.kstar(j1).ne.7)then
         CALL comenv(mass0(j1),mass(j1),massc(j1),aj(j1),jspin(j1),
     &               kstar(j1),mass0(j2),mass(j2),massc(j2),aj(j2),
     &               jspin(j2),kstar(j2),zpars,ecc,sep,jorb,coel)
         com = .true.
      elseif(kstar(j2).ge.2.and.kstar(j2).le.9.and.kstar(j2).ne.7)then
         CALL comenv(mass0(j2),mass(j2),massc(j2),aj(j2),jspin(j2),
     &               kstar(j2),mass0(j1),mass(j1),massc(j1),aj(j1),
     &               jspin(j1),kstar(j1),zpars,ecc,sep,jorb,coel)
         com = .true.
      else
         CALL mix(mass0,mass,aj,kstar,zpars)
      endif
      if(com)then
         jp = MIN(8000000,jp + 1)
         bpp(jp,1) = tphys
         bpp(jp,2) = mass(1)
         if(kstar(1).eq.15) bpp(jp,2) = mass0(1)
         bpp(jp,3) = mass(2)
         if(kstar(2).eq.15) bpp(jp,3) = mass0(2)
         bpp(jp,4) = float(kstar(1))
         bpp(jp,5) = float(kstar(2))
         bpp(jp,6) = sep
         bpp(jp,7) = ecc
         rrl1 = MIN(rrl1,0.99d0)
         rrl2 = MIN(rrl2,0.99d0)
         bpp(jp,8) = rrl1
         bpp(jp,9) = rrl2
         bpp(jp,10) = 7.0
      endif
      epoch(1) = tphys - aj(1)
      epoch(2) = tphys - aj(2)
      if(.not.coel)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Next step should be made without changing the time.                          *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         if(ecc.gt.1.d0)then
            if(kstar(1).ge.13)then
               rc = corerd(kstar(1),mass(1),mass(1),zpars(2))
               ospin(1) = jspin(1)/(k3*rc*rc*mass(1))
            endif
            if(kstar(2).ge.13)then
               rc = corerd(kstar(2),mass(2),mass(2),zpars(2))
               ospin(2) = jspin(2)/(k3*rc*rc*mass(2))
            endif
            goto 135
         endif
         dtm = 0.d0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Reset orbital parameters as separation may have changed.                     *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         tb = (sep/aursun)*SQRT(sep/(aursun*(mass(1)+mass(2))))
         oorb = twopi/tb
         goto 4
      endif
*
 135  continue
*
      sgl = .true.
      if(kstar(1).ne.15.or.kstar(2).ne.15)then
         if(com)then
            com = .false.
         else
            jp = MIN(8000000,jp + 1)
            bpp(jp,1) = tphys
            bpp(jp,2) = mass(1)
            if(kstar(1).eq.15) bpp(jp,2) = mass0(1)
            bpp(jp,3) = mass(2)
            if(kstar(2).eq.15) bpp(jp,3) = mass0(2)
            bpp(jp,4) = float(kstar(1))
            bpp(jp,5) = float(kstar(2))
            bpp(jp,6) = zero
            bpp(jp,7) = zero
            bpp(jp,8) = zero
            bpp(jp,9) = ngtv
            if(coel)then
               bpp(jp,10) = 6.0
            elseif(ecc.gt.1.d0)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Binary dissolved by a supernova or tides.                                    *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               bpp(jp,6) = sep
               bpp(jp,7) = ecc
               bpp(jp,9) = ngtv2
               bpp(jp,10) = 11.0
            else
               bpp(jp,10) = 9.0
            endif
         endif
         if(kstar(2).eq.15)then
            kmax = 1
            rol(2) = -1.d0*rad(2)
            dtmi(2) = tphysf
            if(kstar(1).ge.3.and.kstar(1).le.9) ospin(1) = -1.d0 !KAMLAH 10.02.2021
         elseif(kstar(1).eq.15)then
            kmin = 2
            rol(1) = -1.d0*rad(1)
            dtmi(1) = tphysf
            if(kstar(2).ge.3.and.kstar(2).le.9) ospin(2) = -1.d0 !KAMLAH 10.02.2021
         endif
         ecc = -1.d0
         sep = 0.d0
         dtm = 0.d0
         coel = .false.
         goto 4
      endif
*
 140  continue
*
      if(com)then
         com = .false.
      else
         jp = MIN(8000000,jp + 1)
         bpp(jp,1) = tphys
         bpp(jp,2) = mass(1)
         if(kstar(1).eq.15.and.bpp(jp-1,4).lt.15.0)then
            bpp(jp,2) = mass0(1)
         endif
         bpp(jp,3) = mass(2)
         if(kstar(2).eq.15.and.bpp(jp-1,5).lt.15.0)then
            bpp(jp,3) = mass0(2)
         endif
         bpp(jp,4) = float(kstar(1))
         bpp(jp,5) = float(kstar(2))
         bpp(jp,6) = zero
         bpp(jp,7) = zero
         bpp(jp,8) = zero
         if(coel)then
            bpp(jp,9) = ngtv
            bpp(jp,10) = 6.0
         elseif(kstar(1).eq.15.and.kstar(2).eq.15)then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* Cases of accretion induced supernova or single star supernova.               *
* No remnant is left in either case.                                           *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            bpp(jp,9) = ngtv2
            bpp(jp,10) = 9.0
         else
            bpp(jp,6) = sep
            bpp(jp,7) = ecc
            bpp(jp,8) = rad(1)/rol(1)
            bpp(jp,9) = rad(2)/rol(2)
            bpp(jp,10) = 10.0
         endif
      endif
*
      if((isave.and.tphys.ge.tsave).or.iplot)then
         ip = ip + 1
         bcm(ip,1) = tphys
         bcm(ip,2) = float(kstar(1))
         bcm(ip,3) = mass0(1)
         bcm(ip,4) = mass(1)
         bcm(ip,5) = log10(lumin(1))
         bcm(ip,6) = log10(rad(1))
         teff1 = 1000.d0*((1130.d0*lumin(1)/
     &                    (rad(1)**2.d0))**(1.d0/4.d0))
         bcm(ip,7) = log10(teff1)
         bcm(ip,8) = massc(1)
         bcm(ip,9) = radc(1)
         bcm(ip,10) = menv(1)
         bcm(ip,11) = renv(1)
         bcm(ip,12) = epoch(1)
         bcm(ip,13) = ospin(1)
         bcm(ip,15) = rad(1)/rol(1)
         bcm(ip,16) = float(kstar(2))
         bcm(ip,17) = mass0(2)
         bcm(ip,18) = mass(2)
         bcm(ip,19) = log10(lumin(2))
         bcm(ip,20) = log10(rad(2))
         teff2 = 1000.d0*((1130.d0*lumin(2)/
     &                    (rad(2)**2.d0))**(1.d0/4.d0))
         bcm(ip,21) = log10(teff2)
         bcm(ip,22) = massc(2)
         bcm(ip,23) = radc(2)
         bcm(ip,24) = menv(2)
         bcm(ip,25) = renv(2)
         bcm(ip,26) = epoch(2)
         bcm(ip,27) = ospin(2)
         bcm(ip,29) = rad(2)/rol(2)
         bcm(ip,30) = tb
         bcm(ip,31) = sep
         bcm(ip,32) = ecc
         dt = MAX(dtm,1.0d-12)*1.0d+06
         if(j1.eq.1)then
            bcm(ip,14) = (-1.0*dm1 - dms(1))/dt
            bcm(ip,28) = (dm2 - dms(2))/dt
         else
            bcm(ip,14) = (dm2 - dms(1))/dt
            bcm(ip,28) = (-1.0*dm1 - dms(2))/dt
         endif
         if(isave) tsave = tsave + dtp
         if(tphysf.le.0.d0)then
            ip = ip + 1
            do 145 , k = 1,42
               bcm(ip,k) = bcm(ip-1,k)
 145        continue
         endif
      endif
*
      tphysf = tphys
      if(sgl)then
         if(ecc.ge.0.d0.and.ecc.le.1.d0) ecc = -1.d0
         tb = -1.d0
      endif
      tb = tb*yeardy
      if(jp.ge.80)then
         WRITE(99,*)' EVOLV2 ARRAY ERROR ',mass1i,mass2i,tbi,ecci
         WRITE(*,*)' STOP: EVOLV2 ARRAY ERROR '
*         CALL exit(0)
*         STOP
      elseif(jp.ge.40)then
         WRITE(99,*)' EVOLV2 ARRAY WARNING ',mass1i,mass2i,tbi,ecci,jp
      endif
      bcm(ip+1,1) = -1.0
      bpp(jp+1,1) = -1.0
*
      RETURN
      END
***

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
* SOURCES                                                                      *
* Boffin H. M. J. & Jorissen A. (1988, A&A, 205, 155-)                         *
* Campbell C. G. (1984, MNRAS, 207, 433-)                                      *
* Claeys J. S. W., Pols O. R., Izzard R. G. et al., (2014, A&A, 563, A83)      *
* de Kool M. (1990, ApJ, 358, 189-)                                            *
* Hurley J. R., Tout C. A., Pols O. R. et al. (2002, MNRAS, 329, 897-)         *
* Hut P. (1981, A&A, 99, 126-)                                                 *
* Ivanova N., Justham S., Chen X. et al. (2013, A&A Review, 59-)               *
* Peters P. C. & Mathews J. (1963 Phys. Rev., 131, 435-)                       *
* Peters P. C. (1964 Phys. Rev. 136, B1224)                                    *
* Ulrich R. K. & Burger H. (1976, ApJ, 206, 509-)                              *
* Zahn J. P. (1975, A&A, 41, 329-)                                             *
* Zahn J. P. (1977, A&A, 57, 383-)                                             *
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

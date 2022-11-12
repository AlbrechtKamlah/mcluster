      subroutine input
*
*
*       read input parameters
*       ----------------
*
      include 'common.h'
*
      integer i,ixx1,ixx2,ixx3,verboseTmp, tmpInt
*
      real*8 xxx4, tmp, bodyMin, bodyMax
      real*4 xxx32
      CHARACTER(len=32) :: arg
      CHARACTER(len=100) :: key
      CHARACTER(len=1000) :: conf
      character*1 st1
      character*2 st2

      ixx1 = 0
      ixx2 = 0
      ixx3 = 0
      xxx4 = 0.0d0

      conf = "mcluster.ini"


      ! check input parameters
      DO i = 1, iargc()
        CALL getarg(i, arg)
          print*,"INFO: verbose mode enabled"
        if (arg.eq."--help".OR.arg.eq."-h") then
      print*,"MOCCA = MOnte Carlo Cluster simulAtor"
      print*,""
      print*,"  --help    | -h   prints this help"
      print*,"  --conf    | -c   specify configuration file with the "
      print*,"                   initial parameters. If file is not "
      print*,"                   specified then file mocca.ini is read"
      print*,"                   by default. If mocca.ini file is not "
      print*,"                   present then all default values are "
      print*,"                   used"
          stop
        else if (arg.eq."--conf".OR.arg.eq."-c") then
          CALL getarg(i + 1, conf)
          print*,"INFO: reading configuration from ",conf
        endif
!        WRITE (*,*) arg
      END DO

      ! check if there are unrecognized parameters in the mocca.ini file
      call config_validate(conf)
     
      key = "Mcluster:n"
      call config_getstr(npopchar, "100000", key, conf);
      call char_to_arrayint(npopchar, npop);

      numpop = 0
      do i = 1, 10
        numpop = numpop + 1
        ixx2 = ixx2 + npop(i)
        if (npop(i).EQ.0) then
          numpop = numpop - 1
          exit
        endif
      end do

      print*,"Total number of objects ",ixx2

      key = "Mcluster:fracb"
      call config_getstr(fracbchar, "0.2", key, conf);
      call char_to_arraydouble(fracbchar, fracb);
      do i = 1, numpop
      if (fracb(i).lt.0.0d0.OR.fracb(i).gt.1.0d0) then
      print*,"ERROR Mcluster:fracb ",i,"-th is wrong ",
     &        "file. It should be in the range [0,1]"
        stop
      endif
      end do

      key = "Mcluster:fracb_reference"
      call config_getstr(fracbreferencechar, "0.95, 0.95", key, conf);
      call char_to_arraydouble(fracbreferencechar, fracb_reference);
      do i = 1, numpop
      if (fracb_reference(i).lt.0.0d0.OR.
     &  fracb_reference(i).gt.1.0d0) then
      write(*,*)"ERROR Mcluster:fracb_reference ",i,"-th is wrong ",
     &        "file. It should be in the range [0,1]"
        stop
      endif
      end do

      key = "Mcluster:initialModel"
      call config_getstr(initmodelchar, "1", key, conf);
      call char_to_arrayint(initmodelchar, initmodel);
      do i = 1, numpop
      if (initmodel(i).lt.0.OR.initmodel(i).gt.4) then
      print*,"ERROR Mcluster:initialModel ",i,"-th",
     &     " is wrong. It should be 0, 1, 2, 3 or 4"
         stop
      endif 
      end do

      key = "Mcluster:w0"
      call config_getstr(w0char, "1", key, conf);
      call char_to_arraydouble(w0char, w0);
      do i = 1, numpop
       if (w0(i).lt.1.0d0.OR.w0(i).gt.12.0d0) then
      print*,"ERROR Mcluster:w0 ", i, "th",
     &   " is wrong: ", w0(i) , ". It should be between 1.0 and 12.0"
        stop
      endif
      end do
      
      key = "Mcluster:S"
      call config_getstr(Segchar, "0.0", key, conf);
      call char_to_arraydouble(Segchar, Seg);
      do i = 1, numpop
       if (Seg(i).lt.0.0d0.OR.Seg(i).gt.1.0d0) then
      print*,"ERROR Mcluster:S ",i,"th",
     &   " is wrong. It should be between 0 and 1"
         stop
      endif
      end do
       
      key = "Mcluster:fractal"
      call config_getstr(fractalchar, "3.0", key, conf);
      call char_to_arraydouble(fractalchar, fractal);
      do i = 1, numpop
       if (abs(fractal(i)).lt.1.6d0.OR.
     &     abs(fractal(i)).gt.3.0d0) then
      print*,"ERROR Mcluster:fractal ",i,"th",
     &   " is wrong. It should be between 1.6 and 3.0"
         stop
      endif
      end do

      key = "Mcluster:qvir"
      call config_getdouble(qvir, 0.5d0, key, conf);

      key = "Mcluster:mfunc"
      call config_getstr(imfgchar, "1.0", key, conf);
      call char_to_arrayint(imfgchar, imfg);

      key = "Mcluster:single_mass"
      call config_getstr(equalmasschar, "1.0", key, conf);
      call char_to_arraydouble(equalmasschar, equalmass);

      key = "Mcluster:mlow"
      call config_getstr(mlowchar, "0.08", key, conf);
      call char_to_arraydouble(mlowchar, mlow);

      key = "Mcluster:mup"
      call config_getstr(mupchar, "100.0", key, conf);
      call char_to_arraydouble(mupchar, mup);

      key = "Mcluster:alpha_imf"
      call config_getstr(alphaimfchar,
     &      "-1.35, -2.35, -2.7, 0.0, 0.0", key, conf);

      key = "Mcluster:mlim_imf"
      call config_getstr(mlimimfchar,
     &    "0.08, 0.5, 4.0, 100.0, 0.0, 0.0", key, conf);

      key = "Mcluster:alpha_L3"
      call config_getstr(alpha_L3char, "2.3", key, conf);
      call char_to_arraydouble(alpha_L3char, alpha_L3);

      key = "Mcluster:beta_L3"
      call config_getstr(beta_L3char, "1.4", key, conf);
      call char_to_arraydouble(beta_L3char, beta_L3);

      key = "Mcluster:mu_L3"
      call config_getstr(mu_L3char, "0.2", key, conf);
      call char_to_arraydouble(mu_L3char, mu_L3);

      key = "Mcluster:pairing"
      call config_getstr(pairingchar, "3", key, conf);
      call char_to_arrayint(pairingchar, pairing);
      
      key = "Mcluster:adis"
      call config_getstr(adischar, "3", key, conf);
      call char_to_arrayint(adischar, adis);
      do i = 1, numpop
       if (adis(i).eq.6.AND.fracb(i).gt.0.5d0) then
      print*,"ERROR Mcluster:fracb ",i,"th",
     &   " is wrong. It should not be greater",
     &   " then 0.5 if  adis = 6"
!         stop
      endif
      end do

      key = "Mcluster:eigen"
      call config_getstr(eigenchar, "0", key, conf);
      call char_to_arrayint(eigenchar, eigen);

      key = "Mcluster:amin"
      call config_getstr(aminchar, "-2.0", key, conf);
      call char_to_arraydouble(aminchar, amin);

      key = "Mcluster:amax"
      call config_getstr(amaxchar, "10747.0", key, conf);
      call char_to_arraydouble(amaxchar, amax);

      key = "Mcluster:tf"
      call config_getint(tf, 0, key, conf);

      key = "Mcluster:rbar"
      call config_getdouble(rbar, 35.8d0, key, conf);

      key = "Mcluster:rh_mcl"
      call config_getdouble(rh_mcl, 1.0d0, key, conf);
       if (numpop.eq.1.AND.rh_mcl.lt.0.0d0) then
      print*,"ERROR Mcluster:rh_mcl",
     &   " is wrong. It should be greater",
     &   " then 0.0 for only one population"
         stop
      endif

      key = "Mcluster:conc_pop"
      call config_getstr(conc_popchar, "0.5", key, conf);
      call char_to_arraydouble(conc_popchar, conc_pop);

      key = "Mcluster:potential_energy"
      call config_getint(potential_energy, 1, key, conf);

      key = "Mcluster:epoch"
      call config_getstr(epochchar, "0.0", key, conf);
      call char_to_arraydouble(epochchar, epoch_pop);

      key = "Mcluster:zini"
      call config_getstr(zinichar, "0.001", key, conf);
      call char_to_arraydouble(zinichar, zini_pop);
      do i = 1, numpop
          if (zini_pop(i).le.0.0d0) then
            print*,"ERROR Mcluster:zini ",i,"-th",
     &        " is wrong.  It should be > 0.0"
              stop
       endif
      end do
      zini = zini_pop(1);

      key = "Mcluster:seedmc"
      call config_getint(seedmc, 0, key, conf);

      key = "Mcluster:outputf"
      call config_getint(outputf, 0, key, conf);

      key = "Mcluster:check_en"
      call config_getint(check_en, 1, key, conf);

      key = "Mcluster:BSE"
      call config_getint(BSE, 0, key, conf);
      if (BSE.eq.1.AND.(outputf.eq.0.OR.outputf.eq.2)) then
            print*,"ERROR input for MOCCA cannot be evolved",
     &        ". Please set BSE to 0 or change output "
              stop
       endif
      
      key = "Mcluster:pts1"
      call config_getdouble(pts1, 0.0, key, conf);
          print*,"pts1 =", pts1
          if (BSE.eq.1) then
            print*,"pts1 =", pts1
            if (pts1.eq.0.05) then
        print*,"pts1 = 0.05: Hurley et al. (2000)",
     &    " = Level A,B"
          elseif (pts1.eq.0.001) then
        print*,"pts1 = 0.001: Banerjee et al. (2020)", 
     &     " (Startrack vs. BSE) - Level C"
       endif
      endif
 
      key = "Mcluster:pts2"
      call config_getdouble(pts2, 0.0, key, conf);
      print*,"pts2 =",pts2
      if (BSE.eq.1) then
        print*,"pts2 =", pts2
        if (pts2.eq.0.01) then
        print*,"pts2 = 0.01 [Hurley et al. (2000),",
     &         " Banerjee et al. (2020)]",
     &         " (Startrack vs. BSE) - Level A,B,C"
        endif
      endif

      key = "Mcluster:pts3"
      call config_getdouble(pts3, 0.0, key, conf);
      print*,"pts3 =",pts3
      if (BSE.eq.1) then
        print*,"pts3 =", pts3
        if (pts3.eq.0.02) then
        print*,"pts3 = 0.02 [Hurley et al. (2000),",
     &  " Banerjee et al. (2020)]",
     &  " (Startrack vs. BSE) - Level A,B,C"
        endif
      endif

      key = "Mcluster:mdflag"
      call config_getint(mdflag, 0.0, key, conf);
      if (BSE.eq.1) then
      print*,"mdflag =",mdflag
        if (mdflag.eq.1) then
        print*,"mdflag = 1: original [Hurley et al. (2000)] - Level A"
        elseif (mdflag.eq.2) then
        print*,"mdflag = 2: SSE + LBV added [Hurley et al. (2000),",
     & " Humphreys & Davidson (1994)] - Level B"
        elseif (mdflag.eq.2) then
        print*,"mdflag = 3: updates from [Vink et al. (2001, 2005),",
     & " Belczynski et al. (2010)] - Level C"
        elseif (mdflag.eq.2) then
        print*,"mdflag = 4: 3, but without bi-stability jump,",
     & " [Belczynski et al. (2010)] - Level C"
       endif
      endif

      key = "Mcluster:neta"
      call config_getdouble(neta, 0.0, key, conf);
      print*,"neta =", neta
      if (BSE.eq.1) then
            print*,"neta =", neta
      endif

      key = "Mcluster:bwind"
      call config_getdouble(bwind, 0.0, key, conf);
      print*,"bwind =", bwind
      if (BSE.eq.1) then
         if (bwind.eq.0.0) then
         print*,"bwind = 0.0 [Hurley et al. (2000)] - Level A,B,C"
         endif
      endif
 
      key = "Mcluster:hewind"
      call config_getdouble(hewind, 0.0, key, conf);
      print*, "hewind", hewind
      if (BSE.eq.1) then
        if (hewind.eq.1.0) then
           print*,"= 1.0 [Hurley et al. (2000)] - Level A,B,C"
        endif
      endif

      key = "Mcluster:nsflag"
      call config_getint(nsflag, 0, key, conf);
      print*, "nsflag =", nsflag
      if (BSE.eq.1) then
        if (nsflag.eq.0) then
            print*,"nsflag = 0: Original remnant masses",
     & " [Hurley et al. (2000)] - Level A"
        elseif (nsflag.eq.1) then
        print*,"nsflag = 1: FeNi core mass [Belczynski et al. (2002)]",
     & " - Level B"
        elseif (nsflag.eq.2) then
        print*,"nsflag = 2: FeNi core mass [Belczynski et al. (2008)]",
     & " - Level C"
        elseif (nsflag.eq.3) then
        print*,"mdflag = 3: Rapid SNe [Fryer et al. (2012)]",
     & " - Level C"
        elseif (nsflag.eq.4) then
        print*,"mdflag = 4: Delayed SNe [Fryer et al. (2012)]",
     & " - Level C"
       endif
      endif


      key = "Mcluster:psflag"
      call config_getint(psflag, 0, key, conf);
      print*,"psflag =", psflag
      if (BSE.eq.1) then
         if (psflag.eq.0) then
            print*,"psflag = 0: No (P)PISNe - Level A"
         elseif (psflag.eq.1) then
            print*,"psflag = 1: Strong (P)PISNe",
     &   " [Belczynski et al. (2016)]- Level B"
         elseif (psflag.eq.2) then
            print*,"psflag = 2: Moderate (P)PISNe",
     & " [Leung et al. (2019)],",
     & " similar to [Giacobbo et al. (2018)] - Level C"
         elseif (psflag.eq.3) then
            print*,"psflag = 3: Weak (P)PISNe",
     & " [Leung et al. (2019)],",
     & " similar to [Giacobbo et al. (2018)] - Level C"
         endif
      endif
 
      key = "Mcluster:ecflag"
      call config_getint(ecflag, 0, key, conf);
      print*,"ecflag =", ecflag
      if (BSE.eq.1) then
        if (ecflag.eq.0) then
        print*,"ecflag = 0: No ECSNe - Levels A, B"
        elseif (ecflag.eq.1) then
        print*,"ecflag = 1: Enables ECSNe [Belczynski et al. (2008)]",
     & " with disp0 = ecsig - Level C"
        endif
      endif
  
      key = "Mcluster:ifflag"
      call config_getint(ifflag, 0, key, conf);
      print*,"ifflag =", ifflag
      if (BSE.eq.1) then
        if (ifflag.eq.0) then
         print*,"ifflag = 0: use the mass that results",
     & " from the evolution algorithm",
     & " (core-mass growth vs. envelope mass-loss) - Level A,B"
        elseif (ifflag.eq.1) then
        print*,"ifflag = 1: IFMR following",
     & " [Han et al. (1995)] - Level C"
       endif
      endif

      key = "Mcluster:wdflag"
      call config_getint(wdflag, 0, key, conf);
      print*,"wdflag =", wdflag
      if (BSE.eq.1) then
        if (wdflag.eq.0) then
          print*,"wdflag = 0: standard Mestel",
     &    " cooling [Mestel (1952)] - Level A"
        elseif (wdflag.eq.1) then
          print*,"wdflag = 1: modified Mestel",
     & " cooling [Hurley & Shara (2003)]", 
     & " - Levels B,C"
        endif
      endif

      key = "Mcluster:mxns"
      call config_getdouble(mxns, 0.0, key, conf);
      print*,"mxns =", mxns
      if (BSE.eq.1) then
        if (mxns.eq.3.0) then
        print*,"mxns = 3.0: Follows from causality",
     &  " [Hurley et al. (2000),",
     &  " Lattimer & Prakash (2004)] - Level A"
        elseif (mxns.eq.2.5) then
        print*,"mxns = 2.5: From [Fryer et al. (2012)]",
     &  " and GW detections - Level B,C"
        endif
      endif

      key = "Mcluster:bhflag"
      call config_getint(bhflag, 0, key, conf);
      print*,"bhflag =", bhflag
      if (BSE.eq.1) then
        if (bhflag.eq.0) then
        print*,"bhflag = 0: no BH kick. - Level A"
        elseif (bhflag.eq.1) then
        print*,"bhflag = 1: same as for the regular NSs but scaled",
     & " by fallback - Level B"
        endif
      endif  
 
      key = "Mcluster:kmech"
      call config_getint(kmech, 0, key, conf);
      print*,"kmech =", kmech
      if (BSE.eq.1) then
        if (kmech.eq.1) then
          print*,"kmech = 1: standard momentum-conserving",
     &  " [Banerjee et al. (2020)] (Startrack vs. BSE) - Level B"
        elseif (kmech.eq.2) then
          print*,"kmech = 2: convection-asymmetry-driven",
     &  " [Banerjee et al. (2020)] (Startrack vs. BSE) - Level C"
        elseif (kmech.eq.3) then
          print*,"kmech = 3: collapse-asymmetry-driven",
     &   " [Banerjee et al. (2020)] (Startrack vs. BSE) - Level C"
       elseif (kmech.eq.4) then
          print*,"kmech = 4: neutrino-driven [Banerjee",
     &    " et al. (2020)] (Startrack vs. BSE) - Level C"
       endif
      endif 

      key = "Mcluster:disp"
      call config_getdouble(disp, 0.0, key, conf);
      print*,"disp =", disp
      if (BSE.eq.1) then
        if (disp.eq.190.0) then
          print*,"disp = 190 km/s following [Hurley",
     & " et al. (2000)] - Level A"
      elseif (disp.eq.265.0) then
          print*,"disp = 265 km/s following [Hobbs et al. (2005)",
     & " Belczynski et al. (2008)] - Level B,C"
      endif
      endif

      key = "Mcluster:bhspin"
      call config_getint(bhspin, 0, key, conf);
      print*,"bhspin =", bhspin
      if (BSE.eq.1) then
        if (bhspin.eq.0) then
           print*,"bhspin = 0: Fuller model [Belczinsky",
     & " & Banerjee (2020)] - Level C"
        elseif (bhspin.eq.1) then
           print*,"bhspin = 1: Geneva model [Belczinsky",
     & " & Banerjee (2020)] - Level C"
        elseif (bhspin.eq.2) then
          print*,"bhspin = 2: MESA model [Belczinsky",
     & " & Banerjee (2020)] - Level C"
       endif
      endif 

      key = "Mcluster:beta"
      call config_getdouble(beta, 0.0, key, conf);
      print*,"beta =", beta
      if (BSE.eq.1) then
        if (beta.eq.0.125) then
        print*,"beta = 0.125 (does not depend on stellar",
     & " type - Level A,B,C"
        endif
      endif

      key = "Mcluster:xi"
      call config_getdouble(xi, 0.0, key, conf);
      print*,"xi =", xi
      if (BSE.eq.1) then
        if (xi.eq.1.0) then
          print*,"xi = 1 from [Hurley et al. (2002)]"
        endif
      endif

      key = "Mcluster:acc2"
      call config_getdouble(acc2, 0.0, key, conf);
      print*,"acc2 =", acc2
      if (BSE.eq.1) then
        if (acc2.eq.1.5) then
           print*,"acc2 = 1.5 from [Bondi-Hoyle (1944),",
     &     " [Hurley et al. (2002)]"
        endif
      endif
 
      key = "Mcluster:epsnov"
      call config_getdouble(epsnov, 0.0, key, conf);
      print*,"epsnov =", epsnov
      if (BSE.eq.1) then
        if (epsnov.eq.0.001) then
         print*,"epsnov = 0.001 from [Hurley et al. (2002)]"
        endif
      endif
 
      key = "Mcluster:eddfac"
      call config_getdouble(eddfac, 0.0, key, conf);
      print*,"eddfac =", eddfac
      if (BSE.eq.1) then
         if (eddfac.eq.100) then
             print*,"eddfac = 100 Eddington limit",
     &         " accretion rate [Hurley et al. (2002)]"
        endif
      endif
  
      key = "Mcluster:gamma"
      call config_getdouble(gamma, 0.0, key, conf);
      print*,"gamma =", gamma
      if (BSE.eq.1) then
        if (gamma.gt.0) then
          print*,"> 0.0  takes away a fraction gamma of",
     &         " the orbital angular momentum [Hurley et",
     &         " al. (2002)] - Level A,B"
        elseif (gamma.eq.-1.0) then
          print*,"gamma = = -1.0 then assume the lost material",
     &   " carries with it the specific angular momentum of the",
     &   " primary [Hurley et al. (2002)] - Level C"
        elseif (gamma.eq.-2.0) then
          print*,"gamma = -2.0 assume that material is lost from",
     &   " the system as if a wind from the secondary [Hurley",
     &   " et al. (2002)] - Level C"
        endif
      endif
  
      key = "Mcluster:ceflag"
      call config_getint(ceflag, 0, key, conf);
      print*,"ceflag =", ceflag
      if (BSE.eq.1) then
       if (ceflag.gt.0) then
         print*,"ceflag > 0: activates spin-energy correction",
     &   " in common-envelope - Level A"
       elseif (ceflag.eq.1) then
         print*,"ceflag = 3: de Kool (or Podsiadlowski)",
     & " model [de Kool (1990)] - Level C"
       endif
      endif 

      key = "Mcluster:tflag"
      call config_getint(tflag, 0, key, conf);
      print*,"tflag =", tflag
      if (BSE.eq.1) then
        if (tflag.eq.0) then
           print*,"tflag = 0: no tidal circularisation - Level A,B"
        elseif (tflag.eq.1) then
           print*,"tflag = 1: tidal circularisation",
     & " [Hurley et al. (2002)] - Level C"
        endif
      endif   

      key = "Mcluster:alpha1"
      call config_getdouble(alpha1, 0.0, key, conf);
      print*,"alpha1 =", alpha1
      if (BSE.eq.1) then
        if (alpha1.eq.1.0) then
          print*,"alpha1 = 1.0 from [Hurley et al. (2002)]"
        endif
      endif  
 
      key = "Mcluster:lambda"
      call config_getdouble(lambda, 0.0, key, conf);
      print*,"lambda =", lambda
      if (BSE.eq.1) then
        if (lambda.eq.0.5) then
           print*,"lambda = 0.5 from [Hurley et al. (2002)]"
        endif
      endif
*
      return
*
      end
*

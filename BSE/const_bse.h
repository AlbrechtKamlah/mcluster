*
* const_bse.h
*
      INTEGER idum
      COMMON /VALUE3/ idum
      INTEGER idum2,iy,ir(32)
      COMMON /RAND3/ idum2,iy,ir
      INTEGER ktype(0:15,0:15)
      COMMON /TYPES/ ktype
      INTEGER ceflag,tflag,ifflag,nsflag
      INTEGER wdflag,psflag,ecflag
      COMMON /FLAGS/ ceflag,tflag,ifflag,nsflag,wdflag
      COMMON /FLAGS1/ psflag,ecflag
      INTEGER mdflag
      COMMON /FLAGS2/ mdflag
      INTEGER kmech
      COMMON /FLAGS3/kmech 
      INTEGER bhflag
      REAL*8 disp
*
      REAL*8 neta,bwind,hewind,mxns,alpha1,lambda
      REAL*8 beta,xi,acc2,epsnov,eddfac,gamma,FctorCl
      COMMON /VALUE1/ neta,bwind,hewind,mxns,FctorCl
      COMMON /VALUE2/ alpha1,lambda
      COMMON /VALUE4/ disp,BHFLAG
      COMMON /VALUE5/ beta,xi,acc2,epsnov,eddfac,gamma
      REAL*8 pts1,pts2,pts3
      COMMON /POINTS/ pts1,pts2,pts3
      REAL*8 dmmax,drmax
      COMMON /TSTEPC/ dmmax,drmax
      REAL scm(50000,14),spp(20,3)
      COMMON /SINGLE/ scm,spp
      REAL bcm(50000,34),bpp(80,10)
      COMMON /BINARY/ bcm,bpp
      REAL*8 FBFAC,FBTOT,MCO
      INTEGER ECS
      COMMON /FBACK/ FBFAC,FBTOT,MCO,ECS
      INTEGER BHSPIN
	  COMMON /VALUE6/ BHSPIN
*

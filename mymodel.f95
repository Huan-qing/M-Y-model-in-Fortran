

      PROGRAM MAIN
!**********************************************************************
!   Advised by TE ODU
!
!      THIS IS A 1-D VERSION OF THE MELLOR YAMADA, LEVEL 2.5 MODEL
!   AS APPLIED TO OCEAN BOUNDAY LAYERS. THE MOST RECENT DISCRIPTION
!   OF THE MODEL WILL BE FOUND IN REV. GEOPHYS. AND SPACE PHYS.,
!   VOL.20,PP.851-875,1982.(THE MODEL CAN BE CONVERTED TO AN ATMOS-
!   PHERIC B.L. MODEL WITH DUE REGARD FOR THERMODYNAMICS AND THE DIRECT-
!   ION OF THE GRAVITY VECTOR.)
!      ESSENTIALLY 'THE' MODEL RESIDES IN SUBROUTINE PROFQ WHICH PROVID-
!   ES MIXING COEFFICIENTS, KM FOR VELOCITY AND KH FOR TEMPERATURE AND
!   SALINTY. SUBROUTINES PROFU AND PROFV SOLVE FOR THE U AND V VELOCITY
!   WHEREAS PROFT SOLVES FOR TEMPERATURE AND SALINITY (OR ANY OTHER
!   SCALER). ALL OF THESE "PROF" ROUTINES USE THE INTEGRATION METHOD
!   DESCRIBED ON P.198-201 0F RICHTMEYER AND MORTON, INTERSCIENCE PUB.,1967.
!      THE MODEL IS HERE CONFIGURED TO MATCH THE VELOCITY, TEMPERATURE
!   AND SALINITY OF UNDERLYING WATER. HOWEVER, SEE PROFU,PROFV AND
!   PROFQ FOR SIMPLE MODIFICATIONS TO INCLUDE BOTTOM BOUNDARY LAYERS
!   WHICH, WITH SUFFICIENT VERTICAL RESOLUTION, SHOULD PRODUCE A LOG
!   DISTRIBUTION NEAR THE BOTTOM.
!      THE NUMERICS ARE ADAPTED FROM THE 3-D MODEL OF BLUMBERG AND
!   MELLOR WHICH IS IMPLICIT IN THE VERTICAL DIFFUSION TERMS AND LEAP
!   FROG IN THE ADVECTIVE AND HORIZONTAL DIFFUSION TERMS. SINCE THE
!   LATTER ARE EXCLUDED HERE, THE CODE MAY BE CONVERTED TO A TWO LEVEL
!   TIME STEP ALGORITHM FAIRLY EASILY.
!      THERE ARE NUMERICAL PARAMETERS SUCH AS CHOICE OF GRID
!   AND PHYSICAL PARAMETERS  SUCH AS THE GRAVITY CONSTANT  TO BE CHOSEN
!   BY THE USER. THE MODEL PARAMETERS IN PROFQ MAY BE FIDDLED WITH,
!   OF COURSE, BUT WE PREFER TO THINK OF THEM AS FIXED. THE PARAMETER,
!   UMOL, IS "BACKROUND VISCOSITY" PRESUMABLY RELATED TO INTERNAL WAVES
!   OR WHATEVER. EVIDENCE IN THE OCEAN PLACES IT THE RANGE, 1.E-5 TO
!   1.E-4 MKS.                                                          
!                                   G. MELLOR, MARCH 1983
!
!   This code is not supported in the sense of pom.f. It has been 
!   updated since 1983, but I forget when and how.
!                                   G.M. May 1994
!   It would be best to compare with the 3-D version in pom97.f
!   for the latest updates.
!                                   G.M. March 1999
! This version incorporates alterations to the radiation formulation
! and alterations to the dissipation model as described in Mellor, G. L.
! One-dimensional, ocean surface layer modeling: a problem and a solution:
! 2001, J. Phys. Oceanogr.,31, 790-809. 
!                                   G.M ca. 2000
! This version deletes the above dissipatin alteration and inserts the 
! Craig and Banner wave breaking scheme. See PROFQ for details. If the 
! constants cbcnst and surfl are set to zero, the solution reverts
! to the old scheme where Q2 at the surface was related to the surface
! Reynolds stress.
!                                   G.M. March 2003
!***********************************************************************
        INCLUDE 'obl.c'
      DIMENSION UG(KB),VG(KB),WG(KB)
      DIMENSION WINDMON(12),HEATMON(12),PPMON(12) ! separately define the monthly wind speed and heat when adjust the parameter
!     DIMENSION TM(10000),HSTOR(10000),HSTMOD(10000)
!
! *** for plots KT= number of outputs (now: 0.5/day x 180days) 
! *** Note: must set KT=IEND/IPRINT 
!
      PARAMETER (KT=144)
      DIMENSION TZZ(KT,KB),SZZ(KT,KB),AZZ(KT,KB),VZZ(KT,KB)
 6123  FORMAT(144F8.2)   ! format for output (t,z) data
 6124  FORMAT(6F8.2)    ! format for profile (z)  data
!
      DATA KAPPA/0.40/,EPS/1.E-10/,PI/3.14159/
      DATA ISTART/1/,IEND/481/,IPRINT/96/  ! default values (changed below) 
      REAL :: svd
      DATA WINDMON/19.0,18.0,15.0,11.5,9.0,14.2,17.26,21.30,7.7,21.0,20.0,24.0/
      DATA HEATMON/200.0,100.0,0.0,-90.0,-170.0,-167.0,-200.84,-143.45,-1.0,220.0,170.0,190.0/
      DATA PPMON/100,151.65,122.38,106.81,67.14,66.52,60.13,58.80,54.32,53.32,60.38,70.84/

!
!***  set run length (days), time step (sec) and output interval (timesteps) 
!
      days=1440.0              ! 4-year seasonal run
      DT1=300.
      IPRINT=240*60*60/DT1     ! printout every 10 days
      IEND=days*86400./DT1
!     IPRINT=IEND/24
      write(6,'(''DT1,days,IPRINT ='',2F7.2,I5)') DT1,days,IPRINT 
      iprof=2
      write(6,'(''iprof ='',I3)') iprof
!
!*********************************************************************
!        INSERT USER CHOICE OF PROBLEM PARAMETERS.
!        SUBROUTINE DEPTH ESTABLISHES THE VERTICAL GRID. Z VARIES FROM
!   -1 TO 0 SO THAT H*Z IS THE DIMENSIONAL GRID. H*ZZ IS A STAGGERED
!   GRID ON WHICH ARE LOCATED U,V,T AND S. SEE SUBROUTINE DEPTH FOR
!   MORE DETAILS.
! *********************************************************************
      CALL DEPTH(Z,ZZ,DZ,DZZ,KB)
      N=0
!
!c*** set maximum depth
!
      H=100. 
!
      TIME=0.                   
! ***     Coriolis parameter at 32N
      COR=2.*(2.*PI/86400.)*sin(32.*pi/180.)
!
      GRAV=9.806
      SMOTH=0.05
      UMOL=2.0E-5
!****INSERT USER CHOICE OF INITIAL CONDITIONS*************************
      DO 5 K=1,KB
!
!***  initial velocity profile (zero)
!
      UB(K)=0.
!     UB(K)=2.*EXP(-((-ZZ(K)-.5)/.20)**2)
      U(K)=UB(K)
      UF(K)=U(K)
!     VB(K)=-.50*(1.+ZZ(K))
      VB(K)=0.
      V(K)=VB(K)
      VF(K)=V(K)
      UG(K)=U(K)
      VG(K)=V(K)
!****************** special case of S as tracer *******************
      WT(K)=36.7    !unused var, now for holding S in DENS subroutine  
!
      svd=-5.0    !m/day

!     WG(K)=0.0

!      WG(K)=-0.00008    !particle sinking velocity - 0.72mm/s = 20m/day

!      WG(K)=-0.00006   !particle sinking velocity - 0.12mm/s = 10m/day

!      WG(K)=-0.000024  !particle sinking velocity -             2m/day

!      WG(K)=-0.0000012 !particle sinking velocity -           0.1m/day

       WG(K)=svd/24/60/60   !m/sec

!***  initial temperature profile
!
!---  exponential temperature vs. depth
!     TB(K)=10.+15.*EXP(ZZ(K)*H/50.)
!--- thermocline at 5 m
!     TB(K)=29.+4.*ZZ(K)*H/5.
!     if(tb(k).lt.25.)tb(k)=25.
!--- linear temperature stratification (T=17-27C)
!      TB(K)=27.+10.*ZZ(K)*H/20.
!--- special case: NO stratification and constant KM (Ekman test)
      TB(K)=23.+ 3.*ZZ(K)*H/100.  ! 20C at 100m, 23C at 0m
      T(K)=TB(K)
      TF(K)=T(K)
!
!***  initial salinity profile
!

      IF(ZZ(K)*H.GT.-10.)SB(K)=100.-100.*abs(ZZ(K)*H)/10.
      S(K)=SB(K)
      SF(K)=S(K)
!
!***  turbulent parameters
!
      Q2B(K)=1.E-5
      Q2(K)=Q2B(K)
      Q2LB(K)=1.E-5
      Q2L(K)=Q2LB(K)
      Q2LF(K)=Q2L(K)
      L(K)=Q2L(K)/Q2(K)
!     KM(K)=1.e-2
      KM(K)=1.e-5
      KH(K)=KM(K)
      KQ(K)=KM(K)  
   5  CONTINUE
!
!***  calculate density from temperature and salinity
!
      CALL DENS(T,WT,ZZ,H,RHO,KB)
!
      WUBOT=0.
      WVBOT=0.
!***********************************************************************
      DT2=2.*DT1
      DT4=4.*DT1
      DAYI=DT1/86400.
!
!
      WRITE(6,15) ISTART,IEND,DT1
 15   FORMAT(//' MODEL START UP   ISTART =',I4,'   IEND =',I4, &
        '    DT =',F5.0,' S',/)
!
       ICT=0
!
!***********************************************************************
!*                                                                     *
!*      BEGIN TIME MARCH                                               *
!*                                                                     *
!***********************************************************************
!     IEND=5
!     IPRINT=1

      DO 9000 IINT=ISTART,IEND
!     PER=3600.
      TIME=TIME+  DT1/86400.     ! time (days)
      AMON=   TIME/30.         ! time (month) 30d/mon, 360d/yr
      IYEAR=1.+TIME/360.
      AMON=AMON-12*(IYEAR-1)
!
!
!***********************************************************************
!     INSERT USER CHOICE OF GEOSTROHIC VELOCITY AND SURFACE FORCEING
! AS FUNCTIONS OF TIME.
!     QSURF IS HEAT ABSORBED AT SURFACE (POSITIVE IS WARMING) IN W M-2.
! QSWRAD IS SHORT WAVE RADIATION (CAN ONLY BE POSITIVE) IN W M-2 WHICH
! PENETRATES THE SURFACE AND ATTENUATES ACCORDING TO THE ALGORIITHM
! IN PROFT.
! WUSURF=stress/density (8m/s wind ~ 0.1N/m2=1dyn/cm2 -> WUSURF=0.0001 m2/s2)
!
!***  surface wind stress (ocean stress with opposite sign to wind direction)
!
!     WUSURF=0.0           ! NO wind
!     WUSURF=0.0000125     ! ~1 m/s wind (WESTWARD)
!     WUSURF=-0.0000625     ! ~5 m/s wind (EASTWARD)
! windstress=(RHOair/RHOwater)*CD*WIND**2
!     TSURF=-0.0012*0.0012*( 1.**2)      !  1 m/s wind from west
      IMON=NINT(AMON)
      SSURF=-0.0012*0.0012*(WINDMON(IMON)**2)      ! 10 m/s wind from west
      WVSURF=0.                  
!
!*** wind increased gradually to reduce inertial oscillations
!
      if(TIME.LE.1) WUSURF=SSURF*TIME
! --- WindStress: wind decrease from 10m/s in JAN to 5m/s in end of JUN
                     WUSURF=SSURF    ! constant 10m/s
!      IF(AMON.LE.6.)WUSURF=SSURF-0.5*SSURF*AMON/6.             
! ---             wind increase back to 10m/s JUL to DEC
!      IF(AMON.GT.6.AND.AMON.LE.8.)WUSURF=SSURF-0.5*SSURF*(12-AMON)/6. 
!      IF(AMON.GT.8.AND.AMON.LE.9.)WUSURF=-0.0012*0.0012*(21.0**2)  
!      IF(AMON.GT.9.AND.AMON.LE.10.)WUSURF=-0.0012*0.0012*(0.0**2)  
!      IF(AMON.GT.10.AND.AMON.LE.11.)WUSURF=-0.0012*0.0012*(10.**2)
!      IF(AMON.GT.11.)WUSURF=-0.0012*0.0012*(16.**2)                
!
! --- HeatFlux: cooling 100 W/m2 in DEC to heating 100 w/m2 in end of JUN
      TSURF=180./4.19e6 
      WTSURF=HEATMON(IMON)/4.19e6             
!      IF(AMON.LE.6.)WTSURF=225./4.19e6  -2.*TSURF*AMON/6.             
! ---           heating 100 w/m2 in JUL to cooling in DEC
           
!
!     WTSURF=0.                         
      SWRAD =0.                              
      WSSURF=0.
       TSURF=TB(1)
       SSURF=SB(1)
!     N=N+1
!     TM(N)=TIME
!     HSTOR(N)=HSTOR(N-1)-(WTSURF+SWRAD)*DT1
!
!***  special case: Sal as passive tracer injected at upper 10m (max value=100%).
!
       DELTA = PPMON(IMON)  ! assuming 30 days each
       DO 55 K=1,KB
          IF (ZZ(K)*H .GT. -10.0) THEN
             SB(K)  =  DELTA-DELTA*abs(ZZ(K)*H)/10.
             S(K)=SB(K)
             SF(K)=S(K)
          ENDIF
  55  CONTINUE
!    
!***********************************************************************
!
!***  calculate mixing coefficients KM,KH from Mellor-Yamada scheme
!
        CALL PROFQ(DT2)
!
      DO 325 K=1,KB
!********* special case of constant mixing coeff 20cm2/s (Ekman test) *****
!      KM(k)= 20./10000.
!      KH(k)= 20./10000.
      Q2(K)=Q2(K)+.5*SMOTH*(Q2F(K)+Q2B(K)-2.*Q2(K))
      Q2B(K)=Q2(K)
      Q2(K)=Q2F(K)
      Q2L(K)=Q2L(K)+.5*SMOTH*(Q2LF(K)+Q2LB(K)-2.*Q2L(K))
      Q2LB(K)=Q2L(K)
      Q2L(K)=Q2LF(K)
 325  CONTINUE
!
!***  calculate advection terms for temperature and salinity
!
      CALL PROFT(TF,TB,WTSURF,SWRAD,TSURF,1,DT2,TIME)
      CALL PROFT(SF,SB,WSSURF,0.,SSURF,1,DT2,TIME)

      DO 345 K=1,KB
      T(K)=T(K)+.5*SMOTH*(TF(K)+TB(K)-2.*T(K))
      TB(K)=T(K)
      T(K)=TF(K)
      S(K)=S(K)+.5*SMOTH*(SF(K)+SB(K)-2.*S(K))
      SB(K)=S(K)
      S(K)=SF(K)
 345  CONTINUE
!
!*** SPECIAL CASE: add vertical advection due to particle sinking

 

      DO 346 K=2,KB-1

! because z is more nagative when it goes deeper

      S(K)=S(K)-DT1*WG(K)*(SB(K)-SB(K+1))/(DZZ(K)*H)  
      SB(K)=S(K)

!      S(K)=SF(K)

346  CONTINUE

      S(KB)=S(KB-1)

      SB(KB)=SB(KB-1)

      CALL DENS(T,WT,ZZ,H,RHO,KB)   ! constant salinity (S=tracer)
!
      DO 380 K=1,KB-1
      UF(K)=UB(K)+DT2*COR*(V(K)-VG(K))!-DT2*RDRAG*UB(K)
      VF(K)=VB(K)-DT2*COR*(U(K)-UG(K))!-DT2*RDRAG*VB(K)
 380  CONTINUE
      CALL PROFU(DT2)
      CALL PROFV(DT2)
      DO 382 K=1,KB
      U(K)=U(K)+.5*SMOTH*(UF(K)+UB(K)-2.*U(K))
      V(K)=V(K)+.5*SMOTH*(VF(K)+VB(K)-2.*V(K))
      UB(K)=U(K)
      U(K)=UF(K)
      VB(K)=V(K)
      V(K)=VF(K)
 382  CONTINUE    
!**** Begin diagnostic and print section
! PRCON averages calculated profile data before outputing to unit=20
!     CALL PRCON(TIME,DT1)
!
!     HSTMOD(N)=0.
!     DO K=1,KB-1
!     HSTMOD(N)=HSTMOD(N)+T(K)*(Z(K)-Z(K+1))*H
!     ENDDO
!
!--- kb value still as IC
      KM(KB)=KM(KB-1)
      KH(KB)=KH(KB-1)
      Q2(KB)=Q2(KB-1)
       L(KB)= L(KB-1)
      SM(KB)=SM(KB-1)
      SH(KB)=SH(KB-1)
!
      IF(MOD(IINT,IPRINT).NE.0) GO TO 9000
       ICT=ICT+1
      CALL MLDPTH(ZZ,T,KB-1,ZZMLD)
      ZZDMLD=-ZZMLD*H
      SPROD(1)=(ABS(WUSURF)+EPS)**1.5/KAPPA
      BPROD(1)=0.
      SPROD(KB)=(ABS(WUBOT)+EPS)**1.5/KAPPA
      BPROD(KB)=0.
!
!***  print variables at IPRINT intervals
!
      WRITE(6,500) TIME,ZZDMLD,AMON,IYEAR
 500  FORMAT(/,'TIME=',F7.2,' DAYS MLD =',F7.2,' MON=',F6.2,' YR=',I3)
      WRITE(6,501)
 501  FORMAT(1X,'  DEPTH     U         V         T         S      DEPTH &
         Q2        L         KM        KH        RF        SM       SH &
       ')
      DO 550 K=1,KB
      ZZD=ZZ(K)*H
      ZD=Z(K)*H/10.
!     RF=BPROD(K)/(SPROD(K)+1.E-20)
      RF=BPROD(K)/(SPROD(K)+1.E-10)
      WRITE(6,502) ZZD,U(K),V(K),T(K),S(K),ZD,Q2(K),L(K),KM(K), &
         KH(K),RF,SM(K),SH(K)
 502  FORMAT(1X,F7.2,4(1PE10.2),F7.2,7(1PE10.2))
!
! ----- save for matlab plot of variable(depth,time)
!       (KM in cm2/s, speed VEL in cm/s)
!
       TZZ(ICT,K)=T(K)
       SZZ(ICT,K)=S(K)
       AZZ(ICT,K)=KM(K)*10000.
       VZZ(ICT,K)=SQRT(U(K)**2+V(K)**2)*100.
!
 550  CONTINUE
      WRITE(6,504) 0.,WUSURF,WVSURF,WTSURF,WSSURF
 504  FORMAT(/,1X,F6.0,4(1PE10.2),'     =   SURFACE FLUXES')
      WRITE(6,505) H,WUBOT ,WVBOT,0.,0.
 505  FORMAT(1X,F6.0,4(1PE10.2),'     =   BOTTOM  FLUXES')

 9000 CONTINUE
 9001 WRITE(6,'('' END OF RECORD, TIME ='',F10.2)') TIME

!     WRITE(6,'(''   TIME       HSTOR     HSTORMOD '')')
!     DO I=1,KT
!     HSTMOD(I)=HSTMOD(I)-HSTMOD(1)
!     ENDDO
!     WRITE(6,'(3F10.3)') (TM(I),HSTOR(I),HSTMOD (I),I=1,N,KT)
!
! ----- save for matlab plot of depth vs, time
!
      open(61,file='TZZ_sv50.dat',form='formatted')  ! T(t,z)
      open(62,file='SZZ_sv50.dat',form='formatted')  ! S(t,z)
      open(63,file='KZZ_sv50.dat',form='formatted')  ! Km(t,z)
      open(64,file='VZZ_sv50.dat',form='formatted')  ! Vel spd(t,z)
      open(65,file='UVTS_sv50.dat',form='formatted')  ! profiles at end of run
      do k=1,kb
       write(61,6123)(TZZ(i,k),i=1,KT)
       write(62,6123)(SZZ(i,k),i=1,KT)
       write(63,6123)(AZZ(i,k),i=1,KT)
       write(64,6123)(VZZ(i,k),i=1,KT)
       write(65,6124)ZZ(k)*H,U(k)*100,V(k)*100,T(k),S(k),KM(k)*1e4
!      write(61,'(48f7.3)')(TZZ(i,k),i=1,48)
!      write(62,'(48f7.3)')(SZZ(i,k),i=1,48)
!      write(63,'(48f7.2)')(AZZ(i,k),i=1,48)
      enddo

      STOP
      END
      
      SUBROUTINE PRCON(TIME,DT1)
!  Averages properties over NDAY steps
	INCLUDE 'obl.c'         
	REAL KMAVE
      DIMENSION VAVE(KB),TAVE(KB),QAVE(KB),KMAVE(KB),ZD(KB)
      DATA VAVE/KB*0./,TAVE/KB*0./,QAVE/KB*0./,KMAVE/KB*0./
      DATA N/0/
! NDAY=# of steps in NHOURS  
      NHOURS=6
      NDAY=  NHOURS*3600./DT1
      DO K=1,KB
      VAVE(K)=VAVE(K)+SQRT(V(K)**2+U(K)**2)  
      TAVE(K)=TAVE(K)+T(K)
      QAVE(K)=QAVE(K)+SQRT(Q2(K))
      KMAVE(K)=KMAVE(K)+KM(K)
      ENDDO
      N=N+1
      IF(MOD(IINT,NDAY).EQ.0) THEN     
	CON=1/FLOAT(N)
        DO K=1,KB
	VAVE(K)=VAVE(K)*CON
	TAVE(K)=TAVE(K)*CON
	QAVE(K)=QAVE(K)*CON
	KMAVE(K)=KMAVE(K)*CON
        ZD(K)=H*ZZ(K)
        ENDDO
!       WRITE(6,'('' N = '',I5)') N
        TIM=TIME-0.5
        WRITE(20,'(F12.4)') TIM 
        WRITE(20,'(1x,I5,5F12.4)') &
          (K,ZD(K),VAVE(K),TAVE(K),QAVE(K),KMAVE(K),K=1,KB)                
	N=0
        DO K=1,KB
        VAVE(K)=0.
        TAVE(K)=0.
        QAVE(K)=0.
        KMAVE(K)=0.
        ENDDO
      ENDIF              
      RETURN
      END
!
!*DECK DENS
      SUBROUTINE DENS(t,s,zz,dt,rho,kb)
      implicit double precision (a-h,o-z)
      real t,s,zz,dt,grav,rho
      DIMENSION t(KB),s(KB),rho(KB),zz(KB)
      DATA GRAV/9.806/
!
!   ......... THIS SUBROUTINE COMPUTES DENSITY-1.............
!              T = POTENTIAL TEMPERATURE
!    Mellor,G.L., 1991: An equation of state for numerical models of 
!    Oceans and Estuaries. J.Atmos.and Oceanic Tech. 8, 609-611
!
      DO 1 K=1,KB-1
      TR=T(K)
      SR=S(K)
!         Here, the (approximate) pressure is in units of decibars.
      P=-GRAV*1.025*ZZ(K)*DT*0.1
!
      RHOR = 999.842594 + 6.793952E-2*TR &
           - 9.095290E-3*TR**2 + 1.001685E-4*TR**3 & 
           - 1.120083E-6*TR**4 + 6.536332E-9*TR**5 
!
      RHOR = RHOR + (0.824493 - 4.0899E-3*TR &
           + 7.6438E-5*TR**2 - 8.2467E-7*TR**3 &
           + 5.3875E-9*TR**4) * SR &
           + (-5.72466E-3 + 1.0227E-4*TR &
           - 1.6546E-6*TR**2) *(ABS(SR))**1.5 &
           + 4.8314E-4 * SR**2
!  For shallow water the pressure dependency can be neglected
!  in which case it should also be omitted in PROFQ
!     CR =1449.1+.00821*P+4.55*TR-.045*TR**2+1.34*(SR-35.)
!     RHOR=RHOR + 1.E4*P/CR**2*(1.-0.2*P/CR**2)
!
      RHO(K)=(RHOR-1000.)*1.E-3
    1 CONTINUE

!
      return
      end
!*DECK DEPTH
      SUBROUTINE DEPTH(Z,ZZ,DZ,DZZ,KB)
      DIMENSION Z(KB),ZZ(KB),DZ(KB),DZZ(KB)
      KL1=6        
      KL2=KB-2
!***********************************************************************
!   THIS SUBROUTINE ESTABLISHES THE VERTICAL RESOLUTION WITH LOG
!   DISTRIBUTIONS  AT THE TOP AND BOTTOM AND A LINEAR DISTRIBUTION
!   BETWEEN. CHOOSE KL1 AND KL2. THE DEFAULT KL1 = .3*KB AND KL2 = KB-2
!   YIELDS A LOG DISTRIBUTION AT THE TOP AND NONE AT THE BOTTOM.
!***********************************************************************
!     KL1=.3*KB
      BB=FLOAT(KL2-KL1)+4.
      CC=FLOAT(KL1)-2.
      DEL1=2./BB/EXP(.693147*FLOAT(KL1-2))
      DEL2=2./BB/EXP(.693147*FLOAT(KB-KL2-1))
      Z(1)=0.
      ZZ(1)=-DEL1/2.
      DO 3 K=2,KL1
      Z(K)=-DEL1*EXP(.693147*FLOAT(K-2))
      ZZ(K)=-DEL1*EXP(.693147*(FLOAT(K)-1.5))
    3 CONTINUE
      DO 4 K=KL1,KL2
      Z(K)=-(FLOAT(K)-CC)/BB
      ZZ(K)=-(FLOAT(K)-CC+0.5)/BB
    4 CONTINUE
      DO 5 K=KL2,KB-1
      Z(K)=(1.0-DEL2*EXP(.693147*FLOAT(KB-K-1)))*(-1.)
      ZZ(K)=(1.0-DEL2*EXP(.693147*(FLOAT(KB-K)-1.5)))*(-1.)
    5 CONTINUE
      Z(KB)=-1.0
      ZZ(KB-1)=-1.*(1.-DEL2/2.)
      ZZ(KB)=-1.*(1.+DEL2/2.)
!
!     DO 10 K=1,KB
!  10 Z(K)=-FLOAT(K-1)/FLOAT(KB-1)
!     DO 11 K=1,KB-1
!  11 ZZ(K)=.5*(Z(K)+Z(K+1))
!
      DO 6 K=1,KB-1
      DZ(K)=Z(K)-Z(K+1)
      DZZ(K)=ZZ(K)-ZZ(K+1)
    6 CONTINUE
      WRITE(6,'(''    K       Z       ZZ        DZ       DZZ'')') 
      WRITE(6,'(I5,4F10.5 )') (K,Z(K),ZZ(K),DZ(K),DZZ(K),K=1,KB)
      RETURN
      END
!* DECK MLDPTH
      SUBROUTINE MLDPTH(ZZ,T,KB,ZZMLD)
      DIMENSION ZZ(KB),T(KB)
      DO 100 K=1,KB-1
      IF(T(1).GT.(T(K)+.2)) GO TO 100
      ZZMLD=ZZ(K)-(T(K)+.2-T(1))*(ZZ(K)-ZZ(K+1))/(T(K)-T(K+1))
  100 CONTINUE
      RETURN
      END
!
      SUBROUTINE PROFQ(DT2)
!       Differential length scale version
	INCLUDE 'obl.c'          
      DATA A1,B1,A2,B2,C1/0.92,16.6,0.74,10.1,0.08/
      DATA E1/1.8/,E2/1.33/,E3/1.0/
      DATA SM/KB*0.39/,SH/KB*0.49/,GM/KB*0.154/,GH/KB*0.154/
      DATA FFI/1.0/,GRAV/9.807/,KAPPA/0.4/  
      DATA SMALL/1.E-10/
!
!      cbcnst=0.             ! Craig-Banner constant 
!     surfl=0.              ! Stacey constant
      shiw=0.0              ! Internal wave constant
      cbcnst=100.           ! Craig-Banner constant 
      surfl=2.E5            ! Stacey constant

!
      DH=H
      DO 100 K=2,KBM1
      A(K)=-DT2*(KQ(K+1)+KQ(K)+UMOL )*.5/(DZZ(K-1)*DZ(K) &
          *DH*DH)
      C(K)=-DT2*(KQ(K-1)+KQ(K)+UMOL )*.5/(DZZ(K-1)*DZ(K-1) &
          *DH*DH)
 100  CONTINUE
!***********************************************************************
!                                                                      *
!        THE FOLLOWING SECTION SOLVES THE EQUATION                     *
!        DT2*(KQ*Q2')' - Q2*(2.*DT2*DTEF+1.) = -Q2B                    *
!                                                                      *
!***********************************************************************
      CONST1=(16.6*FFI)**.6666667
      UTAU=(WUSURF**2+WVSURF**2)**.25
! Calculate surface length (Stacey, JPO, 1999, 1363-1365; Terray etal,
! 1999, Proc. Symp. on "The Wind-driven Air-Sea Interface, Sydnery
! Australia)  
      L0=surfl*UTAU**2/GRAV
      EE(1)=1.0
! Craig & Banner (JPO, 1994, 2546-2559) surface Q2 diffusion scheme.
!  DF0 is located at Z(1); DF1 is interpolated diffusion at ZZ(1)
      DF0=2.*cbcnst*utau**3
      DF2=.5*(KQ(2)+KQ(3))*(Q2B(2)-Q2B(3))/(DZ(2)*DH)
      DF1=.66667*DF0+.33333*DF2
      GG(1)=DF1*DZ(1)*H*2./ABS(KQ(1)+KQ(2))                

      Q2F(KB)=.5*SQRT((WUBOT+WUBOT)**2 &
               + (WVBOT+WVBOT)**2)*CONST1
      STF(1)=1.
      DO 101 K=2,KBM1
      BOYGR(K)=GRAV*(RHO(K-1)-RHO(K))/(DZZ(K-1)*DH)
      GH(K)=L(K)**2/ABS(Q2B(K))*BOYGR(K)                                 
      IF(GH(K).GT.0.028) GH(K)=0.028
      STF(K)=1.
      DTEF(K)=SQRT(ABS(Q2B(K)))/(B1*L(K)+1.E-30)*STF(K)
      SPROD(K)=.25*KM(K)*((U(K)+U(K)-U(K-1)- &
         U(K-1))**2 +(V(K)+V(K)-V(K-1)- &
         V(K-1))**2)/(DZZ(K-1)*DH)**2*FFI**2 &
         -shiw*KM(K)*BOYGR(K)
      BPROD(K)=KH(K)*BOYGR(K)
      PROD(K)=SPROD(K)+BPROD(K)
  101 CONTINUE
      DO 102 K=2,KBM1
      GG(K)=1./(A(K)+C(K)*(1.-EE(K-1)) &
          -(2.*DT2*DTEF(K)+1.) )
      EE(K)=A(K)*GG(K)
      GG(K)=(-2.*DT2*PROD(K) &
        +C(K)*GG(K-1)-Q2B(K))*GG(K)
 102  CONTINUE
      DO 103 K=1,KBM1
      KI=KB-K
      Q2F(KI)=EE(KI)*Q2F(KI+1)+GG(KI)
! Code works quite OK with following deleted. Harmless negative
! Q2 are generated of order E-9 to E-11.
      Q2F(KI)=ABS(Q2F(KI))
 103  CONTINUE
!***********************************************************************
!                                                                      *
!        THE FOLLOWING SECTION SOLVES THE EQUATION                     *
!        DT2(KQ*Q2L')' - Q2L*(DT2*DTEF+1.) = -Q2LB                     *
!                                                                      *
!***********************************************************************
      EE(1)=0.
      GG(1)=0.                  
      Q2LF(KB)=0.
      DO 110 K=2,KBM1
      DTEF(K) =DTEF(K)*(1.+E2*((1./ABS(Z(K)-Z(1))+ &
          1./ABS(Z(K)-Z(KB))) *L(K)/(DH*KAPPA))**2)
      GG(K)=1./(A(K)+C(K)*(1.-EE(K-1)) &
          -(DT2*DTEF(K)+1.))
      EE(K)=A(K)*GG(K)
      GG(K)=(DT2*(-(SPROD(K)+E3*BPROD(K)) &
         *L(K)*E1)+C(K)*GG(K-1)-Q2LB(K))*GG(K)
 110  CONTINUE
      DO 111 K=1,KBM1
      KI=KB-K
      Q2LF(KI)=EE(KI)*Q2LF(KI+1)+GG(KI)
 111  CONTINUE
!***********************************************************************
!                                                                      *
!      THE FOLLOWING SECTION SOLVES FOR KM AND KH                      *
!                                                                      *
!***********************************************************************
      L(1)=0.   
      L(KB)=0.
      DO 212 K=1,KBM1
      L(K)=ABS(Q2LF(K)/Q2F(K))
      IF(Z(K).GT.-0.5) &
         L(K)=MAX(L(K),KAPPA*L0)
 212  CONTINUE
!**************************************************************************
! NOTE THAT SM,SH LIMIT TO INFINITY WHEN GH APPROAHES 0.0288
      SM(1)=.39
      SH(1)=.49
      DO 202 K=2,KBM1
      SH(K)=A2*(1.-6.*A1/B1)/(1.-(3.*A2*B2+18.*A1*A2)*GH(K))
      SM(K)=A1*(1.-3.*C1-6.*A1/B1)+SH(K)*(18.*A1*A1+9.*A1*A2)*GH(K)
      SM(K)=SM(K)/(1.-9.*A1*A2*GH(K))
      SH(K)=AMAX1(SH(K),0.)
      SM(K)=AMAX1(SM(K),0.)
  202 CONTINUE
      SM(1)=SM(2)
      SH(1)=SH(2)
!     IF(MOD(IINT,IPRINT).EQ.0) THEN
!     WRITE(6,'(''IINT ='',I7)') IINT
!     DO K=1,5 
!     WRITE (6,'(1X, I5,10E12.3)') K,Q2F(K),L(K),Q2LF(K),GH(K),SM(K),
!    1   SH(K),KN(K),KQ(K),KM(K),KH(K)
!     ENDDO          
!     ENDIF
      DO 204 K=1,KBM1
      KN(K)=L(K)*SQRT(ABS(Q2(K)))
!     KQ(K)=(KN(K)*SQ+KQ(K))*.5
      KQ(K)=(KN(K)*.41*SH(K)+KQ(K))*.5
      KM(K)=(KN(K)*SM(K)+KM(K))*.5
      KH(K)=(KN(K)*SH(K)+KH(K))*.5
  204 CONTINUE
      RETURN
      END


!--------------------------------------------------------------------
      SUBROUTINE PROFT(FF,FB,WFSURF,SWRAD,FSURF,NBC,DT2,TIME)
!
	INCLUDE 'obl.c'              
      DIMENSION FB(KB),FF(KB)
      DIMENSION R(5),AD1(5),AD2(5),RAD(KB)
      DIMENSION ADF(14)
! Irradiance parameters after Paulson and Simpson, JPO, 1977, 952-956.
      NTP=2
!       NTP         =     1      2       3       4       5
!   JERLOV TYPE     =     I      IA      IB      II      III
      DATA R   /        .58 ,   .62  ,  .67  ,  .77  ,  .78   /
      DATA AD1 /        .35 ,   .60  ,  1.0  ,  1.5  ,  1.4   /
      DATA AD2 /        23  ,   20   ,  17   ,  14   ,  7.9   /
!  Special for Flex data.
      DATA ADF / 23.0,20.7,18.2,15.8,13.4,11.0,14.5,19.5, &
             20.9,17.0,13.1,16.5,19.8,22.5/
       NX=1+(TIME-95)/5
       AD2(2)=ADF(NX)+(ADF(NX+1)-ADF(NX))*(TIME-95-NX*5)/5

      UMOLPR=UMOL
!***********************************************************************
!                                                                      *
!        THE FOLLOWING SECTION SOLVES THE EQUATION                     *
!         DTI2*(KH*F')'-F=-FB                                          *
!                                                                      *
!***********************************************************************
      DH=H
      DO 20 K=2,KB-1
      A(K-1)=-DT2*(KH(K)+UMOLPR)/(DZ(K-1)*DZZ(K-1)*DH &
           *DH)
      C(K)=-DT2*(KH(K)+UMOLPR)/(DZ(K)*DZZ(K-1)*DH &
           *DH)
   20 CONTINUE
!-----------------------------------------------------------------------
!   NBC=1: SURF. B.C. IS WFSURF. NO SW RADIATIVE PENETRATION.
!   NBC=2; SURF. B.C. IS WFSURF+SWRAD*(1.-TR). TR*SWRAD PENETRATES WATER COLUMN
!   NBC=3; SURF. B.C. IS TSURF. NO SW RADIATIVE PENETRATION.
!   NBC=4; SURF. B.C. IS TSURF. TR*SWRAD PENETRATES WATER COLUMN
!
! NOTE THAT WTSURF (=WFSURF) AND SWRAD ARE NEG. VALUES WHEN WATER COLUMN IS
!            WARMING.
!-----------------------------------------------------------------------
!------------------------------------------------------------------
!     Penetrative Radiation Calculation. At the bottom any 
!     unattenuated radiation is deposited in the bottom layer.
!------------------------------------------------------------------
        DO 512 K=1,KB
        RAD(K)=0.
  512   CONTINUE
      IF(NBC.EQ.2.OR.NBC.EQ.4.) THEN   
      DO K=1,KB-1
      RAD(K)=SWRAD &
         *(   R(NTP)*EXP(Z(K)*DH/AD1(NTP)) &
                   +(1.-R(NTP))*EXP(Z(K)*DH/AD2(NTP))   )
      ENDDO
      ENDIF
      GO TO (50,51,52,52), NBC
   50 CONTINUE
      EE(1)=A(1)/(A(1)-1.0)
      GG(1)=-DT2*WFSURF/(-DZ(1)*DH)-FB(1)
  500 GG(1)=GG(1)/(A(1)-1.0)
      GO TO 53
!
   51 CONTINUE
      EE(1)=A(1)/(A(1)-1.0)
      GG(1)=DT2*(WFSURF+RAD(1)-RAD(2)) &
           /(DZ(1)*DH)-FB(1)
  510 GG(1)=GG(1)/(A(1)-1.0)
      GO TO 53
!
   52 CONTINUE
      EE(1)=0.
  520 GG(1)=FSURF
!----------------------------------------------------------------------
   53 CONTINUE
!
!----------------------------------------------------------------------
      DO 101 K=2,KB-2
      GG(K)=1./(A(K)+C(K)*(1.-EE(K-1))-1.)
      EE(K)=A(K)*GG(K)
      GG(K)=(C(K)*GG(K-1)-FB(K) &
           +DT2*(RAD(K)-RAD(K+1))/(DH*DZ(K)))*GG(K)
  101 CONTINUE
!-----  BOTTOM ADIABATIC B.C. ------------------------------------------
  102 FF(KB-1)=((C(KB-1)*GG(KB-2)-FB(KB-1) &
              +DT2*(RAD(KB-1)-RAD(KB))/(DH*DZ(KB-1))) &
                /(C(KB-1)*(1.-EE(KB-2))-1.))
!----------------------------------------------------------------------
      DO 105 K=2,KB-1
      KI=KB-K
      FF(KI)=(EE(KI)*FF(KI+1)+GG(KI))
  105 CONTINUE
!
      RETURN
      END
!
      SUBROUTINE PROFU(DT2)
        INCLUDE 'obl.c'            
!***********************************************************************
!                                                                      *
!        THE FOLLOWING SECTION SOLVES THE EQUATION                     *
!         DT2*(KM*U')' - U= -UB                                        *
!                                                                      *
!***********************************************************************
 85   DH=H
      DO 100 K=2,KB-1
      A(K-1)=-DT2*(KM(K)+UMOL  )/(DZ(K-1)*DZZ(K-1)*DH &
           *DH)
      C(K)=-DT2*(KM(K)+UMOL  )/(DZ(K)*DZZ(K-1)*DH &
           *DH)
 100  CONTINUE
      EE(1)=A(1)/(A(1)-1.)
      GG(1)=(-DT2*WUSURF/(-DZ(1)*DH)-UF(1)) &
         /(A(1)-1.)
      DO 101 K=2,KB-2
      GG(K)=1./(A(K)+C(K)*(1.-EE(K-1))-1.)
      EE(K)=A(K)*GG(K)
      GG(K)=(C(K)*GG(K-1)-UF(K))*GG(K)
 101  CONTINUE
      CBC=AMAX1(.0025,.16/ALOG((ZZ(KB-1)-Z(KB))*DH/.01)**2)
      CBC=CBC*SQRT(UB(KB-1)**2+(.25*(VB(KB-1) &
           +VB(KB-1)+VB(KB-1)+VB(KB-1)))**2)
!********************************************************************
!        TO RESTORE BOTTOM B.L. DELETE NEXT LINE
!*********************************************************************
      CBC=0.
      UF(KB-1)=(C(KB-1)*GG(KB-2)-UF(KB-1))/(CBC &
           *DT2/(-DZ(KB-1)*DH)-1.-(EE(KB-2)-1.)*C(KB-1))
      DO 103 K=2,KB-1
      KI=KB-K
      UF(KI)=EE(KI)*UF(KI+1)+GG(KI)
 103  CONTINUE
 92   WUBOT=-CBC*UF(KB-1)
      RETURN
      END
!*DECK PROFV
      SUBROUTINE PROFV(DT2)
        INCLUDE 'obl.c'           
!***********************************************************************
!                                                                      *
!        THE FOLLOWING SECTION SOLVES THE EQUATION                     *
!         DT2*(KM*V')' -V= -VB                                         *
!                                                                      *
!***********************************************************************
      DH=H
      DO 100 K=2,KB-1
      A(K-1)=-DT2*(KM(K)+UMOL  )/(DZ(K-1)*DZZ(K-1)*DH &
           *DH)
      C(K)=-DT2*(KM(K)+UMOL  )/(DZ(K)*DZZ(K-1)*DH &
           *DH)
 100  CONTINUE
      EE(1)=A(1)/(A(1)-1.)
      GG(1)=(-DT2*WVSURF/(-DZ(1)*DH)-VF(1)) &

        /(A(1)-1.)
 98   CONTINUE
      DO 101 K=2,KB-2
      GG(K)=1./(A(K)+C(K)*(1.-EE(K-1))-1.)
      EE(K)=A(K)*GG(K)
      GG(K)=(C(K)*GG(K-1)-VF(K))*GG(K)
 101  CONTINUE
 104  CBC=AMAX1(.0025,.16/ALOG((ZZ(KB-1)-Z(KB))*DH/.01)**2)
      CBC=CBC*SQRT((.25*(UB(KB-1)+UB(KB-1) &
         +UB(KB-1)+UB(KB-1)))**2+VB(KB-1)**2)
!********************************************************************
!        TO RESTORE BOTTOM B.L. DELETE NEXT LINE
!*********************************************************************
      CBC=0.
      VF(KB-1)=(C(KB-1)*GG(KB-2)-VF(KB-1))/(CBC &
           *DT2/(-DZ(KB-1)*DH)-1.-(EE(KB-2)-1.)*C(KB-1))
      DO 103 K=2,KB-1
      KI=KB-K
      VF(KI)=EE(KI)*VF(KI+1)+GG(KI)
 103  CONTINUE
 92   WVBOT=-CBC*VF(KB-1)
      RETURN
      END
!********************************************************************
!                   END OF FILE                       
!*********************************************************************
	    
	    	    
	    

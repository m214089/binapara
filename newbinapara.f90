PROGRAM newbinapara_test_driver

  !    A Fortran 90 test driver for the parameterization
  !
  !    Copyright (C)2018 Määttänen et al. 2018
  !    
  !    anni.maattanen@latmos.ipsl.fr
  !    joonas.merikanto@fmi.fi
  !    hanna.vehkamaki@helsinki.fi
  !
  !    References
  !    A. Määttänen, J. Merikanto, H. Henschel, J. Duplissy, R. Makkonen, 
  !    I. K. Ortega and H. Vehkamäki (2018), New parameterizations for 
  !    neutral and ion-induced sulfuric acid-water particle formation in 
  !    nucleation and kinetic regimes, J. Geophys. Res. Atmos., 122, doi:10.1002/2017JD027429.
  !
  !    Lehtinen, K. E., Dal Maso, M., Kulmala, M., and Kerminen, V. M., 
  !    Journal of Aerosol Science, 38(9), 988-994, 2007

  DOUBLE PRECISION :: t        ! temperature in K 
  DOUBLE PRECISION :: satrat   ! saturatio ratio of water (between zero and 1)
  DOUBLE PRECISION :: rhoa     ! sulfuric acid concentration in 1/cm3
  DOUBLE PRECISION :: csi      ! Inverse of lifetime of ions (s-1)
  DOUBLE PRECISION :: cs       ! H2SO4 condensation sink (h-1)
  DOUBLE PRECISION :: airn     ! Air molecule concentration in (cm-3)
  DOUBLE PRECISION :: ipr      ! Ion pair production rate (cm-3 s-1)
  DOUBLE PRECISION :: jnuc_n   ! Neutral nucleation rate in 1/cm3s (J>10^-7 1/cm3s)
  DOUBLE PRECISION :: ntot_n   ! total number of molecules in the neutral critical cluster
  DOUBLE PRECISION :: jnuc_i   ! Charged nucleation rate in 1/cm3s (J>10^-7 1/cm3s)
  DOUBLE PRECISION :: ntot_i   ! total number of molecules in the charged critical cluster
  DOUBLE PRECISION :: x_n      ! mole fraction of H2SO4 in the neutral critical cluster 
  DOUBLE PRECISION :: x_i      ! mole fraction of H2SO4 in the charged critical cluster 
                               ! (note that x_n=x_i in nucleation regime) 
  DOUBLE PRECISION :: na_n     ! sulfuric acid molecules in the neutral critical cluster
  DOUBLE PRECISION :: na_i     ! sulfuric molecules in the charged critical cluster
  DOUBLE PRECISION :: rc_n     ! radius of the charged critical cluster in meters 
  DOUBLE PRECISION :: rc_i     ! radius of the charged critical cluster in meters 
  DOUBLE PRECISION :: n_i      ! number of ion pairs in air (cm-3) 
  LOGICAL  :: kinetic_n        ! true if kinetic neutral nucleation
  LOGICAL  :: kinetic_i        ! true if kinetic ion-induced nucleation
  DOUBLE PRECISION :: rhoatres ! treshold concentration of h2so4 (1/cm^3) for neutral kinetic nucleation
  DOUBLE PRECISION :: GR       ! Particle growth rate (nm/h) 
  DOUBLE PRECISION :: gamma,m,dx,d1,LKK_n,LKK_i ! Factors for Lehtinen et al. param
!unused  DOUBLE PRECISION :: jnuc_n_17,jnuc_i_17       ! Nucleation rates at 1.7nm mobility diameter

  t=249.0
  satrat=0.4
  rhoa=1.e7  
  csi=1.0/480.  ! Inverse lifetime of ions 
  airn=DBLE(6.023E23)*DBLE(1.E-6)*DBLE(1.013E5)/DBLE(8.31)/t
  ipr=20.

  CALL newbinapara(t,satrat,rhoa,csi,airn,ipr,jnuc_n,ntot_n,jnuc_i,ntot_i,&
                   &   x_n,x_i,na_n,na_i,rc_n,rc_i,n_i,kinetic_n,kinetic_i,rhoatres)


  ! Nucleation rates at 1.7nm using Lehtinen et al.

  m=-1.6  ! Average m-value according to Lehtinen et al. 
  GR=rhoa/(661.1*(satrat*100)**2-1.129E5*(satrat*100)+1.549E7) ! Nieminen et al., 2010
  cs=10  ! Typical CS value in atmosphere in 1/h 
  dx=1.4 ! Target size (in geometric diameter = mobility diameter -0.3nm)
  
  d1=2.*rc_n*1.E9 ! diameter of critical cluster in neutral case
  gamma=MAX(DBLE(0.0),DBLE(1.0)/(m+1)*( (dx/(d1))**(m+1)-1)) ! gamma-factor in Lehtinen et al. (neutral)
  LKK_n=EXP(-gamma*d1*cs/gr)   ! Final scaling factor in Lehtinen et al. for neutral case 

  d1=2.*rc_i*1.E9 ! diameter of critical cluster in charged case 
  gamma=MAX(DBLE(0.0),DBLE(1.0)/(m+1)*( (dx/(d1))**(m+1)-1)) 
  LKK_i=EXP(-gamma*d1*cs/gr)   ! For charged case

  WRITE(*,*) 'Input:'
  WRITE(*,*) ' T(K)    RH(%)      SA(cm-3)     csi (s-1)    airn(cm-3)   ipr(s-1cm-3)'
  WRITE(*,'(F7.1,F11.5,3ES13.4,F8.2)') t,satrat*100,rhoa,csi,airn, ipr
  WRITE(*,*)
  WRITE(*,*) 'Output, neutral nucleation:'
  WRITE(*,*) ' J_n (cm-3s-1)  Ntot_n     x         Na_n   rc_n(nm)       Kinetic  rhoatres(cm-3)'
  WRITE(*,'(ES12.4,3F10.2,ES13.4,L6,ES18.4)') jnuc_n,ntot_n,x_n,na_n,rc_n*1.E9,kinetic_n,rhoatres
  WRITE(*,*)
  WRITE(*,*) 'Output, ion-induced nucleation:'
  WRITE(*,*) ' J_i (cm-3s-1)  Ntot_i     x         Na_i   rc_i(nm)       Kinetic  Ions(cm-3) '
  WRITE(*,'(ES12.4,3F10.2,ES13.4,L6,ES18.4)') jnuc_i,ntot_i,x_i,na_i,rc_i*1.E9,kinetic_i, n_i
  WRITE(*,*)
  WRITE(*,*) 'Output, nucleation rates for 1.7nm mob. diam. particles using Lehtinen et al.(2007)'
  WRITE(*,*) ' J_n17       J_i17'
  WRITE(*,'(4ES12.4)') jnuc_n*LKK_n, jnuc_i*LKK_i


END PROGRAM newbinapara_test_driver


SUBROUTINE newbinapara(t,satrat,rhoa,csi,airn,ipr,jnuc_n,ntot_n,jnuc_i,ntot_i,&
                   &   x_n,x_i,na_n,na_i,rc_n,rc_i,n_i,kinetic_n,kinetic_i,rhoatres)

  !    Fortran 90 subroutine newbinapara
  !
  !    Calculates parametrized values for neutral and ion-induced sulfuric acid-water particle formation rate
  !    of critical clusters, 
  !    number of particles in the critical clusters, the radii of the critical clusters
  !    in H2O-H2SO4-ion system if temperature, saturation ratio of water, sulfuric acid concentration, 
  !    and, optionally, either condensation sink due to pre-existing particles and ion pair production rate,
  !    or atmospheric concentration of negative ions are given. 
  !
  !    The code calculates also the kinetic limit and the particle formation rate
  !    above this limit (in which case we set ntot=1 and na=1)
  !
  !    Copyright (C)2018 Määttänen et al. 2018
  !    
  !    anni.maattanen@latmos.ipsl.fr
  !    joonas.merikanto@fmi.fi
  !    hanna.vehkamaki@helsinki.fi
  !
  !    References
  !    A. Määttänen, J. Merikanto, H. Henschel, J. Duplissy, R. Makkonen, 
  !    I. K. Ortega and H. Vehkamäki (2018), New parameterizations for 
  !    neutral and ion-induced sulfuric acid-water particle formation in 
  !    nucleation and kinetic regimes, J. Geophys. Res. Atmos., 122, doi:10.1002/2017JD027429.
  !
  !    Brasseur, G., and A.  Chatel (1983),  paper  presented  at  the  9th  Annual  Meeting  of  the  
  !    European Geophysical Society, Leeds, Great Britain, August 1982.
  !  
  !    Dunne, Eimear M., et al.(2016), Global atmospheric particle formation from CERN CLOUD measurements,
  !    Science 354.6316, 1119-1124.   
  !

  IMPLICIT NONE 

  !----------------------------------------------------
  
  !Global
  DOUBLE PRECISION,INTENT(in) :: t         ! temperature in K 
  DOUBLE PRECISION,INTENT(in) :: satrat    ! saturatio ratio of water (between zero and 1)
  DOUBLE PRECISION,INTENT(in) :: rhoa      ! sulfuric acid concentration in 1/cm3
  DOUBLE PRECISION,INTENT(in) :: csi        ! Ion condensation sink (s-1)
  DOUBLE PRECISION,INTENT(in) :: airn      ! Air molecule concentration in (cm-3)
  DOUBLE PRECISION,INTENT(in) :: ipr       ! Ion pair production rate (cm-3 s-1)
  DOUBLE PRECISION,INTENT(out) :: jnuc_n   ! Neutral nucleation rate in 1/cm3s (J>10^-7 1/cm3s)
  DOUBLE PRECISION,INTENT(out) :: ntot_n   ! total number of molecules in the neutral critical cluster
  DOUBLE PRECISION,INTENT(out) :: jnuc_i   ! Charged nucleation rate in 1/cm3s (J>10^-7 1/cm3s)
  DOUBLE PRECISION,INTENT(out) :: ntot_i   ! total number of molecules in the charged critical cluster
  DOUBLE PRECISION,INTENT(out) :: x_n      ! mole fraction of H2SO4 in the neutral critical cluster 
  DOUBLE PRECISION,INTENT(out) :: x_i      ! mole fraction of H2SO4 in the charged critical cluster 
                                           ! (note that x_n=x_i in nucleation regime) 
  DOUBLE PRECISION,INTENT(out) :: na_n     ! sulfuric acid molecules in the neutral critical cluster
  DOUBLE PRECISION,INTENT(out) :: na_i     ! sulfuric molecules in the charged critical cluster
  DOUBLE PRECISION,INTENT(out) :: rc_n     ! radius of the charged critical cluster in nm 
  DOUBLE PRECISION,INTENT(out) :: rc_i     ! radius of the charged critical cluster in nm 
  DOUBLE PRECISION,INTENT(out) :: n_i      ! number of ion pairs in air (cm-3) 
  LOGICAL,INTENT(out)  :: kinetic_n        ! true if kinetic neutral nucleation
  LOGICAL,INTENT(out)  :: kinetic_i        ! true if kinetic ion-induced nucleation
  DOUBLE PRECISION,INTENT(out) :: rhoatres ! treshold concentration of h2so4 (1/cm^3) for neutral kinetic nucleation

  ! Local
  DOUBLE PRECISION :: x           ! mole fraction of H2SO4 in the critical cluster 
  DOUBLE PRECISION :: satratln    ! bounded water saturation ratio for neutral case (between 5.E-6 - 1.0)
  DOUBLE PRECISION :: satratli    ! bounded water saturation ratio for ion-induced case (between 1.E-7 - 0.95)
  DOUBLE PRECISION :: rhoaln      ! bounded concentration of h2so4 for neutral case (between 10^10 - 10^19 m-3)
  DOUBLE PRECISION :: rhoali      ! bounded concentration of h2so4 for ion-induced case (between 10^10 - 10^22 m-3)
  DOUBLE PRECISION :: tln         ! bounded temperature for neutral case (between 165-400 K)
  DOUBLE PRECISION :: tli         ! bounded temperature for ion-induced case (195-400 K)
  DOUBLE PRECISION :: kinrhotresn ! threshold sulfuric acid for neutral kinetic nucleation   
  DOUBLE PRECISION :: kinrhotresi ! threshold sulfuric acid for ion-induced kinetic nucleation 
  DOUBLE PRECISION :: jnuc_i1     ! Ion-induced rate for n_i=1 cm-3 
  DOUBLE PRECISION :: xloss       ! Ion loss rate 
  DOUBLE PRECISION :: recomb      ! Ion-ion recombination rate 

  !--- 0) Initializations:

  kinetic_n=.FALSE.
  kinetic_i=.FALSE.
  jnuc_n=0.0
  jnuc_i=0.0
  ntot_n=0.0
  ntot_i=0.0
  na_n=0.0
  na_i=0.0
  rc_n=0.0
  rc_i=0.0
  x=0.0
  x_n=0.0
  x_i=0.0
  satratln=satrat
  satratli=satrat
  rhoaln=rhoa
  rhoali=rhoa
  tln=t
  tli=t
  n_i=0.0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Boundary values according to parameterization limits   

  !Temperature bounds
  IF(t.LE.165.) THEN
     PRINT *,'Warning: temperature < 165.0 K, using 165.0 K in neutral nucleation calculation'
     tln=165.0
  END IF
  IF(t.LE.195.) THEN
     PRINT *,'Warning: temperature < 195.0 K, using 195.0 K in ion-induced nucleation calculation'
     tli=195.0
  END IF
  IF(t.GE.400.) THEN
     PRINT *,'Warning: temperature > 400. K, using 400. K in nucleation calculations'
     tln=400.
     tli=400.
  END IF
  
  ! Saturation ratio bounds
  IF(satrat.LT.1.E-7) THEN
     PRINT *,'Warning: saturation ratio of water < 1.e-7, using 1.e-7 in ion-induced nucleation calculation'
     satratli=1.E-7
  END IF
  IF(satrat.LT.1.E-5) THEN
     PRINT *,'Warning: saturation ratio of water < 1.e-5, using 1.e-5 in neutral nucleation calculation'
     satratln=1.E-5
  END IF
  IF(satrat.GT.0.95) THEN
     PRINT *,'Warning: saturation ratio of water > 0.95, using 0.95 in ion-induced nucleation calculation'
     satratli=0.95
  END IF
  IF(satrat.GT.1.0) THEN
     PRINT *,'Warning: saturation ratio of water > 1 using 1 in neutral nucleation calculation'
     satratln=1.0
  END IF
  
  ! Sulfuric acid concentration bounds
  IF(rhoa.LE.1.e4) THEN
     PRINT *,'Warning: sulfuric acid < 1e4 1/cm3, using 1e4 1/cm3 in nucleation calculation'
     rhoaln=1.e4
     rhoali=1.e4
  END IF
  IF(rhoa.GT.1.e13) THEN
     PRINT *,'Warning: sulfuric acid > 1e13 1/cm3, using 1e13 1/cm3 in neutral nucleation calculation'
     rhoaln=1.e13
  END IF
  IF(rhoa.GT.1.e16) THEN
     PRINT *,'Warning: sulfuric acid concentration > 1e16 1/cm3, using 1e16 1/cm3 in ion-induced nucleation calculation'
     rhoali=1.e16
  END IF
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !Critical cluster composition (valid for both cases, bounds not used here) 
  x_n=  7.9036365428891719e-1 - 2.8414059650092153e-3*tln + 1.4976802556584141e-2*LOG(satratln) &
       & - 2.4511581740839115e-4*tln*LOG(satratln) + 3.4319869471066424e-3 *LOG(satratln)**2     &  
       & - 2.8799393617748428e-5*tln*LOG(satratln)**2 + 3.0174314126331765e-4*LOG(satratln)**3 & 
       & - 2.2673492408841294e-6*tln*LOG(satratln)**3 - 4.3948464567032377e-3*LOG(rhoaln)&
       & + 5.3305314722492146e-5*tln*LOG(rhoaln)
  x_i=  7.9036365428891719e-1 - 2.8414059650092153e-3*tli + 1.4976802556584141e-2*LOG(satratli) &
       & - 2.4511581740839115e-4*tli*LOG(satratli) + 3.4319869471066424e-3 *LOG(satratli)**2     &  
       & - 2.8799393617748428e-5*tli*LOG(satratli)**2 + 3.0174314126331765e-4*LOG(satratli)**3 & 
       & - 2.2673492408841294e-6*tli*LOG(satratli)**3 - 4.3948464567032377e-3*LOG(rhoali)&
       & + 5.3305314722492146e-5*tli*LOG(rhoali)
       
  x_n=MIN(MAX(x_n,DBLE(1.e-30)),DBLE(1.)) 
  x_i=MIN(MAX(x_i,DBLE(1.e-30)),DBLE(1.)) 
  
  !Neutral nucleation
  
  !Kinetic limit check
  IF (satratln .GE. 1.e-2 .AND. satratln .LE. 1.) THEN
     kinrhotresn=EXP(7.8920778706888086e+1 + 7.3665492897447082*satratln - 1.2420166571163805e+4/tln &
          & + (-6.1831234251470971e+2*satratln)/tln - 2.4501159970109945e-2*tln &
          & -1.3463066443605762e-2*satratln*tln + 8.3736373989909194e-06*tln**2   &
          & -1.4673887785408892*LOG(satratln) + (-3.2141890006517094e+1*LOG(satratln))/tln & 
          & + 2.7137429081917556e-3*tln*LOG(satratln)) !1/cm3     
     IF(kinrhotresn.LT.rhoaln) kinetic_n=.TRUE.
  ENDIF
  IF (satratln .GE. 1.e-4  .AND. satratln .LT. 1.e-2) THEN     
     kinrhotresn=EXP(7.9074383049843647e+1 - 2.8746005462158347e+1*satratln - 1.2070272068458380e+4/tln &
          & + (-5.9205040320056632e+3*satratln)/tln - 2.4800372593452726e-2*tln &
          & -4.3983007681295948e-2*satratln*tln + 2.5943854791342071e-5*tln**2   &
          & -2.3141363245211317*LOG(satratln) + (9.9186787997857735e+1*LOG(satratln))/tln & 
          & + 5.6819382556144681e-3*tln*LOG(satratln)) !1/cm3
     IF(kinrhotresn.LT.rhoaln) kinetic_n=.TRUE.
  ENDIF
  IF (satratln .GE. 5.e-6  .AND. satratln .LT. 1.e-4) THEN
     kinrhotresn=EXP(8.5599712000361677e+1 + 2.7335119660796581e+3*satratln - 1.1842350246291651e+4/tln &
          & + (-1.2439843468881438e+6*satratln)/tln - 5.4536964974944230e-2*tln &
          & + 5.0886987425326087*satratln*tln + 7.1964722655507067e-5*tln**2   &
          & -2.4472627526306372*LOG(satratln) + (1.7561478001423779e+2*LOG(satratln))/tln & 
          & + 6.2640132818141811e-3*tln*LOG(satratln)) !1/cm3
     IF(kinrhotresn.LT.rhoaln) kinetic_n=.TRUE. 
  ENDIF
  
  IF(kinetic_n) THEN    
     ! Dimer formation rate
     jnuc_n=1.E6*(2.*0.3E-9)**2.*SQRT(8.*3.141593*1.38E-23*(1./(1.661e-27*98.07)+1./(1.661e-27*98.07)))/2.*SQRT(t)*rhoa**2.
     ntot_n=1. !set to 1 
     na_n=1.   ! The critical cluster contains one molecule, but the produced cluster contains 2 molecules
     x_n=na_n/ntot_n  ! so also set this to 1
     rc_n=0.3E-9
  ELSE
     jnuc_n= 2.1361182605986115e-1 + 3.3827029855551838 *tln -3.2423555796175563e-2*tln**2 +  &
          &  7.0120069477221989e-5*tln**3 +8.0286874752695141/x_n -  &
          &  2.6939840579762231e-1*LOG(satratln) +1.6079879299099518*tln*LOG(satratln) -  &
          &  1.9667486968141933e-2*tln**2*LOG(satratln) +  &
          &  5.5244755979770844e-5*tln**3*LOG(satratln) + (7.8884704837892468*LOG(satratln))/x_n +  &
          &  4.6374659198909596*LOG(satratln)**2 - 8.2002809894792153e-2*tln*LOG(satratln)**2 +  &
          &  8.5077424451172196e-4*tln**2*LOG(satratln)**2 +  &
          &  (-2.6518510168987462e-6)*tln**3*LOG(satratln)**2 +  &
          &  (-1.4625482500575278*LOG(satratln)**2)/x_n - 5.2413002989192037e-1*LOG(satratln)**3 +  &
          &  5.2755117653715865e-3*tln*LOG(satratln)**3 +  &
          &  (-2.9491061332113830e-6)*tln**2*LOG(satratln)**3 +  &
          &  (-2.4815454194486752e-8)*tln**3*LOG(satratln)**3 +  &
          &  (-5.2663760117394626e-2*LOG(satratln)**3)/x_n +  &
          &  1.6496664658266762*LOG(rhoaln) +  &
          &  (-8.0809397859218401e-1)*tln*LOG(rhoaln) +  &
          &  8.9302927091946642e-3*tln**2*LOG(rhoaln) +  &
          &  (-1.9583649496497497e-5)*tln**3*LOG(rhoaln) +  &
          &  (-8.9505572676891685*LOG(rhoaln))/x_n +  &
          &  (-3.0025283601622881e+1)*LOG(satratln)*LOG(rhoaln) +  &
          &  3.0783365644763633e-1*tln*LOG(satratln)*LOG(rhoaln) +  &
          &  (-7.4521756337984706e-4)*tln**2*LOG(satratln)*LOG(rhoaln) +  &
          &  (-5.7651433870681853e-7)*tln**3*LOG(satratln)*LOG(rhoaln) +  &
          &  (1.2872868529673207*LOG(satratln)*LOG(rhoaln))/x_n +  &
          &  (-6.1739867501526535e-1)*LOG(satratln)**2*LOG(rhoaln) +  &
          &  7.2347385705333975e-3*tln*LOG(satratln)**2*LOG(rhoaln) +  &
          &  (-3.0640494530822439e-5)*tln**2*LOG(satratln)**2*LOG(rhoaln) +  &
          &  6.5944609194346214e-8*tln**3*LOG(satratln)**2*LOG(rhoaln) +  &
          &  (-2.8681650332461055e-2*LOG(satratln)**2*LOG(rhoaln))/x_n +  &
          &  6.5213802375160306*LOG(rhoaln)**2 +  &
          &  (-4.7907162004793016e-2)*tln*LOG(rhoaln)**2 +  &
          &  (-1.0727890114215117e-4)*tln**2*LOG(rhoaln)**2 +  &
          &  5.6401818280534507e-7*tln**3*LOG(rhoaln)**2 +  &
          &  (5.4113070888923009e-1*LOG(rhoaln)**2)/x_n +  &
          &  5.2062808476476330e-1*LOG(satratln)*LOG(rhoaln)**2 +  &
          &  (-6.0696882500824584e-3)*tln*LOG(satratln)*LOG(rhoaln)**2 +  &
          &  2.3851383302608477e-5*tln**2*LOG(satratln)*LOG(rhoaln)**2 +  &
          &  (-1.5243837103067096e-8)*tln**3*LOG(satratln)*LOG(rhoaln)**2 +  &
          &  (-5.6543192378015687e-2*LOG(satratln)*LOG(rhoaln)**2)/x_n +  &
          &  (-1.1630806410696815e-1)*LOG(rhoaln)**3 +  &
          &  1.3806404273119610e-3*tln*LOG(rhoaln)**3 +  &
          &  (-2.0199865087650833e-6)*tln**2*LOG(rhoaln)**3 +  &
          &  (-3.0200284885763192e-9)*tln**3*LOG(rhoaln)**3 +  &
          &  (-6.9425267104126316e-3*LOG(rhoaln)**3)/x_n
     jnuc_n=EXP(jnuc_n) 
     
     ntot_n =-3.5863435141979573e-3 - 1.0098670235841110e-1 *tln + 8.9741268319259721e-4 *tln**2 - 1.4855098605195757e-6*tln**3 &
          &   - 1.2080330016937095e-1/x_n + 1.1902674923928015e-3*LOG(satratln) - 1.9211358507172177e-2*tln*LOG(satratln) +  &
          &   2.4648094311204255e-4*tln**2*LOG(satratln) - 7.5641448594711666e-7*tln**3*LOG(satratln) +  &
          &   (-2.0668639384228818e-02*LOG(satratln))/x_n - 3.7593072011595188e-2*LOG(satratln)**2 + &
          &   9.0993182774415718e-4 *tln*LOG(satratln)**2 +&
          &   (-9.5698412164297149e-6)*tln**2*LOG(satratln)**2 + 3.7163166416110421e-8*tln**3*LOG(satratln)**2 +  &
          &   (1.1026579525210847e-2*LOG(satratln)**2)/x_n + 1.1530844115561925e-2 *LOG(satratln)**3 +  &
          &   (-1.8083253906466668e-4) *tln*LOG(satratln)**3 + 8.0213604053330654e-7*tln**2*LOG(satratln)**3 +  &
          &   (-8.5797885383051337e-10)*tln**3*LOG(satratln)**3 + (1.0243693899717402e-3*LOG(satratln)**3)/x_n +  &
          &   (-1.7248695296299649e-2)*LOG(rhoaln) + 1.1294004162437157e-2*tln*LOG(rhoaln) +  &
          &   (-1.2283640163189278e-4)*tln**2*LOG(rhoaln) + 2.7391732258259009e-7*tln**3*LOG(rhoaln) +  &
          &   (6.8505583974029602e-2*LOG(rhoaln))/x_n +2.9750968179523635e-1*LOG(satratln)*LOG(rhoaln) +  &
          &   (-3.6681154503992296e-3) *tln*LOG(satratln)*LOG(rhoaln) + 1.0636473034653114e-5*tln**2*LOG(satratln)*LOG(rhoaln) + &
          &   5.8687098466515866e-9*tln**3*LOG(satratln)*LOG(rhoaln) + (-5.2028866094191509e-3*LOG(satratln)*LOG(rhoaln))/x_n + &
          &   7.6971988880587231e-4*LOG(satratln)**2*LOG(rhoaln) - 2.4605575820433763e-5*tln*LOG(satratln)**2*LOG(rhoaln) +  &
          &   2.3818484400893008e-7*tln**2*LOG(satratln)**2*LOG(rhoaln) +  &
          &   (-8.8474102392445200e-10)*tln**3*LOG(satratln)**2*LOG(rhoaln) +  &
          &   (-1.6640566678168968e-4*LOG(satratln)**2*LOG(rhoaln))/x_n - 7.7390093776705471e-2*LOG(rhoaln)**2 +  &
          &   5.8220163188828482e-4*tln*LOG(rhoaln)**2 + 1.2291679321523287e-6*tln**2*LOG(rhoaln)**2 +  &
          &   (-7.4690997508075749e-9)*tln**3*LOG(rhoaln)**2 + (-5.6357941220497648e-3*LOG(rhoaln)**2)/x_n +  &
          &   (-4.7170109625089768e-3)*LOG(satratln)*LOG(rhoaln)**2 + 6.9828868534370193e-5*tln*LOG(satratln)*LOG(rhoaln)**2 + &
          &   (-3.1738912157036403e-7)*tln**2*LOG(satratln)*LOG(rhoaln)**2 +  &
          &   2.3975538706787416e-10*tln**3*LOG(satratln)*LOG(rhoaln)**2 +  &
          &   (4.2304213386288567e-4*LOG(satratln)*LOG(rhoaln)**2)/x_n + 1.3696520973423231e-3*LOG(rhoaln)**3 +  &
          &   (-1.6863387574788199e-5)*tln*LOG(rhoaln)**3 + 2.7959499278844516e-8*tln**2*LOG(rhoaln)**3 +  &
          &   3.9423927013227455e-11*tln**3*LOG(rhoaln)**3 + (8.6136359966337272e-5*LOG(rhoaln)**3)/x_n
     ntot_n=EXP(ntot_n)
     
     rc_n=EXP(-22.378268374023630+0.44462953606125100*x_n+0.33499495707849131*LOG(ntot_n)) !in meters
     
     na_n=x_n*ntot_n
     IF (na_n .LT. 1.) THEN
        PRINT *, 'Warning: number of acid molecules < 1 in nucleation regime, setting na_n=1'
        na_n=1.0
     ENDIF
  ENDIF
  
  ! Set the neutral nucleation rate to 0.0 if less than 1.0e-7      
  IF(jnuc_n.LT.1.e-7) THEN
     jnuc_n=0.0
  ENDIF
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Threshold neutral nucleation rate (j > 1/cm3s) parameterization (can be commented out if not needed)
  IF (tln .GE. 310.) THEN
     rhoatres=EXP(-2.8220714121794250 + 1.1492362322651116e+1*satratln -3.3034839106184218e+3/tln &
          & + (-7.1828571490168133e+2*satratln)/tln + 1.4649510835204091e-1*tln &
          & -3.0442736551916524e-2*satratln*tln -9.3258567137451497e-5*tln**2   &
          & -1.1583992506895649e+1*LOG(satratln) + (1.5184848765906165e+3*LOG(satratln))/tln & 
          & + 1.8144983916747057e-2*tln*LOG(satratln)) !1/cm3
  ENDIF
  IF (tln .GT. 190. .AND. tln .LT. 310.) THEN
     rhoatres=EXP(-3.1820396091231999e+2 + 7.2451289153199676*satratln + 2.6729355170089486e+4/tln &
          & + (-7.1492506076423069e+2*satratln)/tln + 1.2617291148391978*tln &
          & - 1.6438112080468487e-2*satratln*tln -1.4185518234553220e-3*tln**2   &
          & -9.2864597847386694*LOG(satratln) + (1.2607421852455602e+3*LOG(satratln))/tln & 
          & + 1.3324434472218746e-2*tln*LOG(satratln)) !1/cm3
  ENDIF
  ! to prevent hole in rhoatres calculation set upper limit to 190 K
  ! if (tln .lt. 185. .and. tln .gt. 155.) then   
  IF (tln .LE. 190. .AND. tln .GT. 155.) THEN
     rhoatres=1.1788859232398459e+5 - 1.0244255702550814e+4*satratln + &
          & 4.6815029684321962e+3*satratln**2 -1.6755952338499657e+2*tln
  ENDIF
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  ! Ion-induced nucleation parameterization 
 
  IF(ipr.GT.0.0) THEN ! if the ion production rate is above zero
     
     ! Calculate the ion induced nucleation rate wrt. concentration of 1 ion/cm3
     
     kinrhotresi = 5.3742280876674478e1  - 6.6837931590012266e-3 *LOG(satratli)**(-2) &
          & - 1.0142598385422842e-01 * LOG(satratli)**(-1) - 6.4170597272606873e+00 * LOG(satratli) &
          & - 6.4315798914824518e-01 * LOG(satratli)**2 - 2.4428391714772721e-02 * LOG(satratli)**3 & 
          & - 3.5356658734539019e-04 * LOG(satratli)**4 + 2.5400015099140506e-05 * tli * LOG(satratli)**(-2) & 
          & - 2.7928900816637790e-04 * tli * LOG(satratli)**(-1) + 4.4108573484923690e-02 * tli * LOG(satratli) &
          & + 6.3943789012475532e-03 * tli * LOG(satratli)**(2) + 2.3164296174966580e-04 * tli * LOG(satratli)**(3) &
          & + 3.0372070669934950e-06 * tli * LOG(satratli)**4 + 3.8255873977423475e-06 * tli**2 * LOG(satratli)**(-1) &
          & - 1.2344793083561629e-04 * tli**2 * LOG(satratli) - 1.7959048869810192e-05 * tli**2 * LOG(satratli)**(2) &
          & - 3.2165622558722767e-07 * tli**2 * LOG(satratli)**3 - 4.7136923780988659e-09 * tli**3 * LOG(satratli)**(-1) &
          & + 1.1873317184482216e-07 * tli**3 * LOG(satratli) + 1.5685860354866621e-08 * tli**3 * LOG(satratli)**2 &
          & - 1.4329645891059557e+04 * tli**(-1) + 1.3842599842575321e-01 * tli &
          & - 4.1376265912842938e-04 * tli**(2) + 3.9147639775826004e-07 * tli**3
     
     kinrhotresi=EXP(kinrhotresi) !1/cm3
     
     IF(kinrhotresi.LT.rhoali) kinetic_i=.TRUE.
     
     IF(kinetic_i) THEN    
        jnuc_i1=1.0E6*(0.3E-9 + 0.487E-9)**2.*SQRT(8.*3.141593*1.38E-23*(1./(1.661e-27*98.07)+1./(1.661e-27*98.07)))*&
             &sqrt(tli)*rhoali !1/cm3s  
        ntot_i=1. !set to 1 
        na_i=1.
        x_i=na_i/ntot_i  ! so also set this to 1
        rc_i=0.487E-9
     ELSE
        jnuc_i1 = 3.0108954259038608e+01+tli*6.1176722090512577e+01+(tli**2)*8.7240333618891663e-01+(tli**3)* &
             & (-4.6191788649375719e-03)+(tli**(-1))*8.3537059107024481e-01 + & 
             & (1.5028549216690628e+01+tli*(-1.9310989753720623e-01)+(tli**2)*8.0155514634860480e-04+(tli**3)* &
             & (-1.0832730707799128e-06)+(tli**(-1))*1.7577660457989019)*(LOG(satratli)**(-2)) + &
             & (-2.0487870170216488e-01 +  tli * 1.3263949252910405e-03 +  (tli**2) * (-8.4195688402450274e-06) + &
             & (tli**3) * 1.6154895940993287e-08 + (tli**(-1)) * 3.8734212545203874e+01) * (LOG(satratli)**(-2) * LOG(rhoali)) + &
             & (1.4955918863858371 +  tli * 9.2290004245522454e+01 +  (tli**2) * (-8.9006965195392618e-01) + &
             & (tli**3) * 2.2319123411013099e-03 + (tli**(-1)) * 4.0180079996840852e-03) * &
             & (LOG(satratli)**(-1) * LOG(rhoali)**(-1)) + &
             & (7.9018031228561085 +  tli * (-1.1649433968658949e+01) +  (tli**2) * 1.1400827854910951e-01 + &
             & (tli**3) * (-3.1941526492127755e-04) + (tli**(-1)) * (-3.7662115740271446e-01)) * (LOG(satratli)**(-1)) + &
             & (1.5725237111225979e+02 +  tli * (-1.0051649979836277) +  (tli**2) * 1.1866484014507624e-03 + &
             & (tli**3) * 7.3557614998540389e-06 + (tli**(-1)) * 2.6270197023115189) * (LOG(satratli)**(-1) * LOG(rhoali)) + &
             & (-1.6973840122470968e+01 +  tli * 1.1258423691432135e-01 +  (tli**2) * (-2.9850139351463793e-04) + (tli**3) * &
             & 1.4301286324827064e-07 + (tli**(-1)) * 1.3163389235253725e+01) * (LOG(satratli)**(-1) * LOG(rhoali)**2) + &
             & (-1.0399591631839757 +  tli * 2.7022055588257691e-03 +  (tli**2) * (-2.1507467231330936e-06) + (tli**3) * &
             & 3.8059489037584171e-10 + (tli**(-1)) * 1.5000492788553410e+02) * (LOG(satratli)**(-1) * LOG(rhoali)**3) + &
             & (1.2250990965305315 +  tli * 3.0495946490079444e+01 +  (tli**2) * 2.1051563135187106e+01 + (tli**3) * &
             & (-8.2200682916580878e-02) + (tli**(-1)) * 2.9965871386685029e-02) * (LOG(rhoali)**(-2)) + &
             & (4.8281605955680433 +  tli * 1.7346551710836445e+02 +  (tli**2) * (-1.0113602140796010e+01) + (tli**3) * &
             & 3.7482518458685089e-02 + (tli**(-1)) * (-1.4449998158558205e-01)) * (LOG(rhoali)**(-1)) + &
             & (2.3399230964451237e+02 +  tli * (-2.3099267235261948e+01) +  (tli**2) * 8.0122962140916354e-02 + &
             & (tli**3) * 6.1542576994557088e-05 + (tli**(-1)) * 5.3718413254843007) * (LOG(rhoali)) + &
             & (1.0299715519499360e+02 +  tli * (-6.4663357203364136e-02) +  (tli**2) * (-2.0487150565050316e-03) + &
             & (tli**3) * 8.7935289055530897e-07 + (tli**(-1)) * 3.6013204601215229e+01) * (LOG(rhoali)**2) + &
             & (-3.5452115439584042 +  tli * 1.7083445731159330e-02 +  (tli**2) * (-1.2552625290862626e-05) + (tli**3) * &
             & 1.2968447449182847e-09 + (tli**(-1)) * 1.5748687512056560e+02) * (LOG(rhoali)**3) + &
             & (2.2338490119517975 +  tli * 1.0229410216045540e+02 +  (tli**2) * (-3.2103611955174052) + (tli**3) * &
             & 1.3397152304977591e-02 + (tli**(-1)) * (-2.4155187776460030e-02)) * (LOG(satratli)* LOG(rhoali)**(-2)) + &
             & (3.7592282990713963 +  tli * (-1.5257988769009816e+02) +  (tli**2) * 2.6113805420558802 + (tli**3) * &
             & (-9.0380721653694363e-03) + (tli**(-1)) * (-1.3974197138171082e-01)) * (LOG(satratli)* LOG(rhoali)**(-1)) + &
             & (1.8293600730573988e+01 +  tli * 1.8344728606002992e+01 +  (tli**2) * (-4.0063363221106751e-01) + (tli**3) &
             & * 1.4842749371258522e-03 + (tli**(-1)) * 1.1848846003282287) * (LOG(satratli)) + &
             & (-1.7634531623032314e+02 +  tli * 4.9011762441271278 +  (tli**2) * (-1.3195821562746339e-02) + (tli**3) * &
             & (-2.8668619526430859e-05) + (tli**(-1)) * (-2.9823396976393551e-01)) * (LOG(satratli)* LOG(rhoali)) + &
             & (-3.2944043694275727e+01 +  tli * 1.2517571921051887e-01 +  (tli**2) * 8.3239769771186714e-05 + (tli**3) * &
             & 2.8191859341519507e-07 + (tli**(-1)) * (-2.7352880736682319e+01)) * (LOG(satratli)* LOG(rhoali)**2) + &
             & (-1.1451811137553243 +  tli * 2.0625997485732494e-03 +  (tli**2) * (-3.4225389469233624e-06) + (tli**3) * &
             & 4.4437613496984567e-10 + (tli**(-1)) * 1.8666644332606754e+02) * (LOG(satratli)* LOG(rhoali)**3) + &
             & (3.2270897099493567e+01 +  tli * 7.7898447327513687e-01 +  (tli**2) * (-6.5662738484679626e-03) + (tli**3) * &
             & 3.7899330796456790e-06 + (tli**(-1)) * 7.1106427501756542e-01) * (LOG(satratli)**2 * LOG(rhoali)**(-1)) + &
             & (-2.8901906781697811e+01 +  tli * (-1.5356398793054860) +  (tli**2) * 1.9267271774384788e-02 + (tli**3) * &
             & (-5.3886270475516162e-05) + (tli**(-1)) * 5.0490415975693426e-01) * (LOG(satratli)**2) + &
             & (3.3365683645733924e+01 +  tli * (-3.6114561564894537e-01) +  (tli**2) * 9.2977354471929262e-04 + (tli**3) * &
             & 1.9549769069511355e-07 + (tli**(-1)) * (-8.8865930095112855)) * (LOG(satratli)**2 * LOG(rhoali)) + &
             & (2.4592563042806375 +  tli * (-8.3227071743101084e-03) +  (tli**2) * 8.2563338043447783e-06 + (tli**3) * &
             & (-8.4374976698593496e-09) + (tli**(-1)) * (-2.0938173949893473e+02)) * (LOG(satratli)**2 * LOG(rhoali)**2) + &
             & (4.4099823444352317e+01 +  tli * 2.5915665826835252 +  (tli**2) * (-1.6449091819482634e-02) + (tli**3) * & 
             & 2.6797249816144721e-05 + (tli**(-1)) * 5.5045672663909995e-01)* satratli
        jnuc_i1=EXP(jnuc_i1)
        
        ntot_i = ABS((-4.8324296064013375e+04 +  tli * 5.0469120697428906e+02 +  (tli**2) * (-1.1528940488496042e+00) + &
             & (tli**(-1)) * (-8.6892744676239192e+02) + (tli**(3)) * 4.0030302028120469e-04) + &
             & (-6.7259105232039847e+03 +  tli * 1.9197488157452008e+02 +  (tli**2) * (-1.3602976930126354e+00) + &
             & (tli**(-1)) * (-1.1212637938360332e+02) + (tli**(3)) * 2.8515597265933207e-03) * &
             & LOG(satratli)**(-2) * LOG(rhoali)**(-2) + &
             & (2.6216455217763342e+02 +  tli * (-2.3687553252750821e+00) +  (tli**2) * 7.4074554767517521e-03 + &
             & (tli**(-1)) * (-1.9213956820114927e+03) + (tli**(3)) * (-9.3839114856129453e-06)) * LOG(satratli)**(-2) + &
             & (3.9652478944137344e+00 +  tli * 1.2469375098256536e-02 +  (tli**2) * (-9.9837754694045633e-05) + (tli**(-1)) * &
             & (-5.1919499210175138e+02) + (tli**(3)) * 1.6489001324583862e-07) * LOG(satratli)**(-2) * LOG(rhoali) + &
             & (2.4975714429096206e+02 +  tli * 1.7107594562445172e+02 +  (tli**2) * (-7.8988711365135289e-01) + (tli**(-1)) * &
             & (-2.2243599782483177e+01) + (tli**(3)) * (-1.6291523004095427e-04)) * LOG(satratli)**(-1) * LOG(rhoali)**(-2) + &
             & (-8.9270715592533611e+02 +  tli * 1.2053538883338946e+02 +  (tli**2) * (-1.5490408828541018e+00) + (tli**(-1)) * &
             & (-1.1243275579419826e+01) + (tli**(3)) * 4.8053105606904655e-03) * LOG(satratli)**(-1) * LOG(rhoali)**(-1) + &
             & (7.6426441642091631e+03 +  tli * (-7.1785462414656578e+01) +  (tli**2) * 2.3851864923199523e-01 + (tli**(-1)) * &
             & 8.5591775688708395e+01 + (tli**(3)) * (-3.7000473243342858e-04)) * LOG(satratli)**(-1) + &
             & (-5.1516826398607911e+01 +  tli * 9.1385720811460558e-01 +  (tli**2) * (-3.5477100262158974e-03) + &
             & (tli**(-1)) * 2.7545544507625586e+03 + (tli**(3)) * 5.4708262093640928e-06) * LOG(satratli)**(-1) * LOG(rhoali) + &
             & (-3.0386767129196176e+02 +  tli * (-1.1033438883583569e+04) +  (tli**2) * 8.1296859732896067e+01 + (tli**(-1)) * &
             & 1.2625883141097162e+01 + (tli**(3)) * (-1.2728497822219101e-01)) * LOG(rhoali)**(-2) + &
             & (-3.3763494256461472e+03 +  tli * 3.1916579136391006e+03 +  (tli**2) * (-2.7234339474441143e+01) + (tli**(-1)) * &
             & (-2.1897653262707397e+01) + (tli**(3)) * 5.1788505812259071e-02) * LOG(rhoali)**(-1) + &
             & (-1.8817843873687068e+03 +  tli * 4.3038072285882070e+00 +  (tli**2) * 6.6244087689671860e-03 + (tli**(-1)) * &
             & (-2.7133073605696295e+03) + (tli**(3)) * (-1.7951557394285043e-05)) * LOG(rhoali) + &
             & (-1.7668827539244447e+02 +  tli * 4.8160932330629913e-01 +  (tli**2) * (-6.3133007671100293e-04) + (tli**(-1)) * &
             & 2.5631774669873157e+04 + (tli**(3)) * 4.1534484127873519e-07) * LOG(rhoali)**(2) + &
             & (-1.6661835889222382e+03 +  tli * 1.3708900504682877e+03 +  (tli**2) * (-1.7919060052198969e+01) + (tli**(-1)) * &
             & (-3.5145029804436405e+01) + (tli**(3)) * 5.1047240947371224e-02) * LOG(satratli)* LOG(rhoali)**(-2) + &
             & (1.0843549363030939e+04 +  tli * (-7.3557073636139577e+01) +  (tli**2) * 1.2054625131778862e+00 + (tli**(-1)) * &
             & 1.9358737917864391e+02 + (tli**(3)) * (-4.2871620775911338e-03)) * LOG(satratli)* LOG(rhoali)**(-1) + &
             & (-2.4269802549752835e+03 +  tli * 1.1348265061941714e+01 +  (tli**2) * (-5.0430423939495157e-02) + (tli**(-1)) * &
             & 2.3709874548950634e+03 + (tli**(3)) * 1.4091851828620244e-04) * LOG(satratli) + &
             & (5.2745372575251588e+02 +  tli * (-2.6080675912627314e+00) +  (tli**2) * 5.6902218056670145e-03 + (tli**(-1)) * &
             & (-3.2149319482897838e+04) + (tli**(3)) * (-5.4121996056745853e-06)) * LOG(satratli)* LOG(rhoali) + &
             & (-1.6401959518360403e+01 +  tli * 2.4322962162439640e-01 +  (tli**2) * 1.1744366627725344e-03 + (tli**(-1)) * &
             & (-8.2694427518413195e+03) + (tli**(3)) * (-5.0028379203873102e-06))* LOG(satratli)**(2) + &
             & (-2.7556572017167782e+03 +  tli * 4.9293344495058264e+01 +  (tli**2) * (-2.6503456520676050e-01) + (tli**(-1)) * &
             & 1.2130698030982167e+03 + (tli**(3)) * 4.3530610668042957e-04)* LOG(satratli)**2 * LOG(rhoali)**(-1) + &
             & (-6.3419182228959192e+00 +  tli * 4.0636212834605827e-02 +  (tli**2) * (-1.0450112687842742e-04) + (tli**(-1)) * &
             & 3.1035882189759656e+02 + (tli**(3)) * 9.4328418657873500e-08)* LOG(satratli)**(-3) + &
             & (3.0189213304689042e+03 +  tli * (-2.3804654203861684e+01) +  (tli**2) * 6.8113013411972942e-02 + (tli**(-1)) * & 
             & 6.3112071081188913e+02 + (tli**(3)) * (-9.4460854261685723e-05))* (satratli) * LOG(rhoali) + &
             & (1.1924791930673702e+04 +  tli * (-1.1973824959206000e+02) +  (tli**2) * 1.6888713097971020e-01 + (tli**(-1)) * &
             & 1.8735938211539585e+02 + (tli**(3)) * 5.0974564680442852e-04)* (satratli) + &
             & (3.6409071302482083e+01 +  tli * 1.7919859306449623e-01 +  (tli**2) * (-1.0020116255895206e-03) + (tli**(-1)) * &
             & (-8.3521083354432303e+03) + (tli**(3)) * 1.5879900546795635e-06)* satratli * LOG(rhoali)**(2))
          
          
        rc_i = (-3.6318550637865524e-08 +  tli * 2.1740704135789128e-09   +  (tli**2) * &
             & (-8.5521429066506161e-12) + (tli**3) * (-9.3538647454573390e-15)) + &
             & (2.1366936839394922e-08 +  tli * (-2.4087168827395623e-10) +  (tli**2) * 8.7969869277074319e-13 + &
             & (tli**3) * (-1.0294466881303291e-15))* LOG(satratli)**(-2) * LOG(rhoali)**(-1) + &
             & (-7.7804007761164303e-10 +  tli * 1.0327058173517932e-11 +  (tli**2) * (-4.2557697639692428e-14) + &
             & (tli**3) * 5.4082507061618662e-17)* LOG(satratli)**(-2) + &
             & (3.2628927397420860e-12 +  tli * (-7.6475692919751066e-14) +  (tli**2) * 4.1985816845259788e-16 + &
             & (tli**3) * (-6.2281395889592719e-19))* LOG(satratli)**(-2) * LOG(rhoali) + &
             & (2.0442205540818555e-09 +  tli * 4.0441858911249830e-08 +  (tli**2) * (-3.3423487629482825e-10) + &
             & (tli**3) * 6.8000404742985678e-13)* LOG(satratli)**(-1) * LOG(rhoali)**(-2) + &
             & (1.8381489183824627e-08 +  tli * (-8.9853322951518919e-09) +  (tli**2) * 7.5888799566036185e-11 + &
             & (tli**3) * (-1.5823457864755549e-13))* LOG(satratli)**(-1) * LOG(rhoali)**(-1) + &
             & (1.1795760639695057e-07 +  tli * (-8.1046722896375875e-10) +  (tli**2) * 9.1868604369041857e-14 + &
             & (tli**3) * 4.7882428237444610e-15)* LOG(satratli)**(-1) + &
             & (-4.4028846582545952e-09 +  tli * 4.6541269232626618e-11 +  (tli**2) * (-1.1939929984285194e-13) + &
             & (tli**3) * 2.3602037016614437e-17)* LOG(satratli)**(-1) * LOG(rhoali) + &
             & (2.7885056884209128e-11 +  tli * (-4.5167129624119121e-13) +  (tli**2) * 1.6558404997394422e-15 + &
             & (tli**3) * (-1.2037336621218054e-18))* LOG(satratli)**(-1) * LOG(rhoali)**2 + &
             & (-2.3719627171699983e-09 +  tli * (-1.5260127909292053e-07) +  (tli**2) * 1.7177017944754134e-09 + &
             & (tli**3) * (-4.7031737537526395e-12))* LOG(rhoali)**(-2) + &
             & (-5.6946433724699646e-09 +  tli * 8.4629788237081735e-09 +  (tli**2) * (-1.7674135187061521e-10) + &
             & (tli**3) * 6.6236547903091862e-13)* LOG(rhoali)**(-1) + &
             & (-2.2808617930606012e-08 +  tli * 1.4773376696847775e-10 +  (tli**2) * (-1.3076953119957355e-13) + &
             & (tli**3) * 2.3625301497914000e-16)* LOG(rhoali) + &
             & (1.4014269939947841e-10 +  tli * (-2.3675117757377632e-12) +  (tli**2) * 5.1514033966707879e-15 + &
             & (tli**3) * (-4.8864233454747856e-18))* LOG(rhoali)**2 + &
             & (6.5464943868885886e-11 +  tli * 1.6494354816942769e-08 +  (tli**2) * (-1.7480097393483653e-10) + &
             & (tli**3) * 4.7460075628523984e-13)* LOG(satratli)* LOG(rhoali)**(-2) + &
             & (8.4737893183927871e-09 +  tli * (-6.0243327445597118e-09) +  (tli**2) * 5.8766070529814883e-11 + &
             & (tli**3) * (-1.4926748560042018e-13))* LOG(satratli)* LOG(rhoali)**(-1) + &
             & (1.0761964135701397e-07 +  tli * (-1.0142496009071148e-09) +  (tli**2) * 2.1337312466519190e-12 + &
             & (tli**3) * 1.6376014957685404e-15)* LOG(satratli) + &
             & (-3.5621571395968670e-09 +  tli * 4.1175339587760905e-11 +  (tli**2) * (-1.3535372357998504e-13) + &
             & (tli**3) * 8.9334219536920720e-17)* LOG(satratli)* LOG(rhoali) + &
             & (2.0700482083136289e-11 +  tli * (-3.9238944562717421e-13) +  (tli**2) * 1.5850961422040196e-15 + &
             & (tli**3) * (-1.5336775610911665e-18))* LOG(satratli)* LOG(rhoali)**2 + &
             & (1.8524255464416206e-09 +  tli * (-2.1959816152743264e-11) +  (tli**2) * (-6.4478119501677012e-14) + &
             & (tli**3) * 5.5135243833766056e-16)* LOG(satratli)**2 * LOG(rhoali)**(-1) + &
             & (1.9349488650922679e-09 +  tli * (-2.2647295919976428e-11) +  (tli**2) * 9.2917479748268751e-14 + &
             & (tli**3) * (-1.2741959892173170e-16))* LOG(satratli)**2 + &
             & (2.1484978031650972e-11 +  tli * (-9.3976642475838013e-14) +  (tli**2) * (-4.8892738002751923e-16) + &
             & (tli**3) * 1.4676120441783832e-18)* LOG(satratli)**2 * LOG(rhoali) + &
             & (6.7565715216420310e-13 +  tli * (-3.5421162549480807e-15) +  (tli**2) * (-3.4201196868693569e-18) + &
             & (tli**3) * 2.2260187650412392e-20)* LOG(satratli)**3 * LOG(rhoali)
                    
        na_i=x_i*ntot_i
        IF (na_i .LT. 1.) THEN
           PRINT *, 'Warning: number of acid molecules < 1 in nucleation regime, setting na_n=1'
           na_n=1.0
        ENDIF
     ENDIF
    
     jnuc_i=jnuc_i1 
     ! Ion loss rate (1/s)
     xloss=csi+jnuc_i
     
     ! Recombination (here following Brasseur and Chatel, 1983)   
     recomb=6.0e-8*SQRT(300./tli)+6.0e-26*airn*(300./tli)**4
     
     ! Small ion concentration in air (1/cm3) (following Dunne et al., 2016)
     ! max function is to avoid n_i to go practically zero at very high J_ion 
     n_i=MAX(DBLE(0.01),(SQRT(xloss**2.0+4.0*recomb*ipr)-xloss)/(2.0*recomb))
     
     ! Ion-induced nucleation rate
     ! Min function is to ensure that max function above does not cause J_ion to overshoot 
     jnuc_i=MIN(ipr,n_i*jnuc_i1)
     ! Set the ion-induced nucleation rate to 0.0 if less than 1.0e-7      
     IF(jnuc_i.LT.1.e-7) THEN
        jnuc_i=0.0
     ENDIF

  ENDIF
  
END SUBROUTINE newbinapara

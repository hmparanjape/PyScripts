C ****************************************************************************
C                   VPSC 1-site multiphase  - 06/NOV/2009                    *
C ****************************************************************************
C                        BY R.LEBENSOHN & C.TOME                             *
C ****************************************************************************
C                   MULTIPHASE OR MULTIELEMENT VERSION                       *
C ****************************************************************************
C  RICARDO LEBENSOHN                 |  CARLOS TOME                          *
C  MST-8 - LANL - MS G755            |  MST-8 - LANL - MS G755               *
C  LOS ALAMOS - NM 87545 - USA       |  LOS ALAMOS - NM 87545 - USA          *
C  lebenso@lanl.gov                  |  tome@lanl.gov                        *
C ****************************************************************************
C  GENERAL REFERENCE:                                                        *
C  R.A.Lebensohn & C.N.Tome, Acta metall mater 41, 2611 (1993)               *
C ****************************************************************************
C  COPYRIGHT NOTICE                                                          *
C ****************************************************************************
C Portions of this program were prepared by the Regents of the University of *
C California at Los Alamos National Laboratory (the University) under        *
C Contract No.W-7405-ENG-36 with the U.S. Department of Energy (DOE).        *
C The University has certain rights in the program pursuant to the contract  *
C and the program should not be copied or distributed outside your           *
C organization. All rights in the program are reserved by the DOE and the    *
C University. Neither the U.S. Government nor the University makes any       *
C warranty, express or implied, or assumes any liability or responsibility   *
C for the use of this software.                                              *
C ****************************************************************************
C                              # FEATURES #                                  *
C ****************************************************************************
C  - SINGLE PHASE & MULTIPHASE AGGREGATES.                                   *
C  - SC RELATION FOR MIXED RATE-SENSITIVITY EXPONENTS.                       *
C  - ESHELBY AND PRESSURE TENSORS FOR FULLY INCOMPRESSIBLE INCLUSION & HEM.  *
C  - LOCAL CAUCHY STRESSES                                                   *
C  - MIXED BOUNDARY CONDITIONS ON DISPLACEMENT GRADIENT AND STRESS           *
C  - SYMMETRIC BASIS INTERNAL REPRESENTATION OF TENSORS                      *
C  - OPTION FOR ANALYZING ROLLING FCC COMPONENTS (icubcom=1)                 *
C  - SELF & LATENT HARDENING COEFFICIENTS                                    *
C  - 'VOCE' TYPE HARDENING MODULATION FOR EACH SYSTEM (ihardlaw=0) OR        *
C    'MTS'  HARDENING LAW FOR EACH SYSTEM (ihardlaw=1)                       *
C  - PREDOMINANT TWIN REORIENTATION SCHEME (if ihardlaw=0)                   *
C  - COMPOSITE GRAIN SCHEME (if iharlaw=2)                                   *
C  - DIFFERENT INDIVIDUAL GRAIN SHAPES.                                      *
C  - GENERAL TREATMENT OF GRAIN SHAPE EVOLUTION.                             *
C  - INCREMENTAL CONTROL: VON MISES STRAIN, SINGLE STRAIN COMPONENT OR TIME  *
C  - POSTMORTEM CAPABILITY: RECOVERS ARRAYS FROM PREVIOUS RUN                *
C  - POLYCRYSTAL YIELD SURFACE AND LANKFORD COEFFICIENTS CALCULATION         *
C  - VARIABLE DEFORMATION HISTORY ROUTINE.                                   *
C  - OPTION FOR ROTATION COUPLING BETWEEN ORIENTATIONS                       *
C  - ALL 'CANNED' SUBROUTINES ARE FROM 'NUMERICAL RECIPIES' LIBRARY          *
C  - INTERNAL Bunge CONVENTION (phi1,Phi,phi2 -> phi,theta,omega) FOR        *
C    EULER ANGLES. BUT READS ALSO Kocks and Roe ANGLES AS INPUT.             *
C                                                                            *
C ****************************************************************************
C                         # INPUT FILES #                                    *
C ****************************************************************************
C ALWAYS NEEDED:                                                             *
C                                                                            *
C       VPSC7.IN      - general input and parameters                         *
C       "FILECRYS"    - deformation modes and CRSS (one per phase)           *
C       "FILETEXT"    - initial texture file (one per phase)                 *
C       "FILEPROC"    - deformation process to be simulated                  *
C                                                                            *
C SOMETIMES NEEDED:                                                          *
C                                                                            *
C  if ishape>1     "FILEAXES"  - initial orientation of ellipsoids           *
C                               (morphologic texture, one per phase)         *
C  if irecover=1   POSTMORT.IN - initial state from previous run             *
C  if icubcom=1    CUBCOMP.IN  - ideal orientations - fcc rolling            *
C  if ivgvar=1     "FILEHIST"  - strain history                              *
C  if ihardlaw=1   MTS parameters for the material at the end of "FILECRYS"  *
C                                                                            *
C ****************************************************************************
C                         # OUTPUT FILES #                                   *
C ****************************************************************************
C DEPENDING ON CASE AND DRIVER:                                              *
C        RERR.OUT     - convergence history                                  *
C        STATS.OUT    - statistics                                           *
C        ACT_PHn.OUT  - activity of deformation modes                        *
C        TEX_PHn.OUT  - cryst. texture for each phase/elem at dif steps      *
C        MOR_PHn.OUT  - morph. texture for each phase/elem at dif steps      *
C        STR_STR.OUT  - macroscopic stress-strain components                 *
C        RUN_LOG.OUT  - a log containing the conditions and input to the run *
C        STAT_AX.OUT  - statistics on grain axes                             *
C                                                                            *
C     if isave=1    POSTMORT.OUT - final state in grains and PX              *
C     if icauchy=1  CAUCHY.OUT - final local cauchy stresses                 *
C     if icubcom=1  CUBCOMPn.OUT - rolling components for fcc                *
C                                                                            *
C   PARAMETER & COMMON FILE (incl. during compilation): VPSC6.DIM            *
C                                                                            *
C*****************************************************************************

      !program vpsc_deluxe_ss
      SUBROUTINE VPSC_DELUXE_SS(dbar_, sbar_, eps_tot, 
     $     scauchy_, sdeviat_, dsim_, svm_, dvm_, epsvm_, Np, Ni, M,
     $     jobid)
cc
cc
cc          Np: # of maximum processes
cc          Ni: # of maximum steps in a process
cc          M:  # dimension of the stran and stress (hardwired to be 5)
      INCLUDE 'vpsc7.dim'
cc     -- f2py youngung
c      character*5 foutname
      integer Np, Ni, M
      real*8 dbar_(Np,Ni,M), sbar_(Np,Ni,M)
      real*8 eps_tot(Np,Ni,3,3), scauchy_(Np,Ni,3,3)
      real*8 sdeviat_(Np,Ni,3,3),dsim_(Np,Ni,3,3)
      real*8 svm_(Np, Ni),dvm_(Np, Ni),epsvm_(Np, Ni)
c     real*8 gr_(Nphase, Np, Ni, Ng, 4)
      real*8 aux(3,3)
      character*4 jobid
c                                               ! grain xtal orientation
ccccc Cf2py intent(in) foutname
Cf2py intent(in,out,copy) dbar_,sbar_,eps_tot
Cf2py intent(in,out,copy) scauchy_,sdeviat_,dsim_
Cf2py intent(in,out,copy) svm_, dvm_, epsvm_
Cf2py intent(in) jobid
Cf2py integer intent(hide), depend(dbar_) :: Np     = shape(dbar_,0)
Cf2py integer intent(hide), depend(dbar_) :: Ni     = shape(dbar_,1)
Cf2py integer intent(hide), depend(dbar_) :: M      = shape(dbar_,2)
cccc  #number of grains.
cccCf2__py integer intent(hide), depend(gr_)    :: Ng     = shape(gr_,2)
cccCf2__py integer intent(hide), depend(gr_)    :: Nphase = shape(gr_,0)

      CHARACTER*2 FLABEL
      DIMENSION ROTMAT(3,3)     ! used for rigid sample rotations

ccc   for fortran to python wrapper
      do k=1,nmxprcs
         do i=1,nmxstp
            do j=1,5
c      dbar_=0d0
c      sbar_=0d0
               dbar_(k,i,j)=0.d0
               sbar_(k,i,j)=0.d0
            enddo
            svm_(k,i) = 0.d0
            dvm_(k,i) = 0.d0
            epsvm_(k,i) = 0.d0
         enddo
      enddo

c     modulus file which is written in the subroutine write_stress-strain ----
      open(92, file='modulus_'//jobid,status='unknown')
c     ------------------------------------------------------------------------
      
ccc      

c---------------------------------------------------c
c     OUTPUT FILES are all blocked in order to 
c     import the independency of the wrapped
c     vpsc7.so object
c     All the blocked statement is prepended by 'c##'  --> This is the plan yet to be performed (2011-5-Mar)
c---------------------------------------------------c
c     CALL CPU_TIME (START_TIME)
      i_write_step = 0   !youngung
C ***********************************************************************
C     ASSIGNS # TO UNITS AND OPENS I/O FILES.

      UR0= 0                    ! VPSC7.IN      (OPEN/CLOSE IN MAIN)
      UR1= 1                    ! FILECRYS      (OPEN/READ/CLOSE INSIDE VPSC_INPUT)
      UR2= 2                    ! FILETEXT      (OPEN/READ/CLOSE INSIDE VPSC_INPUT)
      UR3= 3                    ! FILEAXES      (OPEN/READ/CLOSE INSIDE VPSC_INPUT)
      UR4= 4                    ! POSTMORT.IN   (OPEN/CLOSE IN MAIN)
c     UR5= 5     ! CUBCOMP.IN    (OPEN/CLOSE IN MAIN)
c     UR6= 6     ! FILEHIST      (OPEN/CLOSE IN MAIN)
      

C     10     ! RUN_LOG.OUT   (OPEN IN MAIN)
C     11     ! STATS.OUT     (OPEN IN MAIN)
      UW1= 12                   ! RERR.OUT
C     13     ! STR_STR.OUT   (OPEN IN MAIN)
C     14     ! PCYS.OUT      (OPEN IN MAIN)
C     15     ! LANKFORD.OUT  (OPEN IN MAIN)
C     16     ! STAT_AX.OUT   (OPEN IN MAIN)
      UW2= 19                   ! POSTMORT.OUT  (OPEN/CLOSE IN MAIN)
      uw4= 20                   ! for debugging (open and close in the main) -youngung      
C     83     ! FLUCT.OUT     (OPEN IN VPSC_INPUT IF IFLU=1)
C     84     ! FLCUB.OUT     (OPEN IN VPSC_INPUT IF INTERACT=5)
c     92     ! modulus_      (open in main )
c     94     ! crss00000.out - (open and close in write_stress_strain ) - youngung
c     95     ! gr0_000_000 - (open and close in VPSC7.for) -youngung
      uw3= 96                   ! for f2py wrapper out (open and close in the prcs loop) -youngung
C     97              ! SO.OUT        (OPEN IN VPSC_INPUT IF INTERACT=5)
      UR5= 98                   ! CUBCOMP.IN    (OPEN/CLOSE IN MAIN)
      UR6= 99                   ! FILEHIST      (OPEN/CLOSE IN MAIN)

      open(20, file='last_gam.out', status='unknown')
      
      
c##   OPEN(10, FILE='RUN_LOG.OUT',STATUS='UNKNOWN')
      OPEN(11, FILE='STATS'//jobid//'.OUT'  ,STATUS='UNKNOWN')
c##      OPEN(UW1,FILE='RERR.OUT'   ,STATUS='UNKNOWN')
C ***************************************************************************
C     CALL SUBROUTINE VPSC_INPUT FOR READING CRYSTAL, GRAIN AND TEXTURE DATA.
C     INITIALIZES ARRAYS.

      OPEN(UR0,FILE='VPSC7_'//jobid//'.in',STATUS='OLD')

      CALL VPSC_INPUT

      CALL ELSC (0)     ! SC INITIAL ELASTIC MODULI

C ***********************************************************************
C     READS SX & PX STATE FROM PREVIOUS RUN (sav,xmsec,sg,ag,crss,..)
C     WHEN IRECOVER=1.

      IF(IRECOVER.EQ.1) THEN
         OPEN(UR4,file='POSTMORT_'//jobid//'.IN',form='UNFORMATTED',
     $        access='SEQUENTIAL',status='OLD')
         CALL POSTMORTEM (1)    !read the previously saved variables.
         CLOSE(UNIT=UR4)
      ENDIF

C *** SET TO ZERO ACCUMULATED STRAIN COMPONENTS
      EPSACU=0.
      EPSVM =0.
      DO I=1,3
        DO J=1,3
          EPSTOT(I,J)=0.
        ENDDO
      ENDDO

      IF(ICUBCOM.EQ.1) THEN
         OPEN(UR5,file='CUBCOMP_'//jobid//'.IN',status='OLD')
         CALL CUBCOMP(0,0)      ! reads ideal rolling components
         CLOSE(UNIT=UR5)
      ENDIF

C ***********************************************************************
C     OPEN OUTPUT FILES
C ***********************************************************************

      OPEN(13,file='STR_STR'//jobid//'.OUT' ,status='UNKNOWN') ! may be redundant


      IF(ICAUCHY.EQ.1) OPEN(18,file='CAUCHY'//jobid//'.OUT',
     $     status='UNKNOWN')

      IF(NELEM.EQ.1) THEN
        NFILES=NPH
        FLABEL='PH'
      ELSE IF(NELEM.GT.1) THEN
        NFILES=NELEM
        FLABEL='EL'
      ENDIF

      DO I=1,NFILES
        IF(IVGVAR.LE.1) THEN     ! may be redundant
          IUNIT=30+I
          OPEN(IUNIT,
     $         FILE='TEX_'//FLABEL//CHAR(48+I)//'_'//jobid//'.OUT',
     $         STATUS='UNKNOWN')
        ENDIF

        IF(ISHAPE(I).NE.0) THEN
           OPEN(16,FILE='STAT_AX'//'_'//jobid//'.OUT',STATUS='UNKNOWN')
           IUNIT=40+I
           OPEN(IUNIT,FILE='MOR_'//FLABEL//CHAR(48+I)//'_'//
     $          jobid//'.OUT',
     $          STATUS='UNKNOWN')
        ENDIF

        IF(ICUBCOM.EQ.1) THEN
           IUNIT=60+I
           OPEN(IUNIT,FILE='CUBCOMP'//CHAR(48+I)//'_'//jobid//'.OUT',
     $          STATUS='UNKNOWN')
cwx        CALL CUBCOMP (0,1)         ! cubcomp analysis for initial texture
        ENDIF
c
          IUNIT=50+I
          OPEN(IUNIT,
     $         FILE='ACT_'//FLABEL//CHAR(48+I)//'_'//jobid//'.OUT',
     $         STATUS='UNKNOWN')
      ENDDO
c
c     cubcomp analysis for initial texture
c
      IF(ICUBCOM.EQ.1) CALL CUBCOMP (0,1)

C *******************************************************************
C     DO LOOP OVER PROCESSES
C *******************************************************************

C *** READS NUMBER OF PROCESSES TO RUN SEQUENTIALLY. IT CAN BE ANY
C     COMBINATION OF UNIFORM LOAD (IVGVAR=0), VARIABLE LOAD (IVGVAR=1),
C     PCYS (IVGVAR=2) AND LANKFORD (IVGVAR=3).

C *** ONE MORE PROCESS OPTIONS WAS MADE FOR POLYCRYSTAL YIELD SURFACE
C     IN AN APPROXIMATED PLANE STRESS SPACE. (IVGVAR=5)

      READ(UR0,'(A)') PROSA
      READ(UR0,*)     NPROC
      READ(UR0,'(A)') PROSA

      DO 3000 IPROC=1,NPROC
         !open(uw3,file=foutname//char(48+ipro),status='new') !youngung f2py output
C *** READS LOAD CONDITIONS: STRESS, VELOCITY GRADIENT,TEMP, INCR, NSTEPS.
C     CALCULATES SYMMETRIC STRAIN RATE 'DSIM' AND 'DBAR'

      READ(UR0,*) IVGVAR

      IPCYSOPEN=0
      ILANKOPEN=0

      IF(IVGVAR.EQ.0) THEN
        READ(UR0,'(A)') FILEHIST
        OPEN(UR6,file=FILEHIST,status='OLD')
        CALL LOAD_CONDITIONS (UR6)
        CLOSE(UNIT=UR6)
        NSTEPSX=NSTEPS+1      ! EXTRA STEP REQUIRED FOR UPDATING MACRO STRESS
      ELSE IF(IVGVAR.EQ.1) THEN
        READ(UR0,'(A)') FILEHIST
        OPEN(UR6,file=FILEHIST,status='UNKNOWN')
        CALL VAR_VEL_GRAD(0)
        NSTEPSX=NSTEPS
      ELSE IF(IVGVAR.EQ.2) THEN
        if(ipcysopen.eq.0) then
           OPEN(14,file='PCYS'//jobid//'.OUT',status='UNKNOWN')
           ipcysopen=1
        endif
        READ(UR0,*) INDX,INDY
        CALL PCYS (IDUM,INDX,INDY,0)          ! generates PCYS probes
        NSTEPSX=NSTEPS
      ELSE IF(IVGVAR.EQ.3) THEN
        if(ilankopen.eq.0) then
           OPEN(15,file='LANKFORD'//jobid//'.OUT',status='UNKNOWN')
           ilankopen=1
        endif
        READ(UR0,*) DELTALANK
        CALL LANKFORD(ISTEP,DELTALANK,0)     ! initializes arrays
        NSTEPSX=NSTEPS
      ELSE IF(IVGVAR.EQ.4) THEN
         READ(UR0,'(A)') FILEHIST
         OPEN(UR6,file=FILEHIST,status='OLD')
         READ(UR6,'(A)') PROSA
         READ(UR6,*) ((ROTMAT(I,J),J=1,3),I=1,3)
         CLOSE(UNIT=UR6)
         CALL TEXTURE_ROTATION(ROTMAT)
         CALL WRITE_TEXTURE
         
         GO TO 3000

C ***   ADDITIONAL PROCESS FOR PCYS in an approximated PLANE STRESS SPACE
C ***   IT IS ACTUALLY IN A PLANE STRAIN SPACE.
      ELSE IF(IVGVAR.EQ.5) THEN
         if(ipcysopen.eq.0) then
            OPEN(14,file='PCYS_OLD'//jobid//'.OUT',status='UNKNOWN')
            open(141,file='pcys_raw'//jobid//'.out',status='unknown')
            open(142,file='pcys_pl'//jobid//'.out' ,status='unknown')
            ipcysopen=1
        endif
        READ(UR0,*) INDX, INDY, IDUM, DDD
        CALL PCYS_PL(IDUM, INDX, INDY, 0, DDD)
        NSTEPSX=NSTEPS

C ***   additional process for pcys_pl_6d
C ***  
      else if (ivgvar.eq.6) then
         if(ipcysopen.eq.0) then
            open(14, file='PCYS'//jobid//'.OUT', status='unknown')
            open(142,file='pcys_pl_6d'//jobid//'.out', status='unknown')
            ipcysopen=1
         end if
         read(ur0, *) indx, indy, idum, ddd
         call pcys_pl_6d(idum, indx, indy, 0, ddd)
         nstepsx = nsteps
      endif

C *******************************************************************
C     DO LOOP OVER DEFORMATION STEPS
C *******************************************************************

      IF(IVGVAR.LE.1 .AND. NSTEPS.LT.NWRITE) THEN
        WRITE(*,'(/,'' WARNING *** NWRITE='',I3,''  > NSTEPS='',I3,
     $        '' --> WILL NOT WRITE TEXTURE !!!'',/)') NWRITE,NSTEPS
        PAUSE
      ENDIF

      DO 2500 ISTEP=1,NSTEPSX
cw
      if(iflu.eq.1) then
      write(83,'(a,i7)') ' STEP =',istep
      if(icubcom.eq.1) write(84,'(a,i7)') ' STEP =',istep
      endif
cw
      if(interaction.eq.5) then
      write(97,'(a,i7)') ' STEP =',istep
      write(97,'(a)') '    IRS  ITSO  ERRASO      ERRESO'
      endif

cw
c##      WRITE(UW1,'(/,''*******   STEP'',I4)  ') ISTEP
c      WRITE( 10,'(/,''*******   STEP'',I4)  ') ISTEP
c      WRITE(  *,'(/,''*******   STEP'',I4,/)') ISTEP
C      WRITE(*,'(''   ITSGR    SIGGR      SIGAV        DAV    ITTAN'',
c     #               ''     MTAN    ITSEC     MSEC'',/)')

C *** IMPOSE VELOCITY GRADIENT AT EACH STEP:
C       IF IMPOSING MONOTONIC LOADING                 (IVGVAR=0)
C       IF IMPOSING A HISTORY OF DEFORMATION          (IVGVAR=1)
C       IF PROBING THE POLYCRYSTAL YIELD SURFACE      (IVGVAR=2)
C       IF PROBING FOR LANKFORD COEFFICIENTS          (IVGVAR=3)

        IF(IVGVAR.EQ.1) CALL VAR_VEL_GRAD(1)
        IF(IVGVAR.EQ.2) CALL PCYS(ISTEP,INDX,INDY,1)
        IF(IVGVAR.EQ.3) CALL LANKFORD(ISTEP,DELTALANK,1)
C *** POLYCRYSTAL YIELD SURFACE IN AN APPROXIMATED PLANE STRESS SPACE
        IF(IVGVAR.EQ.5) CALL PCYS_PL(ISTEP, INDX,INDY,1, 0.)
        if(ivgvar.eq.6) call pcys_pl_6d(istep,indx,indy,1, 0.)

C *******************************************************************

        CALL UPDATE_SCHMID

C *******************************************************************
C     FOR VOCE NORMALIZE ALL ARRAYS WITH UNITS OF STRESS, STRAIN RATE &
C     COMPLIANCE OR SCALES GAMD0 AND CRSS TO WORK WITH n POWERS OF ORDER ONE.
C     IRATESENS=0 ALLOWS TO ELIMINATE RATE SENSITIVITY.

        IF(IHARDLAW.NE.1) THEN
          IF(IRATESENS.EQ.0) CALL SCALE_3 (ISTEP)   ! rate-insensitive response
C         IF(IRATESENS.EQ.1) CALL NORMALIZE (0,ISTEP)
        ENDIF

C     FOR MTS SCALE GAMD0 TO WORK WITH n POWERS OF ORDER ONE, AND
C     CALCULATE CRRS's ASSOCIATED WITH (POSSIBLE) NEW TEMP AND RATE.

        IF(IHARDLAW.EQ.1) THEN
          CALL SCALE_3 (ISTEP)
          EDOT = SQRT(2./3.)*TNORM(DBAR,5,1)
          CALL UPDATE_CRSS_MTS(EDOT,TEMPERAT,1)
        ENDIF
C *******************************************************************
C     IF IMPOSING STRAIN MAKES A TAYLOR GUESS FOR THE FIRST STEP WHEN
C     IRECOVER=0 OR STARTS FROM RECOVERED STRESS STATE 'SG' WHEN IRECOVER=1.
C     IF IMPOSING STRESS MAKES A SACHS GUESS FOR EVERY STEP.

        IF(ISTEP.EQ.1.AND.IRECOVER.EQ.0 .OR. STRAIN_CONTROL.EQ.0) THEN
          CALL INITIAL_STATE_GUESS
        ENDIF
        IF(ISTEP.EQ.1.AND.IRECOVER.EQ.1) THEN
          KGX=1
          DO IPH=IPHBOT,IPHTOP
            IPHEL=IPH-IPHBOT+1
          DO KKK=NGR(IPH-1)+1,NGR(IPH)
            CALL GRAIN_RATE_AND_MODULI (1,1,KGX,KKK,IPHEL,IPH)
            KGX=KGX+1
          ENDDO
          ENDDO
        ENDIF

c     youngung - store the previous gamdot to check the loading direction change
        if(istep.eq.1.and.irecover.eq.0) then
           do iph = iphbot, iphtop
              do kkk = ngr(iph-1), ngr(iph)
                 do is=1, nsyst(iph)
                    prgamdot(is,iph) = 0.
                 enddo          ! loop over system
              enddo             ! loop over gain
           enddo                ! loop over phase
        else
           do iph = iphbot, iphtop
              do kkk = ngr(iph-1)+1, ngr(iph)
                 do is=1, nsyst(iph)
                    prgamdot(is,kkk)= gamdot(is,kkk)
                 enddo          ! loop over system
              enddo             ! loop over gain
           enddo                ! loop over phase           
        endif
        write(20,'(a,i3)') 'step:', istep
        do iph = iphbot, iphtop
           do kkk= ngr(iph-1)+1, ngr(iph)
              do is=1, nsyst(iph)
                 write(20,'(2e11.3)') prgamdot(is, kkk)
              enddo
           enddo
        enddo
        
c     ------------------------------------------------------------------------        

C *************************************************************************
C     VPSC CALCULATION (OR TAYLOR CALCULATION WHEN INTERACTION=0) OF
C     STRESS AND STRAIN-RATE FOR EVERY GRAIN.

        CALL VPSC (ISTEP)      !!!     THIS IS THE CORE OF THE CODE

C ***************************************************************************
C     VPSC PROVIDES THE MACROSCOPIC STRAIN RATE 'DSIM'.
C     CALCULATES MACROSCOPIC ROTATION RATE AND MACROSCOPIC VELOCITY GRADIENT
C     DEPENDING ON THE BOUNDARY CONDITIONS.

        
      DO I=1,3
        ROTBAR(I,I)=0.
        DO J=1,3
          IF(IUDOT(I,J).EQ.1.AND.IUDOT(J,I).EQ.1) THEN
            ROTBAR(I,J)=(UDOT(I,J)-UDOT(J,I))/2.
          ELSE IF(IUDOT(I,J).EQ.1) THEN
            ROTBAR(I,J)=UDOT(I,J)-DSIM(I,J)
            ROTBAR(J,I)=-ROTBAR(I,J)
          ENDIF
          UDOT(I,J)=DSIM(I,J)+ROTBAR(I,J)
        ENDDO
      ENDDO


      
      

      

C     VON MISES STRAIN-RATE & STRESS.
C     IF INTERACT=0: DBAR IS IMPOSED AND WE DEFINE SBAR=SAV.
C     IF INTERACT>0: DBAR & SBAR FOLLOW FROM THE CALL TO SUBR. STATE6x6.

      SVM=0.
      DVM=0.
      DO I=1,5
        SVM=SVM+SBAR(I)*SBAR(I)
        DVM=DVM+SBAR(I)*DBAR(I)
      ENDDO
      SVM=SQRT(SVM*3./2.)
      DVM=DVM/SVM

C     F2PY BLOCK WRITES -youngung
      do i=1,5
         dbar_(iproc,istep,i) = DBAR(I)
         sbar_(iproc,istep,i) = SBAR(I)
      enddo
      svm_(iproc,istep) = svm
      dvm_(iproc,istep) = dvm
      epsvm_(iproc,istep) = epsvm






      
C     F2PY

cc    F2PY prints to screen
c     write(*,'(a,2i3)') ' Current iproc and step : ', iproc, istep
c     write(*,'(a,5f11.5,3x,a,11f9.5)') 'dbar =', dbar, 'sbar =', sbar
c     write(*,'(a,5f11.5,3x,a,11f9.5)') 'dbar_=', (dbar_(iproc,istep,i),
c     $     i=1,5), 'sbar_=', (sbar_(iproc,istep,i),i=1,5)
ccc   block ends


C *************************************************************************
C     RENORMALIZE ALL ARRAYS WITH UNITS OF STRESS, STRAIN RATE & COMPLIANCE

        IF(IHARDLAW.NE.1) THEN
C         IF(IRATESENS.EQ.1) CALL NORMALIZE (1,ISTEP)
        ENDIF

C *************************************************************************
C     WRITES AND/OR STORES STRESS-STRAIN STATE FOR THE STEP.

      i_write_step = i_write_step + 1  !youngung
c      IF(IVGVAR.LE.1) CALL WRITE_STRESS_STRAIN (ISTEP)
      IF(IVGVAR.LE.1) CALL WRITE_STRESS_STRAIN (i_write_step)
      IF(IVGVAR.EQ.2) CALL PCYS (ISTEP,INDX,INDY,2)
      IF(IVGVAR.EQ.3) CALL LANKFORD (ISTEP,DELTALANK,2)
      IF(IVGVAR.EQ.5) CALL PCYS_PL(ISTEP,INDX,INDY,2,0.)
      if(ivgvar.eq.6) call pcys_pl_6d(istep,indx,indy,2,0.)

cwx      IF(ICUBCOM.EQ.1) CALL CUBCOMP (ISTEP,1)
      IF(ICUBCOM.EQ.1.AND.ISTEP.GT.1) CALL CUBCOMP (ISTEP,1)

C *******************************************************************
C     STATISTICS & OUTPUT
C *******************************************************************

      CALL STAT_SHEAR_ACTIVITY
      IF(IHARDLAW.NE.2) CALL STAT_GRAIN_SHAPE
      IF(IHARDLAW.NE.2) CALL STAT_STRESS_STRAIN

      CALL WRITE_STAT(ISTEP)
      CALL WRITE_ACTIV(ISTEP)

cccc f2py -youngung
      do i=1,3
         do j=1,3
            scauchy_(iproc,istep,i,j) = scauchy(i,j)
            sdeviat_(iproc,istep,i,j) = sdeviat(i,j)
            dsim_(iproc,istep,i,j) = dsim(i,j)
         enddo
      enddo
cccc
      IF(IVGVAR.GE.2) GO TO 2500

      IF(ISTEP.EQ.NSTEPSX) GO TO 2000    ! LAST STEP IS ONLY TO UPDATE STRESS

C ***********************************************************************
C     GIVEN THE OVERALL STRAIN RATE AND STRESS THEN:
C     *  IF ICTRL=0 CALCULATES TIME INCREMENT NECESSARY TO ACHIEVE
C        THE IMPOSED VON MISES STRAIN INCREMENT.
C     *  IF ICTRL=1,6 CALCULATES TIME INCREMENT NECESSARY TO ACHIEVE
C        THE IMPOSED STRAIN COMPONENT INCREMENT.
C     *  IF ICTRL=7 AND STRAIN_CONTROL=1 A TIME INCREMENT IS IMPOSED.
C        --> CALCULATES THE VON MISES STRAIN INCR ASSOCIATED WITH THE
C            IMPOSED STRAIN-RATE.
C     *  IF STRAIN_CONTROL=0 IT CORRESPONDS TO A CREEP SIMULATION
C        --> CALCULATES THE STRAIN INCREMENT ASSOCIATED WITH THE IMPOSED
C            STRESS COMPONENT AND TIME INCREMENT
C ***********************************************************************

      IF(ICTRL.EQ.0) THEN
        TINCR  =EVMINCR/DVM
        EPSACU =EPSACU+EVMINCR
      ELSE IF(ICTRL.GE.1.AND.ICTRL.LE.6 .AND. STRAIN_CONTROL.EQ.1) THEN
        TINCR  =EIJINCR/ABS(DSIMCTRL)
        EPSACU =EPSACU+EIJINCR
      ELSE IF(ICTRL.EQ.7 .OR. STRAIN_CONTROL.EQ.0) THEN
        EVMINCR=DVM*TINCR
        EPSACU =EPSACU+EVMINCR
      ENDIF

      EPSVM=EPSVM+DVM*TINCR
      DO I=1,3
      DO J=1,3
         EPSTOT(I,J) = EPSTOT(I,J) + DSIM(I,J) * TINCR
      ENDDO
      ENDDO

ccc f2py -Youngung
      do i=1,3
         do j=1,3
            eps_tot(iproc,istep,i,j) = epstot(i,j)
         enddo
      enddo
ccc

C ***********************************************************************
C     UPDATES ORIENTATION, HARDENING AND SHAPE OF EVERY GRAIN MAKING
C     A LINEAR EXTRAPOLATION OF THE CALCULATED RATE TIMES 'TINCR'.
C     (SKIPS UPDATE FOR PCYS OR LANKFORD RUN)
C ***********************************************************************

c        WRITE(10,'(/,'' *** UPDATING AT THE END OF STEP'',i3)') ISTEP

        CALL UPDATE_FIJ(0)        ! UPDATES DEFORM TENSOR OF ELEMENT
        CALL UPDATE_SHAPE(0)      ! UPDATES SHAPE OF ELEMENT

        IF(IHARDLAW.LE.1) THEN
           DO IPH=IPHBOT,IPHTOP
              IPHEL=IPH-IPHBOT+1
              IF(NTWMOD(IPHEL).NE.0) CALL UPDATE_TWINNING (IPH)
           ENDDO
           
           IF(IUPDORI.EQ.1) CALL UPDATE_ORIENTATION
           
           DO IPH=IPHBOT,IPHTOP
              CALL UPDATE_FIJ(IPH) ! UPDATES DEF TENSOR OF PHASE & GRAINS
              IF(IUPDSHP.EQ.1) CALL UPDATE_SHAPE (IPH)
           ENDDO
           
           IF(IUPDHAR.EQ.1) THEN
C     *** VOCE HARDENING PLUS PREDOMINANT TWIN REORIENTATION SCHEME
              IF(IHARDLAW.EQ.0)  CALL UPDATE_CRSS_VOCE (2)
C     *** HARDENING FOR IRRADIATED MATERIAL
              IF(IHARDLAW.EQ.-1) CALL UPDATE_CRSS_IRRAD(2)
c     *** hardening by dislocation density evolution
              if(ihardlaw.eq.-2) call update_crss_dislocation(2,istep)
C     *** MECHANICAL THRESHOLD STRESS HARDENING (NO TWINNING, ONLY ONE PHASE)
              IF(IHARDLAW.EQ.1) THEN
                 TEMPX=TEMPERAT
                 EDOT=SQRT(2./3.)*TNORM(DBAR,5,1)
                 CALL UPDATE_CRSS_MTS (EDOT,TEMPX,2)
                 TEMPERAT=TEMPX
              ENDIF
           ENDIF
        ENDIF


c   ***  STRAIN_INDUCED MARTENSITIC TRANSFORMATION IN META-STABLE AUSTENITIC STAINLESS STEEL
c        write(*,*) 'Check itran'
c        write(*,*)itran
c        pause
        IF (ITRAN.eq.1) then
           if (nph.ge.2) then
              write(*,*) ' you are calling phase_trans subr'
              call PHASE_TRANS(ngr(1),istep)
              
           else if (nph.lt.2) then
              write(*,'(a)')"Number of phase must exceed 1"
              write(*,'(a)')"If willing to phase transform"
           endif
        endif



ccc   f2py - Youngung writes texture into gr files
ccc      gr file making
ccc    filename = gr_0_000_000
ccc          1st index: number of ph
ccc          2nd index: number of proc
ccc          3rd index: number of steps
        
        do iph = iphbot, iphtop
           !proc 
           im1 = proc/1000
           im2 = (proc - im1*1000) / 100
           im3 = (proc - im1*1000 - im2*100)/10
           im4 = (proc - im1*1000 - im2*100 - im3*10)
           
           !istep
           in1 = istep/1000
           in2 = (istep - in1*1000)/100
           in3 = (istep - in1*1000 - in2*100)/10
           in4 = (istep - in1*1000 - in2*100 - in3*10)

           !ijob
           
           open(95, file = 'gr_'//jobid//'_'
     $          //char(48+iph)//'_'

     $          //char(48+im1)
     $          //char(48+im2)
     $          //char(48+im3)
     $          //char(48+im4)//'_'
           
     $          //char(48+in1)
     $          //char(48+in2)
     $          //char(48+in3)
     $          //char(48+in4), 
     $          status ='unknown')
           do kkk = ngr(iph-1)+1, ngr(iph)
              do i = 1,3
                 do j = 1,3
                    aux(I,J) = AG(J,I,kkk)
                 enddo
              enddo
              call euler(1, GDA, GDB, GDC, aux)
              write(95, '(3f10.3, f12.7)') GDA,GDB,GDC,wgt(kkk)
c$$$              gr_(iph, iproc, istep, kkk, 1) = GDA
c$$$              gr_(iph, iproc, istep, kkk, 2) = GDB
c$$$              gr_(iph, iproc, istep, kkk, 3) = GDC
c$$$              gr_(iph, iproc, istep, kkk, 4) = wgt(kkk)
           enddo !do -loop over grains : 
           close(95)            ! file closed. 
        enddo  ! do-loop over phase
ccc        

 2000 CONTINUE

C *******************************************************************
C *** END OF FORWARD UPDATE OF STRAIN, SHAPE, ORIENTATION & CRSS
C *******************************************************************

C *** CALCULATES evm/svm/taylor/work FOR EACH GRAIN USING sg/dg/dbar

      CALL GRAIN_INFO

C *******************************************************************
C     WRITE TEXTURE FILES FOR EACH PHASE
C *******************************************************************

      IF(IVGVAR.LE.1) THEN
        IWRITE=0
        IF(NWRITE.EQ.0) THEN
          IF(ISTEP.EQ.NSTEPS) IWRITE=1
        ELSE IF(NWRITE.NE.0) THEN
          IF(MOD(ISTEP,NWRITE).EQ.0) IWRITE=1
        ENDIF
        IF(IWRITE.EQ.1) CALL WRITE_TEXTURE
      ENDIF

C     WRITES GRAIN'S & PX STATES INTO BINARY FILE 'POSTMORT.OUT'

      IF(ISAVE.EQ.ISTEP) THEN
         OPEN(UW2,file='POSTMORT_'//jobid//'.OUT',form='UNFORMATTED',
     $        access='SEQUENTIAL',status='UNKNOWN')
c     write variables into file(uw2)
         CALL POSTMORTEM (2)   
         CLOSE(UNIT=UW2)
         ISAVE=-1               !to avoid overwrite
      ENDIF

      CALL ELSC (1)     ! SC CALCULATION OF ELASTIC MODULI AT EACH STEP

 2500 CONTINUE      ! END OF DO LOOP OVER DEFORMATION INCREMENTS

 3000 CONTINUE      ! END OF DO LOOP OVER PROCESSES

C *******************************************************************

c$$$      CALL CPU_TIME (END_TIME)
c$$$      DELTIME=END_TIME-START_TIME
c$$$      WRITE(*,*)
c$$$      WRITE(*,'('' TOTAL RUN TIME:'',F10.2,'' secs'')') DELTIME
      !STOP
      close(92)
      close(20) !youngung jeong , debugging 
      return 
      END

C*************************************************************************
C     INCLUDE HERE THE FILE 'VPSC7.SUB' WHICH CONTAINS THE SUBROUTINES
C     FOR PERFORMING THE SELF-CONSISTENT CALCULATION AND THE REORIENTATION
C
C     INCLUDE ALSO THE LIBRARY OF CANNED SUBROUTINES 'LIBRARY.SUB'
C
C*************************************************************************

      INCLUDE 'vpsc7.sub'
      INCLUDE 'library7.sub'



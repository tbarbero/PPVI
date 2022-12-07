!! Compiling on analysis:
!    module load netcdf/4.6.1
!    module load intel_compilers/oneapi
!    # Full speed run:
!    ifort -O2 sqinvph.F90 -o sqinvph.exe -lnetcdf -lnetcdff -lhdf5 -lhdf5_hl -lz -limf -lm
!!! Modified by Jan-Huey Chen in May 2021

      PROGRAM QINVERTP25S4
!c  ** modified from qinvertp25.for, adding horizontaL surgery, only leave the PV core  9/1/92  CC WU
!c  ** modified so can surgery on the surface thta too.   cc wu  12/9/92
!c  ** modified so have more choice on nonlinear partition   cc wu  12/9/92

!C************************************************
!C****** Perturbation PV inversion program *******
!C************************************************
!C  Most of the variables used here have meanings similar to those in the total PV inversion program. At the point the variables
!C  are passed to the subroutine, those endine in 'B' are mean quantities (QB for PV , MB for geopotential, SB for stream
!C  function). Those ending in 'P' are perturbation quantities. However, these variables are changed from the way they are read
!C  in, according to some of the options discussed below. There are 4 input files, q and h for the mean and q and h
!C  for either the total fields. Perturbation fields are computed from these.
!C*********************************** OPTIONS ********************************
!C  OMEGAS and OMEGAH are overrelaxation parameters for S and H in subroutine BALP. PRT is the underrelaxation
!C  parameter. These can all have the same values as in tthe total inversion program. THRSH is the dimensional convergence
!C  threshold (meters) and TSCAL and QSCAL are scale factors to multiply the input boundary theta and PV respectively
!C  (generally = 1.0).
!C        
!C  INLIN = option for nonlinear terms. If INLIN = 1, there are no terms dropped in the perturbation inversion equations. 
!C  The equations are still linear, but the nonlinear terms are hidden in the coefficents of the differential operator. The
!C  "mean variables" are redefined to include the actual mean field plus 1/2 the perturbation field. If INLIN=0, nonlinear
!C  terms are dropped altogether. Usually, INLIN =1, otherwise, the sum of piecewise solutions will not equal the total 
!C  perturbation.
!C
!C  IQD = option to make value of q' conditional on some other field. Suppose you wanted the PV perturbation to be only
!C  the perturbation PV in saturated air, you could read in a relative humidity file here and with CRIT, specify the
!C  threshold value of relative humidity below which q'=0. There is also a hardwired option to allow only positive q' in saturated
!C  air: this can be changed.
!C
!C  SUB BALP: These options should be pretty self-explanatory, beginning with the number and list of the output perturbation h and s
!C  levels, the number of inversions to be done and then the number and list of pert. PV levels to be included.
!C  1=lower boundary theta 2 to nl-1 are interior pert. PV levels, bottome to top
!C  NL=upper boundary perturbation theta.
!C
!C  If you don't choose either boundary theta field, homogeneous Neumann conditions are applied at the top and bottom on both
!C  'h' and 's'.
!C
!C  IBC = option for lateral boundary conditions.
!C        
!C  IBC=0 => homogeneous Diriclet conditions.
!C  IBC=1 => 's' and 'h' at the boundaries are equal to the full perturbation
!C  IBC=2 => you will read in a file with an estimate of the interior and boundary values. This option  is useful for "nesting" 
!C           the inversion, i.e., using the boundary conditions from a calculation done on a larger grid.
!C****************************************************************************************
      implicit none
#include <netcdf.inc>

      character*190 :: nmlfile = 'sqinv.nml'
      integer, external ::  iargc
      character*300 :: infile1, infile2, outfile
      character*300 :: maskfiles(40)
      character*90 :: time_unit
      integer :: ntimes(40) = 999
      real, dimension (:,:,:), allocatable :: QP, QB, MP, MB, SIP, SB
      real, dimension (:,:,:), allocatable :: THB, THP, DPF, QPO
      real, dimension (:,:), allocatable :: FM, FC, A, MF, LATT
      real, dimension (:), allocatable :: AP, APP, APM, PI, PR, PIF

      real :: ZHDR(8), HND, SIGM, AA, BETA, CNLINa, CNLINb
      real :: KAP, R, CV, FF, LL, PII, FRC, THO, QCONST, CRIT
      integer :: HDR(8), BRES, NDAYS, LPD
      character*120 ::  DPFIL
      integer :: jsurg1, jsurg2
      integer, dimension (:), allocatable :: IQP
      integer, dimension (:,:), allocatable :: isurg1, isurg2
      integer, parameter :: MAX = 1000, MAXT = 500
      real, parameter :: QMIN = 0.01, DPI=50., CP=1004.5, PO=1.E5, MI=9999.90, GG=9.8066

      real :: OMEGAS, OMEGAH, PRT, THRSH, TSCAL, QSCAL
      integer :: IMAP, INLIN, IQD, NOUT, IBC
      integer :: itsurg = 0

      logical :: mask = .false.
      namelist /nlist/ outfile, infile1, infile2, maskfiles, ntimes, &
                       OMEGAS, OMEGAH, PRT, THRSH, TSCAL, QSCAL, &
                       IMAP, INLIN, IQD, NOUT, itsurg, IBC, mask

      integer :: fid1, fid2, status
      integer :: nx, ny, nl, nt, nnt
      integer :: i, j, k, t, n, nn, tt
      real, dimension (:), allocatable :: lon_in, lat_in, plev_in, time_in
      real, dimension (:,:,:,:), allocatable :: qb_in, psib_in, hb_in
      real, dimension (:,:,:,:), allocatable :: qp_in, psip_in, hp_in
      real, dimension (:,:,:), allocatable :: qp_out3, hp_out3, sp_out3, ur_out3, vr_out3
      real, dimension (:,:,:,:), allocatable :: qp_out, hp_out, sp_out, ur_out, vr_out
      integer :: id, londimid, latdimid, levdimid, timedimid
      integer :: lonid, latid, levid, timeid
      integer :: qpid, hpid, spid, urid, vrid
      integer :: dims(4),start(4), total(4)
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Read in namelist and set up range:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do n=1, 40
         maskfiles(n) = 'ttt'
      enddo

      if (iargc() > 0 ) then
        call getarg(1,nmlfile)
      endif

      open (11, file=trim(nmlfile), status='old')
      read (11, nml=nlist)
      close (11)

      if (mask) then
        nnt = 0
        do n = 1, 40
           if ( ntimes(n)/= 999) nnt = nnt + 1
           print *, 'JHC test, ntimes =', ntimes(n)
           print *, 'JHC test, maskfiles =',maskfiles(n)
        enddo
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! READ dim from netcdf file and setup variables/parameters
!!!!!!!!!!!!!!!!!!!!!!!!!! = 0!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call open_ncfile(infile1,fid1)
      print *, 'Reading ', trim(infile1), fid1

      call get_ncdim1(fid1, 'lon', nx)
      call get_ncdim1(fid1, 'lat', ny)
      call get_ncdim1(fid1, 'lev', nl)
      call get_ncdim1(fid1, 'time', nt)
      print *, 'data nlon, nlat, np, nt = ', nx, ny, nl, nt

      allocate(plev_in(nl), lon_in(nx), lat_in(ny), time_in(nt))
      call get_var1_real (fid1, 'lev', nl, plev_in)
      call get_var1_real (fid1, 'lon', nx, lon_in)
      call get_var1_real (fid1, 'lat', ny, lat_in)
      call get_var1_real (fid1, 'time', nt, time_in)
      call get_att_text (fid1, 'time', 'units', time_unit)
  
      allocate(pr(nl))
      do k = 1, nl
        pr(k) = plev_in(k)/1000.
      enddo
      
      allocate(QP(ny,nx,nl), QB(ny,nx,nl), MP(ny,nx,nl), MB(ny,nx,nl))
      allocate(SIP(ny,nx,nl), SB(ny,nx,nl), DPF(ny,nx,nl), QPO(ny,nx,nl))
      allocate(THB(ny,nx,2), THP(ny,nx,2))
      allocate(FM(ny,nx), FC(ny,nx), MF(ny,nx), LATT(ny,nx))
      allocate(A(ny,5))
      allocate(AP(ny), APP(ny), APM(ny))
      allocate(PI(nl), PIF(nl))
      allocate(IQP(nl))
      allocate(isurg1(100,nl),isurg2(100,nl))

      zhdr(1) = lat_in(1)
      zhdr(2) = lon_in(1)
      zhdr(3) = lat_in(ny)
      zhdr(4) = lon_in(nx)
      zhdr(5) = abs(lon_in(2) - lon_in(1))
      zhdr(6) = abs(lat_in(2) - lat_in(1))

      zhdr(7) = nx
      zhdr(8) = ny

      do i=1,8
        HDR(i)=NINT(ZHDR(i))
      enddo

      PII=4.*ATAN(1.)
      AA=2.E7/PII
      R=2.*CP/7.
      CV=CP-R
      FF=1.E-4
      SIGM=ZHDR(5)/ZHDR(6)
      KAP=R/CP

      ! 2 FOR XY, 1 FOR SPHERICAL
      IF (IMAP.EQ.1) THEN
         LL=AA*ZHDR(5)*PII/180.
         do I=1,HDR(8)
           AP(I)=COS( PII*(ZHDR(3) - (I-1)*ZHDR(6))/180. )
           APM(I)=COS( PII*(ZHDR(3) - (I-1.5)*ZHDR(6))/180. )
           APP(I)=COS( PII*(ZHDR(3) - (I-0.5)*ZHDR(6))/180. )
           A(I,1)=SIGM*SIGM*APM(I)/AP(I)
           A(I,2)=1./( AP(I)*AP(I) )
           A(I,3)=-( 2. + SIGM*SIGM*AP(I)*(APM(I)+APP(I)) )/( AP(I)*AP(I) )
           A(I,4)=1./( AP(I)*AP(I) )
           A(I,5)=SIGM*SIGM*APP(I)/AP(I)
           do J=1,HDR(7)
            FC(I,J)=1.458*SIN( PII*(ZHDR(3) - (I-1)*ZHDR(6))/180. )
            MF(I,J)=1.
            FM(I,J)=FC(I,J)/MF(I,J)
           enddo
         enddo
      ELSE
         LL=ZHDR(5)*1.E3
         OPEN(23,FILE='latlon',STATUS='OLD')
         READ(23,*)(ZHDR(I),I=1,8)
         DO 105 I=1,HDR(8)
           READ(23,*)(MF(I,J),J=1,HDR(7))
 105     CONTINUE        
         DO 103 I=1,HDR(8)
           READ(23,*,END=103)(LATT(I,J),J=1,HDR(7))
 103     CONTINUE        
         DO 102 I=1,HDR(8)
           AP(I)=1.
           A(I,1)=SIGM*SIGM
           A(I,2)=1.
           A(I,3)=-2*( 1. + SIGM*SIGM )
           A(I,4)=1.
           A(I,5)=SIGM*SIGM
           DO 106 J=1,HDR(7)
             FC(I,J)=1.458*SIN( PII*LATT(I,J)/180. )
             MF(I,J)=MF(I,J)*MF(I,J)
             FM(I,J)=FC(I,J)/MF(I,J)
 106       CONTINUE
 102     CONTINUE
      END IF
 
      THO=FF*FF*LL*LL/DPI
      FRC=DPI*THO/(FF*FF*LL*LL)
      QCONST=1.E6*KAP*GG*CP*FF*THO/(DPI*PO)
      BETA=1.458*LL/AA
      HND=DPI*THO/GG
      THRSH=THRSH/HND
      !WRITE(6,*)FRC,QCONST,THO,HND,THRSH
      
      do K=1,NL
       PI(K)=CP*( PR(K)**KAP )/DPI
       PIF(K)=(DPI*PI(K)/CP)**2.5
      enddo

      ! read in data
      allocate(qb_in(nx,ny,nl,nt), psib_in(nx,ny,nl,nt), hb_in(nx,ny,nl,nt))
      allocate(qp_in(nx,ny,nl,nt), psip_in(nx,ny,nl,nt), hp_in(nx,ny,nl,nt))
      call get_var4_real (fid1, 'pv', nx, ny, nl, nt, qb_in)
      call get_var4_real (fid1, 'psi', nx, ny, nl, nt, psib_in)
      call get_var4_real (fid1, 'h', nx, ny, nl, nt, hb_in)

      call open_ncfile(infile2,fid2)
      call get_var4_real (fid2, 'pv', nx, ny, nl, nt, qp_in)
      call get_var4_real (fid2, 'psi', nx, ny, nl, nt, psip_in)
      call get_var4_real (fid2, 'h', nx, ny, nl, nt, hp_in)


      allocate(qp_out3(nx,ny,nl), hp_out3(nx,ny,nl), sp_out3(nx,ny,nl))
      allocate(ur_out3(nx,ny,nl), vr_out3(nx,ny,nl))
      allocate(qp_out(nx,ny,nl,nt), hp_out(nx,ny,nl,nt), sp_out(nx,ny,nl,nt))
      allocate(ur_out(nx,ny,nl,nt), vr_out(nx,ny,nl,nt))
      ! Start to loop with time record

      ! JHC
      !do t = 1, 1
      do t = 1, nt

         do k = 1, nl
            do j = ny,1,-1
               do i = 1, nx
                  qb(j,i,k) = qb_in(i,ny-j+1,k,t)
                  sb(j,i,k) = psib_in(i,ny-j+1,k,t)
                  mb(j,i,k) = hb_in(i,ny-j+1,k,t)
                  qp(j,i,k) = qp_in(i,ny-j+1,k,t)
                  sip(j,i,k) = psip_in(i,ny-j+1,k,t)
                  mp(j,i,k) = hp_in(i,ny-j+1,k,t)
               enddo
             enddo
         enddo

         DO I=1,HDR(8)
           DO J=1,HDR(7)
             THB(I,J,1)=QB(I,J,1)
             THB(I,J,2)=QB(I,J,NL)

             THP(I,J,1)=QP(I,J,1)
             THP(I,J,2)=QP(I,J,NL)

             DO K=1,2
               THB(I,J,K)=THB(I,J,K)/THO
               THP(I,J,K)=TSCAL*THP(I,J,K)/THO - THB(I,J,K)
             ENDDO
           ENDDO
         ENDDO

         !C********* Nondimensionalize variables ******************************
         !c  **  add something new for more choice about partitioning cc wu 12/9/92
         ! LINEARIZATION PARM (0=THROW OUT NONLINEAR TERMS, 1=EQUAL PARTITIONING [FULL LINEAR])
         IF (INLIN.EQ.1) THEN 
           CNLINa=0.5
           CNLINb=0.5
         ELSE IF (INLIN.EQ.2) THEN 
           CNLINa=1.
           CNLINb=0.
         ELSE IF (INLIN.EQ.3) THEN 
           CNLINa=0.
           CNLINb=1.
         ELSE
           CNLINa=0.
           CNLINb=0.
         END IF
         !print *, 'CNLINa= ',cnlina,'CNLINb= ',cnlinb

         do K=1,NL
            do J=1,HDR(7)
            do I=1,HDR(8)
               MP(I,J,K)=MP(I,J,K)/HND
               SIP(I,J,K)=SIP(I,J,K)/HND
               MB(I,J,K)=MB(I,J,K)/HND
               SB(I,J,K)=SB(I,J,K)/HND
               MP(I,J,K)=MP(I,J,K)-MB(I,J,K)
               SIP(I,J,K)=SIP(I,J,K)-SB(I,J,K)
               MB(I,J,K)=MB(I,J,K) + CNLINa*MP(I,J,K)
               SB(I,J,K)=SB(I,J,K) + CNLINb*SIP(I,J,K)
            enddo
            enddo
         enddo
 
         DO 280 K=2,NL-1
            DO 281 J=1,HDR(7)
            DO 282 I=1,HDR(8)
              IF (QP(I,J,K).EQ.MI) THEN
                 QP(I,J,K)=0.
                 GO TO 282
              END IF
              IF (QP(I,J,K).LT.1.E2*QMIN) THEN
                 QP(I,J,K)=PIF(K)*QMIN/QCONST
              ELSE
                 QP(I,J,K)=PIF(K)*QP(I,J,K)/(1.E2*QCONST)
              END IF
              IF (QB(I,J,K).LT.1.E2*QMIN) THEN
                 QB(I,J,K)=PIF(K)*QMIN/QCONST
              ELSE
                 QB(I,J,K)=PIF(K)*QB(I,J,K)/(1.E2*QCONST)
              END IF
              QP(I,J,K)=(QP(I,J,K)-QB(I,J,K))/MF(I,J)
  282       CONTINUE
  281       CONTINUE
  280    CONTINUE

         ! IQD: ENTER "1" IF QPRIME.NE.0 DEPENDS ON ANOTHER FIELD
         IF (IQD.EQ.1) THEN ! JHC: won't come in this part
            PRINT*,'ENTER THE FILE, THE DAY AND NUM LEVELS PER DAY.'
            READ(5,*)DPFIL
            READ(5,*)NDAYS
            READ(5,*)LPD
            OPEN(18,FILE=DPFIL,STATUS='OLD')
            READ(18,*)(ZHDR(I),I=1,8)
            DO 290 NN=1,NDAYS
              DO 291 K=1,LPD
                DO 292 I=1,HDR(8)
                  READ(18,*,END=297)(DPF(I,J,K),J=1,HDR(7))
  292           CONTINUE
  291         CONTINUE
  290       CONTINUE
            OPEN(19,FILE='DT:[CDAVIS]QPR.DAT',STATUS='NEW')
  297       WRITE(19,*)(ZHDR(I),I=1,8)
            DO 300 K=2,NL-1
              PRINT*,'ENTER VALUE OF THIS FIELD BELOW WHICH QPRM=0.'
              READ*,CRIT
              IQP(K)=0
              DO 301 J=2,HDR(7)-1
              DO 302 I=2,HDR(8)-1          
                IF ((DPF(I,J,K).LT.CRIT).OR.(QP(I,J,K).LT.0.)) THEN
                   QP(I,J,K)=0.
                ELSE
                   IQP(K)=IQP(K) + 1
                END IF
 302          CONTINUE
 301          CONTINUE
              IF (IQP(K).GT.0) THEN
                WRITE(6,*)IQP(K),K
                DO 310 I=1,HDR(8)
                DO 312 J=1,HDR(7)
                  QPO(I,J,K)=QP(I,J,K)*1.E2*QCONST/PIF(K)
 312            CONTINUE
                WRITE(19,399)(QPO(I,J,K),J=1,HDR(7))
 310            CONTINUE
              END IF
 300        CONTINUE
 399        FORMAT(10F8.2)
         END IF    

         !C*********** Routine to solve for balanced pert flow by S.O.R. **************

         if (mask) then 
           itsurg = 0
           do tt = 1, nnt
              if (t == ntimes(tt)) then 
                 itsurg = 3
                 exit
              endif
           enddo
         endif

         print *, 'ntime = ', t, 'itsurg = ', itsurg
         print *, 'maskfile = ', maskfiles(t)

         CALL BALP(t, nx, ny, nl, MP, MB, SIP, SB, THP, QP, FC, FM, MF, AP, A, PI, HND, \
                 OMEGAS,OMEGAH,PRT,THRSH,BETA,FRC,MAX,MAXT,SIGM,ZHDR ,pif,tho,qconst,\
                 MI, NOUT, IBC, itsurg, maskfiles(t),                                \
                 qp_out3, hp_out3, sp_out3, ur_out3, vr_out3)

         qp_out(:,:,:,t) = qp_out3(:,:,:)
         hp_out(:,:,:,t) = hp_out3(:,:,:)
         sp_out(:,:,:,t) = sp_out3(:,:,:)
         ur_out(:,:,:,t) = ur_out3(:,:,:)
         vr_out(:,:,:,t) = vr_out3(:,:,:)

      enddo ! t loop

      ! write out netcdf file
      print *, 'outfile = ', trim(outfile)
      status = nf_create(outfile, NF_NOCLOBBER, id)

      status = nf_def_dim(id, 'lon',  nx, londimid)
      status = nf_def_dim(id, 'lat',  ny, latdimid)
      status = nf_def_dim(id, 'lev',  nl, levdimid)
      status = nf_def_dim(id, 'time', nt, timedimid)

      dims (1) = londimid
      status = nf_def_var(id, 'lon', NF_REAL, 1, dims , lonid)
      dims (1) = latdimid
      status = nf_def_var(id, 'lat', NF_REAL, 1, dims , latid)
      dims (1) = levdimid
      status = nf_def_var(id, 'lev', NF_REAL, 1, dims , levid)
      dims (1) = timedimid
      status = nf_def_var(id, 'time', NF_REAL, 1, dims , timeid)

      status = nf_put_att_text(id, lonid, 'units', 12,'degrees_east' )
      status = nf_put_att_text(id, latid, 'units', 13,'degrees_north' )
      status = nf_put_att_text(id, levid, 'units', 3,'hPa' )
      status = nf_put_att_text(id, timeid, 'units', len(trim(time_unit)),trim(time_unit) )

      dims (1) = londimid
      dims (2) = latdimid
      dims (3) = levdimid
      dims (4) = timedimid
      status = nf_def_var(id, 'pv', NF_REAL, 4, dims , qpid)
      status = nf_def_var(id, 'h', NF_REAL, 4, dims , hpid)
      status = nf_def_var(id, 'psi', NF_REAL, 4, dims , spid)
      status = nf_def_var(id, 'u', NF_REAL, 4, dims , urid)
      status = nf_def_var(id, 'v', NF_REAL, 4, dims , vrid)

      status = nf_enddef(id)

      start(1:1)=(/1/)
      total(1:1)=(/nx/)
      status = nf_put_vara_real(id, lonid, start, total, lon_in)
      start(1:1)=(/1/)
      total(1:1)=(/ny/)
      status = nf_put_vara_real(id, latid, start, total, lat_in)
      start(1:1)=(/1/)
      total(1:1)=(/nl/)
      status = nf_put_vara_real(id, levid, start, total, plev_in)
      start(1:1)=(/1/)
      total(1:1)=(/nt/)
      status = nf_put_vara_real(id, timeid, start, total, time_in)

      start(1:4)=(/1,1,1,1/)
      total(1:4)=(/nx,ny,nl,nt/)
      status = nf_put_vara_real(id, qpid, start, total, qp_out)
      status = nf_put_vara_real(id, hpid, start, total, hp_out)
      status = nf_put_vara_real(id, spid, start, total, sp_out)
      status = nf_put_vara_real(id, urid, start, total, ur_out)
      status = nf_put_vara_real(id, vrid, start, total, vr_out)

      status = nf_close(id)


      print*,'end of run'
      END

!C**************************************************************************
       SUBROUTINE BALP( ttt, nx, ny, nl, H,HB,S,SBR,TPR,Q,FCO,FCM,MFC,APS,AC,PE,  \
        HNDM,OMEGS,OMEGH,PART,THRS,BET,FR,MAXX,MAXXT,SIG,HDR, pif,tho,qconst,\
        MI, NOUT, IBC, itsurg, maskfile, qp_out, hp_out, sp_out, ur_out, vr_out)
 
       implicit none
       integer :: nx, ny, nl, NOUT, ttt
       real :: H(NY,NX,NL),HB(NY,NX,NL),HP(NY,NX,NL),SBR(NY,NX,NL)
       real :: S(NY,NX,NL),SP(NY,NX,NL),Q(NY,NX,NL)
       real :: qp_out(nx,ny,nl), hp_out(nx,ny,nl), sp_out(nx,ny,nl)
       real :: ur_out(nx,ny,nl), vr_out(nx,ny,nl)
       real :: MI, PART, SIG, ZM, BET, FR, RS, HNDM, THRS, OMEGS, OMEGH, BI
       real :: SLL(NY,NX,NL),SPP(NY,NX,NL),SLP(NY,NX,NL)
       real :: STB(NY,NX,NL),AVO(NY,NX,NL),RHS(NY,NX,NL),QP(NY,NX,NL)
       real :: ASI(NY,NX,NL),BSI(NY,NX,NL),APHI(NY,NX,NL)
       real :: HRHS(NY,NX,NL),SRHS(NY,NX,NL),OS(NY,NX,NL),OH(NY,NX,NL)
       real :: SISUM(NY,NX,NL),HTSUM(NY,NX,NL),mask(NX,NY,NL)
       real :: TPR(NY,NX,2), TP(NY,NX,2)
       real :: ZNC(NY,NX),FCO(NY,NX),AC(NY,5),FCM(NY,NX),MFC(NY,NX)
       real :: APS(NY),BB(NL),PE(NL),BH(NL),BL(NL),DPI2(NL),PIF(NL)
       real :: XHDR(8), HDR(8)
       real :: UR(NY,NX,NL),VR(NY,NX,NL)
       real :: R1BS, R1BH, R2BS, R2BH, R1PS, R1PH, R2PS, R2PH, ZMRS
       real :: RSA, SXX, SYY, SXY, BETAS, ZSI, RH1, RH2, tho, qconst
       integer :: QLV(NL),SIOUT(NL),HOUT(NL)
       integer :: isurg1(100,NL),isurg2(100,NL)
       integer :: NMLV, GPTS, itsurg, jsurg1, jsurg2, klev1, klev2
       integer :: IBC, IITOT, ITC, ITCC, MAXX, MAXXT
       character*120 :: maskfile
       character*30 IGFIL
       logical :: IT,ICON

       integer :: i, j, k, ivneg, IH, KL 

       SISUM(:,:,:) = 0.
       HTSUM(:,:,:) = 0.

       NMLV = nl
       do k = 1, nl
         QLV(k) = k
       enddo
 
       WRITE(6,*)THRS,OMEGS,OMEGH,PART
       GPTS=NY*NX*NL
       write(6,600)FR,SIG
 600   FORMAT(' FR=',F10.3,' AND SIG=',F10.3)
       WRITE(6,601)GPTS,NY,NX,NL
 601   FORMAT(I6,' gridpoints in domain;',I4,' X',I4,' X',I4)
       !C********** Set coefficients ********************************
       DO I=2,NY-1
         DO J=2,NX-1
           ZNC(I,J)=2.*FR*SIG*SIG*MFC(I,J)/(APS(I)*APS(I))
         enddo 
       enddo 

       DO 202 K=2,NL-1
         BB(K)=-2./( (PE(K+1)-PE(K))*(PE(K)-PE(K-1)) )
         BH(K)=2./( (PE(K+1)-PE(K))*(PE(K+1)-PE(K-1)) )
         BL(K)=2./( (PE(K+1)-PE(K-1))*(PE(K)-PE(K-1)) )
         DPI2(K)=(PE(K+1)-PE(K-1))/2.
         ivneg = 0
         DO 203 J=2,NX-1
         DO 204 I=2,NY-1
           SLL(I,J,K)=ZNC(I,J)*(SBR(I,J+1,K)+SBR(I,J-1,K)-2.*SBR(I,J,K))
           SPP(I,J,K)=ZNC(I,J)*(SBR(I-1,J,K)+SBR(I+1,J,K)-2.*SBR(I,J,K))
           AVO(I,J,K)=FCM(I,J) + FR*( AC(I,1)*SBR(I-1,J,K) +                   \
                                 AC(I,2)*SBR(I,J-1,K) + AC(I,3)*SBR(I,J,K) +   \
                                 AC(I,4)*SBR(I,J+1,K) + AC(I,5)*SBR(I+1,J,K) )
           SLP(I,J,K)=ZNC(I,J)*(SBR(I-1,J+1,K)-SBR(I-1,J-1,K)-                 \
                                SBR(I+1,J+1,K)+SBR(I+1,J-1,K))/2.     !COEFF IS REALLY  2./4.
           STB(I,J,K)=BH(K)*HB(I,J,K+1) + BL(K)*HB(I,J,K-1) + BB(K)*HB(I,J,K)
           IF ( AVO(I,J,K).LT.0.0001 ) THEN
!C	      WRITE(6,96)AVO(I,J,K),I,J,K
             AVO(I,J,K)=0.0001
             ivneg=ivneg + 1
 96          FORMAT(' negative vorticity',f10.3,' at i j k=',3I3)
           END IF
           IF (STB(I,J,K).LT.0.0001) THEN
!c	      WRITE(6,95)stb(i,j,k),I,J,K
             stb(i,j,k)=0.0001
 95          FORMAT(' negative stability',f10.3,' at i j k=',3I3)
           END IF
 204     CONTINUE
 203     CONTINUE
         write(6,*)ivneg,k

         !C**************************************************************
         DO 206 J=2,NX-1
         DO 207 I=2,NY-1
           ASI(I,J,K)=BB(K)*AVO(I,J,K)/(FR*STB(I,J,K)*AC(I,3))
           BSI(I,J,K)=1. + ASI(I,J,K)*FCO(I,J)
           BI=FCO(I,J)*AC(I,3) - 2.*(SLL(I,J,K) + SPP(I,J,K)) 
           IF (BI.GT.0.) THEN
             BI=0.
           END IF
           APHI(I,J,K)=BI/(FR*STB(I,J,K)*AC(I,3))
 207     CONTINUE
 206     CONTINUE
 202   CONTINUE

       !C******** Determine the number of fields to solve for *******
       ! NOUT: THE NUMBER OF INVERSIONS TO BE DONE.'
       DO 210 IH=1,NOUT   !Begin total iteration loop
         !C******* Initialize fields, set boundary conditions *******
         DO J=1,NX
         DO I=1,NY
           TP(I,J,1)=0.
           TP(I,J,2)=0.
           DO K=2,NL-1
             QP(I,J,K)=0.
           enddo
         enddo
         enddo
  
         IF (NMLV.EQ.NL) THEN    !USE ALL LEVELS
           DO J=1,NX
           DO I=1,NY
             TP(I,J,1)=TPR(I,J,1)
             TP(I,J,2)=TPR(I,J,2)
             DO K=2,NL-1
               QP(I,J,K)=Q(I,J,K)
             enddo
           enddo
           enddo
           GO TO 114
         END IF
        !C************************************
         DO 240 KL=1,NMLV
           IF (QLV(KL).EQ.1) THEN
             WRITE(6,602)kl,qlv(kl)
 602         FORMAT(' level ',I4,' field',i4)
             DO 241 J=1,NX
             DO 242 I=1,NY
               TP(I,J,1)=TPR(I,J,1)
 242         CONTINUE
 241         CONTINUE
           ELSE IF (QLV(KL).EQ.NL) THEN
             WRITE(6,602)kl,qlv(kl)
             DO 243 J=1,NX
             DO 244 I=1,NY
               TP(I,J,2)=TPR(I,J,2)
 244         CONTINUE
 243         CONTINUE
           END IF
           DO 245 K=2,NL-1
             IF (K.EQ.QLV(KL)) THEN
               WRITE(6,602)kl,qlv(kl)
               DO 246 J=1,NX
               DO 247 I=1,NY
                  QP(I,J,K)=Q(I,J,K)
 247           CONTINUE
 246           CONTINUE
             END IF
 245       CONTINUE
 240     CONTINUE

         !c  ** read the horizontal surgery domain,  cc wu  , 9/1/92
         !c  ** using mask file -- by cwudpf, 2007/12/06
         ! itsurg: do surgery or not?  0 : no,  1 : surgery inner, 2 : surgery outer, 3 : mask file
  114    print *, 'itsurg = ', itsurg

         if (itsurg.eq.0) go to 1610
                                                       
         !c  ** using mask file
         if (itsurg.eq.3) then
           write (6, *) ' mask file = ', trim(maskfile)

           open (70, file=trim(maskfile), status='old', form='unformatted', access='direct', recl=NX*NY*NL*4)
           read (70, rec=1) (((mask(i,j,k), i=1,NX), j=1,NY), k=1,NL)
           close (70)

           do k = 1, NL
             do j = 1, NX
               do i = 1, NY
                 if (mask(j, NY+1-i, k) .LT. 0.) then
                   Qp(i, j, k) = 0.
                   if (k .EQ. 1) then
                     Tp(i, j, 1) = 0.
                   else if (k .EQ. NL) then
                     Tp(i, j, 2) = 0.
                   end if
                 else
                 end if
               enddo
             enddo
           enddo

         else
         !c  read the vertical levels to have surgery      
           write(6,*) 'vertical surgery levels'        
           read(5,*)  kLev1,kLev2       
           write(6,*) ' kLev1 = ',kLev1,' kLev2= ',kLev2

           do k= klev1 , klev2
             write(6,*) ' k= ', k
             write(6,*) 'surgery in i(y)'         
             !read (5,*)  kindex         
             read (5,*)  (isurg1(j,k),j=1,Nx)         
             read (5,*)  (isurg2(j,k),j=1,Nx)
           end do
                               
           if (itsurg .eq. 1) then

              !c  ** check to see if want to surgery the surface boundary thta
              if (KLev1.eq.1) then
                 print * ,' surgery surface thta ',' KLEV1= ',klev1
                 do j=1,Nx                          
                 do i=1,Ny
                  if((i.ge.isurg1(j,kLEV1)).and.(i.le.isurg2(j,kLEv1)))then        
                          Tp(i,j,1)=0.
                  end if 
                 end do
                 end do
                 klev1=klev1+1
              end if

              !c  ** check to see if want to surgery the top boundary thta
              if (KLev2.eq.NL) then
                 print * ,'surgery top thta ', 'KLEV2= ', klev2
                 do j=1,Nx
                 do i=1,Ny
                  if((i.ge.isurg1(j,kLEV2)).and.(i.le.isurg2(j,kLEv2)))then        
                          Tp(i,j,2)=0.
                  end if 
                 end do
                 end do
                 klev2=klev2-1
              end if

              !c  **  surgery the interior pert. pv
              print *, ' surgery interior pv ', ' KLEV1= ',klev1,' KLEV2= ',klev2
              do k= KLev1, KLev2
                do j=1,Nx
                   do i=1,Ny
                      if( (i.ge.isurg1(j,k)).and.(i.le.isurg2(j,k)) ) then
                         Qp(i,j,k)=0.
                      end if 
                   end do
                end do
              end do

           else if(itsurg .eq. 2) then
             !c  ** check to see if want to surgery the surface boundary thta
             if (KLev1.eq.1) then
                print * ,' surgery surface thta ',' KLEV1= ',klev1
                do j=1,Nx
                do i=1,Ny
                 if((i.lt.isurg1(j,kLEV1)).or.(i.gt.isurg2(j,kLEv1)))then
                         Tp(i,j,1)=0.
                 end if
                end do
                end do
                klev1=klev1+1
             end if

             !c  ** check to see if want to surgery the top boundary thta
             if (KLev2.eq.NL) then
                print * ,'surgery top thta ', 'KLEV2= ', klev2
                do j=1,Nx
                do i=1,Ny
                 if((i.lt.isurg1(j,kLEV2)).or.(i.gt.isurg2(j,kLEv2)))then
                         Tp(i,j,2)=0.
                 end if
                end do
                end do
                klev2=klev2-1
             end if

             !c  **  surgery the interior pert. pv
             print *, ' surgery interior pv ', ' KLEV1= ',klev1,' KLEV2= ',klev2
             do k= KLev1, KLev2
               do j=1,Nx
                  do i=1,Ny
                     if( (i.lt.isurg1(j,k)).or.(i.gt.isurg2(j,k)) ) then
                        Qp(i,j,k)=0.
                     end if
                  end do
               end do
             end do

           end if

         end if
1610     continue
 !c********************************************************************
         ! IBC:  "1" FOR TOTAL PERT B.C., "0" FOR HOMOGEN B.C.
         IF (IBC.EQ.2) THEN
            PRINT*,'ENTER FILENAME.'
            READ(5,*)IGFIL
            OPEN(24,FILE=IGFIL,STATUS='OLD')
            READ(24,*)(XHDR(I),I=1,8)
            DO 250 K=1,NL
              DO 251 I=1,NY
                READ(24,*)(HP(I,J,K),J=1,NX)
                DO 2511 J=1,NX
                 HP(I,J,K)=HP(I,J,K)/HNDM
  2511          CONTINUE
  251         CONTINUE
  250       CONTINUE
            DO 252 K=1,NL
              DO 253 I=1,NY
                READ(24,*,END=253)(SP(I,J,K),J=1,NX)
                DO 2531 J=1,NX
                  SP(I,J,K)=SP(I,J,K)/HNDM
  2531          CONTINUE
  253         CONTINUE
  252       CONTINUE

         ELSE IF (IBC.EQ.1) THEN
            do k=1,nl
              do j=1,nx
                do i=1,ny
                  HP(I,J,K)=0.
                  SP(I,J,K)=0.
                enddo
              enddo
            enddo   
 
            do 254 k=qlv(1),qlv(nmlv) 
              !print *,'k, HTSUM = ',k,  HTSUM(10,10,k)
              DO 255 J=1,NX
              DO 256 I=1,NY
                HP(I,J,K)=H(I,J,K) - HTSUM(I,J,K) !INITIALIZE H1
                SP(I,J,K)=S(I,J,K) - SISUM(I,J,K) !INITIALIZE S1
  256         CONTINUE
  255         CONTINUE
              !print *,'k, HP = ', k,  HP(10,10,k)
  254       CONTINUE
         ELSE
            DO 257 K=1,NL
              DO 258 J=1,NX
              DO 259 I=1,NY
                HP(I,J,K)=0. !INITIALIZE H1
                SP(I,J,K)=0. !INITIALIZE S1
  259         CONTINUE
  258         CONTINUE
  257       CONTINUE
         END IF
         !C********* Calculate upper and lower initial boundary values ***********
         DO 260 J=1,NX
         DO 261 I=1,NY
            HP(I,J,1)=HP(I,J,2) + TP(I,J,1)*(PE(2)-PE(1))
            SP(I,J,1)=SP(I,J,2) + TP(I,J,1)*(PE(2)-PE(1))
            HP(I,J,NL)=HP(I,J,NL-1) - TP(I,J,2)*(PE(NL)-PE(NL-1))
            SP(I,J,NL)=SP(I,J,NL-1) - TP(I,J,2)*(PE(NL)-PE(NL-1))
  261    CONTINUE
  260    CONTINUE

         !C*********** Begin iterations ***************************
         IITOT=0
         ITC=0
  900    CONTINUE     !Total iteration
 
         !C********** Calculate the RHS of the psi equation *********
         DO 270 K=1,NL
           DO 271 J=1,NX
           DO 272 I=1,NY
             OS(I,J,K)=SP(I,J,K)
             OH(I,J,K)=HP(I,J,K)
  272      CONTINUE
  271      CONTINUE
  270    CONTINUE
 
         DO 280 K=2,NL-1
           DO 281 J=2,NX-1
           DO 282 I=2,NY-1
             R1BS=( SBR(I,J+1,K+1)-SBR(I,J-1,K+1)-SBR(I,J+1,K-1)+ SBR(I,J-1,K-1) )/(4.*DPI2(K))
             R1BH=( HB(I,J+1,K+1)-HB(I,J-1,K+1)-HB(I,J+1,K-1)+ HB(I,J-1,K-1) )/(4.*DPI2(K))
             R2BS=( SBR(I-1,J,K+1)-SBR(I+1,J,K+1)-SBR(I-1,J,K-1)+ SBR(I+1,J,K-1) )/(4.*DPI2(K))
             R2BH=( HB(I-1,J,K+1)-HB(I+1,J,K+1)-HB(I-1,J,K-1)+ HB(I+1,J,K-1) )/(4.*DPI2(K))
             R1PS=( SP(I,J+1,K+1)-SP(I,J-1,K+1)-SP(I,J+1,K-1)+ SP(I,J-1,K-1) )/(4.*DPI2(K))
             R1PH=( HP(I,J+1,K+1)-HP(I,J-1,K+1)-HP(I,J+1,K-1)+ HP(I,J-1,K-1) )/(4.*DPI2(K))
             R2PS=( SP(I-1,J,K+1)-SP(I+1,J,K+1)-SP(I-1,J,K-1)+ SP(I+1,J,K-1) )/(4.*DPI2(K))
             R2PH=( HP(I-1,J,K+1)-HP(I+1,J,K+1)-HP(I-1,J,K-1)+ HP(I+1,J,K-1) )/(4.*DPI2(K))
             RHS(I,J,K)=QP(I,J,K) + FR*( (R1BS*R1PH + R1BH*R1PS)/(APS(I)*APS(I)) + SIG*SIG*(R2BS*R2PH + R2BH*R2PS) )
             SRHS(I,J,K)=( RHS(I,J,K) - AVO(I,J,K)*(BH(K)*HP(I,J,K+1)+ \
               BL(K)*HP(I,J,K-1)) )/(FR*STB(I,J,K)) + ASI(I,J,K)*      \
               ( AC(I,1)*HP(I-1,J,K) + AC(I,2)*HP(I,J-1,K) + AC(I,4)*  \
                HP(I,J+1,K) + AC(I,5)*HP(I+1,J,K) )
  282      CONTINUE
  281      CONTINUE
  280    CONTINUE
         !C************* Iteration for psi (2-D) **********************
         ITCC=0
         DO 290 K=2,NL-1
           ITC=0
  800      IT=.TRUE.
           ZMRS=0
           DO 291 J=2,NX-1
           DO 292 I=2,NY-1
              RSA=AC(I,1)*SP(I-1,J,K) + AC(I,2)*SP(I,J-1,K) + AC(I,3)*SP(I,J,K) + AC(I,4)*SP(I,J+1,K) + AC(I,5)*SP(I+1,J,K) 
              SXX=SP(I,J+1,K)+SP(I,J-1,K) -2.*SP(I,J,K)
              SYY=SP(I-1,J,K)+SP(I+1,J,K)-2.*SP(I,J,K)
              SXY=( SP(I-1,J+1,K)-SP(I+1,J+1,K)-SP(I-1,J-1,K)+ SP(I+1,J-1,K) )/4.
              BETAS=SIG*SIG*(FCO(I-1,J)-FCO(I+1,J))* (SP(I-1,J,K)-SP(I+1,J,K))/4. \
                          + (FCO(I,J+1)-FCO(I,J-1))* (SP(I,J+1,K)-SP(I,J-1,K))/4.
              RS=BSI(I,J,K)*RSA + ASI(I,J,K)*( BETAS + SLL(I,J,K)*SYY +    \
                               SPP(I,J,K)*SXX - SLP(I,J,K)*SXY ) - SRHS(I,J,K)
 
              ZSI=SP(I,J,K)
              SP(I,J,K)=ZSI - OMEGS*RS/( BSI(I,J,K)*AC(I,3) - 2.*ASI(I,J,K)*(SLL(I,J,K) + SPP(I,J,K)) )
              ZMRS=ZMRS + ABS(ZSI-SP(I,J,K))/GPTS

              IF (ABS(SP(I,J,K)-ZSI).GT.THRS) THEN
               IT=.FALSE.
              END IF
  292      CONTINUE
  291      CONTINUE
 
           ITC=ITC+1
           IF (AMOD(FLOAT(ITC),10.).EQ.0) THEN
             WRITE(6,*) ITC,ZMRS
           END IF
          
           IF (IT) THEN
             ICON=.TRUE.
             WRITE(6,603)ITC,K
             IF (ITC.GT.1) ITCC=1
  603        FORMAT(I4,' ITERATIONS AT LEVEL',I4)
           ELSE 
             IF (ITC.LT.MAXX) THEN
               GO TO 800
             ELSE
               PRINT*,'TOO MANY ITERATIONS FOR PSI.'
               ICON=.FALSE.
               GO TO 901
             END IF
           END IF
 
  290    CONTINUE
 
         PRINT*,'PSI CONVERGED.'
         IF (IITOT.GT.0) THEN
           DO 300 K=1,NL
             DO 301 J=2,NX-1
             DO 302 I=2,NY-1
               SP(I,J,K)=PART*SP(I,J,K) + (1.-PART)*OS(I,J,K)
  302        CONTINUE
  301        CONTINUE
  300      CONTINUE
         END IF
 
         !C********* Calculate the RHS of the PV+BALANCE equation *******
  700    DO 310 K=2,NL-1
           DO 311 J=2,NX-1
           DO 312 I=2,NY-1
             RH1=(2./AC(I,3))*(SLL(I,J,K) + SPP(I,J,K))* ( AC(I,1)*SP(I-1,J,K) + AC(I,2)*SP(I,J-1,K) + \
                                         AC(I,4)*SP(I,J+1,K) + AC(I,5)*SP(I+1,J,K) )
             BETAS=SIG*SIG*(FCO(I-1,J)-FCO(I+1,J))*                   \
                  (SP(I-1,J,K)-SP(I+1,J,K))/4. + (FCO(I,J+1)-FCO(I,J-1))*(SP(I,J+1,K)-SP(I,J-1,K))/4.
             RH2=BETAS + SLL(I,J,K)*(SP(I-1,J,K)+SP(I+1,J,K)) + SPP(I,J,K)*(SP(I,J-1,K)+SP(I,J+1,K)) - \
                   SLP(I,J,K)*(SP(I-1,J+1,K)-SP(I-1,J-1,K)- SP(I+1,J+1,K)+SP(I+1,J-1,K))/4.
             HRHS(I,J,K)=APHI(I,J,K)*RHS(I,J,K) + RH1 + RH2
  312      CONTINUE
  311      CONTINUE
  310    CONTINUE

         !C************* Solve for phi with 3-D SOR *****************
         ITC=0
         zmrs=0.
  701    IT=.TRUE.
 
         DO 320 K=2,NL-1
         DO 321 J=2,NX-1
         DO 322 I=2,NY-1
           IF (K.EQ.2) THEN
             RS=AC(I,1)*HP(I-1,J,K) + AC(I,2)*HP(I,J-1,K) + \
               (AC(I,3)+ APHI(I,J,K)*(BB(K)+BL(K))*AVO(I,J,K)) * HP(I,J,K) + \
                AC(I,4)*HP(I,J+1,K) + AC(I,5)* HP(I+1,J,K) + \
                APHI(I,J,K)*AVO(I,J,K)*( BH(K)* HP(I,J,K+1) + TP(I,J,1)/DPI2(K) ) - HRHS(I,J,K)
             ZM=HP(I,J,K)
             HP(I,J,K)=ZM - OMEGH*RS/( AC(I,3) + APHI(I,J,K)*(BB(K)+BL(K))*AVO(I,J,K) )
 
           ELSE IF (K.EQ.NL-1) THEN
             RS=AC(I,1)*HP(I-1,J,K) + AC(I,2)*HP(I,J-1,K) + \
               ( AC(I,3) + APHI(I,J,K)*(BB(K)+BH(K))*AVO(I,J,K) )* HP(I,J,K) + \
                AC(I,4)*HP(I,J+1,K) + AC(I,5)*HP(I+1,J,K) + \
                APHI(I,J,K)*AVO(I,J,K)*( BL(K)* HP(I,J,K-1) - TP(I,J,2)/DPI2(K) ) - HRHS(I,J,K)
             ZM=HP(I,J,K)
             HP(I,J,K)=ZM - OMEGH*RS/( AC(I,3) + APHI(I,J,K)* (BB(K)+BH(K))*AVO(I,J,K) )
 
           ELSE
             RS=AC(I,1)*HP(I-1,J,K) + AC(I,2)*HP(I,J-1,K) + \
              ( AC(I,3) + APHI(I,J,K)*BB(K)*AVO(I,J,K) )* HP(I,J,K) + \
                AC(I,4)*HP(I,J+1,K) + AC(I,5)*HP(I+1,J,K) +  \
                APHI(I,J,K)*AVO(I,J,K)*( BH(K)*HP(I,J,K+1) + BL(K)*HP(I,J,K-1) ) - HRHS(I,J,K)
             ZM=HP(I,J,K)
             HP(I,J,K)=ZM-OMEGH*RS/( AC(I,3) + BB(K)*APHI(I,J,K)*AVO(I,J,K) )
           END IF
 
           ZMRS=ZMRS + ABS(ZM-HP(I,J,K))/GPTS
           IF (ABS(ZM-HP(I,J,K)).GT.THRS/2.) THEN
             IT=.FALSE.
           END IF
  322    CONTINUE
  321    CONTINUE
  320    CONTINUE
 
         IF (AMOD(FLOAT(ITC),10.).EQ.0) THEN
           WRITE(6,*)ITC,ZMRS
         END IF
         ZMRS=0.
         
         ITC=ITC+1
         IF (IT) THEN
           PRINT*,'PHI CONVERGED.'
           WRITE(6,606)ITC
  606      FORMAT(' Phi converged in',I4,' iterations.')
           DO 330 J=1,NX
           DO 331 I=1,NY
             HP(I,J,1)=HP(I,J,2) + TP(I,J,1)*(PE(2)-PE(1))
             SP(I,J,1)=SP(I,J,2) + TP(I,J,1)*(PE(2)-PE(1))
             HP(I,J,NL)=HP(I,J,NL-1) - TP(I,J,2)*(PE(NL)-PE(NL-1))
             SP(I,J,NL)=SP(I,J,NL-1) - TP(I,J,2)*(PE(NL)-PE(NL-1))
  331      CONTINUE
  330      CONTINUE
           IF (IITOT.GT.0) THEN
             DO 340 K=1,NL
             DO 341 J=2,NX-1
             DO 342 I=2,NY-1
               HP(I,J,K)=PART*HP(I,J,K) + (1.-PART)*OH(I,J,K)
  342        CONTINUE
  341        CONTINUE
  340        CONTINUE
           END IF
           IF ((ITC.EQ.1).AND.(ITCC.EQ.0)) THEN
             PRINT*,'TOTAL CONVERGENCE.'
           ELSE
             IITOT=IITOT + 1
             WRITE(6,22)IITOT
  22         FORMAT(I4,' TOTAL ITERATION(S).')
             IF (IITOT.GT.MAXXT) THEN
               PRINT*,'TOO MANY TOTAL ITERATIONS.'
               GO TO 901
             ELSE
               GO TO 900
             END IF
           END IF
         ELSE 
           IF (ITC.LT.MAXX) THEN
             GO TO 701
           ELSE
             PRINT*,'TOO MANY ITERATIONS FOR HGHT.'
             ICON=.FALSE.
             GO TO 901
           END IF
         END IF

         !C********** Write out phi and psi fields ***********************
  901    continue
         DO 350 K=1,NL
         DO 351 J=1,NX
         DO 352 I=1,NY
           SISUM(I,J,K)=SISUM(I,J,K) + SP(I,J,K)
           HTSUM(I,J,K)=HTSUM(I,J,K) + HP(I,J,K)
  352    CONTINUE
  351    CONTINUE
  350    CONTINUE
 
  210  CONTINUE

       !********* cwuhuang caculate the balance u,v
       print*,'call ubal'
       print*, 'hdr=',hdr
 
         do k=1,nl
           do i=1,ny
           do j=1,nx
             sp(i,j,k)=sp(i,j,k)*hndm
             hp(i,j,k)=hp(i,j,k)*hndm
           enddo
           enddo
         enddo
 
       call ubal(sp,hdr,ur,vr,nx,ny,nl)
       print*,'call ubal ok!!'
 
       !********* cwuhuang write out data
       do i=1,nx
         do j=1,ny
           qp(j,i,1)=tp(j,i,1)*tho
           qp(j,i,nl)=tp(j,i,2)*tho
           do k=2,nl-1
             qp(j,i,k)=qp(j,i,k)*1.e2*qconst/pif(k)
           enddo
         enddo
       enddo

       !   write(30,rec=(ntt-1)*5+1),(((qp(j,i,k),i=1,nx),j=ny,1,-1),k=1,nl)
       !   write(30,rec=(ntt-1)*5+2),(((hp(j,i,k),i=1,nx),j=ny,1,-1),k=1,nl)
       !   write(30,rec=(ntt-1)*5+3),(((sp(j,i,k),i=1,nx),j=ny,1,-1),k=1,nl)
       !   write(30,rec=(ntt-1)*5+4),(((ur(j,i,k),i=1,nx),j=ny,1,-1),k=1,nl)
       !   write(30,rec=(ntt-1)*5+5),(((vr(j,i,k),i=1,nx),j=ny,1,-1),k=1,nl)

       do k = 1, nl
          do j = ny,1,-1
             do i = 1, nx
                qp_out(i,ny-j+1,k) = qp(j,i,k)
                hp_out(i,ny-j+1,k) = hp(j,i,k)
                sp_out(i,ny-j+1,k) = sp(j,i,k)
                ur_out(i,ny-j+1,k) = ur(j,i,k)
                vr_out(i,ny-j+1,k) = vr(j,i,k)
              enddo
           enddo
       enddo

       RETURN
       end subroutine BALP
!************************************************************************************
       subroutine ubal(psi,hdr,ur,vr,nx,ny,nz)
       !c  **  This program calculated the nondivergent (balanced) wind from the data balp*
 
       integer :: nx, ny, nz
       integer :: kx, ky
       real :: rlat1, rlat2, rlon1, rlon2, dlat ,rlat(ny),ap(ny)
       real :: psi(Ny,Nx,nz),ur(ny,nx,nz),vr(ny,nx,nz)
       real :: hdr(8) 
       real :: pi, aa, dphi, dlam 
 
       kx = Hdr(7)
       ky =  hdr(8)
       dlat= hdr(6)
       dlon= hdr(5)
       rlat1 = hdr(1)
       rlon1 = Hdr(2)
       rlat2 = Hdr(3)                 
       rlon2 = Hdr(4)
 
       pi=4.*atan(1.)       
       aa = 2.e7/pi
       dphi= dlat*pi/180.
       dlam= dlon*pi/180.
  
       !c  **  calculat the wind
       do i=2,ky-1                  
 
         rlat(i)= rlat2-float(i-1)*dlat
         ap(i)= cos(pi*rlat(i)/180.)                          
 
         !c  ** non-divergent wind
         do k=1,nz
            do j=2,kx-1
               ur(i,j,k)= -1.e5 * (psi(i-1,j,k)-psi(i+1,j,k))/ (2.*aa*dphi)
               vr(i,j,k)=  1.e5 * (psi(i,j+1,k)-psi(i,j-1,k))/ (2.*aa*ap(i)*dlam)
            end do
         end do
 
       end do
 
       !C  ** extrapolate the wind at boundary
       do k= 1,nz
          do i=1,ky
 
             ur(i,1,k)= 2.*ur(i,2,k)-ur(i,3,k)
             vr(i,1,k)= 2.*vr(i,2,k)-vr(i,3,k)
 
             ur(i,kx,k)= 2.*ur(i,kx-1,k)-ur(i,kx-2,k)
             vr(i,kx,k)= 2.*vr(i,kx-1,k)-vr(i,kx-2,k)
  
          end do
 
          do j=1,kx
 
             ur(1,j,k)= 2.*ur(2,j,k)-ur(3,j,k)
             vr(1,j,k)= 2.*vr(2,j,k)-vr(3,j,k)
 
             ur(ky,j,k)= 2.*ur(ky-1,j,k)-ur(ky-2,j,k)
             vr(ky,j,k)= 2.*vr(ky-1,j,k)-vr(ky-2,j,k)
  
          end do
 
          ur(1,1,k)=0.5*( ur(1,2,k)+ur(2,1,k) )
          ur(ky,1,k)=0.5*( ur(ky,2,k)+ur(ky-1,1,k) )
          ur(1,kx,k)=0.5*( ur(2,kx,k)+ur(1,kx-1,k) )
          ur(ky,kx,k)=0.5*( ur(ky-1,kx,k)+ur(ky,kx-1,k) )
 
          vr(1,1,k)=0.5*( vr(1,2,k)+vr(2,1,k) )
          vr(ky,1,k)=0.5*( vr(ky,2,k)+vr(ky-1,1,k) )
          vr(1,kx,k)=0.5*( vr(2,kx,k)+vr(1,kx-1,k) )
          vr(ky,kx,k)=0.5*( vr(ky-1,kx,k)+vr(ky,kx-1,k) )
 
        end do
 
        return
 
        end subroutine ubal
!************************************************************************************
     subroutine open_ncfile( iflnm, ncid )
      implicit         none
#include <netcdf.inc>
      character*(*), intent(in)::  iflnm
      integer, intent(out)::      ncid
      integer::  status

      status = nf_open (iflnm, NF_NOWRITE, ncid)
      if (status .ne. NF_NOERR) call handle_err(status)

     end subroutine open_ncfile
!************************************************************************************
     subroutine close_ncfile( iflnm, ncid )
      implicit         none
#include <netcdf.inc>
      character*(*), intent(in)::  iflnm
      integer, intent(in)::      ncid
      integer::  status

      status = nf_close (ncid)
      if (status .ne. NF_NOERR) call handle_err(status)

     end subroutine close_ncfile
!************************************************************************************
     subroutine get_ncdim1( ncid, var1_name, im )
! Get dimension of variable
      implicit         none
#include <netcdf.inc>
      integer, intent(in):: ncid
      character*(*), intent(in)::  var1_name
      integer, intent(out):: im
      integer::  status, var1id

      status = nf_inq_dimid (ncid, var1_name, var1id)
      if (status .ne. NF_NOERR) call handle_err(status)

      status = nf_inq_dimlen (ncid, var1id, im)
      if (status .ne. NF_NOERR) call handle_err(status)

     end subroutine get_ncdim1
!************************************************************************************
     subroutine handle_err(status)
      implicit         none
#     include          <netcdf.inc>
      integer          status

      if (status .ne. nf_noerr) then
        print *, nf_strerror(status)
        stop 'Stopped'
      endif

     end subroutine handle_err
!************************************************************************************
     subroutine get_att_text( ncid, var1_name, att_name, att )
      implicit         none
#include <netcdf.inc>
      integer, intent(in):: ncid
      character*(*), intent(in)::  var1_name, att_name
      character*(*), intent(out)::  att

      integer::  status, var1id

      status = nf_inq_varid (ncid, var1_name, var1id)
      if (status .ne. NF_NOERR) call handle_err(status)
!     write(*,*) 'Got var1id', var1id

      status = nf_get_att_text (ncid, var1id, att_name, att)
      if (status .ne. NF_NOERR) call handle_err(status)

     end subroutine get_att_text
!************************************************************************************
     subroutine get_var1_real( ncid, var1_name, im, var1 )
      implicit         none
#include <netcdf.inc>
      integer, intent(in):: ncid
      character*(*), intent(in)::  var1_name
      integer, intent(in):: im
      real*4, intent(out):: var1(im)

      integer::  status, var1id

      status = nf_inq_varid (ncid, var1_name, var1id)
      if (status .ne. NF_NOERR) call handle_err(status)
!     write(*,*) 'Got var1id', var1id

      status = nf_get_var_real (ncid, var1id, var1)
      if (status .ne. NF_NOERR) call handle_err(status)

     end subroutine get_var1_real
!************************************************************************************
      subroutine get_var4_real( ncid, var4_name, im, jm, km, nt, var4 )
      implicit         none
#include <netcdf.inc>
      integer, intent(in):: ncid
      character*(*), intent(in)::  var4_name
      integer, intent(in):: im, jm, km, nt
      real*4:: wk4(im,jm,km,4)
      real*4, intent(out):: var4(im,jm)
      integer::  status, var4id
      integer:: start(4), icount(4)
      integer:: i,j

      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = 1

      icount(1) = im    ! all range
      icount(2) = jm    ! all range
      icount(3) = km    ! all range
      icount(4) = nt

      status = nf_inq_varid (ncid, var4_name, var4id)
      status = nf_get_vara_real(ncid, var4id, start, icount, var4)

      if (status .ne. NF_NOERR) call handle_err(status)

      end subroutine get_var4_real
!************************************************************************************

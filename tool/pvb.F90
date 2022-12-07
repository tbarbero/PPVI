! Compiling on analysis:
!    module load netcdf/4.6.1
!    module load intel_compilers/oneapi
!    # Full speed run:
!    ifort -O2 pvb.F90 -o pvb.exe -lnetcdf -lnetcdff -lhdf5 -lhdf5_hl -lz -limf -lm

! This program was provided by Prof. CC Wu's group in NTU in Apr. 2021
! Modified and sorted to Fortran 90 by Jan-Huey Chen
! The original comments are lines below starting with "!!"

      PROGRAM pvpiallnu
!!  (the program is just to calculate the theta,pv,h,psi)...schou
!!
!! add the stability check, and correct the Height field to make sure static stability is positive eveywhere. 
!! Also correct Z at lower and upper boundary to make dZ/d(pi)=-theta
!
!! THIS IS THE MODIFICATION OF CHRIS DAVIS''S PROGRAM PVPI.FOR IT CALCULATES PV ON PI SURFACE 
!! AND INVERT THE STRAM FUNCTION FROM RELATIVE VORTICITY.

        implicit none
#include <netcdf.inc>
        character*190 :: nmlfile = 'pvb.nml'
        integer, external ::  iargc
        character*300 :: infile
        character*300 :: outfile
        character*90 :: time_unit
        integer :: nx, ny, nw, nt 
        real, dimension (:,:,:), allocatable :: U, V, TH, Q, VOR, STB
        real, dimension (:,:,:), allocatable :: DU, DV, DTHX, DTHY, H, TT, PSI
        real, dimension (:,:), allocatable :: THB, THT, A
        real, dimension (:), allocatable :: FC, AP, APM, APP, PI 
        real :: PII, DUDY, DVDX, KAP, COEF, DEL, BETA, VORA, DEF, PIB, PIT
        real :: hdr(6) 
        real :: sigm, Lapsi, Rs, thrs, omegs
        real :: hlow, dhlow, dhlowp, htop, dhtop, dhtopp, STB1, DZ1, DZ, FC125
        real :: VL, UPV, dsum, psi11, psidiff, dpsisum, psin, dpsi, dpsiavgnew
        real :: dpsi_ratio, dpsiavgold, ZSHR
!c  ** add to write out the absoloute vorticity  
        real, dimension (:,:,:), allocatable :: AVOR
!c  *********************************************************
        integer :: imax, icount, hdrr(2), ihdr7, ihdr8
!c  ******* The information about the total numbers of halfdays and
!c  ******  the output halfday number we want   
        integer :: nhalfday,nhalfdayo,nhalfout(30),ic
!c  *********************************
        logical :: IT
!C  *******************
        CHARACTER*70 :: F(80),fum(30)
!c  *********************************
        real, parameter :: CP=1004.5, P0=1.E5, MI=9999.90, gg=9.8066 ! gravity
        real :: AA, DL, DP
        integer :: ll, i, j, k, p, n, iloop, iloopc, l

        integer :: fid
        real, dimension (:), allocatable :: lon_in, lat_in, plev_in, time_in
        real, dimension (:), allocatable :: pr
        real, dimension (:,:,:,:), allocatable :: h_in, t_in, u_in, v_in
        real, dimension (:,:,:,:), allocatable :: u_out, v_out, h_out, psi_out, q_out
        real, dimension (:,:,:,:), allocatable :: avor_out
        integer :: status
        integer :: id, londimid, latdimid, levdimid, timedimid
        integer :: lonid, latid, levid, timeid
        integer :: qid, hid, uid, vid, psiid, avorid, tid
        integer :: dims(4),start(4), total(4)

        namelist /nlist/ outfile, infile, imax, omegs, thrs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Read in namelist and set up range:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (iargc() > 0 ) then
          call getarg(1,nmlfile)
        endif

        open (11, file=trim(nmlfile), status='old')
        read (11, nml=nlist)
        close (11)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! READ h,t,u,v from netcdf files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call open_ncfile(infile,fid)
        print *, 'Reading hb, tb, ub, vb from ', trim(infile), fid

        call get_ncdim1(fid, 'lon', nx)
        call get_ncdim1(fid, 'lat', ny)
        call get_ncdim1(fid, 'lev', nw)
        call get_ncdim1(fid, 'time', nt)
        print *, 'data nlon, nlat, np, nt = ', nx, ny, nw, nt

        allocate(plev_in(nw), lon_in(nx), lat_in(ny), time_in(nt))
        call get_var1_real (fid, 'lev', nw, plev_in)
        call get_var1_real (fid, 'lon', nx, lon_in)
        call get_var1_real (fid, 'lat', ny, lat_in)
        call get_var1_real (fid, 'time', nt, time_in)
        call get_att_text (fid, 'time', 'units', time_unit)

        allocate(pr(nw))
        do k = 1, nw
          pr(k) = plev_in(k) * 0.001
        enddo

        allocate(h_in(nx,ny,nw,nt))
        allocate(t_in(nx,ny,nw,nt))
        allocate(u_in(nx,ny,nw,nt))
        allocate(v_in(nx,ny,nw,nt))
        call get_var4_real (fid, 'hb', nx, ny, nw, nt, h_in)
        call get_var4_real (fid, 'tb', nx, ny, nw, nt, t_in)
        call get_var4_real (fid, 'ub', nx, ny, nw, nt, u_in)
        call get_var4_real (fid, 'vb', nx, ny, nw, nt, v_in)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! allocate variables and set up parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate (U(ny,nx,nw), V(ny,nx,nw), TH(ny,nx,nw), H(ny,nx,nw))
        allocate (VOR(ny,nx,nw), STB(ny,nx,nw), TT(ny,nx,nw))
        allocate (DU(ny,nx,nw), DV(ny,nx,nw), DTHX(ny,nx,nw), DTHY(ny,nx,nw))
        allocate (PSI(ny,nx,nw),Q(ny,nx,nw))
        allocate (THB(ny,nx), THT(ny,nx))
        allocate (AVOR(ny,nx,nw))
        allocate (A(ny,5))
        allocate (FC(ny), AP(ny), APM(ny), APP(ny))
        allocate (PI(nw))

        hdr(1) = lat_in(1)
        hdr(2) = lon_in(1)
        hdr(3) = lat_in(ny)
        hdr(4) = lon_in(nx)
        hdr(5) = abs(lon_in(2) - lon_in(1))
        hdr(6) = abs(lat_in(2) - lat_in(1))

        hdrr(1)= nx
        hdrr(2)= ny
        
        ihdr7 = hdrr(1)
        ihdr8 = hdrr(2)
 
        ! the information for the realtive vorticity inversion (JHC: based on the original comments)
        print *, 'imax, omegs, thrs = ', imax, omegs, thrs

        PII=4.*ATAN(1.)
        KAP=2./7.
        AA=2.E7/PII ! earth radius
        SIGM=hdr(5)/hdr(6)
        ! x-dir point length & y-dir point length
        DL=AA*hdr(5)*PII/180.
        DP=AA*hdr(6)*PII/180.

        ! pi coordinate for vertical coord
        DO K=1,NW
         PI(K)= CP *( PR(K)**KAP )
        END DO

        ! lower-boundary( between 1000mb and 925mb)
        ! upper-boundary
        PIB=0.5*(PI(2)+PI(1))
        PIT=0.5*(PI(NW)+PI(NW-1))

        DO I=1,HDRr(2)
        
         AP(I)=COS( PII*( HDR(3) - (I-1)*HDR(6) )/180. )
         APM(I)=COS( PII*( HDR(3) - (I-1.5)*HDR(6) )/180. )
         APP(I)=COS( PII*( HDR(3) - (I-0.5)*HDR(6) )/180. )
         !****above been studied   linly 7/15 *************

         ! ADD THE LAPLACIAN COEFFICIENT HERE    
         A(I,1)=SIGM*SIGM*APM(I)/AP(I)
         A(I,2)=1./( AP(I)*AP(I) )
         A(I,3)=-( 2. + SIGM*SIGM*AP(I)*(APM(I)+APP(I)) )/( AP(I)*AP(I) )
         A(I,4)=1./( AP(I)*AP(I) )
         A(I,5)=SIGM*SIGM*APP(I)/AP(I)

         FC(I)=(1.458E-4)*SIN( PII*( HDR(3) - (I-1)*HDR(6) )/180. )

         !print*,'i, HDR(3) - (I-1)*HDR(6), fc(i)'
         !print*,i, HDR(3) - (I-1)*HDR(6), fc(i)

        END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! BEGIN DO LOOP FOR HALFDAYS
!! mi:the lateral b.c. pv to 9999.90(means pv can't be caculated)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate ( u_out(nx,ny,nw,nt), v_out(nx,ny,nw,nt), psi_out(nx,ny,nw,nt) )
      allocate ( h_out(nx,ny,nw,nt), q_out(nx,ny,nw,nt), avor_out(nx,ny,nw,nt) )
      !open (30,file='pv.dat',status='new',form='unformatted', access='direct',recl=nx*ny*nw*4)
      do 108 ll=1,nt

         ! JHC: for whatever reason, the following code use different index
         ! order (j, i, k) instead of (i, j, k). Also, index j is from north to
         ! south.
         do k = 1, nw
            do j = ny,1,-1
               do i = 1, nx
                  h(j,i,k)  = h_in(i,ny-j+1,k,ll)/gg
                  th(j,i,k) = t_in(i,ny-j+1,k,ll)
                  u(j,i,k)  = u_in(i,ny-j+1,k,ll)
                  v(j,i,k)  = v_in(i,ny-j+1,k,ll)
               enddo
             enddo
         enddo

         DO K=2,NW-1 
          DO I=1,HDRr(2)
           Q(I,1,K)=MI
           Q(I,ihdr7,K)=MI
          END DO
          DO J=1,HDRr(1)
           Q(1,J,K)=MI
           Q(ihdr8,J,K)=MI
          END DO
         END DO

         ! CALCULATE THE POTENTIAL TEMPERATURE, the data is already in K , 
         ! th now becomes theta,tt is temp
         DO K=1,NW
          DO I=1,HDRr(2)
            DO J=1,HDRr(1)
              TT(i,j,k) = TH(i,j,k)
              TH(I,J,K) = TH(I,J,K)*CP/PI(K)
            END DO
         END DO
        END DO

!C***************************************************************
        DO J=1,HDRr(1)
          DO I=1,HDRr(2)
            THB(I,J)=0.5*(TH(I,J,1) + TH(I,J,2)) 
            THT(I,J)=0.5*(TH(I,J,NW-1) + TH(I,J,NW))
          END DO
        END DO

        ! add the stability check, and correct the Height field to make sure static stability is positive eveywhere,  
        ! try to loop 100 times, make sure no negative stability exists **
        do 834 iloop = 1,100

           iloopc = 0

        ! correct Z at top and bottome to let dZ/D(pi) = -1.*theta  ***
           DO J=1,HDRr(1)
              DO I=1,HDRr(2)
                 Hlow= H(i,j,2) - THB(I,J) * (Pi(1)-PI(2)) /9.81
                 dhlow = hlow - H(i,j,1)
                 dhlowp = dhlow/H(i,j,1)
                 H(i,j,1) = hlow
                 htop= H(i,j,NW-1) - THT(I,J)*(Pi(NW)-PI(NW-1))/9.81
                 dhtop = htop - H(i,j,NW)
                 dhtopp = dhtop/H(i,j,NW)
                 H(i,j,NW) = htop
              END DO
           END DO
 
           DO 831 K=2,NW-1
             !WRITE(6,*)' k= ',k
             DO 152 J=1,HDRr(1)
             DO 152 I=1,HDRr(2)
              STB1=( 2.*9.81/(PI(K+1)-PI(K-1)) )* \
              ( (H(I,J,K+1)-H(I,J,K))/(PI(K+1)-PI(K)) - \
                (H(I,J,K)-H(I,J,K-1))/(PI(K)-PI(K-1)) )
              IF(STB1.LT.0.) THEN
               iloopc=iloopc + 1
                DZ1=(PI(K+1)-PI(K))*(H(I,J,K)-H(I,J,K-1))/(PI(K)-PI(K-1))+ 0.000001
                DZ=DZ1 - (H(I,J,K+1) - H(I,J,K))
                H(I,J,K) = H(I,J,K+1) - DZ1 
              END IF
152           continue
831        continue
           !if already no negative stability, then stop looping
           if (iloopc.eq.0) then
           ! write (6,*) 'no negative stability at iloop= ',iloop
             go to 835
           endif
834     continue
835     continue

        ! let initial psi equal H,   
        !FC125=(1.458E-4)*SIN( 12.5*PII/180.)
        FC125=(1.458E-4)*SIN( 5.*PII/180.)

        DO K=1,NW
         DO I=1,HDRr(2)
          DO J=1,HDRr(1)
            if (fc(i).lt. fc125) then
               psi(i,j,k)=H(i,j,k)*gg/fc125
            else
               psi(I,j,k)=H(i,j,k)*gg/fc(i)
            endif
          END DO
         END DO
        END DO

        DO K=2,NW-1
         DO J=1,HDRr(1)
          DO I=1,HDRr(2)        
            DTHY(I,J,K)=(TH(I-1,J,K)-TH(I+1,J,K))/(2.*DP)
            DTHX(I,J,K)=(TH(I,J+1,K)-TH(I,J-1,K))/(2.*AP(I)*DL)
            STB(I,J,K)=(TH(I,J,K+1)-TH(I,J,K-1))/(PI(K+1)-PI(K-1))
          END DO
         END DO
        END DO

        ! ******** Vorticity *******************
        DO K=1,NW
          DO I=2,HDRr(2)-1
          DO J=2,HDRr(1)-1
            VL=(V(I,J+1,K)-V(I,J-1,K))/(2.*DL*AP(I))
            UPV=(AP(I-1)*U(I-1,J,K)-AP(I+1)*U(I+1,J,K))/(2.*DP*AP(I))
            VOR(I,J,K)=VL - UPV
            AVOR(i,j,k) = FC(I)+ VOR(i,j,k) !  to write out the absoloute vorticity
          END DO
          END DO
        END DO

        ! TO INVERT THE STREAM FUNCTION (PSI) FROM RELATIVE VORTICITY(VOR)
        ! solve the lateral psi
        ! calculate the boundary divergence term
        do k=1,NW

           dsum=0.

           do i=1,Hdrr(2)-1
              dsum=dsum -((u(i,1,k)+u(i+1,1,k))/2.)* DP
              dsum=dsum +((u(i,iHdr7,k)+u(i+1,iHdr7,k))/2.)* DP
           End do

           do j=1,Hdrr(1)-1
              dsum=dsum+((v(1,j,k)+v(1,j+1,k))/2.)* DL * AP(1)
              dsum=dsum-((v(iHdr8,j,k)+v(iHdr8,j+1,k))/2.)*DL*AP(iHdr8)
           End do
           !print *, 'ap(1)',ap(1)
           !print *, 'dp',dp
           !print *,'dl',dl
           !print *,'ap(ihdr8)',ap(ihdr8)
           !print *, 'before****dsum=', dsum
           dsum = dsum /( 2.*DP*(Hdrr(2)-1)+ DL*(HDrr(1)-1)*(Ap(1)+Ap(iHdr8))  )
           !print *, 'dsum=', dsum

           ! choose starting point psi(1,1,k)=H(1,1,k)*g/f             
           ! then integrate by Davis (2.40) to get the whole psi      

           do i=1,Hdrr(2)-1
              psi(i+1,1,k)=psi(i,1,k)+(dsum+((u(i  ,1,k) + u(i+1,1,k) )/2.))* DP
           End do

           do j=1,Hdrr(1)-1
              psi(iHDR8,j+1,k)=psi(ihdr8,j,k)+(dsum+((v(iHdr8,j,k)+v(iHdr8,j+1,k))/2.))*DL*AP(iHdr8)
           End do

           do i=Hdrr(2),2,-1
              psi(i-1,iHdr7,k)=psi(i,iHdr7,k)+(dsum-((u(i,iHdr7,k)+u(i-1,iHdr7,k) )/2.))* DP
           End do

           do j=Hdrr(1),3,-1
              psi(1,j-1,k)=psi(1,j,k)+(dsum-((v(1,j,k)+v(1,j-1,k))/2.))*DL*AP(1)
           End do

           ! should correct the difference for original and final psi (1,1,k)
           ! averagedly add the difference on all grids on boundary,   
           psi11 = psi(1,2,k)+ (dsum-((v(1,2,k)+v(1,1,k))/2.))*DL*AP(1)
           psidiff = (psi11 - psi(1,1,k)) / ( 2.*(hdrr(1)+hdrr(2)) - 4. )
           !print *, 'psidiff= ', psidiff

           do i=2,Hdrr(2)
              psi(i,1,k)= psi(i,1,k)+ psidiff
           End do

           do j=2,Hdrr(1)
              psi(iHDR8,j,k)=psi(ihdr8,j,k) + psidiff
           End do

           do i=Hdrr(2)-1,1,-1
              psi(i,iHdr7,k)=psi(i,iHdr7,k) + psidiff
           End do

           do j=Hdrr(1)-1,1,-1
              psi(1,j,k)=psi(1,j,k) + psidiff
           End do
        End do


       !invert to get the interior psi by overrelaxation
       icount=1
155    it=.true.    

       dpsisum=0.

       do k=1,NW

          do i=2,HDRr(2)-1
          do j=2,HDRr(1)-1

            Lapsi=1./(DL*DL)*( A(I,1)*psi(I-1,J,K) + A(I,2)*psi(I,J-1,K) + \
                               A(I,3)*psi(I,J,K)   + A(I,4)*psi(I,J+1,K) + \
                               A(I,5)*psi(I+1,J,K) )
            RS = LApSi - vor(i,j,k)   
            
            !If ((i.eq.5).and.(j.eq.10).and.(k.eq.3)) then
            !     print *,'lapsi = ',lapsi,' vor = ', vor(i,j,k)
            !endif

            psin=psi(i,j,k)
            psi(i,j,k)= psin-omegs*RS/A(i,3)*(DL*DL)

            Dpsi= Psi(i,j,k)-psin

            Dpsisum= Dpsisum + Dpsi

            If (abs(dpsi).GT.thrs) then
               IT= .false.
            end if

          end do
          end do

       end do
       
       dpsiavgnew=Dpsisum/( (hdrr(2)-2)*(hdrr(1)-2)*(float(nw)-2) )
       
       IF (mod(icount,10).eq.0) then
           dpsi_ratio = dpsiavgnew /dpsiavgold
           !print *, ' icount= ', icount, dpsi_ratio
       endif

       dpsiavgold = dpsiavgnew

       icount = icount + 1

       ! check if too many iterations
       if (icount.gt. IMAX) THEN
          PRINT *, 'TOO MANY ITERATION FOR PSI',i,j,k
          GO TO 123
       ENDIF

       ! check if need to do more realxation
       If (IT) then
          print *, 'icount=', icount
          print *, 'psi converged'
       else
          go to  155       
       ENDIF   

123    continue

       !print*,'******************************************'
       !print*,'ok1'

       !** Vertical wind shear *************************************
       DO K=2,NW-1
         DO I=2,HDRr(2)-1
         DO J=2,HDRr(1)-1
           DU(I,J,K)=(U(I,J,K+1)-U(I,J,K-1))/(PI(K+1)-PI(K-1))
           DV(I,J,K)=(V(I,J,K+1)-V(I,J,K-1))/(PI(K+1)-PI(K-1))
         END DO
         END DO
       END DO
       !******* Calculate  ertel PV **************************************
        COEF=1.E2*1.E6*9.81*KAP*(CP**3.5)/P0
        DO L=2,NW-1
          DO J=2,HDRr(1)-1
          DO I=2,HDRr(2)-1
            ZSHR=COEF*(PI(L)**(-2.5))*( DU(I,J,L)*DTHY(I,J,L) - DV(I,J,L)*DTHX(I,J,L) )
            Q(I,J,L)=-COEF*(PI(L)**(-2.5))*( (FC(I)+VOR(i,j,L))* STB(I,J,L) ) - ZSHR

!C**********!** Check for negative values ******
            IF (Q(I,J,L).LE.0.) THEN
              !print *, 'Q is less than zero ', i,' * ', j,' * ', L
            END IF
          END DO
          END DO
        END DO

        ! WRITE OUT "BOUNDARY" THETA AND PV ****************
        do i=1,nx
           do j=1,ny
              q(j,i,nw)=tht(j,i)
              q(j,i,1)=thb(j,i)
           enddo 
        enddo 
 
        ! to scale psi to fit H
        DO K=1,NW
          DO j=1,HDRr(1)
          do i=1,hdrr(2)
             psi(I,J,K)=psi(I,J,K)/1.e5
          enddo
          END DO                          
        END DO
        !call ubal to calculate the balance wind
        print*, 'call ubal'
        call ubal(u,v,psi,nx,ny,nw,hdr)

        ! write out the pv(q)( 1st lev. pv is thb ,last pv is tht)*
        ! pv ,psi(velocity potential),h( geopential height)******** 

        !write(30,rec=(ll-1)*5+1),(((q(j,i,k),  i=1,nx),j=ny,1,-1),k=1,nw)
        !write(30,rec=(ll-1)*5+2),(((psi(j,i,k),i=1,nx),j=ny,1,-1),k=1,nw)
        !write(30,rec=(ll-1)*5+3),(((h(j,i,k),  i=1,nx),j=ny,1,-1),k=1,nw)
        !write(30,rec=(ll-1)*5+4),(((u(j,i,k),  i=1,nx),j=ny,1,-1),k=1,nw)
        !write(30,rec=(ll-1)*5+5),(((v(j,i,k),  i=1,nx),j=ny,1,-1),k=1,nw)

         do k = 1, nw
            do j = ny,1,-1
               do i = 1, nx
                  q_out(i,ny-j+1,k,ll) = q(j,i,k)
                  h_out(i,ny-j+1,k,ll) = h(j,i,k)
                  u_out(i,ny-j+1,k,ll) = u(j,i,k)
                  v_out(i,ny-j+1,k,ll) = v(j,i,k)
                  psi_out(i,ny-j+1,k,ll) = psi(j,i,k)
                  avor_out(i,ny-j+1,k,ll) = avor(j,i,k)
               enddo
             enddo
         enddo

108   continue

      ! Write out the netcdf file
      print *, 'outfile = ', trim(outfile)
      status = nf_create(outfile, NF_NOCLOBBER, id)

      status = nf_def_dim(id, 'lon',  nx, londimid)
      status = nf_def_dim(id, 'lat',  ny, latdimid)
      status = nf_def_dim(id, 'lev',  nw, levdimid)
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
      status = nf_def_var(id, 'pv', NF_REAL, 4, dims , qid)
      status = nf_def_var(id, 'h',  NF_REAL, 4, dims , hid)
      status = nf_def_var(id, 'u',  NF_REAL, 4, dims , uid)
      status = nf_def_var(id, 'v',  NF_REAL, 4, dims , vid)
      status = nf_def_var(id, 'psi',  NF_REAL, 4, dims , psiid)
      status = nf_def_var(id, 'avor',  NF_REAL, 4, dims , avorid)

      status = nf_enddef(id)
 
      start(1:1)=(/1/)
      total(1:1)=(/nx/)
      status = nf_put_vara_real(id, lonid, start, total, lon_in)
      start(1:1)=(/1/)
      total(1:1)=(/ny/)
      status = nf_put_vara_real(id, latid, start, total, lat_in)
      start(1:1)=(/1/)
      total(1:1)=(/nw/)
      status = nf_put_vara_real(id, levid, start, total, plev_in)
      start(1:1)=(/1/)
      total(1:1)=(/nt/)
      status = nf_put_vara_real(id, timeid, start, total, time_in)
     
      start(1:4)=(/1,1,1,1/)
      total(1:4)=(/nx,ny,nw,nt/)
      status = nf_put_vara_real(id, qid, start, total, q_out)
      status = nf_put_vara_real(id, hid, start, total, h_out)
      status = nf_put_vara_real(id, uid, start, total, u_out)
      status = nf_put_vara_real(id, vid, start, total, v_out)
      status = nf_put_vara_real(id, psiid, start, total, psi_out)
      status = nf_put_vara_real(id, avorid, start, total, avor_out)

      status = nf_close(id)

      deallocate ( u_out, v_out, psi_out, h_out, q_out, avor_out )

      print*,'end of run'

      end
!************************************************************************************
      subroutine ubal(ur,vr,psi,nx,ny,nz,hdr)

      real :: rlat1, rlat2, rlon1, rlon2 
      real :: dlat, dlon, rlat, ap
      real :: psi(Ny,Nx,nz),ur(ny,nx,nz),vr(ny,nx,nz)
      real :: hdr(6)
      real :: pi, aa, dphi, dlam

      ! setup domain infomation ***********
      rlat1 = hdr(1)
      rlon1 = hdr(2)
      rlat2 = hdr(3)
      rlon2 = hdr(4)
      dlat= hdr(6)
      dlon= hdr(5)

      ! set parameter ***********************
      pi=4.*atan(1.)
      aa = 2.e7/pi
      dphi= dlat*pi/180.
      dlam= dlon*pi/180.

      ! calculat the wind
      do i=2,ny-1
         rlat= rlat2-float(i-1)*dlat ! the real latitude
         ap= cos(pi*rlat/180.)       ! the factor of latitude

         ! non-divergent wind
         do k=1,nz
           do j=2,nx-1
             ur(i,j,k)= -1.e5 * (psi(i-1,j,k)-psi(i+1,j,k)) / (2.*aa*dphi)
             vr(i,j,k)=  1.e5 * (psi(i,j+1,k)-psi(i,j-1,k)) / (2.*aa*ap*dlam)
           end do
         end do
      end do

      ! extrapolate the wind at boundary
      do k= 1,nz
         do i=1,ny
            ur(i,1,k)= 2.*ur(i,2,k)-ur(i,3,k)
            vr(i,1,k)= 2.*vr(i,2,k)-vr(i,3,k)

            ur(i,nx,k)= 2.*ur(i,nx-1,k)-ur(i,nx-2,k)
            vr(i,nx,k)= 2.*vr(i,nx-1,k)-vr(i,nx-2,k)
         end do

         do j=1,nx
            ur(1,j,k)= 2.*ur(2,j,k)-ur(3,j,k)
            vr(1,j,k)= 2.*vr(2,j,k)-vr(3,j,k)

            ur(ny,j,k)= 2.*ur(ny-1,j,k)-ur(ny-2,j,k)
            vr(ny,j,k)= 2.*vr(ny-1,j,k)-vr(ny-2,j,k)
         end do

         ur(1,1,k)=0.5*( ur(1,2,k)+ur(2,1,k) )
         ur(ny,1,k)=0.5*( ur(ny,2,k)+ur(ny-1,1,k) )
         ur(1,nx,k)=0.5*( ur(2,nx,k)+ur(1,nx-1,k) )
         ur(ny,nx,k)=0.5*( ur(ny-1,nx,k)+ur(ny,nx-1,k) )

         vr(1,1,k)=0.5*( vr(1,2,k)+vr(2,1,k) )
         vr(ny,1,k)=0.5*( vr(ny,2,k)+vr(ny-1,1,k) )
         vr(1,nx,k)=0.5*( vr(2,nx,k)+vr(1,nx-1,k) )
         vr(ny,nx,k)=0.5*( vr(ny-1,nx,k)+vr(ny,nx-1,k) )

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

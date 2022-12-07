!! Compiling on analysis:
!    module load netcdf/4.6.1
!    module load intel_compilers/oneapi
!    ifort -O2 symhtuv1.F90 -o symhtuv.exe -lnetcdf -lnetcdff -lhdf5 -lhdf5_hl -lz -limf -lm

!!! Modified by Jan-Huey Chen in May 2021 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Theis program is to caculate the symmetric u,v and then
!c     Using the Non-linear Balance eq. to get a balance Height
!c     and use Hydro-static eq. to get balance T
!c     so that we have symmetric h, t, u, v 
!c     and asymmetric h, t, u, v
!c     we take this as the basic field for the PV
!c                                               Cwuhuang, 10,30,99"
!c                                               linly,  11,19,99' 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     modified by cwuhuang that we combine the symmetric process
!c     in this program.
!c                                             C.-X. Huang 12,02,99"
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program symbasic
      implicit none
#include <netcdf.inc>
      
      !character*30 :: FMT = "(12X,1x,f5.1,f5.1)"
      !FMT = "(12X,1X,F5.1,1X,F5.1)"

      character*190 :: nmlfile = 'symhtuv.nml'
      integer, external ::  iargc
      character*300 :: infile1, infile2, trackfile
      character*300 :: outfile1, outfile2
      character*90 :: time_unit

      integer :: nx, ny, nz, nnt

      real, dimension (:,:,:), allocatable :: h, t, u, v
      real, dimension (:,:,:), allocatable :: si 
      real, dimension (:,:,:), allocatable :: hb, tb, ub, vb
      real, dimension (:,:,:), allocatable :: sib 
      real, dimension (:,:,:), allocatable :: hp, tp, up, vp
      real, dimension (:,:,:), allocatable :: sip
      real, dimension (:,:,:), allocatable :: rhs, tmp
      real :: hdr(6)
      real :: omegs, thrs, thrs_h, rmax
      integer :: imax, nmx

      real, parameter :: gg = 9.8066
 
      namelist /nlist/ infile1, infile2, trackfile, outfile1, outfile2, &
                       nmx, rmax, imax, omegs, thrs_h

      integer :: i, j, k, nt
      integer :: fid1, fid2
      real, dimension (:), allocatable :: lon_in, lat_in, plev_in, time_in
      real, dimension (:), allocatable :: pr
      real, dimension (:,:,:,:), allocatable :: t_in, h_in, si_in
      real :: hlat, hlon
      real, dimension(:,:,:,:), allocatable :: hb_out, tb_out, ub_out, vb_out
      real, dimension(:,:,:,:), allocatable :: hp_out, tp_out, up_out, vp_out
      real, dimension(:,:,:,:), allocatable :: h_out, t_out, u_out, v_out
      integer :: status
      integer :: id, londimid, latdimid, levdimid, timedimid
      integer :: lonid, latid, levid, timeid
      integer :: hid, tid, uid, vid
      integer :: hpid, tpid, upid, vpid
      integer :: hbid, tbid, ubid, vbid
      integer :: dims(4),start(4), total(4)

      !! Read in namelist
      if (iargc() > 0 ) then
        call getarg(1,nmlfile)
      endif

      open (11, file=trim(nmlfile), status='old')
      read (11, nml=nlist)
      close (11)

      print*, 'Input file for h,t,u,v:'
      print*, trim(infile1)
      print*, 'Input file for h, psi (balance):'
      print*, trim(infile2)
      print*, 'The number of direction to do azimuthal mean = ', nmx
      print*, 'The radius(km) to caculate symmetric u,v = ', rmax
      print*, 'Max. iteration for psi and h (nonlinear balance eq.) = ', imax
      print*, 'Over-relaxation parameter = ', omegs
      print*, 'Threshold of h = ', thrs_h

      !! Read in netcdf input files

      call open_ncfile(infile1,fid1)
    
      call get_ncdim1(fid1, 'lev', nz)
      call get_ncdim1(fid1, 'lon', nx)
      call get_ncdim1(fid1, 'lat', ny)
      call get_ncdim1(fid1, 'time', nnt)

      allocate(plev_in(nz), lon_in(nx), lat_in(ny), time_in(nnt))
      call get_var1_real (fid1, 'lev', nz, plev_in)
      call get_var1_real (fid1, 'lon', nx, lon_in)
      call get_var1_real (fid1, 'lat', ny, lat_in)
      call get_var1_real (fid1, 'time', nnt, time_in)
      call get_att_text (fid1, 'time', 'units', time_unit)

      allocate(pr(nz))
      do k = 1, nz
        pr(k) = plev_in(k) * 0.001
      enddo

      !! range of the data:  south_lat/west_lon/north_lat/east_lon/dlat/dlon :
      hdr(1) = lat_in(1) 
      hdr(2) = lon_in(1)
      hdr(3) = lat_in(ny)
      hdr(4) = lon_in(nx)
      hdr(5) = abs(lon_in(2) - lon_in(1)) 
      hdr(6) = abs(lat_in(2) - lat_in(1))
 
      allocate(t_in(nx,ny,nz,nnt))
      call get_var4_real (fid1, 't', nx, ny, nz, nnt, t_in)

      allocate(h_in(nx,ny,nz,nnt))
      allocate(si_in(nx,ny,nz,nnt))
      call open_ncfile(infile2,fid2)
      call get_var4_real (fid2, 'h', nx, ny, nz, nnt, h_in)
      call get_var4_real (fid2, 'psi', nx, ny, nz, nnt, si_in)
              
      print*, 'before open PV center file' 
      !! the PV center location :
      open(35,file=trim(trackfile),status='old',form='formatted')
      print*, 'after open PV center file'      

      !! allocate all variables which will be used in the loop
      allocate(t(ny,nx,nz), h(ny,nx,nz), si(ny,nx,nz))
      allocate(sib(ny,nx,nz), hb(ny,nx,nz), tb(ny,nx,nz))
      allocate(sip(ny,nx,nz), hp(ny,nx,nz), tp(ny,nx,nz))
      allocate(rhs(ny,nx,nz), tmp(ny,nx,nz))
      allocate(u(ny,nx,nz),v(ny,nx,nz))
      allocate(ub(ny,nx,nz),vb(ny,nx,nz))
      allocate(up(ny,nx,nz),vp(ny,nx,nz))
      !! allocate output variables
      allocate(hb_out(nx,ny,nz,nnt),tb_out(nx,ny,nz,nnt),ub_out(nx,ny,nz,nnt),vb_out(nx,ny,nz,nnt))
      allocate(hp_out(nx,ny,nz,nnt),tp_out(nx,ny,nz,nnt),up_out(nx,ny,nz,nnt),vp_out(nx,ny,nz,nnt))
      allocate(h_out(nx,ny,nz,nnt),t_out(nx,ny,nz,nnt),u_out(nx,ny,nz,nnt),v_out(nx,ny,nz,nnt))

      print*,'*******************Strat to Do loop********************'
      do 214 nt=1,nnt
         print*, 't=',nt

         do k = 1, nz
            do j = ny,1,-1
               do i = 1, nx
                  t(j,i,k)  = t_in(i,ny-j+1,k,nt)
                  h(j,i,k)  = h_in(i,ny-j+1,k,nt)
                  si(j,i,k) = si_in(i,ny-j+1,k,nt)
               enddo
             enddo
         enddo

         ! read in PV center
         print*, 'read in PV center file' 
         !read(35,FMT), hlat, hlon
         !write(*,FMT), hlat, hlo2         
         read(35,215), hlat, hlon
         write(*,215), hlat, hlon
215      format(12x,f5.1,1x,f5.1)
         print*, '######## hlat=',hlat, 'hlon=',hlon, '########'

         do i=1,ny
           do j=1,nx
             do k=1,nz
               h(i,j,k)=h(i,j,k)*gg
               si(i,j,k)=si(i,j,k)*1.e5
             enddo
           enddo
         enddo

         print*, 'call sym to get the symmetric si, h and t'
         call sym(sib,si,hlat,hlon,rmax,nmx,nx,ny,nz,hdr)
         call sym(hb,h,hlat,hlon,rmax,nmx,nx,ny,nz,hdr)
         call sym(tb,t,hlat,hlon,rmax,nmx,nx,ny,nz,hdr)
         do i=1,ny
           do j=1,nx
             do k=1,nz
               hp(i,j,k)=h(i,j,k)-hb(i,j,k)
               sip(i,j,k)=si(i,j,k)-sib(i,j,k)
               tp(i,j,k)=t(i,j,k)-tb(i,j,k)
             enddo
           enddo
         enddo

         print*, 'call right hand side of the non-linear balance eq.'
         call nlbal(rhs,si,sib,nx,ny,nz,hdr)
!cccccccccccccccccccccccccccccccccccc
!c     call relax to solve non-linear balance equation
!c     which we can get the balance height with azimuthal mean psi
!c     note: we use the azimuthal mean h as initial guess filed
!c           and the azimuthal mean h at boundary as b.c.

         thrs=thrs_h
         print*, 'call relaxation to solve perturbation h'
         call relax(hb,tmp,rhs,omegs,thrs,imax,nx,ny,nz,hdr)

         print*, 'call ubal to get the u and v associated psi'
         call ubal(u,v,si,nx,ny,nz,hdr)
         call ubal(ub,vb,sib,nx,ny,nz,hdr)
         call ubal(up,vp,sip,nx,ny,nz,hdr)

         print*, 'Get the basic field of h, t'
         do i=1,ny
           do j=1,nx
             do k=1,nz
               hp(i,j,k)=h(i,j,k)-hb(i,j,k)
             enddo
           enddo
         enddo

         do k = 1, nz
            do j = ny,1,-1
               do i = 1, nx
                 hb_out(i,ny-j+1,k,nt) = hb(j,i,k)
                 tb_out(i,ny-j+1,k,nt) = tb(j,i,k)
                 ub_out(i,ny-j+1,k,nt) = ub(j,i,k)
                 vb_out(i,ny-j+1,k,nt) = vb(j,i,k)

                 hp_out(i,ny-j+1,k,nt) = hp(j,i,k)
                 tp_out(i,ny-j+1,k,nt) = tp(j,i,k)
                 up_out(i,ny-j+1,k,nt) = up(j,i,k)
                 vp_out(i,ny-j+1,k,nt) = vp(j,i,k)

                 h_out(i,ny-j+1,k,nt) = h(j,i,k)
                 t_out(i,ny-j+1,k,nt) = t(j,i,k)
                 u_out(i,ny-j+1,k,nt) = u(j,i,k)
                 v_out(i,ny-j+1,k,nt) = v(j,i,k)
               enddo
             enddo
         enddo

!      write(30,rec=8*(nt-1)+1),(((hb(i,j,k),j=1,nx),i=ny,1,-1),k=1,nz)
!      write(30,rec=8*(nt-1)+2),(((tb(i,j,k),j=1,nx),i=ny,1,-1),k=1,nz)
!      write(30,rec=8*(nt-1)+3),(((ub(i,j,k),j=1,nx),i=ny,1,-1),k=1,nz)
!      write(30,rec=8*(nt-1)+4),(((vb(i,j,k),j=1,nx),i=ny,1,-1),k=1,nz)
!      write(30,rec=8*(nt-1)+5),(((hp(i,j,k),j=1,nx),i=ny,1,-1),k=1,nz)
!      write(30,rec=8*(nt-1)+6),(((tp(i,j,k),j=1,nx),i=ny,1,-1),k=1,nz)
!      write(30,rec=8*(nt-1)+7),(((up(i,j,k),j=1,nx),i=ny,1,-1),k=1,nz)
!      write(30,rec=8*(nt-1)+8),(((vp(i,j,k),j=1,nx),i=ny,1,-1),k=1,nz)
!
!      write(31,rec=4*(nt-1)+1),(((h(i,j,k),j=1,nx),i=ny,1,-1),k=1,nz)
!      write(31,rec=4*(nt-1)+2),(((t(i,j,k),j=1,nx),i=ny,1,-1),k=1,nz)
!      write(31,rec=4*(nt-1)+3),(((u(i,j,k),j=1,nx),i=ny,1,-1),k=1,nz)
!      write(31,rec=4*(nt-1)+4),(((v(i,j,k),j=1,nx),i=ny,1,-1),k=1,nz)

214   continue

      ! Write out the netcdf file
      print*, 'Write out data!!'
      status = nf_create(outfile1, NF_NOCLOBBER, id)
 
      status = nf_def_dim(id, 'lon',  nx, londimid)
      status = nf_def_dim(id, 'lat',  ny, latdimid)
      status = nf_def_dim(id, 'lev',  nz, levdimid)
      status = nf_def_dim(id, 'time', nnt, timedimid)

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
      status = nf_def_var(id, 'hp', NF_REAL, 4, dims , hpid)
      status = nf_def_var(id, 'tp', NF_REAL, 4, dims , tpid)
      status = nf_def_var(id, 'up', NF_REAL, 4, dims , upid)
      status = nf_def_var(id, 'vp', NF_REAL, 4, dims , vpid)
      status = nf_def_var(id, 'hb', NF_REAL, 4, dims , hbid)
      status = nf_def_var(id, 'tb', NF_REAL, 4, dims , tbid)
      status = nf_def_var(id, 'ub', NF_REAL, 4, dims , ubid)
      status = nf_def_var(id, 'vb', NF_REAL, 4, dims , vbid)

      status = nf_enddef(id)

      start(1:1)=(/1/)
      total(1:1)=(/nx/)
      status = nf_put_vara_real(id, lonid, start, total, lon_in)
      start(1:1)=(/1/)
      total(1:1)=(/ny/)
      status = nf_put_vara_real(id, latid, start, total, lat_in)
      start(1:1)=(/1/)
      total(1:1)=(/nz/)
      status = nf_put_vara_real(id, levid, start, total, plev_in)
      start(1:1)=(/1/)
      total(1:1)=(/nnt/)
      status = nf_put_vara_real(id, timeid, start, total, time_in)

      start(1:4)=(/1,1,1,1/)
      total(1:4)=(/nx,ny,nz,nnt/)
      status = nf_put_vara_real(id, hpid, start, total, hp_out)
      status = nf_put_vara_real(id, tpid, start, total, tp_out)
      status = nf_put_vara_real(id, upid, start, total, up_out)
      status = nf_put_vara_real(id, vpid, start, total, vp_out)
      status = nf_put_vara_real(id, hbid, start, total, hb_out)
      status = nf_put_vara_real(id, tbid, start, total, tb_out)
      status = nf_put_vara_real(id, ubid, start, total, ub_out)
      status = nf_put_vara_real(id, vbid, start, total, vb_out)

      status = nf_close(id)

      ! write second file (should consider to combine them)
      status = nf_create(outfile2, NF_NOCLOBBER, id)
 
      status = nf_def_dim(id, 'lon',  nx, londimid)
      status = nf_def_dim(id, 'lat',  ny, latdimid)
      status = nf_def_dim(id, 'lev',  nz, levdimid)
      status = nf_def_dim(id, 'time', nnt, timedimid)

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
      status = nf_def_var(id, 'h', NF_REAL, 4, dims , hid)
      status = nf_def_var(id, 't', NF_REAL, 4, dims , tid)
      status = nf_def_var(id, 'u', NF_REAL, 4, dims , uid)
      status = nf_def_var(id, 'v', NF_REAL, 4, dims , vid)

      status = nf_enddef(id)

      start(1:1)=(/1/)
      total(1:1)=(/nx/)
      status = nf_put_vara_real(id, lonid, start, total, lon_in)
      start(1:1)=(/1/)
      total(1:1)=(/ny/)
      status = nf_put_vara_real(id, latid, start, total, lat_in)
      start(1:1)=(/1/)
      total(1:1)=(/nz/)
      status = nf_put_vara_real(id, levid, start, total, plev_in)
      start(1:1)=(/1/)
      total(1:1)=(/nnt/)
      status = nf_put_vara_real(id, timeid, start, total, time_in)

      start(1:4)=(/1,1,1,1/)
      total(1:4)=(/nx,ny,nz,nnt/)
      status = nf_put_vara_real(id, hid, start, total, h_out)
      status = nf_put_vara_real(id, tid, start, total, t_out)
      status = nf_put_vara_real(id, uid, start, total, u_out)
      status = nf_put_vara_real(id, vid, start, total, v_out)

      status = nf_close(id)

      print*,'End of run'
      end
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     subroutine relax to solve the relaxation
!
!       input rhs(ny,nx,nz) to caculate chi with relaxation method to solve the poission eq.
!
!       lapla(chi(x,y)) + hh(x,y)*chi(x,y) = rhs(x,y)
!
!     for the form
!
!        ( chi(i+1,j,k)+chi(i-1,j,k)+chi(i,j+1,k)+chi(i,j-1,k)-4.*chi(i,j,k) )
!        / (2.*ds*ds) + hh(i,j,k)*chi(i,j,k) = rhs(i,j,k)
!
!     for example to solve stream function 
!
!        lapla(psi(x,y))=vor(x,y)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   note: the data of the first dimension is ny(latitude) from north to south
!         the 2nd dimension is nx(longitude) form west to east
!
!   input rhs, hh, to caculate chi with relaxation method
!
!        hh(ny,nx,nz) : zero in stream function
!        rhs(ny,nx,nz): vorticity in stream function
!
!        omegs: the over-relaxation coefficient
!        imax : the maxium iteration
!        thrs : the threshold
!
!        nx: the number of longitude
!        ny: the number of latitude
!        nz: the number of level
!
!        hdr(1):for initial latitude
!        hdr(2):for initial longitude
!        hdr(3):for final latitude
!        hdr(4):for final longitude
!        hdr(5):for interval of longitude(dx)
!        hdr(6):for interval of latitude(dy)
!
!   output the chi of a poission equation
!
!        chi(ny,nx,nz): psi in stream function
!
!   note that chi(ny,nx,nz) must give the laternal boundary condition first
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine relax(chi,hh,rhs,omegs,thrs,imax,nx,ny,nz,hdr)

      real :: chi(ny,nx,nz), hh(ny,nx,nz)
      real :: rhs(ny,nx,nz)
      real :: hdr(6), a(5)
      real :: sigm, rlat1, rlat2, rlon1, rlon2, dlat, dlon
      real :: omegs, chin, rs, thrs, chisum, chisumold, chi_ratio
      real :: lapchi
      real :: pi, aa, dl, dp, ap, apm, app
      integer :: imax, icount
      logical :: it
      real :: dchisum, dchisumold

!c**** setup domain infomation ***********

      dlat= hdr(6)
      dlon= hdr(5)
      rlat1 = hdr(1)
      rlon1 = Hdr(2)
      rlat2 = Hdr(3)
      rlon2 = Hdr(4)

!c**** set parameter ***********************

      pi=4.*atan(1.)
      aa = 2.e7/pi
      sigm=dlat/dlon

      dl=aa*dlon*pi/180.
      dp=aa*dlat*pi/180.

      do 109 k=1,nz

       do 108 icount=1,imax

         it = .true.
         dchisumold=dchisum
         dchisum=0

         do i=2,ny-1

           ap=cos( pi*(rlat2-float(i-1)*dlat)/180. )
           apm=cos( pi*(rlat2-(float(i)-1.5)*dlat)/180. )
           app=cos( pi*(rlat2-(float(i)-0.5)*dlat)/180. )

!C **       add the laplacian coefficaient

           a(1)=sigm*sigm*apm/ap
           a(2)=1./(ap*ap)
           a(3)=-( 2. + sigm*sigm*ap*(apm+app) )/(ap*ap)
           a(4)=1./(ap*ap)
           a(5)=sigm*sigm*app/ap

           do j=2,nx-1

             lapchi=( 1./(dl*dl) ) *                          \
                   ( a(1)*chi(i-1,j,k) + a(2)*chi(i,j-1,k)+   \
                     a(4)*chi(i,j+1,k) + a(5)*chi(i+1,j,k)+   \
                    (a(3)+hh(i,j,k)*dl*dl )*chi(i,j,k)       )

             rs   = lapchi-rhs(i,j,k)

             chin=chi(i,j,k)

             chi(i,j,k)=chin-omegs*rs*dl*dl / ( a(3)+hh(i,j,k)*dl*dl )

             dchi=chi(i,j,k)-chin

             dchisum=dchisum+dchi/(nx*ny*nz)

             if( abs(dchi) .gt. thrs ) it=.false.

           enddo
         enddo

         if(it) goto 324

107    continue

108    enddo
       print*,'  too many iteration for chi at lev',k
       go to 109

324    print*,'  chi converge at lev',k,', iteration:',icount

109   enddo

      print*,'call relax ok'
      return

      end subroutine relax
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     subroutine nlban to caculate the right hand side of non-linear balance equation that:
!        lap(phi0)= f*lap(psi0)+beta*dpsi0/dy
!                 - 2.*( (ddpsi0/dxdy)**2-(ddpsi0/dxdx)*(ddpsi0/dydy) )
!                 - (   f*lap(psib)+beta*dpsib/dy)
!                     - 2.*( (ddpsib/dxdy)**2-(ddpsib/dxdx)*(ddpsib/dydy)) )
!
!     input si and sib to caculate the rhs of non-linear balance eq.
!
!        si(ny,nx,nz) : the total si
!        sib(ny,nx,nz): the si of basic field
!        u(ny,nx,nz): u-velocity
!        v(ny,nx,nz): v-velocity
!        nx: the number of longitude (west to east)
!        ny: the number of latitude (north to south)
!        nz: the number of level
!        hdr(1):for initial latitude
!        hdr(2):for initial longitude
!        hdr(3):for final latitude
!        hdr(4):for final longitude
!        hdr(5):for interval of longitude(dx)
!        hdr(6):for interval of latitude(dy)
!
!     output is 
!        rhs(ny,nx,nz): rhs of non-linear balance eq.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine nlbal(rhs,si,sib,nx,ny,nz,hdr)

      integer :: nx, ny, nz
      real :: rhs(ny,nx,nz), si(ny,nx,nz), sib(ny,nx,nz)
      real :: hdr(6), dlat, dlon, rlat1, rlat2, rlon1, rlon2
      real :: rlats, rlatn, pi, aa, dphi, dlam, fc
      real :: rhs1, rhs2, rhs3, rhs4, rhs5, rhs6, rhs7, rhs8
      real :: omega
      integer :: i, j, k

!c******************************************
      dlat= hdr(6)
      dlon= hdr(5)
      rlat1 = hdr(1)
      rlon1 = Hdr(2)
      rlat2 = Hdr(3)
      rlon2 = Hdr(4)
!c******************************************
      pi=4.*atan(1.)
      aa = 2.e7/pi
      dphi= dlat*pi/180.
      dlam= dlon*pi/180.
      omega= 2*pi/(24.*60.*60.)

      do k=1,nz
        do i=2,ny-1

!c ****************************
!c         rlat: the real latitude
!c         ap  : the factor of latitude
!c         fc  : Coriolis parameter (vary with latitude)

          rlat= rlat2-float(i-1)*dlat
          rlatn=rlat2-float(i-2)*dlat
          rlats=rlat2-float(i)*dlat

          ap= cos(pi*rlat/180.)
          fc= 2.*omega*sin(pi*rlat/180.)
          beta=2.*omega*( sin(pi*rlatn/180.)-sin(pi*rlats/180.) ) \
               /( 2.*aa*dphi )

          do j=2,nx-1

             rhs1=fc*(   (si(i,j+1,k)-2.*si(i,j,k)+si(i,j-1,k)) \
                         / ((aa*dlam*ap)*(aa*dlam*ap))          \
                      +  (si(i-1,j,k)-2.*si(i,j,k)+si(i+1,j,k)) \
                         / ((aa*dphi)*(aa*dphi))              )

             rhs2=beta*( ( si(i-1,j,k)-si(i+1,j,k) )/(2.*aa*dphi) )

             rhs3=-2.*( ( si(i-1,j+1,k)-si(i-1,j-1,k)-           \
                          si(i+1,j+1,k)+si(i+1,j-1,k) )          \
                         /((2.*aa*dphi)*(2.*aa*dlam*ap)) )**2

             rhs4=2.*(  (si(i,j+1,k)-2.*si(i,j,k)+si(i,j-1,k))   \
                         / ((aa*dlam*ap)*(aa*dlam*ap))       )   \
                    *(  (si(i-1,j,k)-2.*si(i,j,k)+si(i+1,j,k))   \
                         / ((aa*dphi)*(aa*dphi))             )

             rhs5=fc*(   (sib(i,j+1,k)-2.*sib(i,j,k)+sib(i,j-1,k)) \
                         / ((aa*dlam*ap)*(aa*dlam*ap))             \
                      +  (sib(i-1,j,k)-2.*sib(i,j,k)+sib(i+1,j,k)) \
                         / ((aa*dphi)*(aa*dphi))              )

             rhs6=beta*( ( sib(i-1,j,k)-sib(i+1,j,k) )/(2.*aa*dphi) )

             rhs7=-2.*( (sib(i-1,j+1,k)-sib(i-1,j-1,k)-            \
                         sib(i+1,j+1,k)+sib(i+1,j-1,k))            \
                       /((2.*aa*dphi)*(2.*aa*dlam*ap))   )**2

             rhs8=2.*(  (sib(i,j+1,k)-2.*sib(i,j,k)+sib(i,j-1,k))  \
                         / ((aa*dlam*ap)*(aa*dlam*ap))          )  \
                    *(  (sib(i-1,j,k)-2.*sib(i,j,k)+sib(i+1,j,k))  \
                         / ((aa*dphi)*(aa*dphi))                )

!c            rhs(i,j,k)=(rhs1+rhs2+rhs3+rhs4)-(rhs5+rhs6+rhs7+rhs8)
!ccccccccccccccccccccccccccccccc                 cwuhuang 2000,03,31
             rhs(i,j,k)=(rhs5+rhs6+rhs7+rhs8)

          enddo
        enddo
!c       print*,'rhs1,rhs2,rhs3,rhs4,rhs5,rhs6,rhs7,rhs8='
!c       print*,rhs1,rhs2,rhs3,rhs4,rhs5,rhs6,rhs7,rhs8
      enddo

      return

      end subroutine nlbal

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  **  subroutine ubal to  calculated the nondivergent (balanced) wind from the balance psi data
!c
!c   note: the data of the first dimension is ny(latitude) from north to south
!c         the 2nd dimension is nx(longitude) form west to east
!c
!c   input psi to caculate balance u,v
!c        psi(ny,nx,nz): stream function
!c        nx: the number of longitude
!c        ny: the number of latitude
!c        nz: the number of level
!c        hdr(1):for initial latitude
!c        hdr(2):for initial longitude
!c        hdr(3):for final latitude
!c        hdr(4):for final longitude
!c        hdr(5):for interval of longitude(dx)
!c        hdr(6):for interval of latitude(dy)
!c   output is
!c        ur(ny,nx,nz): balance u
!c        vr(ny,nx,nz): balance v
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ubal(ur,vr,psi,nx,ny,nz,hdr)

      integer :: nx, ny, nz
      REAL :: rlat1, rlat2, rlon1, rlon2, dlat, dlon, rlat, ap
      REAL :: psi(Ny,Nx,nz),ur(ny,nx,nz),vr(ny,nx,nz)
      REAL :: hdr(6)  
      real :: pi, aa, dphi, dlam

!c**** setup domain infomation ***********  
 
      dlat= hdr(6)
      dlon= hdr(5)
      rlat1 = hdr(1)
      rlon1 = Hdr(2)
      rlat2 = Hdr(3)                 
      rlon2 = Hdr(4)

!c**** set parameter ***********************

      pi=4.*atan(1.)       
      aa = 2.e7/pi
      dphi= dlat*pi/180.
      dlam= dlon*pi/180.
 
!c****** calculat the wind

       do i=2,ny-1                  

!c ****************************
!c        rlat: the real latitude
!c        ap  : the factor of latitude

         rlat= rlat2-float(i-1)*dlat
         ap= cos(pi*rlat/180.)

!c  ** non-divergent wind
         do k=1,nz
           do j=2,nx-1

             ur(i,j,k)= -(psi(i-1,j,k)-psi(i+1,j,k)) / (2.*aa*dphi)
             vr(i,j,k)=  (psi(i,j+1,k)-psi(i,j-1,k)) / (2.*aa*ap*dlam)
           end do
         end do

       end do

!C  ** extrapolate the wind at boundary
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
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       program sym
!c  **** this program is modified by cwuhuang
!c  **** which calculate the axi-symmetric mean
!c                                             modified by cwuhuang 99'11'30
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     input xi(ny,nx,nz) to caculate symmetric xs(ny,nx,nz)
!c
!c     input : xi(ny,nx,nz)
!c             hlat, hlon : the center of hurricane
!c             rmax : the max. radius to do symmetric
!c             nmx  : # of the directions for a circular
!c
!c     output: xs(ny,nx,nz)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sym(xs,xi,hlat,hlon,rmax,nmx,nx,ny,nz,hdr)
      real :: xs(ny,nx,nz), xi(ny,nx,nz)
      real :: del, tha, hlat, hlon, rmax
      real :: rpr(1000), xms(1000), xmm
      real :: hdr(6), rd(4)
      integer :: nmx, ird, ir
      real :: pi, dlat, dlon
      real :: rlat1, rlat2, rlon1, rlon2
      real :: drr, drmax, pi180, xc, yc, xcp, ycp, rd2
      real :: theta, pp, qq, dist
      real :: xx, yy, ix, iy, dz, c1, c2, r
      integer :: i, j, k, ii

!cccccccccccccccccccccccccccc
!c     set parameter
!c
      pi = 4.*atan(1.)
      dlat = hdr(6)
      dlon = hdr(5)
      rlat1 = hdr(1)
      rlon1 = hdr(2)
      rlat2 = hdr(3)
      rlon2 = hdr(4)
!cccccccccccccccccccccccccccc
!c     we take the drr as 10
!c     so that dr increase 0.1 degree every step
!c
      drr=10.

      pi180 = 4.*atan(1.)/180.
      drmax = dlat*110.e3
      xc = 360. + hlon
      yc = hlat
      xcp = xc*pi180
      ycp = yc*pi180
      rd2=rmax/drmax

      print*, 'nx,ny,dlat,dlon,rlat1,rlon1,rlat2,rlon2='
      print*, nx,ny,dlat,dlon,rlat1,rlon1,rlat2,rlon2
      print*, 'hlat, hlon =', hlat, hlon
!c***********************************************
!c     check the radius of azimuthal mean field
!c     which can not be larger than the distancd 
!c     between center and boundary
!c
      rd (1) = hlon - rlon1
      rd (2) = rlon2 - hlon
      rd (3) = rlat2 - hlat
      rd (4) = hlat - rlat1
      !print*,'rd(i),rd2=',(rd(k),k=1,4),rd2

      do ii = 1,4
         if (rd2.gt.rd(ii))  rd2 = rd(ii)
      end do

      ird = int(drr*rd2/1.)
!c************************************************
!c     loop start

      do 108 k=1,nz
      !print*,'z=',k

      r = 0.
      do ir = 1, ird
        r=r+1./drr

!c**********************************************************
!c       interpolate the grid point to (r, theta) coordinate
!c       and caculate the azimuthal mean of radius r
!c
        xms(ir) = 0.
        do i=1, nmx
          theta = 2.*pi*float(i-1)/float(nmx)
          xx = r*cos(theta) + hlon
          yy = r*sin(theta) + hlat

          ix = int( ( xx-rlon1)/dlon )+1
          iy = int( ( rlat2-yy)/dlat )+1
         !print*,'ix,iy=',ix,iy

          pp = (xx -(rlon1+(ix-1)*dlon) )/dlon
          qq = ((rlat2-(iy-1)*dlat) -yy )/dlat
         !print*,'pp,qq=',pp,qq

          xmm  =   (1.-pp)*(1.-qq)*xi(iy  ,ix  ,k) \
                  +(1.-pp)*    qq *xi(iy+1,ix  ,k) \
                  +    pp *(1.-qq)*xi(iy  ,ix+1,k) \
                  +    pp *    qq *xi(iy+1,ix+1,k) 

          rpr(ir) = r
          xms(ir) = xms(ir)+xmm/float(nmx)

        enddo
      enddo
!     print*,'xms(ir)=',xms
!c****************************************************************
!c     inverse the (r, theta) coordinate to the lat-lon coordinate
!c     and we will get the azimuthal mean at grid point

      do i=1,ny
        do j=1,nx

!c*********
!c  ** get the position of each grid

          del = (360. + rlon1 + float (j-1) * dlon) * pi180
          tha = (rlat2 - float (i-1) * dlat) * pi180

          dist = (del-xcp)*(del-xcp) + (tha-ycp)*(tha-ycp)
          dist = sqrt(dist)/pi180

!c********
!c         do interpolation
!c
          do ir=1, ird-1


            if(dist .ge. rpr(ir) .and. dist .le. rpr(ir+1)) then
              dz=rpr(ir+1) - rpr(ir)
              if( dz .le. 0.0001) print*,'dz = 0.'
              c1 = (rpr(ir+1) -dist)/dz
              c2 = 1. - c1
              xs(i,j,k) = c2*xms(ir+1) + c1*xms(ir)
              goto 324
            endif

          enddo
!c********
!c         do extrapolation
          if (dist .lt. rpr(1) ) then

            c1= rpr(1) - dist
            c2= rpr(2) - rpr(1) 
            xs(i,j,k) = xms(1) + c1/c2*(xms(1)-xms(2))

          elseif (dist .gt. rpr(ird) ) then 
            xs(i,j,k) = xms(ird)
          endif

324     enddo
      enddo

108   continue
      print*, '****hlat, hlon =', hlat, hlon

      return

      end subroutine sym
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

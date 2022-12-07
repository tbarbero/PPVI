!! Compiling on analysrs:
!    module load netcdf/4.6.1
!    module load intel_compilers/oneapi
!    # Full speed run:
!    ifort -O2 steersq.F90 -o steersq.exe -lnetcdf -lnetcdff -lhdf5 -lhdf5_hl -lz -limf -lm
!!! Modified by Jan-Huey Chen in May 2021
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c     Theis program is to caculate the steering of typhoon
!c
!c                                             C.-X. Huang 12,12,99"
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c                                 -- modified by cwudpf, 2008/01/02
      program symbasic

      implicit none
#include <netcdf.inc>
      character*190 :: nmlfile = 'steersq.nml'
      integer, external ::  iargc
      character*300 :: infile, trackfile
      character*300 :: outfile
      character*90 :: time_unit
      integer, parameter :: nmx=48
      integer :: nx, ny, nz, nnt
      real, parameter :: zerowind=0
      real, parameter :: rmax=1000.e5
      real, dimension(:,:,:), allocatable :: u, v, ub, vb, up, vp
      real, dimension(:,:,:), allocatable :: tmp
      real, dimension(:,:), allocatable :: udlm, vdlm, stu, stv, uu, vv
      real, dimension(:), allocatable :: pr, hlat, hlon, ubt, vbt
      real, dimension(:), allocatable :: rlat, rlon
      character*4, dimension(:), allocatable :: tdate
      character*2, dimension(:), allocatable :: dd, hh
      real :: hdr(6), xlat, xlon, pi, ae, dt
      integer :: imax, ideg, nrec, iplow, iptop
      character*4 :: starttime
      logical :: isFindTime

      namelist /nlist/ infile, outfile, trackfile, starttime, ideg
  
      integer :: i, j, k, nt, kk, ii
      integer :: fid, totline, count
      real, dimension(:), allocatable :: lon_in, lat_in, plev_in, time_in
      real, dimension(:,:,:,:), allocatable :: u_in, v_in
      integer :: stas
      character*1 :: lable

      integer :: id, status, levdimid, timedimid
      integer :: levid, timeid
      integer :: uid, vid, hlatid, hlonid, ubtid, vbtid
      integer :: dims(4),start(4), total(4)

! ******************************************************************************
      !! Read in namelist
      if (iargc() > 0 ) then
        call getarg(1,nmlfile)
      endif

      open (11, file=trim(nmlfile), status='old')
      read (11, nml=nlist)
      close (11)

      !! Read in netcdf input files

      call open_ncfile(infile,fid)

      call get_ncdim1(fid, 'lev', nz)
      call get_ncdim1(fid, 'lon', nx)
      call get_ncdim1(fid, 'lat', ny)
      call get_ncdim1(fid, 'time', nnt)

      allocate(plev_in(nz), lon_in(nx), lat_in(ny), time_in(nnt))
      call get_var1_real (fid, 'lev', nz, plev_in)
      call get_var1_real (fid, 'lon', nx, lon_in)
      call get_var1_real (fid, 'lat', ny, lat_in)
      call get_var1_real (fid, 'time', nnt, time_in)
      call get_att_text (fid, 'time', 'units', time_unit)

      allocate(pr(nz))
      do k = 1, nz
        pr(k) = plev_in(k) * 0.001
      enddo

      allocate(u_in(nx,ny,nz,nnt), v_in(nx,ny,nz,nnt))
      call get_var4_real (fid, 'u', nx, ny, nz, nnt, u_in)
      call get_var4_real (fid, 'v', nx, ny, nz, nnt, v_in)

      !! range of the data:  south_lat/west_lon/north_lat/east_lon/dlat/dlon :
      hdr(1) = lat_in(1)
      hdr(2) = lon_in(1)
      hdr(3) = lat_in(ny)
      hdr(4) = lon_in(nx)
      hdr(5) = abs(lon_in(2) - lon_in(1))
      hdr(6) = abs(lat_in(2) - lat_in(1))

      ! set parameter
      pi=4.*atan(1.)
      ae=2.e7/pi

      ! read TC center data:
      open(31,file=trackfile,status='old',form='formatted')
 
      count = 0
      do while(.true.)
        read(31,100,iostat=stas) lable
        if (stas/=0) then
           exit
        else
           if (lable=="Z") then
             count=count+1
           endif
        endif
      enddo
100   format(2x,a1)

      totline = count

      allocate(u(ny,nx,nz), v(ny,nx,nz))
      allocate(ub(ny,nx,nz), vb(ny,nx,nz))
      allocate(up(ny,nx,nz), vp(ny,nx,nz))
      !allocate(udlm(nz,nz), vdlm(nz,nz))
      allocate( stu(0:10,nz), stv(0:10,nz))
      allocate(tmp(ny,nx,nz))
      allocate(uu(nz,nnt), vv(nz,nnt)) 
      allocate(hlat(nnt), hlon(nnt), ubt(nnt), vbt(nnt))
      allocate(tdate(totline), rlat(totline), rlon(totline), dd(totline), hh(totline))

      rewind(31)
      do ii=1,totline
        read(31,120) hh(ii), dd(ii), rlat(ii), rlon(ii)
        tdate(ii) = dd(ii)//hh(ii)
      enddo
120   format(a2, 1x, a2, 7x, f5.1, 1x, f5.1)

  
!      !open(20,file=infile,status='old',form='unformatted', access='direct',recl=4*nx*ny*nz)
!      !open(40,file=outfile,status='new',form='unformatted', access='direct', recl=56+4*(2+nz*(nz+2)))

      kk = totline

      nrec=1
      isFindTime = .false.

      do ii=1,kk
        if (tdate(ii) .eq. starttime) then
          isFindTime = .true.
        endif
        if (isFindTime .and. (nrec .le. nnt)) then
           hlat(nrec)=rlat(ii)
           hlon(nrec)=rlon(ii)
          
          !c caculate the best track
          if(ii .eq. 1) then
            dt=6.*60.*60.
            ubt(nrec) = ae*(rlon(2)-rlon(1))*pi/180.*cos(rlat(1)*pi/180.)/dt
            vbt(nrec) = ae*(rlat(2)-rlat(1))*pi/180./dt
          else if (ii .eq. kk) then
            dt=6.*60.*60.
            ubt(nrec) = ae*(rlon(ii)-rlon(ii-1))*pi/180.*cos(rlat(ii)*pi/180.)/dt
            vbt(nrec) = ae*(rlat(ii)-rlat(ii-1))*pi/180./dt
          else
            dt=12.*60.*60.
            ubt(nrec) = ae*(rlon(ii+1)-rlon(ii-1))*pi/180.*cos(rlat(ii)*pi/180.)/dt
            vbt(nrec) = ae*(rlat(ii+1)-rlat(ii-1))*pi/180./dt
          endif

          print*,'ubt,vbt = ',ubt(nrec),vbt(nrec)
          nrec=nrec+1
        endif
      enddo

      do nt=1,nnt
        !read(20,rec=5*(nt-1)+4)(((u(i,j,k),j=1,nx),i=ny,1,-1),k=1,nz)
        !read(20,rec=5*(nt-1)+5)(((v(i,j,k),j=1,nx),i=ny,1,-1),k=1,nz) 
        do k = 1, nz
           do j = ny,1,-1
              do i = 1, nx
                 u(j,i,k)  = u_in(i,ny-j+1,k,nt)
                 v(j,i,k)  = v_in(i,ny-j+1,k,nt)
              enddo
            enddo
        enddo

        if (nt == 1) then
          print *, 'k = 1', nt, u(:,:,1)
          print *, 'k = 2', nt, u(:,:,2)
        endif  

        call sym(stu,u,rlat(nt),rlon(nt),rmax,nmx,nx,ny,nz,hdr)
        call sym(stv,v,rlat(nt),rlon(nt),rmax,nmx,nx,ny,nz,hdr)

   
        if (nt == 1 .or. nt == 2) then
          print *, 't = ', nt, stu(3,1)
        endif  
!!!
       ! if (hlat(nt) == 0 .and. hlon(nt) == 0) then
       !   stu(ideg,:) = 0
       !   stv(ideg,:) = 0
       ! endif
!!!  
        do i=1,nz
          if (rlat(nt) == 0 .and. rlon(nt) == 0) then
            print *, 'No BT data for storm, setting DLM_SF = 0'
            uu(i,nt)=zerowind 
            vv(i,nt)=zerowind
          else 
            uu(i,nt)=stu(ideg,i)
            vv(i,nt)=stv(ideg,i)
          endif
        enddo

        ! move compute Deep layer mean to python script
        !do iplow=1,nz-1
        !  do iptop=iplow+1, nz
        !    call dlm(udlm(iplow,iptop),uu,pr,iptop,iplow,nz)
        !    call dlm(vdlm(iplow,iptop),vv,pr,iptop,iplow,nz)
        !  if (( iplow .EQ. 2 ).and. ( iptop .EQ. 7 )) then
        !    print*,'u,v at levs 2-7 = ',udlm(iplow,iptop),vdlm(iplow,iptop)
        !  endif
        !  enddo
        !enddo

      enddo ! nt loop

      ! write out ubt(nnt), vbt(nnt), hlat(nnt), hlon(nnt), uu(nz,nnt), vv(nz,nnt)
      print *, 'outfile = ', trim(outfile)
      status = nf_create(outfile, NF_NOCLOBBER, id)

      status = nf_def_dim(id, 'lev',  nz, levdimid)
      status = nf_def_dim(id, 'time', nnt, timedimid)

      dims (1) = levdimid
      status = nf_def_var(id, 'lev', NF_REAL, 1, dims , levid)
      dims (1) = timedimid
      status = nf_def_var(id, 'time', NF_REAL, 1, dims , timeid)

      status = nf_put_att_text(id, levid, 'units', 3,'hPa' )
      status = nf_put_att_text(id, timeid, 'units', len(trim(time_unit)),trim(time_unit) )

      dims(1) = levdimid
      dims(2) = timedimid
      status = nf_def_var(id, 'u_ave', NF_REAL, 2, dims , uid)
      status = nf_def_var(id, 'v_ave', NF_REAL, 2, dims , vid)

      dims(1) = timedimid
      status = nf_def_var(id, 'ubt', NF_REAL, 1, dims , ubtid)
      status = nf_def_var(id, 'vbt', NF_REAL, 1, dims , vbtid)
      status = nf_def_var(id, 'lat_TC_center', NF_REAL, 1, dims , hlatid)
      status = nf_def_var(id, 'lon_TC_center', NF_REAL, 1, dims , hlonid)

      status = nf_enddef(id)

      start(1:1)=(/1/)
      total(1:1)=(/nz/)
      status = nf_put_vara_real(id, levid, start, total, plev_in)
      start(1:1)=(/1/)
      total(1:1)=(/nnt/)
      status = nf_put_vara_real(id, timeid, start, total, time_in)

      start(1:2)=(/1,1/)
      total(1:2)=(/nz,nnt/)
      status = nf_put_vara_real(id, uid, start, total, uu)
      status = nf_put_vara_real(id, vid, start, total, vv)

      start(1:1)=(/1/)
      total(1:1)=(/nnt/)
      status = nf_put_vara_real(id, ubtid, start, total, ubt)
      status = nf_put_vara_real(id, vbtid, start, total, vbt)
      status = nf_put_vara_real(id, hlatid, start, total, hlat)
      status = nf_put_vara_real(id, hlonid, start, total, hlon)

      status = nf_close(id)


      do nt = 1, nnt
        print *, nt, uu(:,nt), vv(:,nt)
      enddo

      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       program steering
!c  **** this program is modified by cwuhuang
!c  **** which calculate the steering axi-symmetric mean
!c
!c                                             modified by cwuhuang 99'11'30
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     input xi(ny,nx,nz) to caculate symmetric steering flow
!c
!c     input : xi(ny,nx,nz)
!c             hlat, hlon : the center of hurricane
!c             rmax : the max. radius to do symmetric
!c             nmx  : # of the directions for a circular
!c
!c     output: stu(0:10, nz)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sym(stu,xi,hlat,hlon,rmax,nmx,nx,ny,nz,hdr)
      real :: xi(ny,nx,nz)
      real :: del, tha, hlat, hlon, rmax
      real :: rpr(1000), xms(1000), xmm
      real :: hdr(6), rd(4), stu(0:10,nz)
      integer nmx, ird, ir
      real :: pi, dlat, dlon
      real :: rlat1, rlat2, rlon1, rlon2
      real :: drr, drmax, pi180, xc, yc, xcp, ycp, rd2
      real :: theta, pp, qq, dist
      real :: xx, yy, ix, iy, dz, c1, c2, r
      integer :: i, j, k, ii
   
      !c set parameter
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
      drr=10.

      pi180 = 4.*atan(1.)/180.
      drmax = dlat*110.e3
      xc = 360. + hlon
      yc = hlat
      xcp = xc*pi180
      ycp = yc*pi180
      rd2=rmax/drmax
!
!c***********************************************
!c     check the radius of azimuthal mean field
!c     which can not be larger than the distancd 
!c     between center and boundary
!c
      rd (1) = hlon - rlon1
      rd (2) = rlon2 - hlon
      rd (3) = rlat2 - hlat
      rd (4) = hlat - rlat1
      !c  print*,'rd(i),rd2=',(rd(k),k=1,4),rd2

      do ii = 1,4
         if (rd2.gt.rd(ii))  rd2 = rd(ii)
      end do

      ird = int(drr*rd2/1.)
!c************************************************
!c     loop start

      do 108 k=1,nz
!c     print*,'z=',k

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

!cccccccccccccccccccccccccccccccccc
!c     caculat the steering flow

      stu(0,k)=0.
      stu(0,k)=2.*xms(1)-xms(2)

      do i=1,10
        stu(i,k)=0.
        nmean=0
        do j=1,i*drr
          nmean=nmean+j
        enddo
        do j=1,i*drr
          stu(i,k)=stu(i,k)+j*xms(j)/float(nmean)
        enddo
      enddo
!c     print*,'steering at lev ',k
!c     print*,(stu(i,k)-xms(i*drr),i=0,5)

108   continue

       return

       end subroutine sym
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Subroutine dlm to calculate the deep layer mean
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dlm(ud,u,pr,iptop,iplow,nz)

      real ud, u(nz), pr(nz)
      integer nz, iptop, iplow

      dp=pr(iptop)-pr(iplow)
      ud=0.

      do i=iplow, iptop-1
        ud=ud + ( u(i+1)+u(i) )/2.*( pr(i+1)-pr(i) )/dp
      enddo
!c     print*,'ud=',ud
      return
      end subroutine dlm
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

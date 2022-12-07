#!/bin/csh -f
source $MODULESHOME/init/tcsh
source parameters.csh
module unload netcdf
module load netcdf

set outfile_dir = '/archive/twb/pv_inversion/2017/NAtl/'${num}'/'${model}'' 
set infile_dir = '/archive/twb/pv_inversion/2017/NAtl/'${num}'/'${model}''

foreach date ($ft)
  foreach ini_hour ( 00 12 )
   
    # GFS/SHiELD_gfsIC/SHiELD_ifsIC
    set infile1 = $infile_dir/${date}${ini_hour}/${date}_${ini_hour}Z_htuv.nc
    set infile2 = $infile_dir/${date}${ini_hour}/${date}_${ini_hour}Z_pv.nc
    set trackfile = $infile_dir/${date}${ini_hour}/center.txt
    set outfile1 = $outfile_dir/${date}${ini_hour}/${date}_${ini_hour}Z_symhtuv.nc
    set outfile2 = $outfile_dir/${date}${ini_hour}/${date}_${ini_hour}Z_htuv2.nc
    
cat >! symhtuv.nml << EOF
&nlist     ! -*- F90 -*-
   infile1 = '$infile1'
   infile2 = '$infile2'
   trackfile = '$trackfile'
   outfile1 = '$outfile1'
   outfile2 = '$outfile2'
   nmx = 48        ! The number of direction to do azimuthal mean
   rmax = 2000000. ! The radius(m) to caculate symmetric u,v
   imax = 3000     ! Max. iteration for psi and h (nonlinear balance eq.)
   omegs = 1.6     ! Over-relaxation parameter
   thrs_h = 1.     ! Thres hold of h
&end
EOF

#/home/Jan-Huey.Chen/NGGPS/TC/pv_inversion/new/symhtuv.exe
./symhtuv.exe
  end
end 

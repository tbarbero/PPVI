#!/bin/csh -f
source $MODULESHOME/init/tcsh
source parameters.csh
module unload netcdf
module load netcdf

set outfile_dir = '/archive/twb/pv_inversion/2017/NAtl/'${num}'/'${model}''
set infile_dir = '/archive/twb/pv_inversion/2017/NAtl/'${num}'/'${model}''

foreach date ($ft)
  foreach ini_hour ( 00 12 )
    
    set infile = $infile_dir/${date}${ini_hour}/${date}_${ini_hour}Z_symhtuv.nc
    set outfile = $outfile_dir/${date}${ini_hour}/${date}_${ini_hour}Z_pvb.nc
    
cat >! pvb.nml << EOF
&nlist     ! -*- F90 -*-
   infile = '$infile'
   outfile = '$outfile'
   imax = 3000
   omegs = 1.85
   thrs = 1.e5
&end
EOF

#/home/Jan-Huey.Chen/NGGPS/TC/pv_inversion/new/pvb.exe
./pvb.exe
  end
end 

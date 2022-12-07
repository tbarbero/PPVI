#!/bin/csh -f
source $MODULESHOME/init/tcsh
source parameters.csh
module unload netcdf
module load netcdf

set outfile_dir = '/archive/twb/pv_inversion/2017/NAtl/'${num}'/'${model}''
set infile_dir = '/archive/twb/pv_inversion/2017/NAtl/'${num}'/'${model}''

#foreach date ( `seq 20170914 1 20170930` )
foreach date ($ft)
  foreach ini_hour ( 00 12 )
    
    # GFS/SHiELD_gfsIC/SHiELD_ifsIC
    set outfile = $outfile_dir/${date}${ini_hour}/${date}_${ini_hour}Z_pv2.nc
    set infile = ${infile_dir}/${date}${ini_hour}/${date}_${ini_hour}Z_htuv2.nc
    
    # IFS 
   # set outfile = $outfile_dir/${date}_${ini_hour}/${date}_${ini_hour}Z_pv2.nc
   # set infile = ${infile_dir}/${date}_${ini_hour}/${date}_${ini_hour}Z_htuv2.nc

cat >! pv2.nml << EOF
&nlist     ! -*- F90 -*-
   outfile = '$outfile'
   infile = '$infile'
   imax = 3000
   omegs = 1.85
   thrs = 1.e5
&end
EOF

#/home/Jan-Huey.Chen/NGGPS/TC/pv_inversion/new/pv2.exe
./pv2.exe
  end
end 

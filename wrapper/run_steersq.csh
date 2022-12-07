#!/bin/csh -f
source $MODULESHOME/init/tcsh
source parameters.csh
module unload netcdf
module load netcdf

# start time for 1700 file was '1706'
# just every file date (ddhh) + 06 hours?

set outfile_dir = '/archive/twb/pv_inversion/2017/NAtl/'${num}'/'${model}''
set infile_dir = '/archive/twb/pv_inversion/2017/NAtl/'${num}'/'${model}''

#foreach date ( 20170924 20170925 )

foreach date ($ft)
  set d = `echo $date | cut -c 7-`

  foreach ini_hour ( 00 12 )
    if ($ini_hour == 00) then
      set h = 06
    else if ($ini_hour == 12) then
      set h = 18 
    endif
    set starttime = $d$h

    foreach tag ('ne' 'se' 'sw' 'nw' 'bhigh' 'chigh' 'sqinvph')
   


    set infile = $infile_dir/${date}${ini_hour}/${date}_${ini_hour}Z_${tag}.nc
    #set infile = $infile_dir/${date}${ini_hour}/${date}_${ini_hour}Z_sqinvph.nc
    set trackfile = $infile_dir/${date}${ini_hour}/center.txt

    # ideg is the radius to compute axisymmetric mean steering flow in degree, from 3 to 10
    set ideg = 3 
    set outfile = ${outfile_dir}/${date}${ini_hour}/${date}_${ini_hour}Z_steering_${ideg}degree_${tag}.nc
    
  cat >! steersq.nml << EOF
&nlist     ! -*- F90 -*-
   infile = '$infile'
   trackfile = '$trackfile'
   outfile = '$outfile'
   starttime = '$starttime'
   ideg = $ideg  ! radius to compute axisymmetric mean steering flow in degree, from 3 to 10
&end
EOF

#/home/Jan-Huey.Chen/NGGPS/TC/pv_inversion/new/steersq.exe
./steersq.exe 
    end
    end
end
 

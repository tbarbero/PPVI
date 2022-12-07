#!/bin/csh -f

set echo
source parameters.csh
source $MODULESHOME/init/tcsh
module unload nco
module load nco
module load netcdf

foreach hr ($num)
#foreach hr (`seq -f "%02g" 04 04`)
  set root_dir = '/archive/twb/pv_inversion/2017/NAtl/'${hr}'/IFS'
  cd ${root_dir}
  set dates = `ls` 
  echo $dates

  foreach date ($dates)
    cd ${root_dir}/$date
    foreach file (`ls *pv.nc`)
      echo $hr $file
      ncks -O -d time,1,40 $file $file 
    end
    foreach file (`ls *htuv.nc`)
      echo $hr $file
      ncks -O -d time,1,40 $file $file
    end

    
  end
end


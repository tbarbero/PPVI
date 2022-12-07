#!/bin/csh -f
source $MODULESHOME/init/tcsh
source parameters.csh
module unload netcdf
module load netcdf

set outfile_dir = '/archive/twb/pv_inversion/2017/NAtl/'${num}'/'${model}''
set infile_dir = '/archive/twb/pv_inversion/2017/NAtl/'${num}'/'${model}''

#set system = 'nw'

foreach date ($ft)
#foreach date (20170916)
  foreach ini_hour ( 00 12 )
    foreach system ('ne' 'se' 'sw' 'nw' 'chigh' 'bhigh') 
    #set system = 'chigh'

#    set system = ${sys}
   # GFS/SHiELD models
    set outfile = ${outfile_dir}/${date}${ini_hour}/${date}_${ini_hour}Z_${system}.nc
    set infile1 = ${infile_dir}/${date}${ini_hour}/${date}_${ini_hour}Z_pvb.nc
    set infile2 = ${infile_dir}/${date}${ini_hour}/${date}_${ini_hour}Z_pv2.nc

    # IFS
#    set outfile = $outfile_path/${date}_${ini_hour}/${date}_${ini_hour}Z_${system}.nc
#    set infile1 = ${infile_dir}/${date}_${ini_hour}/${date}_${ini_hour}Z_pvb.nc
#    set infile2 = ${infile_dir}/${date}_${ini_hour}/${date}_${ini_hour}Z_pv2.nc


    set ntimes = '1'
    foreach mask_file (`ls -v ${infile_dir}/${date}${ini_hour}/mask/${system}*`)
      set nt = `echo $mask_file | cut -f2 -d"." `
      if ($nt != '1') then 
        set ntimes = `echo $ntimes $nt | tr  -s ' ' '\n' | tr -s '\n' ',' `
      endif
    end

    set maskfiles = `ls -v ${infile_dir}/${date}${ini_hour}/mask/${system}* | tr ' ' '\n'  |  sed 's/^\(.*\)$/"\1"/' | tr '\n' ',' | tr ' ' ','`
    

cat >! sqinv.nml << EOF
&nlist     ! -*- F90 -*-
   outfile = '$outfile'
   infile1 = '$infile1'
   infile2 = '$infile2'
   mask = .true.      ! if mask = T, isturg will be set up as "3" in the code
   maskfiles = $maskfiles
   ntimes = $ntimes
   OMEGAS = 1.8        ! relaxation parameter 
   OMEGAH = 1.6        ! relaxation parameter
   PRT = 0.5           ! underrelaxation parameter
   THRSH = 0.1         ! CONVERGENCE THRESHOLD
   TSCAL = 1.          ! SCALE FACTOR FOR THTA
   QSCAL = 1.          ! SCALE FOR Q
   IMAP = 1            ! 2 FOR XY, 1 FOR SPHERICAL
   INLIN = 1           ! LINEARIZATION PARM (0=THROW OUT NONLINEAR TERMS, 1=EQUAL PARTITIONING)
   IQD = 0             ! QPRM DEPENDS ON ANOTHER FILE? (O = NO)
   NOUT = 1            ! NUMBERS OF INVERSIONS
   IBC = 1             ! BOUNDARY CONDITIONS (0=HOMOGENEOUS)
     
&end
EOF

#/home/Jan-Huey.Chen/NGGPS/TC/pv_inversion/new/sqinvph.exe 
./sqinvph.exe
        end
    end
end     

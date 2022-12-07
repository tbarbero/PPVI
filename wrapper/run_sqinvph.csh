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
    set outfile = ${outfile_dir}/${date}${ini_hour}/${date}_${ini_hour}Z_sqinvph.nc
    set infile1 = ${infile_dir}/${date}${ini_hour}/${date}_${ini_hour}Z_pvb.nc
    set infile2 = ${infile_dir}/${date}${ini_hour}/${date}_${ini_hour}Z_pv2.nc
    
    # IFS
#    set outfile = $outfile_path/${date}_${ini_hour}/${date}_${ini_hour}Z_sqinvph.nc
#    set infile1 = ${infile_dir}/${date}_${ini_hour}/${date}_${ini_hour}Z_pvb.nc
#    set infile2 = ${infile_dir}/${date}_${ini_hour}/${date}_${ini_hour}Z_pv2.nc    

cat >! sqinv.nml << EOF
&nlist     ! -*- F90 -*-
   outfile = '$outfile'
   infile1 = '$infile1'
   infile2 = '$infile2'
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
   itsurg = 0          ! do surgery or not: 0: no; 1: surge inner; 2: surge outer
   IBC = 1             ! BOUNDARY CONDITIONS (0=HOMOGENEOUS)
     
&end
EOF

#/home/Jan-Huey.Chen/NGGPS/TC/pv_inversion/new/sqinvph.exe
./sqinvph.exe
  end
end 

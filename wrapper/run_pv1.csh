#!/bin/csh -f
source $MODULESHOME/init/tcsh
source /home/Tyler.Barbero/scripts/pv_inversion_tool/parameters.csh
module unload netcdf
module load netcdf
set echo
echo $num
echo $ft
echo $model
echo $geopotential

if ( ! -d /archive/twb/pv_inversion/2017/NAtl/$num ) then
    mkdir /archive/twb/pv_inversion/2017/NAtl/$num
endif  

set outfile_path = '/archive/Tyler.Barbero/pv_inversion/2017/NAtl/'${num}'/'${model}''
set infile_dir = '/archive/Tyler.Barbero/'${model}''

if ( $model == "SHiELD_ifsIC" ) then
    echo 'model: ${model}'
    set x_name = 'grid_xt'
    set y_name = 'grid_yt'
    set level_name = 'plev'
    set time_name = 'time'
    set h_name = 'h_plev'
    set t_name = 't_plev'
    set u_name = 'u_plev'
    set v_name = 'v_plev'  
else if ( ${model} == "SHiELD_gfsIC" ) then
    echo 'model: ${model}'
    set x_name = 'grid_xt'
    set y_name = 'grid_yt'
    set level_name = 'plev'
    set time_name = 'time'
    set h_name = 'h_plev'
    set t_name = 't_plev'
    set u_name = 'u_plev'
    set v_name = 'v_plev'    
else if ( ${model} == "GFS" ) then
    echo 'model: ${model}'
    set x_name = 'lon'
    set y_name = 'lat'
    set level_name = 'lev'
    set time_name = 'time'
    set h_name = 'hgt'
    set t_name = 'temp'
    set u_name = 'uwind'
    set v_name = 'vwind'
else if ( ${model}  == "GFS_analysis" ) then
    echo 'model: ${model}'
    set infile_dir = '/archive/Tyler.Barbero/GFS/analysis'
    set x_name = 'lon'
    set y_name = 'lat'
    set level_name = 'lev'
    set time_name = 'time'
    set h_name = 'hgt'
    set t_name = 'temp'
    set u_name = 'uwind'
    set v_name = 'vwind'
else if ( ${model} == "IFS" ) then
    echo 'model: ${model}'
    set x_name = 'longitude'
    set y_name = 'latitude'
    set level_name = 'level'
    set time_name = 'time'
    set h_name = 'z'
    set t_name = 't'
    set u_name = 'u'
    set v_name = 'v'
else if ( ${model} == "IFS_analysis" ) then
    echo 'model: ${model}'
    set infile_dir = '/archive/Tyler.Barbero/IFS/analysis'
    set x_name = 'longitude'
    set y_name = 'latitude'
    set level_name = 'level'
    set time_name = 'time'
    set h_name = 'z'
    set t_name = 't'
    set u_name = 'u'
    set v_name = 'v'
else
    break
endif 


foreach date ($ft)
  foreach ini_hour ( 00 12 )
    
    # make directory
    mkdir -p ${outfile_path}/${date}${ini_hour}
    
    # GFS: lon, lat, lev, time, hgt, temp, uwind, vwind   
    # IFS: longitude, latitude, level, time, z, t, u, v
    # IFS_analysis:
#    set outfile1 = ${outfile_path}/${date}_${ini_hour}/${date}_${ini_hour}Z_htuv.nc
#    set outfile2 = ${outfile_path}/${date}_${ini_hour}/${date}_${ini_hour}Z_pv.nc
#    set h_file = ${infile_dir}/1_AN0_${date}_10day_analysis.z.nc
#    set t_file = ${infile_dir}/1_AN0_${date}_10day_analysis.t.nc
#    set u_file = ${infile_dir}/1_AN0_${date}_10day_analysis.u.nc
#    set v_file = ${infile_dir}/1_AN0_${date}_10day_analysis.v.nc 
    # SHiELD_gfsIC: grid_xt, grid_yt, plev, time, h_plev, t_plev, u_plev, v_plev
    # SHiELD_ifsIC: grid_xt, grid_yt, plev, time, h_plev, t_plev, u_plev, v_plev
    if ($model == "SHiELD_ifsIC") then 
        set outfile1 = ${outfile_path}/${date}${ini_hour}/${date}_${ini_hour}Z_htuv.nc
        set outfile2 = ${outfile_path}/${date}${ini_hour}/${date}_${ini_hour}Z_pv.nc
        set h_file = ${infile_dir}/${date}_${ini_hour}Z_SHiELD_ifsIC_4xdaily_h_1deg.nc
        set t_file = ${infile_dir}/${date}_${ini_hour}Z_SHiELD_ifsIC_4xdaily_t_1deg.nc
        set u_file = ${infile_dir}/${date}_${ini_hour}Z_SHiELD_ifsIC_4xdaily_u_1deg.nc
        set v_file = ${infile_dir}/${date}_${ini_hour}Z_SHiELD_ifsIC_4xdaily_v_1deg.nc
    else if ($model == "SHiELD_gfsIC") then
        set outfile1 = ${outfile_path}/${date}${ini_hour}/${date}_${ini_hour}Z_htuv.nc
        set outfile2 = ${outfile_path}/${date}${ini_hour}/${date}_${ini_hour}Z_pv.nc
        set h_file = ${infile_dir}/${date}_${ini_hour}Z_SHiELD_gfsIC_4xdaily_h_1deg.nc
        set t_file = ${infile_dir}/${date}_${ini_hour}Z_SHiELD_gfsIC_4xdaily_t_1deg.nc
        set u_file = ${infile_dir}/${date}_${ini_hour}Z_SHiELD_gfsIC_4xdaily_u_1deg.nc
        set v_file = ${infile_dir}/${date}_${ini_hour}Z_SHiELD_gfsIC_4xdaily_v_1deg.nc
        # these are ifs files, but i accidently named them 'SHiELD_gfsIC' 
    else if ($model == "GFS") then
        set outfile1 = $outfile_path/${date}${ini_hour}/${date}_${ini_hour}Z_htuv.nc
        set outfile2 = $outfile_path/${date}${ini_hour}/${date}_${ini_hour}Z_pv.nc
        set h_file = ${infile_dir}/${date}${ini_hour}/${date}_${ini_hour}Z_h_360x180.nc
        set t_file = ${infile_dir}/${date}${ini_hour}/${date}_${ini_hour}Z_t_360x180.nc
        set u_file = ${infile_dir}/${date}${ini_hour}/${date}_${ini_hour}Z_u_360x180.nc
        set v_file = ${infile_dir}/${date}${ini_hour}/${date}_${ini_hour}Z_v_360x180.nc
    else if ($model == "GFS_analysis") then
        set outfile1 = $outfile_path/${date}${ini_hour}/${date}_${ini_hour}Z_htuv.nc
        set outfile2 = $outfile_path/${date}${ini_hour}/${date}_${ini_hour}Z_pv.nc
        set h_file = ${infile_dir}/${date}_${ini_hour}Z00_h_360x180.10day.nc
        set t_file = ${infile_dir}/${date}_${ini_hour}Z00_t_360x180.10day.nc
        set u_file = ${infile_dir}/${date}_${ini_hour}Z00_u_360x180.10day.nc
        set v_file = ${infile_dir}/${date}_${ini_hour}Z00_v_360x180.10day.nc
    else if ($model == 'IFS') then
        set outfile1 = ${outfile_path}/${date}${ini_hour}/${date}_${ini_hour}Z_htuv.nc
        set outfile2 = ${outfile_path}/${date}${ini_hour}/${date}_${ini_hour}Z_pv.nc
        if ($ini_hour == 00) then
          #set ini_hour = `echo ${ini_hour} | sed 's/.$//'`
          set ini_hour = 0
        else if ($ini_hour == 12) then
          set ini_hour = `echo ${ini_hour}00`
        endif
        set h_file = ${infile_dir}/1_FC0_${date}_${ini_hour}.z.nc
        set t_file = ${infile_dir}/1_FC0_${date}_${ini_hour}.t.nc
        set u_file = ${infile_dir}/1_FC0_${date}_${ini_hour}.u.nc
        set v_file = ${infile_dir}/1_FC0_${date}_${ini_hour}.v.nc
    else if ($model == 'IFS_analysis') then
        set outfile1 = ${outfile_path}/${date}${ini_hour}/${date}_${ini_hour}Z_htuv.nc
        set outfile2 = ${outfile_path}/${date}${ini_hour}/${date}_${ini_hour}Z_pv.nc
        if ($ini_hour == 00) then
          set ini_hour = 0
        else if ($ini_hour == 12) then
          set ini_hour = `echo ${ini_hour}00`
        endif
        set h_file = ${infile_dir}/1_AN0_${date}_${ini_hour}.z.10day.nc
        set t_file = ${infile_dir}/1_AN0_${date}_${ini_hour}.t.10day.nc
        set u_file = ${infile_dir}/1_AN0_${date}_${ini_hour}.u.10day.nc
        set v_file = ${infile_dir}/1_AN0_${date}_${ini_hour}.v.10day.nc  
    endif

cat >! pv1.nml << EOF
&nlist     ! -*- F90 -*-
   outfile1 = '$outfile1'
   outfile2 = '$outfile2'
   h_infile = '$h_file'
   t_infile = '$t_file'
   u_infile = '$u_file'
   v_infile = '$v_file'
   plev = 1000., 925., 850., 700., 500., 400., 300., 250., 200., 150., 100., 70.
   dx = 1.
   dy = 1.
   west = 250.
   east = 355.
   south = 0.
   north = 55.
   imax = 3000
   omegs = 1.85
   thrs = 1.e5
   ! change parameters for each model below
   geopotential = $geopotential
   x_name = '$x_name'
   y_name = '$y_name'
   level_name = '$level_name'
   time_name = '$time_name' 
   h_name = '$h_name'
   t_name = '$t_name'
   u_name = '$u_name'
   v_name = '$v_name'
&end
EOF

./pv1.exe
#/home/Jan-Huey.Chen/NGGPS/TC/pv_inversion/new/pv1.exe

  end
end 

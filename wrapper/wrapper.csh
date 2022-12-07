#!/bin/csh -f 

# turn on source with python scripts
source parameters.csh
set models = ( "GFS_analysis" )
#set numbers = ( "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" )
#set numbers = ( "02" "03" "05" "06" "07" "08" "09" "10" "11" "12" "13")
set numbers = ("04")
# need 04

#set numbers = ("01")
foreach model ($models)
  foreach num ($numbers)
    ### get forecast time list for each storm
    if ( $model == 'IFS_analysis' || $model == 'IFS' ) then
      set line = `head -1 /archive/twb/pv_inversion/2017/NAtl/best_track/IFS/2017.NAtl.${num}.txt`
    else if ( $model == 'GFS_analysis' || $model == 'GFS' ) then
      set line = `head -1 /archive/twb/pv_inversion/2017/NAtl/best_track/GFS/2017.NAtl.${num}.txt`
    else
      set line = `head -1 /archive/twb/pv_inversion/2017/NAtl/best_track/${model}/2017.NAtl.${num}.txt`
    endif
    ### set geopotential in nml 
    if ( $model == 'GFS_analysis' || $model == 'SHiELD_gfsIC' || $model == 'SHiELD_ifsIC' || $model == 'GFS' ) then
      set geopotential = .false.
    else
      set geopotential = .true.
    endif
    
    ### set forecast times (ft)
    set t1 = `echo $line | cut -d" " -f9 | sed 's/..$//'`
    set t2 = `echo $line | sed 's/.* //'| sed 's/..$//' `
    set t11 = `echo $t1 | head -c 6 | tail -c 1`
    set t22 = `echo $t2 | head -c 6 | tail -c 1`
    set yyyymm1 = `echo $t1 | head -c 6`
    set yyyymm2 = `echo $t2 | head -c 6`
    
    if ($t11 != $t22) then
      echo 'not equal'
      set interval1 = `seq $t1 1 ${yyyymm1}31`
      set interval2 = `seq ${yyyymm2}01 1 $t2`
      set ft = `echo $interval1 $interval2`
    else if ($t11 == $t22) then
      set ft = `seq $t1 1 $t2`
    endif
    echo 'forecast_times' $ft 
    ### set forecast times (ft) 

    ###### TESTING ########
    # get first run time of each storm to test
    #set ft = `echo $ft | head -n1 | cut -d " " -f1`
    #echo $ft
    ####### TESTING ########
    #set ft = "20170821" 
    
    ### send variables to parameter file for PPVI tool to read in
    echo 'set model = "'$model'"' > parameters.csh # > important to rewrite file
    echo 'set num = "'$num'"' >> parameters.csh
    echo 'set ft = "'$ft'"' >> parameters.csh 
    echo 'set geopotential = '${geopotential}'' >> parameters.csh
    
    ./run_pv1.csh
    #./correct_IFS_time.csh
     python find_TC_position.py "$model" "$num" 
     ./run_symhtuv.csh
     ./run_pv2.csh
     ./run_pvb.csh
     ./run_sqinvph.csh
     python mask_region.py "$model" "$num" "$ft"
     ./run_sqinvmask.csh
     ./run_steersq.csh

  end
end

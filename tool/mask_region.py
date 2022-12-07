#!/usr/bin/env python
import numpy as np
import xarray as xr
from array import array
import sys
import os

########################################
def tic():
    import time
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    import time
    if 'startTime_for_tictoc' in globals():
        print("Elapsed time is " + str(time.time() - startTime_for_tictoc) + " seconds.")
    else:
        print ("Toc: start time not set")

def read_csh_variables():
    model = sys.argv[1]
    num = sys.argv[2]
    ft = sys.argv[3]
    return model, num, ft 

def reset_grid():
    lon_2d, lat_2d = np.meshgrid(ds.lon.values, ds.lat.values)
    mask_region = np.zeros_like(lon_2d, dtype="bool")
    return mask_region, lon_2d, lat_2d

def read_bt_file(bt_file_path):
    file1 = open(bt_file_path)
    lines = file1.readlines()
    n = len(lines)
    lat_bt = np.zeros(n)
    lon_bt = np.zeros(n)

    for i in np.arange(0,n):
        lat_bt[i] = lines[i].split()[1]
        lon_bt[i] = lines[i].split()[2]
    return lat_bt, lon_bt

def read_pv_file(pv_file_path):
    ds = xr.open_dataset(pv_file_path)
    return ds

def save_binary_file(region_name,t):
    # inputs:
    # region name (str) i.e., "chigh"
    # t: time index (int)
    
    # convert numpy array to binary flat GrAds file
    binary_array = ds[region_name+"_mask"][t,:,:,:].data.flatten(order='C').astype('float32')

    # save file
    filename = maskdir + region_name + '.' + str(t+1)
    with open(filename, 'wb') as wf:
        binary_array.tofile(wf)
    wf.close()
    print(filename, ' saved!')

#####################################
### Main Function

model,num,ft = read_csh_variables()
ft = ft.split()
root_dir = '/archive/twb/pv_inversion/2017/NAtl/'+num+'/'+model+'/'

for day in ft:
    for ini_hour in ['00','12']:
        maskdir = root_dir+day+ini_hour+'/mask/'
        print('Saving file',day+ini_hour)
        
        if os.path.isdir(maskdir):
            print('directory already exists')
        else:
            os.mkdir(maskdir)
        
        try:
            # get lat and lon best track
            lat_bt, lon_bt = read_bt_file(root_dir+day+ini_hour+'/center.txt')
            # read in q' file once
            ds = read_pv_file(root_dir+day+ini_hour+'/'+day+'_'+ini_hour+'Z_sqinvph.nc')
        except FileNotFoundError:
            print('No data for: ',day+ini_hour)
            continue

        #########
        tic()
        levs_ignore = [0,7,8,9,10,11]
        ### calculate high pressure mask [static boundaries]
        # BHIGH
        mask_region, lon_2d, lat_2d = reset_grid()
        mask_region[:,ds.lon.values >= 280.5] = True # slab
        mask_region[(lon_2d - 280.5) <= (lat_2d - 30.5)] = False # triangle
        ds["bhigh_mask"] = ds["pv"].copy()
        ds["bhigh_mask"].data = -1 * np.ones_like(ds["bhigh_mask"].data)
        ds["bhigh_mask"].data[(ds["pv"].data < 0) & mask_region[None,None,:,:]] = 1
        ds["bhigh_mask"][:,levs_ignore,:,:] = -1
        
        # CHIGH
        mask_region, lon_2d, lat_2d = reset_grid()
        mask_region[:, ds.lon.values <= 280.5] = True # slab
        mask_region[(lon_2d - 280.5) <= (lat_2d - 30.5)] = True # triangle
        ds["chigh_mask"] = ds["pv"].copy()
        ds["chigh_mask"].data = -1 * np.ones_like(ds["chigh_mask"].data)
        ds["chigh_mask"].data[(ds["pv"].data < 0) & mask_region[None, None, :, :]] = 1
        ds["chigh_mask"][:,levs_ignore,:,:] = -1
        
        for t in np.arange(0,len(lon_bt)):
            if (lat_bt[t]==0) & (lon_bt[t] == 0):
                # set mask at all levs to zero
                print('No best track position for time t=', str(t))
                ds["ne_mask"] = ds["pv"].copy()
                ds["se_mask"] = ds["pv"].copy()
                ds["sw_mask"] = ds["pv"].copy()
                ds["nw_mask"] = ds["pv"].copy()
                ds["ne_mask"].data[t,:,:,:] = -1
                ds["se_mask"].data[t,:,:,:] = -1
                ds["sw_mask"].data[t,:,:,:] = -1
                ds["nw_mask"].data[t,:,:,:] = -1 
                save_binary_file(region_name="ne",t=t)
                save_binary_file(region_name="se",t=t)
                save_binary_file(region_name="sw",t=t)
                save_binary_file(region_name="nw",t=t)
                ### HIGH PRESSURES
                # BHIGH
                save_binary_file(region_name="bhigh",t=t) 
                # CHIGH
                save_binary_file(region_name="chigh",t=t)     
            else:
                               
                ### HIGH PRESSURES
                # BHIGH
                save_binary_file(region_name="bhigh",t=t)
        
                # CHIGH
                save_binary_file(region_name="chigh",t=t)

                ### QUADRANTS
                # NE 
                mask_region, lon_2d, lat_2d = reset_grid()
                mask_region[:,:] = True
                mask_region[ds.lat.values <= lat_bt[t],:] = False
                mask_region[:,ds.lon.values <= lon_bt[t]] = False
                ds["ne_mask"] = ds["pv"].copy()
                ds["ne_mask"].data = -1 * np.ones_like(ds["ne_mask"].data)
                ds["ne_mask"].data[(ds["pv"].data > 0) & mask_region[None,None,:,:]] = 1
                ds["ne_mask"][:,levs_ignore,:,:] = -1
                save_binary_file(region_name="ne",t=t)
        
                # SE 
                mask_region, lon_2d, lat_2d = reset_grid()
                mask_region[:,:] = True
                mask_region[ds.lat.values >= lat_bt[t],:] = False
                mask_region[:,ds.lon.values <= lon_bt[t]] = False
                ds["se_mask"] = ds["pv"].copy()
                ds["se_mask"].data = -1 * np.ones_like(ds["se_mask"].data)
                ds["se_mask"].data[(ds["pv"].data > 0) & mask_region[None,None,:,:]] = 1
                ds["se_mask"][:,levs_ignore,:,:] = -1
        
                save_binary_file(region_name="se",t=t)
        
                # SW
                mask_region, lon_2d, lat_2d = reset_grid()
                mask_region[:,:] = True
                mask_region[ds.lat.values >= lat_bt[t],:] = False
                mask_region[:,ds.lon.values >= lon_bt[t]] = False
                ds["sw_mask"] = ds["pv"].copy()
                ds["sw_mask"].data = -1 * np.ones_like(ds["sw_mask"].data)
                ds["sw_mask"].data[(ds["pv"].data > 0) & mask_region[None,None,:,:]] = 1
                ds["sw_mask"][:,levs_ignore,:,:] = -1
        
                save_binary_file(region_name="sw",t=t)
        
                # NW
                mask_region, lon_2d, lat_2d = reset_grid()
                mask_region[:,:] = True
                mask_region[ds.lat.values <= lat_bt[t],:] = False
                mask_region[:,ds.lon.values >= lon_bt[t]] = False
                ds["nw_mask"] = ds["pv"].copy()
                ds["nw_mask"].data = -1 * np.ones_like(ds["nw_mask"].data)
                ds["nw_mask"].data[(ds["pv"].data > 0) & mask_region[None,None,:,:]] = 1
                ds["nw_mask"][:,levs_ignore,:,:] = -1
        
                save_binary_file(region_name="nw",t=t)
        toc() 

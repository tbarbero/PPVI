#!/usr/bin/env python
import sys
import numpy as np
import xarray as xr
import scipy
from scipy.interpolate import RectBivariateSpline
import os.path
import datetime


# Functions
def vert_avg(array,axis):
    # array: [time,levels,lat,lon]
    # axis = choose levels index
    # returns array of averaged values [time,lat,lon]
    
    array_interp = np.mean(array[:,int(np.array(np.where(ds.lev==925))):int(np.array(np.where(ds.lev==300))),:,:],axis=1);    
    return array_interp

def lines_to_da(ft_lines):
    # https://xarray.pydata.org/en/stable/user-guide/combining.html
    array = np.zeros([5,len(ft_lines)])
    time = []
    intensity_list = []

    for i in np.arange(0,len(ft_lines)):

        line = ft_lines[i]
        data = line.split()

        t = int(data[0])
        time.append(t)
        lon = float(data[1])
        lat = float(data[2])
        pressure = float(data[3])
        wspd = float(data[4])
        wdir = float(data[5])
        intensity = str(data[6])
        intensity_list.append(intensity)

        array[0,i] = lon
        array[1,i] = lat
        array[2,i] = pressure
        array[3,i] = wspd
        array[4,i] = wdir

    lon = xr.DataArray(name='longitude',data=array[0,:],dims=['time'],coords=dict(time=time))
    lat = xr.DataArray(name='latitude',data=array[1,:],dims=['time'],coords=dict(time=time))
    pressure = xr.DataArray(name='pressure',data=array[2,:],dims=['time'],coords=dict(time=time))
    wspd = xr.DataArray(name='wspd',data=array[3,:],dims=['time'],coords=dict(time=time))
    wdir = xr.DataArray(name='wdir',data=array[4,:],dims=['time'],coords=dict(time=time))
    inten = xr.DataArray(name='intensity',data=intensity_list,dims=['time'],coords=dict(time=time))
    ds_grid = [[lon,lat,pressure,wspd,wdir,inten]]
    da = xr.combine_nested(ds_grid, concat_dim=["time", None])

    return da

def return_pos_ind(ds,lat,lon):
    # objective: return index of lat/lon position of boundaries
    # ds = xr.dataset
    # ds.lon, ds.lat are coordinates
    # lon, lat from forecast BT files to get lat/lon 5x5 gride
    

#    # IF IFS DATA
    if (float(ds.lat[0])).is_integer():
        try:
            lon_lower = np.where(ds.lon == int(lon)-5)[0][0]
            lon_upper = np.where(ds.lon == int(lon)+5)[0][0]
            lat_lower = np.where(ds.lat == int(lat)-5)[0][0]
            lat_upper = np.where(ds.lat == int(lat)+5)[0][0]

        except IndexError:
            try:
                lon_lower = np.where(ds.lon == int(lon)-4)[0][0]
                lon_upper = np.where(ds.lon == int(lon)+4)[0][0]
                lat_lower = np.where(ds.lat == int(lat)-4)[0][0]
                lat_upper = np.where(ds.lat == int(lat)+4)[0][0]

            except IndexError:
                try:
                    lon_lower = np.where(ds.lon == int(lon)-3)[0][0]
                    lon_upper = np.where(ds.lon == int(lon)+3)[0][0]
                    lat_lower = np.where(ds.lat == int(lat)-3)[0][0]
                    lat_upper = np.where(ds.lat == int(lat)+3)[0][0]

                except IndexError:
                    try:
                        lon_lower = np.where(ds.lon == int(lon)-2)[0][0]
                        lon_upper = np.where(ds.lon == int(lon)+2)[0][0]
                        lat_lower = np.where(ds.lat == int(lat)-2)[0][0]
                        lat_upper = np.where(ds.lat == int(lat)+2)[0][0]

                    except IndexError:
                        try:
                            lon_lower = np.where(ds.lon == int(lon)-1)[0][0]
                            lon_upper = np.where(ds.lon == int(lon)+1)[0][0]
                            lat_lower = np.where(ds.lat == int(lat)-1)[0][0]
                            lat_upper = np.where(ds.lat == int(lat)+1)[0][0]

                        except IndexError:
                            try:
                                print('Using original forecast model BT as storm center')
                                lon_lower = 0
                                lon_upper = 0
                                lat_lower = 0
                                lat_upper = 0
                                lat = lat
                                lon = lon
                            except IndexError:
                                print('Exiting script')
    
    # ELSE SHIELD_ifsIC/SHiELD_gfsIC/GFS
    else:
        try:
            lon_lower = np.where(ds.lon == int(lon)-5.5)[0][0]
            lon_upper = np.where(ds.lon == int(lon)+5.5)[0][0]
            lat_lower = np.where(ds.lat == int(lat)-5.5)[0][0]
            lat_upper = np.where(ds.lat == int(lat)+5.5)[0][0]

        except IndexError:
            try:
                lon_lower = np.where(ds.lon == int(lon)-4.5)[0][0]
                lon_upper = np.where(ds.lon == int(lon)+4.5)[0][0]
                lat_lower = np.where(ds.lat == int(lat)-4.5)[0][0]
                lat_upper = np.where(ds.lat == int(lat)+4.5)[0][0]

            except IndexError:
                try:
                    lon_lower = np.where(ds.lon == int(lon)-3.5)[0][0]
                    lon_upper = np.where(ds.lon == int(lon)+3.5)[0][0]
                    lat_lower = np.where(ds.lat == int(lat)-3.5)[0][0]
                    lat_upper = np.where(ds.lat == int(lat)+3.5)[0][0]

                except IndexError:
                    try:
                        lon_lower = np.where(ds.lon == int(lon)-2.5)[0][0]
                        lon_upper = np.where(ds.lon == int(lon)+2.5)[0][0]
                        lat_lower = np.where(ds.lat == int(lat)-2.5)[0][0]
                        lat_upper = np.where(ds.lat == int(lat)+2.5)[0][0]

                    except IndexError:
                        try:
                            lon_lower = np.where(ds.lon == int(lon)-1.5)[0][0]
                            lon_upper = np.where(ds.lon == int(lon)+1.5)[0][0]
                            lat_lower = np.where(ds.lat == int(lat)-1.5)[0][0]
                            lat_upper = np.where(ds.lat == int(lat)+1.5)[0][0]

                        except IndexError:
                            try:
                                print('Using original forecast model BT as storm center')
                                lon_lower = 0
                                lon_upper = 0
                                lat_lower = 0
                                lat_upper = 0
                                lat = lat
                                lon = lon
                            except IndexError:
                                print('Exiting script')

    # returns indices of lat/lon boundaries
    return lon_lower, lon_upper, lat_lower, lat_upper, lat, lon 

def find_min_ind(array):
    ind = np.where(array == np.max(array))
    return ind

def read_csh_variables():
    # model and storm number
    model = sys.argv[1]
    num = sys.argv[2]
    return model, num


# Main Function
# read in pushed variables from csh script 
model,num = read_csh_variables()
print('model:',model,'storm_num:',num)

# open TC center forecast data
f = open('/archive/twb/pv_inversion/2017/NAtl/best_track/'+model+'/2017.NAtl.'+num+'.txt','r')
textfile = f.readlines()
f.close()
obs_times = textfile[0]
initialization_times = obs_times.split()[8:]
#print(initialization_times)
lines = textfile[1:]
ind = 0
ind_list = []
it_list = []
'forecast_time_list:', 
for line in lines:
    if 'overlap_BT_for' in line:
        ind_list.append(ind)
        tmp = line.split()[2]
        it_list.append(tmp)
    if 'forecast' in line:
        ind_list.append(ind)
        tmp = line.split()[2]
        it_list.append(tmp)
        
    ind+=1
#test
#initialization_times = ['2017091500']
#print(it_list)

rootdir = '/archive/twb/pv_inversion/2017/NAtl/'+num+'/'+model+'/'

for t in initialization_times:
    print('Opening PV file for forecast time: ',t) 
    date = t[0:8]
    hr = t[8:]
    try:
        ds = xr.open_dataset(rootdir+t+'/'+date+'_'+hr+'Z_pv.nc')    
    except FileNotFoundError:
	print('File not found, continue to next file')
        continue

    pv_interp = vert_avg(array=ds.pv,axis=1)

    # create 40 forecast times starting at each init time to point to times in model forecast BT data
    forecast_time_list = []
#    m = 0
    ini_hour = t 
    while len(forecast_time_list) < 40: 
        dt_obj = datetime.datetime.strptime(ini_hour,'%Y%m%d%H')
        dt_obj = dt_obj + datetime.timedelta(hours=6)
        dt_str = datetime.datetime.strftime(dt_obj, '%Y%m%d%H')
        ini_hour = dt_str
        #print(ini_hour)
        forecast_time_list.append(dt_str)
        #i+=1       
    try:
        init = it_list.index(t) + 1 # init forecast ind 
        fin = it_list.index(t) + 2 # end of forecast ind
        it_lines = lines[ind_list[init]:ind_list[fin]][1:] # ignore header
    except IndexError:
        init = it_list.index(t) + 1 # init forecast ind
        fin = it_list.index(t) + 2 # end of forecast ind
        it_lines = lines[ind_list[init]:][1:] # ignore header 
    except ValueError:
        continue
    #except ValueError:
    #    print('end of BT location file')
    #    break

#    if t and 'forecast' not in lines:
 #       breakm#    try: 
#        lines[ind_list[fin]];
#    except IndexError:
#        print('end of BT location file')
#        break
#      
    print(*it_lines)   
    da = lines_to_da(it_lines)
    datimes = np.array(da.time)
    datimes = [str(datime) for datime in datimes]
    
    # open file to write
    abs_fpath = os.path.join(rootdir+t+'/center.txt')
    #print(abs_fpath)
#    try:
#        os.remove(abs_fpath)
#    except FileNotFoundError:
#        print('Center.txt already exists!')
#        pass
#        
    file1 = open(abs_fpath,'w')
    # compare each forecast time to time in datimes
    ii=0


#    print('forecast_time_list:',forecast_time_list)
    for ft in forecast_time_list:

        ft_str = datetime.datetime.strptime(ft,'%Y%m%d%H')
        ft_str = datetime.datetime.strftime(ft_str,'%HZ%d%^b%Y')

        try:
            ind = datimes.index(ft)
        except:
            # if datime not in forecast time, skip and set as 0
            lat = 0
            lon = 0
            linestr = str(ft_str)+' '+str(lat)+' '+str(lon)+'\n'
            file1.write(linestr)
            print('Writing location data to file:')
            print(linestr)
            continue # skip to next ft
  
        lat_bt = float(da.latitude[ind])
        lon_bt = float(da.longitude[ind])
      
        # get actual max PV value with 5x5 box around lat_bt and lon_bt
        lon_lower, lon_upper, lat_lower, lat_upper, lat, lon  = return_pos_ind(ds,lat_bt,lon_bt)
        
        if lon_lower==lon_upper==lat_lower==lat_upper==0:
            # if boundaries == 0, use original forecast BT location
            linestr = str(ft_str)+' '+str(round(lat_bt,1))+' '+str(round(lon_bt,1))+'\n'
            file1.write(linestr)
            print('Writing location data to file:')
            print(linestr)
            continue # skip to next ft    
        else:

            tmp = pv_interp.isel({'lon':np.arange(lon_lower,lon_upper,1),'lat':np.arange(lat_lower,lat_upper,1)})
            # create higher-res grid for spline inter
            lat_fine = np.linspace(ds.lat[lat_lower],ds.lat[lat_upper],len(ds.lat[lat_lower:lat_upper])*10)
            lon_fine = np.linspace(ds.lon[lon_lower],ds.lon[lon_upper],len(ds.lon[lon_lower:lon_upper])*10)

           # create spline object using original coordinates (x,y) and gridded data (z)
            try:
                spl = RectBivariateSpline(ds.lat[lat_lower:lat_upper],ds.lon[lon_lower:lon_upper],tmp[ii,:,:])
            except:
                # cannot create spline, too little data (position close to edge domain), use forecast bt position 
                linestr = str(ft_str)+' '+str(round(lat_bt,1))+' '+str(round(lon_bt,1))+'\n'
                file1.write(linestr)
                print('Writing location data to file:')
                print(linestr)
                continue # skip to next ft
  
            
            pv_fine = spl(lat_fine,lon_fine)
            
            # find lat/lon position of min pv val
            pos = find_min_ind(pv_fine);
            lat_ind = pos[0][0]
            lon_ind = pos[1][0]
            lon = round(lon_fine[lon_ind],1)
            lat = round(lat_fine[lat_ind],1)

        # write lat/lon to text file
        linestr = str(ft_str)+' '+str(lat)+' '+str(lon)+'\n'
        file1.write(linestr)
        print('Writing location data to file:')
        print(linestr)
        
        ii+=1
    
        
        # iterate through pv time file
    file1.close()

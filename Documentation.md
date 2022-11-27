The Piecewise Potential Vorticity Inverison (PPVI) is computed through a series of scripts, in order below:
1. run_pv1.csh
2. run_pv1.csh
3. find_TC_center.py
4. run_symhtuv.csh
5. run_pv2.csh
6. run_pvb.csh
7. run_sqinvph.csh
8. mask_region.py
9. run_sqinvmask.csh
10. run_steersq.csh 

Step 1. run_pv1.csh <-- pv1.exe <-- pv1.F90
* Input files:  global h, t, u, v data at 1 by 1 degree
* Output files: yyyymmdd_hhZ_htuv.nc  yyyymmdd_hhZ_pv.nc
    * yyyymmdd_hhZ_htuv.nc: regional and smoothed h,t,u,v
	* h: height [m]
	* t: temperature [K]
	* u,v: wind components [m/s]
    * yyyymmdd_hhZ_pv.nc
        * pv: potential vorticity [0.01PVU; 1e-2 x 1e-6 m^2 1/s K 1/kg]
        * h: adjusted h [m]
        * psi: streamfunction [m^2/s]
        * u,v: non-divergent wind [m/s]
        * avor: absolute vorticity [1/s]

* *Objective:*
    * Crop global data grid into basin-wide region given by the namelist variables 'west', 'east', 'south', 'north' in run_pv1.csh 
    * 'east' and 'west' $\epsilon$ (0,360)
    * 'south' and 'north' $\epsilon$ (-90,90)
    * Apply 9-point local smoothing function within new domain
    * Save regional, smoothed h,t,u,v in 'yyyymmdd_hhZ_htuv.nc'
    * Adjust geopoential height, h, to make sure static stability is positive everywhere (dh/dp>0)
    * Compute relative vorticity from u,v, then compute the streamfunction from relative vorticity 
    * Compute Ertel Potential Vorticity (PV)
    * Compute non-divergent wind from streamfunction 
    * Save regional variables: Ertel PV, adjusted h, non-divergent u and v, psi, absolute vorticity in 'yyyymmdd_hhZ_pv.nc'
 
#```
#``` 

### Math


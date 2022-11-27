# Documentation
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

Moreover, the csh/python scripts above, is wrapped in .csh scripts, to run seamlessly all together

### Step 1. run_pv1.csh

Executable: pv1.exe, pv1.F90

Input files:  global h,t,u,v data; 1deg x 1deg

Output files: *yyyymmdd_hhZ_htuv.nc*, *yyyymmdd_hhZ_pv.nc*

* *yyymmdd_hhZ_htuv.nc*: regional and smoothed h,t,u,v
    * h: height [m] 
    * t: temperature [K] 
    * u,v: wind components [m/s]
* *yyyymmdd_hhZ_pv.nc*
    * pv: potential vorticity [0.01PVU; 1e-2 x 1e-6 m^2 1/s K 1/kg]
    * h: adjusted h [m] 
    * psi: streamfunction [m^2/s]
    * u,v: non-divergent wind [m/s]
    * avor: absolute vorticity [1/s]

Objective:
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

### Step 2. find_TC_center.py

Input files: *2017.NAtl.num.txt*, *yyyymmdd_hhZ_pv.nc*

Output files: *center.txt*

1. 06Z26AUG2017 lat.x lon.x
2. 12Z26AUG2017 lat.x lon.x

...

40. 00Z05SEP2017 lat.x lon.x

Objective:
* For all initialized runs of a storm, compute the TC storm center out to 240Z every 6 hours
    * Load in regional PV file, vertically average PV from 925hPa to 300hPa [forecast_times,lat,lon]
    * Grab initial TC storm center guess (lat.xx, lon.xx) and make a 5deg x 5deg box around center on PV data
    * Compute rectangular bivariate spline interpolation for PV box out to the hundreths place (.xx) for a more accurate center
    * Find (lat,lon) coordinate of maximum PV value and define as TC center
    * Write out 6-hr TC center forecast out to 240Z for each initialized run in 'center.txt'
 
### Step 3. run_symhtuv.csh

Executable: symhtuv.exe, symhtuv1.F90

Input files: *yyyymmdd_hhZ_htuv.nc*, *yyyymmdd_hhZ_pv.nc*, *center.txt*

Output files: *yyyymmdd_hhZ_htuv2.nc*,  *yyyymmdd_hhZ_symhtuv.nc*
* *yyyymmdd_hhZ_htuv2.nc*
    * h: geopotential height (h in *yyyymmdd_hhZ_htuv.nc* multiplied by gravity coefficient 9.8 m/s^2), [m^2/s^2]
    * t: temperature (identical to t in *yyyymmdd_hhZ_htuv.nc*) [K]
    * u,v: non-divergent wind (similar to u,v, in *yyyymmdd_hhZ_pv.nc*)

* *yyyymmdd_hhZ_symhtuv.nc*
    * hb: basic geopotential height 
    * tb: basic temperature
    * ub, vb: non-divergent wind associated with the basic stream function
    * hp: perturbation geopotential height 
    * tp: perturbation temperature 
    * up, vb: non-divergent wind associated with the perturbation stream function

Objective:
* Read in adjusted height and streamfunction from *yyyymmdd_hhZ_pv.nc*, and temperature from  *yyyymmdd_hhZ_htuv.nc*
* Compute basic fields: symmetric height(h), streamfunction(psi), and temperature(t) using TC center
* Compute perturbation fields: subtract basic fields from total fields of h, psi, and t
* Compute basic and perturbation non-divergent winds from associated stream functions, individually
 
### Step 4. run_pv2.csh

Executable: pv2.exe, pv2.F90

Input files: *yyyymmdd_hhZ_htuv2.nc*

Output files: *yyyymmdd_hhZ_pv2.nc*

* *yyyymmdd_hhZ_pv2.nc*
    * pv: potential vorticity [0.01PVU; 1e-2 x 1e-6 m^2 1/s K 1/kg]
    * h: adjusted height [m] 
    * psi: stream function  [m^2/s]
    * u,v: non-divergent wind [m/s]
    * avor: absolute vorticity [1/s]
    * Note: these are the **TOTAL** fields

Objective:
* Repeat computing parts in pv1.F90, but use the non-divergent wind

### Step 5. run_pvb.csh

Executable: pvb.exe, pvb.F90 

Input files: *yyyymmdd_hhZ_symhtuv.nc*

Output files: *yyyymmdd_hhZ_pvb.nc*

* *yyyymmdd_hhZ_pvb.nc*
    * pv: potential vorticity [0.01PVU; 1e-2 x 1e-6 m^2 1/s K 1/kg]
    * h: adjusted height [m] 
    * psi: stream function  [m^2/s]
    * u,v: non-divergent wind [m/s]
    * avor: absolute vorticity [1/s]
    * Note: these are the **BASIC** fields

Objective:
* Repeat computing parts in pv1.F90, but use the basic fields 

### Step 6. run_sqinvph.csh

Executable: sqinvph.exe, sqinvph.F90

Input files:  *yyyymmdd_hhZ_pvb.nc*, *yyyyyyymmdd_hhZ_pv2.nc*

Output files: *yyyymmdd_hhZ_sqinvph.nc*

* *yyyymmdd_hhZ_sqinvph.nc*
    * pv: potential vorticity [0.01PVU; 1e-2 x 1e-6 m^2 1/s K 1/kg]
    * h: adjusted height [m]
    * psi: stream function  [m^2/s]
    * u,v: non-divergent wind [m/s]

Objective:
* Compute perturbation by subtracting basic fields from total field
* Compute the non-divergent wind associated with the total PV perturbation 

### Step 7: mask_region.py

Input files: *yyyymmdd_hhZ_sqinvph.nc*

Output files: Binary flat masked PV regions i.e.,
* 'bhigh.1','bhigh.2',...,'bhigh.40'
* 'chigh.1','chigh.2',...,'chigh.40'
* 'nw.1', 'nw.2', ...,'nw.40'
* 'ne.1', 'ne.2', ...,'ne.40'
* 'se.1', 'se.2', ...,'se.40'
* 'sw.1', 'sw.2', ...,'sw.40'

Obective:
* Read in total PV perturbation field and TC center location
* Cut up the total PV perturbation field into different regions:
    * Note: the decomposition of the total PV perturbation field is subjective to what you systems you want to investigate
    * Cut the negative PV perturbation into static Bermuda High (BH) and Continental High (CH)
    * Cut the positive PV pertubation into four ordinal (NW,NE,SE,SW) qudrants relative to the storm center at each forecast time
* Save the cut PV perturbation for each region for each forecast time as a binary flat file
### Reminder: put in picture here about decomposition of PV reigons

### Step 8. run_sqinvmask.csh

Executable: sqinvph.exe, sqinvph.F90

Input files:  *yyyymmdd_hhZ_pvb.nc*,  *yyyyyyymmdd_hhZ_pv2.nc*, mask files for one targeted system (i.e., BH)

Output files: *yyyymmdd_hhZ_$system.nc*

* *yyyymmdd_hhZ_$system.nc*
    * pv: potential vorticity [0.01PVU; 1e-2 x 1e-6 m^2 1/s K 1/kg]
    * h: adjusted height [m]
    * psi: stream function  [m^2/s]
    * u,v: non-divergent wind [m/s]

Objective: Compute the non-divergent wind associated with with the piecewise pv perturbation for a targeted system

### Step 9. run_steersq.csh

Executable: *steersq.exe*, *steersq.F90*

Input files:
* *yyyymmdd_hhZ_sqinvph.nc* (u, v associated with total pv perturbation) _**or**_ 
* *yyyymmdd_hhZ_$system.nc* (u, v associate with pv perturbations of particular system)

* *center.txt* TC track file


Output files: *yyyymmdd_hhZ_steering_3degree_$system*

* *yyyymmdd_hhZ_steering_3degree_$system*

    * u_ave,v_ave (forecast_time,level): average wind in circle of 3-degree radius based on the TC center at each level
    * ubt,vbt (time): TC movement vector at each time step, computed based on the TC locations in the best track data
    * lat_TC_center,lon_TC_center (time): TC center at each time step. They are the same as in center.txt. 

Objective: Calculate the average u-,v-wind in a circle of 3deg radius around the TC center for each level at each forecast time

### Step 10. run_dlmsf

Input files:

Output files:

Objective: Calculate the deep layer mean steering flow (DLMSF) from 925hPa to 300hPa, weighted on the width of each layer.

<br>
<br>
<br>

# Equations and Notation

Ertels Potential Vorticity

$$\\Large q = \\frac{g\\kappa\\pi}{p}[(f+\\nabla^2\\Phi)\\frac{\\partial^2\\Phi}{\\partial\\pi^2}-\\frac{1}{a^2cos^2\\phi}\\frac{\\partial^2\\Phi}{\\partial\\lambda\\partial\\pi}\\frac{\\partial^2\\Phi}{\\partial\\lambda\\partial\\pi}-\\frac{1}{a^2}\\frac{\\partial^2\\Psi}{\\partial\\phi\\partial\\pi}\\frac{\\partial^2\\Phi}{\\partial\\phi\\partial\\pi}]$<br>$$

Charney's Nonlinear Balanced Equation

$$\\Large\\nabla^2 \\Phi = \\nabla \\cdot ( f \\nabla \\Psi ) + \\frac{2}{a^4cos^2\\phi} \\frac{\\partial(\\partial\\Phi/\\partial\\lambda,\\partial\\Phi/\\partial\\phi)}{\\partial(\\lambda,\\phi)}$$

Nondivergent Wind relation

$$\\Large V = \\hat{k} \\times \\nabla\\Psi$$

Hydrostatic Balance

$$\\Large\\theta = -\\frac{\\partial\\phi}{\\partial\\pi}$$

Relative vorticity and streamfunction relation (Elliptical 2nd order PDE)

$$\\Large\\zeta = \\frac{\\partial^2\\Psi}{\\partial x^2} + \\frac{\\partial ^2\\Psi}{\\partial y^2}$$

**Boundary Conditions**

Lateral boundadry of $\\theta$ and $\\psi$ and $\\theta$ on the upper and lower boundaries

$$\\Large \\psi' =  \\psi - \\overline{\\psi}$$

$$\\Large\\phi' =  \\phi - \\hat{\\phi}$$

$$\\Large q' = q - \\hat{q}$$

* $x'$ denotes perturbation fields,
* $x$ denotes total fields,
* $\\hat{x}$ denotes basic fields

* $Global h,t,u,v$


* $q$: potential vorticity<br>,
* $\\Phi$: geopotential height<br>,
* $\\Psi$: streamfunction<br>,
* $f$: coriolis parameter<br>,
* $a$: Earth's radius<br>,
* $\\kappa$:$\\frac{R_d}{C_p}$<br>",
* $\\pi$: $C_p(p/p_o)^\\kappa$ - vertical coordinate<br>",
* $\\lambda$: latitude<br>",
* $\\phi$: longitude"

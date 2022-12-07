#!/usr/bin/env python

import numpy as np

def dlm(u,plev):
    '''
     u: 3-degree averaged u-, v-component of wind at each pressure level fora single time (levels)
     plev: pressure level array (levels)
     ud: u or v DLMSF component (t)
    '''
    # calculate dlmsf for single forecast time from 925 to 300 hPa 
    ud = 0 
    ind1 = int(np.where(ds.lev==925)[0])
    ind2 = int(np.where(ds.lev==300)[0])
    dp = int(plev[ind2])-int(plev[ind1])

    for i in range(ind1,ind2):
        ud = ud + ((u[i+1]+u[i])/2)*np.abs((plev[i+1]-plev[i])/(dp))
    
    return ud

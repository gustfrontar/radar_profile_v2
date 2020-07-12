#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 19:33:33 2020

@author: juan
"""

import glob
import scipy.io as sio
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import rvd_read as rr 

def grid_profile( ref , z , z_min=0.0 , z_max=15000.0 , delta_z = 500.0 , undef = -32.0 ) :

   zp = np.arange( z_min + delta_z / 2.0 , z_max + delta_z / 2.0 , delta_z )
   zp_min = np.arange( z_min , z_max , delta_z )
   zp_max = zp_min + delta_z
   nbin = np.size( zp_min ) 

   meanp=np.zeros(nbin) + np.nan
   stdp = np.copy(meanp)
   minp = np.copy(meanp)
   maxp = np.copy(meanp)
   nump = np.zeros(nbin)   

   z[z==undef] = np.nan

   for ii in range( nbin ) :
       my_mask = ( z <= zp_max[ii] ) & ( z >= zp_min[ii] ) & ( ~ np.isnan( z ) )
       tmp_num = np.sum( my_mask )
       if  tmp_num > 30 :
          meanp[ii] = np.mean( ref[my_mask] )
          stdp[ii]  = np.std ( ref[my_mask] )
          minp[ii]  = np.min ( ref[my_mask] )
          maxp[ii]  = np.max ( ref[my_mask] )
          nump[ii]  = tmp_num

   return zp , meanp , stdp , minp , maxp , nump

def extract_profile_data( filename , radius , lonp , latp , lonradar , latradar , altradar ) :
    
    radar = rr.rvd_read( filename , lonradar , latradar , altradar )

    # Buscamos los puntos que estan en el cilindro
    dlon =  np.cos(radar.gate_latitude['data']*np.pi/180.0)*( radar.gate_longitude['data'] - lonp )    
    dlat =  ( radar.gate_latitude['data'] - latp )
    distancia = np.sqrt( ( dlon * 111000.0 )**2 + ( dlat * 111000.0 )**2 )

    mascara = distancia <= radius 

    elevations= np.tile( radar.elevation['data'] , (radar.fields['reflectivity']['data'].shape[1] ,1 )).transpose()

    ref = radar.fields['reflectivity']['data'][ mascara ]

    elev = elevations[ mascara ]

    alt  = radar.gate_z['data'][ mascara ]

    return ref , alt , elev 



DataPath = './data/'
FilelistRVD = np.sort(glob.glob(DataPath+'*.z.rvd'))
nfiles = len( FilelistRVD )


latradar = -34.787778
lonradar = -58.536667
altradar = 30          # m

# Coordenadas del centro del area de interes
latp = -34.590187
lonp = -58.483944

# Radio del cilindro en metros
radius = 10000.0
    

# Obtenemos los puntos que estan dnetro del cilindro y obtenemos los valores de 
# reflectividad, altura, elevacion  para dichos puntos.

for ifile,my_file in enumerate( FilelistRVD ) :
   print('Reading ',my_file)

   [ref , alt , elev ] = extract_profile_data( my_file , radius , lonp , latp , lonradar , latradar , altradar ) 

   [zp , meanp , stdp , minp , maxp , nump] = grid_profile( ref , alt )
   
   if ifile == 0 :
      nz = np.size( zp ) 
      z_th_profile = np.zeros( (nz , nfiles) )
      meanref_th_profile = np.zeros( (nz , nfiles) )
      stdref_th_profile  = np.zeros( (nz , nfiles) )
      maxref_th_profile  = np.zeros( (nz , nfiles) )  
      minref_th_profile  = np.zeros( (nz , nfiles) )
      num_th_profile  = np.zeros( (nz , nfiles) )

   z_th_profile[:,ifile] = zp
   meanref_th_profile[:,ifile] = meanp
   stdref_th_profile[:,ifile]  = stdp
   maxref_th_profile[:,ifile]  = maxp  
   minref_th_profile[:,ifile]  = minp
   num_th_profile[:,ifile]   = nump

plt.figure()
plt.pcolormesh( meanref_th_profile )
plt.colorbar()
plt.show()

#plt.plot(ref_cilindro, alt_cilindro , '.')
#plt.plot( meanp , zp , 'r-')
#plt.plot( meanp - stdp , zp , 'k--')
#plt.plot( meanp + stdp , zp , 'k--')
#plt.show()



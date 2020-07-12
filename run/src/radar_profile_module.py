#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 19:33:33 2020

@author: juan
"""
import os
import glob
import scipy.io as sio
import numpy as np
from src import rvd_read as rr 
from datetime import datetime as dt
from datetime import timedelta
import gc 
#from tqdm import tqdm


def get_stations( database_file ) :

    f = open( database_file , "r" )
    station_dict=dict()
    cont = True
    while( cont ) :
       tmp=f.readline()
       if len(tmp) > 1 :
          newline=tmp.split(';')
          station_dict[newline[0]]=dict()
          station_dict[newline[0]]['lon']=newline[1]
          station_dict[newline[0]]['lat']=newline[2]
          station_dict[newline[0]]['alt']=newline[3]
          station_dict[newline[0]]['name']=newline[4].rstrip("\n")
       else :
          cont=False

    return station_dict


def get_cases( database_file )  :
    
    f = open( database_file , "r" )
    cases_dict=dict()
    cases_dict['ID']=list()
    cases_dict['date']=list()
    cases_dict['rank']=list()
    cases_dict['pp6h']=list()
    cont = True
    while( cont ) :
       tmp=f.readline()
       if len(tmp) > 1 :
          newline=tmp.split(',')          
          cases_dict['ID'].append( newline[0] )
          cases_dict['date'].append( dt.strptime( newline[2] + '/' + newline[1] ,'%H/%d/%m/%Y' ) )
          cases_dict['rank'].append( newline[3] )
          cases_dict['pp6h'].append( newline[4].rstrip("\n") )
       else :
          cont=False

    return cases_dict

def get_filelist( ini_date , end_date , data_path , name_filter=None , sufix='.z.rvd' , data_order='RVD' )  :
    #Inidate esta en el formato YYYYMMDDHHmm
    #Enddate esta en el formato YYYYMMDDHHmm

    if data_order == 'RVD' :

       #Generamos la lista de dias que abarcan el periodo solicitado.
       ini_date_dt = dt.strptime( ini_date[0:8] , '%Y%m%d' )
       end_date_dt = dt.strptime( end_date[0:8] , '%Y%m%d' )

       c_date = ini_date_dt
       date_folder_list = list() 
       while ( c_date <= end_date_dt ) :

          date_folder_list.append( c_date ) 
          c_date = c_date + timedelta( days=1 )
          
       #Generamos la lista de dias que abarcan el periodo solicitado.
       ini_date_dt = dt.strptime( ini_date[0:10] , '%Y%m%d%H' )
       end_date_dt = dt.strptime( end_date[0:10] , '%Y%m%d%H' )   
          
       my_filelist = list()
       for my_date_folder in date_folder_list :

          my_path = data_path + '/' + str( my_date_folder.year ) + '/' + dt.strftime( my_date_folder , '%Y%m%d' ) + '/rvd/' 
          tmp_filelist = glob.glob(my_path + '*' + sufix )
          tmp_filelist_cp = tmp_filelist.copy()
          
          #Chequeamos si los archivos que encontramos en esta carpeta estan dentro del rango. 
          for my_file in tmp_filelist :
              my_filename = os.path.basename(my_file )
                         
              if name_filter is not None :
                if not any(my_filter in my_filename for my_filter in name_filter ) :
                    tmp_filelist_cp.remove( my_file )
                    continue 
              #print(my_filename)

              tmp_date = dt.strptime( my_filename[12:20] + my_filename[21:25] , '%Y%m%d%H%M' )
              if tmp_date > end_date_dt or tmp_date < ini_date_dt :
                 tmp_filelist_cp.remove( my_file )

          my_filelist = my_filelist + tmp_filelist_cp 

       my_filelist.sort()
       #print(my_filelist)

    else :
       print('Data order not recognized :( ')
       my_filelist=list()

    return my_filelist


def get_profiles( filelist , lonradar , latradar , altradar , lonp , latp , radius ) :

   my_profile = dict()
   
   my_profile['lon'] = lonp
   my_profile['lat'] = latp
   my_profile['radius'] = radius
   my_profile['lonradar'] = lonradar
   my_profile['latradar'] = latradar
   my_profile['altradar'] = lonradar
   my_profile['file_list'] = filelist   
   
   
   nfiles = len( filelist )

   for ifile,my_file in enumerate( filelist ) :
      #print('Reading ',my_file)

      [ref , alt , elev , date ] = extract_profile_data( my_file , radius , lonp , latp , lonradar , latradar , altradar ) 

      [zp , meanp , stdp , minp , maxp , nump] = grid_profile( ref , alt )
   
      if ifile == 0 :
         nz = np.size( zp ) 
         my_profile['z_th_profile']       = np.zeros( (nz , nfiles) )
         my_profile['meanref_th_profile'] = np.zeros( (nz , nfiles) )
         my_profile['stdref_th_profile']  = np.zeros( (nz , nfiles) )
         my_profile['maxref_th_profile']  = np.zeros( (nz , nfiles) )  
         my_profile['minref_th_profile']  = np.zeros( (nz , nfiles) )
         my_profile['num_th_profile']     = np.zeros( (nz , nfiles) )
         my_profile['date']               = list()

      my_profile['z_th_profile'][:,ifile]       = zp
      my_profile['meanref_th_profile'][:,ifile] = meanp
      my_profile['stdref_th_profile'][:,ifile]  = stdp
      my_profile['maxref_th_profile'][:,ifile]  = maxp  
      my_profile['minref_th_profile'][:,ifile]  = minp
      my_profile['num_th_profile'][:,ifile]     = nump
      my_profile['date'].append( date )

   return my_profile


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
    date = radar.metadata['start_datetime']

    mascara = distancia <= radius 

    elevations= np.tile( radar.elevation['data'] , (radar.fields['reflectivity']['data'].shape[1] ,1 )).transpose()

    ref = radar.fields['reflectivity']['data'][ mascara ]

    elev = elevations[ mascara ]

    alt  = radar.gate_z['data'][ mascara ]

    return ref , alt , elev , date 









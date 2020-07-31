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

def get_filelist( ini_date , end_date , data_path , name_filter=None , sufix='.z.rmvd' , data_order='RVD' )  :
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

          my_path = data_path + '/' + str( my_date_folder.year ) + '/' + dt.strftime( my_date_folder , '%Y%m%d' ) + '/*/' 
          tmp_filelist = glob.glob(my_path + '*' + sufix )
          #For some dates there is no /rvd/ folder 
          my_path = data_path + '/' + str( my_date_folder.year ) + '/' + dt.strftime( my_date_folder , '%Y%m%d' ) + '/'
          tmp_filelist = tmp_filelist +  glob.glob(my_path + '*' + sufix ) 

          tmp_filelist_cp = tmp_filelist.copy()

          #Chequeamos si los archivos que encontramos en esta carpeta estan dentro del rango. 
          for my_file in tmp_filelist :

              my_filename = os.path.basename( my_file )
                         
              if name_filter is not None :
                if not any(my_filter in my_filename for my_filter in name_filter ) :
                    tmp_filelist_cp.remove( my_file )
                    continue 

              tmp_date = dt.strptime( my_filename[12:20] + my_filename[21:25] , '%Y%m%d%H%M' )
              if tmp_date > end_date_dt or tmp_date < ini_date_dt :
                 tmp_filelist_cp.remove( my_file )

          my_filelist = my_filelist + tmp_filelist_cp 

       my_filelist.sort()

    else :
       arint('Data order not recognized :( ')
       my_filelist=list()

    return my_filelist

def get_profiles( conf ) : #filelist , lonradar , latradar , altradar , lonp , latp , radius ) :
   #Esta funcion ingresa la lon y lat de un punto en el area del radar.
   #La funcion toma un cilindro de radio "radius" alrededor del punto lonp,latp y 
   #genera un perfil de reflectividad promedio (de los puntos en el interior del cilindro)
   #El perfil se genera con una cierta resolucion espacial (dz=500 m por defecto) y 
   #por defecto para todas las alturas entre 0 y 15000 metros.
   #El perfil se guarda en un diccionario cuyas entradas incluyen:
   #z_th_profile , la altura sobre el terreno
   #meanref_th_profile , la reflectividad media a cada altura
   #stdref_th_profile , la desviacion estandard de la reflectividad a cada altura
   #calculada sobre todos los puntos que fueron usados para obtener el valor de reflectividad a una altura determinada
   #maxref_th_profile , la maxima reflectividad detectada a una determinada altura
   #minref_th_profile , la minimia reflectvidad detectada
   #num_th_profile , la cantidad de valores de reflectividad promediados para cada altura.
   #date , un objeto fecha de python indicando la hora a la que fue tomado el volumen.

   my_profile = dict()

   my_profiles = conf.copy
   
   nfiles = len( conf['filelist'] )

   for ifile,my_file in enumerate( conf['filelist'] ) :
      #print('Reading ',my_file)

      #Esta funcion devuelve todos los puntos que estan dentro del cilindro (ref,alt,elev)
      #Tambien devuelve el perfil de vecinos mas cercanos nn_ref , nn_alt y nn_elev
      [ref , alt , elev , date , nn_ref , nn_alt , nn_elev ] = extract_profile_data( my_file , conf ) #radius , lonp , latp , lonradar , latradar , altradar ) 

      [zp , meanp , stdp , minp , maxp , nump, etop , vil , vild] = grid_profile( ref , alt , conf )
   
      if ifile == 0 :
         nz = np.size( zp ) 
         my_profile['z_th_profile']       = np.zeros( (nz , nfiles) )
         my_profile['meanref_th_profile'] = np.zeros( (nz , nfiles) )
         my_profile['stdref_th_profile']  = np.zeros( (nz , nfiles) )
         my_profile['maxref_th_profile']  = np.zeros( (nz , nfiles) )  
         my_profile['minref_th_profile']  = np.zeros( (nz , nfiles) )
         my_profile['num_th_profile']     = np.zeros( (nz , nfiles) )
         my_profile['z_raw_profile']      = list()
         my_profile['alt_raw_profile']    = list()
         my_profile['elev_raw_profile']   = list()
         my_profile['z_nn_profile']       = list()
         my_profile['ref_nn_profile']     = list()
         my_profile['elev_nn_profile']    = list()
         my_profile['date']               = list()
         my_profile['etop']               = np.zeros( nfiles )
         my_profile['vil']                = np.zeros( nfiles )
         my_profile['vild']               = np.zeros( nfiles )


      my_profile['z_th_profile'][:,ifile]       = zp
      my_profile['meanref_th_profile'][:,ifile] = meanp
      my_profile['stdref_th_profile'][:,ifile]  = stdp
      my_profile['maxref_th_profile'][:,ifile]  = maxp  
      my_profile['minref_th_profile'][:,ifile]  = minp
      my_profile['num_th_profile'][:,ifile]     = nump
      my_profile['date'].append( date )
      my_profile['z_raw_profile'].append( ref )
      my_profile['alt_raw_profile'].append( alt )
      my_profile['elev_raw_profile'].append( elev )
      my_profile['z_nn_profile'].append( nn_alt )
      my_profile['ref_nn_profile'].append( nn_ref )
      my_profile['elev_nn_profile'].append( nn_elev )
      my_profile['etop'][ifile] = etop
      my_profile['vil'][ifile]  = vil
      my_profile['vild'][ifile] = vild

   return my_profile


def grid_profile( ref , z , conf ) : #z_min=0.0 , z_max=15000.0 , delta_z = 500.0 , undef = -32.0 ) :


   z_min = conf['z_min']
   z_max = conf['z_max']
   delta_z = conf['delta_z']
   undef = conf['undef']
   min_ref = conf['min_ref_profile']
   min_ref_etop = conf['min_ref_etop']

   zp = np.arange( z_min + delta_z / 2.0 , z_max + delta_z / 2.0 , delta_z )
   zp_min = np.arange( z_min , z_max , delta_z )
   zp_max = zp_min + delta_z
   nbin = np.size( zp_min ) 

   meanp=np.zeros(nbin) + np.nan
   stdp = np.copy(meanp)
   minp = np.copy(meanp)
   maxp = np.copy(meanp)
   nump = np.zeros(nbin)   

   #ref[ref==undef] = np.nan

   for ii in range( nbin ) :
       my_mask = ( z <= zp_max[ii] ) & ( z >= zp_min[ii] ) & ( ~ np.isnan( ref ) )  & ( ref > min_ref ) & ( ref != undef )
       tmp_num = np.sum( my_mask )
       if  tmp_num > 30 :
          meanp[ii] = np.mean( ref[my_mask] )
          stdp[ii]  = np.std ( ref[my_mask] )
          minp[ii]  = np.min ( ref[my_mask] )
          maxp[ii]  = np.max ( ref[my_mask] )
          nump[ii]  = tmp_num

   #Detect the highest echo top
   etop = np.nan  #Echo top
   vil  = np.nan  #Vertically integrated liquid
   vild = np.nan  #VIL density
   for ii in range( nbin -1 )  :
       if ( meanp[ii+1] < min_ref_etop ) & ( meanp[ii] >= min_ref_etop ) :
          etop = zp[ii] 

       if ~ np.isnan( meanp[ii+1] ) & ~ np.isnan( meanp[ii] ) :
          tmp_z_mean = 0.5 * (meanp[ii+1] + meanp[ii])
          if tmp_z_mean > 56.0 :
             tmp_z_mean = 56.0
          vil_inc = 3.44e-6 * ( ( tmp_z_mean )**(4.0/7.0) ) * delta_z
          if np.isnan( vil ) :
             vil = vil_inc
          else               :
             vil = vil + vil_inc
   if ~ np.isnan( etop ) & ~ np.isnan( vil ) :
       vild = 1000.0 * vil / etop 

   

   return zp , meanp , stdp , minp , maxp , nump , etop , vil , vild

def extract_profile_data( filename , conf ) : #radius , lonp , latp , lonradar , latradar , altradar ) :
    
    radar = rr.rvd_read( filename , conf['lonradar'] , conf['latradar'] , conf['altradar'] )

    # Buscamos los puntos que estan en el cilindro y calculamos el perfil medio sobre el cilindro.
    dlon =  np.cos(radar.gate_latitude['data']*np.pi/180.0)*( radar.gate_longitude['data'] - conf['lon'] )    
    dlat =  ( radar.gate_latitude['data'] - conf['lat'] )
    distancia = np.sqrt( ( dlon * 111000.0 )**2 + ( dlat * 111000.0 )**2 )
    date = radar.metadata['start_datetime']

    mascara = distancia <= conf['radius'] 

    elevations= np.tile( radar.elevation['data'] , (radar.fields['reflectivity']['data'].shape[1] ,1 )).transpose()

    ref = radar.fields['reflectivity']['data'][ mascara ]

    elev = elevations[ mascara ]

    alt  = radar.gate_z['data'][ mascara ]

    #Buscamos los puntos que son el vecino mas cercano (nn) del punto seleccionado en cada ppi para
    #construir el perfil de vecinos mas cercanos. 

    elevations = radar.elevation['data']
    unique_elevs = np.unique( elevations ) 
    nelevs = np.size( unique_elevs )
    
    ref_nn = np.zeros( nelevs )
    alt_nn = np.zeros( nelevs )

    for ie,my_elev in enumerate( unique_elevs ) :
        my_mask = elevations == my_elev
        my_ref  = radar.fields['reflectivity']['data'][my_mask,:] 
        my_z    = radar.gate_z['data'][my_mask,:]
        my_dist = distancia[my_mask,:]
        [minx , miny] = np.where( my_dist == np.min( my_dist ) )
        ref_nn[ie] = np.copy( my_ref[minx,miny] )
        alt_nn[ie] = np.copy( my_z[minx,miny] )       

    return ref , alt , elev , date , ref_nn , alt_nn , unique_elevs 





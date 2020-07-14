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

      #Esta funcion devuelve todos los puntos que estan dentro del cilindro (ref,alt,elev)
      #Tambien devuelve el perfil de vecinos mas cercanos nn_ref , nn_alt y nn_elev
      [ref , alt , elev , date , nn_ref , nn_alt , nn_elev ] = extract_profile_data( my_file , radius , lonp , latp , lonradar , latradar , altradar ) 

      [zp , meanp , stdp , minp , maxp , nump] = grid_profile( ref , alt )
   
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

    # Buscamos los puntos que estan en el cilindro y calculamos el perfil medio sobre el cilindro.
    dlon =  np.cos(radar.gate_latitude['data']*np.pi/180.0)*( radar.gate_longitude['data'] - lonp )    
    dlat =  ( radar.gate_latitude['data'] - latp )
    distancia = np.sqrt( ( dlon * 111000.0 )**2 + ( dlat * 111000.0 )**2 )
    date = radar.metadata['start_datetime']

    mascara = distancia <= radius 

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


def get_missing_periods( data , delta_t_max ) :
   #Hacemos un resumen de los datos faltantes.
   dt_max = timedelta( minutes = delta_t_max )    #Este es el maximo dt que toleramos entre 2 volumenes.

   data['datedt'] = list()
   for idate , my_date in enumerate( data['date'] ) :
       data['datedt'].append( dt.strptime( data['date'][idate] , '%Y-%m-%dT%H:%M:%SZ' ) )

   missing_list = list()  
   missing_index_list = list()
   for idate , my_date in enumerate( data['datedt'] ) : 
      if idate == 0 :
         dtt = data['datedt'][idate] - dt.strptime( data['ini_date'] , '%Y%m%d%H' )
      else          :
         dtt = data['datedt'][idate] - data['datedt'][idate-1] 
      if dtt > dt_max :
         missing_index_list.append( idate )
         missing_list.append( dtt.total_seconds() )
   dtt = dt.strptime( data['end_date'] , '%Y%m%d%H' ) - data['datedt'][idate-1]  
   #Para el ultimo elemento.
   if dtt > dt_max :
      missing_index_list.append( idate+1 )
      missing_list.append( dtt.total_seconds() ) 

   #La funcion agrega el resultado en el mismo diccionario data. Es decir el diccionario
   #data que devuelve la funcion tiene todo lo que tenia el original + estas 3 variables que indican la presencia de periodos faltantes.

   data['missing_index_list'] = missing_index_list
   data['missing_list']        = missing_list
   data['dt_max']             = dt_max 

   return data 

def add_missing_periods( data , delta_t_max ) :

   data = get_missing_periods( data , delta_t_max )

   #A partir de las variable originales, generamos una nueva variable y su correspondiente escala de tiempos
   #que tiene nans en las zonas donde hay datos faltantes (huecos en los datos mayores a dt_max )
   nz = data['meanref_th_profile'].shape[0]
   nt = data['meanref_th_profile'].shape[1]
   dt_add = timedelta( minutes = 1.0 )

   data['meanref_th_profile_ext'] = np.zeros( ( nz , nt + 2*len(data['missing_list']) ) ) + np.nan
   data['maxnref_th_profile_ext'] = np.zeros( ( nz , nt + 2*len(data['missing_list']) ) ) + np.nan
   data['minnref_th_profile_ext'] = np.zeros( ( nz , nt + 2*len(data['missing_list']) ) ) + np.nan
   data['stdref_th_profile_ext'] = np.zeros( ( nz , nt + 2*len(data['missing_list']) ) ) + np.nan
   data['num_th_profile_ext'] = np.zeros( ( nz , nt + 2*len(data['missing_list']) ) ) + np.nan
   
   data['date_ext']  = list()

   dindex = 0
   for it in range( nt ) :

      if ( it in data['missing_index_list'] ) :
          dindex = dindex + 2
          if  it == 0 :
             data['date_ext'].append( dt.strptime( data['ini_date'] , '%Y%m%d%H' ) )
             data['date_ext'].append( data['datedt'][0] - dt_add )           
          if  it > 0 & it < nt-1 :
             data['date_ext'].append( data['datedt'][it-1] + dt_add )
             data['date_ext'].append( data['datedt'][it] - dt_add )

      data['meanref_th_profile_ext'][:,it+dindex] = data['meanref_th_profile'][:,it]
      data['maxnref_th_profile_ext'][:,it+dindex] = data['maxref_th_profile'][:,it]
      data['minnref_th_profile_ext'][:,it+dindex] = data['minref_th_profile'][:,it]
      data['stdref_th_profile_ext'][:,it+dindex] = data['stdref_th_profile'][:,it]  
      data['num_th_profile_ext'][:,it+dindex] = data['num_th_profile'][:,it]

      data['date_ext'].append( data['datedt'][it] )

   if nt in data['missing_index_list'] :
      data['date_ext'].append( data['datedt'][-1] + dt_add )  
      data['date_ext'].append( dt.strptime( data['end_date'] , '%Y%m%d%H' ) )


   return data










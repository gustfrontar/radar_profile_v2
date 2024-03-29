#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 19:33:33 2020

@author: juan
"""
import numpy as np
from datetime import datetime as dt
from datetime import timedelta

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
   data['maxref_th_profile_ext'] = np.zeros( ( nz , nt + 2*len(data['missing_list']) ) ) + np.nan
   data['minref_th_profile_ext'] = np.zeros( ( nz , nt + 2*len(data['missing_list']) ) ) + np.nan
   data['stdref_th_profile_ext'] = np.zeros( ( nz , nt + 2*len(data['missing_list']) ) ) + np.nan
   data['vil_ext'] = np.zeros( ( nt + 2*len(data['missing_list']) ) ) + np.nan
   data['etop_ext'] = np.zeros( ( nt + 2*len(data['missing_list']) ) ) + np.nan
   data['num_th_profile_ext'] = np.zeros( ( nz , nt + 2*len(data['missing_list']) ) ) + np.nan
   
   data['z_raw_profile_ext'] = list()
   data['alt_raw_profile_ext'] = list()

   data['date_ext']  = list()

   dindex = 0
   for it in range( nt ) :

      if ( it in data['missing_index_list'] ) :
          dindex = dindex + 2
          data['z_raw_profile_ext'].append( [] )
          data['alt_raw_profile_ext'].append( [] )
          data['z_raw_profile_ext'].append( [] )
          data['alt_raw_profile_ext'].append( [] )
          
          if  it == 0 :
             data['date_ext'].append( dt.strptime( data['ini_date'] , '%Y%m%d%H' ) )
             data['date_ext'].append( data['datedt'][0] - dt_add )
          if  it > 0 & it < nt-1 :
             data['date_ext'].append( data['datedt'][it-1] + dt_add )
             data['date_ext'].append( data['datedt'][it] - dt_add )

      data['meanref_th_profile_ext'][:,it+dindex] = data['meanref_th_profile'][:,it]
      data['maxref_th_profile_ext'][:,it+dindex] = data['maxref_th_profile'][:,it]
      data['minref_th_profile_ext'][:,it+dindex] = data['minref_th_profile'][:,it]
      data['stdref_th_profile_ext'][:,it+dindex] = data['stdref_th_profile'][:,it]
      data['num_th_profile_ext'][:,it+dindex] = data['num_th_profile'][:,it]
      data['vil_ext'][it+dindex] = data['vil'][it]
      data['etop_ext'][it+dindex] = data['etop'][it]
      data['z_raw_profile_ext'].append( data['z_raw_profile'][it] )
      data['alt_raw_profile_ext'].append( data['alt_raw_profile'][it] )

      data['date_ext'].append( data['datedt'][it] )

   if nt in data['missing_index_list'] :
      data['date_ext'].append( data['datedt'][-1] + dt_add )
      data['date_ext'].append( dt.strptime( data['end_date'] , '%Y%m%d%H' ) )
      data['z_raw_profile_ext'].append( [] )
      data['alt_raw_profile_ext'].append( [] )
      data['z_raw_profile_ext'].append( [] )
      data['alt_raw_profile_ext'].append( [] )



   return data

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 18:59:51 2020

@author: jruiz
"""

#Importamos todas las librerias necesarias
import os
import numpy as np
import pickle as pkl
import glob
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime as dt
from datetime import timedelta
import matplotlib.dates as mdates
#Estas dos lineas actualizan el repositorio (en caso que hayamos hecho cambios durante la sesion en colab)
import sys
from src import missing_period_module as mpm 


pkl_path = '../pkl/'
#Definimos el evento que vamos a estudiar
event='87585_2013040212'

#Los casos posibles son: 87553_2013040212 ; 87585_2013040212 ; 87593_2013040300

#Seteamos el rango de variacion de los parametros radio y dz (en m). Estos rangos no pueden ser 
#cualquier cosa, tienen que coincidir con los que generamos oportunamente. Podemos
#sacar valores pero no agregar nuevos ya que eso requeriria hacer los calculos en el servidor 
#para esos valores. Si del analisis preliminar se desprende que estaria bueno provar con mas 
#valores en algun rango vamos al servidor y generamos esos experimentos para analizarlos aqui.

radius = [ 1000.0 , 2500.0 , 5000.0 , 10000.0 ]
dz = [ 1000.0 , 2000.0 ]

for my_radius in radius :
  for my_dz in dz :
      #Genero el nombre del pkl que guarda los datos para el evento seleccionado y para la combinacion de radio y dz correspondiente.
      my_file = pkl_path + '/event_profiles_' + event + '_sens_r' + str(my_radius) + '_dz_' + str(my_dz) + '.pkl' 

      my_f = open( my_file , 'rb' )  #my_f es un objeto que refiere al archivo y su contenido (similar a cuando leemos un archivo en fortran y le asignamos un numero)  

      data = pkl.load( my_f )        #Aca hacemos la lectura

      #Le agrego a data las entradas correspondientes a las series completadas con los datos faltantes (para generar las figuras con los blancos en donde faltan datos)
      #add_missing_periods es un par de funciones que busca los periodos faltantes y completa las series temporales de los perfiles de manera acorde.

      #TODO: Ahora agregamos muchos NaNs debido a que solo promediamos valores por encima de 5dBZ esto hace que graficamente los tiempos que no hay lluvia
      #sean indistinguibles de los faltantes. Hay que ver como solucionar eso. 
      data = mpm.add_missing_periods( data , 20 ) 


      #==============================================================================#
      # COMIENZA EL GRAFICADO

      n_times = len(data['date_ext'])
      times = np.empty(n_times)

    
      for ii in range(n_times):
         dateobj =   dt.strftime( data['date_ext'][ii] , '%Y-%m-%d-%H:%M:%S' ) 
         times[ii] = mdates.datestr2num(dateobj)


      #==============================================================================#
      # GRAFICO REFLECTIVIDAD MEDIA EN FUNCION DE LA ALTURA Y EL TIEMPO
      #==============================================================================#

      plt.figure(figsize=[24,8])

      plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
      plt.gca().xaxis.set_major_locator(mdates.HourLocator())

      plt.pcolormesh(times,
               data['z_th_profile'][:,0]/1000,
               data['meanref_th_profile_ext'],
               cmap='gist_ncar', 
               vmin=5,
               vmax=55)     
      plt.colorbar(shrink=0.7, pad=0.01)
      plt.gcf().autofmt_xdate()
      plt.yticks(np.arange(0,16))

      plt.title(data['ID']+' '+data['name']+'  Periodo: '+data['ini_date']+'Z-'+data['end_date']+'Z' + ' R=' + str(my_radius) + ' DZ=' + str(my_dz), fontsize=18)
      plt.xlabel('Tiempo [UTC]', fontsize=14)
      plt.ylabel('Altura [km]', fontsize=14)

      plt.show()


      #==============================================================================#
      # PARA EL PUNTO DONDE LA REFLECTIVIDAD ES MAS INTENSA (maxima suma de la ref en la vertical)
      # GRAFICO LA NUBE ORIGINAL DE PUNTOS, EL PERFIL MEDIO, EL PERFIL DE LOS VECINOS MAS CERCANOS
      # EL PERFIL MAXIMO, EL MINIMO Y LA DESVIACION ESTANDARD.
      #==============================================================================#
      sum_ref = np.nansum( data['meanref_th_profile_ext'] , axis=0 )
      sum_ref[np.isnan( sum_ref )] = 0.0
      #Selecciono el indice temporal donde la suma vertical de la reflectividad es maxima
      #(esto es solo un criterio arbitario, se puede fijar un tiempo que sea interesante para el 
      #evento en base a cualquier otro criterio)
      time_profile = np.argmax( sum_ref )

      plt.figure(figsize=[24,8])

      plt.plot( data['z_raw_profile'][time_profile] , data['alt_raw_profile'][time_profile]/1000 ,'ko',alpha=0.3)
      plt.plot( data['meanref_th_profile'][:,time_profile] , data['z_th_profile'][:,time_profile]/1000 ,'-bo',label='mean')

      plt.plot( data['maxref_th_profile'][:,time_profile]  , data['z_th_profile'][:,time_profile]/1000 ,'-go',label='Max')
      plt.plot( data['minref_th_profile'][:,time_profile]  , data['z_th_profile'][:,time_profile]/1000 ,'-go',label='Min')
      plt.plot( data['meanref_th_profile'][:,time_profile] + data['stdref_th_profile'][:,time_profile] , data['z_th_profile'][:,time_profile]/1000 ,'--b')
      plt.plot( data['meanref_th_profile'][:,time_profile] - data['stdref_th_profile'][:,time_profile] , data['z_th_profile'][:,time_profile]/1000 ,'--b')
      plt.plot( data['ref_nn_profile'][time_profile]       , data['z_nn_profile'][time_profile]/1000 ,'-ro')

      plt.legend()
      plt.title('Datos de reflectividad crudos para el tiempo ' + str(time_profile) + ' R=' + str(my_radius) + ' DZ=' + str(my_dz), fontsize=18)
      plt.xlabel('Reflectividad', fontsize=14)
      plt.ylabel('Altura [km]', fontsize=14)
      
      plt.ylim([0,15])
      plt.xlim([-10,60])

      plt.show()

      #==============================================================================#
      # GRAFICO REFLECTIVIDAD MEDIA EN FUNCION DE LA ALTURA Y EL TIEMPO
      #==============================================================================#
      n_times=len( data['z_raw_profile'] )
      drefdz = np.zeros( n_times )  #La pendiente de reflectividad con la altura (primeros 5 km)
      drefdzs = np.zeros( n_times ) #Igual a la anterior pero con un suavizado temporal.
      maxref = np.zeros( n_times )  #La maxima reflectividad (primeros 5 km)
      vil    = np.zeros( n_times )  #El vil
      etop   = np.zeros( n_times )  #El echo top a partir del perfil medio.
      reflow = np.zeros( n_times )  #La reflectividad en el nivel mas bajo

      for mytime in range( n_times ) :
         
          tmp_ref = np.copy( data['z_raw_profile'][mytime] )
          tmp_z   = np.copy( data['alt_raw_profile'][mytime]/1000 )
          #Me quedo con los datos por debajo de 5 km.
          tmp_ref = tmp_ref[ tmp_z < 5.0 ]
          tmp_z   = tmp_z[ tmp_z < 5.0 ]
          if tmp_ref.size > 0 :
             pars = np.polyfit( tmp_z , tmp_ref , 1 )
             drefdz[mytime] = pars[0]
             maxref[mytime] = np.max( tmp_ref )
          else  :
             maxref[mytime] = np.nan
             drefdz[mytime] = np.nan
          vil[mytime] = data['vil'][mytime]
          reflow[mytime] = data['meanref_th_profile'][0,mytime]
       

      for mytime in range( n_times ) :
          mini = mytime - 3
          maxi = mytime + 3
          if mini < 0 : 
             mini=0
          if maxi > n_times-1 : 
             maxi= n_times-1
          drefdzs[mytime] = np.nanmean( drefdz[mini:maxi] )
          
      plt.figure(figsize=[24,8])
      plt.plot( reflow , '-b' , label='RefLow')
      plt.plot( maxref , '-r', label='MaxRef' )
      plt.plot( vil * 100 , '-g' , label='Vil * 1e2')
      plt.plot( drefdzs * 10 , '-k' , label='DZr/Dz * 10.0')
      plt.ylim([-50,60])
      plt.grid()
      plt.legend()
      plt.title('Series de tiempo para ' + str(time_profile) + ' R=' + str(my_radius) + ' DZ=' + str(my_dz), fontsize=18)
      plt.show()



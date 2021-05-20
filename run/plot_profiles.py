# -*- coding: utf-8 -*-
"""
La idea de este script es tomar la lista de eventos y para cada evento graficar el pseudo-vertical profile de radar.
"""


#Importamos todas las librerias necesarias
import os
import numpy as np
import pickle as pkl
import glob
import matplotlib 
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


from datetime import datetime as dt
from datetime import timedelta
import matplotlib.dates as mdates
import sys
sys.path.insert(1,'../run/src/')
import missing_period_module as mpm
import radar_profile_module  as rpm
import warnings 
warnings.filterwarnings('ignore')
#================== CONFIGURACION ========================
station_database_file = './csv/estaciones.csv'
event_list_file = '../run/csv/lista_casos_con_radar_sondeo.csv'
fig_path='../fig/'
pkl_path='../pkl/'
rayos_path='../rayos/pkl/'

radius = 5000.0      #Radio del cilindro en el pseudo vertical profile.
dz = 1000.0          #Resolucion vertical del pseudo vertical profile.
vil_tr = 0.5         #Threshold de VIL que se usa para calcular el valor medio
                     #de diferentes parametros sobre el periodo de mayor precipitacion.

#==========================================================


#LEO LA LISTA DE EVENTOS. 
#sounding_data = True lee las isotermas de 0, -10 y -20 del sondeo mas cercano en espacio y tiempo.
cases_dict = rpm.get_cases( event_list_file , sounding_data = True )

ncases = len( cases_dict['date'] )  #Obtengo el numero de casos a procesar.

#cases_dict['H0degC'] = np.ones( ncases ) * 3000.0 

#Cargo la base de datos de las estaciones.
station_dict = rpm.get_stations( station_database_file )

#==============================================================================#
# COMIENZA EL GRAFICADO

for ievent in range( ncases ) :

   pkl_file = pkl_path + '/event_profiles_' + cases_dict['ID'][ievent] + '_' + dt.strftime( cases_dict['date'][ievent] , '%Y%m%d%H' ) + '.pkl'
   data = pkl.load( open( pkl_file , 'rb' ) )   #Aca hacemos la lectura

   #Le agrego a data las entradas correspondientes a las series completadas con los datos faltantes (para generar las figuras con los blancos en donde faltan datos)
   #add_missing_periods es un par de funciones que busca los periodos faltantes y completa las series temporales de los perfiles de manera acorde.

   data['meanref_th_profile'][ np.isnan(data['meanref_th_profile']) ] = 0.0
   data = mpm.add_missing_periods( data , 20 )
 
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
   plt.plot([times[0],times[-1]],[cases_dict['H0degC'][ievent]/1000.0,cases_dict['H0degC'][ievent]/1000.0],'--k',linewidth=4) 
   plt.plot([times[0],times[-1]],[cases_dict['Hm20degC'][ievent]/1000.0,cases_dict['Hm20degC'][ievent]/1000.0],'--k',linewidth=4)
    
   plt.colorbar(shrink=0.7, pad=0.01)
   plt.gcf().autofmt_xdate()
   plt.yticks(np.arange(0,16))
   plt.title(data['ID']+' '+data['name']+'  Periodo: '+data['ini_date']+'Z-'+data['end_date']+'Z' + ' R=' + str(radius) + ' DZ=' + str(dz), fontsize=18)
   plt.xlabel('Tiempo [UTC]', fontsize=14)
   plt.ylabel('Altura [km]', fontsize=14)
   plt.savefig( fig_path + '/MeanRefZvsTime_' + cases_dict['ID'][ievent] + '_' + dt.strftime( cases_dict['date'][ievent] , '%Y%m%d%H' ) + '.png' )
   plt.close()

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
   plt.title('Datos de reflectividad crudos para el tiempo ' + str(time_profile) + ' R=' + str(radius) + ' DZ=' + str(dz), fontsize=18)
   plt.xlabel('Reflectividad', fontsize=14)
   plt.ylabel('Altura [km]', fontsize=14)
   
   plt.ylim([0,15])
   plt.xlim([-10,60])
   plt.savefig( fig_path + '/MeanRefProfile_atMaxRef_' + cases_dict['ID'][ievent] + '_' + dt.strftime( cases_dict['date'][ievent] , '%Y%m%d%H' ) + '.png' )

   plt.close()

   #==============================================================================#
   # GRAFICO GRADIENTE VERTICAL DE REFLECTIVIDAD Y REFLECTIVIDAD MAXIMA EN DIFERENTES CAPAS
   #==============================================================================#
   #Check if we have valid sounding data.
   #If we dont have a valid height for any temperatur, then we replace the missing value by 
   #the average height of that temperature. 
   if cases_dict['H0degC'][ievent] == -999 :
      tmp_array = np.array( cases_dict['H0degC'] )
      cases_dict['H0degC'][ievent] = np.mean( tmp_array[tmp_array > 0] )
   if cases_dict['Hm10degC'][ievent] == -999 :
      tmp_array = np.array( cases_dict['Hm10degC'] )
      cases_dict['Hm10degC'][ievent] = np.mean( tmp_array[tmp_array > 0] )
   if cases_dict['Hm20degC'][ievent] == -999 :
      tmp_array = np.array( cases_dict['Hm20degC'] )
      cases_dict['Hm20degC'][ievent] = np.mean( tmp_array[tmp_array > 0] )


   n_times=len( data['z_raw_profile_ext'] )
   #Altura maxima y minima para el calculo de la pendiente de reflectividad por debajo de la banda brillante.
   drefdzlow            = np.zeros( n_times )  #La pendiente de reflectividad con la altura (entre la sup y la banda brillante).
   drefdzlows           = np.zeros( n_times )  #Igual a la anterior pero con un suavizado temporal.
   #Altura maximia y minima para el calculo de la pendiente de reflectividad entre la banda brillante y el echo top
   drefdzhi            = np.zeros( n_times )  #La pendiente de reflectividad con la altura (entre la sup y la banda brillante).
   drefdzhis           = np.zeros( n_times )  #Igual a la anterior pero con un suavizado temporal.
   #Altura maximia y minima para el calculo de la pendiente de reflectividad entre la banda brillante el nivel de -20
   drefdzm20            = np.zeros( n_times )  #La pendiente de reflectividad con la altura (entre la sup y la banda brillante).
   drefdzm20s           = np.zeros( n_times )  #Igual a la anterior pero con un suavizado temporal.

   #Reflectivity at -10 and -20 C
   refm10 = np.zeros( n_times ) + np.nan
   refm20 = np.zeros( n_times ) + np.nan
   maxreflow = np.zeros( n_times )  #La maxima reflectividad por debajo de la banda brillante
   maxrefhi  = np.zeros( n_times )  #La maxima reflectividad por encima de la banda brillante
   maxrefm20 = np.zeros( n_times )  #La reflectividad maxima entre la banda brillante y los -20C
   vil    = np.zeros( n_times )  #El vil
   etop   = np.zeros( n_times )  #El echo top a partir del perfil medio.
   reflow = np.zeros( n_times )  #La reflectividad en el nivel mas bajo

   for mytime in range( n_times ) :

       etop[mytime] = data['etop_ext'][mytime]
         
       myref = np.copy( data['z_raw_profile_ext'][mytime] )
       myz   = np.copy( data['alt_raw_profile_ext'][mytime] )
       #Obtengo los puntos de reflectividad por debajo de la banda brillante
       mymask = np.logical_and( 1000.0  < myz ,  cases_dict['H0degC'][ievent] - 500.0 > myz)
       tmp_ref = myref[ mymask ]
       tmp_z   = myz[ mymask ]/1000.0
       if tmp_ref.size > 0 :
          pars = np.polyfit( tmp_z , tmp_ref , 1 )
          drefdzlow[mytime] = pars[0]
          maxreflow[mytime] = np.max( tmp_ref )
       else  :
          maxreflow[mytime] = np.nan
          drefdzlow[mytime] = np.nan

       #Obtengo los puntos de reflectividad por encima de la banda brillante y por debajo del echotop
       mymask = np.logical_and( cases_dict['H0degC'][ievent] + 500.0 < myz ,  data['etop_ext'][mytime] > myz)
       tmp_ref = myref[ mymask ]
       tmp_z   = myz[ mymask ]/1000.0
       if tmp_ref.size > 0 :
          pars = np.polyfit( tmp_z , tmp_ref , 1 )
          drefdzhi[mytime] = pars[0]
          maxrefhi[mytime] = np.max( tmp_ref )
       else  :
          maxrefhi[mytime] = np.nan
          drefdzhi[mytime] = np.nan

       #Obtengo los puntos de reflectividad por encima de la banda brillante y por de la temperatura de -20C
       mymask = np.logical_and( cases_dict['H0degC'][ievent] + 500.0 < myz , cases_dict['Hm20degC'][ievent] > myz)
       tmp_ref = myref[ mymask ]
       tmp_z   = myz[ mymask ]/1000.0
       if tmp_ref.size > 0 :
          pars = np.polyfit( tmp_z , tmp_ref , 1 )
          drefdzm20[mytime] = pars[0]
          maxrefm20[mytime] = np.max( tmp_ref )
       else  :
          maxrefm20[mytime] = np.nan
          drefdzm20[mytime] = np.nan

       #Obtengo la reflectividad media en un kilometro alrededor de -20 C
       mymask = np.logical_and( cases_dict['Hm20degC'][ievent] - 500.0  < myz ,  cases_dict['Hm20degC'][ievent] + 500.0 > myz )
       tmp_ref = myref[ mymask ]
       refm20[mytime] = np.nanmean( tmp_ref )

       #Obtengo la reflectividad media en un kilometro alrededor de -10 C
       mymask = np.logical_and( cases_dict['Hm10degC'][ievent] - 500.0  < myz ,  cases_dict['Hm10degC'][ievent] + 500.0 > myz )
       tmp_ref = myref[ mymask ]
       refm10[mytime] = np.nanmean( tmp_ref )

       vil[mytime] = data['vil_ext'][mytime]
       reflow[mytime] = data['meanref_th_profile_ext'][0,mytime]
       
   smooth_window_size = 3 
   for mytime in range( n_times ) :
       mini = mytime - smooth_window_size  # 3
       maxi = mytime + smooth_window_size  # 3
       if mini < 0 : 
          mini=0
       if maxi > n_times-1 : 
          maxi= n_times-1
       drefdzlows[mytime] = np.nanmean( drefdzlow[mini:maxi] )
       drefdzhis[mytime] = np.nanmean( drefdzhi[mini:maxi] )
       drefdzm20s[mytime] = np.nanmean( drefdzm20[mini:maxi] )

   vil_mask = vil > vil_tr 
   #Generate a report with the mean values for different parameters during period of high VIL.
   print('Propiedades medias en la region que mas llueve')
   print('Evento: ' + cases_dict['ID'][ievent] + ' ' + dt.strftime( cases_dict['date'][ievent] , '%Y%m%d%H' ))
   print('Pendiente ref (sup-0) : ' + str( np.round( np.nanmean( drefdzlow[vil_mask] ),2))) 	 
   print('Pendiente ref (0- -20): ' + str( np.round( np.nanmean( drefdzm20[vil_mask] ),2)))
   print('Pendiente ref (0-etop): ' + str( np.round( np.nanmean( drefdzhis[vil_mask] ),2)))
   print('Ref en el 1er nivel   : ' + str( np.round( np.nanmean( reflow[vil_mask] ),2))) 
   print('Ref en -10            : ' + str( np.round( np.nanmean( refm10[vil_mask] ),2)))
   print('Ref en -20            : ' + str( np.round( np.nanmean( refm20[vil_mask] ),2)))
   print('MaxRef (sup-0)        : ' + str( np.round( np.nanmean( maxreflow[vil_mask] ),2)))
   print('MaxRef (0 - -20)      : ' + str( np.round( np.nanmean( maxrefm20[vil_mask] ),2)))
   print('MaxRef (0-etop )      : ' + str( np.round( np.nanmean( maxrefhi[vil_mask] ),2)))
   print('Etop                  : ' + str( np.round( np.nanmean( etop[vil_mask] ),2)))  
   print()
   print()
   print('===========================================================================')

   #Significant reflectivities           
   plt.figure(figsize=[24,8])
   plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
   plt.gca().xaxis.set_major_locator(mdates.HourLocator())
   plt.plot( times , reflow ,'--r' ,linewidth=3,label='Lowest level ref')
   plt.plot( times , maxreflow ,'-r' ,linewidth=3,label='RefMax (Surface 0C)')
   plt.plot( times , maxrefhi ,'-g',linewidth=3, label='RefMax (0C Echo Top)' )
   plt.plot( times , maxrefm20 ,'-b',linewidth=3, label='RefMax (0C -20C)' )

   plt.ylim([-30,70])
   plt.grid()
   #plt.gcf().autofmt_xdate()
   plt.legend()
   plt.title('Series de tiempo para ' + cases_dict['ID'][ievent] + '_' + dt.strftime( cases_dict['date'][ievent] , '%Y%m%d%H') + ' R=' + str(radius) + '  DZ=' + str(dz) , fontsize=18)
   plt.savefig( fig_path + '/RefTimeSeries_' + cases_dict['ID'][ievent] + '_' + dt.strftime( cases_dict['date'][ievent] , '%Y%m%d%H' ) + '.png' )
   plt.close()


   #Reflectivity slopes
   plt.figure(figsize=[24,8])
   plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
   plt.gca().xaxis.set_major_locator(mdates.HourLocator())

   drefdzlows[np.isnan(drefdzlows)] = 0.0
   drefdzhis[np.isnan(drefdzhis)] = 0.0
   drefdzm20s[np.isnan(drefdzm20s)] = 0.0   
   plt.plot( times , -drefdzlows , '-r' ,linewidth=3, label='-DZr/Dz  (Surface- 0C')   # se puede graficar con o sin suavizado
   plt.plot( times , -drefdzhis  , '-g' ,linewidth=3, label='-DZr/Dz  (0C-Echo Top)')   # se puede graficar con o sin suavizado
   plt.plot( times , -drefdzm20s  , '-b' ,linewidth=3, label='-DZr/Dz  (0C- -20C')   # se puede graficar con o sin suavizado
   plt.ylim([-3,20])
   plt.grid()
   #plt.gcf().autofmt_xdate()
   plt.legend()
   plt.title('Series de tiempo para ' + cases_dict['ID'][ievent] + '_' + dt.strftime( cases_dict['date'][ievent] , '%Y%m%d%H') + ' R=' + str(radius) + '  DZ=' + str(dz) , fontsize=18)
   plt.savefig( fig_path + '/DRefDzTimeSeries_' + cases_dict['ID'][ievent] + '_' + dt.strftime( cases_dict['date'][ievent] , '%Y%m%d%H' ) + '.png' )
   plt.close()

   #VIL          
   plt.figure(figsize=[24,8])
   plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
   plt.gca().xaxis.set_major_locator(mdates.HourLocator())
   vil[np.isnan(vil)]=0.0
   plt.plot( times , vil , '-g' ,linewidth=3, label='Vil ')
   plt.ylim([0,5.0])
   plt.grid()
   #plt.gcf().autofmt_xdate()
   plt.legend()
   plt.title('Series de tiempo para ' + cases_dict['ID'][ievent] + '_' + dt.strftime( cases_dict['date'][ievent] , '%Y%m%d%H') + ' R=' + str(radius) + '  DZ=' + str(dz) , fontsize=18)
   plt.savefig( fig_path + '/VIL_' + cases_dict['ID'][ievent] + '_' + dt.strftime( cases_dict['date'][ievent] , '%Y%m%d%H' ) + '.png' )
   plt.close()










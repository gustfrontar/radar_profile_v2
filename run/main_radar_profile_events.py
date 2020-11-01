from src import radar_profile_module as rpm 
from datetime import datetime as dt
from datetime import timedelta
import pickle as pkl

event_list_file       = './csv/lista_casos.csv'
station_database_file = './csv/estaciones.csv'
radar_data_path       = '/ms270/raddata/base_datos_radar/eze/'
output_path           = '../pkl/'

dt_ini = timedelta(hours=15.0)  #ancho del evento en horas. Un evento lo vamos a definir como las 9 horas anteriores y las 3 horas
dt_end = timedelta(hours=9.0)  #posteriores a la hora del acumulado (-3horas)

conf=dict()
conf['latradar']=-34.787778
conf['lonradar']=-58.536667
conf['altradar']=30
conf['radius']=10000.0          #Este es el radio del cilindro que se toma alrededor de la estacion y sobre el cual se promedia la reflectividad.
conf['min_ref_profile']=5.0     #Los valores de reflecividad por debajo de este valor no van a ser considerados en el calculo del perfil vertical de reflectividad.
conf['min_ref_etop']=10.0       #Este es el valor de reflectividad usado para el calculo del echotop (y del VILD)
conf['z_min']=0.0               #La altura minima donde arranca el perfil.
conf['z_max']=15000.0           #La altura maxima donde termina el perfil.
conf['delta_z']=500.0           #La resolucion vertical del perfil.
conf['undef']=-32.0             #Valor de reflectividad que indica la falta de dato. 
conf['vil_threshold']=5.0       #Valor de VIL a partir del cual los perfiles son incorporados. 


#Cargo la base de datos de las estaciones.
station_dict = rpm.get_stations( station_database_file )

#Carglo la base de datos de eventos.
cases_dict = rpm.get_cases( event_list_file )

print('ID, DATE , NFILES')
for ievent , my_event in enumerate( cases_dict['date'] ) :

   #Para cada evento definimos la fecha de inicio y de fin.    
   event_date = cases_dict['date'][ievent]
   event_ini_date = dt.strftime( event_date - dt_ini , '%Y%m%d%H' )
   event_end_date = dt.strftime( event_date + dt_end , '%Y%m%d%H' )
   
   #Busco la lista de archivos que necesito leer para este evento.
   conf['filelist'] = rpm.get_filelist( event_ini_date , event_end_date , radar_data_path , name_filter=['cz240p1','cz240p0'] , sufix='.z.rvd' , data_order='RVD' ) 

   print( cases_dict['ID'][ievent] + ' , ' + dt.strftime( cases_dict['date'][ievent] , '%Y%m%d%H' ) + ' , ' + str(len( conf['filelist']  )) )

   if len( conf['filelist'] ) > 0 :
      #Identifico la longitud y latitud del evento.
      conf['lon'] = float( station_dict[ cases_dict['ID'][ievent] ]['lon'] )
      conf['lat'] = float( station_dict[ cases_dict['ID'][ievent] ]['lat'] )

      #Caclulo el perfil vertical de la reflectividad en base a la lista que encontre.
      event_profiles = rpm.get_profiles( conf )
   
      #Agregamos algo de metadata al diccionario que contiene el perfil promediado sobre el cilindro.
      event_profiles['ini_date'] = event_ini_date
      event_profiles['end_date'] = event_end_date
      event_profiles['pp6h_obs_date'] = event_date
      event_profiles['pp6h'] = cases_dict['pp6h'][ievent]
      event_profiles['ID'] = cases_dict['ID'][ievent]
      event_profiles['name'] = station_dict[ cases_dict['ID'][ievent] ]['name']

      #Para ahorrar espacio en disco.
      del event_profiles['z_raw_profile']
      del event_profiles['alt_raw_profile']
      del event_profiles['elev_raw_profile']

      #Guardamos el perfil medio sobre el cilindro. 
      pickle_file = output_path + '/event_profiles_' + cases_dict['ID'][ievent] + '_' + dt.strftime( cases_dict['date'][ievent] , '%Y%m%d%H' ) + '.pkl'
      pkl.dump( event_profiles , open( pickle_file , "wb" ) )


   #event_th_nn_profile = ..... #Este seria el perfil del vecino mas cercano (nn).    

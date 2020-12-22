
# =============================================================================

import numpy as np
import zlib
import io
import pyart
#import matplotlib.pyplot as plt

shift_az=7.4999   #No usar 7.5 porque genera un problema con la founcion np.round()
shift_date = '20150420'

def rvd_read( input_file , radar_lon , radar_lat , radar_alt , correct_az = True , correct_elev = True ) :

   my_file = open(input_file,'rb')
   header_bufr=io.BytesIO( my_file.read(512*4) )
   data_bufr = my_file.read()
   my_file.close()

   range_res = 0.5  #Los datos de salida van a tener esta resolucion.
   max_nrange = 480  #Maxima cantidad de pixeles.

   #Get the data
   data = np.frombuffer( zlib.decompress( data_bufr , 16+zlib.MAX_WBITS ) , dtype=np.uint8) / 2 -32
   #Get the header
   header=dict()
   header['VERS']=np.array(header_bufr.readline()[6:]).astype(float)
   header['RADAR']=header_bufr.readline().decode('utf-8')[7:].rstrip(' \n')
   header['DATE']=header_bufr.readline().decode('utf-8')[6:].rstrip(' \n')

   header['MOMENT']=header_bufr.readline().decode('utf-8')[8:].rstrip(' \n')
   header['JOB_NAME']=header_bufr.readline().decode('utf-8')[10:].rstrip(' \n')
   header['NELEVS']=np.int(header_bufr.readline()[8:])
   header_bufr.readline()
   header['ELEVS']=np.zeros(header['NELEVS'])
   header['NRANGE']=np.zeros(header['NELEVS']).astype(int)
   header['RANGE_RES']=np.zeros(header['NELEVS'])
   header['NAZIMUTH']=np.zeros(header['NELEVS']).astype(int)
   header['V_MAX']=np.zeros(header['NELEVS'])
   for ielev in range( header['NELEVS']) :
      tmp_par = header_bufr.readline().decode('utf-8').rstrip('\n').split(':')
      #print(tmp_par)
      header['ELEVS'][ielev] = np.array( tmp_par[1] ).astype(float)
      header['NRANGE'][ielev] = np.array( tmp_par[2] ).astype(int)
      header['RANGE_RES'][ielev] = np.array( tmp_par[3] ).astype(float)
      header['NAZIMUTH'][ielev] = np.array( tmp_par[4] ).astype(int)
      header['V_MAX'][ielev] = np.array( tmp_par[5] ).astype(float)

   nazimuth = np.unique( header['NAZIMUTH'])
   if ( nazimuth.size > 1 ):
      print('Error: Tenemos diferentes numeros de azimuths!')


   # Make a empty radar with the dimensions of the dataset.
   radar = pyart.testing.make_empty_ppi_radar( max_nrange , nazimuth , header['NELEVS'])

   if correct_elev :
      if header['NELEVS']== 12 :
         header['ELEVS'] = np.array([0.3,0.5,0.9,1.3,1.8,2.3,3.1,4.0,5.1,6.3,8.0,10.0])
      elif header['NELEVS']== 14 :
         header['ELEVS'] = np.array([0.3,0.5,0.9,1.3,1.8,2.3,3.1,4.0,5.1,6.3,8.0,10.0,13.0,19.0])
      elif header['NELEVS']== 15 :
         header['ELEVS'] = np.array([0.4,0.9,1.3,1.7,2.2,2.8,3.5,4.4,5.5,6.8,8.4,10.3,13.4,19.2,34.3])
      else  :
         if header['ELEVS'][0] == header['ELEVS'][1] :
            header['ELEVS'][0]=0.3 

   # Start filling the radar attributes with variables in the dataset.

   radar.latitude['data']     = np.array([radar_lat])
   radar.longitude['data']    = np.array([radar_lon])
   radar.altitude['data']     = np.array([radar_alt])

   radar.range['data']        = range_res*np.arange(0,max_nrange)*1e3
   radar.fixed_angle['data']  = header['ELEVS']
   radar.sweep_number['data'] = header['NELEVS']

   radar.azimuth['data']=np.zeros( nazimuth * header['NELEVS'] )
   radar.elevation['data']=np.zeros( nazimuth * header['NELEVS'] )
   radar.time['data']=np.zeros( nazimuth * header['NELEVS'] )
   radar.ray_angle_res=dict()
   radar.ray_angle_res['data'] = 360.0 / nazimuth.astype(float)

   ref=np.nan * np.zeros((nazimuth[0]*header['NELEVS'],max_nrange))
   position=0
   ray=0

   for el in range(header['NELEVS']) :
      if  header['RANGE_RES'][el] == range_res    :
         my_range= int( min([ header['NRANGE'][el] , max_nrange ]) )
         for az in range(nazimuth[0]) :
            ref[ray,0:my_range]=data[ position:position+my_range ]
            position = position + header['NRANGE'][el]
            radar.azimuth['data'][ray] = az * radar.ray_angle_res['data']
            radar.elevation['data'][ray] = radar.fixed_angle['data'][el]
            ray=ray+1
      if  header['RANGE_RES'][el] == range_res / 2.0   :
         my_range= int( min([ header['NRANGE'][el] , max_nrange * 2 ]) / 2 )
         for az in range(nazimuth[0]) :
            ref[ray,0:my_range]=running_mean(data[ position:position+ (my_range*2) ], 2)[::2]
            position = position + header['NRANGE'][el]
            radar.azimuth['data'][ray] = az * radar.ray_angle_res['data']
            radar.elevation['data'][ray] = radar.fixed_angle['data'][el]
            ray=ray+1
            
            
   #Rotate azimuth for dates before shift_date 
   if correct_az :
      cdate = header['DATE']
      if float(cdate[0:8]) < float(shift_date)   :
         radar.azimuth['data'] = radar.azimuth['data'] + shift_az 
         radar.azimuth['data'][radar.azimuth['data']>360.0]=radar.azimuth['data'][radar.azimuth['data']>360.0]-360.0

   # Let's work on the field data, we will just do reflectivity for now, but any of the
   # other fields can be done the same way and added as a key pair in the fields dict.
   ref_dict = pyart.config.get_metadata('reflectivity')
   ref_dict['data'] = np.array(ref)
   radar.fields = {'reflectivity': ref_dict}
   
   radar.metadata['instrument_name']='DWSR-2500C'
   
   radar.metadata['start_datetime']=cdate[0:4]+'-'+cdate[4:6]+'-'+cdate[6:8]+'T'+cdate[8:10]+':'+cdate[10:12]+':00Z'
   radar.metadata['start_time']=cdate[0:4]+'-'+cdate[4:6]+'-'+cdate[6:8]+' '+cdate[8:10]+':'+cdate[10:12]+':00.000'
   radar.metadata['end_datetime']=radar.metadata['start_datetime']
   radar.metadata['end_time']=radar.metadata['start_time']

   return radar


def order_variable ( radar , var_name , undef )  :

   import numpy as np
   import numpy.ma as ma

   #From azimuth , range -> azimuth , range , elevation

   if radar.ray_angle_res != None   :
      #print( radar.ray_angle_res , radar.ray_angle_res == None )
      ray_angle_res = np.unique( radar.ray_angle_res['data'] )
   else                             :
      print('Warning: ray_angle_res no esta definido, estimo la resolucion en radio como la diferencia entre los primeros angulos')
      ray_angle_res = np.min( np.abs( radar.azimuth['data'][1:] - radar.azimuth['data'][0:-1] ) )
      print('La resolucion en rango estimada es: ',ray_angle_res)


   if( np.size( ray_angle_res ) >= 2 )  :
      print('Warning: La resolucion en azimuth no es uniforme en los diferentes angulos de elevacion ')
      print('Warning: El codigo no esta preparado para considerar este caso y puede producir efectos indeseados ')
   ray_angle_res=np.nanmean( ray_angle_res )

   levels=np.sort( np.unique(radar.elevation['data']) )
   nb=radar.azimuth['data'].shape[0]

   order_azimuth=np.arange(0.0,360.0,ray_angle_res) #Asuming a regular azimuth grid

   na=np.size(order_azimuth)
   ne=np.size(levels)
   nr=np.size(radar.range['data'].data)


   var = np.ones( (nb,nr) )

   if ( var_name == 'altitude' ) :
       var[:]=radar.gate_altitude['data']
   elif( var_name == 'longitude' ) :
       var[:]=radar.gate_longitude['data']
   elif( var_name == 'latitude'  ) :
       var[:]=radar.gate_latitude['data']
   elif( var_name == 'x' )         :
       var[:]=radar.gate_x['data']
   elif( var_name == 'y' )         :
       var[:]=radar.gate_y['data']
   else  :
       var[:]=radar.fields[var_name]['data'].data


   #Allocate arrays
   order_var    =np.zeros((na,nr,ne))
   order_time   =np.zeros((na,ne))
   azimuth_exact=np.zeros((na,ne))
   order_n      =np.zeros((na,nr,ne),dtype='int')

   current_lev = radar.elevation['data'][0]
   ilev = np.where( levels == current_lev )[0]

   for iray in range( 0 , nb )  :   #Loop over all the rays

     #Check if we are in the same elevation.
     if  radar.elevation['data'][iray] != current_lev  :
         current_lev = radar.elevation['data'][iray]
         ilev=np.where( levels == current_lev  )[0]

     #Compute the corresponding azimuth index.
     az_index = np.round( radar.azimuth['data'][iray] / ray_angle_res ).astype(int)
     #Consider the case when azimuth is larger than na*ray_angle_res-(ray_angle_res/2)
     if az_index >= na   :
        az_index = 0

     tmp_var = var[iray,:]
     undef_mask = tmp_var == undef
     tmp_var[ undef_mask ] = 0.0

     order_var [ az_index , : , ilev ] = order_var [ az_index , : , ilev ] + tmp_var
     order_n   [ az_index , : , ilev ] = order_n   [ az_index , : , ilev ] + np.logical_not(undef_mask).astype(int)

     order_time[ az_index , ilev ] = order_time[ az_index , ilev ] + radar.time['data'][iray]
     azimuth_exact[ az_index , ilev ] = azimuth_exact[ az_index , ilev ] + radar.azimuth['data'][ iray ]

   order_var[ order_n > 0 ] = order_var[ order_n > 0 ] / order_n[ order_n > 0 ]
   order_var[ order_n == 0] = undef

   return order_var , order_azimuth , levels , order_time , azimuth_exact

def running_mean(x, N):
    tmp=np.zeros(np.shape(x))
    cumsum = np.cumsum(np.insert(x, 0, 0))
    tmp[0:-1]=(cumsum[N:] - cumsum[:-N]) / float(N)
    tmp[-1]=tmp[-2]
    return tmp

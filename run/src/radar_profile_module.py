#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 19:33:33 2020

@author: juan
"""
import os
import glob
import numpy as np
from src import rvd_read as rr 
from datetime import datetime as dt
from datetime import timedelta

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
      print('Reading ',my_file)

      #Esta funcion devuelve todos los puntos que estan dentro del cilindro (ref,alt,elev)
      #Tambien devuelve el perfil de vecinos mas cercanos nn_ref , nn_alt y nn_elev
      #[ref , alt , elev , date , nn_ref , nn_alt , nn_elev ] = extract_profile_data( my_file , conf ) #radius , lonp , latp , lonradar , latradar , altradar ) 

      [ref , alt , elev , nn_ref , nn_alt , nn_vil , nn_elev , date ] = extract_profile_data_interp( my_file , conf ) #radius , lonp , latp , lonradar , latradar , altradar ) 

      [zp , meanp , stdp , minp , maxp , nump, etop , vil , vild ] = grid_profile( ref , alt , conf )
      
   
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
      
      # import matplotlib.pyplot as plt
      # plt.figure(figsize=[24,8])

      # plt.plot( my_profile['z_raw_profile'][-1] , my_profile['alt_raw_profile'][-1]/1000 ,'ko',alpha=0.3)
      # plt.plot( my_profile['meanref_th_profile'][:,ifile] , my_profile['z_th_profile'][:,ifile]/1000 ,'-bo',label='mean')
      # plt.plot( my_profile['ref_nn_profile'][-1]       , my_profile['z_nn_profile'][-1]/1000 ,'-ro',label='NeNe')
      # plt.plot( my_profile['maxref_th_profile'][:,ifile]  , my_profile['z_th_profile'][:,ifile]/1000 ,'-go',label='Max')
      # plt.plot( my_profile['minref_th_profile'][:,ifile]  , my_profile['z_th_profile'][:,ifile]/1000 ,'-go',label='Min')
      # plt.plot( my_profile['meanref_th_profile'][:,ifile] + my_profile['stdref_th_profile'][:,ifile] , my_profile['z_th_profile'][:,ifile]/1000 ,'--b')
      # plt.plot( my_profile['meanref_th_profile'][:,ifile] - my_profile['stdref_th_profile'][:,ifile] , my_profile['z_th_profile'][:,ifile]/1000 ,'--b')
      # plt.show()

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
       if  tmp_num > 0 :
          meanp[ii] = 10.0*np.log10( np.mean( 10.0**(ref[my_mask]/10.0) ) )
          stdp[ii]  = np.std ( ref[my_mask] )
          minp[ii]  = np.min ( ref[my_mask] )
          maxp[ii]  = np.max ( ref[my_mask] )
          nump[ii]  = tmp_num

   #Detect the highest echo top
   etop = np.nan  #Echo top
   vil  = np.nan  #Vertically integrated liquid
   vild = np.nan  #VIL density
   meanpower = 10.0 ** ( meanp / 10.0 )
   for ii in range( nbin -1 )  :
       if ( meanp[ii+1] < min_ref_etop ) & ( meanp[ii] >= min_ref_etop ) :
          etop = zp[ii] 

       if ~ np.isnan( meanp[ii+1] ) & ~ np.isnan( meanp[ii] ) :
          tmp_z_mean = 0.5 * (meanpower[ii+1] + meanpower[ii])
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


def extract_profile_data_interp( filename , conf ) : #radius , lonp , latp , lonradar , latradar , altradar ) :
    #Esta version de la funcion primero interpola todos los datos a la posicion de los 
    #pixeles en el nivel mas bajo. Luego define ahi el punto mas cercano y los puntos que caen en el area del 
    #cilindro. Tambien usa el VIL para determinar que perfiles se retienen y cuales no. 
    undef = conf['undef']
    
    radar = rr.rvd_read( filename , conf['lonradar'] , conf['latradar'] , conf['altradar'] )
        
    [dbz3d , azimuth , levels , time , azimuthe , npixdbz ] = order_variable ( radar , 'reflectivity' , undef )  
    [lon3d , azimuth , levels , time , azimuthe , npix ] = order_variable ( radar , 'longitude' , undef )  
    [lat3d , azimuth , levels , time , azimuthe , npix ] = order_variable ( radar , 'latitude' , undef ) 
    [z3d , azimuth , levels , time , azimuthe , npix ]   = order_variable ( radar , 'altitude' , undef )
    
    #Reemplazamos los nans que puedan haber en la reflectividad por undef.
    dbz3d[np.isnan(dbz3d)] = undef 
    dbz3d[dbz3d == 0.0 ]   = undef
        
    #Dbz_int y z_int es la reflectividad y la altura interpoladas a la reticula del angulo de elevacion
    #mas bajo. De manera que ahora los puntos correspondientes a las diferentes elevaciones son coincidentes
    #en la vertical. 
    tmp_lon = lon3d - conf['lonradar'] #Estas variables sirven para tener un sistema de referencia centrado en el radar.
    tmp_lat = lat3d - conf['latradar']
    [vil_int , etop_int , dbz_int , z_int ] = calcula_vil( dbz3d , tmp_lon , tmp_lat , z3d , undef )  
    
    #print( np.sum( radar.fields['reflectivity']['data'] == 0 ),np.sum( dbz3d==0.0 ) )
    #import matplotlib.pyplot as plt
    #print( np.max(vil_int[:,0:200]) , np.min(vil_int) )
    #plt.pcolor(vil_int);plt.colorbar();plt.show()
    #print( np.max(etop_int) , np.min(etop_int) )
    #plt.pcolor(etop_int[:,0:200]);plt.colorbar();plt.show()
    #plt.pcolor(dbz_int[210,0:200,:]);plt.contour(z_int[210,0:200,:],levels=[5000,10000,20000]);plt.colorbar();plt.show()
    #plt.pcolor(dbz3d[210,0:200,:]);plt.contour(z3d[210,0:200,:],levels=[5000,10000,20000]);plt.colorbar();plt.show()
    #plt.pcolor(z_int[210,0:200,:]);plt.colorbar();plt.show()

    #print(np.sum(np.abs( dbz3d ) < 9.99997940e-10) , np.sum(np.abs( radar.fields['reflectivity']['data'] ) < 9.99997940e-10) ) 


    # Buscamos los puntos que estan en el cilindro y calculamos el perfil medio sobre el cilindro.    
    dlon =  np.cos( lat3d[:,:,0]*np.pi/180.0)*( lon3d[:,:,0] - conf['lon'] )
    dlat =  ( lat3d[:,:,0] - conf['lat'] )
    distancia = np.sqrt( ( dlon * 111000.0 )**2 + ( dlat * 111000.0 )**2 )
    date = radar.metadata['start_datetime']

    mascara = np.logical_and( distancia <= conf['radius'] , vil_int > conf['vil_threshold'] ) 
    mascara = np.logical_and( mascara , etop_int > conf['etop_threshold'] )
    mascara = np.logical_and( mascara , dbz_int[:,:,0] > conf['lowref_threshold'] )
    #mascara =  distancia <= conf['radius'] 

    mascara3d = np.repeat( mascara[:,:,np.newaxis],np.shape(dbz_int)[2],axis=2)
    lev3d   = np.repeat( np.repeat( levels[np.newaxis,:],np.shape(dbz_int)[0],axis=0)[:,np.newaxis,:] ,np.shape(dbz_int)[1],axis=1)

                    
    ref  = dbz_int[ mascara3d ]
    elev = lev3d[ mascara3d ]
    alt  = z_int[ mascara3d ]
    #vil  = vil_int[ mascara ]

    #Buscamos los puntos que son el vecino mas cercano (nn) del punto seleccionado en cada ppi para
    #construir el perfil de vecinos mas cercanos. 

    [minx , miny] = np.where( distancia == np.min( distancia ) )
    
    ref_nn = np.copy( dbz_int[minx,miny,:] )
    alt_nn = np.copy( z_int[minx,miny,:] )
    vil_nn = np.copy( vil_int[minx,miny] )

    return ref , alt , elev , ref_nn , alt_nn , vil_nn , levels , date


def local_mean( array , kernel_x , kernel_y , undef ) :
    #Asumimos que hay condiciones ciclicas en el axis 0 pero no en el 1.
    #array es el array de datos de entrada
    #kernel_x es cual es el desplazamiento maximo (hacia cada lado) en la direccion de x
    #kernel_y es cual es el desplazamiento maximo (hacia cada lado) en la direccion de y
    #undef son los valores invalidos en el array de entrada.

    [nx,ny]=np.shape(array)
    arraym = np.zeros( np.shape(array) )
    countm = np.zeros( np.shape(array) )
    for ix in range(-kernel_x,kernel_x+1) :
        for iy in range(-kernel_y,kernel_y +1) :
          tmp_array = np.zeros( np.shape(array) )
          if iy > 0 :
             tmp_array[:,0+iy:] = array[:,0:-iy]
          if iy == 0 :
             tmp_array = np.copy(array)
          if iy < 0 :
             tmp_array[:,0:iy] = array[:,-iy:]
          tmp_array=np.roll( tmp_array , ix , axis=0 )
          mask = tmp_array != undef
          arraym[ mask ] = arraym[mask] + tmp_array[mask]
          countm[ mask ] = countm[mask] + 1
    mask = countm > 0
    arraym[mask] = arraym[mask] / countm[mask]
    arraym[~mask] = undef 

    return arraym


def calcula_vil( dbz_in , x_in , y_in , z_in , undef , etop_thresh=5.0 )  :
  from scipy.interpolate import interp1d
  dbz = np.copy(dbz_in)
  x   = np.copy(x_in)
  y   = np.copy(y_in)
  z   = np.copy(z_in)
  
  fill_value_power = 1.0e-10

  [na,nr,ne] = dbz.shape

  dbz_int = np.zeros( dbz.shape )
  z_int   = np.zeros( dbz.shape )
  vil_int     = np.zeros( (na , nr) )
  etop_int    = np.zeros( (na , nr) )
  
  etop_init_mask = np.zeros( (na , nr) ).astype(bool)
  etop_detected_mask = np.zeros( (na , nr) ).astype(bool)
  
  
  

  ranger = ( x**2 + y**2 )**0.5
  ranger0 = ranger[:,:,0]
  
  #Calculo el VIL en la reticula x0 , y0
  for ie in range(ne)   :
    dbz2d = np.copy( dbz[:,:,ie] )
    dbz2d_mean = local_mean( dbz2d , 1 , 1 , undef )
    #Intento salvar algunos agujeros que pueda haber en el campo de reflectividad.
    mask = np.logical_or( dbz2d == undef , dbz2d < 0.0 )
    dbz2d[ mask ] = dbz2d_mean[mask]
    #Los undef que quedaron pasan a ser 0 para el calcuo del VIL 
    #dbz2d[dbz2d == undef ] = fill_value_power  
    dbz2d = 10.0 ** (  dbz2d / 10.0 )
    dbz2d[ dbz[:,:,ie] == undef ] = fill_value_power
              
    for ia in range(na)   :
      
      interpolator = interp1d(ranger[ia,:,ie] , dbz2d[ia,:] , kind='linear' , bounds_error = False , fill_value = fill_value_power )
      dbz_int[ia,:,ie] = interpolator(ranger0[ia,:])
      interpolator = interp1d(ranger[ia,:,ie] , z[ia,:,ie] , kind='linear' , bounds_error = False , fill_value = np.nan)
      z_int[ia,:,ie] = interpolator(ranger0[ia,:])
      #Completo algunos niveles repitiendo el ultimo valor para hacer mas robusto el calculo del VIL
    if ie > 0 :
      dz = z_int[:,:,ie] - z_int[:,:,ie-1] ; dz[dz==0] = np.nan
      vil_inc = 3.44e-6 * ( ( 0.5*(dbz_int[:,:,ie] + dbz_int[:,:,ie-1]) ) ** (4.0/7.0) ) * ( dz )
      vil_inc[np.isnan(vil_inc)] = 0.0 
      vil_int = vil_int + vil_inc

    #Echo top computation      
    etop_init_mask[ np.logical_and( dbz[:,:,ie] > etop_thresh , np.logical_not( etop_detected_mask ) ) ]=True
    
    if ie > 0 :
        etop_detection_mask=np.logical_and( etop_init_mask , np.logical_or( dbz[:,:,ie] < etop_thresh , np.isnan(dbz[:,:,ie] ) ) )
        etop_detected_mask[ etop_detection_mask ] = True
        etop_int[ etop_detection_mask ] = ( z_int[:,:,ie] )[ etop_detection_mask ]
        etop_init_mask[ etop_detected_mask ] = False
    if ie == ne-1 :
        etop_detection_mask=np.logical_and( etop_init_mask , dbz[:,:,ie] > etop_thresh )
        etop_int[ etop_detection_mask ] = ( z_int[:,:,ie] )[ etop_detection_mask ]
      
      
      
  #Hasta aca tenemos vil_int que es el vil en la reticula x0 y0. Para las cuentas en general nos puede venir bien
  #tener el vil interpolado a la reticula x,y (es decir un vil definido para todoas las elevaciones del radar)
  vil_int[ np.isnan(vil_int) ] = 0.0
  dbz_int_out = 10.0*np.log10(dbz_int)
  dbz_int_out[ dbz_int == fill_value_power ] = undef 
  
  
     
  return vil_int , etop_int , dbz_int_out , z_int  


def var_int( var_in , x_in , y_in , int_lev = 0 , fill_value = 0.0 ) :
    from scipy.interpolate import interp1d
    var = np.copy(var_in)
    x   = np.copy(x_in)
    y   = np.copy(y_in)

    [na,nr,ne] = var.shape

    var_int = np.zeros( var.shape )
    ranger = ( x**2 + y**2 )**0.5
    ranger0 = ranger[:,:,int_lev]

    for ie in range(ne)   :
       for ia in range(na)   :
          interpolator = interp1d(ranger[ia,:,ie] , var[ia,:,ie] , kind='linear' , bounds_error = False , fill_value = 0.0 )
          var_int[ia,:,ie] = interpolator(ranger0[ia,:])

    return var_int 


def order_variable ( radar , var_name , undef )  :  

   import numpy as np
   #import warnings 
   #import matplotlib.pyplot as plt

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

     tmp_var = np.copy(var[iray,:])
     undef_mask = tmp_var == undef 
     tmp_var[ undef_mask ] = 0.0
    
     order_var [ az_index , : , ilev ] = order_var [ az_index , : , ilev ] + tmp_var
     order_n   [ az_index , : , ilev ] = order_n   [ az_index , : , ilev ] + np.logical_not(undef_mask).astype(int)

     order_time[ az_index , ilev ] = order_time[ az_index , ilev ] + radar.time['data'][iray]
     azimuth_exact[ az_index , ilev ] = azimuth_exact[ az_index , ilev ] + radar.azimuth['data'][ iray ]

   order_var[ order_n > 0 ] = order_var[ order_n > 0 ] / order_n[ order_n > 0 ]
   order_var[ order_n == 0] = undef

   return order_var , order_azimuth , levels , order_time , azimuth_exact , order_n 

def order_variable_inv (  radar , var , undef )  :

   import numpy as np
   
   #From azimuth , range , elevation -> azimuth , range

   na=var.shape[0]
   nr=var.shape[1]
   ne=var.shape[2]

   nb=radar.azimuth['data'].shape[0]

   levels=np.sort( np.unique(radar.elevation['data']) )

   if radar.ray_angle_res != None   :
      #print( radar.ray_angle_res , radar.ray_angle_res == None )
      ray_angle_res = np.unique( radar.ray_angle_res['data'] )
   else                             :
      print('Warning: ray_angle_res no esta definido, estimo la resolucion en radio como la diferencia entre los primeros angulos')
      ray_angle_res = np.min( np.abs( radar.azimuth['data'][1:] - radar.azimuth['data'][0:-1] ) )
      print('La resolucion en rango estimada es: ',ray_angle_res)

   if( np.size( ray_angle_res ) >= 2 )  :
      print('Warning: La resolucion en azimuth no es uniforme en los diferentes angulos de elevacion ')
      print('Warning: El codigo no esta preparado para considerar este caso y puede producir efectos indesaedos ')
   ray_angle_res=np.nanmean( ray_angle_res )

   current_lev = radar.elevation['data'][0]
   ilev = np.where( levels == current_lev  )[0]

   output_var = np.zeros((nb,nr) )
   output_var[:] = undef

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

      output_var[ iray , : ] = var[ az_index , : , ilev ]

   return output_var





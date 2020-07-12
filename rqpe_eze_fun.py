import glob
import os
import numpy as np
import rvd_read as rr
from datetime import datetime , timedelta
import gc 
import pickle as pkl

#TODO LIST:
# -Hacer una pequeña estadistica de la covertura temporal de los datos en un periodo (dt medio entre imagenes, maximo periodo sin imagenes)
# -Una funcion para acumular precipitacion estimada por radar en periodos de 6, 12, etc para comparar con el pluviometro.
# -La funcion que extrae el perfil vertical de la reflectividad para un punto sobre un periodo determinado de tiempo. 
# -Una funcion que permita calibrar lo que ve el radar con los acumulados a 24 horas.


def rqpe_eze_fun(lat_est,lon_est,id_est,fecha_ini,fecha_fin,path,ppi_pp=0,save=True,out_path=None):


    
    # Radar EZE
    radar_lat = -34.787778
    radar_lon = -58.536667
    radar_alt = 30 # m

    if out_path is None :
        out_path = './'

    fecha_tmp=datetime.strptime( fecha_ini[0:8] , '%Y%m%d')
    fecha_fin_tmp=datetime.strptime( fecha_fin[0:8] , '%Y%m%d')

    dt = timedelta( days = 1.0 )

    FilelistRVD=[]
    while( fecha_tmp <= fecha_fin_tmp ) : 
       
       DataPath = path + '/' + fecha_tmp.strftime('%Y%m%d') + '/'
       FilelistRVD +=  glob.glob(DataPath+'ar1.cz240*.z.rvd') 

       fecha_tmp = fecha_tmp + dt
    FilelistRVD.sort()

    fecha_ini_tmp=datetime.strptime( fecha_ini , '%Y%m%d%H%M')
    fecha_fin_tmp=datetime.strptime( fecha_fin , '%Y%m%d%H%M')
    to_be_removed = []
    for ifile , my_file in enumerate( FilelistRVD ) :
        filename = os.path.basename( my_file )
        filedate = datetime.strptime( filename[12:20] + filename[21:25] ,'%Y%m%d%H%M') 
        if( filedate < fecha_ini_tmp or filedate > fecha_fin_tmp ) :
            to_be_removed.append( my_file )
    for my_file in to_be_removed :
        FilelistRVD.remove( my_file )

    vector_hora = []

    pac = 0
    
    #Initialize output dictionary.
    nfiles = len( FilelistRVD )
    output=dict()
    output['Date']=list()
    output['Time']=list()
    output['MRef']=np.zeros( nfiles ) + np.nan
    output['Ref']=np.zeros( nfiles ) + np.nan
    output['RR']=np.zeros( nfiles ) + np.nan
    output['PAC']=np.zeros( nfiles ) + np.nan
    
    for ifile,fname in enumerate(FilelistRVD) :

        #print(fname)

        radar = rr.rvd_read( fname , radar_lon , radar_lat , radar_alt )

        #====================================================================================#

        t = datetime.strptime(fname[len(DataPath)+21:-6], "%H%M")
        hour = int(t.strftime('%H%M')[0:2])
        mins = int(t.strftime('%H%M')[2:4])/60
        vector_hora.append(hour+mins)
        
        elevs=np.unique( radar.elevation['data'] )
        if ppi_pp <= np.size(elevs) :
           clev = elevs[ppi_pp]
        else                       :
           clev = elevs[-1]
           print('Warning: PPI_PP > NELEVS')
        data_mask = radar.elevation['data'] == clev 

        # Genero un subset de los datos del PPI mas bajo:
        lat_ppi = radar.gate_latitude['data'][data_mask,:]
        lon_ppi = radar.gate_longitude['data'][data_mask,:]
        ref_ppi = radar.fields['reflectivity']['data'][data_mask,:]
        #ref_ppi = np.roll(radar.fields['reflectivity']['data'][0:359],7,axis=0)

        # Calculo la reflectividad media alrededor del radar para identificar momentos de atenuacion
        ref_ppi_linear = np.power(10.0,ref_ppi/10.0)
        ref_ppi_linear_mean = np.mean(ref_ppi_linear)
        ref_ppi_mean = 10*np.log10(ref_ppi_linear_mean)

        # Buscamos el pixel mas cercano a la estacion
        dlon = np.cos(lat_ppi*np.pi/180.0)*( lon_ppi - lon_est )
        dlat = ( lat_ppi - lat_est )
        dist = np.sqrt( ( dlon * 111000.0 )**2 + ( dlat * 111000.0 )**2 )

        mascara = dist == dist.min()

        ref_est = ref_ppi[ mascara ][0]

        dBZthrMin=5.0
        dBZthrHail=55.0
        if ref_est>dBZthrHail:
            ref_est=dBZthrHail

        dBZbaselinear = np.power(10.0,ref_est/10.0)
        
        if ref_est < dBZthrMin :
            dBZbaselinear = 0.0 
        
        # Calculo tasa instantanea de precipitacion usando parámetro de la relación Marshall-Palmer:
        a=200.0
        b=1.6
        rr_est = np.power(dBZbaselinear/a,1/b) # En mm/hr
        pac_est = rr_est*(10/60)
        pac = pac + pac_est

        print( fname[len(DataPath)+12:-6]+'Z', str(ref_ppi_mean)+'dBZ', str(ref_est)+'dBZ', str(rr_est)+'mm/h')
        
        output['Date'].append(fname[len(DataPath)+12:-11])
        output['Time'].append(fname[len(DataPath)+21:-6])
        output['MRef'][ifile] = ref_ppi_mean
        output['Ref'][ifile] = ref_est
        output['RR'][ifile] = rr_est
        output['PAC'][ifile] = pac

        #====================================================================================#
        gc.collect()


    out_file = out_path + '/EZE_' + str(id_est) + '_' + fecha_ini + '_' + fecha_fin +  '.pkl'
    pkl.dump( output , open( out_file , "wb" ) )

    return output, vector_hora

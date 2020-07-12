import glob
import numpy as np
import matplotlib.pyplot as plt
import rqpe_eze_fun as rqpe

# ARGUMENTOS DE ENTRADA DE LA FUNCION
lat_est = -34.59022 # >>> Latitud del punto donde queremos la serie temporal
lon_est = -58.48397 # >>> Longitud del punto donde queremos la serie temporal
id_est  = '87SARASA' 
fecha_ini = '201309070400'  # >>> Fecha
fecha_fin = '201309071200'
path='/media/juan/Seagate Backup Plus Drive/data/radar/EZE/'

# Calculamos la serie temporal de precipitacion
output, vector_hora = rqpe.rqpe_eze_fun(lat_est,lon_est,id_est,fecha_ini,fecha_fin,path)

#==============================================================================#

fig=plt.figure(figsize=[15,5])
ax = fig.add_subplot(111)

ax.plot(vector_hora,output['RR'],'-bo')
ax.set_xlabel('Time [UTC]', fontsize=14)
ax.set_ylabel('Rain Rate [mm/h]', color='b', fontsize=14)
#ax.set_xlim([0,24])
ax.set_ylim([0,20])
#ax.set_xticks(np.arange(0, 25, step=1))
ax.set_yticks(np.arange(0, 20.1, step=2))

ax2 = ax.twinx()
#ax2.plot(vector_hora,output['MRef'],'--r')
ax2.plot(vector_hora,output['PAC'],'--r')
ax2.set_xlabel('Time [UTC]', fontsize=14)
ax2.set_ylabel('PAC [mm]', color='r', fontsize=14)
#ax2.set_ylabel('Mean Reflectivity [dBZ]', fontsize=14)
#ax2.set_xlim([0,24])
ax2.set_ylim([0,30])
ax2.set_yticks(np.arange(0, 30.1, step=2))

plt.title('RQPE-OCBA   Date: '+fecha_ini, fontsize=20)

fig.savefig('RQPE_OCBA_'+fecha_ini+'.png', bbox_inches='tight', dpi=100)

plt.show()

#==============================================================================#

from rediscluster import RedisCluster
#import redis
import json
import numpy as np
from datetime import datetime,timedelta
import cal_indicator as indicator
import testplotmulti
import os
import time
import configparser
from mpi4py import MPI 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.pyplot import MultipleLocator
from mpl_toolkits.axisartist.parasite_axes import HostAxes, ParasiteAxes
import sys




comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
##-----------------------------------

pic_path=sys.argv[1]

##################start_change #########start_change #############start_change ###########start_change#########



cf = configparser.ConfigParser()
cf.read("../namelist/initial_for_obstat.input")
casename = eval(cf.get("obstat_ini", "case_name")) #'case_yz'
#print(casename)


#-------------------------------set redis-------------------------------------


redisaddress = eval(cf.get("obstat_ini", "redis_address")).split(',')
redisaddress = [x.split(':') for x in redisaddress]
redisaddress = np.array(redisaddress)

startup_nodes = []
address_id = 0
for x in redisaddress:
   host = redisaddress[address_id,0]
   port = redisaddress[address_id,1]
   redis_dictionary = {"host":host, "port":6379}
   startup_nodes.append(redis_dictionary)
   address_id = address_id+1

#print(startup_nodes)

#startup_nodes = [                           #集群节点及端口
    #{"host":"192.168.1.1", "port":6379},
    #{"host":"192.168.1.2", "port":6379},
    #{"host":"192.168.1.3", "port":6379}
#]

#print(startup_nodes)

conn = RedisCluster(startup_nodes=startup_nodes, decode_responses=True)



#------------------------------set time-------------------------------------

redis_da_start_time = conn.hget(casename,'da_start_time')
da_start0,da_start1,da_start2,da_start3,da_start4 = redis_da_start_time.split(',',4)
da_start_time = datetime(int(da_start0),int(da_start1),int(da_start2),int(da_start3),int(da_start4))

delta_seconds = int(conn.hget(casename,'da_in_seconds')) #5*60

start_time_get = eval(cf.get("obstat_ini", "obsstat_start_time"))
start0,start1,start2,start3,start4 = start_time_get.split(',',4)
start_time = datetime(int(start0),int(start1),int(start2),int(start3),int(start4))

end_time_get = eval(cf.get("obstat_ini", "obsstat_end_time"))
end0,end1,end2,end3,end4 = end_time_get.split(',',4)
end_time = datetime(int(end0),int(end1),int(end2),int(end3),int(end4))
run_end_time = end_time.strftime("%Y%m%d%H%M")

Numtime_to_start = int((start_time-da_start_time).total_seconds() // delta_seconds )

if (start_time-da_start_time).total_seconds() % delta_seconds == 0:
    run_start_time = start_time
else:
    run_start_time = da_start_time
    for i in range(0,Numtime_to_start+1):
        run_start_time = run_start_time + timedelta(seconds = delta_seconds)

Numtime = int(((end_time-run_start_time).total_seconds()) // delta_seconds + 1) 


print('start time:',run_start_time)
print('numtime:',Numtime)

#------------------------------set space------------------------------------

NumYGrid = int(conn.hget(casename,'num_lat')) #180
NumXGrid = int(conn.hget(casename,'num_lon')) #360
NumZGrid = int(conn.hget(casename,'mpas_num_lev')) #26
delta_lev_id = int(conn.hget(casename,'obs_stat_interval_vertical')) #4

vari_all = eval(cf.get("obstat_ini", "obsstat_var_name")).split(',')

prob_vari_surf = ['ps']

vari_surf = [x for x in vari_all if x in prob_vari_surf]
vari_norm = [x for x in vari_all if x not in prob_vari_surf]

#vari_norm = [vari_all[0],vari_all[1],vari_all[2]] #['u','v','t']
#vari_surf = [vari_all[3]] #['ps']
 
#print(vari_norm)

#print(NumZGrid)
#print(-delta_lev_id)

vari_norm_lev_id_ave = [x for x in np.arange(NumZGrid,0,-delta_lev_id)]
#vari_norm_lev_id_ave.reverse()

print(vari_norm_lev_id_ave)

#print(vari_norm_lev_id_ave)
#vari_norm_lev_id_add = [1] 
#隔几层画隔几层写,可追加指定层
vari_norm_lev_id = vari_norm_lev_id_ave#+vari_norm_lev_id_add
vari_surf_lev_id = [0]
#vari_norm_lev_id = [1,32,40,59]
#vari_surf_lev_id = [60]

onepoint_y = int(conn.hget(casename,'onepoint_y')) 
onepoint_x = int(conn.hget(casename,'onepoint_x'))
onepoint_z = int(conn.hget(casename,'onepoint_z'))

#onepoint_func = 'scatter' 
onepoint_func = 'imshow' 


#-------------------------------set hash-key name---------------------------

res_in_key = str(NumXGrid)+'x'+str(NumYGrid)
deg= eval(conn.hget(casename,'atm_mpas_sceneid'))#int(conn.hget(casename,'atm_mpas_sceneid'))
test_case = eval(conn.hget(casename,'test_case'))

cc_case = deg.split('-',2)[2]

#print(deg)

hashkey_lat = 'lat:'+res_in_key #'lat:1deg-201902'
hashkey_lon = 'lon:'+res_in_key

hashkey_tr = 'realfield:'+deg+':'+res_in_key #'realfield:1deg-201902'
hashkey_bg = 'axbfield:'+casename+':'+deg+':'+res_in_key
hashkey_amb = 'ambfield:'+casename+':'+deg+':'+res_in_key

#print(hashkey_lat)


##################end_change #########end_change #############end_change ###########end_change################


#------------------------------prepare---------------------------------------


Numlevel = len(vari_norm_lev_id)
Numproc = len(vari_norm)+len(vari_surf)

#print(Numproc)

vari = vari_norm + vari_surf

vari_lev_group_id = []
for x in vari_norm:
   vari_lev_group_id.append(vari_norm_lev_id)
for x in vari_surf:
   vari_lev_group_id.append(vari_surf_lev_id)
vari_lev_id = np.array(vari_lev_group_id)


#----------------------------------rank-------------------------------------


nl_sec = Numproc // size
nl_res = Numproc % size
if(rank < nl_res) :
    nl_start = rank * (nl_sec + 1)
    nl_end = nl_start + nl_sec + 1
else:
    nl_start = rank * nl_sec + nl_res
    nl_end = nl_start + nl_sec

print("myrank:", rank, nl_start, nl_end)


#######################################取经纬度信息################################



lat_key_get = []
for lat_in in range(1,NumYGrid+1):
        lat_key_str = str(lat_in)+':'+'lat'
        lat_key_get.append(lat_key_str)

lon_key_get = []
for lon_in in range(1,NumXGrid+1):
        lon_key_str = str(lon_in)+':'+'lon'
        lon_key_get.append(lon_key_str)

lat_get = conn.hmget(hashkey_lat, lat_key_get)
lon_get = conn.hmget(hashkey_lon, lon_key_get)

lat_get = [float(x) for x in lat_get]
lon_get = [float(x) for x in lon_get]

lat = np.array(lat_get)
lon = np.array(lon_get)



#########################################取背景场，真实场，分析场，诊断绘图##########################################

runtime = []
curr_time_from_start = run_start_time
curr_time = curr_time_from_start.strftime("%Y%m%d%H%M")

for time_in in range(0,Numtime):
    runtime.append(time_in)


    for vari_in in range(nl_start, nl_end):
        level_in_key = onepoint_z
        time_in_key = curr_time
        variable_in_key = vari[vari_in]
        level_in_key_list = vari_lev_id[vari_in]
        print(variable_in_key)

        if variable_in_key == 'ps':
            level_in_range = 1
        else:
            level_in_range = Numlevel



#---------------------------------------------剖面：刀经向---------------------------------------------------------------

        key_str = []
        key_get = []

        key_str_x = []
        key_get_x = []

        key_get_y = []

        lon_in_key = onepoint_x

        for level_in in range(0,level_in_range):
            level_in_key = level_in_key_list[level_in]

            for lat_in in range(1,NumYGrid+1):
                lat_in_key = lat_in
                key_str = str(lon_in_key)+':'+str(lat_in_key)+':'+str(level_in_key)+':'+str(variable_in_key)
                key_get.append(key_str)
                key_str_x= str(lat_in_key)+':'+'lat'
                key_get_x.append(key_str_x)
                key_get_y=key_get_y+[level_in_key]



                #print(key_get[2]) 
                #print(hashkey_tr+':'+str(time_in_key))


        lev_lat_get = conn.hmget(hashkey_amb+':'+str(time_in_key), key_get)
        lev_lat_x = conn.hmget(hashkey_lat, key_get_x)
        lev_lat_y = key_get_y


#--------------------------------------------------剖面：刀纬向--------------------------------------------------------


        key_str = []
        key_get = []

        key_str_x = []
        key_get_x = []

        key_get_y = []

        lat_in_key = onepoint_y

        for level_in in range(0,level_in_range):
            level_in_key = level_in_key_list[level_in]

            for lon_in in range(1,NumXGrid+1):
                lon_in_key = lon_in
                key_str = str(lon_in_key)+':'+str(lat_in_key)+':'+str(level_in_key)+':'+str(variable_in_key)
                key_get.append(key_str)
                key_str_x = str(lon_in_key)+':'+'lon'
                key_get_x.append(key_str_x)
                key_get_y=key_get_y+[level_in_key]
                #print(key_get[2]) 
                #print(hashkey_tr+':'+str(time_in_key))


        lev_lon_get = conn.hmget(hashkey_amb+':'+str(time_in_key), key_get)
        lev_lon_x = conn.hmget(hashkey_lon, key_get_x)
        lev_lon_y = key_get_y


#---------------------------------------------------剖面：刀水平------------------------------------------------------------


        key_str = []
        key_get = []

        if variable_in_key == 'ps':
            level_in_key = 0
        else:
            level_in_key = onepoint_z

        for lat_in in range(1,NumYGrid+1):
            lat_in_key = lat_in

            for lon_in in range(1,NumXGrid+1):
                lon_in_key = lon_in
                key_str = str(lon_in_key)+':'+str(lat_in_key)+':'+str(level_in_key)+':'+str(variable_in_key)
                key_get.append(key_str)
                #print(key_get[2]) 
                #print(hashkey_tr+':'+str(time_in_key))


        lat_lon_get = conn.hmget(hashkey_amb+':'+str(time_in_key), key_get)


#---------------------------------------------------准备剖面数据-----------------------------------------------------------


        #print(tr_get)

        lev_lat_get = [float(x) for x in lev_lat_get]
        lev_lon_get = [float(x) for x in lev_lon_get]
        lat_lon_get = [float(x) for x in lat_lon_get]

        lev_lat_x = [float(x) for x in lev_lat_x]
        lev_lat_y = [float(x) for x in lev_lat_y]
        lev_lon_x = [float(x) for x in lev_lon_x]
        lev_lon_y = [float(x) for x in lev_lon_y]

        if onepoint_func == 'scatter':
           lev_lat = np.array(lev_lat_get)
           lev_lon = np.array(lev_lon_get)
           lat_lon = np.array(lat_lon_get).reshape(NumYGrid,NumXGrid)

           lev_lat_xx = np.array(lev_lat_x)
           lev_lat_yy = np.array(lev_lat_y)
           lev_lon_xx = np.array(lev_lon_x)
           lev_lon_yy = np.array(lev_lon_y)

           if variable_in_key == 'ps':
               level_in = 0
               onepoint_z_key = 0
           else:
               level_in = vari_norm_lev_id_ave.index(onepoint_z)
               onepoint_z_key = onepoint_z

           pic_param = [curr_time, variable_in_key, level_in, onepoint_z_key, pic_path, cc_case]

           testplotmulti.plots(lev_lat_xx, lev_lat_yy, lev_lat, lev_lon_xx, lev_lon_yy, lev_lon, lat_lon, pic_param)


        elif onepoint_func == 'imshow':
           lev_lat = np.array(lev_lat_get).reshape(level_in_range,NumYGrid)
           lev_lon = np.array(lev_lon_get).reshape(level_in_range,NumXGrid)
           lat_lon = np.array(lat_lon_get).reshape(NumYGrid,NumXGrid)

           if variable_in_key == 'ps':
               level_in = 0
               onepoint_z_key = 0
           else:
               level_in = vari_norm_lev_id_ave.index(onepoint_z)
               onepoint_z_key = onepoint_z

           extent_top = NumZGrid+delta_lev_id/2

           pic_param = [curr_time, variable_in_key, level_in, onepoint_z_key, extent_top, pic_path, cc_case]

           testplotmulti.ploti(lev_lat, lev_lon, lon, lat, lat_lon, pic_param)


        else:
           print('break')


    curr_time_from_start = curr_time_from_start + timedelta(seconds = delta_seconds)
    curr_time = curr_time_from_start.strftime("%Y%m%d%H%M")



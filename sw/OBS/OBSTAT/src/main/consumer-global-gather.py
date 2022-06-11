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

obs_open = int(eval(cf.get("obstat_ini", "obs_open")))

NumYGrid = int(conn.hget(casename,'num_lat')) #180
NumXGrid = int(conn.hget(casename,'num_lon')) #360
NumZGrid = int(conn.hget(casename,'mpas_num_lev')) #26
delta_lev_id = int(conn.hget(casename,'obs_stat_interval_vertical')) #4

vari_all = eval(cf.get("obstat_ini", "obsstat_var_name")).split(',')

prob_vari_surf = ['ps']

vari_surf = [x for x in vari_all if x in prob_vari_surf]
vari_norm = [x for x in vari_all if x not in prob_vari_surf]

#print(vari_norm)

#print(NumZGrid)
#print(-delta_lev_id)

vari_norm_lev_id_ave = [x for x in np.arange(NumZGrid,0,-delta_lev_id)]
#vari_norm_lev_id_ave.reverse()

#print(vari_norm_lev_id_ave)

#print(vari_norm_lev_id_ave)
#vari_norm_lev_id_add = [1] 
#隔几层画隔几层写,可追加指定层
vari_norm_lev_id = vari_norm_lev_id_ave#+vari_norm_lev_id_add
vari_surf_lev_id = [0]
#vari_norm_lev_id = [1,32,40,59]
#vari_surf_lev_id = [60]

#-------------------------------set hash-key name---------------------------

res_in_key = str(NumXGrid)+'x'+str(NumYGrid)
deg= eval(conn.hget(casename,'atm_mpas_sceneid')) #int(180/NumYGrid) #写表名deg改成用小数

#print(deg)

hashkey_lat = 'lat:'+res_in_key #'lat:1deg-201902'
hashkey_lon = 'lon:'+res_in_key


hashkey_tr = 'realfield:'+deg+':'+res_in_key #'realfield:1deg-201902'
hashkey_bg = 'axbfield:'+casename+':'+deg+':'+res_in_key
hashkey_amb = 'ambfield:'+casename+':'+deg+':'+res_in_key

if obs_open == 1:

    obs_hashkey_lat = 'obslat:'+casename+':'+res_in_key
    obs_hashkey_lon = 'obslon:'+casename+':'+res_in_key
    hashkey_obs = 'obs:'+casename+':'+res_in_key

print(hashkey_obs)


##################end_change #########end_change #############end_change ###########end_change################


#------------------------------prepare---------------------------------------


Numlevel = len(vari_norm_lev_id)
Numproc = Numlevel*len(vari_norm)+len(vari_surf)

#print(Numproc)

vari = vari_norm + vari_surf

vari_lev_group_id = []
for x in vari_norm:
   vari_lev_group_id.append(vari_norm_lev_id)
for x in vari_surf:
   vari_lev_group_id.append(vari_surf_lev_id)
vari_lev_id = np.array(vari_lev_group_id)

draw_norm_root = range(0,Numlevel*len(vari_norm),Numlevel)
draw_surf_root = range(Numlevel*len(vari_norm),Numproc,1)
draw_root = list(draw_norm_root)+list(draw_surf_root)
print(draw_root)


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
lat_id=[]
for lat_in in range(1,NumYGrid+1):
        lat_key_str = str(lat_in)+':'+'lat'
        lat_key_get.append(lat_key_str)
        lat_id.append(lat_in)

lon_key_get = []
lon_id=[]
for lon_in in range(1,NumXGrid+1):
        lon_key_str = str(lon_in)+':'+'lon'
        lon_key_get.append(lon_key_str)
        lon_id.append(lon_in)

lat_get = conn.hmget(hashkey_lat, lat_key_get)
lon_get = conn.hmget(hashkey_lon, lon_key_get)

lat_get = [float(x) for x in lat_get]
lon_get = [float(x) for x in lon_get]

lat_theta_get = [[np.cos(x*np.pi/180)] for x in lat_get]

lat = np.array(lat_get)
lon = np.array(lon_get)

lat_theta = np.array(lat_theta_get)

if obs_open == 1:
    
    NumYGrid_obs=int(conn.hlen(obs_hashkey_lat))
    NumXGrid_obs=int(conn.hlen(obs_hashkey_lon))
    Numobs = NumYGrid_obs * NumXGrid_obs
    print(Numobs)

    obs_lat_key_get=[]
    obs_lat_id=[]
    for lat_in in range(1,NumYGrid_obs+1):
            obs_lat_key_str = str(lat_in)+':'+'lat'
            obs_lat_key_get.append(obs_lat_key_str)
            obs_lat_id.append(lat_in)

    obs_lon_key_get=[]
    obs_lon_id=[]
    for lon_in in range(1,NumXGrid_obs+1):
            obs_lon_key_str = str(lon_in)+':'+'lon'
            obs_lon_key_get.append(obs_lon_key_str)
            obs_lon_id.append(lon_in)

    print(obs_lat_key_get)

    obs_lat_get = conn.hmget(obs_hashkey_lat, obs_lat_key_get)
    obs_lon_get = conn.hmget(obs_hashkey_lon, obs_lon_key_get)

    obs_lat_get = [float(x) for x in obs_lat_get]
    obs_lon_get = [float(x) for x in obs_lon_get]


    #obs_lon_id = [1,3,5,7,9,11,13,15]
    #obs_lat_id = [1,3,5,7]

    #obs_lon_get = [0.5, 2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5]
    #obs_lat_get = [-89.5, -87.5, -85.5, -83.5]


    selected_l_id=[]
    selected_r_id=[]
    selected_bot_id=[]
    selected_top_id=[]

    obs_get_latlon_interval=360/NumXGrid

    for y in obs_lon_get:
        m = [lon_id[lon_get.index(x)] for x in lon_get if (x > y-obs_get_latlon_interval and x < y)]
        n = [lon_id[lon_get.index(x)] for x in lon_get if (x > y and x < y+obs_get_latlon_interval)]
        selected_l_id=selected_l_id+m
        selected_r_id=selected_r_id+n

        #print('left',lon_get[m[0]-1])
        #print('right',lon_get[n[0]-1])
        #print('obs-lon',y)


    for y in obs_lat_get:
        m = [lat_id[lat_get.index(x)] for x in lat_get if (x > y-obs_get_latlon_interval and x < y)]
        n = [lat_id[lat_get.index(x)] for x in lat_get if (x > y and x < y+obs_get_latlon_interval)]
        selected_bot_id=selected_bot_id+m
        selected_top_id=selected_top_id+n

        #print('bot',lat_get[m[0]-1])
        #print('top',lat_get[n[0]-1])
        #print('obs-lat',y)

    #print(selected_l_id,selected_r_id,selected_bot_id,selected_top_id)



##########################################取背景场，真实场，分析场，诊断绘图##########################################

va = []
runtime = []
curr_time_from_start = run_start_time
curr_time = curr_time_from_start.strftime("%Y%m%d%H%M")

##--------------------------------------------------obs get---------------------------------------------------------------------

for time_in in range(0,Numtime):
    #runtime.append(curr_time_from_start.strftime("%Y-%m-%d-%H-%M"))
    runtime.append(time_in)


    for proc_in in range(nl_start, nl_end):

        proc_for_norm = Numlevel*(len(vari_norm))
        
        if nl_start < proc_for_norm:
           level_in = nl_start % Numlevel
           vari_in = nl_start // Numlevel
        else:
           level_in = 0
           vari_in = len(vari_norm)-1 + (nl_end - proc_for_norm)


        level_in_key_list = vari_lev_id[vari_in]
        level_in_key = level_in_key_list[level_in]

        print(level_in_key)

        variable_in_key = vari[vari_in]

        key_str = []
        key_get = []
        time_in_key = curr_time

        print(variable_in_key)


        for lat_in in range(1,NumYGrid+1):
            lat_in_key = lat_in
            for lon_in in range(1,NumXGrid+1):
                lon_in_key = lon_in
                key_str = str(lon_in_key)+':'+str(lat_in_key)+':'+str(level_in_key)+':'+str(variable_in_key)
                key_get.append(key_str)


        tr_get_global = conn.hmget(hashkey_tr+':'+str(time_in_key), key_get)
        amb_get_global = conn.hmget(hashkey_amb+':'+str(time_in_key), key_get)
        bg_get_global = conn.hmget(hashkey_bg+':'+str(time_in_key), key_get)

        tr_get_global = [float(x) for x in tr_get_global]
        amb_get_global = [float(x) for x in amb_get_global]
        bg_get_global = [float(x) for x in bg_get_global]

        tr_global_array = np.array(tr_get_global).reshape(NumYGrid,NumXGrid)
        amb_global_array = np.array(amb_get_global).reshape(NumYGrid,NumXGrid)
        bg_global_array = np.array(bg_get_global).reshape(NumYGrid,NumXGrid)
        an_global_array = bg_global_array+amb_global_array

        tr_global_theta = tr_global_array * lat_theta
        bg_global_theta = bg_global_array * lat_theta
        an_global_theta = an_global_array * lat_theta

        tr_global_gather = np.array(comm.allgather(tr_global_array))
        bg_global_gather = np.array(comm.allgather(bg_global_array))
        an_global_gather = np.array(comm.allgather(an_global_array))

        tr_global_theta_gather = np.array(comm.allgather(tr_global_theta))
        bg_global_theta_gather = np.array(comm.allgather(bg_global_theta))
        an_global_theta_gather = np.array(comm.allgather(an_global_theta))

        interval_i = len(level_in_key_list)
        istart = vari_in*Numlevel
        iend = istart+interval_i

        tr_g = tr_global_gather[istart:iend]
        bg_g = bg_global_gather[istart:iend]
        an_g = an_global_gather[istart:iend]

        tr_g_theta = tr_global_theta_gather[istart:iend]
        bg_g_theta = bg_global_theta_gather[istart:iend]
        an_g_theta = an_global_theta_gather[istart:iend]


        grmse = indicator.cal_rmse(an_g, tr_g)
        grmse_theta = indicator.cal_rmse(an_g_theta, tr_g_theta)

        if obs_open == 1:


            if variable_in_key == 'ps':
                level_in_key_obs = 0
            else:
                level_in_key_obs = Numlevel-level_in


            key_get_l_top = []
            key_get_r_top = []
            key_get_l_bot = []
            key_get_r_bot = []
            key_get_obs = []


            for lat_in in range(0,NumYGrid_obs):
                lat_in_key_bot = selected_bot_id[lat_in]
                lat_in_key_top = selected_top_id[lat_in]
                lat_in_key_obs = obs_lat_id[lat_in]

                for lon_in in range(0,NumXGrid_obs):
                    lon_in_key_l = selected_l_id[lon_in]
                    lon_in_key_r = selected_r_id[lon_in]
                    lon_in_key_obs = obs_lon_id[lon_in]

                    key_str_l_top = str(lon_in_key_l)+':'+str(lat_in_key_top)+':'+str(level_in_key)+':'+str(variable_in_key)
                    key_str_r_top = str(lon_in_key_r)+':'+str(lat_in_key_top)+':'+str(level_in_key)+':'+str(variable_in_key)
                    key_str_l_bot = str(lon_in_key_l)+':'+str(lat_in_key_bot)+':'+str(level_in_key)+':'+str(variable_in_key)
                    key_str_r_bot = str(lon_in_key_r)+':'+str(lat_in_key_bot)+':'+str(level_in_key)+':'+str(variable_in_key)
                    key_str_obs = str(lon_in_key_obs)+':'+str(lat_in_key_obs)+':'+str(level_in_key_obs)+':'+str(variable_in_key)

                        #print(key_str_l_top)


                    key_get_l_top.append(key_str_l_top)
                    key_get_r_top.append(key_str_r_top)
                    key_get_l_bot.append(key_str_l_bot)
                    key_get_r_bot.append(key_str_r_bot)
                    key_get_obs.append(key_str_obs)

                    #print(key_get_obs)

            #key_get_4 = [key_get_l_top, key_get_r_top, key_get_l_bot, key_get_r_bot]
            key_get_4 = key_get_l_top + key_get_r_top + key_get_l_bot + key_get_r_bot

            #print(key_get_obs)

            tr_get = conn.hmget(hashkey_tr+':'+str(time_in_key), key_get_4)
            amb_get = conn.hmget(hashkey_amb+':'+str(time_in_key), key_get_4)
            bg_get = conn.hmget(hashkey_bg+':'+str(time_in_key), key_get_4)
            obs_get = conn.hmget(hashkey_obs+':'+str(time_in_key), key_get_obs)

            #print(obs_get)

            tr_get = [float(x) for x in tr_get]
            amb_get = [float(x) for x in amb_get]
            bg_get = [float(x) for x in bg_get]
            obs_get = [float(x) for x in obs_get]

            tr_array = np.array(tr_get)
            amb_array = np.array(amb_get)
            bg_array = np.array(bg_get)
            an_array = bg_array+amb_array
            obs_array = np.array(obs_get)

            tr_obs_1lev = (tr_array[0:Numobs]+tr_array[Numobs:2*Numobs]+tr_array[2*Numobs:3*Numobs]+tr_array[3*Numobs:4*Numobs])/4
            bg_obs_1lev = (bg_array[0:Numobs]+bg_array[Numobs:2*Numobs]+bg_array[2*Numobs:3*Numobs]+bg_array[3*Numobs:4*Numobs])/4
            an_obs_1lev = (an_array[0:Numobs]+an_array[Numobs:2*Numobs]+an_array[2*Numobs:3*Numobs]+an_array[3*Numobs:4*Numobs])/4
            obs_1lev = obs_array

            tr_obs_gather = np.array(comm.allgather(tr_obs_1lev))
            bg_obs_gather = np.array(comm.allgather(bg_obs_1lev))
            an_obs_gather = np.array(comm.allgather(an_obs_1lev))
            obs_gather = np.array(comm.allgather(obs_1lev))

            tr_1d = tr_obs_gather[istart:iend].flatten()
            bg_1d = bg_obs_gather[istart:iend].flatten()
            an_1d = an_obs_gather[istart:iend].flatten()
            ob_1d = obs_gather[istart:iend].flatten()


            osd, asd, bsd, asd_osd, asd_bsd, ormse, brmse, armse = indicator.cal_theta(bg_1d, an_1d, tr_1d, ob_1d)
            #ormse, brmse, armse = indicator.cal_theta(bg_1d, an_1d, tr_1d, ob_1d)
            va.append([osd, asd, bsd, asd_osd, asd_bsd, ormse, brmse, armse, grmse, grmse_theta])
            #va.append([ormse, brmse, armse, grmse, grmse_theta])


        else:


            va.append([grmse, grmse_theta])


        if nl_start in draw_root:

            testplotmulti.ploth(va, runtime, curr_time, variable_in_key, obs_open, pic_path)

    curr_time_from_start = curr_time_from_start + timedelta(seconds = delta_seconds)
    curr_time = curr_time_from_start.strftime("%Y%m%d%H%M")









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
import sys
import netCDF4


t_start=time.time()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
#size = comm.Get_size()


pic_path=sys.argv[1]
#-----------------------------------

cf = configparser.ConfigParser()
cf.read("../namelist/initial_for_obstat.input")

#-------------------------------set redis-------------------------------------

path_nc=eval(cf.get("obstat_ini", "path_nc"))

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

##################start_change #########start_change #############start_change ###########start_change#########
casename = eval(cf.get("obstat_ini", "case_name")) #'case_yz'
#print(casename)

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

#print(Numtime_to_start)

if (start_time-da_start_time).total_seconds() % delta_seconds == 0:
    run_start_time = start_time
else:
    run_start_time = da_start_time
    for i in range(0,Numtime_to_start+1):
        run_start_time = run_start_time + timedelta(seconds = delta_seconds)

Numtime = int(((end_time-run_start_time).total_seconds()) // delta_seconds + 1) 


#print('start time:',run_start_time)
#print('numtime:',Numtime)

#------------------------------set space------------------------------------

NumYGrid = 1260 #int(conn.hget(casename,'num_lat')) #180
NumXGrid = 2560 #int(conn.hget(casename,'num_lon')) #360
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
#print(vari_norm_lev_id_ave)
#vari_norm_lev_id_add = [1] 
#隔几层画隔几层写,可追加指定层
vari_norm_lev_id = vari_norm_lev_id_ave#+vari_norm_lev_id_add
vari_surf_lev_id = [0]
#vari_norm_lev_id = [1,32,40,59]
#vari_surf_lev_id = [60]

Numblock=eval(cf.get("obstat_ini", "numblock"))
NumYblock=NumYGrid//Numblock
NumXblock=NumXGrid//Numblock

#-------------------------------set hash-key name---------------------------

res_in_key = str(NumXGrid)+'x'+str(NumYGrid)
#res_in_key = str(360)+'x'+str(180)
deg = eval(conn.hget(casename,'atm_mpas_sceneid')) 
test_case = eval(conn.hget(casename,'test_case'))

#print(deg)
cc_case = deg.split('-',2)[2]

#print('cc_case',cc_case)

#print(deg)

hashkey_lat = 'lat:'+res_in_key #'lat:1deg-201902'
hashkey_lon = 'lon:'+res_in_key

#print(hashkey_lat)

#hashkey_tr = 'realfield:'+deg+':'+res_in_key #'realfield:1deg-201902'
#hashkey_bg = 'axbfield:'+casename+':'+deg+':'+res_in_key
#hashkey_amb = 'ambfield:'+casename+':'+deg+':'+res_in_key

hashkey_tr = 'yuz:2580x1260' #'realfield:1deg-201902'
#hashkey_bg = 'yuz:'+res_in_key
#hashkey_amb = 'yuz:'+res_in_key
bg_name_nc = 'test_15km'
amb_name_nc = 'test_15km'

#print(hashkey_tr)
#print(hashkey_bg)


##################end_change #########end_change #############end_change ###########end_change################
#------------------------------prepare---------------------------------------


Numlevel = len(vari_norm_lev_id)
#Numproc = Numlevel*(len(vari_norm))+len(vari_surf)
Numproc = (Numlevel*(len(vari_norm))+len(vari_surf))*Numblock
rank_gather = 0

#print(Numproc)

vari = vari_norm + vari_surf

vari_lev_group_id = []
for x in vari_norm:
   vari_lev_group_id.append(vari_norm_lev_id)
for x in vari_surf:
   vari_lev_group_id.append(vari_surf_lev_id)
vari_lev_id = np.array(vari_lev_group_id)


#----------------------------------sspadd-------------------------------------
#nl_sec = Numproc // size
#nl_res = Numproc % size
nl_sec = 1
nl_res = 0
if(rank < nl_res) :
    nl_start = rank * (nl_sec + 1)
    nl_end = nl_start + nl_sec + 1
else:
    nl_start = rank * nl_sec + nl_res
    nl_end = nl_start + nl_sec

print("myrank:", rank, nl_start, nl_end)
########################################取经纬度信息################################



lat_key_get = []
for lat_in in range(1,NumYGrid+1):
        lat_key_str = str(lat_in)+':'+'lat'
        lat_key_get.append(lat_key_str)

lon_key_get = []
for lon_in in range(1,NumXGrid+1):
        lon_key_str = str(lon_in)+':'+'lon'
        lon_key_get.append(lon_key_str)

#lat_get = conn.hmget(hashkey_lat, lat_key_get)
#lon_get = conn.hmget(hashkey_lon, lon_key_get)

lat_get = np.arange(-90,90,180/1260)
lon_get = np.arange(0,360,360/2560)

lat_get = [float(x) for x in lat_get]
lon_get = [float(x) for x in lon_get]


lat_theta_get = [[np.cos(x*np.pi/180)] for x in lat_get]

lat = np.array(lat_get)
lon = np.array(lon_get)

lat_theta = np.array(lat_theta_get)
#print(lat_theta)

id_interval = len(lat)//6
id_lat = np.arange(0,len(lat)+1,id_interval)


#########################################取背景场，真实场，分析场，诊断绘图##########################################

va = []
va_theta = []
runtime = []
curr_time_from_start = run_start_time
curr_time = curr_time_from_start.strftime("%Y%m%d%H%M")

#for time_in in range(0,Numtime):
for time_in in range(0,1):
    key_str = []
    key_get = []
    runtime.append(curr_time)


    for proc_in in range(nl_start, nl_end):

        #proc_for_norm = Numlevel*(len(vari_norm))
        proc_for_norm = Numlevel*(len(vari_norm))*Numblock
        
        if nl_start < proc_for_norm:
           level_in = nl_start // Numblock % Numlevel
           vari_in = nl_start // Numblock // Numlevel
        else:
           level_in = 0
           #vari_in = len(vari_norm)-1 + (nl_end - proc_for_norm)
           vari_in = len(vari_norm)-1 + ((nl_start//Numblock+1) - int(proc_for_norm//Numblock))


        level_in_key_list = vari_lev_id[vari_in]
        level_in_key = level_in_key_list[level_in]

        #print(level_in_key)

        variable_in_key = vari[vari_in]

        time_in_key = curr_time

        t0=time.time()
        #print(variable_in_key)


        #for lat_in in range(1,NumYGrid+1):
        #lat_block_start=nl_start%Numblock*NumYblock+1
        #lat_block_end=(nl_start%Numblock+1)*NumYblock+1
        lon_block_start=nl_start%Numblock*NumXblock+1
        lon_block_end=(nl_start%Numblock+1)*NumXblock+1

        pipe1=conn.pipeline()
        #pipe2=conn.pipeline()
        #pipe3=conn.pipeline()

        #for lat_in in range(lat_block_start,lat_block_end):
        for lon_in in range(lon_block_start,lon_block_end): 

        #for lat_in in range(1,5+1):
            #for lon_in in range(1,NumXGrid+1): 
            key_get = []
            for lat_in in range(1,NumYGrid+1):
            #for lon_in in range(1,5+1): 

                lat_in_key = lat_in
                lon_in_key = lon_in 
                key_str = str(lat_in_key)+':'+str(level_in_key)+':'+str(variable_in_key)
                key_get.append(key_str)

            pipe1.hmget(hashkey_tr+':'+str(time_in_key)+":"+str(lon_in_key), key_get)
            #pipe2.hmget(hashkey_bg+':'+str(time_in_key)+":"+str(lon_in_key), key_get)
            #pipe3.hmget(hashkey_amb+':'+str(time_in_key)+":"+str(lon_in_key), key_get)
                #print(key_get[2])
                #print(hashkey_tr+':'+str(time_in_key))

        t1=time.time()
        print("拼接：",t1-t0)


        #print(key_get)
        #print(hashkey_tr+':'+str(time_in_key))
        tr_get = np.array(pipe1.execute()).flatten()
        #bg_get = np.array(pipe2.execute()).flatten()
        #amb_get = np.array(pipe3.execute()).flatten()

        t2=time.time()
        print("hmget：",t2-t1)
        #print("len_get:", len(tr_get[0]))

        tnc0=time.time()

        if nl_start < proc_for_norm:
            #bg_get=netCDF4.Dataset(path_nc+bg_name_nc+':'+str(time_in_key)+".nc").variables[variable_in_key][(level_in_key-1),(lat_block_start-1):(lat_block_end-1),:]
            #amb_get=netCDF4.Dataset(path_nc+amb_name_nc+':'+str(time_in_key)+".nc").variables[variable_in_key][(level_in_key-1),(lat_block_start-1):(lat_block_end-1),:]
            bg_get=netCDF4.Dataset(path_nc+bg_name_nc+".nc").variables[variable_in_key][(level_in_key-1),0:NumYGrid,(lon_block_start-1):(lon_block_end-1)]
            amb_get=netCDF4.Dataset(path_nc+amb_name_nc+".nc").variables[variable_in_key][(level_in_key-1),0:NumYGrid,(lon_block_start-1):(lon_block_end-1)]
        else:
            #bg_get=netCDF4.Dataset(path_nc+bg_name_nc+':'+str(time_in_key)+".nc").variables[variable_in_key][(lat_block_start-1):(lat_block_end-1),:]
            #amb_get=netCDF4.Dataset(path_nc+amb_name_nc+':'+str(time_in_key)+".nc").variables[variable_in_key][(lat_block_start-1):(lat_block_end-1),:]
            bg_get=netCDF4.Dataset(path_nc+bg_name_nc+".nc").variables['s'][0:NumYGrid,(lon_block_start-1):(lon_block_end-1)]
            amb_get=netCDF4.Dataset(path_nc+amb_name_nc+".nc").variables['s'][0:NumYGrid,(lon_block_start-1):(lon_block_end-1)]

        tnc1=time.time()
        print("nc_get:",tnc1-tnc0)



        #print(amb_get)

        tr_get = [float(x) for x in tr_get]
        #bg_get = [float(x) for x in bg_get]
        #amb_get = [float(x) for x in amb_get]
        bg_get = bg_get.T
        amb_get = amb_get.T

        tr_array = np.array(tr_get)
        bg_array = np.array(bg_get)
        amb_array = np.array(amb_get)

        t3=time.time()


        t_comm_block0=time.time()

        color_block = rank / Numblock
        block_comm = comm.Split(color_block, rank)
        block_comm_rank = block_comm.Get_rank()

        t_comm_block1=time.time()
        print("comm_block:", t_comm_block1 - t_comm_block0)

        #tr=np.array(block_comm.gather(tr_array, root=0))
        t222=time.time()
        tr_block_gather = block_comm.gather(tr_array, root=0)
        bg_block_gather = block_comm.gather(bg_array, root=0)
        amb_block_gather = block_comm.gather(amb_array, root=0)
        t333=time.time()
        print("gather",t333-t222)

        if block_comm_rank == 0:
           tr = ((np.array(tr_block_gather)).flatten()).reshape(NumXGrid,NumYGrid).T
           bg = ((np.array(bg_block_gather)).flatten()).reshape(NumXGrid,NumYGrid).T
           amb = ((np.array(amb_block_gather)).flatten()).reshape(NumXGrid,NumYGrid).T
           an = bg + amb
           print(len(tr))
           print(len(tr[0]))
           tr_theta = tr * lat_theta
           bg_theta = bg * lat_theta
           amb_theta = amb * lat_theta
           an_theta = bg_theta + amb_theta

           t_lev0=time.time()

           rmse_theta_6090S = indicator.cal_rmse(an_theta[id_lat[0]:id_lat[1],:], tr_theta[id_lat[0]:id_lat[1],:])
           rmse_theta_3060S = indicator.cal_rmse(an_theta[id_lat[1]:id_lat[2],:], tr_theta[id_lat[1]:id_lat[2],:])
           rmse_theta_0030S = indicator.cal_rmse(an_theta[id_lat[2]:id_lat[3],:], tr_theta[id_lat[2]:id_lat[3],:])
           rmse_theta_0030N = indicator.cal_rmse(an_theta[id_lat[3]:id_lat[4],:], tr_theta[id_lat[3]:id_lat[4],:])
           rmse_theta_3060N = indicator.cal_rmse(an_theta[id_lat[4]:id_lat[5],:], tr_theta[id_lat[4]:id_lat[5],:])
           rmse_theta_6090N = indicator.cal_rmse(an_theta[id_lat[5]:id_lat[6],:], tr_theta[id_lat[5]:id_lat[6],:])
           rmse_theta_level = indicator.cal_rmse(an_theta[:,:], tr_theta[:,:])

           t_lev1=time.time()
           print("rmse_lev:",t_lev1-t_lev0)

        t_comm_group0=time.time()

        group_draw_tag = np.arange(0,Numproc,Numblock)
        group_draw = comm.Get_group()
        group_draw = group_draw.Incl(group_draw_tag)
        group_draw_comm = comm.Create(group_draw)

        t_comm_group1=time.time()
        print("comm_group:",t_comm_group1 - t_comm_group0)

        if group_draw_comm != MPI.COMM_NULL:

            t_comm_draw0=time.time()

            group_draw_comm_rank = group_draw_comm.Get_rank()
            color_draw = group_draw_comm_rank / Numlevel
            draw_comm = group_draw_comm.Split(color_draw, group_draw_comm_rank)
            draw_comm_rank = draw_comm.Get_rank()

            t22=time.time()
            print("comm_draw:",t22 - t_comm_draw0)

            tr_theta_allgather = draw_comm.allgather(tr_theta)
            bg_theta_allgather = draw_comm.allgather(bg_theta)
            an_theta_allgather = draw_comm.allgather(an_theta)
            t33=time.time()
            print("allgather",t33-t22)
            rmse_global_theta = indicator.cal_rmse(np.array(an_theta_allgather), np.array(tr_theta_allgather))
            t_glo_1=time.time()
            print("rmse_global",t_glo_1-t33)


            va_theta.append([rmse_theta_6090N, rmse_theta_3060N, rmse_theta_0030N, rmse_theta_0030S, rmse_theta_3060S, rmse_theta_6090S, rmse_theta_level, rmse_global_theta])

            print(len(tr_theta_allgather))
            print(len(tr))


            pic_param = [runtime, curr_time, level_in_key, level_in, variable_in_key, pic_path, cc_case]
            testplotmulti.plotlev(bg, an, tr, va_theta, pic_param)

            draw_comm.Free()
            group_draw_comm.Free()
            block_comm.Free()




        t4=time.time()
        print("gather_draw",t4-t3)




        #va_print_gather = np.array(comm.gather(rmse_global_theta, root=0))


        #if rank == 0:
            #va_print_id = [int(x*Numlevel) for x in range(0,len(vari_norm))]+[Numproc-1]
            #print(va_print_id)
            #va_print = str([round(va_print_gather[x],3) for x in va_print_id])
            #print(va_print)
            #vari_f = str(vari)

            #if time_in == 0:
                #f = open('./vv.dat', 'w')
                #f.writelines(vari_f)
                #f.write('\n')

                #f.writelines(va_print)
                #f.write('\n')
                #f.close()
            #else:
                #f = open('./vv.dat', 'a')
                #f.writelines(va_print)
                #f.write('\n')
                #f.close()



    curr_time_from_start = curr_time_from_start + timedelta(seconds = delta_seconds)
    curr_time = curr_time_from_start.strftime("%Y%m%d%H%M")

t_end=time.time()
print("end-start:",t_end - t_start)



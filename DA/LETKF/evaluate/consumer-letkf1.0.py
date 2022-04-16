from rediscluster import RedisCluster
#import redis
import json
import numpy as np
from datetime import datetime
#from colorrange import colorrange
import cal_indicator as indicator
import testplot
import os
import time
#sspadd------------------------------
from mpi4py import MPI 

#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#size = comm.Get_size()
#-----------------------------------
t1 = time.time()
Numtime  = 20
Numlevel = 1
NumYGrid = 180
NumXGrid = 360
va = []
lat_key_get = []
lon_key_get = []
key_get = []
key_tr_get = []
runtime = []

startup_nodes = [                           #集群节点及端口
    {"host":"192.168.1.1", "port":6379},
    {"host":"192.168.1.2", "port":6379},
    {"host":"192.168.1.3", "port":6379}
]

#conn = redis.StrictRedis(host='localhost', port=6379, decode_responses=True)
conn = RedisCluster(startup_nodes=startup_nodes, decode_responses=True)
t2 = time.time()
print("connect:", t2-t1)

##--------------sspadd-----------------
#nl_sec = Numlevel // size
#nl_res = Numlevel % size
#if(rank < nl_res) :
    #nl_start = rank * (nl_sec + 1)
    #nl_end = nl_start + nl_sec + 1
#else:
    #nl_start = rank * nl_sec + nl_res
    #nl_end = nl_start + nl_sec
###------------------------------------
#print("myrank:", rank, nl_start, nl_end)
#########################################取经纬度信息################################

#for lat_in in range(0,NumYGrid):
#        lat_key_str = str(lat_in)
#        lat_key_get.append(lat_key_str)
#
#for lon_in in range(0,NumXGrid):
#        lon_key_str = str(lon_in)
#        lon_key_get.append(lon_key_str)
#
#lat_get = conn.hmget('lat_key', lat_key_get)
#lon_get = conn.hmget('lon_key', lon_key_get)
#
#lat_get = [float(x) for x in lat_get]
#lon_get = [float(x) for x in lon_get]
#
#lat = np.array(lat_get)
#lon = np.array(lon_get)
#
#id_interval = int(len(lat)/6)
#id_lat = np.arange(0,len(lat),id_interval)
##print(id_lat)
#
t3 = time.time()
print("hmget: ", t3-t2)
#
##print(lat)
##print(lon)

#id_interval=30
#id_lat = np.arange(NumYGrid,1,-id_interval)
lat = range(1,180)
lon = range(1,360)
####################################################################################

for level_in in range(0,Numlevel):
    for lat_in in range(1,NumYGrid+1):
        for lon_in in range(1,NumXGrid+1): 
            key_str = str(lat_in)+'.'+str(lon_in)
            key_tr_get.append(key_str)

h_tr_get = conn.hmget('tr-letkf1.0', key_tr_get)
h_tr_get = [float(x) for x in h_tr_get]
h_tr_array = np.array(h_tr_get)
h_tr = h_tr_array.reshape(NumYGrid,NumXGrid)


for time_in in range(1,Numtime+1):
    key_str = []
    key_get = []
    runtime.append(time_in)
    curr_time=str(time_in)

    for level_in in range(0,Numlevel):
    #for level_in in range(nl_start, nl_end):
        for lat_in in range(1,NumYGrid+1):
            for lon_in in range(1,NumXGrid+1): 
                key_str = str(time_in)+'.'+str(lat_in)+'.'+str(lon_in)
                key_get.append(key_str)
#               print(get)

        t4 = time.time()
        print("拼接:", t4-t3)

    #    print(key_get[0])

        #h_tr_get = conn.hmget('tr-key', key_get)
        h_bg_get = conn.hmget('bg-letkf1.0', key_get)        
        h_an_get = conn.hmget('an-letkf1.0', key_get)

        t5 = time.time()
        print("hmget", t5-t4)
    #    print(h_tr_get[0])

        #h_tr_get = [float(x) for x in h_tr_get]
        h_bg_get = [float(x) for x in h_bg_get]
        h_an_get = [float(x) for x in h_an_get]

        #h_tr_array = np.array(h_tr_get)
        h_bg_array = np.array(h_bg_get)
        h_an_array = np.array(h_an_get)

        
        #h_tr = h_tr_array.reshape(NumYGrid,NumXGrid)
        h_bg = h_bg_array.reshape(NumYGrid,NumXGrid)
        h_an = h_an_array.reshape(NumYGrid,NumXGrid)

#        print(h_tr[0,2])

        t9 = time.time()

        h_rmse_6090N = indicator.cal_rmse(h_an[150:180,:], h_tr[150:180,:])
        h_rmse_3060N = indicator.cal_rmse(h_an[120:150,:], h_tr[120:150,:])
        h_rmse_0030N = indicator.cal_rmse(h_an[90:120,:], h_tr[90:120,:])
        h_rmse_0030S = indicator.cal_rmse(h_an[60:90,:], h_tr[60:90,:])
        h_rmse_3060S = indicator.cal_rmse(h_an[30:60,:], h_tr[30:60,:])
        h_rmse_6090S = indicator.cal_rmse(h_an[0:30,:], h_tr[0:30,:])
        t10 = time.time() 
        print("h_rmse:", t10-t9)



#        time.append(time_in)
        t90 = time.time() 
        
        f = open('./vv.dat', 'a')
        f.writelines('{} {} {} {} {} {} {} {}'.format(time_in, '%.2f' % level_in, '%.2f' % h_rmse_6090N, '%.2f' % h_rmse_3060N, '%.2f' % h_rmse_0030N, '%.2f' % h_rmse_0030S, '%.2f' % h_rmse_3060S, '%.2f' % h_rmse_6090S))
        f.write('\n')
        f.close()
        
        t100 = time.time() 
        print("h_rmse2:", t100-t90)
        
        va.append([h_rmse_6090N, h_rmse_3060N, h_rmse_0030N, h_rmse_0030S, h_rmse_3060S, h_rmse_6090S])

        colorrange = np.array([[0, 200], [-100, 100], [-65, 65], [-110, 110], [-0.55, 0.55], [-0.55, 0.55]])

#        print(h_bg)
        
        va_test=np.array(va)
#        print(va_test[:,0])

        t101 = time.time() 
        testplot.ploth(lon, lat, h_bg, h_an, h_tr, va, runtime, curr_time, colorrange, level_in, Numlevel) 
        t102 = time.time() 
        print("testplot:", t102-t101)

    t6 = time.time()
    print("ploth: ", t6-t5)
#    print(va[0::5,0])



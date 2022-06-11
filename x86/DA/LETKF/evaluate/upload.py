from rediscluster import RedisCluster
import numpy as np
import os


startup_nodes = [                           #集群节点及端口
     {"host":"192.168.1.1", "port":6379},
     {"host":"192.168.1.2", "port":6379},
     {"host":"192.168.1.3", "port":6379}
]

conn = RedisCluster(startup_nodes=startup_nodes, decode_responses=True)

path_tr = os.getcwd() + "/../result/actual/"
path_bg = os.getcwd() + "/../result/back/"
path_an = os.getcwd() + "/../result/ana/"

file_tr=os.listdir(path_tr)
file_bg=os.listdir(path_bg)
file_an=os.listdir(path_an)

print(file_tr)

#lat_tr=[]
#lon_tr=[]
#tr=[]


#lat_bg=[]
#lon_bg=[]
#bg=[]

#lat_an=[]
#lon_an=[]
#an=[]

dict_tr={}
dict_bg={}
dict_an={}

pipe=conn.pipeline()

for file in file_tr:
   position_tr = path_tr+file
   f_tr = open(position_tr)
   data_tr = f_tr.readlines()

   for string in data_tr:
       string = string.split()
       lat_get=int(float(string[0]))
       lon_get=int(float(string[1]))
       tr_get=float(string[2])
       key_tr=str(lat_get)+'.'+str(lon_get)
       append_dict_tr={key_tr:tr_get}
       pipe.hmset('tr-letkf1.0',append_dict_tr)
       #dict_tr.update(append_dict_tr)

       
       #lat_tr.append(lat_get)
       #lon_tr.append(lon_get)
       #tr.append(tr_get)



for file in file_bg:
  position_bg = path_bg+file
  f_bg = open(position_bg)
  data_bg = f_bg.readlines()

  for string in data_bg:
      string = string.split()
      nt_get=int(float(string[0]))
      lat_get=int(float(string[1]))
      lon_get=int(float(string[2]))
      bg_get=float(string[3])
      key_bg=str(nt_get)+'.'+str(lat_get)+'.'+str(lon_get)
      append_dict_bg={key_bg:bg_get}
      pipe.hmset('bg-letkf1.0',append_dict_bg)
      dict_bg.update(append_dict_bg)


      #lat_bg.append(lat_get)
      #lon_bg.append(lon_get)
      #bg.append(bg_get)




for file in file_an:
  position_an = path_an+file
  f_an = open(position_an)
  data_an = f_an.readlines()

  for string in data_an:
      string = string.split()
      nt_get=int(float(string[0]))
      lat_get=int(float(string[1]))
      lon_get=int(float(string[2]))
      an_get=float(string[3])
      key_an=str(nt_get)+'.'+str(lat_get)+'.'+str(lon_get)
      append_dict_an={key_an:an_get}
      pipe.hmset('an-letkf1.0',append_dict_an)
      #dict_an.update(append_dict_an)


print('start-redis')


      #lat_an.append(lat_get)
      #lon_an.append(lon_get)
      #an.append(an_get)

  #data_tr = np.loadtxt(position_tr)
  #print(data_tr)
  #tr.append(data_tr)

#f_test=np.loadtxt("./actual/actual.out1")



#ff=np.array(f_test)
#print(lat)

#print(dict_bg)

#conn.hmset('tr-letkf1.0',dict_tr)
#conn.hmset('bg-letkf1.0',dict_bg)
#conn.hmset('an-letkf1.0',dict_an)
pipe.execute()









import mpi4py.MPI as MPI
import numpy as np
import sys
import time
import PostProcessor
import TrainData 
import Model
import config as cf
import sys
import scipy as sp
import scipy.ndimage
# saved_model_path = '/home/xuewei/ddn/share/sin_cnnModel/'
# saved_model_path = '/home/xuewei/ddn/share/16_cnnModel/'
saved_model_path = '/home/xuewei/ddn/share/test_cnn/'
# saved_model_path = '/home/xuewei/ddn/share/test_cnn_50km/'
# data_path = '/home/xuewei/ddn/share/svr_data/svr_data_100km_12month_ACGrid/'
# saved_model_path = data_path +'model/'

cf.init()
inter_num  = cf.get_inter_num()
halo = 8
x_num = 2
y_num = 2
comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()
comm_size = comm.Get_size()
bias_rank = 91
group_num = 1
layer_num = cf.get_layer_num()
model_num = cf.get_model_num()
all_size = cf.get_all_size()

group_rank = (comm_rank - bias_rank) % group_num
group_root = group_rank + bias_rank
group_bloc = int((comm_rank - bias_rank) / group_num)
data_sz = np.empty(10,dtype=np.int32)*0
newGroup = comm.group.Excl([x for x in range(bias_rank)])
new_comm = comm.Create_group(newGroup)
new_comm = new_comm.Split(group_rank)
new_comm_rank = new_comm.Get_rank()
new_comm_size = new_comm.Get_size()
train_sz = int(data_sz[1] * data_sz[2] * data_sz[3]*data_sz[8] / new_comm_size)
train_data = np.empty(train_sz)

proxy = PostProcessor.PostProcessor()


def run(iternum = 0):
	#************
	# get config color proc
	#************
	if comm_rank == group_root:
		# !!!!!!!!!!1
		comm.Recv([data_sz, 10, MPI.INT], comm_rank - group_num, 0)
		# DNN: recv info [ 64  97 180 360   1 180   1 360   2   8]
		print('DNN: recv info',data_sz,flush = True)
	new_comm.Bcast([data_sz, 10, MPI.INT],root = int(group_rank/group_num))
	new_comm.barrier()
	inter_num = data_sz[9]
	total_num = int(sys.argv[1])

	base_num = data_sz[8]
	gen_num =int((total_num - base_num)/base_num)
	cf.set_div_range(gen_num)
	div_range = cf.get_div_range()

	# return
			
	# 	#************
	# 	#init model
	# 	#************
	if iternum == 0:
		svr = Model.MulityCNN(div_range = cf.get_div_range(), comm = new_comm)
		svr.data_sz = data_sz
		svr.halo = halo
		proxy.model = svr
		proxy.load(saved_model_path)

	if comm_rank == group_root:
		recv_sz = data_sz[1] * data_sz[2] * data_sz[3]*data_sz[8]
		recv_data = np.empty(recv_sz)

		# !!!!!!!!!!2
		comm.Recv([recv_data, recv_sz, MPI.DOUBLE], comm_rank - group_num, 0)
		print('DNN: recv info',recv_data.shape,flush = True)

		# 数据处理
		# [ 32  107  90 180   1  90   1 180   2]
		recv_data = recv_data.reshape(data_sz[3],data_sz[2],data_sz[1],data_sz[8])
		recv_data = recv_data.transpose(3,2,1,0)
		print('recv_data:',recv_data.shape,flush = True)

		tmp_recv = recv_data
		data_sz_lr = recv_data.shape[0]*halo*recv_data.shape[1]*recv_data.shape[2]

		item_recv_right =  recv_data[:,:,:,: halo]
		item_recv_left = recv_data[:,:,:,-halo:]

		recv_data = np.append(item_recv_left,recv_data,axis =3)
		recv_data = np.append(recv_data,item_recv_right,axis =3)
		print('recv_data:',recv_data.shape,flush = True)

		item_recv_up = np.repeat(recv_data[:,:,[0]], repeats=halo, axis=2)
		item_recv_down = np.repeat(recv_data[:,:,[-1]], repeats=halo, axis=2)

		print('item_recv_up:',item_recv_up.shape,flush = True)
		print('item_recv_down:',item_recv_down.shape,flush = True)


		recv_data = np.append(item_recv_up,recv_data,axis =2)
		recv_data = np.append(recv_data,item_recv_down,axis =2)

		# item_recv_up = np.repeat(recv_data[:,:,[0]], repeats=halo, axis=2)
		# item_recv_down = np.repeat(recv_data[:,:,[-1]], repeats=halo, axis=2)

		# recv_data = np.append(item_recv_up,recv_data,axis =2)
		# recv_data = np.append(recv_data,item_recv_down,axis =2)

		print('recv_data:',recv_data.shape,flush = True)
	else:
		recv_data = train_data = np.empty([data_sz[8],data_sz[1],data_sz[2]+halo*2,data_sz[3]+halo*2])
	# new_comm.barrier()
	new_comm.Bcast(recv_data,0)
	proxy.data = recv_data

	if comm_rank == group_root:
		recv_sz =  data_sz[1] * data_sz[2] * data_sz[3] * data_sz[8]
		result_data = np.empty((data_sz[8] * data_sz[2] * data_sz[3] * div_range*new_comm.Get_size()))
	else:
		result_data = None

	div_begin = cf.get_div_begin()
	new_comm.barrier()

	result = proxy.predict()
	result = result.reshape(div_range * data_sz[8] *  data_sz[2] * data_sz[3] )
	new_comm.barrier()
	new_comm.Gather(result, result_data, root=0)
	new_comm.barrier()



	if comm_rank == group_root:
		result_data = result_data.reshape(new_comm.Get_size(),div_range*data_sz[8],data_sz[2],data_sz[3])
		# 180 * 90 * 97 * 2
		result_data = result_data.transpose(1,0,2,3)
		result_data = np.append(tmp_recv,result_data,axis=0)
		if iternum == 0:
			np.save(str(iternum)+'_'+str(comm_rank)+"result.npy",result_data)
			
		result_data = result_data.transpose(3,2,1,0)
		total_size = (div_range+1) * recv_sz
		send_data = tmp_recv.reshape(recv_sz,1)
		result_data = result_data.reshape(total_size)
		comm.Send([result_data, total_size, MPI.DOUBLE], comm_rank - group_num, 0)
		
	del result_data
	del recv_data
	return inter_num


# 程序执行多少次
inter_num = run(0) 
for x in range(1,inter_num):
	run(x)





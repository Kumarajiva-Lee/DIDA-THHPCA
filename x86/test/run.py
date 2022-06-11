import mpi4py.MPI as MPI
import numpy as np
import sys
import PostProcessor
from vae import Model
import config as cf
import sys

saved_model_path = '/ddnstor/xuewei/zzh/VAE/4/models/'

cf.init()
inter_num  = 40
comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()
bias_rank = 91
group_num = 1

group_rank = (comm_rank - bias_rank) % group_num
group_root = group_rank + bias_rank
data_sz = np.empty(10,dtype=np.int32)*0
newGroup = comm.group.Excl([x for x in range(bias_rank)])
new_comm = comm.Create_group(newGroup)
new_comm = new_comm.Split(group_rank)

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
		svr = Model(comm = new_comm)#todo
		proxy.model = svr
		proxy.load(saved_model_path)#todo

	if comm_rank == group_root:
		recv_sz = data_sz[1] * data_sz[2] * data_sz[3]*data_sz[8]
		recv_data = np.empty(recv_sz)

		# !!!!!!!!!!2
		comm.Recv([recv_data, recv_sz, MPI.DOUBLE], comm_rank - group_num, 0)
		print('DNN: recv data',recv_data.shape,flush = True)

		# [ 32  107  90 180   1  90   1 180   2]
		recv_data = recv_data.reshape(data_sz[3],data_sz[2],data_sz[1],data_sz[8])
		recv_data = recv_data.transpose(3,2,1,0)
		print('recv_data:',recv_data.shape,flush = True)
		if iternum == 0:
			np.save(str(iternum)+'_'+str(comm_rank)+"original.npy",recv_data)

		tmp_recv = recv_data
	else:
		recv_data = np.empty([data_sz[8],data_sz[1],data_sz[2],data_sz[3]])
	new_comm.barrier()
	new_comm.Bcast(recv_data,root=0)
	proxy.data = recv_data

	if comm_rank == group_root:
		recv_sz =  data_sz[1] * data_sz[2] * data_sz[3] * data_sz[8]
		result_data = np.empty((data_sz[8] * data_sz[2] * data_sz[3] * div_range*new_comm.Get_size()))
	else:
		result_data = None

	new_comm.barrier()
	result = proxy.predict()#todo
	result = result.reshape(div_range * data_sz[8] *  data_sz[2] * data_sz[3] )
	new_comm.barrier()
	new_comm.Gather(result, result_data, root=0)
	new_comm.barrier()



	if comm_rank == group_root:
		result_data = result_data.reshape(new_comm.Get_size(),div_range*data_sz[8],data_sz[2],data_sz[3])
		# 180 * 90 * 97 * 2
		result_data = result_data.transpose(1,0,2,3)
		result_data = np.append(tmp_recv,result_data,axis=0)
		print(result_data.shape,flush = True)
		if iternum == 0:
			np.save(str(iternum)+'_'+str(comm_rank)+"result.npy",result_data)
			
		result_data = result_data.transpose(3,2,1,0)
		total_size = (div_range+1) * recv_sz
		result_data = result_data.reshape(total_size)
		comm.Send([result_data, total_size, MPI.DOUBLE], comm_rank - group_num, 0)
		
	del result_data
	del recv_data
	return inter_num


inter_num = run(0) 
for x in range(1,inter_num):
	run(x)

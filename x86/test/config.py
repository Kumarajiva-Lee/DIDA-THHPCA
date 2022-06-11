def init():
	global bias_rank
	bias_rank = 64
	global group_num
	group_num = 4
	global layer_num
	layer_num = 97
	global basic_num
	basic_num = 4
	global model_num
	model_num = 1
	global all_size 
	# all_size = 1940
	all_size = 1552 
	global lev_l
	lev_l ={}
	lev_l['phs'] = 1
	lev_l['pt'] = 32
	lev_l['u_Agrid'] = 32
	lev_l['v_Agrid'] = 32

	global tar_l
	tar_l = ['phs','pt','u_Agrid','v_Agrid']
	global ord_
	ord_ = 0

	global div_range
	# svr 生成的member数量
	div_range = 31

	global div_begin
	# svr 生成的member数量
	div_begin = 2

	global inter_num
	# python 与 fortran 通信的次数
	inter_num = 40

def set_div_range(num):
	global div_range
	div_range = num

def get_bias_rank():
	return bias_rank

def get_bias_rank():
	return bias_rank

def get_group_num():
	return group_num

def get_layer_num():
	return layer_num

def get_basic_num():
	return basic_num

def get_model_num():
	return model_num

def get_all_size ():
	return all_size 

def get_tar_l():
	return tar_l

def get_ord_():
	return ord_ 

def get_div_range():
	return div_range 

def get_div_begin():
	return div_begin 

def get_inter_num():
	return inter_num  

import os

chrd = '/home/xuewei/ddn/liminyan/sere_test/dida_cnn_merge/test_50km/test_m'
exe = '/home/xuewei/ddn/liminyan/sere_test/dida_cnn_merge/alldida/dida_'
os.chdir('/home/xuewei/ddn/liminyan/sere_test/dida_cnn_merge/test')


print(os.getcwd())
keyl = ['4']

# ,'1','2','4']
member = ['128']

# ,'32','64','128']
os.system('pwd')

for tar,fix_member in zip(keyl,member):
	# fix_member = str(int(tar)*32)
	print(tar,fix_member)
	# fix_member = '32'
	work_path = chrd + tar + '_cnn'+fix_member
	namelist = 'namelists_'+tar+'.input'
	tar_exe = exe +tar+'.exe'
	# python_runable = ' ../python_test.py '
	python_runable = ' /home/xuewei/ddn/liminyan/sere_test/dida_cnn_merge/test/python_test50km.py '

	os.chdir(work_path)
	print(os.getcwd())
	os.system('ls')

	cmd = 'salloc --exclusive -N 113 -n 452  -p cnCPU  --ntasks-per-node=4  mpirun -n 64  ' \
	      + tar_exe + ' ' + namelist + ' : -n 388  ~/.conda/envs/py36/bin/python3.6 '+ python_runable + fix_member + ' | tee -a log.out'
	print(cmd)
	os.system(cmd)
	# os.chdir('../../test/')









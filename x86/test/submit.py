import os

chrd = 'test_m'
exe = '../../alldida/'
print(os.getcwd())
keyl = ['4']
member = [128]
# keyl = ['2']
# member = [64]
keyl = ['1']
member = [32]


for tar,fix_member in zip(keyl,member):
	work_path = chrd + tar + '_cnn'+str(fix_member)
	namelist = 'namelists_'+tar+'.input'
	tar_exe = exe + 'dida_'+tar+'.exe'
	python_runable = ' ../python_test.py '

	os.chdir(work_path)
	print(os.getcwd())

	cmd = 'salloc --exclusive -N 113 -n 452  -p cnCPU  --ntasks-per-node=4  mpirun -n 64  ' \
	      + tar_exe + ' ' + namelist + ' : -n 388   ~/.conda/envs/py36/bin/python3.6 '+ python_runable + str(fix_member) + ' | tee -a log.out'
	print(cmd)
	#os.system(cmd)
	os.chdir('../')









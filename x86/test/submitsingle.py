import os

chrd = 'single_point_m'
exe = '../single_didaexe/'
python_runable = ' ../run_single.py '
# python_runable = ' ../python_test.py '


print(os.getcwd())
keyl = ['1']
# ,'1']
# ,'2','4']
member = [32]

# fix_member = 32
for tar,fix_member in zip(keyl,member):
	work_path = chrd + tar + '_cnn'+str(fix_member)
	namelist = 'namelists'+'.input'
	tar_exe = exe + 'dida_'+tar+'.exe'

	os.chdir(work_path)
	print(os.getcwd())

	cmd = 'salloc --exclusive -N 47 -n 188  -p cnCPU  --ntasks-per-node=4  mpirun -n 91  ' \
	      + tar_exe + ' ' + namelist + ' : -n 97  ~/.conda/envs/py36/bin/python3.6 '+ python_runable + str(fix_member) + ' | tee -a log.out'
	print(cmd)
	os.system(cmd)
	os.chdir('../')



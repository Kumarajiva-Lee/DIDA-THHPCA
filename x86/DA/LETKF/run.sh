rm -r result
mkdir result
cd result
mkdir actual
mkdir back
mkdir ana
mkdir param
cd ..

srun -N 1 -n 8 -w cn008 ./build/test.exe
#srun -N 1 -n 8 ./build/test.exe

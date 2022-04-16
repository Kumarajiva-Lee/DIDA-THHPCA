mpif90 -O3 -g -o letkf_core_fc.o -c letkf_core_half.F90
mpicc -O3 -g -I/usr/sw/yyzlib/hgemm/include -I../../../../LIB/gptl/include -mieee -mhost -o half_master.o -c letkf_core_half_master.c
mpicc -O3 -g -mieee -mslave -msimd -mftz -o half_slave.o -c letkf_core_half_slave.c
ar rc libletkf_half.a letkf_core_fc.o half_master.o half_slave.o

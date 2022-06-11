FC = mpif90
CC = sw9gcc

GPTL_DIR = /home/export/online1/mdt00/shisuan/swthhpca/DIDA-SHARE/3rdParty/gptl_sw/lib
XMATH_PATH = /home/export/online1/mdt00/shisuan/swthhpca/DIDA-SHARE/3rdParty/xMath-new
HGEMM_DIR = /usr/sw/yyzlib/hgemm

test_sw: master.o letkf_core.o letkf_core_fc.o half_master.o half_slave.o
	$(FC) $^ -g -o $@ -L$(XMATH_PATH) $(HGEMM_DIR)/libswhgemm.a $(GPTL_DIR)/libgptl.a -lxMath -lswperf -lm_slave
master.o: master.F90
	$(FC) -O3 -g -o $@ -c $<
letkf_core.o: letkf_core.F90
	$(FC) -g -mieee -o $@ -c $<
letkf_core_fc.o: letkf_core_half.F90
	$(FC) -O3 -g -o $@ -c $<
half_master.o: letkf_core_half_master.c
	$(CC) -O3 -g -I$(HGEMM_DIR)/include -O3 -mieee -mhost -o $@ -c $^
half_slave.o: letkf_core_half_slave.c
	$(CC) -O3 -g -mieee -mslave -msimd -mftz -o $@ -c $^

run: test_sw
	bsub -b -I -q q_test_2cab -n 1 -cgsp 64 -share_size 10000 -host_stack 256 ./test_sw


.PHONY: clean
clean :
	rm *.o test_sw

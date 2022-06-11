#This makefile is for acc optimization for letkf
ACC=swafort

MLAPA_DIR=/usr/sw/yyzlib/lapack-3.8.0
XMATH_DIR=/usr/sw/yyzlib/xMath-SACA
GPTL_DIR=/home/export/online1/mdt00/shisuan/swthhpca/DIDA-SHARE/3rdParty/gptl_sw/lib

letkf_core.o: letkf_core.F90
	$(ACC) $^ -cpp -mieee -msimd -c -o letkf_core.o

.PHONY: clean

clean:
	rm -f *.a *.o

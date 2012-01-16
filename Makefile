C = gfortran
O = -O2 
D = 
L = ./liblapack.a /usr/lib/atlas/libblas.a  
#L = -framework Accelerate
#I = -I ~/myprogs/include
all: test dgemm_test

clean:
	rm *.o test dgemm_test doublechksum 

ftdpotrf.o: ftdpotrf.f90
	$C $O $D -c ftdpotrf.f90
mydpotf2.o: mydpotf2.f90
	$C $O $D -c mydpotf2.f90
ftdgemm.o: ftdgemm.f90
	$C $O $D -c ftdgemm.f90

test: dpotrf_test.f90 ftdpotrf.o mydpotf2.o ftdgemm.o
	$C $O $D $I dpotrf_test.f90 ftdpotrf.o mydpotf2.o ftdgemm.o $L -o test

doublechksum: doblechksum.f90
	$C $O $D doblechksum.f90 mydpotf2.o $L -o doublechksum	 

dgemm_test: ftdgemm.o ftdgemm_test.f90
	$C $O $D ftdgemm_test.f90 ftdgemm.o $L -o dgemm_test


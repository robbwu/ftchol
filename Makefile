C = gfortran
O = -O3 
D =
L = ./liblapack.a /usr/lib/atlas/libblas.a  ~/myprogs/lib/libpapi.a
I = -I ~/myprogs/include
all: test

ftdpotrf.o: ftdpotrf.f90
	$C $O $D -c ftdpotrf.f90
mydpotf2.o: mydpotf2.f90
	$C $O $D -c mydpotf2.f90

test: dpotrf_test.f90 ftdpotrf.o mydpotf2.o
	$C $O $D $I dpotrf_test.f90 ftdpotrf.o mydpotf2.o $L -o test

doublechksum: doblechksum.f90
	$C $O $D doblechksum.f90 mydpotf2.o $L -o doublechksum	 


all: test

ftdpotrf.o: ftdpotrf.f90
	gfortran -c ftdpotrf.f90
mydpotf2.o: mydpotf2.f
	gfortran -c mydpotf2.f

test: dpotrf_test.f90 ftdpotrf.o mydpotf2.o
	gfortran dpotrf_test.f90 ftdpotrf.o mydpotf2.o ../lapack-3.4.0/liblapack.a /usr/lib/atlas/libblas.a -o test

doublechksum: doblechksum.f90
	gfortran doblechksum.f90 mydpotf2.o ../lapack-3.4.0/liblapack.a /usr/lib/atlas/libblas.a -o doublechksum	 


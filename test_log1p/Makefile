CC=gcc
FC=gfortran

test_log1p: test_log1p.o
	$(FC) -o test_log1p test_log1p.o

test_log1p.o:
	$(FC) -c -g -fPIC test_log1p.f90

clean:
	-rm -f *.o *.mod test_log1p

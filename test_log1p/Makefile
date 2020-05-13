test_log1p: test_log1p.o
	gfortran -o test_log1p test_log1p.o

test_log1p.o:
	gfortran -c -fPIC test_log1p.f90

clean:
	-rm -f *.o *.mod test_log1p

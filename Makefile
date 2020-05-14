libsymmetry.so: symmetry.o
	icc -shared -o libsymmetry.so symmetry.o -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

symmetry.o: symmetry.c
	icc -c -fPIC -I${MKLROOT}/include -o symmetry.o symmetry.c

clean:
	-rm -vf *.so *.o *.pyc



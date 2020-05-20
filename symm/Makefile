libsymm.so: equiv_map.o
	gcc -shared -o libsymm.so equiv_map.o 

equiv_map.o: equiv_map.c
	gcc -g -c -fPIC -o equiv_map.o equiv_map.c

clean:
	-rm -vf *.so *.o



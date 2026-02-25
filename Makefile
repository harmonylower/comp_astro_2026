FC=gfortran
FFLAGS=-c
SRC=say_hello.f90 learning.f90
OBJ=$(SRC:.f90=.o)

%.o: %.f90
	$(FC) $(FFLAGS) $<	-o $@

test: say_hello
	$(FC) -o $@ learning.f90 say_hello.o

say_hello:
	$(FC) $(FFLAGS) say_hello.f90



clean:
	rm -f *.o *.mod test

.PHONY: say_hello
FC = gfortran
FFLAGS = -O3 -Wextra -fdefault-real-8 -finit-real=nan
SRC = say_hello.f90 write_output.f90 setup.f90 density.f90 boundary.f90 SPH.f90 
OBJ = $(SRC:.f90=.o)


%.o: %.f90
	$(FC) $(FFLAGS)	-o $@ -c $<

sph: ${OBJ}
	$(FC) $(FFLAGS) -o $@ ${OBJ}

clean:
	rm -f *.o *.mod test


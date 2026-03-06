FC = gfortran
FFLAGS = -O3 -Wextra -fdefault-real-8 -finit-real=nan -Wall
SRC = say_hello.f90 eos.f90 kernal.f90 calcs.f90 setup.f90 boundary.f90 write_output.f90 leapfrog.f90 SPH.f90 
OBJ = $(SRC:.f90=.o)


%.o: %.f90
	$(FC) $(FFLAGS)	-o $@ -c $<

sph: ${OBJ}
	$(FC) $(FFLAGS) -o $@ ${OBJ}

clean:
	rm -f *.o *.mod test



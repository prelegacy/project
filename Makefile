FC=gfortran
FFLAGS=-O3 -fdefault-real-8 -ffpe-trap=invalid,zero,overflow -O -fcheck=all -g -fbacktrace
ifeq ($(DEBUG), yes)
  FFLAGS+=-fcheck=all -ffpe-trap=invalid,zero,overflow
endif
SRC=  main.f90 setup.f90 func.f90 heateq.f90
OBJ=${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

main: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

main.o: setup.o func.o heateq.o

heateq.o: func.o setup.o

clean:
	rm main *.o *.mod *.dat *.png *.txt

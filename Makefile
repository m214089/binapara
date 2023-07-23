FC     = gfortran
FFLAGS = -g -fbacktrace -ffpe-trap=invalid,zero,overflow -std=f2018 -Wpedantic -Wall -Wextra -Wunderflow -freal-4-real-8

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

all: newbinapara

newbinapara: newbinapara.o
	$(FC) $(FFLAGS) -o $@ $<

newbinapara.o: newbinapara.f90

clean:
	rm -rf newbinapara newbinapara.o

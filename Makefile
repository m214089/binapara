FC     = gfortran
FFLAGS = -std=f2018 -Wpedantic -Wall -Wextra -Wunderflow -ffpe-trap=invalid,zero,overflow,inexact

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

all: newbinapara

newbinapara: newbinapara.o
	$(FC) $(FFLAGS) -o $@ $<

newbinapara.o: newbinapara.f90

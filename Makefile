FC = gfortran
FFLAGS = -w -ffixed-line-length-none -fd-lines-as-comments -fno-automatic -fPIC -freal-4-real-8

all: bunema

bunema: bunema.o
	$(FC) $(FFLAGS) -o bunema bunema.o

bunema.o: bunema.f
	$(FC) $(FFLAGS) -c bunema.f

clean:
	rm -f *.o bunema

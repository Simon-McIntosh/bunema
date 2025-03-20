FC = gfortran
FFLAGS = -w -ffixed-line-length-none -fd-lines-as-comments -fno-automatic -fPIC -freal-4-real-8

all: bunema

bunema: main.o bunema.o
	$(FC) $(FFLAGS) -o bunema main.o bunema.o

main.o: main.f
	$(FC) $(FFLAGS) -c main.f

bunema.o: bunema.f
	$(FC) $(FFLAGS) -c bunema.f

clean:
	rm -f *.o bunema

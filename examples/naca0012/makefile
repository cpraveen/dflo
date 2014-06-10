FC=gfortran

grid:
	$(FC) -o struct2gmsh struct2gmsh.f
	./struct2gmsh

run:
	../../src-entropy/dflo m063.prm

all: grid run

clean:
	rm -f *.plt

allclean:
	rm -f *.plt *.xml naca.msh struct2gmsh

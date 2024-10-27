dflo
====

> This code is rather old, for more recent examples of DG codes, see
> https://github.com/cpraveen/fem/tree/master/dg1d
> https://github.com/cpraveen/fem/tree/master/dg2d

Discontinuous Galerkin solver for compressible flows. Some features of the code are

* cartesian and quadrilateral cells
* Qk basis: nodal Lagrange polynomials on Gauss points
* Pk basis: Legendre polynomials
* TVB limiter
* Positivity preserving limiter
* Flux functions: Lax-Friedrichs, Roe, HLLC, KFVS

###Getting dflo

```git clone https://github.com/cpraveen/dflo```

### To compile code in "src"
You have to first compile deal.II with Trilinos library, UMFPACK, threads. Serial versions of these libraries are sufficient since the code in "src" is a serial code.

* ```cd dflo/src```
* You must set the variable ```DEAL_II_DIR``` to your deal.II installation directory. It is good to set this in your bashrc file.
* ```cmake .```
* ```make```

By default, this compiles a DEBUG version of the code which is good for code development but runs very slowly. In order to compile optimized version, do

```make release```

### To compile code in "src-mpi"
This is the mpi version which uses p4est. So you must install p4est and compile deal.II with p4est support. It does not require trilinos.

To install p4est, see this page

https://www.dealii.org/developer/external-libs/p4est.html

Obtain latest version of deal.II from github

```git clone https://github.com/dealii/dealii```

Change into dealii directory

```
cd dealii
mkdir build
cd build
```

Configure deal.II. A basic setup is given in the file ```dealii_mpi.sh```. Just run this inside ```build``` directory

```
./dealii_mpi.sh
```

Now compile deal.II and install it

```
make all
make install
```

After this you can compile dflo.

### Running the code
Many test cases are included in the ```examples``` directory. In most examples, a geo file is provided to generate the mesh. You need [Gmsh](http://geuz.org/gmsh) to generate the mesh as follows

```gmsh -2 grid.geo```

Then run dflo

```dflo input.prm```

Some parts of dflo are parallelized with threads. You can specify number of threads to use while starting dflo

```dflo input.prm 4```

which will run with 4 threads.

To run the mpi version, do

```mpirun -np 4 dflo input.prm```

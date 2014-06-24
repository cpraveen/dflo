dflo
====

Discontinuous Galerkin solver for compressible flows. Some features of the code are

* cartesian and quadrilateral cells
* nodal Lagrange basis on Gauss points
* TVB limiter
* Positivity preserving limiter
* Flux functions: Lax-Friedrichs, Roe, HLLC, KFVS

### To compile
* ```cd dflo/src```
* You must set the variable ```DEAL_II_DIR``` to your deal.II installation directory. It is good to set this in your bashrc file.
* ```cmake .```
* ```make```

By default, this compiles a DEBUG version of the code which is good for code development but runs very slowly. In order to compile optimized version, do

```make release```

### Running the code
Many test cases are included in the ```examples``` directory. In most examples, a geo file is provided to generate the mesh. You need [Gmsh](http://geuz.org/gmsh) to generate the mesh as follows

```gmsh -2 grid.geo```

Then run dflo

```dflo input.prm```

Some parts of dflo are parallelized with threads. You can specify number of threads to use while starting dflo

```dflo input.prm 4```

which will run with 4 threads.

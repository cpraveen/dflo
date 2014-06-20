dflo
====

Discontinuous Galerkin solver for compressible flows

### To compile ###
* cd dflo/src
* You must set the path to your deal.II installation. It is good to set this in your bashrc file.
* cmake .
* make

By default, this compiles a DEBUG version of the code which is good for code development but runs slowly. In order to compile optimized version, do

make release

dflo
====

Discontinuous Galerkin solver for compressible flows

### To compile ###
* cd dflo/src
* You must set the variable DEAL_II_DIR to your deal.II installation directory. It is good to set this in your bashrc file.
* cmake .
* make

By default, this compiles a DEBUG version of the code which is good for code development but runs very slowly. In order to compile optimized version, do

make release

ifort -Vaxlib -openmp -extend_source -mcmodel=large -i-dynamic -g -r8 -i8 -c galtogad.f90 createcool.f90 makehalo.f90 calculate_h0_fun.f90 routines.f90
ifort -Vaxlib -openmp -extend_source -mcmodel=large -i-dynamic -g -r8 -i8 galtogad.o createcool.o makehalo.o calculate_h0_fun.o routines.o -o galtogadrun
chmod u+x galtogadrun

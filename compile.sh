gfortran -c var.f90
gfortran main.f90 var.o -llapack

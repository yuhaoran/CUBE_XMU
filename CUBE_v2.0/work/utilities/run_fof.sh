gfortran -O3 -cpp -fopenmp -fcoarray=single -DPID ../main/parameters.f90 fof.f90 && time ./a.out

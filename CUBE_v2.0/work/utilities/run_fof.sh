gfortran -O3 -cpp -fopenmp -fcoarray=single -mcmodel=medium -DPID ../main/parameters.f90 CUBE_fof_halofinder.f90 -o fof.x && time ./fof.x

export FC='ifort'
export XFLAG_NO_OMP='-O3 -fpp -qopenmp -coarray=distributed -mcmodel=medium'
export XFLAG='-O3 -fpp -qopenmp -coarray=distributed -mcmodel=medium'
export OFLAG_NO_OMP=$XFLAG_NO_OMP' -c'
export OFLAG=$XFLAG' -c'
export FFTFLAG='-I/opt/intel/mkl/include/fftw/ -mkl'

export OMP_STACKSIZE=16000M
#export KMP_STACKSIZE=16000M
export OMP_NUM_THREADS=64
#export OMP_THREAD_LIMIT=4
export FOR_COARRAY_NUM_IMAGES=1
ulimit
ulimit -s unlimited

# run executable by ./a.out and set number of images by FOR_COARRAY_NUM_IMAGES
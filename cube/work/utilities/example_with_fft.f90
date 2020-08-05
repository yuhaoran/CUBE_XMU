#define macbook
program spinfield
  use omp_lib
  use parameters
  !use pencil_fft
  use iso_fortran_env, only : int64
  implicit none

  integer,parameter :: ngrid=16
  integer ii,jj,kk
  integer(8) plan_fft_fine,plan_ifft_fine
  real rho_f(ngrid+2,ngrid,ngrid)
  real vel(ngrid,ngrid,ngrid,3)
  real rvec(3)

  call create_cubefft_plan

  do kk=1,ngrid
  do jj=1,ngrid
  do ii=1,ngrid
    rvec=modulo([ii,jj,kk]+ngrid/2-1,ngrid)-ngrid/2;
    vel(ii,jj,kk,1)= rvec(2)*exp(-sum(rvec(1:2)**2)/2);
    vel(ii,jj,kk,2)=-rvec(1)*exp(-sum(rvec(1:2)**2)/2);
    vel(ii,jj,kk,3)=0;
  enddo
  enddo
  enddo
  vel(1,1,1,:)=0;













call destroy_cubefft_plan

contains

  subroutine create_cubefft_plan
    use,intrinsic :: ISO_C_BINDING
    use omp_lib
    implicit none
    save
    !include 'fftw3.f'
    include 'fftw3.f03'
    integer istat,icore

#ifndef macbook
    call sfftw_init_threads(istat)
    print*, 'sfftw_init_threads status',istat
    icore=omp_get_max_threads()
    print*, 'omp_get_max_threads() =',icore
    !call sfftw_plan_with_nthreads(icore)

    call sfftw_plan_with_nthreads(64)
#endif

    call sfftw_plan_dft_r2c_3d(plan_fft_fine,ngrid,ngrid,ngrid,rho_f,rho_f,FFTW_MEASURE)
    call sfftw_plan_dft_c2r_3d(plan_ifft_fine,ngrid,ngrid,ngrid,rho_f,rho_f,FFTW_MEASURE)
  endsubroutine create_cubefft_plan

  subroutine destroy_cubefft_plan
    use,intrinsic :: ISO_C_BINDING
    implicit none
    save
    !include 'fftw3.f'
    include 'fftw3.f03'
    call sfftw_destroy_plan(plan_fft_fine)
    call sfftw_destroy_plan(plan_ifft_fine)
#ifndef macbook
    call fftw_cleanup_threads()
#endif
  endsubroutine destroy_cubefft_plan

endprogram

!gfortran -O3 -cpp -fopenmp -fcoarray=single -DFFTFINE -Dsigma8 -DPID ../main/parameters.f90 gwspin_reco.f90 -I/usr/local/Cellar/fftw/3.3.8_1/include/ -L/usr/local/Cellar/fftw/3.3.8_1/lib/ -lfftw3f -lm -ldl
#define macbook
program gwspin_reco
  use omp_lib
  use parameters
  use iso_fortran_env, only : int64
  implicit none

  integer,parameter :: ngrid=nf
  integer t1,t2,tt1,tt2,ttt1,ttt2,t_rate,kg,jg,ig,ii,jj,cur_checkpoint
  integer i,j,k,n_rsmall,n_ratio,nmassbin,itemp,l
  integer(8) plan_fft_fine,plan_ifft_fine
  real kr,kx,ky,kz,pow,r_filter,spin_u(3),spin_v(3),spin_r(3),spin_e(3,3),hpos1(3),hpos0(3)
  real t11,t22,t33,t12,t23,t31,tsmall(3,3),tlarge(3,3),torque(3,3)
  real spin(3,ngrid,ngrid,ngrid),rho_f(ngrid+2,ngrid,ngrid)
  complex crho_f(ngrid/2+1,ngrid,ngrid)
  complex,allocatable :: phi_k(:,:,:)
  real,allocatable :: phi(:,:,:),phi_large(:,:,:)
  real,allocatable :: r_small(:)
  
  equivalence(rho_f,crho_f)

  !call omp_set_num_threads(ncore)
  !call geometry
  image=1
  open(16,file='../main/z_checkpoint.txt',status='old')
  do i=1,nmax_redshift
    read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
  enddo
  71 n_checkpoint=i-1
  close(16)
  cur_checkpoint=1
  sim%cur_checkpoint=cur_checkpoint
  print*, ''
  print*, 'program gwspin_reco on single node'
  print*, 'on',ncore,' cores'
  print*, 'Resolution ngrid =', ngrid
  print*, 'at redshift=',z_checkpoint(cur_checkpoint)
  print*, 'Box size', box
  print*, 'output: ', opath
  print*, '-----------------------------------------'
  sync all
  print*,''
  print*,'Creating FFT plans'
  call system_clock(t1,t_rate)
  call create_cubefft_plan
  call system_clock(t2,t_rate)
  print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';
  sync all

  ! construct mass bins
  n_rsmall=10
  n_ratio=1

  ! construct scale bins
  allocate(r_small(n_rsmall))
  do i=1,n_rsmall
    r_small(i)=1.0+0.4*(i-1)
  enddo

  open(11,file=output_name('phik'),status='old',action='read',access='stream')
  read(11) crho_f
  close(11)
  allocate(phi_k(ngrid/2+1,ngrid,ngrid))
  allocate(phi(0:ngrid+1,0:ngrid+1,0:ngrid+1))
  allocate(phi_large(0:ngrid+1,0:ngrid+1,0:ngrid+1))
  phi_k=crho_f
  print*,'phi_k(k=0)=',phi_k(1,1,1)

  do ii=1,n_rsmall
    print*, jj,'/',n_ratio,' rs, ratio'
    call system_clock(tt1,t_rate)
    call correlate_spin(r_small(ii),1.1) ! loop over mass bins
    call system_clock(tt2,t_rate)
    print*, '  elapsed time =',real(tt2-tt1)/t_rate,'secs';
  enddo
  deallocate(r_small)
  call destroy_cubefft_plan

contains


  subroutine correlate_spin(r_small,ratio_scale)
    implicit none
    save
    real r_small,ratio_scale

    print*,'r [Mpc/h], r_ratio =', r_small, ratio_scale
    call gaussian_fourier_filter(phi_k,r_small)
    call sfftw_execute(plan_ifft_fine)
    phi(1:ngrid,1:ngrid,1:ngrid)=rho_f(:ngrid,:,:)
    call buffer_1layer(phi)

    call gaussian_fourier_filter(phi_k,r_small*ratio_scale)
    call sfftw_execute(plan_ifft_fine)
    phi_large(1:ngrid,1:ngrid,1:ngrid)=rho_f(:ngrid,:,:)
    call buffer_1layer(phi_large)

    call spinfield !(phi,phi_large,spin)
    print*,'b1,b2,b3 ='
    print*,sum(1d0*spin(1,:,:,:)**2)/ngrid/ngrid/ngrid
    print*,sum(1d0*spin(2,:,:,:)**2)/ngrid/ngrid/ngrid
    print*,sum(1d0*spin(3,:,:,:)**2)/ngrid/ngrid/ngrid
    print*,''
  endsubroutine

  real function ccc(vec_i,vec2)
    implicit none
    real vec_i(3),vec2(3)
    vec_i=vec_i/norm2(vec_i)
    vec2=vec2/norm2(vec2)
    ccc=sum(vec_i*vec2)
  endfunction

  subroutine buffer_1layer(phi)
    implicit none
    real phi(0:ngrid+1,0:ngrid+1,0:ngrid+1)
    phi(0,:,:)=phi(ngrid,:,:); sync all
    phi(ngrid+1,:,:)=phi(1,:,:); sync all
    phi(:,0,:)=phi(:,ngrid,:); sync all
    phi(:,ngrid+1,:)=phi(:,1,:); sync all
    phi(:,:,0)=phi(:,:,ngrid); sync all
    phi(:,:,ngrid+1)=phi(:,:,1); sync all
  endsubroutine

  subroutine spinfield !(phi,phi_large,spin)
    implicit none
    save
    do k=1,ngrid
    do j=1,ngrid
    do i=1,ngrid
      tsmall(1,1)=phi(i+1,j,k)-2*phi(i,j,k)+phi(i-1,j,k)
      tsmall(2,2)=phi(i,j+1,k)-2*phi(i,j,k)+phi(i,j-1,k)
      tsmall(3,3)=phi(i,j,k+1)-2*phi(i,j,k)+phi(i,j,k-1)
      tsmall(1,2)=(phi(i+1,j+1,k)+phi(i-1,j-1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k))/4
      tsmall(2,3)=(phi(i,j+1,k+1)+phi(i,j-1,k-1)-phi(i,j+1,k-1)-phi(i,j-1,k+1))/4
      tsmall(3,1)=(phi(i+1,j,k+1)+phi(i-1,j,k-1)-phi(i+1,j,k-1)-phi(i-1,j,k+1))/4
      tsmall(2,1)=tsmall(1,2)
      tsmall(3,2)=tsmall(2,3)
      tsmall(1,3)=tsmall(3,1)

      tlarge(1,1)=phi_large(i+1,j,k)-2*phi_large(i,j,k)+phi_large(i-1,j,k)
      tlarge(2,2)=phi_large(i,j+1,k)-2*phi_large(i,j,k)+phi_large(i,j-1,k)
      tlarge(3,3)=phi_large(i,j,k+1)-2*phi_large(i,j,k)+phi_large(i,j,k-1)
      tlarge(1,2)=(phi_large(i+1,j+1,k)+phi_large(i-1,j-1,k)-phi_large(i+1,j-1,k)-phi_large(i-1,j+1,k))/4
      tlarge(2,3)=(phi_large(i,j+1,k+1)+phi_large(i,j-1,k-1)-phi_large(i,j+1,k-1)-phi_large(i,j-1,k+1))/4
      tlarge(3,1)=(phi_large(i+1,j,k+1)+phi_large(i-1,j,k-1)-phi_large(i+1,j,k-1)-phi_large(i-1,j,k+1))/4
      tlarge(2,1)=tlarge(1,2)
      tlarge(3,2)=tlarge(2,3)
      tlarge(1,3)=tlarge(3,1)

      torque=-matmul(tsmall,tlarge)
      spin(1,i,j,k)=-torque(2,3)+torque(3,2)
      spin(2,i,j,k)=-torque(3,1)+torque(1,3)
      spin(3,i,j,k)=-torque(1,2)+torque(2,1)
      spin(:,i,j,k)=spin(:,i,j,k)/norm2(spin(:,i,j,k))
    enddo
    enddo
    enddo

  endsubroutine

  subroutine gaussian_fourier_filter(phi_k,r_filter)
    ! apply Gaussian window function with radius r_filter to Fourier space phi_k
    ! returns Fourier space crho_f
    implicit none
    complex phi_k(ngrid/2+1,ngrid,ngrid)
    real r_filter
    crho_f=0
    do k=1,ngrid
    do j=1,ngrid
    do i=1,ngrid/2+1
      kg=k
      jg=j
      ig=i
      kz=mod(kg+ngrid/2-1,ngrid)-ngrid/2
      ky=mod(jg+ngrid/2-1,ngrid)-ngrid/2
      kx=ig-1
      kr=sqrt(kx**2+ky**2+kz**2)
      kr=max(kr,1.0)
      kr=2*pi*kr/box
      pow=exp(-kr**2*r_filter**2/2)**0.25 ! apply E-mode window function
      crho_f(i,j,k)=phi_k(i,j,k)*pow
    enddo
    enddo
    enddo
    crho_f(1,1,1)=0 ! DC frequency
  endsubroutine

  subroutine tophat_fourier_filter(phi_k,r_filter)
    ! apply tophat window function with radius r_filter to Fourier space phi_k
    ! returns Fourier space crho_f
    implicit none
    integer,parameter :: n_dist=ngrid/2*sqrt(3.)+1
    integer i_dist(3),ibin
    real r_small,ratio_scale,r_dist,r_filter
    complex phi_k(ngrid/2+1,ngrid,ngrid)

    ! construct real space tophat
    do k=1,ngrid
    do j=1,ngrid
    do i=1,ngrid
      i_dist=mod((/i,j,k/)+ngrid/2-1,ng)-ngrid/2
      r_dist=norm2(real(i_dist))*box/ngrid
      rho_f(i,j,k)=merge(1.0,0.0,r_dist<=r_filter)
    enddo
    enddo
    enddo
    !open(11,file=output_name('tophat'),status='replace',access='stream')
    !write(11) r3
    !close(11)
    !call pencil_fft_forward ! get complex tophat
    call sfftw_execute(plan_fft_fine)
    !open(11,file=output_name('ctophat'),status='replace',access='stream')
    !write(11) crho_f
    !close(11)
    crho_f=phi_k*crho_f
  endsubroutine

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

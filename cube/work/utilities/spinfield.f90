#define macbook
program helicityfield
  use omp_lib
  use parameters
  !use pencil_fft
  use iso_fortran_env, only : int64
  implicit none

  integer,parameter :: ngrid=nf
  integer,parameter :: ny=ngrid/2 ! Nyquist frequency
  integer,parameter :: nbin=nint(ny*sqrt(3.))

  integer ii,jj,kk
  integer(8) plan_fft_fine,plan_ifft_fine
  real rho_f(ngrid+2,ngrid,ngrid)
  complex crho_f(ngrid/2+1,ngrid,ngrid),phi_k(ngrid/2+1,ngrid,ngrid)
  real phi(0:ngrid+1,0:ngrid+1,0:ngrid+1)
  real phi_large(0:ngrid+1,0:ngrid+1,0:ngrid+1)
  real tsmall(3,3),tlarge(3,3),torque(3,3)
  real spin(ngrid,ngrid,ngrid,3)
  real,dimension(ngrid,ngrid,ngrid,3) :: vel,velE,velL,velR
  real,dimension(ngrid+2,ngrid,ngrid,3) :: velk_real,velkE_real,velkL_real,velkR_real
  complex,dimension(ngrid/2+1,ngrid,ngrid,3) :: velk,velkE,velkL,velkR
  real rvec(3),kvec(3),k1,k2,k3,k12,k23,k31,k123
  complex proj_E(3,3),proj_L(3,3),proj_R(3,3)

  real xip(10,nbin)[*]

  equivalence(rho_f,crho_f)
  equivalence(velk_real,velk)
  equivalence(velkE_real,velkE)
  equivalence(velkL_real,velkL)
  equivalence(velkR_real,velkR)


  print*,'helicityfield'
  call geometry
  image=1
  open(16,file='../main/z_checkpoint.txt',status='old')
  do ii=1,nmax_redshift
    read(16,end=71,fmt='(f8.4)') z_checkpoint(ii)
  enddo
  71 n_checkpoint=ii-1
  close(16)
  sim%cur_checkpoint=1 ! read halos first

  print*,'call create_cubefft_plan'
  call create_cubefft_plan

! simple test field
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

! LPT velocity field
print*,'read potential',output_name('phi1')
open(11,file=output_name('phi1'),access='stream')
read(11) rho_f(:ngrid,:,:)
close(11)

call sfftw_execute(plan_fft_fine)
phi_k=crho_f

call gaussian_fourier_filter(phi_k,2.0)
call sfftw_execute(plan_ifft_fine); rho_f=rho_f/real(ngrid)/real(ngrid)/real(ngrid);
phi(1:ngrid,1:ngrid,1:ngrid)=rho_f(:ngrid,:,:)
call buffer_1layer(phi)

call gaussian_fourier_filter(phi_k,4.0)
call sfftw_execute(plan_ifft_fine); rho_f=rho_f/real(ngrid)/real(ngrid)/real(ngrid);
phi_large(1:ngrid,1:ngrid,1:ngrid)=rho_f(:ngrid,:,:)
call buffer_1layer(phi_large)

!open(11,file=output_name('phi'),status='replace',access='stream')
!write(11) phi(1:ngrid,1:ngrid,1:ngrid)
!close(11)
!open(11,file=output_name('phi_large'),status='replace',access='stream')
!write(11) phi_large(1:ngrid,1:ngrid,1:ngrid)
!close(11)
print*,'call spinfield'
call spinfield !(phi,phi_large,spin)
open(11,file=output_name('spin_2_4'),status='replace',access='stream')
write(11) spin
close(11)


  do ii=1,3
    rho_f(:ngrid,:,:)=spin(:,:,:,ii)
    call sfftw_execute(plan_fft_fine)
    velk_real(:,:,:,ii)=rho_f
  enddo

  velkE=0; velkL=0; velkR=0;
  ! Fourier space
  do kk=1,ngrid
  do jj=1,ngrid
  do ii=1,ngrid/2+1
    kvec=modulo([ii,jj,kk]+ngrid/2-1,ngrid)-ngrid/2;
    k1=kvec(1); k2=kvec(2); k3=kvec(3);
    k12=sqrt(k1**2+k2**2); k23=sqrt(k2**2+k3**2); k31=sqrt(k3**2+k1**2);
    k123=sqrt(k1**2+k2**2+k3**2);

    if (k123==0) cycle
    proj_E(1,1)=k1**2/k123**2; proj_E(2,2)=k2**2/k123**2; proj_E(3,3)=k3**2/k123**2;
    proj_E(1,2)=k1*k2/k123**2; proj_E(2,1)=proj_E(1,2);
    proj_E(2,3)=k2*k3/k123**2; proj_E(3,2)=proj_E(2,3);
    proj_E(3,1)=k3*k1/k123**2; proj_E(1,3)=proj_E(3,1);
    proj_L(1,1)=k23**2/k123**2/2; proj_L(2,2)=k31**2/k123**2/2; proj_L(3,3)=k12**2/k123**2/2;
    proj_L(1,2)=cmplx(-k1*k2/k123**2/2,k3/k123/2); proj_L(2,1)=conjg(proj_L(1,2));
    proj_L(2,3)=cmplx(-k2*k3/k123**2/2,k1/k123/2); proj_L(3,2)=conjg(proj_L(2,3));
    proj_L(3,1)=cmplx(-k3*k1/k123**2/2,k2/k123/2); proj_L(1,3)=conjg(proj_L(3,1));

    proj_R=transpose(proj_L);

    velkE(ii,jj,kk,:)=matmul(proj_E,velk(ii,jj,kk,:));
    velkL(ii,jj,kk,:)=matmul(proj_L,velk(ii,jj,kk,:));
    velkR(ii,jj,kk,:)=matmul(proj_R,velk(ii,jj,kk,:));

  enddo
  enddo
  enddo

  ! inverse transform
  do ii=1,3
    rho_f=velkE_real(:,:,:,ii)
    call sfftw_execute(plan_ifft_fine); rho_f=rho_f/real(ngrid)/real(ngrid)/real(ngrid);
    velE(:,:,:,ii)=rho_f(:ngrid,:,:)

    rho_f=velkL_real(:,:,:,ii)
    call sfftw_execute(plan_ifft_fine); rho_f=rho_f/real(ngrid)/real(ngrid)/real(ngrid);
    velL(:,:,:,ii)=rho_f(:ngrid,:,:)

    rho_f=velkR_real(:,:,:,ii)
    call sfftw_execute(plan_ifft_fine); rho_f=rho_f/real(ngrid)/real(ngrid)/real(ngrid);
    velR(:,:,:,ii)=rho_f(:ngrid,:,:)
  enddo

  open(11,file=output_name('spinfields'),status='replace',access='stream')
  write(11) spin,velE,velL,velR
  close(11)

  print*,'call vec_power'
  call vec_power(xip,spin)
  if (head) then
    open(15,file=output_name('velpower'),status='replace',access='stream')
    write(15) xip
    close(15)
  endif

  print*,'call inv_curl'
  call inv_curl
  open(11,file=output_name('velofspin'),status='replace',access='stream')
  write(11) spin
  close(11)




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

  subroutine vec_power(xip,vec1)
    ! decompose vec1 into E, L, R components
    ! compute helicity power spectrum
    implicit none

    integer ibin,idim,ii,jj,kk
    real sincx,sincy,sincz,sinc
    real ampvv,ampEE,ampLL,ampRR
    real vec1(ngrid,ngrid,ngrid,3),vk_real(ngrid+2,ngrid,ngrid,3)
    complex vk(ngrid/2+1,ngrid,ngrid,3),vkE(3),vkL(3),vkR(3)
    real xi(10,0:nbin),xip(10,nbin)[*]

    equivalence(vk_real,vk)
    xi=0
    do idim=1,3
      rho_f(:ngrid,:,:)=vec1(:,:,:,idim);
      call sfftw_execute(plan_fft_fine)
      vk_real(:,:,:,idim)=rho_f
    enddo

    do kk=1,ngrid
    do jj=1,ngrid
    do ii=1,ngrid/2+1
      !print*,ii,jj,kk
      kvec=modulo([ii,jj,kk]+ngrid/2-1,ngrid)-ngrid/2;
      k1=kvec(1); k2=kvec(2); k3=kvec(3);
      k12=sqrt(k1**2+k2**2); k23=sqrt(k2**2+k3**2); k31=sqrt(k3**2+k1**2);
      k123=sqrt(k1**2+k2**2+k3**2);
      sincx=merge(1.0,sin(pi*k1/ngrid)/(pi*k1/ngrid),k1==0.0)
      sincy=merge(1.0,sin(pi*k2/ngrid)/(pi*k2/ngrid),k2==0.0)
      sincz=merge(1.0,sin(pi*k3/ngrid)/(pi*k3/ngrid),k3==0.0)
      sinc=sincx*sincy*sincz
      if (k123/=0) then
        proj_E(1,1)=k1**2/k123**2; proj_E(2,2)=k2**2/k123**2; proj_E(3,3)=k3**2/k123**2;
        proj_E(1,2)=k1*k2/k123**2; proj_E(2,1)=proj_E(1,2);
        proj_E(2,3)=k2*k3/k123**2; proj_E(3,2)=proj_E(2,3);
        proj_E(3,1)=k3*k1/k123**2; proj_E(1,3)=proj_E(3,1);
        proj_L(1,1)=k23**2/k123**2/2; proj_L(2,2)=k31**2/k123**2/2; proj_L(3,3)=k12**2/k123**2/2;
        proj_L(1,2)=cmplx(-k1*k2/k123**2/2,k3/k123/2); proj_L(2,1)=conjg(proj_L(1,2));
        proj_L(2,3)=cmplx(-k2*k3/k123**2/2,k1/k123/2); proj_L(3,2)=conjg(proj_L(2,3));
        proj_L(3,1)=cmplx(-k3*k1/k123**2/2,k2/k123/2); proj_L(1,3)=conjg(proj_L(3,1));
      endif
      proj_R=transpose(proj_L);

      vkE=matmul(proj_E,vk(ii,jj,kk,:));
      vkL=matmul(proj_L,vk(ii,jj,kk,:));
      vkR=matmul(proj_R,vk(ii,jj,kk,:));
      ampvv=sum(vk(ii,jj,kk,:)*conjg(vk(ii,jj,kk,:)))
      ampEE=sum(vkE*conjg(vkE))
      ampLL=sum(vkL*conjg(vkL))
      ampRR=sum(vkR*conjg(vkR))
      ibin=nint(k123)

      xi(1,ibin)=xi(1,ibin)+1 ! number count
      xi(2,ibin)=xi(2,ibin)+k123 ! k count
      xi(3,ibin)=xi(3,ibin)+ampvv/(sinc**4.0)*4*pi*k123**3
      xi(4,ibin)=xi(4,ibin)+ampEE/(sinc**4.0)*4*pi*k123**3
      xi(5,ibin)=xi(5,ibin)+ampLL/(sinc**4.0)*4*pi*k123**3
      xi(6,ibin)=xi(6,ibin)+ampRR/(sinc**4.0)*4*pi*k123**3

    enddo
    enddo
    enddo
    xip=xi(:,1:);
    sync all

    ! co_sum
    if (head) then
      do ii=2,nn**3
        xip=xip+xip(:,:)[ii]
      enddo
    endif
    sync all

    ! broadcast
    xip=xip(:,:)[1]
    sync all

    xip(2,:)=xip(2,:)/xip(1,:)*(2*pi)/box;
    xip(3,:)=xip(3,:)/xip(1,:);
    xip(4,:)=xip(4,:)/xip(1,:);
    xip(5,:)=xip(5,:)/xip(1,:);
    xip(6,:)=xip(6,:)/xip(1,:);

  endsubroutine

  subroutine inv_curl
    ! decompose vec1 into E, L, R components
    ! compute helicity power spectrum
    implicit none

    integer ii,jj,kk
    real vk_real(ngrid+2,ngrid,ngrid,3)
    complex vk(ngrid/2+1,ngrid,ngrid,3),vkE(3),vkL(3),vkR(3)
    equivalence(vk_real,vk)

    rho_f=0;
    do ii=1,3
      rho_f(:ngrid,:,:)=spin(:,:,:,ii);
      call sfftw_execute(plan_fft_fine)
      vk_real(:,:,:,ii)=rho_f
    enddo

    do kk=1,ngrid
    do jj=1,ngrid
    do ii=1,ngrid/2+1
      !print*,ii,jj,kk
      if (ii*jj*kk==1) cycle
      kvec=modulo([ii,jj,kk]+ngrid/2-1,ngrid)-ngrid/2;
      k1=kvec(1); k2=kvec(2); k3=kvec(3);
      k12=sqrt(k1**2+k2**2); k23=sqrt(k2**2+k3**2); k31=sqrt(k3**2+k1**2);
      k123=sqrt(k1**2+k2**2+k3**2);
      proj_L(1,1)=k23**2/k123**2/2; proj_L(2,2)=k31**2/k123**2/2; proj_L(3,3)=k12**2/k123**2/2;
      proj_L(1,2)=cmplx(-k1*k2/k123**2/2,k3/k123/2); proj_L(2,1)=conjg(proj_L(1,2));
      proj_L(2,3)=cmplx(-k2*k3/k123**2/2,k1/k123/2); proj_L(3,2)=conjg(proj_L(2,3));
      proj_L(3,1)=cmplx(-k3*k1/k123**2/2,k2/k123/2); proj_L(1,3)=conjg(proj_L(3,1));
      proj_R=transpose(proj_L);
      !vk(ii,jj,kk,:) = (matmul(proj_R,vk(ii,jj,kk,:))-matmul(proj_L,vk(ii,jj,kk,:)))/k123;
      !vk(ii,jj,kk,:) = -matmul(proj_L,vk(ii,jj,kk,:))/k123;
      vk(ii,jj,kk,:) =  matmul(proj_R,vk(ii,jj,kk,:))/k123;
    enddo
    enddo
    enddo
    vk(1,1,1,:)=0
    do ii=1,3
      crho_f=vk(:,:,:,ii);
      call sfftw_execute(plan_ifft_fine)
      spin(:,:,:,ii)=rho_f(:ngrid,:,:)
    enddo
  endsubroutine

  subroutine gaussian_fourier_filter(phi_k,r_filter)
    ! apply Gaussian window function with radius r_filter to phi_k
    ! returns Fourier space crho_f
    implicit none
    complex phi_k(ngrid/2+1,ngrid,ngrid)
    real r_filter,pow,kz,ky,kx,kr
    integer i,j,k,kg,jg,ig
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
    integer i,j,k
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
      spin(i,j,k,1)=-torque(2,3)+torque(3,2)
      spin(i,j,k,2)=-torque(3,1)+torque(1,3)
      spin(i,j,k,3)=-torque(1,2)+torque(2,1)
      !spin(:,i,j,k)=spin(:,i,j,k)/norm2(spin(:,i,j,k))

    enddo
    enddo
    enddo

  endsubroutine

endprogram

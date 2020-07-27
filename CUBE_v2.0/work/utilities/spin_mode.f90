#define macbook
program lpt
  use omp_lib
  use parameters
  use halo_output
  !use pencil_fft
  use iso_fortran_env, only : int64
  implicit none
  save

  integer,parameter :: ngrid=ng
  real,parameter :: b_link=0.2 ! linking length
  character(4) b_link_string
  character(300) fn_phi,fn_out,fn_raw
  integer(4) t1,t2,tt1,tt2,t_rate,nhalo,ihalo,hgrid(3),ninfo
  integer(8) kg,jg,ig,ir,jj,imassbin,cur_checkpoint
  real kr,kx,ky,kz,pow,r_filter,rm,masstemp
  real spin_u(3),spin_v(3),spin_r(3),spin_e(3,3),hpos1(3),hpos0(3)
  integer(8) plan_fft_fine,plan_ifft_fine
  !real,allocatable :: rho_f(:,:,:)
  !complex,allocatable :: crho_f(:,:,:)
  real rho_f(ngrid+2,ngrid,ngrid)
  complex crho_f(ngrid/2+1,ngrid,ngrid)

  real,allocatable :: phi(:,:,:),phi_large(:,:,:),idsp(:,:,:,:)
  complex,allocatable :: phi_k(:,:,:)
  integer i,j,k,n_rsmall,nmassbin,itemp,l
  real t11,t22,t33,t12,t23,t31
  real tsmall(3,3),tlarge(3,3),torque(3,3)
  real spin(3,ngrid,ngrid,ngrid)
  real,allocatable :: theta(:,:),imass_info(:,:)
  real,allocatable :: corr_tqx(:,:,:),r_small(:)
  integer,allocatable :: ind(:,:),isort_mass(:),i1(:),i2(:),ii1,ii2
  equivalence(rho_f,crho_f)

  type(type_halo_catalog_header) fof_header
  type(type_halo_catalog_array),allocatable :: hcat(:)

  !call omp_set_num_threads(ncore)
  !call geometry
  write(b_link_string,'(f4.2)') b_link
  image=1
  open(16,file='../main/z_checkpoint.txt',status='old')
  do i=1,nmax_redshift
    read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
  enddo
  71 n_checkpoint=i-1
  close(16)
  cur_checkpoint=n_checkpoint ! read halos first
  sim%cur_checkpoint=cur_checkpoint
  print*, ''
  print*, 'spin mode correlation code on single node'
  print*, 'on',ncore,' cores'
  print*, 'Resolution ngrid =', ngrid
  print*, 'at redshift=',z_checkpoint(cur_checkpoint)
  print*, 'Box size', box
  print*, 'output: ', opath
  print*, '-----------------------------------------'

  sim%cur_checkpoint=1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  fn_phi=output_name('phik')
  sim%cur_checkpoint=n_checkpoint
  fn_out=output_name('spin_muq') !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  fn_raw=output_name('spin_muq_raw')

  sync all
  print*,''
  print*,'Creating FFT plans'
  call system_clock(t1,t_rate)
  call create_cubefft_plan
  call system_clock(t2,t_rate)
  print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';
  sync all
  ! open halo catalog
  print*,''
  sim%cur_checkpoint=n_checkpoint
  print*,'read halo catalog ',output_name('FoF_b'//b_link_string)
  open(11,file=output_name('FoF_b'//b_link_string),status='old',access='stream')
  read(11) fof_header
  nhalo=fof_header%nhalo; ninfo=fof_header%ninfo
  allocate(hcat(nhalo))
  read(11) hcat
  close(11)

  allocate(ind(3,nhalo),theta(3,nhalo))
  allocate(phi(0:ngrid+1,0:ngrid+1,0:ngrid+1))
  allocate(phi_large(0:ngrid+1,0:ngrid+1,0:ngrid+1))
  allocate(phi_k(ngrid/2+1,ngrid,ngrid))

  do ihalo=1,nhalo
    !ind(:,ihalo)=ceiling(hcat(ihalo)%q/real(ng_global)*real(ngrid)) ! indeces in phi
    ind(:,ihalo)=ceiling(hcat(ihalo)%q/real(ng_global)*real(ngrid)) ! indeces in phi !!!!!!!!!!!!
  enddo

 ! construct halo mass bins
  rm=2.0 ! mass bin ratio
  n_rsmall=25
  nmassbin=ceiling(log(hcat(1)%hmass/hcat(nhalo)%hmass)/log(rm))
  print*,'  FoF linking parameter =',fof_header%linking_parameter
  print*,'  nhalo =',nhalo
  print*,'  max,min of hcat%hmass =',hcat(1)%hmass,hcat(nhalo)%hmass
  print*,'  mass bin, rm =',rm
  print*,'  nmassbin =',nmassbin
  allocate(imass_info(6,nmassbin),i1(nmassbin),i2(nmassbin))
  !allocate(n_opt(n_rsmall))
  allocate(corr_tqx(n_rsmall,nmassbin,6))
  allocate(r_small(n_rsmall))
  !n_opt=[25,25,25,24,23,22,21,19,18,17,16,15];
  imassbin=nmassbin
  i2(nmassbin)=nhalo
  i1(1)=1
  itemp=1
  do ihalo=nhalo,1,-1 ! construct mass bin boundaries
    if (hcat(ihalo)%hmass<hcat(nhalo)%hmass*rm**itemp .or. imassbin==1) then
      i1(imassbin)=ihalo ! still in the bin
    else
      imassbin=imassbin-1 ! (must) assign to previous bin
      i2(imassbin)=ihalo
      i1(imassbin)=ihalo
      itemp=itemp+1
    endif
  enddo
  print*,'  boundaries ='
  print*, i1
  print*, i2
  ! construct scale bins in logspace
  do i=1,n_rsmall
    r_small(i)=0.1*1.2**(i-1)
  enddo
  print*,'  smoothing scale r ='
  print*,r_small

  sim%cur_checkpoint=1 ! read IC Fourier space gravitational potential
  open(11,file=fn_phi,status='old',access='stream')
  read(11) phi_k
  close(11)

  cur_checkpoint=n_checkpoint
  sim%cur_checkpoint=cur_checkpoint
  imass_info=0 ! to be computed later


  !call correlate_spin(3.0)  
  !open(12,file=fn_raw,status='replace',access='stream')
  !write(12) 3, nhalo
  !write(12) theta
  !close(12)
  !stop

  print*, jj,'/',' rs, t, q, x'
  do ir=1,n_rsmall
    call system_clock(tt1,t_rate)
    call correlate_spin(r_small(ir)) ! loop over mass bins
    call system_clock(tt2,t_rate)
    print*, '  elapsed time =',real(tt2-tt1)/t_rate,'secs';
  enddo
  print*,''

  do imassbin=1,nmassbin
    ii1=i1(imassbin)
    ii2=i2(imassbin)
    imass_info(1,imassbin)=minval(hcat(ii1:ii2)%hmass)
    imass_info(2,imassbin)=maxval(hcat(ii1:ii2)%hmass)
    imass_info(3,imassbin)=sum(hcat(ii1:ii2)%hmass)/(ii2-ii1+1)
    imass_info(4,imassbin)=ii2-ii1+1
    imass_info(5,imassbin)=ii1
    imass_info(6,imassbin)=ii2
  enddo
  open(11,file=fn_out,status='replace',access='stream')
  write(11) nmassbin,n_rsmall,imass_info(:,:),r_small,corr_tqx
  close(11)
!call destroy_penfft_plan
call destroy_cubefft_plan


contains


  subroutine correlate_spin(r_small)
    implicit none
    save
    real r_small

    call gaussian_fourier_filter(phi_k,r_small)
    call sfftw_execute(plan_ifft_fine)
    !print*,rho_f(1:3,1,1); stop
    phi(1:ngrid,1:ngrid,1:ngrid)=rho_f(:ngrid,:,:)
    call buffer_1layer(phi)
    !call gaussian_fourier_filter(phi_k,r_small*1.05)
    call gaussian_fourier_filter_k2(phi_k,r_small)
    call sfftw_execute(plan_ifft_fine)
    phi_large(1:ngrid,1:ngrid,1:ngrid)=rho_f(:ngrid,:,:)
    call buffer_1layer(phi_large)

    call spinfield !(phi,phi_large,spin)
    do ihalo=1,nhalo
      !print*,ind(:,ihalo)
      !print*,spin(:,ind(1,ihalo),ind(2,ihalo),ind(3,ihalo))
      !stop
      theta(1,ihalo)=ccc(hcat(ihalo)%jt,spin(:,ind(1,ihalo),ind(2,ihalo),ind(3,ihalo)))
      theta(2,ihalo)=ccc(hcat(ihalo)%jl,spin(:,ind(1,ihalo),ind(2,ihalo),ind(3,ihalo)))
      theta(3,ihalo)=ccc(hcat(ihalo)%je,spin(:,ind(1,ihalo),ind(2,ihalo),ind(3,ihalo)))
    enddo
    
    do imassbin=1,nmassbin
      ii1=i1(imassbin)
      ii2=i2(imassbin)
      corr_tqx(ir,imassbin,1)=sum(theta(1,ii1:ii2))/(ii2-ii1+1)
      corr_tqx(ir,imassbin,2)=sum(theta(2,ii1:ii2))/(ii2-ii1+1)
      corr_tqx(ir,imassbin,3)=sum(theta(3,ii1:ii2))/(ii2-ii1+1)
      corr_tqx(ir,imassbin,4)=sum(theta(1,ii1:ii2)**2)/(ii2-ii1+1)
      corr_tqx(ir,imassbin,5)=sum(theta(2,ii1:ii2)**2)/(ii2-ii1+1)
      corr_tqx(ir,imassbin,6)=sum(theta(3,ii1:ii2)**2)/(ii2-ii1+1)
    enddo
    print*, r_small, sum(theta,2)/nhalo
  endsubroutine

  real function ccc(vec_i,vec2)
    implicit none
    real vec_i(3),vec2(3)
    vec_i=vec_i/norm2(vec_i)
    vec2=vec2/norm2(vec2)
    ccc=sum(vec_i*vec2)
    !ccc=ccc**2 !!! unoriented correlation
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

    subroutine gaussian_fourier_filter_k2(phi_k,r_filter)
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
      crho_f(i,j,k)=-phi_k(i,j,k)*pow * kr**2 ! converted to density
    enddo
    enddo
    enddo
    crho_f(1,1,1)=0 ! DC frequency
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

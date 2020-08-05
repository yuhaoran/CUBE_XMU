#define macbook
!#define ELUCID
!#define reco_idsp
program lpt
  use omp_lib
  use parameters
  use halo_output
  !use pencil_fft
  use iso_fortran_env, only : int64
  implicit none
  save
#ifdef ELUCID
  integer,parameter :: ngrid=500 ! ELUCID grid number
#else
  integer,parameter :: ngrid=ng
#endif
  character(300) fn_phi,fn_out
  integer(4) t1,t2,tt1,tt2,t_rate,nhalo,ihalo,hgrid(3),ninfo
  integer(8) kg,jg,ig,ii,jj,imassbin,cur_checkpoint
  real kr,kx,ky,kz,pow,r_filter,rm,masstemp
  real spin_u(3),spin_v(3),spin_r(3),spin_e(3,3),hpos1(3),hpos0(3)
  integer(8) plan_fft_fine,plan_ifft_fine
  !real,allocatable :: rho_f(:,:,:)
  !complex,allocatable :: crho_f(:,:,:)
  real rho_f(ngrid+2,ngrid,ngrid)
  complex crho_f(ngrid/2+1,ngrid,ngrid)

  real,allocatable :: phi(:,:,:),phi_large(:,:,:),idsp(:,:,:,:)
  complex,allocatable :: phi_k(:,:,:)
  integer i,j,k,n_rsmall,n_ratio,nmassbin,itemp,l
  real t11,t22,t33,t12,t23,t31
  real tsmall(3,3),tlarge(3,3),torque(3,3)
  real spin(3,ngrid,ngrid,ngrid)
  real,allocatable :: theta(:,:),imass_info(:,:)
  real,allocatable :: corr_x(:,:,:),corr_q(:,:,:),corr_t(:,:,:),r_small(:),ratio_scale(:)
  integer,allocatable :: ind(:,:),isort_mass(:),i1(:),i2(:),ii1,ii2
  equivalence(rho_f,crho_f)

  type type_halo_catalog
    integer nhalo_tot,nhalo
    real den_odc
  endtype

  type type_halo_info
    real hpos(3)
    real mass_odc,radius_odc,v_disp
    real x_mean(3),v_mean(3),ang_mom(3),var_x(3),inertia(3,3)
    real q_mean(3),inertia_q(3,3),s_mean(3)
  endtype

  type(type_halo_info) halo_info
  type(type_halo_catalog) halo_catalog
  type(type_halo_catalog_header) halo_header
  type(type_halo_catalog_array),allocatable :: hcat(:)

  !call omp_set_num_threads(ncore)
  !call geometry
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
  if (.true.) then
    sim%cur_checkpoint=n_checkpoint
    !fn_phi=output_name('phik_E')
    fn_phi='../../output/universe3/image1/0.000_phik_E_1.bin'
    fn_out=output_name('lptcorr_e')
  else
    sim%cur_checkpoint=1
    fn_phi=output_name('phik')
    !fn_phi='../../output/universe3/image1/100.000_phik_1.bin'
    fn_out=output_name('lptcorr_i')
  endif

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
  print*,'read halo catalog ',output_name('fof')
  open(11,file=output_name('fof'),status='old',access='stream')
  read(11) halo_header
  nhalo=halo_header%nhalo; ninfo=halo_header%ninfo
  allocate(hcat(nhalo))
  read(11) hcat
  close(11)

  allocate(ind(3,nhalo),theta(3,nhalo))
  allocate(phi(0:ngrid+1,0:ngrid+1,0:ngrid+1))
  allocate(phi_large(0:ngrid+1,0:ngrid+1,0:ngrid+1))
  allocate(phi_k(ngrid/2+1,ngrid,ngrid))
#ifdef reco_idsp
  allocate(idsp(3,ngrid,ngrid,ngrid))
  open(10,file=output_name('idsp_c'),status='old',access='stream')
  read(10) idsp ! map from Eulerian space to Lagrangian space
  close(10)
#endif
  do ihalo=1,nhalo
#   ifdef reco_idsp
      hpos1 = hcat(ihalo)%s/real(ng_global)*real(ngrid)
      hgrid = ceiling(hpos1)
      hpos0 = hpos1 + idsp(:,hgrid(1),hgrid(2),hgrid(3))
      hpos0 = modulo(hpos0,real(ngrid))
      ind(:,ihalo) = ceiling(hpos0)
#   else
      ind(:,ihalo)=ceiling(hcat(ihalo)%q/real(ng_global)*real(ngrid)) ! indeces in phi
#   endif
  enddo

 ! construct halo mass bins
  rm=2.0 ! mass bin ratio
  n_rsmall=100
  n_ratio=1
  nmassbin=ceiling(log(hcat(1)%hmass/hcat(nhalo)%hmass)/log(rm))
  print*,'  FoF linking parameter =',halo_header%linking_parameter
  print*,'  nhalo =',nhalo
  print*,'  max,min of hcat%hmass =',hcat(1)%hmass,hcat(nhalo)%hmass
  print*,'  mass bin ratio, rm =',rm
  print*,'  nmassbin =',nmassbin
  allocate(imass_info(4,nmassbin),i1(nmassbin),i2(nmassbin))
  allocate(corr_x(n_rsmall,n_ratio,nmassbin),corr_q(n_rsmall,n_ratio,nmassbin),corr_t(n_rsmall,n_ratio,nmassbin))
  allocate(r_small(n_rsmall),ratio_scale(n_ratio))
  imassbin=nmassbin
  i2(nmassbin)=nhalo
  i1(1)=1
  itemp=1
  do ihalo=nhalo,1,-1 ! construct mass bin boundaries
    if (hcat(ihalo)%hmass<hcat(nhalo)%hmass*rm**itemp .or. imassbin==1) then
      i1(imassbin)=ihalo ! still in the bin
    else
      imassbin=imassbin-1 ! assign to previous bin
      i2(imassbin)=ihalo
      itemp=itemp+1
    endif
  enddo
  print*,'  boundaries ='
  print*, i1
  print*, i2
  ! construct scale bins
  do i=1,n_rsmall
    r_small(i)=0.2+0.2*(i-1)
  enddo
  do i=1,n_ratio
    ratio_scale(i)=1.0+0.2*(i-1)
  enddo
  print*,'  smoothing scale r ='
  print*,r_small
  print*,'  R/r'
  print*,ratio_scale

#ifdef ELUCID
  open(11,file='../../S500_5001/cxyz_251_500_500.bin',access='stream')
  read(11) rho_f
  close(11)
  ! covert to potential
  do k=1,ngrid
  do j=1,ngrid
  do i=1,ngrid/2+1
    kg=k
    jg=j
    ig=i
    kz=mod(kg+ngrid/2-1,ngrid)-ngrid/2
    ky=mod(jg+ngrid/2-1,ngrid)-ngrid/2
    kx=ig-1
    kz=2*sin(pi*kz/ngrid)
    ky=2*sin(pi*ky/ngrid)
    kx=2*sin(pi*kx/ngrid)

    kr=kx**2+ky**2+kz**2
    kr=max(kr,1.0/ngrid**2)
    !kr=2*pi*kr/ngrid
    pow=-4*pi/kr
    crho_f(i,j,k)=crho_f(i,j,k)*pow
  enddo
  enddo
  enddo
  crho_f(1,1,1)=0
  phi_k=crho_f
#else
  sim%cur_checkpoint=1 ! read IC Fourier space gravitational potential
  open(11,file=fn_phi,status='old',access='stream')
  read(11) phi_k
  close(11)
#endif

  cur_checkpoint=n_checkpoint
  sim%cur_checkpoint=cur_checkpoint
  imass_info=0 ! to be computed later
  open(11,file=fn_out,status='replace',access='stream')
  write(11) nmassbin,n_rsmall,n_ratio,imass_info(:,:),r_small(:),ratio_scale(:)
  do jj=1,n_ratio
    print*, jj,'/',n_ratio,' rs, ratio, t, q, x'
    do ii=1,n_rsmall
      call system_clock(tt1,t_rate)
      call correlate_spin(r_small(ii),ratio_scale(jj)) ! loop over mass bins
      call system_clock(tt2,t_rate)
      print*, '  elapsed time =',real(tt2-tt1)/t_rate,'secs';
    enddo
    print*,''
  enddo

  do imassbin=1,nmassbin
    ii1=i1(imassbin)
    ii2=i2(imassbin)
    imass_info(1,imassbin)=minval(hcat(ii1:ii2)%hmass)
    imass_info(2,imassbin)=maxval(hcat(ii1:ii2)%hmass)
    imass_info(3,imassbin)=sum(hcat(ii1:ii2)%hmass)/(ii2-ii1+1)
    imass_info(4,imassbin)=ii2-ii1+1
    write(11) corr_t(:,:,imassbin)
    write(11) corr_q(:,:,imassbin)
    write(11) corr_x(:,:,imassbin)
  enddo
  rewind(11)
  write(11) nmassbin,n_rsmall,n_ratio,imass_info(:,:)
  close(11)
!call destroy_penfft_plan
call destroy_cubefft_plan


contains


  subroutine correlate_spin(r_small,ratio_scale)
    implicit none
    save
    real r_small,ratio_scale

    call gaussian_fourier_filter(phi_k,r_small)
    call sfftw_execute(plan_ifft_fine)
    !print*,rho_f(1:3,1,1); stop
    phi(1:ngrid,1:ngrid,1:ngrid)=rho_f(:ngrid,:,:)
    call buffer_1layer(phi)
    call gaussian_fourier_filter(phi_k,r_small*1.05)
    !call gaussian_fourier_filter_k2(phi_k,r_small)
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
    !print*,'a'
    do imassbin=1,nmassbin
      ii1=i1(imassbin)
      ii2=i2(imassbin)
      corr_t(ii,jj,imassbin)=sum(theta(1,ii1:ii2))/(ii2-ii1+1)
      corr_q(ii,jj,imassbin)=sum(theta(2,ii1:ii2))/(ii2-ii1+1)
      corr_x(ii,jj,imassbin)=sum(theta(3,ii1:ii2))/(ii2-ii1+1)
    enddo
    print*, r_small, ratio_scale, sum(theta,2)/nhalo
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

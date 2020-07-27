module pp_force
  use variables
  implicit none
  save


contains

subroutine ext_pp_force
  use omp_lib
  use variables
  implicit none
  save

  integer ilayer
  integer(8) ntest,ip_offset,ipll1,ipll2,it1,it2
  integer,allocatable :: ll(:),hoc(:,:,:)
  integer idxf,ipp1,ipp2
  integer,allocatable :: ipf(:),ift1(:)
  real,allocatable :: xf(:,:),vf(:,:),af(:,:)
  integer :: ipll1_list(10240*14*3)
  real :: x_tmp(10240*14*8,3),force_tmp(10240*14*8,3)
  integer :: ipll1_size, ipll2_size


  integer icellpp,icellpp2,ncellpp,ijk(3,62)

  if (head) then
    print*, ''
    print*, 'ext_pp_force'
    print*, '  pp_range=',int(pp_range,1)
    print*, '  ext_pp_force over',int(nnt**3,2),'tiles'
    print*, '  np_pp_max=',np_pp_max
    call system_clock(tt1,t_rate)
  endif
  nth=omp_get_max_threads()
  print*,'  max num_threads =',nth
  allocate(ll(np_pp_max),hoc(1-ncell:nft+ncell,1-ncell:nft+ncell,1-ncell:nft+ncell))
  allocate(ipf(np_pp_max),ift1(np_pp_max/800),xf(3,np_pp_max),vf(3,np_pp_max),af(3,np_pp_max))
  if (pp_range==1) then
    ncellpp=13
    ijk(:,1)=(/-1,-1,-1/) ! 9 cells in z-
    ijk(:,2)=(/ 0,-1,-1/)
    ijk(:,3)=(/ 1,-1,-1/)
    ijk(:,4)=(/-1, 0,-1/)
    ijk(:,5)=(/ 0, 0,-1/)
    ijk(:,6)=(/ 1, 0,-1/)
    ijk(:,7)=(/-1, 1,-1/)
    ijk(:,8)=(/ 0, 1,-1/)
    ijk(:,9)=(/ 1, 1,-1/)
    ijk(:,10)=(/-1,-1, 0/) ! 3 cells in y-
    ijk(:,11)=(/ 0,-1, 0/)
    ijk(:,12)=(/ 1,-1, 0/)
    ijk(:,13)=(/-1, 0, 0/) ! 1 cell in x-
  elseif (pp_range==2) then
    ncellpp=62
    ! tbd
  endif

  itest1=0
  f2_max_pp=0
  print*, '  sim%nplocal =',sum(rhoc(1:nt,1:nt,1:nt,:,:,:))
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    call system_clock(t1,t_rate)
    call system_clock(t_ll1,t_rate)
    ip_offset=idx_b_l(0,0,itx,ity,itz)+sum(rhoc(:-1,0,0,itx,ity,itz))

    hoc=0; ll=0; idxf=0;
    !!$omp paralleldo &
    !!$omp& default(shared) schedule(dynamic,1)&
    !!$omp& private(igz,igy,igx,np,nzero,lp,ip1,xvec) &
    !!$omp& private()
    do igz=0,nt+1
    do igy=0,nt+1
    do igx=0,nt+1
      np=rhoc(igx,igy,igz,itx,ity,itz)
      nzero=idx_b_r(igy,igz,itx,ity,itz)-sum(rhoc(igx:,igy,igz,itx,ity,itz))
      do lp=1,np ! loop over cdm particles
        ip1=nzero+lp
        xvec1=ncell*((/igx,igy,igz/)-1)+ncell*(int(xp(:,ip1)+ishift,izipx)+rshift)*x_resolution
        vreal=tan(pi*real(vp(:,ip1))/real(nvbin-1))/(sqrt(pi/2)/(sigma_vi*vrel_boost))
        idxf=idxf+1 ! index in the xf
        xf(:,idxf)=xvec1 ! make temp list xf
        vf(:,idxf)=vreal ! make temp list vf
        ipf(idxf)=ip1-ip_offset ! record index of ip in xf order, not consecutive due to buffer depth
        ivec1=floor(xvec1)+1
        ll(idxf)=hoc(ivec1(1),ivec1(2),ivec1(3))
        hoc(ivec1(1),ivec1(2),ivec1(3))=idxf
      enddo
    enddo
    enddo
    enddo
    !!$omp endparalleldo
    call system_clock(t_ll2,t_rate)

    call system_clock(t_af1,t_rate)
    af=0
    do ilayer=0,pp_range
    !$omp paralleldo &
    !$omp& default(shared) schedule(dynamic,1)&
    !$omp& private(igz,igy,igx,ipll1,ipll2,xvec21,rmag,rcut,pcut) &
    !$omp& private(force_pp,icellpp,icellpp2,ipll1_size,ipll2_size,ipll1_list,x_tmp,force_tmp)
    do igz=1+ilayer,nft+pp_range,pp_range+1
    do igy=1-pp_range,nft+pp_range
    do igx=1-pp_range,nft+pp_range
      ipll1_size=0
      ipll1=hoc(igx,igy,igz)
      do while (ipll1/=0)
        ipll1_size=ipll1_size+1
        ipll1_list(ipll1_size)=ipll1
        ipll1=ll(ipll1)
      enddo
      if( ipll1_size==0 ) then
        cycle
      endif

      ipll2_size=ipll1_size

      do icellpp=1,ncellpp
        ipll2=hoc(igx+ijk(1,icellpp),igy+ijk(2,icellpp),igz+ijk(3,icellpp))
        do while (ipll2/=0)
          ipll2_size=ipll2_size+1
          ipll1_list(ipll2_size)=ipll2
          ipll2=ll(ipll2)
        enddo
      enddo
      
      do icellpp=1,ipll2_size
        x_tmp(icellpp,:)=xf(:,ipll1_list(icellpp))
      enddo
      force_tmp(:ipll2_size+1,:)=0

      !kernel
      do icellpp=1,ipll1_size
        do icellpp2=icellpp+1,ipll2_size
          xvec21=x_tmp(icellpp2,:)-x_tmp(icellpp,:)
          rmag=sum(xvec21**2)+rsoft*rsoft
          rcut=rmag*rmag*rmag
          pcut=sim%mass_p_cdm/sqrt(rcut)
          force_pp=xvec21*pcut
          force_tmp(icellpp,:)=force_tmp(icellpp,:)+force_pp
          force_tmp(icellpp2,:)=force_tmp(icellpp2,:)-force_pp
        enddo
      enddo

      do icellpp=1,ipll2_size
        af(:,ipll1_list(icellpp))=af(:,ipll1_list(icellpp))+force_tmp(icellpp,:)
      enddo
    enddo
    enddo
    enddo
    !$omp endparalleldo
    enddo
    call system_clock(t_af2,t_rate)

    call system_clock(t_max1,t_rate)
    do igz=1,nft ! maxval on hoc(1:nft,1:nft,1:nft)
    do igy=1,nft
    do igx=1,nft
      ipll1=hoc(igx,igy,igz)
      do while (ipll1/=0)
        f2_max_pp=max(f2_max_pp,sum(af(:,ipll1)**2))
        ipll1=ll(ipll1)
      enddo
    enddo
    enddo
    enddo
    call system_clock(t_max2,t_rate)

    !print*, '    update vp'
    call system_clock(t_vp1,t_rate)
    !$omp paralleldo &
    !$omp& default(shared) schedule(dynamic,64)&
    !$omp& private(ipp1,ip1,vreal)
    do ipp1=1,idxf
      ip1=ipf(ipp1)+ip_offset
      vreal=vf(:,ipp1)+af(:,ipp1)*a_mid*dt/6/pi
      vp(:,ip1)=nint(real(nvbin-1)*atan(sqrt(pi/2)/(sigma_vi*vrel_boost)*vreal)/pi,kind=izipv)
    enddo
    !$omp endparalleldo
    call system_clock(t_vp2,t_rate)

    print*,'  f2_max_pp',f2_max_pp
    call system_clock(t2,t_rate)
    if (head) print*, '  tile:',int(itx,1),int(ity,1),int(itz,1),' elapsed time =',real(t2-t1)/t_rate,'secs'


    t0_ll=t0_ll+real(t_ll2-t_ll1)/t_rate
    t0_af=t0_af+real(t_af2-t_af1)/t_rate
    t0_max=t0_max+real(t_max2-t_max1)/t_rate
    t0_vp=t0_vp+real(t_vp2-t_vp1)/t_rate

  enddo
  enddo
  enddo !! itz
  !stop
  ! for reference, in pm.f90
  !dt_fine=sqrt( 1.0 / (sqrt(maxval(f2_max_fine))*a_mid*GG) )
  !dt_coarse=sqrt( real(ncell) / (sqrt(f2_max_coarse)*a_mid*GG) )
  !dt_pp=sqrt(0.1*rsoft) / max(sqrt(maxval(f2_max_pp))*a_mid*GG,1e-3)
  sim%dt_pp=1.0*sqrt(1.0) / (sqrt(f2_max_pp)*a_mid*GG)
  !sim%dt_pp=max(sim%dt_pp,0.01)
  sync all
  do i=1,nn**3
    sim%dt_pp=min(sim%dt_pp,sim[i]%dt_pp)
  enddo
  sync all

  deallocate(ll,hoc,ipf,ift1,xf,vf,af)

  if (head) then
    call system_clock(tt2,t_rate)
    print*,'  max of f_pp', sqrt(f2_max_pp)
    print*,'  dt_pp',sim%dt_pp
    print*,'  elapsed time =',real(tt2-tt1)/t_rate,'secs'
    print*,'  cumulative time: ll, af, max, vp'
    print*,'  ',t0_ll,t0_af,t0_max,t0_vp
    print*,''
  endif
  sync all

endsubroutine


endmodule
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

  integer(8) ntest,ip_offset,ipll1,ipll2,it1,it2,ilayer
  integer hoc(1-ncell:nft+ncell,1-ncell:nft+ncell,1-ncell:nft+ncell)
  integer ll(np_pp_max),idxf,ipf(np_pp_max),ift1(np_pp_max/8),ipp1,ipp2
  real xf(3,np_pp_max),vf(3,np_pp_max),af(3,np_pp_max)
  real ftemp

  integer icellpp,ncellpp,ijk(3,62)

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
    !if (np_pp_max<sum(rhoc(0:nt+1,0:nt+1,0:nt+1,itx,ity,itz))) then
    !  print*, 'np_pp_max too small'
    !  print*, np_pp_max, sum(rhoc(0:nt+1,0:nt+1,0:nt+1,itx,ity,itz))
    !  stop
    !endif
    ! the offset to use tile-based linked-list
    !print*,'  ip_offset',ip_offset
    !print*,'           ',idx_b_l(0,0,itx,ity,itz)+sum(rhoc(:-1,0,0,itx,ity,itz))
    !ip_offset=idx_b_l(0,0,itx,ity,itz)+sum(rhoc(:-1,0,0,itx,ity,itz))



    !print*,'    linked list'
    hoc=0; ll=0; idxf=0;
    do igz=0,nt+1
    do igy=0,nt+1
    do igx=0,nt+1
      np=rhoc(igx,igy,igz,itx,ity,itz)
      nzero=idx_b_r(igy,igz,itx,ity,itz)-sum(rhoc(igx:,igy,igz,itx,ity,itz))
      do lp=1,np ! loop over cdm particles
        ip1=nzero+lp
        xvec1=ncell*((/igx,igy,igz/)-1)+ncell*(int(xp(:,ip1)+ishift,izipx)+rshift)*x_resolution
        vreal=tan(pi*real(vp(:,ip1))/real(nvbin-1))/(sqrt(pi/2)/(sigma_vi*vrel_boost))
        !if (igx==1.and.igy==1.and.igz==1) then
        !  print*,'ip1',ipll1,ip_offset,ip1
        !  print*,'xf1',xvec1
        !  print*,'vf1',vreal
        !  print*,'vp1',vp(:,ip1)
          !stop
        !endif
        idxf=idxf+1 ! index in the xf
        xf(:,idxf)=xvec1 ! make temp list xf
        vf(:,idxf)=vreal ! make temp list vf
        ipf(idxf)=ip1 ! record index of ip in xf order, not consecutive due to buffer depth
        ivec1=floor(xvec1)+1
        !ipll1=ip1-ip_offset
        !ll(ipll1)=hoc(ivec1(1),ivec1(2),ivec1(3))
        !hoc(ivec1(1),ivec1(2),ivec1(3))=ipll1
        ll(idxf)=hoc(ivec1(1),ivec1(2),ivec1(3))
        hoc(ivec1(1),ivec1(2),ivec1(3))=idxf
        !ip1=ipll1+ip_offset
      enddo
    enddo
    enddo
    enddo
    !print*, 'ntemp',sum(rhoc(0:nt+1,0:nt+1,0:nt+1,itx,ity,itz)),idxf
    !print*, 'xf1',xf(:,:1)
    !print*, 'vf1',vf(:,:1)
    !print*, 'ipf1',ipf(1)
    !print*, hoc(1,1,1),ipf(hoc(1,1,1))
    !print*,xf(:,hoc(1,1,1))
    !print*,vf(:,hoc(1,1,1))
    !print*,ipf(1),ipf(idxf); stop


    !print*,'    force'
    af=0 ! account acceleration
    do ilayer=0,pp_range
    !$omp paralleldo &
    !$omp& default(shared) schedule(dynamic,1)&
    !$omp& private(igz,igy,igx,ipll1,it1,ift1,it2,icellpp,ipll2,ipp1,ipp2) &
    !$omp& private(xvec21,rmag,rcut,pcut,force_pp)
    do igz=1+ilayer,nft+pp_range,pp_range+1
    !do igz=1,nft+pp_range    do igy=0,nft+1
    do igy=1-pp_range,nft+pp_range
    do igx=0,nft+1
      ipll1=hoc(igx,igy,igz)
      it1=0 ! temp list in cell
      do while (ipll1/=0) ! central cell
        it1=it1+1
        ift1(it1)=ipll1 ! ipll1 is the index of xp, ift1 is a temp index list
        ipll1=ll(ipll1)
      enddo
      if (it1==0) cycle
      !print*, 'cell central=',igx,igy,igz
      !print*, 'it1=',it1
      it2=it1

      do icellpp=1,ncellpp ! temp list around cell
        ipll2=hoc(igx+ijk(1,icellpp),igy+ijk(2,icellpp),igz+ijk(3,icellpp))
        !print*,'considering cell',igx+ijk(1,icellpp),igy+ijk(2,icellpp),igz+ijk(3,icellpp)
        do while (ipll2/=0)
          it2=it2+1
          ift1(it2)=ipll2
          ipll2=ll(ipll2)
        enddo
      enddo
      if (it2==1) cycle

      !print*, 'cell surround=',ii,jj,kk
      !print*, 'it2=',it2
      !print*, xf(:,ift1(:it2))

      do ipp1=1,it1 ! update acceleration list
      do ipp2=ipp1+1,it2
        xvec21=xf(:,ift1(ipp2))-xf(:,ift1(ipp1))
        rmag=sqrt(sum(xvec21**2))
        rmag=merge(1d0,rmag,rmag==0)
        rcut=rmag/nf_cutoff
        pcut=1-(7./4*rcut**3)+(3./4*rcut**5)
        force_pp=sim%mass_p_cdm*(xvec21/rmag**3)*pcut
        force_pp=merge(force_pp,force_pp*0,rmag>rsoft)
        af(:,ift1(ipp1))=af(:,ift1(ipp1))+force_pp
        af(:,ift1(ipp2))=af(:,ift1(ipp2))-force_pp
      enddo
      enddo
    enddo
    enddo
    enddo
    !$omp endparalleldo
    enddo

    do igz=1,nft ! maxval on hoc(1:nft,1:nft,1:nft)
    do igy=1,nft
    do igx=1,nft
      ipll1=hoc(igx,igy,igz)
      f2_max_pp=max(f2_max_pp,sum(af(:,ipll1)**2))
      ipll1=ll(ipll1)
    enddo
    enddo
    enddo

    !print*, '    update vp'
    !$omp paralleldo &
    !$omp& default(shared) schedule(dynamic,64)&
    !$omp& private(ipp1,ip1,vreal)
    do ipp1=1,idxf
      ip1=ipf(ipp1)
      vreal=vf(:,ipp1)+af(:,ipp1)*a_mid*dt/6/pi
      vp(:,ip1)=nint(real(nvbin-1)*atan(sqrt(pi/2)/(sigma_vi*vrel_boost)*vreal)/pi,kind=izipv)
      !if (ip1==78049) then
      !  print*,'ipp1,ip1',ipp1,ip1
      !  print*,'vf1',vf(:,ipp1)
      !  print*,'af1',sum(af(:,ipp1)**2)
      !  print*,'vreal',vreal
      !  print*,vp(:,ip1)
      !  stop
      !endif
    enddo
    !$omp endparalleldo

    !print*,sum(af(:,5489)**2) ; stop
    !print*, size(af(:,:100)**2), size(sum(af(:,:100)**2,1)) ;stop
    !print*,f2_max_pp
    call system_clock(t2,t_rate)
    if (head) print*, '  tile:',int(itx,1),int(ity,1),int(itz,1),' elapsed time =',real(t2-t1)/t_rate,'secs'
    !! print*, '  itest1 =',itest1
  enddo
  enddo
  enddo !! itz
  !stop
  ! for reference, in pm.f90
  !dt_fine=sqrt( 1.0 / (sqrt(maxval(f2_max_fine))*a_mid*GG) )
  !dt_coarse=sqrt( real(ncell) / (sqrt(f2_max_coarse)*a_mid*GG) )
  !dt_pp=sqrt(0.1*rsoft) / max(sqrt(maxval(f2_max_pp))*a_mid*GG,1e-3)
  sim%dt_pp=1.0*sqrt(1.0) / (sqrt(f2_max_pp)*a_mid*GG)
  sim%dt_pp=max(sim%dt_pp,0.0003)
  sync all
  do i=1,nn**3
    sim%dt_pp=min(sim%dt_pp,sim[i]%dt_pp)
  enddo
  sync all

  if (head) then
    call system_clock(tt2,t_rate)
    print*,'  max of f_pp', sqrt(f2_max_pp)
    print*,'  dt_pp',sim%dt_pp
    print*,'  updated',itest1,'particles'
    print*,'  average pairs',npairs/real(sim%nplocal)
    print*,'  elapsed time =',real(tt2-tt1)/t_rate,'secs'
    print*,''
  endif
  sync all
!stop
endsubroutine


endmodule

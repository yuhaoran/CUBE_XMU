module halo_output
  implicit none
  type type_halo_catalog_header
    integer nhalo_tot,nhalo,ninfo
    real linking_parameter
  endtype
  type type_halo_catalog_array
    real hmass ! number of particles ! 1:1
    real,dimension(3) :: x,v,q,u,jl,je,jt,s ! pos & vel in E,L,S spaces ! 2:25
    real,dimension(3,3) :: tide,xx,qq,vv,uu ! 2nd order stat ! 26:70
    real lambda_qq(3),vec_qq(3,3) ! eigenvalue/vector of qq ! 71:82
    real lambda_tide(3),vec_tide(3,3) ! eigenvalue/vector of tide ! 83:94
    real spin_ratio_e,spin_ratio_l ! kinematic spin parameter 95:96
    real lam_e,lam_l ! 97:98 ! Bullock spin parameter
  endtype
endmodule

subroutine halofind_FoF
  use parameters
  use halo_output
  use variables, only: xp,vp,pid,rhoc,vfield
  use variables, only: i,j,k,l,itx,ity,itz,nlast,ip,np,vreal
  use variables, only: nzero,idx_b_r
  implicit none

  !! smaller nfof is memory-lite but time-consuming
  integer,parameter :: nfof=nf_global ! nfof is the resolution to do percolation
  integer,parameter :: ninfo=98 ! number of real numbers per halo in the halo catalog
  real,parameter :: b_link=0.2 ! linking length
  character(4) b_link_string
  real,parameter :: np_halo_min=100 ! minimum number of particles to be a halo

  type(type_halo_catalog_header) halo_header
  type(type_halo_catalog_array),allocatable :: hcat(:)
  
  integer(8) nlist,nplocal
  integer idx(3),iq1,iq2,iq3,jp,jq(3),is1,is2,pbc(3),np_central,np_candidate,n_friend
  integer np_iso,np_mem,np_head,ngroup,nhalo,nhalo_tot,ihalo
  real msim2phys,pos1(3),rp,rp2,rsq,xf_hoc(3),r_equiv
  
  integer(4),allocatable :: hoc(:,:,:),ll(:),ip_party(:),ip_friend(:),llgp(:),hcgp(:),ecgp(:),isort_mass(:)
  integer(4),allocatable :: iph_halo_all(:),iph_halo(:)
  real,allocatable :: xf(:,:),vf(:,:),xf_party(:,:),dx(:,:),dv(:,:),dq(:,:),du(:,:),pos_fof(:,:),mu_quj(:),mu_xvj(:)
  real,allocatable :: x_mean_all(:,:),x_mean(:,:),v_mean_all(:,:),v_mean(:,:)
  real,allocatable :: q_mean_all(:,:),q_mean(:,:),u_mean_all(:,:),u_mean(:,:)
  real,allocatable :: mu(:,:)

  real,parameter :: lsim2pc=1e6 * box/nf_global/h0
  real,parameter :: G_phys=0.0043 ! pc * M_solar^-1 * (km/s)^2
  real,parameter :: rho_crit=2.7756e11 ! M_solar Mpc^-3 h^-2
  real,parameter :: grid_per_p=(ncell/np_nc)**3 * merge(2,1,body_centered_cubic)

  write(b_link_string,'(f4.2)') b_link
  if (head) then
    print*,'Runtime FoF halofinder'
    print*,'  nc,nnt,nt,ng,nf='
    print*,'  ',int(nc,2),int(nnt,2),int(nt,2),int(ng,2),int(nf,2)
    print*,'  nfof=',int(nfof,2)
    print*,'  np_halo_min=',np_halo_min
    print*,'  output file =', output_name_halo('FoF_b'//b_link_string)
  endif

  nplocal=sim%nplocal;    print*, '  nplocal=',nplocal
  msim2phys=rho_crit*sim%omega_m*sim%box**3/sim%npglobal/sim%h0;  
  print*,'  msim2phys',msim2phys
  nlist=nplocal/8

  print*, '  convert zip format to float format'
  ! xf(1:3,1:np_local) is the position list for all particles
  ! this consumes a lot memory and needs optimization
  allocate(xf(3,nplocal),vf(3,nplocal),xf_party(3,nlist))
  nlast=0
  do itz=1,nnt
  do ity=1,nnt
  do itx=1,nnt
    do k=1,nt
    do j=1,nt
    do i=1,nt
    np=rhoc(i,j,k,itx,ity,itz)
    nzero=idx_b_r(j,k,itx,ity,itz)-sum(rhoc(i:,j,k,itx,ity,itz))
    do l=1,np
      ip=nzero+l
      pos1=nt*((/itx,ity,itz/)-1)+ ((/i,j,k/)-1) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
      xf(:,nlast+l)=pos1/real(nc) ! xf is in range [0,1]
      vreal=tan((pi*real(vp(:,ip)))/real(nvbin-1)) / (sqrt(pi/2)/(sim%sigma_vi*vrel_boost))
      vreal=vreal+vfield(:,i,j,k,itx,ity,itz)
      vf(:,nlast+l)=vreal
    enddo
    nlast=nlast+np
    enddo
    enddo
    enddo
  enddo
  enddo
  enddo
  ! check the range of data
  print*, '  minval xf',minval(xf,dim=2)
  print*, '  maxval xf',maxval(xf,dim=2)
  print*, '  minval vf',minval(vf,dim=2)
  print*, '  maxval vf',maxval(vf,dim=2)
  
  print*,'  create hoc ll'
  allocate(hoc(nfof,nfof,nfof),ll(nplocal),ip_party(nlist),ip_friend(nlist),llgp(nplocal),hcgp(nplocal),ecgp(nplocal))
  hoc=0; ll=0
  do ip=1,nplocal
    idx=ceiling(xf(:,ip)*nfof) ! index of the grid
    ll(ip)=hoc(idx(1),idx(2),idx(3))
    hoc(idx(1),idx(2),idx(3))=ip
    hcgp(ip)=ip ! initialize hcgp(ip)=ip for isolated particles
  enddo
  llgp=0; ecgp=0; ! initialize group link list

  print*,'  percolation'
  rp=b_link*(float(nplocal))**(-1./3.) ! distance threshold in unit of box size
  rp2=rp**2
  print*,'  rp2=',rp2

  ! loop over fof cells
  do iq3=1,nfof
  do iq2=1,nfof
  do iq1=1,nfof
    ! np_central: particle number in the central cell
    ! np_candidate: particle number in central and surrounding cells
    np_candidate=0; np_central=0;
    ! make the temporary interaction list
    ! xf_party is a temporary position list
    ! ip_party is a temporary linked-list
    ! pbc is {-1,0,1} for {left,current,right} periodic boundary condition (PBC)

    ! central cell
    ip=hoc(iq1,iq2,iq3)
    do while (ip/=0)
      np_central=np_central+1
      xf_party(:,np_central)=xf(:,ip) ! position list in the party
      ip_party(np_central)=ip ! particle indext list in the party
      ip=ll(ip) ! find next particle in the chain
    enddo
    if (np_central==0) cycle ! for empty cell, cycle to the next cell
    np_candidate=np_central

    ! 9 cells in -z
    jq(3)=iq3-1
    do is2=-1,1
    do is1=-1,1
      jq(:2)=[iq1+is1,iq2+is2]
      pbc=floor(real(jq-1)/nfof)
      jq=jq-pbc*nfof ! apply PBC
      ip=hoc(jq(1),jq(2),jq(3))
      do while (ip/=0)
        np_candidate=np_candidate+1
        xf_party(:,np_candidate)=xf(:,ip)+pbc ! positions can be <0 or >1
        ip_party(np_candidate)=ip ! continue writing the linked list
        ip=ll(ip) ! find next particle in the chain
      enddo
    enddo
    enddo

    ! 3 cells in -y
    jq(3)=iq3; jq(2)=iq2-1
    do is1=-1,1
      jq(1)=iq1+is1
      pbc=floor(real(jq-1)/nfof)
      jq=jq-pbc*nfof
      ip=hoc(jq(1),jq(2),jq(3))
      do while (ip/=0)
        np_candidate=np_candidate+1
        xf_party(:,np_candidate)=xf(:,ip)+pbc
        ip_party(np_candidate)=ip
        ip=ll(ip)
      enddo
    enddo

    ! 1 cell in -x
    jq=[iq1-1,iq2,iq3]
    pbc=floor(real(jq-1)/nfof)
    jq=jq-pbc*nfof
    ip=hoc(jq(1),jq(2),jq(3))
    do while (ip/=0)
      np_candidate=np_candidate+1
      xf_party(:,np_candidate)=xf(:,ip)+pbc
      ip_party(np_candidate)=ip
      ip=ll(ip)
    enddo

    if (np_candidate<=1) cycle ! for no surrounding particle case, cycle to the next cell
    if (np_candidate>nlist) stop 'np_candidate>nlist' ! too many candidates in the party

    ! now, loop over all pairs of particles, but one particle needs to be in the central cell
    do i=1,np_central ! loop over central cell
      ! n_friend is the number of pairs linked
      n_friend=0
      do j=i+1,np_candidate ! loop over central and surroundings, but loop only to the right
        rsq=sum((xf_party(:,i)-xf_party(:,j))**2)
        if (rsq<=rp2) then
          n_friend=n_friend+1 ! increment n_friend
          ip_friend(n_friend)=ip_party(j) ! ip_friend
        endif
      enddo
      ! for friends, merge the friendship of the two with other particles
      do j=1,n_friend
        call merge_chain(ip_party(i),ip_friend(j))
      enddo
    enddo

  enddo
  enddo
  enddo

  ! count isolated, member and leader particles
  np_iso=0; np_mem=0; np_head=0;
  do i=1,nplocal
    if (hcgp(i)==i) then
      np_iso=np_iso+1
    elseif (hcgp(i)==0) then
      np_mem=np_mem+1
    else
      np_head=np_head+1
    endif
  enddo
  ngroup=np_head ! there are ngroup groups with at least 2 particles
  print*,'  np_iso, np_mem, np_head',np_iso, np_mem, np_head
  print*,'  percentage',[np_iso, np_mem, np_head]/real(nplocal)*100,' %'

  print*,'  fof particles'
  ! loop over groups and select halos
  ! find the particle member of each group/halo
  ! group := at least 2 particles
  ! halo := groups with at least np_halo_min particles
  allocate(x_mean_all(3,nlist),v_mean_all(3,nlist),iph_halo_all(nlist),dx(3,nlist),dv(3,nlist))
  allocate(q_mean_all(3,nlist),u_mean_all(3,nlist),dq(3,nlist),du(3,nlist))
  allocate(mu_quj(nlist),mu_xvj(nlist))

  ! read phi from IC and compute Lagrangian properties
  ! skipped in runtime version

  nhalo=0
  do ip=1,nplocal
    if (hcgp(ip)/=ip .and. hcgp(ip)/=0) then
      np=0
      jp=hcgp(ip)
      xf_hoc=xf(:,jp)
      !qi_hoc=qgrid(pid(jp))
      do while (jp/=0)
        np=np+1
        dx(:,np)=modulo(xf(:,jp)-xf_hoc+0.5,1.)-0.5
        dv(:,np)=vf(:,jp)
        jp=llgp(jp)
      enddo
      if (np<np_halo_min) cycle
      nhalo=nhalo+1
      iph_halo_all(nhalo)=hcgp(ip)
      x_mean_all(:,nhalo)=modulo(xf_hoc+sum(dx(:,:np),2)/np,1.) ! center of mass, \in [0,1)
      v_mean_all(:,nhalo)=sum(dv(:,:np),2)/np ! mean velocity
      endif
  enddo
  print*,'  nplocal =',nplocal
  print*,'  ngroup  =',ngroup
  print*,'  nhalo   =',nhalo

  ! transfer data to smaller arrays
  allocate(x_mean(3,nhalo),v_mean(3,nhalo),q_mean(3,nhalo),u_mean(3,nhalo),iph_halo(nhalo))
  x_mean=x_mean_all(:,:nhalo)*nf_global ! convert into unit of fine cell
  v_mean=v_mean_all(:,:nhalo)
  !q_mean=q_mean_all(:,:nhalo)
  !u_mean=u_mean_all(:,:nhalo)
  iph_halo=iph_halo_all(:nhalo)
  deallocate(x_mean_all,v_mean_all,q_mean_all,u_mean_all,iph_halo_all)

  allocate(hcat(nhalo),pos_fof(3,nlist))
  do i=1,3 ! transfer data to hcat
    hcat%x(i)=x_mean(i,:)
    hcat%v(i)=v_mean(i,:)
    !hcat%q(i)=q_mean(i,:)
    !hcat%u(i)=u_mean(i,:)
  enddo
  deallocate(x_mean,v_mean,q_mean,u_mean)

  print*,'  calculate additional halo properties'
    do ihalo=1,nhalo !! omp parallel do
      jp=iph_halo(ihalo) ! start from ip-header
      np=0; !tide=0;
      do while (jp/=0) ! loop over halo members
        np=np+1
        dx(:,np)=modulo(xf(:,jp)-hcat(ihalo)%x+0.5,1.)-0.5 ! \in [-0.5,0.5) ! topological distince
        dv(:,np)=vf(:,jp)-hcat(ihalo)%v ! peculiar velocity wrt com
        jp=llgp(jp)
      enddo
      hcat(ihalo)%hmass=np
      ! angular momentum vector
      hcat(ihalo)%je(1)=sum(dx(2,:np)*dv(3,:np)-dx(3,:np)*dv(2,:np))
      hcat(ihalo)%je(2)=sum(dx(3,:np)*dv(1,:np)-dx(1,:np)*dv(3,:np))
      hcat(ihalo)%je(3)=sum(dx(1,:np)*dv(2,:np)-dx(2,:np)*dv(1,:np))
      ! spin parameter
      r_equiv=(hcat(ihalo)%hmass * grid_per_p * 3/800/pi)**(1./3.) * lsim2pc
      hcat(ihalo)%lam_e = norm2(hcat(ihalo)%je)/hcat(ihalo)%hmass * (sim%vsim2phys*lsim2pc*nf_global)
      hcat(ihalo)%lam_e = hcat(ihalo)%lam_e / sqrt(2*G_phys*hcat(ihalo)%hmass*msim2phys* r_equiv)

      hcat(ihalo)%je=hcat(ihalo)%je/norm2(hcat(ihalo)%je)
      hcat(ihalo)%q=0 ! below properties see offline version
      hcat(ihalo)%u=0
      hcat(ihalo)%s=0
      hcat(ihalo)%jl=0; hcat(ihalo)%jt=0;
      hcat(ihalo)%tide=0; hcat(ihalo)%xx=0; hcat(ihalo)%qq=0
      hcat(ihalo)%vv=0; hcat(ihalo)%uu=0
      hcat(ihalo)%lambda_qq=0; hcat(ihalo)%vec_qq=0
      hcat(ihalo)%lambda_tide=0; hcat(ihalo)%vec_tide=0
      hcat(ihalo)%spin_ratio_e=0; hcat(ihalo)%spin_ratio_l=0
      hcat(ihalo)%lam_l=0
    enddo !! do ihalo=1,nhalo

    ! sort halos according to mass
    if (nhalo>0) then
      print*,'fof_sort'
      allocate(isort_mass(nhalo))
      isort_mass=0
      call indexx(nhalo,hcat(1:nhalo)%hmass,isort_mass)
      hcat=hcat(isort_mass(nhalo:1:-1)) ! sort halo info
      iph_halo=iph_halo(isort_mass(nhalo:1:-1)) ! sort halo hoc
      print*,'  hmass top-3'
      print*, hcat(:3)%hmass
      deallocate(isort_mass)
    endif

    ! write halo catalog
    halo_header%nhalo_tot=nhalo
    halo_header%nhalo=nhalo
    halo_header%ninfo=ninfo
    halo_header%linking_parameter=b_link
    open(21,file=output_name_halo('FoFruntime_b'//b_link_string),status='replace',access='stream')
    write(21) halo_header,hcat
    do ihalo=1,nhalo
      jp=iph_halo(ihalo)
      do while (jp/=0)
        write(21) pid(jp) ! write PID
        jp=llgp(jp)
      enddo
      jp=iph_halo(ihalo)
      do while (jp/=0)
        write(21) xf(:,jp),vf(:,jp) ! write Eulerian position and velocity
        jp=llgp(jp)
      enddo
    enddo
    close(21)
    deallocate(iph_halo)
    deallocate(xf,vf,xf_party,hoc,ll,ip_party,ip_friend,llgp,hcgp,ecgp,pos_fof)
    deallocate(dx,dv,dq,du,mu_quj,mu_xvj)
    deallocate(hcat)

  contains

  subroutine merge_chain(i,j)
    ! llgp is a linked list: llgp(ip) means ip->llgp
    ! ip1->ip2->...->ipn
    ! ip1 is the head of chain (hoc)
    ! ipn is the end of chain (eoc)
    integer i,j,ihead,jhead,iend,jend,ipart
    jend=merge(j,ecgp(j),ecgp(j)==0)
    iend=merge(i,ecgp(i),ecgp(i)==0)
    if (iend==jend) return ! same chain
    ihead=max(hcgp(i),hcgp(iend))
    jhead=max(hcgp(j),hcgp(jend))
    ipart=jhead ! change eoc of the chain-j
    do while (ipart/=0)
      ecgp(ipart)=iend ! set chain-j's eoc to be iend
      ipart=llgp(ipart) ! next particle
    enddo
    llgp(jend)=ihead ! link j group to i group
    ecgp(i)=iend
    hcgp(iend)=jhead ! change hoc
    hcgp(jend)=0 ! set jend as a member
  endsubroutine

  subroutine indexx(N,ARRIN,INDX)
    implicit none
    integer N
    integer INDX(n),INDXT,IR
    real ARRIN(N),Q

    DO 11 J=1,N
      INDX(J)=J
    11    CONTINUE
    L=N/2+1
    IR=N
    10    CONTINUE
      IF(L.GT.1)THEN
        L=L-1
        INDXT=INDX(L)
        Q=ARRIN(INDXT)
      ELSE
        INDXT=INDX(IR)
        Q=ARRIN(INDXT)
        INDX(IR)=INDX(1)
        IR=IR-1
        IF(IR.EQ.1)THEN
          INDX(1)=INDXT
          RETURN
        ENDIF
      ENDIF
      I=L
      J=L+L
    20      IF(J.LE.IR)THEN
        IF(J.LT.IR)THEN
          IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
        ENDIF
        IF(Q.LT.ARRIN(INDX(J)))THEN
          INDX(I)=INDX(J)
          I=J
          J=J+J
        ELSE
          J=IR+1
        ENDIF
      GO TO 20
      ENDIF
      INDX(I)=INDXT
    GO TO 10
  endsubroutine indexx

  
endsubroutine

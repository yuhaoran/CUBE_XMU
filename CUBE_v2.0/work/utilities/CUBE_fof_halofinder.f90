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
  endtype
endmodule

program CUBE_FoF
  use parameters
  use halo_output
  implicit none

  !! smaller nfof is memory-lite but time-consuming
  integer,parameter :: nfof=nf_global ! nfof is the resolution to do percolation
  integer,parameter :: ninfo=94 ! number of real numbers per halo in the halo catalog
  real,parameter :: b_link=0.20 ! linking length
  real,parameter :: np_halo_min=100 ! minimum number of particles to be a halo

  type(type_halo_catalog_header) halo_header
  type(type_halo_catalog_array),allocatable :: hcat(:)

  integer np_iso,np_mem,np_head,cur_checkpoint
  integer i,j,k,l,itx,ity,itz,nlast,ip,jp,np,nplocal,nlist,idx(3),n_friend,ngroup,ihalo,nhalo
  integer iq1,iq2,iq3,np_candidate,np_central,pbc(3),jq(3),is1,is2,qi(3),qi_hoc(3)
  integer(4) t1,t2,tt1,tt2,t_rate
  real pos1(3),rp,rp2,rsq,xf_hoc(3),vreal(3),torque(3,3),tide(6),jeig(3),t_opt(3),a_opt

  integer(4),allocatable :: rhoc(:,:,:,:,:,:)
  real(4),allocatable :: vfield(:,:,:,:,:,:,:),phi(:,:,:)[:]
  integer(izipx),allocatable :: xp(:,:)
  integer(izipv),allocatable :: vp(:,:)
  integer(4),allocatable :: pid(:),iph_halo_all(:),iph_halo(:)

  integer(4),allocatable :: hoc(:,:,:),ll(:),ip_party(:),ip_friend(:),llgp(:),hcgp(:),ecgp(:),isort_mass(:)
  real,allocatable :: xf(:,:),vf(:,:),xf_party(:,:),dx(:,:),dv(:,:),dq(:,:),du(:,:),pos_fof(:,:)
  real,allocatable :: x_mean_all(:,:),x_mean(:,:),v_mean_all(:,:),v_mean(:,:)
  real,allocatable :: q_mean_all(:,:),q_mean(:,:),u_mean_all(:,:),u_mean(:,:)
  real,allocatable :: mu(:,:)

  real,dimension(3,3) :: qq_halo,tide_halo
  real,dimension(3,3,3) :: dd,qq

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call geometry
  if (head) then
    print*,'External CUBE_FoF halofinder'
    print*,'nc,nnt,nt,ng,nf='
    print*,int(nc,2),int(nnt,2),int(nt,2),int(ng,2),int(nf,2)
    print*,'nfof='
    print*,int(nfof,2)
    print*,'np_halo_min='
    print*,np_halo_min

    print*, 'checkpoints:'
    open(16,file='../main/z_checkpoint.txt',status='old')
    do i=1,nmax_redshift
      read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
      print*, z_checkpoint(i)
    enddo
    71 n_checkpoint=i-1
    close(16)
    print*,''
  endif
  sync all
  n_checkpoint=n_checkpoint[1]
  z_checkpoint(:)=z_checkpoint(:)[1]
  sync all

  do cur_checkpoint= 2,n_checkpoint
    sim%cur_checkpoint=cur_checkpoint
    print*, 'start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))
    print*, '  read checkpoint header',output_name('info')
    ! read CUBE checkpoint
    open(11,file=output_name('info'),status='old',action='read',access='stream')
    read(11) sim
    close(11)
    nplocal=sim%nplocal;    print*, '  nplocal=',nplocal
    nlist=nplocal/8

    ! allocate checkpoint arrays
    allocate(xp(3,nplocal),vp(3,nplocal),rhoc(nt,nt,nt,nnt,nnt,nnt),vfield(3,nt,nt,nt,nnt,nnt,nnt),pid(nplocal))
    allocate(xf(3,nplocal),vf(3,nplocal),xf_party(3,nlist)) ! traditional checkpoint arrays
    open(11,file=output_name('xp'),status='old',action='read',access='stream')
    read(11) xp
    close(11)
    open(11,file=output_name('vp'),status='old',action='read',access='stream')
    read(11) vp
    close(11)
    open(11,file=output_name('np'),status='old',action='read',access='stream')
    read(11) rhoc
    close(11)
    open(11,file=output_name('vc'),status='old',action='read',access='stream')
    read(11) vfield
    close(11)
    open(11,file=output_name('id'),status='old',action='read',access='stream')
    read(11) pid
    close(11)
    pid=pid-1

    print*, '  convert zip format to float format'
    ! xf(1:3,1:np_local) is the position list for all particles
    ! this consumes a lot memory and needs optimization
    nlast=0
    do itz=1,nnt
    do ity=1,nnt
    do itx=1,nnt
      do k=1,nt
      do j=1,nt
      do i=1,nt
        np=rhoc(i,j,k,itx,ity,itz)
        do l=1,np
          ip=nlast+l
          pos1=nt*((/itx,ity,itz/)-1)+ ((/i,j,k/)-1) + (int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution
          xf(:,ip)=pos1/real(nc) ! xf is in range [0,1]
          vreal=tan((pi*real(vp(:,ip)))/real(nvbin-1)) / (sqrt(pi/2)/(sim%sigma_vi*vrel_boost))
          vreal=vreal+vfield(:,i,j,k,itx,ity,itz)
          vf(:,ip)=vreal
        enddo
        nlast=nlast+np
      enddo
      enddo
      enddo
    enddo
    enddo
    enddo
    deallocate(xp,vp,rhoc,vfield)

    ! check the range of data
    print*, '  minval xf',minval(xf,dim=2)
    print*, '  maxval xf',maxval(xf,dim=2)
    print*, '  minval vf',minval(vf,dim=2)
    print*, '  maxval vf',maxval(vf,dim=2)

    call system_clock(tt1,t_rate)
    ! create head-of-chain and linked-list for all particles
    print*, 'create hoc ll'
    ! the following consume a lot memory, needs to do tile by tile
    allocate(hoc(nfof,nfof,nfof),ll(nplocal),ip_party(nlist),ip_friend(nlist),llgp(nplocal),hcgp(nplocal),ecgp(nplocal))
    hoc=0; ll=0
    do ip=1,nplocal
      idx=ceiling(xf(:,ip)*nfof) ! index of the grid
      ll(ip)=hoc(idx(1),idx(2),idx(3))
      hoc(idx(1),idx(2),idx(3))=ip
      hcgp(ip)=ip ! initialize hcgp(ip)=ip for isolated particles
    enddo
    llgp=0; ecgp=0; ! initialize group link list

    ! use percolation method to find FoF groups
    print*,'percolation'
    call system_clock(t1,t_rate)
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

    call system_clock(t2,t_rate)
    print*, '  percolation time =',real(t2-t1)/t_rate,'secs';



    print*,'fof_particles'
    ! loop over groups and select halos
    ! find the particle member of each group/halo
    ! group := at least 2 particles
    ! halo := groups with at least np_halo_min particles
    allocate(x_mean_all(3,nlist),v_mean_all(3,nlist),iph_halo_all(nlist),dx(3,nlist),dv(3,nlist))
    allocate(q_mean_all(3,nlist),u_mean_all(3,nlist),dq(3,nlist),du(3,nlist))

    sim%cur_checkpoint=1
    print*,'  read initial potential field',output_name('phi1')
    allocate(phi(0:nf+1,0:nf+1,0:nf+1)[*])
    open(11,file=output_name('phi1'),status='old',action='read',access='stream')
    read(11) phi(1:nf,1:nf,1:nf)
    close(11)
    sim%cur_checkpoint=cur_checkpoint
    print*, '  buffer phi'
    phi(0,:,:)=phi(nf,:,:)[image1d(inx,icy,icz)]
    phi(nf+1,:,:)=phi(1,:,:)[image1d(ipx,icy,icz)]
    sync all
    phi(:,0,:)=phi(:,nf,:)[image1d(icx,iny,icz)]
    phi(:,nf+1,:)=phi(:,1,:)[image1d(icx,ipy,icz)]
    sync all
    phi(:,:,0)=phi(:,:,nf)[image1d(icx,icy,inz)]
    phi(:,:,nf+1)=phi(:,:,1)[image1d(icx,icy,ipz)]

    nhalo=0;
    do ip=1,nplocal
      if (hcgp(ip)/=ip .and. hcgp(ip)/=0) then
        np=0
        jp=hcgp(ip) ! hoc of the group
        xf_hoc=xf(:,jp) ! position of hoc
        qi_hoc=qgrid(pid(jp)) ! q-index of hoc
        !print*,pid(jp),qi_hoc; stop
        do while (jp/=0) ! loop over group members
          np=np+1
          ! position wrt hoc ! \in [-0.5,0.5) ! topological distance
          dx(:,np)=modulo(xf(:,jp)-xf_hoc+0.5,1.)-0.5
          ! peculiar velocity
          dv(:,np)=vf(:,jp)
          ! initial position
          qi=qgrid(pid(jp)) ! Lagrangian grid number ! \in {1,2,...,nf_global}
          dq(:,np)=modulo(qi-qi_hoc+nf_global/2,nf_global)-nf_global/2 ! qpos wrt qhoc
          ! initial velocity
          !print*,'qi',qi
          du(1,np)=-phi(qi(1)+1,qi(2),qi(3))+phi(qi(1)-1,qi(2),qi(3))
          du(2,np)=-phi(qi(1),qi(2)+1,qi(3))+phi(qi(1),qi(2)-1,qi(3))
          du(3,np)=-phi(qi(1),qi(2),qi(3)+1)+phi(qi(1),qi(2),qi(3)-1)
          jp=llgp(jp) ! next particle in chain
        enddo
        if (np<np_halo_min) cycle ! skip the current do-loop
        nhalo=nhalo+1 ! found one halo
        iph_halo_all(nhalo)=hcgp(ip) ! hoc of the halo "ip-header"
        x_mean_all(:,nhalo)=modulo(xf_hoc+sum(dx(:,:np),2)/np,1.) ! center of mass, \in [0,1)
        v_mean_all(:,nhalo)=sum(dv(:,:np),2)/np ! mean velocity
        q_mean_all(:,nhalo)=modulo(qi_hoc+sum(dq(:,:np),2)/real(np),real(nf_global)) ! q center of mass
        u_mean_all(:,nhalo)=sum(du(:,:np),2)/np
        !if (nhalo==1) then
        !  print*,'q_mean_all',q_mean_all(:,nhalo)
        !endif
      endif
    enddo
    print*,'  nplocal =',nplocal
    print*,'  ngroup  =',ngroup
    print*,'  nhalo   =',nhalo
    print*,'  calculated mean of x,v,q,u'

    ! transfer data to smaller arrays
    allocate(x_mean(3,nhalo),v_mean(3,nhalo),q_mean(3,nhalo),u_mean(3,nhalo),iph_halo(nhalo))
    x_mean=x_mean_all(:,:nhalo)*nf_global ! convert into unit of fine cell
    v_mean=v_mean_all(:,:nhalo)
    q_mean=q_mean_all(:,:nhalo)
    u_mean=u_mean_all(:,:nhalo)
    iph_halo=iph_halo_all(:nhalo)
    deallocate(x_mean_all,v_mean_all,q_mean_all,u_mean_all,iph_halo_all)

    allocate(hcat(nhalo),pos_fof(3,nlist))
    do i=1,3 ! transfer data to hcat
      hcat%x(i)=x_mean(i,:)
      hcat%v(i)=v_mean(i,:)
      hcat%q(i)=q_mean(i,:)
      hcat%u(i)=u_mean(i,:)
    enddo
    deallocate(x_mean,v_mean,q_mean,u_mean)

    print*,'calculate additional halo properties'
    do ihalo=1,nhalo !! omp parallel do
      jp=iph_halo(ihalo) ! start from ip-header
      np=0; tide=0;
      do while (jp/=0) ! loop over halo members
        np=np+1

        dx(:,np)=modulo(xf(:,jp)-hcat(ihalo)%x+0.5,1.)-0.5 ! \in [-0.5,0.5) ! topological distince

        dv(:,np)=vf(:,jp)-hcat(ihalo)%v ! peculiar velocity wrt com

        qi=qgrid(pid(jp)) ! Lagrangian grid number ! \in {1,2,...,nf_global}
        dq(:,np)=modulo(qi-hcat(ihalo)%q+nf_global/2,nf_global*1.)-nf_global/2 ! in unit of fine grid

        du(1,np)=-phi(qi(1)+1,qi(2),qi(3))+phi(qi(1)-1,qi(2),qi(3))
        du(2,np)=-phi(qi(1),qi(2)+1,qi(3))+phi(qi(1),qi(2)-1,qi(3))
        du(3,np)=-phi(qi(1),qi(2),qi(3)+1)+phi(qi(1),qi(2),qi(3)-1)
        du(:,np)=du(:,np)-hcat(ihalo)%u

        tide(1)=tide(1)-phi(qi(1)-1,qi(2),qi(3))+2*phi(qi(1),qi(2),qi(3))-phi(qi(1)+1,qi(2),qi(3));
        tide(2)=tide(2)-phi(qi(1),qi(2)-1,qi(3))+2*phi(qi(1),qi(2),qi(3))-phi(qi(1),qi(2)+1,qi(3));
        tide(3)=tide(3)-phi(qi(1),qi(2),qi(3)-1)+2*phi(qi(1),qi(2),qi(3))-phi(qi(1),qi(2),qi(3)+1);
        tide(4)=tide(4)-(phi(qi(1)-1,qi(2)-1,qi(3))+phi(qi(1)+1,qi(2)+1,qi(3))&
                        -phi(qi(1)-1,qi(2)+1,qi(3))-phi(qi(1)+1,qi(2)-1,qi(3)))/4
        tide(5)=tide(5)-(phi(qi(1),qi(2)-1,qi(3)-1)+phi(qi(1),qi(2)+1,qi(3)+1)&
                        -phi(qi(1),qi(2)-1,qi(3)+1)-phi(qi(1),qi(2)+1,qi(3)-1))/4
        tide(6)=tide(6)-(phi(qi(1)-1,qi(2),qi(3)-1)+phi(qi(1)+1,qi(2),qi(3)+1)&
                        -phi(qi(1)+1,qi(2),qi(3)-1)-phi(qi(1)-1,qi(2),qi(3)+1))/4
        jp=llgp(jp)
      enddo
      hcat(ihalo)%tide(1,1)=tide(1)/np; hcat(ihalo)%tide(2,2)=tide(2)/np; hcat(ihalo)%tide(3,3)=tide(3)/np;
      hcat(ihalo)%tide(1,2)=tide(4)/np; hcat(ihalo)%tide(2,1)=tide(4)/np;
      hcat(ihalo)%tide(2,3)=tide(5)/np; hcat(ihalo)%tide(3,2)=tide(5)/np;
      hcat(ihalo)%tide(3,1)=tide(6)/np; hcat(ihalo)%tide(1,3)=tide(6)/np;
      ! halo mass in unit of particle number
      hcat(ihalo)%hmass=np
      ! angular momentum vector
      hcat(ihalo)%je(1)=sum(dx(2,:np)*dv(3,:np)-dx(3,:np)*dv(2,:np))
      hcat(ihalo)%je(2)=sum(dx(3,:np)*dv(1,:np)-dx(1,:np)*dv(3,:np))
      hcat(ihalo)%je(3)=sum(dx(1,:np)*dv(2,:np)-dx(2,:np)*dv(1,:np))
      hcat(ihalo)%je=hcat(ihalo)%je/norm2(hcat(ihalo)%je)
      ! initial angular momentum vector
      hcat(ihalo)%jl(1)=sum(dq(2,:np)*du(3,:np)-dq(3,:np)*du(2,:np))
      hcat(ihalo)%jl(2)=sum(dq(3,:np)*du(1,:np)-dq(1,:np)*du(3,:np))
      hcat(ihalo)%jl(3)=sum(dq(1,:np)*du(2,:np)-dq(2,:np)*du(1,:np))
      hcat(ihalo)%jl=hcat(ihalo)%jl/norm2(hcat(ihalo)%jl)
      ! moment of inertia tensor
      hcat(ihalo)%xx(1,1)=sum(dx(1,:np)**2)
      hcat(ihalo)%xx(2,2)=sum(dx(2,:np)**2)
      hcat(ihalo)%xx(3,3)=sum(dx(3,:np)**2)
      hcat(ihalo)%xx(1,2)=sum(dx(1,:np)*dx(2,:np))
      hcat(ihalo)%xx(2,3)=sum(dx(2,:np)*dx(3,:np))
      hcat(ihalo)%xx(3,1)=sum(dx(3,:np)*dx(1,:np))
      ! Lagrangian inertia tensor
      hcat(ihalo)%qq(1,1)=sum(dq(1,:np)**2)/np
      hcat(ihalo)%qq(2,2)=sum(dq(2,:np)**2)/np
      hcat(ihalo)%qq(3,3)=sum(dq(3,:np)**2)/np
      hcat(ihalo)%qq(1,2)=sum(dq(1,:np)*dq(2,:np))/np; hcat(ihalo)%qq(2,1)=hcat(ihalo)%qq(1,2)
      hcat(ihalo)%qq(2,3)=sum(dq(2,:np)*dq(3,:np))/np; hcat(ihalo)%qq(3,2)=hcat(ihalo)%qq(2,3)
      hcat(ihalo)%qq(3,1)=sum(dq(3,:np)*dq(1,:np))/np; hcat(ihalo)%qq(1,3)=hcat(ihalo)%qq(3,1)
      ! TTT
      torque=matmul(hcat(ihalo)%qq,hcat(ihalo)%tide)
      hcat(ihalo)%jt(1)=torque(2,3)-torque(3,2)
      hcat(ihalo)%jt(2)=torque(3,1)-torque(1,3)
      hcat(ihalo)%jt(3)=torque(1,2)-torque(2,1)
      hcat(ihalo)%jt=hcat(ihalo)%jt/norm2(hcat(ihalo)%jt)
      ! eigen decomposition
      !print*,'eigenvalue'
      !print*,'hcat(ihalo)%hmass=',hcat(ihalo)%hmass
      !print*,'hcat(ihalo)%qq'
      !print*,hcat(ihalo)%qq
      
      qq_halo=hcat(ihalo)%qq
      call dsyevj3(qq_halo,hcat(ihalo)%lambda_qq,hcat(ihalo)%vec_qq)
      qq_halo=hcat(ihalo)%qq/sum(hcat(ihalo)%lambda_qq) ! rescale such that sum(lambda_qq)=1
      hcat(ihalo)%tide=hcat(ihalo)%tide*sum(hcat(ihalo)%lambda_qq)/norm2(hcat(ihalo)%jt) ! rescale such that |jt|=1
      call dsyevj3(qq_halo,hcat(ihalo)%lambda_qq,hcat(ihalo)%vec_qq)
      call eigsrt(hcat(ihalo)%lambda_qq,hcat(ihalo)%vec_qq)
      !print*,'hcat(ihalo)%lambda_qq'
      !print*,hcat(ihalo)%lambda_qq
      !print*,'hcat(ihalo)%vec_qq'
      !print*,hcat(ihalo)%vec_qq(:,1)
      !print*,hcat(ihalo)%vec_qq(:,2)
      !print*,hcat(ihalo)%vec_qq(:,3)

            ! to compare with qq, we use negative tide tensor
      tide_halo=-hcat(ihalo)%tide
      call dsyevj3(tide_halo,hcat(ihalo)%lambda_tide,hcat(ihalo)%vec_tide)
      call eigsrt(hcat(ihalo)%lambda_tide,hcat(ihalo)%vec_tide)
      !print*,'hcat(ihalo)%lambda_tide'
      !print*,hcat(ihalo)%lambda_tide
      !print*,'hcat(ihalo)%vec_tide'
      !print*,hcat(ihalo)%vec_tide(:,1)
      !print*,hcat(ihalo)%vec_tide(:,2)
      !print*,hcat(ihalo)%vec_tide(:,3)
      !stop

      ! vv
      hcat(ihalo)%vv=0
      ! uu
      hcat(ihalo)%uu=0
      ! redshift space center of mass
      hcat(ihalo)%s=0
    enddo !! do ihalo=1,nhalo

    ! symmetrize 3x3 matrices
    do j=1,3
      i=modulo(j,3)+1
      hcat%xx(j,i)=hcat%xx(i,j)
      hcat%vv(j,i)=hcat%vv(i,j)
      hcat%qq(j,i)=hcat%qq(i,j)
      hcat%uu(j,i)=hcat%uu(i,j)
    enddo
    deallocate(phi)

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

    !! write halo catalog consistent with CUBE_SO halofiner
    halo_header%nhalo_tot=nhalo
    halo_header%nhalo=nhalo
    halo_header%ninfo=ninfo
    halo_header%linking_parameter=b_link
    print*,'  write', output_name('fof')
    open(11,file=output_name('fof'),status='replace',access='stream')
    write(11) halo_header,hcat
    close(11)

    !! write halo PIDs into file
    print*,'  write', output_name('fofpid')
    open(11,file=output_name('fofpid'),status='replace',access='stream')
    write(11) nhalo ! write nhalo as header
    do ihalo=1,nhalo
      jp=iph_halo(ihalo)
      write(11) nint(hcat(ihalo)%hmass) ! write particle number of halo
      do while (jp/=0)
        write(11) pid(jp) ! write PID
        jp=llgp(jp)
      enddo
      flush(11)
    enddo
    close(11)

    deallocate(iph_halo)
    deallocate(xf,vf,xf_party,hoc,ll,ip_party,ip_friend,llgp,hcgp,ecgp,pos_fof,pid)
    deallocate(dx,dv,dq,du)
    deallocate(hcat)
    
    call system_clock(tt2,t_rate)
    print*, 'total time =',real(tt2-tt1)/t_rate,'secs';
  enddo
  print*,'CUBE_FoF done'

  
  stop
  print*,''
  print*,'spin correlation analysis'
  allocate(mu(6,nhalo))
  do ihalo=6,nhalo
    !! eigen-optimization
    !print*,'lambda qq',hcat(ihalo)%lambda_qq
    t_opt=[1.0, 1.0, 1.0] ! weights before normalization
    a_opt=sum(hcat(ihalo)%lambda_qq)/sum(t_opt*hcat(ihalo)%lambda_qq)/norm2(hcat(ihalo)%jt) ! normalization factor
    !print*,'norm2(hcat(ihalo)%jt)',norm2(hcat(ihalo)%jt)
    !print*,'a_opt',a_opt

    hcat(ihalo)%jt=hcat(ihalo)%jt/norm2(hcat(ihalo)%jt)
    mu(1,ihalo)=sum( hcat(ihalo)%jl * hcat(ihalo)%je )
    mu(2,ihalo)=sum( hcat(ihalo)%jl * hcat(ihalo)%jt )
    mu(3,ihalo)=sum( hcat(ihalo)%jt * hcat(ihalo)%je )
    dd=0; qq=0;

    !print*,mu(3,ihalo)
    
    do i=1,3
      dd(i,i,i)=hcat(ihalo)%lambda_qq(i)*t_opt(i)!*a_opt
      qq(:,:,i)=matmul(matmul(hcat(ihalo)%vec_qq,dd(:,:,i)),transpose(hcat(ihalo)%vec_qq))
      torque=matmul(qq(:,:,i),hcat(ihalo)%tide)
      !torque=matmul(sum(qq,3),hcat(ihalo)%tide)
      jeig(1)=torque(2,3)-torque(3,2)
      jeig(2)=torque(3,1)-torque(1,3)
      jeig(3)=torque(1,2)-torque(2,1)
      !print*,'norm2 jeig',norm2(jeig)
      !print*,'norm',norm2(hcat(ihalo)%je),norm2(jeig)
      mu(3+i,ihalo)=sum(hcat(ihalo)%je*jeig)/norm2(jeig)
      !print*,'mu',mu(3+i,ihalo)
    enddo
    !print*,'sum of qq' ! should equal to hcat(ihalo)%qq
    !print*,sum(qq,3)
    !stop

  enddo
  print*,'mu_le',sum(mu(1,:))/nhalo
  print*,'mu_lt',sum(mu(2,:))/nhalo
  print*,'mu_te',sum(mu(3,:))/nhalo
  !print*,'mu(I1T,e)',sum(mu(4,:))/nhalo
  !print*,'mu(I2T,e)',sum(mu(5,:))/nhalo
  !print*,'mu(I3T,e)',sum(mu(6,:))/nhalo
  deallocate(mu)

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

    function qgrid(pid0)
      implicit none
      integer pid0,nfg2,qgrid(3)
      nfg2=nf_global**2
      qgrid(3)=1+pid0/nfg2
      qgrid(2)=1+(pid0-(pid0/nfg2)*nfg2)/nf_global
      qgrid(1)=1+modulo(pid0,nf_global)
    endfunction

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

    SUBROUTINE dsyevj3(A,W,Q)
    REAL A(3,3)
    REAL Q(3,3)
    REAL W(3)
    INTEGER          N
    PARAMETER        ( N = 3 )
    REAL SD, SO
    REAL S, C, T
    REAL G, H, Z, THETA
    REAL THRESH
    INTEGER I, X, Y, R

    DO 10 X = 1, N
      Q(X,X) = 1.0D0
      DO 11, Y = 1, X-1
        Q(X, Y) = 0.0D0
        Q(Y, X) = 0.0D0
 11   CONTINUE
 10 CONTINUE

    DO 20 X = 1, N
      W(X) = A(X, X)
 20 CONTINUE

    SD = 0.0D0
    DO 30 X = 1, N
      SD = SD + ABS(W(X))
 30 CONTINUE
    SD = SD**2

    DO 40 I = 1, 50
      SO = 0.0D0
      DO 50 X = 1, N
        DO 51 Y = X+1, N
          SO = SO + ABS(A(X, Y))
 51     CONTINUE
 50   CONTINUE
      IF (SO .EQ. 0.0D0) THEN
        RETURN
      END IF

      IF (I .LT. 4) THEN
        THRESH = 0.2D0 * SO / N**2
      ELSE
        THRESH = 0.0D0
      END IF

      DO 60 X = 1, N
        DO 61 Y = X+1, N
          G = 100.0D0 * ( ABS(A(X, Y)) )
          IF ( I .GT. 4 .AND. ABS(W(X)) + G .EQ. ABS(W(X)) &
                       .AND. ABS(W(Y)) + G .EQ. ABS(W(Y)) ) THEN
            A(X, Y) = 0.0D0
          ELSE IF (ABS(A(X, Y)) .GT. THRESH) THEN
            H = W(Y) - W(X)
            IF ( ABS(H) + G .EQ. ABS(H) ) THEN
              T = A(X, Y) / H
            ELSE
              THETA = 0.5D0 * H / A(X, Y)
              IF (THETA .LT. 0.0D0) THEN
                T = -1.0D0 / (SQRT(1.0D0 + THETA**2) - THETA)
              ELSE
                T = 1.0D0 / (SQRT(1.0D0 + THETA**2) + THETA)
              END IF
            END IF

            C = 1.0D0 / SQRT( 1.0D0 + T**2 )
            S = T * C
            Z = T * A(X, Y)

            A(X, Y) = 0.0D0
            W(X)    = W(X) - Z
            W(Y)    = W(Y) + Z
            DO 70 R = 1, X-1
              T       = A(R, X)
              A(R, X) = C * T - S * A(R, Y)
              A(R, Y) = S * T + C * A(R, Y)
 70         CONTINUE
            DO 80, R = X+1, Y-1
              T       = A(X, R)
              A(X, R) = C * T - S * A(R, Y)
              A(R, Y) = S * T + C * A(R, Y)
 80         CONTINUE
            DO 90, R = Y+1, N
              T       = A(X, R)
              A(X, R) = C * T - S * A(Y, R)
              A(Y, R) = S * T + C * A(Y, R)
 90         CONTINUE

            DO 100, R = 1, N
              T       = Q(R, X)
              Q(R, X) = C * T - S * Q(R, Y)
              Q(R, Y) = S * T + C * Q(R, Y)
100         CONTINUE
          END IF
 61     CONTINUE
 60   CONTINUE
 40 CONTINUE
    PRINT *, "DSYEVJ3: No convergence."
    END SUBROUTINE

    SUBROUTINE eigsrt(d,v)
    INTEGER,PARAMETER:: n=3
    REAL d(n),v(n,n)
    INTEGER i,j,k
    REAL p
    do 13 i=1,n-1
      k=i
      p=d(i)
      do 11 j=i+1,n
        if(d(j).ge.p)then
          k=j
          p=d(j)
        endif
11      continue
      if(k.ne.i)then
        d(k)=d(i)
        d(i)=p
        do 12 j=1,n
          p=v(j,i)
          v(j,i)=v(j,k)
          v(j,k)=p
12        continue
      endif
13    continue
    return
    END SUBROUTINE

    pure function matinv3(A) result(B)
      !! Performs a direct calculation of the inverse of a 3×3 matrix.
      real(4), intent(in) :: A(3,3)   !! Matrix
      real(4)             :: B(3,3)   !! Inverse matrix
      real(4)             :: detinv
      ! Calculate the inverse determinant of the matrix
      detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
                - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
                + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))
      ! Calculate the inverse of the matrix
      B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
      B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
      B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
      B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
      B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
      B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
      B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
      B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
      B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    end function

end

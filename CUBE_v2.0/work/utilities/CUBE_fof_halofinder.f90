module halo_output
  implicit none
  type type_halo_catalog_header
    integer nhalo_tot,nhalo
    real linking_parameter
  endtype
  type type_halo_cat_array
    real hpos(3)
    real mass_fof,radius_fof,v_disp
    real x_mean(3),v_mean(3),ang_mom(3),var_x(3),inertia(3,3)
    real q_mean(3),inertia_q(3,3),s_mean(3)
  endtype
endmodule

program CUBE_FoF
  !! written by Haoran Yu. haoran@cita.utoronto.ca
  !! test program for FoF halofinder
  !! suitable for CUBE's zip checkpoints
  !! optimized percolation algorithm
  !! todo list: OpenMP and MPI parallelization, memory optimization
  !! acknowledgement: Prof. Yipeng Jing's data structure for percolation
  use parameters
  use halo_output
  implicit none

  !! smaller nfof is memory-lite but time-consuming
  integer,parameter :: nfof=nc*4 ! nfof is the resolution to do percolation
  real,parameter :: b_link=0.2 ! linking length
  real,parameter :: np_halo_min=30 ! minimum number of particles to be a halo

  type(type_halo_catalog_header) halo_header !! halo catalog header
  type(type_halo_cat_array),allocatable :: hcat(:) !! halo catalog array

  integer np_iso,np_mem,np_head,cur_checkpoint
  integer i,j,k,l,itx,ity,itz,nlast,ip,jp,np,nplocal,nlist,idx(3),n_friend,ngroup,ihalo,nhalo
  integer iq1,iq2,iq3,np_candidate,np_central,pbc(3),jq(3),is1,is2,is3
  integer(4) t1,t2,tt1,tt2,ttt1,ttt2,t_rate
  !integer(4) rhoc(nt,nt,nt,nnt,nnt,nnt)

  real pos1(3),rp,rp2,rsq,xf_hoc(3),vreal(3)

  integer(4),allocatable :: rhoc(:,:,:,:,:,:)
  real(4),allocatable :: vfield(:,:,:,:,:,:,:)
  integer(izipx),allocatable :: xp(:,:)
  integer(izipv),allocatable :: vp(:,:)
  integer(4),allocatable :: pid(:),iph_halo_all(:),iph_halo(:)

  integer(4),allocatable :: hoc(:,:,:),ll(:),ip_party(:),ip_friend(:),llgp(:),hcgp(:),ecgp(:),isort_mass(:),iwsp(:)
  real,allocatable :: xf(:,:),vf(:,:),xf_party(:,:),dx(:,:),dv(:,:),vi(:,:),pos_fof(:,:),hm(:)
  real,allocatable :: x_mean_all(:,:),x_mean(:,:),v_mean_all(:,:),v_mean(:,:)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call geometry
  if (head) then
    print*, 'External CUBE_FoF halofinder'
    print*, 'nc,nnt,nt,ng,nf='
    print*, int(nc,2),int(nnt,2),int(nt,2),int(ng,2),int(nf,2)
    print*, 'nfof='
    print*, int(nfof,2)

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

  do cur_checkpoint= n_checkpoint,n_checkpoint
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
    nhalo=0;
    do ip=1,nplocal
      if (hcgp(ip)/=ip .and. hcgp(ip)/=0) then
        np=0
        jp=hcgp(ip) ! hoc of the group
        xf_hoc=xf(:,jp) ! position of hoc
        do while (jp/=0) ! loop over group members
          np=np+1
          dx(:,np)=modulo(xf(:,jp)-xf_hoc+0.5,1.)-0.5 ! position wrt hoc ! \in [-0.5,0.5) ! topological distance
          dv(:,np)=vf(:,jp) ! peculiar velocity
          jp=llgp(jp) ! next particle in chain
        enddo
        if (np<np_halo_min) cycle ! skip the current do-loop
        !print*,np
        nhalo=nhalo+1 ! found one halo
        !print*,nhalo
        iph_halo_all(nhalo)=hcgp(ip) ! hoc of the halo "ip-header"
        !print*,ip
        x_mean_all(:,nhalo)=modulo(xf_hoc+sum(dx(:,:np),2)/np,1.) ! center of mass, \in [0,1)
        v_mean_all(:,nhalo)=sum(dv(:,:np),2)/np ! mean velocity
      endif
    enddo
    print*,'  nplocal =',nplocal
    print*,'  ngroup  =',ngroup
    print*,'  nhalo   =',nhalo

    ! transfer data to smaller arrays
    allocate(x_mean(3,nhalo),v_mean(3,nhalo),iph_halo(nhalo))
    x_mean=x_mean_all(:,:nhalo)
    v_mean=v_mean_all(:,:nhalo)
    iph_halo=iph_halo_all(:nhalo)
    deallocate(x_mean_all,v_mean_all,iph_halo_all)

    allocate(hcat(nhalo),pos_fof(3,nlist))
    do i=1,3 ! transfer data to hcat
      hcat%hpos(i)=x_mean(i,:); hcat%x_mean(i)=x_mean(i,:)
      hcat%v_mean(i)=v_mean(i,:)
    enddo
    deallocate(x_mean,v_mean)

    print*,'calculate additional halo properties'
    do ihalo=1,nhalo !! omp parallel do
      jp=iph_halo(ihalo) ! start from ip-header
      np=0
      do while (jp/=0) ! loop over halo members
        np=np+1
        dx(:,np)=modulo(xf(:,jp)-hcat(ihalo)%hpos+0.5,1.)-0.5 ! \in [-0.5,0.5) ! topological distince
        dv(:,np)=vf(:,jp)-hcat(ihalo)%v_mean ! peculiar velocity
        jp=llgp(jp)
      enddo
      ! halo mass in unit of particle number
      hcat(ihalo)%mass_fof=np
      ! angular momentum vector
      hcat(ihalo)%ang_mom(1)=sum(dx(2,:np)*dv(3,:np)-dx(3,:np)*dv(2,:np))
      hcat(ihalo)%ang_mom(2)=sum(dx(3,:np)*dv(1,:np)-dx(1,:np)*dv(3,:np))
      hcat(ihalo)%ang_mom(3)=sum(dx(1,:np)*dv(2,:np)-dx(2,:np)*dv(1,:np))
      ! moment of inertia tensor
      hcat(ihalo)%inertia(1,1)=sum(dx(1,:np)**2)
      hcat(ihalo)%inertia(2,2)=sum(dx(2,:np)**2)
      hcat(ihalo)%inertia(3,3)=sum(dx(3,:np)**2)
      hcat(ihalo)%inertia(1,2)=sum(dx(1,:np)*dx(2,:np))
      hcat(ihalo)%inertia(2,3)=sum(dx(2,:np)*dx(3,:np))
      hcat(ihalo)%inertia(3,1)=sum(dx(3,:np)*dx(1,:np))
      ! radius
      hcat(ihalo)%radius_fof=0
      ! velocity dispersion
      hcat(ihalo)%v_disp=0
      ! position variance
      hcat(ihalo)%var_x=0
      ! Lagrangian center of mass
      hcat(ihalo)%q_mean=0
      ! Lagrangian moment intertia tensor
      hcat(ihalo)%inertia_q=0
      ! redshift space center of mass
      hcat(ihalo)%s_mean=0
    enddo
    ! symmetrize tensors
    hcat%inertia(2,1)=hcat%inertia(1,2)
    hcat%inertia(3,2)=hcat%inertia(2,3)
    hcat%inertia(1,3)=hcat%inertia(3,1)
    hcat%inertia_q(2,1)=hcat%inertia_q(1,2)
    hcat%inertia_q(3,2)=hcat%inertia_q(2,3)
    hcat%inertia_q(1,3)=hcat%inertia_q(3,1)

    if (nhalo>0) then
      print*,'fof_sort'
      allocate(isort_mass(nhalo))
      isort_mass=0
      call indexx(nhalo,hcat(1:nhalo)%mass_fof,isort_mass)
      hcat=hcat(isort_mass(nhalo:1:-1)) ! sort halo info
      iph_halo=iph_halo(isort_mass(nhalo:1:-1)) ! sort halo hoc
      print*,'  mass_fof top-10'
      print*, hcat(:10)%mass_fof
      deallocate(isort_mass)
    endif

    !! write halo catalog consistent with CUBE_SO halofiner
    halo_header%nhalo_tot=nhalo
    halo_header%nhalo=nhalo
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
      write(11) nint(hcat(ihalo)%mass_fof) ! write particle number of halo
      do while (jp/=0)
        write(11) pid(jp) ! write PID
        jp=llgp(jp)
      enddo
      flush(11)
    enddo
    close(11)

    deallocate(iph_halo)
    deallocate(xf,xf_party,hoc,ll,ip_party,ip_friend,llgp,hcgp,ecgp,dx,pos_fof,pid,hcat)

    call system_clock(tt2,t_rate)
    print*, 'total time =',real(tt2-tt1)/t_rate,'secs';
    print*, ''
  enddo
  print*,'CUBE_FoF done'

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

end

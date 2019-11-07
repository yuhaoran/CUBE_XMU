!#define READ_PLIST

program fof
  use parameters
  implicit none

  integer,parameter :: nfof=nft*2 ! nfof is the resolution to do percolation
  real,parameter :: b_link=0.2
  real,parameter :: massll=2

  integer nptest,ipart_iso,ipart_mem,ipart_head
  integer i,j,k,l,itx,ity,itz,nlast,ip,np,npart,idx(3),ntp,ngr,ngr0,grm
  integer iq1,iq2,iq3,nt1,nt2,peri(3),jq1,jq2,jq3,is1,is2,is3
  integer(4) t1,t2,tt1,tt2,ttt1,ttt2,t_rate
  integer(4) rhoc(nt,nt,nt,nnt,nnt,nnt)

  real pos1(3),rp,rp2,rsq,x0(3),x1(3)

  integer(izipx),allocatable :: xp(:,:)
  integer(4),allocatable :: hoc(:,:,:),ll(:),ll1(:),ll2(:),llgp(:),hcgp(:),ecgp(:),nab(:),iwsp(:)
  real,allocatable :: xf(:,:),xf1(:,:),xi(:,:),xgr(:,:),cm(:),wsp(:)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call geometry
  if (head) then
    print*, 'fof'
    print*, 'nc,nnt,nt,ng='
    print*, int(nc,2),int(nnt,2),int(nt,2),int(ng,2)
    print*, 'nfof='
    print*, int(nfof,2)
  endif

  if (head) then
    print*, 'checkpoint at:'
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

  do cur_checkpoint= 1,n_checkpoint
    print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))
    print*, 'read checkpoint header',output_name('info')
    open(11,file=output_name('info'),status='old',action='read',access='stream')
    read(11) sim
    close(11)
    npart=sim%nplocal
    print*, '  npart=',npart
    allocate(xp(3,npart),xf(3,npart),xf1(3,npart/8))
    open(11,file=output_name('xp'),status='old',action='read',access='stream')
    read(11) xp
    close(11)
    open(11,file=output_name('np'),status='old',action='read',access='stream')
    read(11) rhoc
    close(11)

    print*, 'convert zip format to float format'
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
          xf(:,ip)=pos1/real(nc)
        enddo
        nlast=nlast+np
      enddo
      enddo
      enddo
    enddo
    enddo
    enddo
    deallocate(xp)

    print*, 'write float format into file'
    open(11,file=output_name('xf'),status='replace',access='stream')
    write(11) xf
    close(11)

#   ifdef READ_PLIST
      npart=2*128**3
      allocate(xp(3,npart),xf(3,npart),xf1(3,npart/8))
      print*, 'read float format from file'
      open(11,file=output_name('xf_pp2'),status='old',access='stream')
      read (11) xf
      close(11)
#   endif

    ! check the range of data
    print*, '  minval',minval(xf,dim=2)
    print*, '  maxval',maxval(xf,dim=2)

    call system_clock(tt1,t_rate)
    print*, 'create hoc ll'
    allocate(hoc(nfof,nfof,nfof),ll(npart),ll1(npart),ll2(npart),llgp(npart),hcgp(npart),ecgp(npart))
    hoc=0; ll=0
    do ip=1,npart
      idx=ceiling(xf(:,ip)*nfof)
      ll(ip)=hoc(idx(1),idx(2),idx(3))
      hoc(idx(1),idx(2),idx(3))=ip
    enddo

    print*,'percolation'
    call system_clock(t1,t_rate)
    rp=b_link*(float(npart))**(-1./3.)
    rp2=rp**2
    print*,'  rp2=',rp2
    llgp=0; ecgp=0
    do ip=1,npart
      hcgp(ip)=ip ! hcgp(ip)=ip for isolated particles
    enddo
    !nptest=0
    do iq3=1,nfof
    do iq2=1,nfof
    do iq1=1,nfof
      nt1=0; nt2=0
      ! make the temporary interaction list

      ! central cell
      ip=hoc(iq1,iq2,iq3)
      do while (ip/=0)
        nt1=nt1+1
        xf1(:,nt1)=xf(:,ip) ! temp xf
        ll1(nt1)=ip ! temp ll
        ip=ll(ip)
      enddo
      if (nt1==0) cycle
      nt2=nt1 ! np_central

      ! 9 cells in -z
      jq3=iq3-1
      do is2=-1,1
      do is1=-1,1
        jq1=iq1+is1
        jq2=iq2+is2
        !jq1=0;jq2=1;jq3=nfof+1
        peri(:)=floor(real([jq1,jq2,jq3]-1)/nfof)
        jq1=jq1-peri(1)*nfof
        jq2=jq2-peri(2)*nfof
        jq3=jq3-peri(3)*nfof
        !print*,jq1,jq2,jq3; stop
        ip=hoc(jq1,jq2,jq3)
        do while (ip/=0)
          nt1=nt1+1
          xf1(:,nt1)=xf(:,ip)+peri
          ll1(nt1)=ip
          ip=ll(ip)
        enddo
      enddo
      enddo

      ! 3 cells in -y
      jq3=iq3
      jq2=iq2-1
      do is1=-1,1
        jq1=iq1+is1
        peri(:)=floor(real([jq1,jq2,jq3]-1)/nfof)
        jq1=jq1-peri(1)*nfof
        jq2=jq2-peri(2)*nfof
        jq3=jq3-peri(3)*nfof
        ip=hoc(jq1,jq2,jq3)
        do while (ip/=0)
          nt1=nt1+1
          xf1(:,nt1)=xf(:,ip)+peri
          ll1(nt1)=ip
          ip=ll(ip)
        enddo
      enddo

      ! 1 cell in -x
      jq1=iq1-1
      jq2=iq2
      jq3=iq3
      peri(:)=floor(real([jq1,jq2,jq3]-1)/nfof)
      jq1=jq1-peri(1)*nfof
      jq2=jq2-peri(2)*nfof
      jq3=jq3-peri(3)*nfof
      ip=hoc(jq1,jq2,jq3)
      do while (ip/=0)
        nt1=nt1+1
        xf1(:,nt1)=xf(:,ip)+peri
        ll1(nt1)=ip
        ip=ll(ip)
      enddo

      if (nt1<=1) cycle
      if (nt1>npart/8) stop 'nt1>npart/8'

      !ll2=0
      do i=1,nt2 ! central cell
        ntp=0
        do j=i+1,nt1
          rsq=sum((xf1(:,i)-xf1(:,j))**2)
          if (rsq<=rp2) then
            ntp=ntp+1
            ll2(ntp)=ll1(j)
          endif
        enddo

        do j=1,ntp
          call merge_group(hcgp,llgp,ecgp,npart,ll1(i),ll2(j))
        enddo
      enddo

    enddo
    enddo
    enddo

    ipart_iso=0; ipart_mem=0; ipart_head=0;
    do i=1,npart
      if (hcgp(i)==i) then
        ipart_iso=ipart_iso+1
      elseif (hcgp(i)==0) then
        ipart_mem=ipart_mem+1
      else
        ipart_head=ipart_head+1
      endif
    enddo
    print*,'  ipart_iso, ipart_mem, ipart_head'
    print*,ipart_iso, ipart_mem, ipart_head

    call system_clock(t2,t_rate)
    print*, '  percolation time =',real(t2-t1)/t_rate,'secs';

    print*,'fof_member'
    allocate(xi(3,npart/8),xgr(3,npart/2),cm(npart/2))
    cm=0; ngr=0
    do ip=1,npart
      if (hcgp(ip)/=ip .and. hcgp(ip)/=0) then ! not ordered
        grm=0; x1=0
        j=hcgp(ip)
        do while (j/=0)
          grm=grm+1
          xi(:,grm)=xf(:,j)
          !print*,j,xf(:,j)
          j=llgp(j)
        enddo

        x0=xi(:,1) ! set to the first particle
        do j=1,grm
          xi(:,j)=xi(:,j)-x0 ! relative dist
          !xi(:,j)=[0.4,0.51,0.99]
          xi(:,j)=modulo(xi(:,j)+0.5,1.)-0.5 ! topo dist
          !print*,xi(:,j); stop
          x1=x1+xi(:,j)
        enddo
        x1=x1/grm+x0
        x1=modulo(x1,1.)
        ngr=ngr+1
        cm(ngr)=grm
        xgr(:,ngr)=x1
        !print*, ngr,grm,x1
        hcgp(ngr)=hcgp(ip) ! why?
      endif
    enddo
    print*,'  found',ngr,' halos'

    if (ngr>0) then
      print*,'fof_sort'
      allocate(nab(npart/2),wsp(npart/2),iwsp(npart/2))
      nab=0; wsp=0; iwsp=0
      call indexx(ngr,cm,nab)

      cm(1:ngr)=cm(nab(ngr:1:-1))
      hcgp(1:ngr)=hcgp(nab(ngr:1:-1))
      xgr(:,1:ngr)=xgr(:,nab(ngr:1:-1))
      print*,'  most massive cm'
      print*, cm(:10)
      !print*,xgr(:,1:10)

      ngr0=0
      do i=1,ngr
        if (cm(i)>=massll) ngr0=ngr0+1
      enddo
    endif

    print*,'write into file',output_name('halofof')
    print*,'  ngr0 =',ngr0
    open(11,file=output_name('halofof'),status='replace',access='stream')
    write(11) ngr0,b_link ! total number of halos, linking parameter
    write(11) xgr(:,:ngr0) ! write mean location
    write(11) cm(:ngr0) ! write particle number
    close(11)

    open(11,file=output_name('halofof_xfidx'),status='replace',access='stream')
    !write(11) will add stuff soon
    close(11)

    

    deallocate(xf,xf1,hoc,ll,ll1,ll2,llgp,hcgp,ecgp,xi,xgr,cm)
    if (ngr>0) deallocate(nab,wsp,iwsp)

    call system_clock(tt2,t_rate)
    print*, 'total time =',real(tt2-tt1)/t_rate,'secs';
    print*, ''
    print*, ''

  enddo



  print*,'fof done'












  contains


    subroutine merge_group(hcgp,llgp,ecgp,npart,i,j)
      integer npart,i,j,iend,jend,ihead,jhead
      integer hcgp(npart),llgp(npart),ecgp(npart)
      if(hcgp(i)==i)then      !i has been isolated
         if(hcgp(j)==j)then   !j has been isolated
            hcgp(i)=j
            llgp(j)=i
            hcgp(j)=0
            ecgp(i)=i
            ecgp(j)=i
         else if(hcgp(j).eq.0)then !j is in group, but not the
            llgp(i)=llgp(j)     !end of chain
            llgp(j)=i
            hcgp(i)=0
            ecgp(i)=ecgp(j)
         else
            llgp(i)=hcgp(j)     !j is in group and the end
            hcgp(j)=i           !of chain
            hcgp(i)=0
            ecgp(i)=ecgp(j)
         endif
      elseif(hcgp(i).eq.0)then !i is in group but not the end of Ch
         if(hcgp(j).eq.j)then   !j has been isolated
            llgp(j)=llgp(i)
            llgp(i)=j
            hcgp(j)=0
            ecgp(j)=ecgp(i)
         elseif(hcgp(j).eq.0)then !j is in group, but not the
            jend=ecgp(j)        !end of chain
            iend=ecgp(i)
            if(iend.ne.jend)then !if not the same one
               !print*,'not the same one'
               jhead=hcgp(jend) !j group merged into i group
               call chend(llgp,ecgp,npart,jhead,iend)
               llgp(jend)=llgp(i)
               llgp(i)=hcgp(jend)
               hcgp(jend)=0
            endif
         else                   !j is in group and end of chain
            iend=ecgp(i)        !j group merged into i group
            if(iend.ne.j)then   !if not the same one
               jhead=hcgp(j)
               call chend(llgp,ecgp,npart,jhead,iend)
               llgp(j)=llgp(i)
               llgp(i)=hcgp(j)
               hcgp(j)=0
            endif
         endif
      else                      !i is in group and the end of chain
         if(hcgp(j).eq.j)then   !j has been isolated
            llgp(j)=hcgp(i)
            hcgp(i)=j
            hcgp(j)=0
            ecgp(j)=i
         elseif(hcgp(j).eq.0)then !j is in group, but not the
            jend=ecgp(j)        !end of chain
            if(jend.ne.i)then   !i group merged into j group
               ihead=hcgp(i)    !if not the same one
               call chend(llgp,ecgp,npart,ihead,jend)
               llgp(i)=llgp(j)
               llgp(j)=hcgp(i)
               hcgp(i)=0
            endif
         else
            call chend(llgp,ecgp,npart,hcgp(j),i)
            llgp(j)=hcgp(i)     !j is in group and end
            hcgp(i)=hcgp(j)     !of chain
            hcgp(j)=0           !j group merged into i group
         endif
      endif
    endsubroutine merge_group

    subroutine chend(llgp,ecgp,npart,jhead,iend)
      integer llgp(npart),ecgp(npart),j
      integer npart,jhead,iend
      j=jhead
      do while (j/=0)
        ecgp(j)=iend
        j=llgp(j)
      enddo
    endsubroutine chend

    subroutine indexx(n,arrin,indx)
      implicit none
      integer n
      integer indx(n),indxt,ir
      real arrin(n),q

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

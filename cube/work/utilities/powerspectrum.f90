!! module for power spectrum analysis
! uses pencil_fft however work for single image only
! cx1,cx2 can be memory optimized
! auto-power can be memory optimized
! check nexp frequently
#define linear_kbin
!#define pl2d

module powerspectrum
use pencil_fft

#ifdef linear_kbin
  integer(8),parameter :: nbin=nint(nyquist*sqrt(3.))
#else
  integer(8),parameter :: nbin=floor(4*log(nyquist*sqrt(3.)/0.95)/log(2.))
#endif
#ifdef pl2d
  integer kp,kl
  integer nmode(nint(nyquist*sqrt(2.))+1,nyquist+1)
  real pow2d(nint(nyquist*sqrt(2.))+1,nyquist+1)
  real pow2drsd(nint(nyquist*sqrt(2.))+1,nyquist+1)
#endif
complex cx1(ng*nn/2+1,ng,npen),cx2(ng*nn/2+1,ng,npen)


contains

subroutine cross_power(xip,cube1,cube2)
  implicit none

  integer i,j,k,ig,jg,kg,ibin
  real kr,kx(3),sincx,sincy,sincz,sinc,rbin

  real cube1(ng,ng,ng),cube2(ng,ng,ng)
  real xi(10,0:nbin),xip(10,nbin)[*]
  real amp11,amp12,amp22
  complex cx1(ng*nn/2+1,ng,npen),cx2(ng*nn/2+1,ng,npen)

  real,parameter :: nexp=4.0 ! CIC kernel
print*,'in cross_power'
  xi=0

  r3=cube1
  call pencil_fft_forward
  cx1=cxyz/ng_global/ng_global/ng_global

  r3=cube2
  call pencil_fft_forward
  cx2=cxyz/ng_global/ng_global/ng_global
print*,'a'
  xi=0
  sync all
#ifdef pl2d
  print*, 'size of pow2d',nint(nyquist*sqrt(2.))+1,nyquist+1
  pow2d=0
  nmode=0
#endif

  do k=1,npen
  do j=1,ng
  do i=1,nyquist+1
    kg=(nn*(icz-1)+icy-1)*npen+k
    jg=(icx-1)*ng+j
    ig=i
    kx=mod((/ig,jg,kg/)+nyquist-1,ng_global)-nyquist
    if (ig==1.and.jg==1.and.kg==1) cycle ! zero frequency
    if ((ig==1.or.ig==ng*nn/2+1) .and. jg>ng*nn/2+1) cycle
    if ((ig==1.or.ig==ng*nn/2+1) .and. (jg==1.or.jg==ng*nn/2+1) .and. kg>ng*nn/2+1) cycle
    kr=sqrt(kx(1)**2+kx(2)**2+kx(3)**2)
    sincx=merge(1.0,sin(pi*kx(1)/ng_global)/(pi*kx(1)/ng_global),kx(1)==0.0)
    sincy=merge(1.0,sin(pi*kx(2)/ng_global)/(pi*kx(2)/ng_global),kx(2)==0.0)
    sincz=merge(1.0,sin(pi*kx(3)/ng_global)/(pi*kx(3)/ng_global),kx(3)==0.0)
    sinc=sincx*sincy*sincz
#   ifdef linear_kbin
      ibin=nint(kr)
#   else
      rbin=4.0/log(2.)*log(kr/0.95)
      ibin=merge(ceiling(rbin),floor(rbin),rbin<1)
#   endif
    xi(1,ibin)=xi(1,ibin)+1 ! number count
    xi(2,ibin)=xi(2,ibin)+kr ! k count
    amp11=real(cx1(i,j,k)*conjg(cx1(i,j,k)))/(sinc**4.0)*4*pi*kr**3
    amp22=real(cx2(i,j,k)*conjg(cx2(i,j,k)))/(sinc**4.0)*4*pi*kr**3
    amp12=real(cx1(i,j,k)*conjg(cx2(i,j,k)))/(sinc**4.0)*4*pi*kr**3
#ifdef pl2d
    kp=nint(sqrt(kx(1)**2+kx(3)**2))+1
    kl=abs(kx(2))+1
    nmode(kp,kl)=nmode(kp,kl)+1
    pow2d(kp,kl)=pow2d(kp,kl)+amp11
    pow2drsd(kp,kl)=pow2drsd(kp,kl)+amp22
#endif
    xi(3,ibin)=xi(3,ibin)+amp11 ! auto power 1
    xi(4,ibin)=xi(4,ibin)+amp22 ! auto power 2
    xi(5,ibin)=xi(5,ibin)+amp12 ! cross power
    xi(6,ibin)=xi(6,ibin)+1/sinc**2.0 ! kernel 1
    xi(7,ibin)=xi(7,ibin)+1/sinc**4.0 ! kernel 2

  enddo
  enddo
  enddo
  sync all
#ifdef pl2d
  nmode=max(1,nmode)
  pow2d=pow2d/nmode
  pow2drsd=pow2drsd/nmode
  open(55,file=output_name('pow2d'),status='replace',access='stream')
  write(55) pow2d
  write(55) pow2drsd
  close(55)
#endif
sync all
xip=xi(:,1:)

  ! co_sum
  if (head) then
    do i=2,nn**3
      xip=xip+xip(:,:)[i]
    enddo
  endif
  sync all

  ! broadcast
  xip=xip(:,:)[1]
  sync all

  ! divide and normalize
  xip(2,:)=xip(2,:)/xip(1,:)*(2*pi)/box ! k_phy
  xip(3,:)=xip(3,:)/xip(1,:) ! Delta_LL
  xip(4,:)=xip(4,:)/xip(1,:) ! Delta_RR
  xip(5,:)=xip(5,:)/xip(1,:) ! Delta_LR ! cross power
  xip(6,:)=xip(6,:)/xip(1,:) ! kernel
  xip(7,:)=xip(7,:)/xip(1,:) ! kernel
  xip(8,:)=xip(5,:)/sqrt(xip(3,:)*xip(4,:)) ! r
  xip(9,:)=sqrt(xip(4,:)/xip(3,:)) ! b
  xip(10,:)=xip(8,:)**4/xip(9,:)**2 * xip(4,:) ! P_RR*r^4/b^2 reco power

  sync all
endsubroutine

subroutine auto_power(xi,cube1)
  implicit none

  integer i,j,k,ig,jg,kg,ibin
  real kr,kx(3),sincx,sincy,sincz,sinc,rbin

  real cube1(ng,ng,ng)
  real xi(10,nbin)[*]
  complex C_0(nbin)
  real amp11

  real,parameter :: nexp=4.0 ! CIC kernel
  print*,'in auto_power'

  xi=0
  r3=cube1
  call pencil_fft_forward
  cxyz=cxyz/ng_global/ng_global/ng_global
  sync all

  C_0=0
  do k=1,npen
  do j=1,ng
  do i=1,nyquist+1
    kg=(nn*(icz-1)+icy-1)*npen+k
    jg=(icx-1)*ng+j
    ig=i
    kx=mod((/ig,jg,kg/)+nyquist-1,ng_global)-nyquist
    if (ig==1.and.jg==1.and.kg==1) cycle ! zero frequency
    if ((ig==1.or.ig==ng*nn/2+1) .and. jg>ng*nn/2+1) cycle
    if ((ig==1.or.ig==ng*nn/2+1) .and. (jg==1.or.jg==ng*nn/2+1) .and. kg>ng*nn/2+1) cycle
    kr=sqrt(kx(1)**2+kx(2)**2+kx(3)**2)
    sincx=merge(1.0,sin(pi*kx(1)/ng_global)/(pi*kx(1)/ng_global),kx(1)==0.0)
    sincy=merge(1.0,sin(pi*kx(2)/ng_global)/(pi*kx(2)/ng_global),kx(2)==0.0)
    sincz=merge(1.0,sin(pi*kx(3)/ng_global)/(pi*kx(3)/ng_global),kx(3)==0.0)
    sinc=sincx*sincy*sincz
    ibin=nint(kr)
    xi(1,ibin)=xi(1,ibin)+1 ! number count
    xi(2,ibin)=xi(2,ibin)+kr ! k count
    amp11=real(cxyz(i,j,k)*conjg(cxyz(i,j,k)))/(sinc**4.0)*4*pi*kr**3
    xi(3,ibin)=xi(3,ibin)+amp11
    C_0(ibin)=C_0(ibin)+cxyz(i,j,k)/(sinc**2.0)*sqrt(4*pi*kr**3)
    !if (ibin==4) print*, kx,C_0(ibin)
  enddo
  enddo
  enddo
  sync all

  ! co_sum
  if (head) then
    do i=2,nn**3
      xi=xi+xi(:,:)[i]
    enddo
  endif
  sync all

  ! broadcast
  xi=xi(:,:)[1]
  sync all

  ! divide and normalize
  xi(2,:)=xi(2,:)/xi(1,:)*(2*pi)/box ! k_phy
  xi(3,:)=xi(3,:)/xi(1,:)
  xi(4,:)=real(C_0*conjg(C_0))/xi(1,:)
  sync all
  !print*,'        ',C_0(4)

endsubroutine


subroutine density_to_potential(cube1)
  implicit none
  integer i,j,k,ig,jg,kg
  real kr,kx(3),cube1(ng,ng,ng)
  print*,'convert density to potential'
  r3=cube1
  call pencil_fft_forward
  cxyz=cxyz/ng_global/ng_global/ng_global
  sync all

  do k=1,npen
  do j=1,ng
  do i=1,nyquist+1
    kg=(nn*(icz-1)+icy-1)*npen+k
    jg=(icx-1)*ng+j
    ig=i
    kx=mod((/ig,jg,kg/)+nyquist-1,ng_global)-nyquist
    kx=2*sin(pi*kx/ng_global)
    kr=kx(1)**2+kx(2)**2+kx(3)**2
    kr=max(kr,1.0/ng_global**2)
    cxyz(i,j,k) = -4*pi/kr * cxyz(i,j,k)
  enddo
  enddo
  enddo
  cxyz(1,1,1)=0
  sync all

  open(11,file=output_name('phik'),status='replace',access='stream')
    write(11) cxyz
  close(11)
  sync all

  call pencil_fft_backward
  sync all

  open(11,file=output_name('phi'),status='replace',access='stream')
    write(11) r3
  close(11)
  sync all

endsubroutine


endmodule

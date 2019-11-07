use omp_lib
use pencil_fft
implicit none
save

integer i,j,k,ig,jg,kg
real kx,ky,kz,kr

call omp_set_num_threads(ncore)
call geometry
call create_penfft_plan


open(11,file='input.bin',access='stream',status='old')
read(11) r3
close(11)

!! this transfroms r3 into Fourier space, stored as cxyz
call pencil_fft_forward

!! i,j,k are the coordinates in the node
do k=1,npen
do j=1,nf
do i=1,nyquest+1
  !! ig,jg,kg are the global coordinates, labeled from the first node
  kg=(nn*(icz-1)+icy-1)*npen+k
  jg=(icx-1)*nf+j
  ig=i
  !! kx,ky,kz are Fourier frequencies
  kz=mod(kg+nyquest-1,nf_global)-nyquest
  ky=mod(jg+nyquest-1,nf_global)-nyquest
  kx=ig-1
  !! can optionally turn off the sinc function
  kz=2*sin(pi*kz/nf_global)
  ky=2*sin(pi*ky/nf_global)
  kx=2*sin(pi*kx/nf_global)
  !! kr is the k^2 we are using
  kr=kx**2+ky**2+kz**2
  kr=max(kr,1.0/nf_global**2) ! avoid kr being 0
  cxyz(i,j,k)=-4*pi/kr
enddo
enddo
enddo

if (head) cxyz(1,1,1)=0 ! DC frequency

sync all

!! this transforms cxyz into real space, stored as r3
call pencil_fft_backward

open(11,file='output.bin',status='replace',access='stream')
write(11) r3
close(11)

end

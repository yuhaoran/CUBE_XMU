use omp_lib
  integer,parameter :: n=20000
  integer,parameter :: m=2
  integer,parameter :: l=2
  real vfield(3,n,n,n,m,m,m)[l,l,*]
  real,allocatable :: phi(:,:,:,:,:,:)[:,:,:]


  allocate(phi(n,n,n,m,m,m)[l,l,*])
  vfield=0;
  print*, vfield(:,1,1,1,1,1,1)
  deallocate(phi)
end

program read
implicit none

   integer, parameter :: FILES = 1          ! number of files per snapshot

   character*200 filename, fnumber

   character*5, parameter ::       zz='3.000' 
   character*200, parameter :: input_xp=zz//'_xp_1.bin'
   character*200, parameter :: input_vp=zz//'_vp_1.bin'
   character*200, parameter :: input_np=zz//'_np_1.bin'
   character*200, parameter :: input_vc=zz//'_vc_1.bin'
   character*200, parameter :: input_id=zz//'_id_1.bin'
   character*200, parameter :: input_info=zz//'_info_1.bin'
   character*200, parameter :: output=zz//'_tmp.bin'
   integer(8) itx,ity,itz,i,j,k
   integer(8) nti,ntj,ntk
   integer(8) nlast,ip,np,l

   ! For pos:
   integer(8),parameter :: izipx = 2
   real(8),parameter :: x_resolution=1.0/(int(2,8)**(izipx*8))
   integer(8),parameter :: ishift=-(int(2,8)**(izipx*8-1))
   real(8),parameter :: rshift=0.5-ishift
   integer(izipx),allocatable :: xp(:,:)
   real(4),allocatable :: pos(:,:)

   ! For vel:
   integer(8),parameter :: izipv = 2
   integer(8),parameter :: nvbin=int(2,8)**(8*izipv)
   real(8),parameter :: vrel_boost=2.5
   real,parameter :: pi=4*atan(1.)
   integer(izipv),allocatable :: vp(:,:)
   real(4),allocatable :: vfield(:,:,:,:,:,:,:)
   real*4,allocatable    :: vel(:,:)
 
   ! For density on mesh:
   integer(4),allocatable :: rhoc(:,:,:,:,:,:)

   type sim_header
    integer(8) nplocal,npglobal,nplocal_nu,npglobal_nu
    integer(8) izipx,izipv,izipx_nu,izipv_nu
    integer(8) image
    integer(8) nn,nnt,nt,ncell,ncb
    integer(8) timestep
    integer(8) cur_checkpoint
    integer(8) cur_halofind

    real a, t, tau
    real dt_pp, dt_fine, dt_coarse, dt_vmax, dt_vmax_nu
    real mass_p_cdm,mass_p_nu
    real box

    real h0
    real omega_m
    real omega_l
    real s8
    real vsim2phys
    real sigma_vres
    real sigma_vi
    real sigma_vi_nu
    real z_i,z_i_nu
    real vz_max
  endtype
  type(sim_header) sim

  integer(8),parameter :: nc=128 ! nc/image/dim, in physical volume, >=24
  integer(8) ng, nf_global

! ----For Gadget ---
   integer*4 npart(0:5) ! In total, 6 types: 0, Gas/ 1, Halo/ 2, Disk/ 3, Bulge/ 4, Stars/ 5, Bndry
   real*8    massarr(0:5)
   real*8    a
   real*8    redshift
   integer*4 flag_sfr,flag_feedback
   integer*4 nall(0:5)
   integer*4 flag_cooling,num_files
   real*8    BoxSize,Omega0,OmegaLambda,HubbleParam
   character*(256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8) fill
   integer*8 I8
   integer*4 I4
   integer*4 fn
   fn=0
  ! ------ To convert into global coord.
  !integer*8, allocatable   :: ID(:)

  open(unit=11, file=input_info,status='old',access='stream')
  read(11) sim
  close(11)
  print*,'box        = ',sim%box, ' Mpc/h'
  print*,'h0         = ',sim%h0,' km/s/Mpc' 
  print*,'a          = ',sim%a
  print*,'vsim2phys  = ',sim%vsim2phys, ' (km/s)/(1.0)'
  print*,'omega_m    = ',sim%omega_m
  print*,'omega_l    = ',sim%omega_l
  print*,'z_i        = ',sim%z_i
  print*,'z_i_nu     = ',sim%z_i_nu
  print*,'nplocal    = ',sim%nplocal
  print*,'npglobal   = ',sim%npglobal
  print*,'nn         = ',sim%nn
  print*,'nnt        = ',sim%nnt
  print*,'ncell      = ',sim%ncell
  print*,'ncb        = ',sim%ncb
  print*,'nt         = ',sim%nt
  print*,'izipx      = ',sim%izipx
  print*,'izipv      = ',sim%izipv
  print*,'izipx_nu   = ',sim%izipx_nu
  print*,'izipv_nu   = ',sim%izipv_nu
  print*,'sigma_vi   =',sim%sigma_vi, '(simulation unit)'
  print*,'sigma_vi_nu=',sim%sigma_vi_nu, '(simulation unit)'
  nf_global=nc*sim%ncell*sim%nn
  BoxSize = dble(sim%box)*1000.d0 ! in kpc/h for Gadget
  Omega0 = dble(sim%omega_m)
  OmegaLambda = dble(sim%omega_l)
  HubbleParam = dble(sim%h0)
  a = dble(sim%a)
  redshift = 1.d0/a - 1.d0

  do i=0,5
     massarr(i)=0.d0
  enddo
  num_files = INT(FILES,KIND(I4))
  flag_sfr=0
  flag_feedback=0
  flag_cooling=0

  nall(0)=0
  nall(1)=INT(sim%npglobal,KIND(I8))
  nall(2)=0
  nall(3)=0
  nall(4)=0
  nall(5)=0

  npart(0)=0
  npart(1)=INT(sim%nplocal,KIND(I8))
  npart(2)=0
  npart(3)=0
  npart(4)=0
  npart(5)=0

  do i=0,5
     print *,'Type=',i,'/    Particles (nall)=', nall(i)
  end do

  print *,'HubbleParam =',HubbleParam
  print *,'Omega0      =',Omega0
  print *,'OmegaLambda =',OmegaLambda
  print *,'num_files   =',num_files
  print *,'BoxSize     =',BoxSize

  print *, ' '
  print *,'*** For pos ***'
  print *, ' '
  allocate(rhoc(sim%nt,sim%nt,sim%nt,sim%nnt,sim%nnt,sim%nnt))

  open(unit=13,file=input_np,status='old',access='stream')
  read(13)rhoc
  close(13)
  print*, ' The # of ptl in First 5 coarse grid number: ',rhoc(1,1,1:5,1,1,1)
  print*, 'The # of ptl in Middle 5 coarse grid number: ',rhoc(sim%nt,sim%nt,sim%nt-4:sim%nt,1,1,2)
  print*, '  The # of ptl in Last 5 coarse grid number: ',rhoc(sim%nt,sim%nt,sim%nt-4:sim%nt,2,2,2)
  print*, 'sum(rhoc) = ',sum(rhoc)


  allocate(xp(3,sim%nplocal))
  open(unit=12, file=input_xp,status='old',access='stream')!By default, form='unformatted'
  read(12) xp
  close(12)
  print*, 'Before conversion......'
  print*, 'First 50 xp in x-aixs: ',xp(1,1:50)
  print*, ' Last 50 xp in x-aixs: ',xp(1,sim%nplocal-49:sim%nplocal)
  print*, 'max xp',maxval(xp)
  print*, 'min xp',minval(xp)

  print*, 'To convert xp into global coordinate...'
  allocate(pos(1:3,1:sim%nplocal))

  print*, 'Start to find global coordinate, in each image...!'
  nlast=0
  do itx=1,sim%nnt ! loop over tile. In case you have many images, you need another loop for image
  do ity=1,sim%nnt
  do itz=1,sim%nnt
     do nti=1,sim%nt ! loop over coarse cell in each tile
     do ntj=1,sim%nt
     do ntk=1,sim%nt
        np=rhoc(nti,ntj,ntk,itx,ity,itz)
        do l=1,np
          ip=nlast+l
          pos(1,ip)=sim%nt*(itx-1) + (nti-1) + (int(xp(1,ip)+ishift,izipx)+rshift)*x_resolution
          pos(2,ip)=sim%nt*(ity-1) + (ntj-1) + (int(xp(2,ip)+ishift,izipx)+rshift)*x_resolution
          pos(3,ip)=sim%nt*(itz-1) + (ntk-1) + (int(xp(3,ip)+ishift,izipx)+rshift)*x_resolution
        enddo
        nlast=nlast+np
     enddo
     enddo
     enddo
  enddo
  enddo
  enddo
  !deallocate(rhoc)
  deallocate(xp)

  print*, 'After conversion ......'
  print*, 'First 50 pos in x-aixs: ',pos(1,1:50)
  print*, ' Last 50 pos in x-aixs: ',pos(1,sim%nplocal-49:sim%nplocal)

  print*, 'max readin x',maxval(pos)
  print*, 'min readin x',minval(pos)
  do i=1,sim%nplocal   ! Mpc -> kpc
     pos(1,i) = pos(1,i)*sim%box/sim%nt/sim%nnt*1000.d0
     pos(2,i) = pos(2,i)*sim%box/sim%nt/sim%nnt*1000.d0
     pos(3,i) = pos(3,i)*sim%box/sim%nt/sim%nnt*1000.d0
  enddo
  print*, 'First 50 pos in x-aixs: ',pos(1,1:50)
  print*, ' Last 50 pos in x-aixs: ',pos(1,sim%nplocal-49:sim%nplocal)

  print*, 'max readin x',maxval(pos)
  print*, 'min readin x',minval(pos)

  print *,' '
  print *,' *** pos Done *** '
  print *, ' '
  deallocate(rhoc)


  print*, 'Start to save Gadget format!'
  write(fnumber,'(I3)') fn
  filename= output(1:len_trim(output)) // '.' // fnumber(verify(fnumber,' '):3)
  print*,filename
  open(unit=12, file=filename,status='replace',form='unformatted')
  write(unit=12) npart, massarr, a, redshift, flag_sfr, flag_feedback, nall, flag_cooling, num_files, BoxSize, Omega0, OmegaLambda, HubbleParam, fill
  write(unit=12)pos
  close(12)
  deallocate(pos)
  print *,' '
  print *,' *** All Done ***'
  print *, ' '



!================================================================================
end program read

program cube_write_snapshot
implicit none

! Number of particles, snapshot Id 
integer(4) np,ips

! current redshift, current omega_m, current Omega_lambda, boxsize in unit of Mpc/h, xscale, vscale （xscale和vscale给定任意数都可以，比如1.0）     
real(4) ztp,omgt,lbdt,boxsize,xscale,vscale

！以上为头文件需要写入的参数

!particles position
real(4) pos(np, 3)
!particles velocity
real(4) vel(np, 3)

! 文件命名规则，比如8213模拟：pos8213.0780     vel8213.4333
open(11,file=SnapshotPath+'pos'+'SimulationName'+'.SnapshotId',access='stream')
open(12,file=SnapshotPath+'vel'+'SimulationName'+'.SnapshotId',access='stream')

write(11) np,ips,ztp,omgt,lbdt,boxsize,xscale,vscale
write(12) np,ips,ztp,omgt,lbdt,boxsize,xscale,vscale


! 粒子的位置要统一成相对位置（x/BoxSize), 粒子的速度比较奇怪，
! 要改写为：V (km/s) / (H0*L*(H/H0)*a*a*R0),
! 其中L: BoxSize, a: 1/(1+z), R0: 1+z_Ini (1+初始红移)
write(11) pos
write(12) vel

close(11)
close(12)

end


program cube_write_fof
implicit none

! 应该是fof找halo的参数，在jing的格式中 b=0.2
real(4) b

！group的数目
integer(4) ngroup

！以上为头文件需要写入的参数

! 一个group中的粒子数
integer(8) NpInOneGroup

! 一个group中的粒子ID
integer(4) Pid(NpInOneGroup)

! 文件命名规则，比如8213模拟：fof.b20.8213.5000
open(13,file=HaloPath+'fof.b20.SimulationName.SnapshotId',access='stream')

! 写头文件
write(13) b, ngroup

! 接下来，HBT在读group文件的时候跳过了前三行，所以这三行应该可以随便写，但我还没弄清楚原因。
do i=1, 3
 	write(13)
endo

！接下来，按照每个group中先写粒子数，后写粒子ID的方式，循环写入

do i=1, ngroup
	write(13) NpInOneGroup
	write(13) Pid(NpInOneGroup)
endo

close(13)

end
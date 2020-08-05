function [ sim ] = get_sim_info( prefix )
%% read info
fid=fopen([prefix,'info_1.bin']);
  disp([prefix,'info_1.bin'])
  sim.nplocal=fread(fid,1,'integer*8');
  sim.npglobal=fread(fid,1,'integer*8');
  sim.nplocal_nu=fread(fid,1,'integer*8');
  sim.npglobal_nu=fread(fid,1,'integer*8');

  sim.izipx=fread(fid,1,'integer*8');
  sim.izipv=fread(fid,1,'integer*8');
  sim.izipx_nu=fread(fid,1,'integer*8');
  sim.izipv_nu=fread(fid,1,'integer*8');

  sim.image=fread(fid,1,'integer*8');

  sim.nn=fread(fid,1,'integer*8');
  sim.nnt=fread(fid,1,'integer*8');
  sim.nt=fread(fid,1,'integer*8');
  sim.ncell=fread(fid,1,'integer*8');
  sim.ncb=fread(fid,1,'integer*8');

  sim.istep=fread(fid,1,'integer*8');
  sim.cur_checkpoint=fread(fid,1,'integer*8');
  sim.cur_halofind=fread(fid,1,'integer*8');

  sim.a=fread(fid,1,'real*4');
  sim.t=fread(fid,1,'real*4');
  sim.tau=fread(fid,1,'real*4');
  sim.dt=fread(fid,5,'real*4');
  sim.mass_p=fread(fid,1,'real*4');
  sim.mass_p_nu=fread(fid,1,'real*4');
  sim.box=fread(fid,1,'real*4');

  sim.h0=fread(fid,1,'real*4');
  sim.omega_m=fread(fid,1,'real*4');
  sim.omega_l=fread(fid,1,'real*4');
  sim.s8=fread(fid,1,'real*4');
  sim.vsim2phys=fread(fid,1,'real*4');
  sim.sigma_vres=fread(fid,1,'real*4');
  sim.sigma_vi=fread(fid,1,'real*4');
  sim.sigma_vi_nu=fread(fid,1,'real*4');
  sim.z_i=fread(fid,1,'real*4');
  sim.z_i_nu=fread(fid,1,'real*4');
  sim.vz_max=fread(fid,1,'real*4');
fclose(fid);

sim.nc=sim.nt*sim.nnt;
sim.nf=sim.nc*sim.ncell;
G=6.67408e-11; % in m^3 kg^-1 s^-2
rho_crit=30000/8/pi/G/(3.086e19)^2; % in kg m^-3 h^-2
rho_crit=rho_crit/1.98855e30*(3.086e22)^3; % in M_solar Mpc^-3 h^-2
sim.mass_p_solar=rho_crit*sim.omega_m*sim.box^3/sim.npglobal/sim.h0;

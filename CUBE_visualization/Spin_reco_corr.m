clear; SetDefault
format short
global Redshift Redshift_i Universe Path Dir
Universe='1';
Redshift_i='100.000';
Redshift='0.000';
%Path='/Users/haoran/cloud/CUBE_spin/CUBE_v2.0/output/';
Path='../CUBE_v2.0/output/';
Dir=['universe',Universe,'/image1/'];
sim=get_sim_info([Path,Dir,Redshift_i,'_']);
ng=sim.nf;
disp('-------------------------------------------------------------------')
disp('nf ='); disp(sim.nf)
disp('box ='); disp(sim.box)
disp('mass_p ='); disp(sim.mass_p_solar)
%%
fid=fopen([Path,Dir,Redshift,'_lptcorr_i_1.bin'],'r');
nmassbin=fread(fid,1,'integer*4');
n_rsmall=fread(fid,1,'integer*4');
n_ratio=fread(fid,1,'integer*4');
imass_info=fread(fid,[4,nmassbin],'real*4');
rsmall=fread(fid,n_rsmall,'real*4');
ratio_scale=fread(fid,n_ratio,'real*4');
mass_info=imass_info(1:3,:)*sim.mass_p_solar;
%%
xx=[rsmall(1),rsmall(n_rsmall)];
yy=[ratio_scale(1),ratio_scale(n_ratio)];

for imassbin=1:nmassbin
  corr_t=fread(fid,[n_rsmall,n_ratio],'real*4');
  corr_q=fread(fid,[n_rsmall,n_ratio],'real*4');
  corr_x=fread(fid,[n_rsmall,n_ratio],'real*4');
  
  figure(imassbin)
  
  subplot(3,1,1)
  imagesc(xx,yy,corr_t'); axis xy; colorbar; label('$r\, [{\rm Mpc}/h]$','$R/r$'); 
  %colormap(1-gray);caxis([-0.01,0.01])
  h=title(['$j_T$ for $M/M_\odot\simeq$ ',num2str(mass_info(3,imassbin),'%5.1E')]); hold on
  h.FontSize=16;
  hold on; [corr_max,loc]=maxloc(corr_t);
  plot(rsmall(loc(1)),ratio_scale(loc(2)),'r*');
  text(rsmall(loc(1)),ratio_scale(loc(2)),num2str(corr_max)); %caxis([0,0.2])
  
  subplot(3,1,2)
  imagesc(xx,yy,corr_q'); axis xy; colorbar; label('$r\, [{\rm Mpc}/h]$','$R/r$'); 
  %colormap(1-gray);caxis([-0.01,0.01])
  h=title('$j_q$'); hold on
  h.FontSize=16;
  hold on; [corr_max,loc]=maxloc(corr_q);
  plot(rsmall(loc(1)),ratio_scale(loc(2)),'r*'); 
  text(rsmall(loc(1)),ratio_scale(loc(2)),num2str(corr_max)); %caxis([0,0.2])
  
  subplot(3,1,3)
  imagesc(xx,yy,corr_x'); axis xy; colorbar; label('$r\, [{\rm Mpc}/h]$','$R/r$'); 
  %colormap(1-gray);caxis([-0.01,0.01])
  h=title('$j_0$'); hold on
  h.FontSize=16;
  hold on; [corr_max,loc]=maxloc(corr_x);
  plot(rsmall(loc(1)),ratio_scale(loc(2)),'r*'); 
  text(rsmall(loc(1)),ratio_scale(loc(2)),num2str(corr_max)); %caxis([0,0.2])
end
fclose(fid);


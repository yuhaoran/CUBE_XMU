clear; SetDefault
format short
global Redshift Redshift_i Universe Path Dir
Universe='2';
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
fid=fopen([Path,Dir,Redshift,'_lptcorr_e_1.bin'],'r');
nmassbin=fread(fid,1,'integer*4');
n_rsmall=fread(fid,1,'integer*4');
n_ratio=fread(fid,1,'integer*4');
imass_info=fread(fid,[4,nmassbin],'real*4');
rsmall=fread(fid,n_rsmall,'real*4');
ratio_scale=fread(fid,n_ratio,'real*4');
mass_info=imass_info(1:3,:)*sim.mass_p_solar;
%%
if 1 % plot r_opt vs mass
  figure
  xx=[rsmall(1),rsmall(n_rsmall)];
  yy=[log10(mass_info(3,2)),log10(mass_info(3,nmassbin))];
  corr_x=zeros(n_rsmall,nmassbin); corr_t=corr_x; corr_q=corr_x;
  corrmax=zeros(nmassbin,1);loc=corrmax;
  for imassbin=1:nmassbin
    temp=fread(fid,[n_rsmall,n_ratio],'real*4');
    corr_t(:,imassbin)=temp(:,1);
    temp=fread(fid,[n_rsmall,n_ratio],'real*4');
    corr_q(:,imassbin)=temp(:,1);
    temp=fread(fid,[n_rsmall,n_ratio],'real*4');
    corr_x(:,imassbin)=temp(:,1);
  end
  fclose(fid);
  ngp=sim.nf^3/sim.nplocal;
  rq=sim.box/sim.nf*(3*ngp*imass_info(3,1:nmassbin)/4/pi).^(1/3.);
  
  subplot(3,2,1)
  imagesc(xx,yy,corr_t(:,2:nmassbin)'); hold on; axis xy; colorbar; label('$r\, [{\rm Mpc}/h]$','$\log_{10}(M/M_\odot)$'); colormap(1-gray);caxis([0,0.6])
  h=title('$j_T$'); h.FontSize=16;
  for imassbin=1:nmassbin
    [corrmax(imassbin),loc(imassbin)]=max(corr_t(:,imassbin));
  end
  plot(rsmall(loc(2:nmassbin)),log10(mass_info(3,2:nmassbin)),'r--',rq(2:nmassbin),log10(mass_info(3,2:nmassbin)),'y--')
  subplot(3,2,2)
  [hx,hy1,hy2]=plotyy(rq(2:nmassbin),rsmall(loc(2:nmassbin)),rq(2:nmassbin),corrmax(2:nmassbin)); hold on
  hx(2).YLim=[0 1];
  %hx(1).XLim=[0,max(rsmall)]; hx(2).XLim=[0,max(rsmall)];
  %ax=gca; ax.YLim=[0 1];
  hy1.LineStyle='-'; hy1.Marker='.'; hy1.MarkerSize=10;
  hy2.LineStyle='none'; hy2.Marker='o';
  ylabel(hx(1),'$r\, [{\rm Mpc}/h]$')
  ylabel(hx(2),'maximum correlation')
  xlabel('$r_q\, [{\rm Mpc}/h]$')
  plot([0,max(rq)],[0,max(rq)],'k--')
 
 
 
  subplot(3,2,3)
  imagesc(xx,yy,corr_q(:,2:nmassbin)'); hold on; axis xy; colorbar; label('$r\, [{\rm Mpc}/h]$','$\log_{10}(M/M_\odot)$'); colormap(1-gray);caxis([0,0.6])
  h=title('$j_q$'); h.FontSize=16;
  for imassbin=1:nmassbin
    [corrmax(imassbin),loc(imassbin)]=max(corr_q(:,imassbin));
  end
  plot(rsmall(loc(2:nmassbin)),log10(mass_info(3,2:nmassbin)),'r--',rq(2:nmassbin),log10(mass_info(3,2:nmassbin)),'y--')
  subplot(3,2,4)
  [hx,hy1,hy2]=plotyy(rq(2:nmassbin),rsmall(loc(2:nmassbin)),rq(2:nmassbin),corrmax(2:nmassbin)); hold on
  hx(2).YLim=[0 1];
  %hx(1).XLim=[0,max(rsmall)]; hx(2).XLim=[0,max(rsmall)];
  %ax=gca; ax.YLim=[0 1];
  hy1.LineStyle='-'; hy1.Marker='.'; hy1.MarkerSize=10;
  hy2.LineStyle='none'; hy2.Marker='o';
  ylabel(hx(1),'$r\, [{\rm Mpc}/h]$')
  ylabel(hx(2),'maximum correlation')
  xlabel('$r_q\, [{\rm Mpc}/h]$')
  plot([0,max(rq)],[0,max(rq)],'k--')
  
  subplot(3,2,5)
  imagesc(xx,yy,corr_x(:,2:nmassbin)'); hold on; axis xy; colorbar; label('$r\, [{\rm Mpc}/h]$','$\log_{10}(M/M_\odot)$'); colormap(1-gray);caxis([0,0.6])
  h=title('$j_0$'); h.FontSize=16;
  for imassbin=1:nmassbin
    [corrmax(imassbin),loc(imassbin)]=max(corr_x(:,imassbin));
  end
  plot(rsmall(loc(2:nmassbin)),log10(mass_info(3,2:nmassbin)),'r--',rq(2:nmassbin),log10(mass_info(3,2:nmassbin)),'y--')
  subplot(3,2,6)
  [hx,hy1,hy2]=plotyy(rq(2:nmassbin),rsmall(loc(2:nmassbin)),rq(2:nmassbin),corrmax(2:nmassbin)); hold on
  hx(2).YLim=[0 1];
  %hx(1).XLim=[0,max(rsmall)]; hx(2).XLim=[0,max(rsmall)]; 
  %ax=gca; ax.YLim=[0 1];
  hy1.LineStyle='-'; hy1.Marker='.'; hy1.MarkerSize=10;
  hy2.LineStyle='none'; hy2.Marker='o';
  ylabel(hx(1),'$r\, [{\rm Mpc}/h]$')
  ylabel(hx(2),'maximum correlation')
  xlabel('$r_q\, [{\rm Mpc}/h]$')
  plot([0,max(rq)],[0,max(rq)],'k--')
else % plot for each mass bin
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
end

clear; SetDefault
format short
global Redshift Redshift_i Universe Path Dir
Universe='12';
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
%% autopower
n_row_xi=10;
fid=fopen([Path,Dir,Redshift,'_cicpower_1.bin']);
  xi=fread(fid,'real*4')';
fclose(fid);
xi=reshape(xi,n_row_xi,numel(xi)/n_row_xi)';
camb=load('CAMB/pkz0.txt');
figure(2); loglog(camb(:,1),camb(:,2).*camb(:,1).^3/(2*pi^2),'--',xi(:,2),xi(:,3))
grid on; hold on
%% phi
if 1
  phi1=loadfield3d([Path,Dir,Redshift_i,'_phi1_1.bin']);
  xgrid=[0.5,sim.box-0.5];
  figure; imagesc(xgrid,xgrid,reshape(-mean(phi1(:,:,:),3),ng,ng)'); hold on
  axis xy square; colorbar; title('$-\phi_1$')
end
%% delta_L
if 1
  delta_L=loadfield3d([Path,Dir,'delta_L_1.bin']);
  xgrid=[0.5,sim.box-0.5];
  figure; imagesc(xgrid,xgrid,reshape(mean(delta_L(:,:,:),3),ng,ng)'); hold on
  axis xy square; colorbar; caxis([-3,3]); title('$\delta_L$')
  %colormap(1-gray)
end
%% delta_c
if 1
  Redshift='0.000';
  delta_c=loadfield3d([Path,Dir,Redshift,'_delta_c_1.bin']);
  xgrid=[0.5,sim.box-0.5];
  figure; imagesc(xgrid,xgrid,reshape(mean(delta_c(:,:,:),3),ng,ng)'); hold on
  axis xy square; colorbar; caxis([-1,3]); title('$\delta_c$'); colormap(1-gray);
end
%% delta_E
if 1
  delta_E=loadfield3d([Path,Dir,Redshift,'_delta_E_1.bin']);
  xgrid=[0.5,sim.box-0.5];
  figure; imagesc(xgrid,xgrid,reshape(mean(delta_E(:,:,:),3),ng,ng)'); hold on
  axis xy square; colorbar; caxis([-3,3]); title('$\delta_E$')
end
%colormap(1-gray);
%% FoF halo catalog
Redshift='0.000';
fid=fopen([Path,Dir,Redshift,'_fof_1.bin']);
  nhalo_tot=fread(fid,1,'integer*4')';
  nhalo=fread(fid,1,'integer*4')';
  ninfo=fread(fid,1,'integer*4');
  linking_parameter=fread(fid,1,'real*4')';
  hcat=fread(fid,[ninfo,nhalo],'real*4');
fclose(fid);
n_plot=100;
if 1
  figure
  h1=plot3(hcat(2,:),hcat(3,:),hcat(4,:),'k.'); hold on
  h1.MarkerSize=0.1;
  plot3(hcat(2,1:n_plot),hcat(3,1:n_plot),hcat(4,1:n_plot),'ro')
  view([-30 30]);  xlabel('X'); ylabel('Y'); zlabel('Z'); axis equal
  box on
end
% alignment statistics
n1=1; n2=nhalo;
disp('=================')
disp('correlations')
mean(dot(hcat(14:16,n1:n2),hcat(17:19,n1:n2)))
mean(dot(hcat(14:16,n1:n2),hcat(20:22,n1:n2)))
mean(dot(hcat(17:19,n1:n2),hcat(20:22,n1:n2)))
%return

% b
n1=1; n2=nhalo;
disp('Lagrangian b_1, b_2, b_3')
disp(mean(hcat(14:16,n1:n2).^2,2)-1/3)
disp('\pm')
disp(1/3/sqrt(n2-n1+1))

disp('Eulerian b_1, b_2, b_3')
disp(mean(hcat(17:19,n1:n2).^2,2)-1/3)
disp('\pm')
disp(1/3/sqrt(n2-n1+1))

disp('TTT b_1, b_2, b_3')
disp(mean(hcat(20:22,n1:n2).^2,2)-1/3)
disp('\pm')
disp(1/3/sqrt(n2-n1+1))

return

%% plot histogram of theta(j_q,j_0); test rebin first
hmass_min=hcat(1,nhalo)*sim.mass_p_solar;
hmass_max=hcat(1,1)*sim.mass_p_solar;
%hmass_lim=[hmass_min,1e12,2e12,4e12,8e12,1.6e13,3e13,1e14,2e14,hmass_max];
hmass_lim=[hmass_min,1e11,2e11,5e11,1e12,2e12,4e12,8e12,2e13,hmass_max];
nhbin=numel(hmass_lim)-1;
nlim=zeros(nhbin+1,1); 
for ihbin=1:nhbin+1
  nlim(ihbin)=sum(hcat(1,:)*sim.mass_p_solar>=hmass_lim(ihbin));
end
nlim
%% plot hist
%nlim(1)=1; nlim(nhbin+1)=nhalo;
nmu=30;
for ihbin=1:nhbin
  i1=nlim(ihbin+1);
  i2=nlim(ihbin);
  mu=dot(hcat(14:16,i1:i2),hcat(17:19,i1:i2));
  mu_mean=mean(mu);
  [pdfy,pdfx]=hist(mu,nmu);
  pdfy=pdfy/(nlim(ihbin)-nlim(ihbin+1)+1)/(2/nmu);
  h=plot(pdfx,pdfy); hold on
  hcolor=[1 1 1]-[0.06 0.1 0.06]*ihbin;
  h.Color=hcolor; h.LineWidth=2.5;
  lname=['$M/M_\odot\in\,$[',num2str(hmass_lim(ihbin),'%6.1e'),',',num2str(hmass_lim(ihbin+1),'%6.1e'),'], $\langle\mu\rangle=$',num2str(mu_mean,'%6.2f')];
  eval(['t',num2str(ihbin),'=lname;'])
end
h=legend(t1,t2,t3,t4,t5,t6,t7,t8,t9); h.FontSize=14; h.Location='northwest';
h=xlabel('$\mu$'); h.FontSize=14;
h=ylabel('PDF'); h.FontSize=14;

simname=['Simulation: $L=$',num2str(sim.box),'Mpc/$h$, $N_p=$',num2str(sim.nf),'$^3$; standard'];
h=title(simname); h.FontSize=14;
hold off
%% protohalo shape
tau=sum(hcat(71:73,:),1);
ellip=(hcat(71,:)-hcat(73,:))./(tau*2);
prol=(hcat(71,:)-2*hcat(72,:)+hcat(73,:))./(tau*2);
figure; plot(ellip,prol,'.')
mu=zeros(3,nhalo);
mu(1,:)=dot(hcat(14:16,:),hcat(17:19,:));
mu(2,:)=dot(hcat(14:16,:),hcat(20:22,:));
mu(3,:)=dot(hcat(17:19,:),hcat(20:22,:));

figure; plot3(ellip,prol,mu(1,:),'.'); box on; 
xlabel('$e$'),ylabel('$p$'),zlabel('$\mu$')
%% FoF PID
n_plot=1;
fid=fopen([Path,Dir,Redshift,'_fofpid_1.bin']);
  nhalo_fof=fread(fid,1,'integer*4')';
  for ihalo=1:n_plot
    nphalo=fread(fid,1,'integer*4')';
    pidhalo=fread(fid,nphalo,'integer*4')';
    qpos=zeros(3,nphalo); % create Lagrangian region array
    pidhalo=pidhalo-1; % qid will start from 0
    dx_mean=0;
    for ip=1:nphalo
      qpos(3,ip)=floor(pidhalo(ip)/ng^2);
      qpos(2,ip)=floor((pidhalo(ip)-qpos(3,ip)*ng^2)/ng);
      qpos(1,ip)=mod(pidhalo(ip),ng);
      qpos(:,ip)=(qpos(:,ip)+0.5)*sim.box/ng;
      %idx=ceil(qpos(1,ip))*2-1;
      %idy=ceil(qpos(2,ip))*2-1;
      %if ihalo==n_halo_plot
      %  rhoq(idx:idx+1,idy:idy+1)=rhoq(idx:idx+1,idy:idy+1)+1;
      %end
      %dx=qpos(:,ip)-haloinfo(7:9,ihalo);
      %dx=mod(dx+ng/2,ng)-ng/2;
      %dx_mean=dx_mean+dx;
    end
    plot3(qpos(1,:),qpos(2,:),qpos(3,:),'.')
  end
fclose(fid);
return
%% SO halo catalog
ninfo=42;
fid=fopen([Path,Dir,Redshift,'_halo_1.bin']);
  nhalo_tot=fread(fid,1,'integer*4')';
  nhalo=fread(fid,1,'integer*4')';
  halo_odc=fread(fid,1,'real*4')';
  haloinfo=fread(fid,[ninfo,nhalo],'real*4');
fclose(fid);
plot(haloinfo(7,:)*sim.box/ng,haloinfo(8,:)*sim.box/ng,'r.')
% hpos 1:3; mass 4; r 5; v_disp 6;
% x_mean 7:9; v_mean 10:12; ang_mom 13:15; var_x 16:18; inertia 19:27;
% q_mean 28:30; inertia_q 31:33; s_mean 

for ihalo=1:nhalo
  haloinfo(13:15,ihalo)=haloinfo(13:15,ihalo)/norm(haloinfo(13:15,ihalo));
end
disp(mean(haloinfo(13:15,:).^2,2))

%% HMF
nbin=10;
subplot(2,1,2)

fid=fopen('angulo12.bin','r');
hmf=fread(fid,'real*4'); fclose(fid); hmf=reshape(hmf,[numel(hmf)/2,2]);
loglog(hmf(:,1),hmf(:,2));hold on

fid=fopen('press74.bin','r');
hmf=fread(fid,'real*4'); fclose(fid); hmf=reshape(hmf,[numel(hmf)/2,2]);
loglog(hmf(:,1),hmf(:,2));hold on

fid=fopen('bocquet16.bin','r');
hmf=fread(fid,'real*4'); fclose(fid); hmf=reshape(hmf,[numel(hmf)/2,2]);
loglog(hmf(:,1),hmf(:,2));hold on

legend('Angulo12 FoF','Press74 FoF','Bocquet16 SO')

xlabel('$M/M_\odot$'); ylabel('$dN/d\log_{10}M$')
hmass=haloinfo_fof(4,:)*sim.mass_p_solar;

mbin_cube=logspace(log10(min(hmass)),log10(max(hmass)),nbin);

[pdfy,pdfx]=hist(log10(hmass),50);
pdfy=pdfy/(pdfx(2)-pdfx(1))/sim.box^3*sim.h0^3;
loglog(10.^pdfx,pdfy,'.','MarkerSize',15);hold on
xlabel('$M/M_\odot$'); ylabel('$dN/d\log_{10}M$')
grid on
ax1=gca;
ax1.XLim=[1e10,2e15];
ax1.YLim=[1e-8,2];
return
%% spin field from halos
spin=loadvectorfield3d([Path,Dir,Redshift,'_spin_halo_1.bin']);

for idim=1:3
  figure(100+idim)
  imagesc(xgrid,xgrid,reshape(sum(spin(:,:,:,idim),3),ng,ng)');
  axis xy square; colorbar; caxis([-1,1])
end


%% autopower
n_row_xi=10;
fid=fopen([Path,Dir,Redshift,'_cicpower_1.bin']);
  xi=fread(fid,'real*4')';
fclose(fid);
xi=reshape(xi,n_row_xi,numel(xi)/n_row_xi)';
camb=load('CAMB/pkz0.txt');
figure(2); loglog(camb(:,1),camb(:,2).*camb(:,1).^3/(2*pi^2),'--',xi(:,2),xi(:,3))
grid on; hold on
%% plot vector decomposition power spectrum
fid=fopen('../../CUBE_spin/CUBE_v2.0/output/universe1/image1/20.000_velpower_1.bin','r');
  a=fread(fid,'real*4'); xi=reshape(a,[10,numel(a)/10]);
fclose(fid);  
figure; loglog(xi(2,:),xi(3,:),xi(2,:),xi(4,:),xi(2,:),xi(5,:),xi(2,:),xi(6,:))
legend('vv','EE','LL','RR'); title('R')
%% plot here !!!!!!!!!!
Redshift='0.000';
fid=fopen([Path,Dir,Redshift,'_spinpower_1.bin']);
  a=fread(fid,'real*4'); xi=reshape(a,[10,numel(a)/10]);
fclose(fid);  
figure; semilogx(xi(2,:),xi(4,:)./xi(3,:), xi(2,:),xi(5,:)./xi(3,:), xi(2,:),xi(6,:)./xi(3,:))
legend('E-mode','L-mode','R-mode'); title('spin power')
ax1=gca; ax1.YLim=[0.1,0.6];
%%
figure; loglog(xi(2,:),xi(3,:), xi(2,:),xi(4,:), xi(2,:),xi(5,:), xi(2,:),xi(6,:))
legend('Total power','E-mode','L-mode','R-mode'); title('spin power')

return

%%
Redshift='20.000';
fid=fopen([Path,Dir,Redshift,'_spin_2_4_1.bin'],'r');
  fx=fread(fid,ng^3,'real*4'); fx=reshape(fx,[ng,ng,ng]);
  fy=fread(fid,ng^3,'real*4'); fy=reshape(fy,[ng,ng,ng]);
  fz=fread(fid,ng^3,'real*4'); fz=reshape(fz,[ng,ng,ng]);
fclose(fid);
figure; imagesc(xgrid,xgrid,reshape(mean(fx(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar
figure; imagesc(xgrid,xgrid,reshape(mean(fy(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar
figure; imagesc(xgrid,xgrid,reshape(mean(fz(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar
%%
fid=fopen([Dir,int2str(universe),'/image1/','100.000','_spin_5_10_1.bin'],'r');
  fx=fread(fid,ng^3,'real*4'); fx=reshape(fx,[ng,ng,ng]);
  fy=fread(fid,ng^3,'real*4'); fy=reshape(fy,[ng,ng,ng]);
  fz=fread(fid,ng^3,'real*4'); fz=reshape(fz,[ng,ng,ng]);
fclose(fid);
figure; imagesc(xgrid,xgrid,reshape(mean(fx(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar
figure; imagesc(xgrid,xgrid,reshape(mean(fy(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar
figure; imagesc(xgrid,xgrid,reshape(mean(fz(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar

%% helicity
% decomposision
fid=fopen([Dir,int2str(universe),'/image1/','100.000','_spinfields_1.bin'],'r');
  fx=fread(fid,ng^3,'real*4'); fx=reshape(fx,[ng,ng,ng]);
  fy=fread(fid,ng^3,'real*4'); fy=reshape(fy,[ng,ng,ng]);
  fz=fread(fid,ng^3,'real*4'); fz=reshape(fz,[ng,ng,ng]);
  fEx=fread(fid,ng^3,'real*4'); fEx=reshape(fEx,[ng,ng,ng]);
  fEy=fread(fid,ng^3,'real*4'); fEy=reshape(fEy,[ng,ng,ng]);
  fEz=fread(fid,ng^3,'real*4'); fEz=reshape(fEz,[ng,ng,ng]);
  fLx=fread(fid,ng^3,'real*4'); fLx=reshape(fLx,[ng,ng,ng]);
  fLy=fread(fid,ng^3,'real*4'); fLy=reshape(fLy,[ng,ng,ng]);
  fLz=fread(fid,ng^3,'real*4'); fLz=reshape(fLz,[ng,ng,ng]);
  fRx=fread(fid,ng^3,'real*4'); fRx=reshape(fRx,[ng,ng,ng]);
  fRy=fread(fid,ng^3,'real*4'); fRy=reshape(fRy,[ng,ng,ng]);
  fRz=fread(fid,ng^3,'real*4'); fRz=reshape(fRz,[ng,ng,ng]);
fclose(fid);
figure; imagesc(xgrid,xgrid,reshape(mean(fx(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar
figure; imagesc(xgrid,xgrid,reshape(mean(fy(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar
figure; imagesc(xgrid,xgrid,reshape(mean(fz(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar
figure; imagesc(xgrid,xgrid,reshape(mean(fEx(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar
figure; imagesc(xgrid,xgrid,reshape(mean(fEy(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar
figure; imagesc(xgrid,xgrid,reshape(mean(fEz(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar
figure; imagesc(xgrid,xgrid,reshape(mean(fLx(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar
figure; imagesc(xgrid,xgrid,reshape(mean(fLy(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar
figure; imagesc(xgrid,xgrid,reshape(mean(fLz(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar
figure; imagesc(xgrid,xgrid,reshape(mean(fRx(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar
figure; imagesc(xgrid,xgrid,reshape(mean(fRy(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar
figure; imagesc(xgrid,xgrid,reshape(mean(fRz(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar
%%
fid=fopen([Dir,int2str(universe),'/image1/','100.000','_velofspin_1.bin'],'r');
  fx=fread(fid,ng^3,'real*4'); fx=reshape(fx,[ng,ng,ng]);
  fy=fread(fid,ng^3,'real*4'); fy=reshape(fy,[ng,ng,ng]);
  fz=fread(fid,ng^3,'real*4'); fz=reshape(fz,[ng,ng,ng]);
fclose(fid);
figure; imagesc(xgrid,xgrid,reshape(mean(fx(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar
figure; imagesc(xgrid,xgrid,reshape(mean(fy(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar
figure; imagesc(xgrid,xgrid,reshape(mean(fz(:,:,1:20),3),ng,ng)'); hold on; axis xy square; colorbar


%% RSD E-mode 2D power
fid=fopen([Dir,int2str(universe),'/image1/',redshift,'_pow2d_1.bin']);
xi=fread(fid,'real*4')';
fclose(fid);

%%  RSD 1D
n_row_xi=10;
fid=fopen([Dir,int2str(universe),'/image1/',redshift,'_cicpower_1.bin']);
xi=fread(fid,'real*4')';
fclose(fid);
xi=reshape(xi,n_row_xi,numel(xi)/n_row_xi)';
%figure;

semilogx(xi(:,2),xi(:,8))
grid on; hold on
label('$k$/($h$ Mpc$^{-1}$)','$\xi_{NL}$')
ax1=gca;
  ax1.XLim=[xi(1,2),xi(length(xi),2)];
  ax1.YLim=[0,1];

%% subplots:
figure(universe+plot_offset)

n_row_xi=10;
fid=fopen([Dir,int2str(universe),'/image1/',redshift,'_cicpower_1.bin']);
xi=fread(fid,'real*4')';
fclose(fid);
xi=reshape(xi,n_row_xi,numel(xi)/n_row_xi)';

subplot(2,2,1)
camb=load('CAMB/pkz0.txt');
%loglog(a(:,1),a(:,2).*a(:,1).^-0.04*0.88,'--'); hold on
k3=xi(:,2).^3/(2*pi^2);
Deltak_shot=Pk_shot*k3;
%loglog(xi(:,2),xi(:,3)./k3-Pk_shot,xi(:,2),xi(:,4)./k3-Pk_shot,xi(:,2),xi(:,5)./k3-Pk_shot)
loglog(camb(:,1),camb(:,2).*camb(:,1).^3/(2*pi^2),'--',...
       xi(:,2),(xi(:,3)),...
       xi(:,2),(xi(:,4)),...
       xi(:,2),Deltak_shot,'--')
grid on; hold on
label('$k$/($h$ Mpc$^{-1}$)','$\Delta^2(k)$')
legend('CLASS','$cc-$shot','linear','shot','Location','northwest')
ax1=gca;
  ax1.XLim=[xi(1,2),xi(length(xi),2)];
  ax1.YLim=[1e-3,1e3];
  
subplot(2,2,3)
semilogx(xi(:,2),xi(:,8))
grid on; hold on
label('$k$/($h$ Mpc$^{-1}$)','$\xi_{NL}$')
ax1=gca;
  ax1.XLim=[xi(1,2),xi(length(xi),2)];
  ax1.YLim=[0,1];
  
subplot(2,2,2)
delta_c=loadfield2d([Dir,int2str(universe),'/image1/',redshift,'_delta_c_proj_1.bin']);
%delta_rsd=loadfield2d([dir,int2str(universe),'/image1/',redshift,'_delta_rsd_proj_1.bin']);
imagesc(xgrid,xgrid,delta_c'); hold on
axis xy square; colorbar; colormap(1-gray); caxis([-1,5]); title('$\delta_c$')
%plot(haloinfo(1,:),haloinfo(2,:),'r.')
%% recon
camb=loadfield3d([Dir,int2str(universe),'/image1/','delta_L_1.bin']);
  figure(1); imagesc(xgrid,xgrid,reshape(mean(camb(:,:,z1:z2),3),ng,ng)'); hold on
  axis xy square; colorbar; title('$L$')

  %%
camb=loadfield3d([Dir,int2str(universe),'/image1/',redshift,'_delta_E_1.bin']);
  figure(3); imagesc(xgrid,xgrid,reshape(mean(camb(:,:,z1:z2),3),ng,ng)'); hold on
  axis xy square; colorbar; title('$R$')

%%  
camb=loadfield3d([Dir,int2str(universe),'/image1/',redshift,'_recon_1.bin']);
  figure(3); imagesc(xgrid,xgrid,reshape(mean(camb(:,:,z1:z2),3),ng,ng)'); hold on
  axis xy square; colorbar; title('$R$')

camb=loadfield3d([Dir,int2str(38),'/image1/','delta_L_1.bin']);
  figure(4); imagesc(xgrid,xgrid,reshape(mean(camb(:,:,z1:z2),3),ng,ng)'); hold on
  axis xy square; colorbar; title('$R$')  
  
camb=loadfield3d([Dir,int2str(38),'/image1/',redshift,'_delta_E_1.bin']);
  figure(5); imagesc(xgrid,xgrid,reshape(mean(camb(:,:,z1:z2),3),ng,ng)'); hold on
  axis xy square; colorbar; title('$R$')




%%
subplot(2,2,4)
delta_L=loadfield2d([Dir,int2str(universe),'/image1/delta_L_proj_1.bin']);
imagesc(xgrid,xgrid,delta_L'); axis xy square; colorbar; colormap(1-gray); caxis([-3,3])
title('$\delta_L$')

% overplot (2,2,3) for 2LPT divergence delta2
figure;
delta2=loadfield2d([Dir,int2str(universe),'/image1/0.000_recon_1.bin']);
imagesc(xgrid,xgrid,delta2'); axis xy square; colorbar; colormap(1-gray); caxis([-3,3])
title('2LPT $\delta^{(2)}$')

%% projection and halos

ninfo=39;

fid=fopen([Dir,int2str(universe),'/image1/',redshift,'_halo_1.bin']);
nhalo_tot=fread(fid,1,'integer*4')';
nhalo=fread(fid,1,'integer*4')';
den_odc=fread(fid,1,'real*4')';
haloinfo=fread(fid,[ninfo,nhalo],'real*4');
fclose(fid);

figure(100+universe+plot_offset)
%imagesc(x,x,log10(1+delta_c)'); hold on
%axis xy square; colorbar; colormap(1-gray); title('$\delta_c$')
imagesc(xgrid,xgrid,delta_c'); hold on
axis xy square; colorbar; colormap(1-gray); caxis([-1,5]); title('$\delta_c$')
plot(haloinfo(1,:),haloinfo(2,:),'r.')

% redshift space density field
figure(150+universe+plot_offset)
imagesc(xgrid,xgrid,delta_rsd'); hold on
axis xy square; colorbar; colormap(1-gray); caxis([-1,5]); title('$\delta_{\rm rsd}$')

figure(200+universe+plot_offset)
imagesc(xgrid,xgrid,delta_L'); hold on
axis xy square; colorbar; colormap(1-gray); caxis([-3,3]); title('$\delta_L$')
plot(haloinfo(28,:),haloinfo(29,:),'r.')
%plot(pfail(1,:),pfail(2,:),'g.')
return
%%
figure(300+universe+plot_offset)
imagesc(xgrid,xgrid,log10(1+delta_c)'); hold on
axis xy square; colorbar; colormap(1-gray); title('$\delta_c$')
%imagesc(x,x,delta_L'); hold on
%axis xy square; colorbar; colormap(1-gray); caxis([-3,3]); title('$\delta_c$')
for i=1:nhalo
  r_qx=sqrt((haloinfo(1,i)-haloinfo(28,i))^2+(haloinfo(2,i)-haloinfo(29,i))^2);
  if r_qx<ng/2
    plot([haloinfo(1,i),haloinfo(28,i)],[haloinfo(2,i),haloinfo(29,i)],'b'); hold on
  end
end
plot(haloinfo(1,:),haloinfo(2,:),'r.')
axis xy square;

%% halo Voronoi
den_vor=loadfield3d('~/cloud/cafproject/CUBEnu/work/utilities/Voronoi/den_interp.dat');
%den_vor=loadfield3d([dir,int2str(universe),'/image1/',redshift,'_protohaloden_vor_1.bin']);
figure
ratio=length(den_vor)/ng;
xgrid=[0.5,length(den_vor)-0.5];
imagesc(xgrid,xgrid,mean(den_vor,3)'); hold on
axis xy square; colorbar; colormap(1-gray); title('$\rho_v$')
%plot(haloinfo(1,:)*ratio,haloinfo(2,:)*ratio,'r.')
plot(haloinfo(28,:)*ratio,haloinfo(29,:)*ratio,'r.')

%% produce my own colormap
load bbr_color.mat
ntop=158; n=ntop/2;
cosmo=bbr_color(2:2:ntop,:);
for i=1:n
  cosmo(i,:)=cosmo(i,:)*i/n;
end
cosmo(:,1)=cosmo(:,1)*2;
cosmo(:,2)=cosmo(:,2).^2*2;
cosmo(:,3)=cosmo(:,3).^1.5;
cosmo=min(1,cosmo);cosmo=max(0,cosmo);

figure(100)
x=1:n;
plot(x,cosmo(1:n,1),'r',x,cosmo(1:n,2),'g',x,cosmo(1:n,3),'b')

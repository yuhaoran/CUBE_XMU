clear; SetDefault
global Redshift Redshift_i Universe Path Dir
Universe='1';
Redshift_i='100.000';
Redshift='0.000';
Path='../CUBE_v2.0/output/';
Dir=['universe',Universe,'/image1/'];
sim=get_sim_info([Path,Dir,Redshift_i,'_']);
ng=sim.nf;

%% potential
phi1=loadfield3d([Path,Dir,Redshift_i,'_phi1_1.bin']);
xgrid=[0.5,sim.box-0.5];
figure; imagesc(xgrid,xgrid,reshape(-mean(phi1(:,:,:),3),ng,ng)'); hold on
axis xy square; colorbar; title('$-\phi_1$')
%% linear density field
phi1=loadfield3d([Path,Dir,'delta_L_1.bin']);
xgrid=[0.5,sim.box-0.5];
figure; imagesc(xgrid,xgrid,reshape(mean(phi1(:,:,:),3),ng,ng)'); hold on
axis xy square; colorbar; title('$\delta_L$')
%% nonlinear density field
delta_c=loadfield3d([Path,Dir,Redshift,'_delta_c_1.bin']);
xgrid=[0.5,sim.box-0.5];
figure; imagesc(xgrid,xgrid,reshape(mean(delta_c(:,:,:),3),ng,ng)'); hold on
axis xy square; colorbar; caxis([-1,3]); title('$\delta_c$')
%colormap(1-gray);
%% FoF halo catalog
ninfo=42;
fid=fopen([Path,Dir,Redshift,'_fof_1.bin']);
  nhalo_fof_tot=fread(fid,1,'integer*4')';
  nhalo_fof=fread(fid,1,'integer*4')';
  linking_parameter=fread(fid,1,'real*4')';
  hcat=fread(fid,[ninfo,nhalo_fof],'real*4');
fclose(fid);
hcat(1:3,:)=hcat(1:3,:)*sim.box;
n_plot=10;
plot3(hcat(1,:),hcat(2,:),hcat(3,:),'k.')
plot3(hcat(1,1:n_plot),hcat(2,1:n_plot),hcat(3,1:n_plot),'ro')
%% FoF PID plot test
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
    end
    plot3(qpos(1,:),qpos(2,:),qpos(3,:),'.')
  end
fclose(fid);
for ihalo=1:nhalo_fof
  hcat(13:15,ihalo)=hcat(13:15,ihalo)/norm(hcat(13:15,ihalo)); % normalize
end
disp(mean(hcat(13:15,:).^2,2))
disp('expected dispersion')
disp(1/3/sqrt(nhalo_fof))

return
%% SO halo catalog % to be deprecated soon
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

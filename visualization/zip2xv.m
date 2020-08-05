function [ sim,xv ] = zip2xv( prefix )
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
nc=sim.nt*sim.nnt;
nf=nc*sim.ncell;
%% np
fid=fopen([prefix,'np_1.bin']);
rhoc=fread(fid,nc^3,'integer*4');
rhoc=reshape(rhoc,sim.nt,sim.nt,sim.nt,sim.nnt,sim.nnt,sim.nnt);

%% read xp
fid=fopen([prefix,'xp_1.bin']);
xzip=eval(['fread(fid,[3,sim.nplocal],''uint',int2str(8*sim.izipx),''')']);
fclose(fid);
%fid=fopen([prefix,'zip1_0.dat']);
%vzip=eval(['fread(fid,[3,sim.nplocal],''integer*',int2str(sim.izipv),''')']);
%fclose(fid);
%% xv
xv=zeros(3,sim.nplocal);
ip=0;

for itz=1:sim.nnt
for ity=1:sim.nnt
for itx=1:sim.nnt
  for k=1:sim.nt
  for j=1:sim.nt
  for i=1:sim.nt
    np=rhoc(i,j,k,itx,ity,itz);
    for l=1:np
      ip=ip+1;
      sim.xv(1:3,ip)=([itx;ity;itz]-1)*sim.nt*4 + ([i;j;k]-1)*4 + (xzip(:,ip)+0.5)/2^(sim.izipx*8)*4;
    end
  end
  end
  end
end
end
end
sim.xv=sim.xv*sim.box/nf;
sim.delta_c=loadfield3d([prefix,'delta_c_1.bin']);
return
%% read zipid
%fid=fopen([prefix,'id_0.bin']);
%sim.iczip=eval(['fread(fid,[4,sim.nplocal],''integer*',int2str(2),''')']);
%fclose(fid);

%sim.ic=nf*(sim.iczip(2:4,:)+32768.5)/65536;
%sim.disp=mod(sim.xv(1:3,:)-sim.ic+nf/2,nf)-nf/2;
%return
%if 1
%    figure
%    plot3(sim.disp(1,:),sim.disp(2,:),sim.disp(3,:),'k.')
%    grid on
%end
%% density
n=nf*2;
sim.rho_f=zeros(n+2,n+2,n+2);
mass_p=sim.mass_p;
for ip=1:sim.nplocal
    tempx=sim.xv(1:3,ip)*n/nf+0.5;
    idx1=floor(tempx)+1;
    idx2=idx1+1;
    dx1=idx1-tempx;
    dx2=1-dx1;
    sim.rho_f(idx1(1),idx1(2),idx1(3))=sim.rho_f(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*mass_p;
    sim.rho_f(idx2(1),idx1(2),idx1(3))=sim.rho_f(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*mass_p;
    sim.rho_f(idx1(1),idx2(2),idx1(3))=sim.rho_f(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*mass_p;
    sim.rho_f(idx1(1),idx1(2),idx2(3))=sim.rho_f(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*mass_p;
    sim.rho_f(idx1(1),idx2(2),idx2(3))=sim.rho_f(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*mass_p;
    sim.rho_f(idx2(1),idx1(2),idx2(3))=sim.rho_f(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*mass_p;
    sim.rho_f(idx2(1),idx2(2),idx1(3))=sim.rho_f(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*mass_p;
    sim.rho_f(idx2(1),idx2(2),idx2(3))=sim.rho_f(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*mass_p;
end
sim.rho_f(2,:,:)=sim.rho_f(2,:,:)+sim.rho_f(n+2,:,:);
sim.rho_f(n+1,:,:)=sim.rho_f(n+1,:,:)+sim.rho_f(1,:,:);
sim.rho_f(:,2,:)=sim.rho_f(:,2,:)+sim.rho_f(:,n+2,:);
sim.rho_f(:,n+1,:)=sim.rho_f(:,n+1,:)+sim.rho_f(:,1,:);
sim.rho_f(:,:,2)=sim.rho_f(:,:,2)+sim.rho_f(:,:,n+2);
sim.rho_f(:,:,n+1)=sim.rho_f(:,:,n+1)+sim.rho_f(:,:,1);

sim.rho_f=sim.rho_f(2:n+1,2:n+1,2:n+1);
if 1
    figure;
    imagesc(reshape(sum(sim.rho_f,3),n,n)'); axis xy square; colorbar
end
return
%% disp field
n=nf/2;
sim.disp=zeros(3,n+2,n+2,n+2);
for ip=1:sim.nplocal
    tempx=sim.ic(1:3,ip)*n/nf+0.5;
    idx1=floor(tempx)+1;
    idx2=idx1+1;
    dx1=idx1-tempx;
    dx2=1-dx1;
    dr=mod(sim.xv(1:3,ip)-sim.ic(1:3,ip)+nf/2,nf)-nf/2;
    sim.disp(:,idx1(1),idx1(2),idx1(3))=sim.disp(:,idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*dr;
    sim.disp(:,idx2(1),idx1(2),idx1(3))=sim.disp(:,idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*dr;
    sim.disp(:,idx1(1),idx2(2),idx1(3))=sim.disp(:,idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*dr;
    sim.disp(:,idx1(1),idx1(2),idx2(3))=sim.disp(:,idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*dr;
    sim.disp(:,idx1(1),idx2(2),idx2(3))=sim.disp(:,idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*dr;
    sim.disp(:,idx2(1),idx1(2),idx2(3))=sim.disp(:,idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*dr;
    sim.disp(:,idx2(1),idx2(2),idx1(3))=sim.disp(:,idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*dr;
    sim.disp(:,idx2(1),idx2(2),idx2(3))=sim.disp(:,idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*dr;
end
sim.disp(:,2,:,:)=sim.disp(:,2,:,:)+sim.disp(:,n+2,:,:);
sim.disp(:,n+1,:,:)=sim.disp(:,n+1,:,:)+sim.disp(:,1,:,:);
sim.disp(:,:,2,:)=sim.disp(:,:,2,:)+sim.disp(:,:,n+2,:);
sim.disp(:,:,n+1,:)=sim.disp(:,:,n+1,:)+sim.disp(:,:,1,:);
sim.disp(:,:,:,2)=sim.disp(:,:,:,2)+sim.disp(:,:,:,n+2);
sim.disp(:,:,:,n+1)=sim.disp(:,:,:,n+1)+sim.disp(:,:,:,1);

sim.disp=sim.disp(:,2:n+1,2:n+1,2:n+1);
if 1
    figure; imagesc(reshape(mean(sim.disp(1,:,:,:),4),n,n)'); axis xy square; colorbar
    figure; imagesc(reshape(mean(sim.disp(2,:,:,:),4),n,n)'); axis xy square; colorbar
    figure; imagesc(reshape(mean(sim.disp(3,:,:,:),4),n,n)'); axis xy square; colorbar    
end
sim.div=divergence(reshape(sim.disp(1,:,:,:),[n,n,n]),...
                   reshape(sim.disp(2,:,:,:),[n,n,n]),...
                   reshape(sim.disp(3,:,:,:),[n,n,n]));
[sim.curlx,sim.curly,sim.curlz,sim.cav]=curl(reshape(sim.disp(1,:,:,:),[n,n,n]),...
                   reshape(sim.disp(2,:,:,:),[n,n,n]),...
                   reshape(sim.disp(3,:,:,:),[n,n,n]));
if 1
    figure; imagesc(-reshape(sum(sim.div(:,:,:),3),n,n)'); axis xy square; colorbar
end



end

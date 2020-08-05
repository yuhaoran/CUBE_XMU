%% Reads CUBE zip format output
SetDefault
clear; clc; close all
%sim=zip2xv('/Users/haoran/cloud/CUBE_XMU/CUBE_v2.0/output/universe1/image1/0.000_');
sim=zip2xv('../CUBE_v2.0/output/universe1/image1/0.000_');
disp(sim)
%% Visualization
figure
% this plots a density projection over the 3rd direction
xgrid=[0.5,sim.box-0.5];
nf=sim.nt*sim.nnt*sim.ncell; % fine grid number per dimension
imagesc(xgrid,xgrid,reshape(mean(sim.delta_c(:,:,:),3),nf,nf)'); hold on
axis xy square; colorbar; %colormap(1-gray); 
caxis([-1,3]); title('$\delta_c$')
%% this plots all the particles in 3D
% if you have a lot (>64^3) particles this visualization might slow down your computer
h1=plot3(sim.xv(1,:),sim.xv(2,:),sim.xv(3,:),'k.'); axis square xy
ax1=gca;
ax1.XLim=[0,sim.box];
ax1.YLim=[0,sim.box];
ax1.ZLim=[0,sim.box];
h1.MarkerSize=0.01;
xlabel('$x/({\rm Mpc}/h)$')
ylabel('$y/({\rm Mpc}/h)$')
zlabel('$z/({\rm Mpc}/h)$')
view(20,30)
return

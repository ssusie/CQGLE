clear all;close all; clc;

load A1.mat
load A3.mat
load A5.mat

figure(1)
% subplot(3,1,1)
% waterfall(abs(A1))
% subplot(3,1,2)
% waterfall(abs(A3))
% subplot(3,1,3)
% waterfall(abs(A5))

L=80; n=1024;
t2=linspace(-L/2,L/2,n+1); t=t2(1:n);
p=512-100:512+100;
tp=t(p);
z=linspace(0,80,161);
%%
figure(1)
subplot(3,1,1),% waterfall(t,zt2,abs(usol(1:2:41,:))); shading interp, map=[0 0 0]; colormap(map)
waterfall(t(p),z(40:2:80), abs(A1(40:2:80,p)));
shading interp, map=[0 0 0]; colormap(map)
xlim([-10 10])
ylim([z(40) z(80)])
% surfl(t(p),z,abs(usol(:,p))), shading interp, colormap(gray)
set(gcf,'Position',[100 100 300 300]);
view(24,44)
%%
subplot(3,1,2)
waterfall(t(p),z(40:3:100), abs(A3(40:3:100,p)));
shading interp, map=[0 0 0]; colormap(map)
xlim([-10 10])
ylim([z(40) z(100)])
% surfl(t(p),z,abs(usol(:,p))), shading interp, colormap(gray)
set(gcf,'Position',[100 100 300 300]);
view(24,44)
%%
subplot(3,1,3)
waterfall(t(p),z(40:2:80), abs(A5(40:2:80,p)));
shading interp, map=[0 0 0]; colormap(map)
xlim([-10 10])
ylim([z(40) z(80)])
% surfl(t(p),z,abs(usol(:,p))), shading interp, colormap(gray)
set(gcf,'Position',[100 100 300 300]);
view(24,44)








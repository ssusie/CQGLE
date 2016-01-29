
clear all; close all; clc

n=1024; L=80;
zspan=linspace(0,100,41);
t2=linspace(-L/2,L/2,n+1); t=t2(1:n);
k=(2*pi/L)*[0:n/2-1 -n/2:-1].';

u=sech(1.0*t).'; ut=fft(u); umaster=[];

for jl=1:6

if jl==1
% parameters: three bump
tau=-0.3;
kappa=-0.05;
beta=1.45;
nu=0;
sigma=-0.1;
gamma=-0.5;

elseif jl==2
% parameters: breather 1 L=20
tau=-0.3;
kappa=-0.05;
beta=1.4;
nu=0;
sigma=-0.1;
gamma=-0.5;


elseif jl==3
% parameters: breather 1 L=20
tau=0.08;
kappa=0;
beta=0.66;
nu=-0.1;
sigma=-0.1;
gamma=-0.1;

elseif j==4
tau=0.125;
kappa=0;
beta=1;
nu=-0.6;
sigma=-0.1;
gamma=-0.1;

elseif jl==5
% parameters: fat bump L=20
tau=0.08;
kappa=-0.05;
beta=0.6;
nu=-0.1;
sigma=-0.1;
gamma=-0.1;


else jl==6
% parameters: fat bump L=20
tau=0.08;
kappa=-0.05;
beta=0.5;
nu=-0.1;
sigma=-0.1;
gamma=-0.1;
end

[z,utsol]=ode45('cqgle_rhs',zspan,ut,[],k,tau,kappa,beta,nu,sigma,gamma);

for j=1:length(z)
  usol(j,:)=ifft(utsol(j,:));
end
ut=utsol(end,:).';

umaster=[umaster; usol];

end

save master

%%
z=linspace(0,600,6*41);
surfl(t,z,abs(umaster)); shading interp, colormap(gray)
view(15,65)
set(gca,'Zlim',[0 4],'Ztick',[0 2 4],'Xlim',[-10 10],'Xtick',[-10 0 10],'Ylim',[0 300],'Ytick',[0 100 200 300],'Fontsize',[15])
xlabel('x','Fontsize',[15]), ylabel('t','Fontsize',[15])


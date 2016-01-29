clear all; close all; clc

load master

p=512-128:512+128;
tp=t(p);
z=linspace(0,300,123);

hold on


% n1=513; n2=522; n3=531;  % works with 3
% n1=513; n2=517; n3=519;    
n1=513; n2=519; n3=526; % locations from DEIM (NL_all quintic)

phi=zeros(3,1024);
phi(1,n1)=1;
phi(2,n2)=1;
phi(3,n3)=1;

nz1=11; nz2=51; nz3=90;

tt=ones(1,length(p));

load NL
A=phi*abs(nl_Psi);


b1=[(umaster(nz1,n1)); (umaster(nz1,n2)); (umaster(nz1,n3))]; %abs(umaster(nz1,n4))];
b2=[(umaster(nz2,n1)); (umaster(nz2,n2)); (umaster(nz2,n3))]; %abs(umaster(nz1,n4))];
b3=[(umaster(nz3,n1)); (umaster(nz3,n2)); (umaster(nz3,n3))]; %abs(umaster(nz1,n4))];
   
u1=(abs(b1)).^2.*b1;
u2=(abs(b2)).^2.*b2;
u3=(abs(b3)).^2.*b3;



%%
m=24;
cvx_begin;
variable x2(m); 
   minimize( norm(x2,1) ); 
   subject to
    A*x2 == abs(u1);
cvx_end;

figure(1)
subplot(3,1,1)
plot(1,x2(1),'mo','Linewidth',[3]), hold on
plot(2,x2(2),'ko','Linewidth',[3])
plot(3:8,x2(3:8),'bo','Linewidth',[3])
plot(9:22,x2(9:22),'co','Linewidth',[3])
plot(23,x2(23),'ro','Linewidth',[3])
plot(24,x2(24),'go','Linewidth',[3])
title('nl for beta1 ')
clear x2
% break

%%

m=24;
cvx_begin;
variable x2(m); 
   minimize( norm(x2,1) ); 
   subject to
    A*x2 ==abs(u2);
cvx_end;

% figure(2)
subplot(3,1,2)
plot(1,x2(1),'mo','Linewidth',[3]), hold on
plot(2,x2(2),'ko','Linewidth',[3])
plot(3:8,x2(3:8),'bo','Linewidth',[3])
plot(9:22,x2(9:22),'co','Linewidth',[3])
plot(23,x2(23),'ro','Linewidth',[3])
plot(24,x2(24),'go','Linewidth',[3])
title('nl for beta3 ')
clear x2

%%

m=24;
cvx_begin;
variable x3(m); 
   minimize( norm(x3,1) ); 
   subject to
    A*x3 == abs(u3);
cvx_end;

% figure(3)
subplot(3,1,3)
plot(1,x3(1),'mo','Linewidth',[3]), hold on
plot(2,x3(2),'ko','Linewidth',[3])
plot(3:8,x3(3:8),'bo','Linewidth',[3])
plot(9:22,x3(9:22),'co','Linewidth',[3])
plot(23,x3(23),'ro','Linewidth',[3])
plot(24,x3(24),'go','Linewidth',[3])
title('nl for beta5 ')






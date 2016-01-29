clear all; close all; clc

load master

p=512-128:512+128;
tp=t(p);

z=linspace(0,300,123);

hold on


   
n1=513; n2=519; n3=526; % locations from DEIM (average)

phi=zeros(3,1024);
phi(1,n1)=1;
phi(2,n2)=1;
phi(3,n3)=1;

nz1=11; nz2=51; nz3=90;

tt=ones(1,length(p));


load nl2_Psi
load nl4_Psi


% psi2=abs(NL_Psi);

A2=phi*abs(nl2_Psi);
A4=phi*abs(nl4_Psi);


b1=[abs(umaster(nz1,n1)); abs(umaster(nz1,n2)); abs(umaster(nz1,n3))]; %abs(umaster(nz1,n4))];
b2=[abs(umaster(nz2,n1)); abs(umaster(nz2,n2)); abs(umaster(nz2,n3))]; %abs(umaster(nz1,n4))];
b3=[abs(umaster(nz3,n1)); abs(umaster(nz3,n2)); abs(umaster(nz3,n3))]; %abs(umaster(nz1,n4))];

u12=(abs(b1)).^2.*b1;
u14=(abs(b1)).^4.*b1;

u22=(abs(b2)).^2.*b2;
u24=(abs(b2)).^4.*b2;

u32=(abs(b3)).^2.*b3;
u34=(abs(b3)).^4.*b3;


%%
m=24;
cvx_begin;
variable x2(m); 
   minimize( norm(x2,1) ); 
   subject to
    A2*x2 == u12;
cvx_end;




m=24;
cvx_begin;
variable x4(m); 
   minimize( norm(x4,1) ); 
   subject to
    A4*x4 == u14;
cvx_end;
x2=abs(x2);
figure(1)
subplot(3,1,1)
plot(1,x2(1),'mo','Linewidth',[3]), hold on
plot(2,x2(2),'ko','Linewidth',[3])
plot(3:8,x2(3:8),'bo','Linewidth',[3])
plot(9:22,x2(9:22),'co','Linewidth',[3])
plot(23,x2(23),'ro','Linewidth',[3])
plot(24,x2(24),'go','Linewidth',[3])
title('cubic')
ylim([-25 max(x2)+50])
set(gcf,'Position',[100 100 300 300]); 

x4=abs(x4);
figure(2)
subplot(3,1,1)
plot(1,x4(1),'mo','Linewidth',[3]), hold on
plot(2,x4(2),'ko','Linewidth',[3])
plot(3:8,x4(3:8),'bo','Linewidth',[3])
plot(9:22,x4(9:22),'co','Linewidth',[3])
plot(23,x4(23),'ro','Linewidth',[3])
plot(24,x4(24),'go','Linewidth',[3])
ylim([-300 max(x4)+500])
title('quintic ')
set(gcf,'Position',[100 100 300 300]); 

% break
%%
%%%%%%%% compute errors %%%%%%

er2(1)=norm(u12-A2(:,1)*x2(1))/norm(u12);
er2(2)=norm(u12-A2(:,2)*x2(2))/norm(u12);
er2(3)=norm(u12-A2(:,3:8)*x2(3:8))/norm(u12);
er2(4)=norm(u12-A2(:,9:22)*x2(9:22))/norm(u12);
er2(5)=norm(u12-A2(:,23)*x2(23))/norm(u12);
er2(6)=norm(u12-A2(:,24)*x2(24))/norm(u12);


er4(1)=norm(u14-A4(:,1)*x4(1))/norm(u14);
er4(2)=norm(u14-A4(:,2)*x4(2))/norm(u14);
er4(3)=norm(u14-A4(:,3:8)*x4(3:8))/norm(u14);
er4(4)=norm(u14-A4(:,9:22)*x4(9:22))/norm(u14);
er4(5)=norm(u14-A4(:,23)*x4(23))/norm(u14);
er4(6)=norm(u14-A4(:,24)*x4(24))/norm(u14);

figure(11)
subplot(3,1,1)
plot(er2,'ko','Linewidth',[3])
hold on
plot(er4,'rv','Linewidth',[3])
% title('E')
legend('E cubic', 'E quintic')
set(gcf,'Position',[100 100 300 300]); 
% break
%%

m=24;
cvx_begin;
variable x2(m); 
   minimize( norm(x2,1) ); 
   subject to
    A2*x2 == u22;
cvx_end;



m=24;
cvx_begin;
variable x4(m); 
   minimize( norm(x4,1) ); 
   subject to
    A4*x4 == u24;
cvx_end;

figure(1)
subplot(3,1,2)
plot(1,x2(1),'mo','Linewidth',[3]), hold on
plot(2,x2(2),'ko','Linewidth',[3])
plot(3:8,x2(3:8),'bo','Linewidth',[3])
plot(9:22,x2(9:22),'co','Linewidth',[3])
plot(23,x2(23),'ro','Linewidth',[3])
plot(24,x2(24),'go','Linewidth',[3])
ylim([-5 max(x2)+10])

figure(2)
subplot(3,1,2)
plot(1,x4(1),'mo','Linewidth',[3]), hold on
plot(2,x4(2),'ko','Linewidth',[3])
plot(3:8,x4(3:8),'bo','Linewidth',[3])
plot(9:22,x4(9:22),'co','Linewidth',[3])
plot(23,x4(23),'ro','Linewidth',[3])
plot(24,x4(24),'go','Linewidth',[3])
ylim([-50 max(x4)+100])

%%
er2(1)=norm(u22-A2(:,1)*x2(1))/norm(u22);
er2(2)=norm(u22-A2(:,2)*x2(2))/norm(u22);
er2(3)=norm(u22-A2(:,3:8)*x2(3:8))/norm(u22);
er2(4)=norm(u22-A2(:,9:22)*x2(9:22))/norm(u22);
er2(5)=norm(u22-A2(:,23)*x2(23))/norm(u22);
er2(6)=norm(u22-A2(:,24)*x2(24))/norm(u22);


er4(1)=norm(u24-A4(:,1)*x4(1))/norm(u24);
er4(2)=norm(u24-A4(:,2)*x4(2))/norm(u24);
er4(3)=norm(u24-A4(:,3:8)*x4(3:8))/norm(u24);
er4(4)=norm(u24-A4(:,9:22)*x4(9:22))/norm(u24);
er4(5)=norm(u24-A4(:,23)*x4(23))/norm(u24);
er4(6)=norm(u24-A4(:,24)*x4(24))/norm(u24);

figure(11)
subplot(3,1,2)
plot(er2,'ko','Linewidth',[3])
hold on
plot(er4,'rv','Linewidth',[3])


%%

m=24;
cvx_begin;
variable x2(m); 
   minimize( norm(x2,1) ); 
   subject to
    A2*x2 == u32;
cvx_end;


m=24;
cvx_begin;
variable x4(m); 
   minimize( norm(x4,1) ); 
   subject to
    A4*x4 == u34;
cvx_end;

figure(1)
subplot(3,1,3)
plot(1,x2(1),'mo','Linewidth',[3]), hold on
plot(2,x2(2),'ko','Linewidth',[3])
plot(3:8,x2(3:8),'bo','Linewidth',[3])
plot(9:22,x2(9:22),'co','Linewidth',[3])
plot(23,x2(23),'ro','Linewidth',[3])
plot(24,x2(24),'go','Linewidth',[3])
ylim([-10 max(x2)+20])

figure(2)
subplot(3,1,3)
plot(1,x4(1),'mo','Linewidth',[3]), hold on
plot(2,x4(2),'ko','Linewidth',[3])
plot(3:8,x4(3:8),'bo','Linewidth',[3])
plot(9:22,x4(9:22),'co','Linewidth',[3])
plot(23,x4(23),'ro','Linewidth',[3])
plot(24,x4(24),'go','Linewidth',[3])
ylim([-50 max(x4)+100])


%%
er2(1)=norm(u32-A2(:,1)*x2(1))/norm(u32);
er2(2)=norm(u32-A2(:,2)*x2(2))/norm(u32);
er2(3)=norm(u32-A2(:,3:8)*x2(3:8))/norm(u32);
er2(4)=norm(u32-A2(:,9:22)*x2(9:22))/norm(u32);
er2(5)=norm(u32-A2(:,23)*x2(23))/norm(u32);
er2(6)=norm(u32-A2(:,24)*x2(24))/norm(u32);


er4(1)=norm(u34-A4(:,1)*x4(1))/norm(u34);
er4(2)=norm(u34-A4(:,2)*x4(2))/norm(u34);
er4(3)=norm(u34-A4(:,3:8)*x4(3:8))/norm(u34);
er4(4)=norm(u34-A4(:,9:22)*x4(9:22))/norm(u34);
er4(5)=norm(u34-A4(:,23)*x4(23))/norm(u34);
er4(6)=norm(u34-A4(:,24)*x4(24))/norm(u34);


figure(11)
subplot(3,1,3)
plot(er2,'ko','Linewidth',[3])
hold on
plot(er4,'rv','Linewidth',[3])












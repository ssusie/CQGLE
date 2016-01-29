clear all; close all; clc

load master

p=512-128:512+128;
tp=t(p);


%%%%% locations from DEIM (average cubic and quintic)
% n1=513; n2=518; n3=528; %3 sensors 0 5 15 
n1=513; n2=519; n3=526;%locs from full nl
%  n1=513; n2=518; n3=526; n4=530;  %4 sensors 0 5 13 17
%  n1=513; n2=518; n3=523; n4=528; n5=531;  % 5 sensors 0 5 10 15 18
%  n1=513; n2=518; n3=521; n4=525; n5=528; n6=531;  % 6 sensors 
%0 5 8 12 15 18

phi_n=zeros(3,1024);
phi_n(1,n1)=1;
phi_n(2,n2)=1;
phi_n(3,n3)=1;
% phi_n(4,n4)=1;
% phi_n(5,n5)=1;
% phi_n(6,n6)=1;

nz1=11; nz2=51; nz3=90;

tt=ones(1,length(p));

load nl2_Psi
load nl4_Psi


noise=0.2;

A2=phi_n*abs(nl2_Psi);
A4=phi_n*abs(nl4_Psi);

N2=zeros(1,6);
N4=zeros(1,6);
K2=zeros(1,6);
K4=zeros(1,6);
J2=zeros(1,6);
J4=zeros(1,6);
nloop=400;
m=24;

for j=1:nloop

    b1_n=[abs(umaster(nz1,n1)); abs(umaster(nz1,n2)); abs(umaster(nz1,n3))]+noise*randn(3,1);
    b2_n=[abs(umaster(nz2,n1)); abs(umaster(nz2,n2)); abs(umaster(nz2,n3))]+noise*randn(3,1);
    b3_n=[abs(umaster(nz3,n1)); abs(umaster(nz3,n2)); abs(umaster(nz3,n3))]+noise*randn(3,1);

    u12=(abs(b1_n)).^2.*b1_n;
    u14=(abs(b1_n)).^4.*b1_n;

    u22=(abs(b2_n)).^2.*b2_n;
    u24=(abs(b2_n)).^4.*b2_n;

    u32=(abs(b3_n)).^2.*b3_n;
    u34=(abs(b3_n)).^4.*b3_n;



    cvx_begin;
    variable x2(m); 
       minimize( norm(x2,1) ); 
       subject to
        A2*x2 == u12;
    cvx_end;


    cvx_begin;
    variable x4(m); 
       minimize( norm(x4,1) ); 
       subject to
        A4*x4 == u14;
    cvx_end;


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
    
    [M2, I2]=min(abs(er2));
    [M4, I4]=min(abs(er4));
      
    N2(I2)=N2(I2)+1;
    N4(I4)=N4(I4)+1;
      
     clear x2 x4 M2 M4 I2 I4 er2 er4;
      
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


    %%%%%%%% compute errors %%%%%%

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

      
      [M2, I2]=min(abs(er2));
      [M4, I4]=min(abs(er4));
      
      K2(I2)=K2(I2)+1;
      K4(I4)=K4(I4)+1;
      
      clear x2 x4 M2 M4 I2 I4  er2 er4;
      
          cvx_begin;
        variable x2(m); 
           minimize( norm(x2,1) ); 
           subject to
            A2*x2 == u32;
        cvx_end;


        cvx_begin;
        variable x4(m); 
           minimize( norm(x4,1) ); 
           subject to
            A4*x4 == u34;
        cvx_end;


    %%%%%%%% compute errors %%%%%%

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

      
      [M2, I2]=min(abs(er2));
      [M4, I4]=min(abs(er4));
      
      J2(I2)=J2(I2)+1;
      J4(I4)=J4(I4)+1;

      clear  x2 x4 M2 M4 I2 I4 er2 er4;
       

end


nn2=N2*100/nloop;
nn4=N4*100/nloop;
k2=K2*100/nloop;
k4=K4*100/nloop;
j2=J2*100/nloop;
j4=J4*100/nloop;



save sens3_400 nn2 nn4 k2 k4 j2 j4 

%%
figure(1)
B1=[nn2; nn4];
B3=[k2; k4];
B5=[j2;j4];
figure(1) 
subplot(1,3,1)
bar(B1')
ylim([0 100])
xlim([0.5 6.5])
subplot(1,3,2)
bar(B3')
ylim([0 100])
xlim([0.5 6.5])
subplot(1,3,3)
bar(B5')
ylim([0 100])
xlim([0.5 6.5])

% 
% % % set(gcf,'Position',[300 300 100 100]); 
% % set(gcf,'Position',[100 100 700 150])
% 
% 

clear all; close all; clc

load master

p=512-128:512+128;
tp=t(p);

n1=513; n2=519; n3=526;  % DEIM NL ALL


phi_n=zeros(3,1024);
phi_n(1,n1)=1;
phi_n(2,n2)=1;
phi_n(3,n3)=1;

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
nloop=10;
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


    cvx_begin quiet ;
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
      
      cvx_begin quiet;
    variable x2(m); 
       minimize( norm(x2,1) ); 
       subject to
        A2*x2 == u22;
    cvx_end;


    cvx_begin quiet;
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
      
      clear x2 x4 M2 M4 I2 I4 er2 er4;
      
          cvx_begin quiet;
        variable x2(m); 
           minimize( norm(x2,1) ); 
           subject to
            A2*x2 == u32;
        cvx_end;


        cvx_begin quiet;
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

      clear er2 er4 x2 x4 M2 M4 I2 I4;
       

end


% n2=[N2(1) N2(2) sum(N2(3:8)) sum(N2(9:22)) N2(23) N2(24)];
% n4=[N4(1) N4(2) sum(N4(3:8)) sum(N4(9:22)) N4(23) N4(24)]; 
% 
% k2=[K2(1) K2(2) sum(K2(3:8)) sum(K2(9:22)) K2(23) K2(24)];
% k4=[K4(1) K4(2) sum(K4(3:8)) sum(K4(9:22)) K4(23) K4(24)];
% 
% j2=[J2(1) J2(2) sum(J2(3:8)) sum(J2(9:22)) J2(23) J2(24)];
% j4=[J4(1) J4(2) sum(J4(3:8)) sum(J4(9:22)) J4(23) J4(24)];

nn2=N2*100/nloop;
nn4=N4*100/nloop;
k2=K2*100/nloop;
k4=K4*100/nloop;
j2=J2*100/nloop;
j4=J4*100/nloop;

%%
figure(1)
subplot(3,1,1)
bar(nn2)
title('cubic library')
ylim([0 100])
subplot(3,1,2)
bar(k2)
ylim([0 100])
subplot(3,1,3)
bar(j2)
ylim([0 100])
xlabel('regime')
ylabel('accuracy')
set(gcf,'Position',[100 100 300 300]); 


%%
figure(2)
subplot(3,1,1)
bar(nn4)
title('quartic library')
ylim([0 100])
subplot(3,1,2)
bar(k4)
ylim([0 100])
subplot(3,1,3)
bar(j4)
ylim([0 100])
xlabel('regime')
ylabel('accuracy')
set(gcf,'Position',[100 100 300 300]); 
%%

Y2=[nn2' k2' j2' ];
Y4=[nn4' k4' j4'];

figure(3)
subplot(2,1,1)
bar3(Y2)
xlim([0.5 3.5])
ylim([0.5 6.5])
% view(-60, 14)
title('cubic library')
subplot(2,1,2)
bar3(Y4)
xlim([0.5 3.5])
ylim([0.5 6.5])
% view(-60, 14)
title('quintic library')
set(gcf,'Position',[100 100 300 300]); 






break







% % break
% 
% M2=mean(er2,2);
% M4=mean(er4,2);
% 
% figure(2)
% bar(M2)
% figure(2)
% bar(M4)
% 




%%

cvx_begin;
variable x2n(m); 
   minimize( norm(x2n,1) ); 
   subject to
    A*x2n == b2_n;
cvx_end;








%%
cvx_begin;
variable x3n(m); 
   minimize( norm(x3n,1) ); 
   subject to
    A*x3n == b3_n;
cvx_end;




break
%  mode reconstruction

z=[25:2:125]; n=1024; L=80;
k=(2*pi/L)*[0:n/2-1 -n/2:-1].';

psiA=psiL(:,1);
for j=1:1
  psiA2x(:,j)=ifft(-(k.^2).*fft(psiA(:,j)));
  psiA4x(:,j)=ifft( (k.^4).*fft(psiA(:,j)));
end
lhs=inv((psiA.')*psiA);

A_test1=phi*psiA;
a_ic=A_test1\b1;

tau=-0.3;
ka=-0.05;
be=1.45;
nu=0;
si=-0.1;
ga=-0.5;


[z,a_sol]=ode45('a_rhs',z,a_ic,[],psiA,psiA2x,psiA4x,lhs,tau,ka,be,nu,si,ga);
for j=1:length(z)
    u_sol(j,:)=( psiA*(a_sol(j,:).') ).';
end
u_lr_master=u_sol;


z=[125:2:225]; 
psiA=psiL(:,3:8);

for j=1:6
  psiA2x(:,j)=ifft(-(k.^2).*fft(psiA(:,j)));
  psiA4x(:,j)=ifft( (k.^4).*fft(psiA(:,j)));
end
lhs=inv((psiA.')*psiA);

A_test2=phi*psiA;
a_ic=A_test2\b2;

tau=0.08;
ka=0;
be=0.66;
nu=-0.1;
si=-0.1;
ga=-0.1;

[z,a_sol]=ode45('a_rhs',z,a_ic,[],psiA,psiA2x,psiA4x,lhs,tau,ka,be,nu,si,ga);
for j=1:length(z)
    u_sol(j,:)=( psiA*(a_sol(j,:).') ).';
end
u_lr_master=[u_lr_master; u_sol];

clear psiA2x
clear psiA4x

z=[225:2:325]; 
psiA=psiL(:,23);

for j=1:1
  psiA2x(:,j)=ifft(-(k.^2).*fft(psiA(:,j)));
  psiA4x(:,j)=ifft( (k.^4).*fft(psiA(:,j)));
end
lhs=inv((psiA.')*psiA);

A_test3=phi*psiA;
a_ic=A_test3\b3;

tau=0.08;
ka=-0.05;
be=0.6;
nu=-0.1;
si=-0.1;
ga=-0.1;

[z,a_sol]=ode45('a_rhs',z,a_ic,[],psiA,psiA2x,psiA4x,lhs,tau,ka,be,nu,si,ga);
for j=1:length(z)
    u_sol(j,:)=( psiA*(a_sol(j,:).') ).';
end
u_lr_master=[u_lr_master; u_sol];

figure(1)
z=linspace(25,325,153);
subplot(2,2,2), surfl(tp,z,abs(u_lr_master(:,p))); shading interp, colormap(gray)
view(15,65)
set(gca,'Zlim',[0 4],'Ztick',[0 2 4],'Xlim',[-10 10],'Xtick',[-10 0 10],'Ylim',[0 325],'Ytick',[0 100 200 300],'Fontsize',[15])
xlabel('x','Fontsize',[15]), ylabel('t','Fontsize',[15])
text(-16.5,0,3,'|U|','Fontsize',[15])

break


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% identify & reconstruct with noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_lr_mastern=[];
for jn=1:3;

    if jn==1
        z=[25:2:125]; 
        E(1)=abs(x1n(1)).^2;
        E(2)=abs(x1n(2)).^2;
        E(3)=sum(abs(x1n(3:8)).^2);
        E(4)=sum(abs(x1n(9:22)).^2);
        E(5)=abs(x1n(23)).^2;
        E(6)=abs(x1n(24)).^2;
        b_n=b1_n;
    elseif jn==2
        z=[125:2:225]; 
        E(1)=abs(x2n(1)).^2;
        E(2)=abs(x2n(2)).^2;
        E(3)=sum(abs(x2n(3:8)).^2);
        E(4)=sum(abs(x2n(9:22)).^2);
        E(5)=abs(x2n(23)).^2;
        E(6)=abs(x2n(24)).^2;
        b_n=b2_n;
    else
        z=[225:2:325]; 
        E(1)=abs(x3n(1)).^2;
        E(2)=abs(x3n(2)).^2;
        E(3)=sum(abs(x3n(3:8)).^2);
        E(4)=sum(abs(x3n(9:22)).^2);
        E(5)=abs(x3n(23)).^2;
        E(6)=abs(x3n(24)).^2;
        b_n=b3_n;
    end

[Emax,jloop]=max(E);  % identify regime
clear psiA2x
clear psiA4x
clear psiA

%jloop

if jloop==1
tau=-0.3;
ka=-0.05;
be=1.45;
nu=0;
si=-0.1;
ga=-0.5;

psiA=psiL(:,1);
for j=1:1
  psiA2x(:,j)=ifft(-(k.^2).*fft(psiA(:,j)));
  psiA4x(:,j)=ifft( (k.^4).*fft(psiA(:,j)));
end
lhs=inv((psiA.')*psiA);


elseif jloop==2
tau=-0.3;
ka=-0.05;
be=1.4;
nu=0;
si=-0.1;
ga=-0.5;

psiA=psiL(:,2);
for j=1:1
  psiA2x(:,j)=ifft(-(k.^2).*fft(psiA(:,j)));
  psiA4x(:,j)=ifft( (k.^4).*fft(psiA(:,j)));
end
lhs=inv((psiA.')*psiA);


elseif jloop==3
tau=0.08;
ka=0;
be=0.66;
nu=-0.1;
si=-0.1;
ga=-0.1;

psiA=psiL(:,3:8);
for j=1:6
  psiA2x(:,j)=ifft(-(k.^2).*fft(psiA(:,j)));
  psiA4x(:,j)=ifft( (k.^4).*fft(psiA(:,j)));
end
lhs=inv((psiA.')*psiA);


elseif jloop==4
tau=0.125;
ka=0;
be=1;
nu=-0.6;
si=-0.1;
ga=-0.1;

psiA=psiL(:,9:22);
for j=1:14
  psiA2x(:,j)=ifft(-(k.^2).*fft(psiA(:,j)));
  psiA4x(:,j)=ifft( (k.^4).*fft(psiA(:,j)));
end
lhs=inv((psiA.')*psiA);


elseif jloop==5
tau=0.08;
ka=-0.05;
be=0.6;
nu=-0.1;
si=-0.1;
ga=-0.1;

psiA=psiL(:,23);
for j=1:1
  psiA2x(:,j)=ifft(-(k.^2).*fft(psiA(:,j)));
  psiA4x(:,j)=ifft( (k.^4).*fft(psiA(:,j)));
end
lhs=inv((psiA.')*psiA);


elseif jloop==6
tau=0.08;
ka=-0.05;
be=0.5;
nu=-0.1;
si=-0.1;
ga=-0.1;

psiA=psiL(:,24);

for j=1:1
  psiA2x(:,j)=ifft(-(k.^2).*fft(psiA(:,j)));
  psiA4x(:,j)=ifft( (k.^4).*fft(psiA(:,j)));
end
lhs=inv((psiA.')*psiA);


end

A_test1=phi*psiA;
a_icn=A_test1\b_n;
    
[z,a_soln]=ode45('a_rhs',z,a_icn,[],psiA,psiA2x,psiA4x,lhs,tau,ka,be,nu,si,ga);
for j=1:length(z)
    u_soln(j,:)=( psiA*(a_soln(j,:).') ).';
end
u_lr_mastern=[u_lr_mastern; u_soln];
      
end    


figure(1)
z=linspace(25,325,153);
subplot(2,2,3), surfl(tp,z,abs(u_lr_mastern(:,p))); shading interp, colormap(gray)
view(15,65)
set(gca,'Zlim',[0 4],'Ztick',[0 2 4],'Xlim',[-10 10],'Xtick',[-10 0 10],'Ylim',[0 325],'Ytick',[0 100 200 300],'Fontsize',[15])
xlabel('x','Fontsize',[15]), ylabel('t','Fontsize',[15])
text(-16.5,0,3,'|U|','Fontsize',[15])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% identify & reconstruct with noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear u_lr_mastern;
u_lr_mastern=[];
for jn=1:3;

    if jn==1
        z=[25:2:125]; 
        E(1)=abs(x1_5(1)).^2;
        E(2)=abs(x1_5(2)).^2;
        E(3)=sum(abs(x1_5(3:8)).^2);
        E(4)=sum(abs(x1_5(9:22)).^2);
        E(5)=abs(x1_5(23)).^2;
        E(6)=abs(x1_5(24)).^2;
        b_n=b1_5;
    elseif jn==2
        z=[125:2:225]; 
        E(1)=abs(x2_5(1)).^2;
        E(2)=abs(x2_5(2)).^2;
        E(3)=sum(abs(x2_5(3:8)).^2);
        E(4)=sum(abs(x2_5(9:22)).^2);
        E(5)=abs(x2_5(23)).^2;
        E(6)=abs(x2_5(24)).^2;
        b_n=b2_5;
    else
        z=[225:2:325]; 
        E(1)=abs(x3_5(1)).^2;
        E(2)=abs(x3_5(2)).^2;
        E(3)=sum(abs(x3_5(3:8)).^2);
        E(4)=sum(abs(x3_5(9:22)).^2);
        E(5)=abs(x3_5(23)).^2;
        E(6)=abs(x3_5(24)).^2;
        b_n=b3_5;
    end

[Emax,jloop]=max(E);  % identify regime
clear psiA2x 
clear psiA4x
clear psiA
%jloop

if jloop==1
tau=-0.3;
ka=-0.05;
be=1.45;
nu=0;
si=-0.1;
ga=-0.5;

psiA=psiL(:,1);
for j=1:1
  psiA2x(:,j)=ifft(-(k.^2).*fft(psiA(:,j)));
  psiA4x(:,j)=ifft( (k.^4).*fft(psiA(:,j)));
end
lhs=inv((psiA.')*psiA);


elseif jloop==2
tau=-0.3;
ka=-0.05;
be=1.4;
nu=0;
si=-0.1;
ga=-0.5;

psiA=psiL(:,2);
for j=1:1
  psiA2x(:,j)=ifft(-(k.^2).*fft(psiA(:,j)));
  psiA4x(:,j)=ifft( (k.^4).*fft(psiA(:,j)));
end
lhs=inv((psiA.')*psiA);


elseif jloop==3
tau=0.08;
ka=0;
be=0.66;
nu=-0.1;
si=-0.1;
ga=-0.1;

psiA=psiL(:,3:8);
for j=1:6
  psiA2x(:,j)=ifft(-(k.^2).*fft(psiA(:,j)));
  psiA4x(:,j)=ifft( (k.^4).*fft(psiA(:,j)));
end
lhs=inv((psiA.')*psiA);


elseif jloop==4
tau=0.125;
ka=0;
be=1;
nu=-0.6;
si=-0.1;
ga=-0.1;

psiA=psiL(:,9:22);
for j=1:14
  psiA2x(:,j)=ifft(-(k.^2).*fft(psiA(:,j)));
  psiA4x(:,j)=ifft( (k.^4).*fft(psiA(:,j)));
end
lhs=inv((psiA.')*psiA);


elseif jloop==5
tau=0.08;
ka=-0.05;
be=0.6;
nu=-0.1;
si=-0.1;
ga=-0.1;

psiA=psiL(:,23);
for j=1:1
  psiA2x(:,j)=ifft(-(k.^2).*fft(psiA(:,j)));
  psiA4x(:,j)=ifft( (k.^4).*fft(psiA(:,j)));
end
lhs=inv((psiA.')*psiA);


elseif jloop==6
tau=0.08;
ka=-0.05;
be=0.5;
nu=-0.1;
si=-0.1;
ga=-0.1;

psiA=psiL(:,24);
for j=1:1
  psiA2x(:,j)=ifft(-(k.^2).*fft(psiA(:,j)));
  psiA4x(:,j)=ifft( (k.^4).*fft(psiA(:,j)));
end
lhs=inv((psiA.')*psiA);


end

A_test1=phi_5*psiA;
a_icn=(A_test1\b_n);
    
[z,a_soln]=ode45('a_rhs',z,a_icn,[],psiA,psiA2x,psiA4x,lhs,tau,ka,be,nu,si,ga);
for j=1:length(z)
    u_soln(j,:)=( psiA*(a_soln(j,:).') ).';
end
u_lr_mastern=[u_lr_mastern; u_soln];
      
end    


figure(1)
z=linspace(25,325,153);
subplot(2,2,4), surfl(tp,z,abs(u_lr_mastern(:,p))); shading interp, colormap(gray)
view(15,65)
set(gca,'Zlim',[0 4],'Ztick',[0 2 4],'Xlim',[-10 10],'Xtick',[-10 0 10],'Ylim',[0 325],'Ytick',[0 100 200 300],'Fontsize',[15])
xlabel('x','Fontsize',[15]), ylabel('t','Fontsize',[15])
text(-16.5,0,3,'|U|','Fontsize',[15])




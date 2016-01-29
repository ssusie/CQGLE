clear all; close all; clc;

load psi_L.mat
load nl2_Psi.mat
load nl4_Psi.mat

L=80; n=1024;
t2=linspace(-L/2,L/2,n+1); t=t2(1:n);
sli=ones(1,1024);

p=512-50:512+50;
tp=t(p);
slip=sli(p);

% break
%%
figure(1)

subplot(3,1,1)
plot3(tp,1*slip,abs(psiL(p,1)),'m','Linewidth',[1]), hold on
plot3(tp,3*slip,abs(psiL(p,2)),'k','Linewidth',[1])


plot3(tp,5*slip,abs(psiL(p,3)),'b','Linewidth',[1])
plot3(tp,7*slip,abs(psiL(p,4)),'b','Linewidth',[1])
plot3(tp,9*slip,abs(psiL(p,5)),'b','Linewidth',[1])
plot3(tp,11*slip,abs(psiL(p,6)),'b','Linewidth',[1])
plot3(tp,13*slip,abs(psiL(p,7)),'b','Linewidth',[1])
plot3(tp,15*slip,abs(psiL(p,8)),'b','Linewidth',[1])


plot3(tp,18*slip,abs(psiL(p,23)),'r','Linewidth',[1])
plot3(tp,20*slip,abs(psiL(p,24)),'g','Linewidth',[1])
grid on
view(98, 24)
xlim([-4 4])
title('modes of the full system')

%  break

%%

subplot(3,1,2)

plot3(tp,1*slip,abs(nl2_Psi(p,1)),'m','Linewidth',[1]), hold on
plot3(tp,3*slip,abs(nl2_Psi(p,2)),'k','Linewidth',[1])


plot3(tp,5*slip,abs(nl2_Psi(p,3)),'b','Linewidth',[1])
plot3(tp,7*slip,abs(nl2_Psi(p,4)),'b','Linewidth',[1])
plot3(tp,9*slip,abs(nl2_Psi(p,5)),'b','Linewidth',[1])
plot3(tp,11*slip,abs(nl2_Psi(p,6)),'b','Linewidth',[1])
plot3(tp,13*slip,abs(nl2_Psi(p,7)),'b','Linewidth',[1])
plot3(tp,15*slip,abs(nl2_Psi(p,8)),'b','Linewidth',[1])

plot3(tp,18*slip,abs(nl2_Psi(p,23)),'r','Linewidth',[1])
plot3(tp,20*slip,abs(nl2_Psi(p,24)),'g','Linewidth',[1])
grid on
view(98,24)
xlim([-4 4])
title('nonlinear cubic modes')


%%

subplot(3,1,3)

plot3(tp,1*slip,abs(nl4_Psi(p,1)),'m','Linewidth',[1]), hold on
plot3(tp,3*slip,abs(nl4_Psi(p,2)),'k','Linewidth',[1])


plot3(tp,5*slip,abs(nl4_Psi(p,3)),'b','Linewidth',[1])
plot3(tp,7*slip,abs(nl4_Psi(p,4)),'b','Linewidth',[1])
plot3(tp,9*slip,abs(nl4_Psi(p,5)),'b','Linewidth',[1])
plot3(tp,11*slip,abs(nl4_Psi(p,6)),'b','Linewidth',[1])
plot3(tp,13*slip,abs(nl4_Psi(p,7)),'b','Linewidth',[1])
plot3(tp,15*slip,abs(nl4_Psi(p,8)),'b','Linewidth',[1])


plot3(tp,18*slip,abs(nl4_Psi(p,23)),'r','Linewidth',[1])
plot3(tp,20*slip,abs(nl4_Psi(p,24)),'g','Linewidth',[1])
grid on

view(98, 24)
xlim([-4 4])
zlim([0 0.4])
title('nonlinear quintic modes')
set(gcf,'Position',[100 100 300 300]); 

 break
%%


figure(2)

subplot(3,1,1)
plot(t,abs(psiL(:,1)), 'm','Linewidth',[2])
hold on
plot(t,abs(psiL(:,2)), 'k','Linewidth',[2])
plot(t,abs(psiL(:,3)), 'b','Linewidth',[2])
plot(t,abs(psiL(:,23)), 'r','Linewidth',[2])
plot(t,abs(psiL(:,24)), 'g','Linewidth',[2])
xlim([0 5])
ylim([0 0.5])
title('nonlinear modes')


subplot(3,1,2)
plot(t,abs(nl2_Psi(:,1)), 'm','Linewidth',[2])
hold on
plot(t,abs(nl2_Psi(:,2)), 'k','Linewidth',[2])
plot(t,abs(nl2_Psi(:,3)), 'b','Linewidth',[2])
plot(t,abs(nl2_Psi(:,23)), 'r','Linewidth',[2])
plot(t,abs(nl2_Psi(:,24)), 'g','Linewidth',[2])
xlim([0 5])
ylim([0 0.5])
title('nonlinear cubic modes')


subplot(3,1,3)
plot(t,abs(nl4_Psi(:,1)), 'm','Linewidth',[2])
hold on
plot(t,abs(nl4_Psi(:,2)), 'k','Linewidth',[2])
plot(t,abs(nl4_Psi(:,3)), 'b','Linewidth',[2])
plot(t,abs(nl4_Psi(:,23)), 'r','Linewidth',[2])
plot(t,abs(nl4_Psi(:,24)), 'g','Linewidth',[2])
xlim([0 5])
ylim([0 0.5])
title('nonlinear quintic modes')

set(gcf,'Position',[100 100 300 300]); 





















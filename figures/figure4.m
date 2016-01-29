clear all; close all; clc;

load nl2_Psi.mat
load nl4_Psi.mat

L=80; n=1024;
t2=linspace(-L/2,L/2,n+1); t=t2(1:n);
sli=ones(1,1024);

p=513:512+50;
tp=t(p);
slip=sli(p);

%%
figure(1)
 subplot(2,1,1)
hold on
 plot3(tp,1*slip,abs(nl2_Psi(p,1)),'m','Linewidth',[1.5]), hold on
% plot3(tp,4*slip,abs(nl2_Psi(p,2)),'c','Linewidth',[1])

plot3(tp,5*slip,abs(nl2_Psi(p,3)),'b','Linewidth',[1.5])
plot3(tp,9*slip,abs(nl2_Psi(p,4)),'b','Linewidth',[1.5])
plot3(tp,13*slip,abs(nl2_Psi(p,5)),'b','Linewidth',[1.5])
plot3(tp,17*slip,abs(nl2_Psi(p,6)),'b','Linewidth',[1.5])
plot3(tp,21*slip,abs(nl2_Psi(p,7)),'b','Linewidth',[1.5])
plot3(tp,25*slip,abs(nl2_Psi(p,8)),'b','Linewidth',[1.5])

plot3(tp,29*slip,abs(nl2_Psi(p,23)),'r','Linewidth',[1.5])
% plot3(tp,28*slip,abs(nl2_Psi(p,24)),'g','Linewidth',[1])
 grid on


%%
n1=513; n2=519; n3=526;

hold on
a=0:0.1:30;
for j=[n1, n2, n3]
    
    plot3(t(j)+0*a, a, 0*a, 'k', 'Linewidth',[1.5])
    
    b1=[0 abs(nl2_Psi(j,1))];
    b3=[0 abs(nl2_Psi(j,3))];
    b4=[0 abs(nl2_Psi(j,4))];
    b5=[0 abs(nl2_Psi(j,5))];
    b6=[0 abs(nl2_Psi(j,6))];
    b7=[0 abs(nl2_Psi(j,7))];
    b8=[0 abs(nl2_Psi(j,8))];
    b23=[0 abs(nl2_Psi(j,23))];

    plot3(t(j)+0*b1, 1+0*b1, b1, 'k', 'Linewidth',[1.5])
    plot3(t(j)+0*b3, 5+0*b3, b3, 'k', 'Linewidth',[1.5])
    plot3(t(j)+0*b4, 9+0*b4, b4, 'k', 'Linewidth',[1.5])
    plot3(t(j)+0*b5, 13+0*b5, b5, 'k', 'Linewidth',[1.5])
    plot3(t(j)+0*b6, 17+0*b6, b6, 'k', 'Linewidth',[1.5])
    plot3(t(j)+0*b7, 21+0*b7, b7, 'k', 'Linewidth',[1.5])
    plot3(t(j)+0*b8, 25+0*b8, b8, 'k', 'Linewidth',[1.5])
    plot3(t(j)+0*b23, 29+0*b23, b23, 'k', 'Linewidth',[1.5])

    view(54, 62)
    % view([180 72])
    % xlim([0 4])
     set(gcf,'Position',[100 100 300 300]); 
    %  hold on
    % plot3(tp,1*slip,abs(nl2_Psi(p,1)),'k','Linewidth',[1]) 

end
title('nonlinear cubic modes')



subplot(2,1,2)
%figure(2)
hold on
 plot3(tp,1*slip,abs(nl4_Psi(p,1)),'m','Linewidth',[1.5]), hold on
% plot3(tp,4*slip,abs(nl2_Psi(p,2)),'c','Linewidth',[1])

plot3(tp,5*slip,abs(nl4_Psi(p,3)),'b','Linewidth',[1.5])
plot3(tp,9*slip,abs(nl4_Psi(p,4)),'b','Linewidth',[1.5])
plot3(tp,13*slip,abs(nl4_Psi(p,5)),'b','Linewidth',[1.5])
plot3(tp,17*slip,abs(nl4_Psi(p,6)),'b','Linewidth',[1.5])
plot3(tp,21*slip,abs(nl4_Psi(p,7)),'b','Linewidth',[1.5])
plot3(tp,25*slip,abs(nl4_Psi(p,8)),'b','Linewidth',[1.5])

plot3(tp,29*slip,abs(nl4_Psi(p,23)),'r','Linewidth',[1.5])
% plot3(tp,28*slip,abs(nl2_Psi(p,24)),'g','Linewidth',[1])
 grid on


%%
n1=513; n2=519; n3=526;

hold on
a=0:0.1:30;
for j=[n1, n2, n3]
    
    plot3(t(j)+0*a, a, 0*a, 'k', 'Linewidth',[1.5])
    
    b1=[0 abs(nl4_Psi(j,1))];
    b3=[0 abs(nl4_Psi(j,3))];
    b4=[0 abs(nl4_Psi(j,4))];
    b5=[0 abs(nl4_Psi(j,5))];
    b6=[0 abs(nl4_Psi(j,6))];
    b7=[0 abs(nl4_Psi(j,7))];
    b8=[0 abs(nl4_Psi(j,8))];
    b23=[0 abs(nl4_Psi(j,23))];

    plot3(t(j)+0*b1, 1+0*b1, b1, 'k', 'Linewidth',[1.5])
    plot3(t(j)+0*b3, 5+0*b3, b3, 'k', 'Linewidth',[1.5])
    plot3(t(j)+0*b4, 9+0*b4, b4, 'k', 'Linewidth',[1.5])
    plot3(t(j)+0*b5, 13+0*b5, b5, 'k', 'Linewidth',[1.5])
    plot3(t(j)+0*b6, 17+0*b6, b6, 'k', 'Linewidth',[1.5])
    plot3(t(j)+0*b7, 21+0*b7, b7, 'k', 'Linewidth',[1.5])
    plot3(t(j)+0*b8, 25+0*b8, b8, 'k', 'Linewidth',[1.5])
    plot3(t(j)+0*b23, 29+0*b23, b23, 'k', 'Linewidth',[1.5])

    view(54, 62)
    % view([180 72])
    % xlim([0 4])
     set(gcf,'Position',[100 100 300 300]); 
    %  hold on
    % plot3(tp,1*slip,abs(nl2_Psi(p,1)),'k','Linewidth',[1]) 

end
title('nonlinear quintic modes')





% 
% 513   530   507 
% 
% 
% 522   480   491   
% 
% 
% 506   531   536
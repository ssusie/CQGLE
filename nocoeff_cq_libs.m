clear all; close all; clc;

load A1.mat
A1=abs(A1.');

% waterfall(abs(A1(300:600,100:161)'))
% break
% 

st=100;
fin=161;
X=A1(:,st:end);

%compute the nonlinear part
n=1024;
NL=zeros(n,fin-st+1);

% tau=-0.03;
% kappa=-0.05;
% beta=1.45; %mu
% nu=0;
% sigma=-0.1; %eps
% gamma=-0.5;

for j=1:fin-st+1
    
    NL2(:,j)=(abs(X(:,j))).^2.*X(:,j);
    NL4(:,j)=(abs(X(:,j))).^4.*X(:,j);
    
end

[u2,s2,v2]=svd(abs(NL2),0);
[u4,s4,v4]=svd(abs(NL4),0);
nl2_Psi=[u2(:,1)];
nl4_Psi=[u4(:,1)];

figure(11)
plot(diag(s2)/sum(diag(s2)), 'ko')
xlim([0 10])
hold on
plot(diag(s4)/sum(diag(s4)), 'ro')
xlim([0 10])

title('svals nl2 nl4 beta1')
legend('nl2','nl4')

figure(12)

plot(cumsum(diag(s2)/sum(diag(s2))), 'ko')
xlim([0 10])
hold on
plot(cumsum(diag(s4)/sum(diag(s4))), 'ro')
xlim([0 10])

title('cumsum svals beta 2')
legend('nl2','nl4')
% break
clear NL2 NL4 X u2 s2 v2 u4 s4 v4
 
%%

load A2.mat
A2=abs(A2.');
st=100;
fin=161;
X=A2(:,st:end);
%compute the nonlinear part
n=1024;
NL=zeros(n,fin-st+1);

% tau=-0.03;
% kappa=-0.05;
% beta=1.4; %mu
% nu=0;
% sigma=-0.1; %eps
% gamma=-0.5;

for j=1:fin-st+1
    
    NL2(:,j)=(abs(X(:,j))).^2.*X(:,j);
    NL4(:,j)=(abs(X(:,j))).^4.*X(:,j);
    
end

[u2,s2,v2]=svd(abs(NL2),0);
[u4,s4,v4]=svd(abs(NL4),0);
nl2_Psi=[nl2_Psi u2(:,1)];
nl4_Psi=[nl4_Psi u4(:,1)];

figure(21)
plot(diag(s2)/sum(diag(s2)), 'ko')
xlim([0 10])
hold on
plot(diag(s4)/sum(diag(s4)), 'ro')
xlim([0 10])

title('svals nl2 nl4 beta2 ')
legend('nl2','nl4')

figure(22)

plot(cumsum(diag(s2)/sum(diag(s2))), 'ko')
xlim([0 10])
hold on
plot(cumsum(diag(s4)/sum(diag(s4))), 'ro')
xlim([0 10])

title('cumsum svals beta2')
legend('nl2','nl4')
% break
clear NL2 NL4 X u2 s2 v2 u4 s4 v4


%%
load A3.mat
A3=abs(A3.');
st=100;
fin=161;
X=A3(:,st:end);

%compute the nonlinear part
n=1024;
NL=zeros(n,fin-st+1);

% tau=0.08;
% kappa=0;
% beta=0.66; %mu
% nu=-0.1;
% sigma=-0.1; %eps
% gamma=-0.1;

for j=1:fin-st+1
    
    NL2(:,j)=(abs(X(:,j))).^2.*X(:,j);
    NL4(:,j)=(abs(X(:,j))).^4.*X(:,j);
    
end

[u2,s2,v2]=svd(abs(NL2),0);
[u4,s4,v4]=svd(abs(NL4),0);
nl2_Psi=[nl2_Psi u2(:,1:6)];
nl4_Psi=[nl4_Psi u4(:,1:6)];

figure(31)
plot(diag(s2)/sum(diag(s2)), 'ko')
xlim([0 10])
hold on
plot(diag(s4)/sum(diag(s4)), 'ro')
xlim([0 10])

title('svals nl2 nl4 beta3 ')
legend('nl2','nl4')

figure(32)

plot(cumsum(diag(s2)/sum(diag(s2))), 'ko')
xlim([0 10])
hold on
plot(cumsum(diag(s4)/sum(diag(s4))), 'ro')
xlim([0 10])

title('cumsum svals beta3')
legend('nl2','nl4')
% break
clear NL2 NL4 X u2 s2 v2 u4 s4 v4


%%
load A4.mat
A4=abs(A4.');
st=100;
fin=161;
X=A4(:,st:end);

%compute the nonlinear part
n=1024;
NL=zeros(n,fin-st+1);

% tau=0.125;
% kappa=0;
% beta=1; %mu
% nu=-0.6;
% sigma=-0.1; %eps
% gamma=-0.1;

for j=1:fin-st+1
    
    NL2(:,j)=(abs(X(:,j))).^2.*X(:,j);
    NL4(:,j)=(abs(X(:,j))).^4.*X(:,j);
    
end

[u2,s2,v2]=svd(abs(NL2),0);
[u4,s4,v4]=svd(abs(NL4),0);
%%%%%%11 modes are neough though
nl2_Psi=[nl2_Psi u2(:,1:14)];
nl4_Psi=[nl4_Psi u4(:,1:14)];

figure(41)
plot(diag(s2)/sum(diag(s2)), 'ko')
xlim([0 20])
hold on
plot(diag(s4)/sum(diag(s4)), 'ro')
xlim([0 20])

title('svals nl2 nl4 beta4 ')
legend('nl2','nl4')

figure(42)

plot(cumsum(diag(s2)/sum(diag(s2))), 'ko')
xlim([0 20])
hold on
plot(cumsum(diag(s4)/sum(diag(s4))), 'ro')
xlim([0 20])
ylim([0.98 1])
title('cumsum svals beta4')
legend('nl2','nl4')

clear NL2 NL4 X u2 s2 v2 u4 s4 v4



%%

load A5.mat
A5=abs(A5.');
st=100;
fin=161;
X=A5(:,st:end);


%compute the nonlinear part
n=1024;
NL=zeros(n,fin-st+1);

% tau=0.08;
% kappa=-0.05;
% beta=0.6; %mu
% nu=-0.1;
% sigma=-0.1; %eps
% gamma=-0.1;

for j=1:fin-st+1
    
    NL2(:,j)=(abs(X(:,j))).^2.*X(:,j);
    NL4(:,j)=(abs(X(:,j))).^4.*X(:,j);
    
end

[u2,s2,v2]=svd(abs(NL2),0);
[u4,s4,v4]=svd(abs(NL4),0);

nl2_Psi=[nl2_Psi u2(:,1)];
nl4_Psi=[nl4_Psi u4(:,1)];

figure(51)
plot(diag(s2)/sum(diag(s2)), 'ko')
xlim([0 10])
hold on
plot(diag(s4)/sum(diag(s4)), 'ro')
xlim([0 10])

title('svals nl2 nl4 beta5 ')
legend('nl2','nl4')

figure(52)

plot(cumsum(diag(s2)/sum(diag(s2))), 'ko')
xlim([0 10])
hold on
plot(cumsum(diag(s4)/sum(diag(s4))), 'ro')
xlim([0 10])

title('cumsum svals beta5')
legend('nl2','nl4')

clear NL2 NL4 X u2 s2 v2 u4 s4 v4




%%
load A6.mat
A6=abs(A6.');
st=100;
fin=161;
X=A6(:,st:end);

[VX,SX,WX]=svd(X,0);
figure(61)
plot(diag(SX)/sum(diag(SX)), 'ko')
figure(62)
plot(cumsum(diag(SX)/sum(diag(SX))), 'ko')


%compute the nonlinear part
n=1024;
NL=zeros(n,fin-st+1);

% tau=0.08;
% kappa=-0.05;
% beta=0.5; %mu
% nu=-0.1;
% sigma=-0.1; %eps
% gamma=-0.1;

for j=1:fin-st+1
    
    NL2(:,j)=(abs(X(:,j))).^2.*X(:,j);
    NL4(:,j)=(abs(X(:,j))).^4.*X(:,j);
    
end

[u2,s2,v2]=svd(abs(NL2),0);
[u4,s4,v4]=svd(abs(NL4),0);

nl2_Psi=[nl2_Psi u2(:,1)];
nl4_Psi=[nl4_Psi u4(:,1)];

figure(61)
plot(diag(s2)/sum(diag(s2)), 'ko')
xlim([0 10])
hold on
plot(diag(s4)/sum(diag(s4)), 'ro')
xlim([0 10])

title('svals nl2 nl4 beta6 ')
legend('nl2','nl4')

figure(62)

plot(cumsum(diag(s2)/sum(diag(s2))), 'ko')
xlim([0 10])
hold on
plot(cumsum(diag(s4)/sum(diag(s4))), 'ro')
xlim([0 10])

title('cumsum svals beta6')
legend('nl2','nl4')

clear NL2 NL4 X u2 s2 v2 u4 s4 v4


save nl2_Psi nl2_Psi
save nl4_Psi nl4_Psi


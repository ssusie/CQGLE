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

beta=1.45; %mu
nu=0;
sigma=-0.1; %eps

for j=1:fin-st+1
    
    NL(:,j)=(i+beta)*(abs(X(:,j))).^2.*X(:,j)+...
            (i*nu+sigma)*(abs(X(:,j))).^4.*X(:,j);
    
end

[u,s,v]=svd(abs(NL),0);

nl_Psi=[u(:,1)];

figure(11)
plot(diag(s)/sum(diag(s)), 'ko')
xlim([0 10])


title('svals nl beta1')


figure(12)

plot(cumsum(diag(s)/sum(diag(s))), 'ko')
xlim([0 10])


title('cumsum svals beta1')
clear NL X u s v 
 
%%

load A2.mat
A2=abs(A2.');
st=100;
fin=161;
X=A2(:,st:end);
%compute the nonlinear part
n=1024;
NL=zeros(n,fin-st+1);

beta=1.4; %mu
nu=0;
sigma=-0.1; %eps

for j=1:fin-st+1
    
    NL(:,j)=(i+beta)*(abs(X(:,j))).^2.*X(:,j)+...
            (i*nu+sigma)*(abs(X(:,j))).^4.*X(:,j);
    
end

[u,s,v]=svd(abs(NL),0);

nl_Psi=[nl_Psi u(:,1)];


figure(21)
plot(diag(s)/sum(diag(s)), 'ko')
xlim([0 10])

title('svals nl beta2 ')

figure(22)

plot(cumsum(diag(s)/sum(diag(s))), 'ko')
xlim([0 10])


title('cumsum svals beta2')

clear NL X u s v 

%%
load A3.mat
A3=abs(A3.');
st=100;
fin=161;
X=A3(:,st:end);

%compute the nonlinear part
n=1024;
NL=zeros(n,fin-st+1);

beta=0.66; %mu
nu=-0.1;
sigma=-0.1; %eps

for j=1:fin-st+1
    
    NL(:,j)=(i+beta)*(abs(X(:,j))).^2.*X(:,j)+...
            (i*nu+sigma)*(abs(X(:,j))).^4.*X(:,j);
    
end

[u,s,v]=svd(abs(NL),0);
nl_Psi=[nl_Psi u(:,1:6)];

figure(31)
plot(diag(s)/sum(diag(s)), 'ko')
xlim([0 10])

title('svals nl beta3 ')


figure(32)

plot(cumsum(diag(s)/sum(diag(s))), 'ko')
xlim([0 10])

title('cumsum svals beta3')

% break
clear NL  X u s v 


%%
load A4.mat
A4=abs(A4.');
st=100;
fin=161;
X=A4(:,st:end);

%compute the nonlinear part
n=1024;
NL=zeros(n,fin-st+1);

beta=1; %mu
nu=-0.6;
sigma=-0.1; %eps

for j=1:fin-st+1
    
    NL(:,j)=(i+beta)*(abs(X(:,j))).^2.*X(:,j)+...
            (i*nu+sigma)*(abs(X(:,j))).^4.*X(:,j);
    
end

[u,s,v]=svd(abs(NL),0);

%%%%%% 11 modes are neough though
nl_Psi=[nl_Psi u(:,1:14)];

figure(41)
plot(diag(s)/sum(diag(s)), 'ko')
xlim([0 20])

title('svals nl beta4 ')


figure(42)

plot(cumsum(diag(s)/sum(diag(s))), 'ko')
xlim([0 20])
ylim([0.98 1])
title('cumsum svals beta4')

clear NL X u s v 


%%

load A5.mat
A5=abs(A5.');
st=100;
fin=161;
X=A5(:,st:end);


%compute the nonlinear part
n=1024;
NL=zeros(n,fin-st+1);

beta=0.6; %mu
nu=-0.1;
sigma=-0.1; %eps

for j=1:fin-st+1
    
    NL(:,j)=(i+beta)*(abs(X(:,j))).^2.*X(:,j)+...
             (i*nu+sigma)*(abs(X(:,j))).^4.*X(:,j);
    
end

[u,s,v]=svd(abs(NL),0);

nl_Psi=[nl_Psi u(:,1)];

figure(51)
plot(diag(s)/sum(diag(s)), 'ko')
xlim([0 10])

title('svals nl beta5 ')


figure(52)

plot(cumsum(diag(s)/sum(diag(s))), 'ko')
xlim([0 10])
title('cumsum svals beta5')

clear NL X u s v 



%%
load A6.mat
A6=abs(A6.');
st=100;
fin=161;
X=A6(:,st:end);

% [VX,SX,WX]=svd(X,0);
% figure(61)
% plot(diag(SX)/sum(diag(SX)), 'ko')
% figure(62)
% plot(cumsum(diag(SX)/sum(diag(SX))), 'ko')


%compute the nonlinear part
n=1024;
NL=zeros(n,fin-st+1);

beta=0.5; %mu
nu=-0.1;
sigma=-0.1; %eps

for j=1:fin-st+1
    
    NL(:,j)=(i+beta)*(abs(X(:,j))).^2.*X(:,j)+...
            (i*nu+sigma)*(abs(X(:,j))).^4.*X(:,j);
    
end

[u,s,v]=svd(abs(NL),0);

nl_Psi=[nl_Psi u(:,1)];

figure(61)
plot(diag(s)/sum(diag(s)), 'ko')
xlim([0 10])
title('svals nl beta6 ')


figure(62)

plot(cumsum(diag(s)/sum(diag(s))), 'ko')
xlim([0 10])

title('cumsum svals beta6')
clear NL X u s v 

save NL nl_Psi



clear all; close all; clc;

load master
%  break
%%%%%%%%% betta_1 regime %%%%%%%%%

st=27;
fin=41;
% sample=[st:fin st+41:fin+41 st+82:fin+82];
% break

%%
%%%%%% b_1
X1=umaster(st:fin,:)';

%compute the nonlinear part
NL1=zeros(n,fin-st+1);
beta=1.45; %mu
nu=0;
sigma=-0.1; %eps

for j=1:fin-st+1
    
   NL1(:,j)=(i+beta)*(abs(X1(:,j))).^2.*X1(:,j)+...
        (i*nu+sigma)*(abs(X1(:,j))).^4.*X1(:,j);
    
end


%%
%%%%%% b_3
X3=umaster(st+41:fin+41,:)';
beta=0.66; %mu
nu=-0.1;
sigma=-0.1; %eps

for j=1:fin-st+1
    
   NL2(:,j)=(i+beta)*(abs(X3(:,j))).^2.*X3(:,j)+...
        (i*nu+sigma)*(abs(X3(:,j))).^4.*X3(:,j);
    
end



%%
%%%%%%% b_5
X5=umaster(st+82:fin+82,:)';
beta=0.6; %mu
nu=-0.1;
sigma=-0.1; %eps

for j=1:fin-st+1
    
   NL3(:,j)=(i+beta)*(abs(X5(:,j))).^2.*X5(:,j)+...
        (i*nu+sigma)*(abs(X5(:,j))).^4.*X5(:,j);
    
end

%%

NL=[NL1 NL2 NL3];
[u,s,v]=svd(abs(NL),0);

% figure(1)
% plot(cumsum(diag(s2)/sum(diag(s2))), 'ko')
% hold on
% plot(cumsum(diag(s4)/sum(diag(s4))), 'rv')
% break

m=3; %number of sensors
[ro,g(1)]=max(abs(u(:,1)));
U=[u(:,1)]; 
z=zeros(n,1);
P=z; P(g(1))=1;

for l=2:m
    c=(P'*U)\(P'*u(:,l));
    r=u(:,l)-U*c;
    [ro,g(l)]=max(abs(r));
    U=[U,u(:,l)]; 
    P=[P,z]; P(g(l),l)=1;
   
end


abs(g-513)



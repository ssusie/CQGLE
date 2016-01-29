clear all; close all;clc;

load A1
% break
%%%%%%%%% betta_1 regime %%%%%%%%%

st=1;
fin=161;
n=1024;
sample=[st:fin];
% break

% X=abs(A1');
X=A1.';
m=3;

%compute the nonlinear part

for j=1:fin-st+1
    
    NL2(:,j)= (abs(X(:,j))).^2.*X(:,j);
    NL4(:,j)= (abs(X(:,j))).^4.*X(:,j);
    
end

[u2,s2,v2]=svd(abs(NL2),0);
[u4,s4,v4]=svd(abs(NL4),0);
figure(1)
plot(cumsum(diag(s2)/sum(diag(s2))), 'ko')
hold on
plot(cumsum(diag(s4)/sum(diag(s4))), 'rv')
% break

[ro2,g2(1)]=max(abs(u2(:,1)));
U2=[u2(:,1)]; 
z=zeros(n,1);
P2=z; P2(g2(1))=1;

for l=2:m
    c=(P2'*U2)\(P2'*u2(:,l));
    r2=u2(:,l)-U2*c;
    [ro2,g2(l)]=max(abs(r2));
    U2=[U2,u2(:,l)]; 
    P2=[P2,z]; P2(g2(l),l)=1;
   
end

[ro4,g4(1)]=max(abs(u4(:,1)));
U4=[u4(:,1)]; 
z=zeros(n,1);
P4=z; P4(g4(1))=1;

for l=2:m
    c=(P4'*U4)\(P4'*u4(:,l));
    r=u4(:,l)-U4*c;
    [ro4,g4(l)]=max(abs(r));
    U4=[U4,u4(:,l)]; 
    P4=[P4,z]; P4(g4(l),l)=1;
   
end

abs(513-g2)
abs(513-g4)









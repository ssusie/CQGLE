clear all; close all; clc;

x=0:12; 
y=0:34;
z=ones(1, length(x));
w=ones(1, length(y));
 
P=[0  0  0  0  0  0  0  0  0  9  0  0;
   5  15 12 6  4  13 10 6  6  21 6 6;
   13 26 17 22 13 23 15 20 13 32 15 13];
P=flipud(P);
 B=zeros(33, 12);
%%
% % Create a blue rectangle
figure(1)
subplot(1,3,1)
h1=axes;
for i=1:size(P,1)
    for j=1:size(P,2)
        left = j-1;
        right = left + 1 ;%j;
        bottom = P(i,j);
        top = bottom + 1;
        xx = [left left right right];
        yy = [bottom top top bottom];
         if (j~=12)%if (j~=4) &(j~=8)&(j~=12)
            fill(xx, yy, 'k');hold on;
        else 
            fill(xx, yy, 'r');hold on;
        end
        

         B(P(i,j)+1, j)=1;
    end
end
 

% break

for j=1:length(y)

    plot(x, y(j)*z,'k')
    hold on
    
end


for j=1:length(x)
    
    plot( x(j)*w,y, 'k')
    hold on
end

xlim([0 12])
ylim([0 34])
% set(gcf,'Position',[100 100 700 150])

set(h1, 'Ydir', 'reverse')
set(h1, 'YAxisLocation', 'Right') 

%%
[u,s,v]=svd(B,0);

% figure(2)
% subplot(2,1,1)
% plot(abs(diag(s))/sum(abs(diag(s))), 'ko')
% title('svd')
% subplot(2,1,2)
% plot(cumsum(diag(s))/sum(diag(s)), 'b*')
% title('cumsum')

% figure(3)
% for j=1:5%11
%     plot(abs(u(:,j)))
%     hold on
% end
% 
%%

load nl2_Psi
M=abs(nl2_Psi(513:546,1));
B1=[4,5, 9, 10, 12,15,17,20, 21,22, 23, 26, 32];
B2=[0, 6, 13];
% figure(5)
% plot(abs(nl2_Psi(513:512+33,1)), 'm')
figure(4)
A=subplot(1,3,1)
plot(0:33,M, 'm')
xlim([-1 33])
e=[1,1];
 % hold on
%plot(B(3)*e, [0, M(B(3)+1)] , '--')
% h1=axes;
% set(h1, 'xdir', 'reverse')
% set(h1, 'XAxisLocation', 'Right') 
% hold on 
for k=1: length(B1)
    hold on
    plot(B1(k)*e, [0, M(B1(k)+1)] , 'k--*')
end

for k=1: length(B2)
    hold on
    plot(B2(k)*e, [0, M(B2(k)+1)] , 'r--*')
end

 view(90, 90)
% set(gcf,'Position',[100 100 300 150])
set(A,'YTick',[]) 
set(A,'XTick',[]) 

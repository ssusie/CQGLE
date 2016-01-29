clear all; close all; clc
load sens5_400

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


close all
clear all
clc

[xx, yy] = textread('hpts.in','%f %f','headerlines',1); 

tf = length(xx);
tt = tf/2;
ss1 = zeros(1:2);

plot(xx(1:tt),yy(1:tt),'b*','markersize',4)
hold on
plot(xx(tt+1:tf),yy(tt+1:tf),'ro','markersize',4)
plot(ss1,[-1 1],'g')
plot([-1 1],ss1,'g')
daspect([1 1 1])
%% Homework 1 - Francisco Castillo
clear all; close all; clc
format long
param.a = 2.12;
param.b = -0.3;
param.maxiter = 100;
param.L = 0.01;

xmin=-3;
xmax=3;
ymin=-3;
ymax=3;
n=1000
x=linspace(xmin,xmax,n);
y=linspace(ymin,ymax,n);
basin = basin_henon(n,x,y,param);
% x = linspace(xmin,xmax,n)';
% y = linspace(ymax,ymin,n)';
% [X,Y]  = meshgrid(x,y)
% v(:,1) = reshape(X,[n^2,1])
% v(:,2) = reshape(Y,[n^2,1])
% x
% y
% basin = reshape((v(:,1)>0 | v(:,2)>0),[n,n])
h=surf(x,y,basin)
view(0,90)
set(h,'edgecolor','none')
saveas(gcf,'henon','epsc')
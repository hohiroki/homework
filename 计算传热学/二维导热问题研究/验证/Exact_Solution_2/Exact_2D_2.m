clc;
clear all;
close all;

x=0:0.01:pi;
y=0:0.01:pi;
[xx,yy]=meshgrid(x,y);

t=0.5;

T=600.*sin(xx).*sin(yy).*exp(-2*t);

surf(xx,yy,0*xx,T);view(2);axis  equal ;set(gcf,'renderer','zbuffer');
shading   interp,  
colorbar
figure
contour(xx,yy,T);
fid=fopen('E:\Exact_2D.dat','wb');
csvwrite('E:\Exact_2D.dat',[xx(:),yy(:),T(:)])
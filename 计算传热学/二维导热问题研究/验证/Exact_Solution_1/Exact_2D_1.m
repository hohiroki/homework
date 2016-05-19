clc;
clear all;
close all;

x=0:0.01:1;
y=0:0.01:1;
[xx,yy]=meshgrid(x,y);
n=0:1:200;
T=zeros(size(xx));
for i=1:length(n)
   T=T+4*(-1)^n(i)/(2*n(i)+1)/pi*cos((2*n(i)+1)*pi/2*xx).*sinh((2*n(i)+1)*pi/2*yy)/sinh((2*n(i)+1)*pi/2);
end
T=300*(T+1);
surf(xx,yy,0*xx,T);view(2);axis  equal ;set(gcf,'renderer','zbuffer');
shading   interp,  
colorbar
!figure
!contour(xx,yy,T,300:10:max(T(:)));
fid=fopen('Exact_2D.dat','wb')
csvwrite('Exact_2D.dat',[xx(:),yy(:),T(:)])
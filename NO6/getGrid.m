function [X,Y,xnum,ynum] = getGrid()
%描述喷管形状
Rt = 5*1e-3;
Re = Rt*2^0.5;
Rb = 1.5*Rt;
theta = 15*pi/180;
%自己定义一个Ld
Ld = 2*1e-3;
%afa角度与fai角度
afa = 20*pi/180;
fai = 60*pi/180;
Le = (Re-Rt)/tan(theta);
L = (2*Rt+Le*tan(fai))/(tan(fai)-2*tan(theta));
Rout = L*tan(theta)+Rt;
Lx = Ld+(Rb-Rt)/tan(afa)+L;
Ly = 2*Rout;
%画网格
%x方向分为3部分
xnum1 = 4; %入口段
xnum2 = 7; %收缩段
xnum3 = 15; %扩张段
xnum = xnum1+xnum2+xnum3;
ynum = 8;
X = ones(xnum,ynum);
Y = ones(xnum,ynum);
for i = 1:1:xnum1
    for j = 1:1:ynum
        X(i,j) = (i-1)*Ld/(xnum1-1);
        Y(i,j) = (j-1)*2*Rb/(ynum-1)-Rb;
    end
end
for i = xnum1+1:1:xnum1+xnum2
    for j = 1:1:ynum
        x=Ld+(i-xnum1)*(Rb-Rt)/tan(afa)/(xnum2);
        x1 = Ld;
        y1 = (j-1)*2*Rb/(ynum-1)-Rb;
        x2 = Ld+(Rb-Rt)/tan(afa);
        y2 = (j-1)*2*Rt/(ynum-1)-Rt;
        y = y1+(y2-y1)/(x2-x1)*(x-x1);
        X(i,j) = x;
        Y(i,j) = y;
    end
end
%控制斜率
for i = xnum1+xnum2+1:1:xnum
    for j = 1:1:ynum
        x1 = Ld+(Rb-Rt)/tan(afa);
        y1 = (j-1)*2*Rt/(ynum-1)-Rt;
        x2 = Ld+(Rb-Rt)/tan(afa)+L-cos(fai)*(j-1)*2*Rout/sin(fai)/(ynum-1);
        y2 = (j-1)*2*Rout/(ynum-1)-Rout;
        x3 = x1+(i-xnum1-xnum2)*Le/xnum3;
        y3 = Rt+(x3-x1)*tan(theta);
        x4 = x1+(i-xnum1-xnum2)*L/xnum3;
        y4 = -Rt-(x4-x1)*tan(theta);
        k = (y2-y1)/(x2-x1)-(y4-y3)/(x4-x3);
        x = (y3-y1+x1*(y2-y1)/(x2-x1)-x3*(y4-y3)/(x4-x3))/k;
        y = (y2-y1)/(x2-x1)*(x-x1)+y1;
        X(i,j) = x;
        Y(i,j) = y;
    end
end
end
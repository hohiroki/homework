clear all;
clc;
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
xnum1 = 10; %入口段
xnum2 = 20; %收缩段
xnum3 = 60; %扩张段
xnum = xnum1+xnum2+xnum3;
ynum = 30;
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
%画网格图
% for j=1:1:ynum    
%    plot(X(:,j),Y(:,j));hold on; 
% end
% for i = 1:1:xnum
%     plot(X(i,:),Y(i,:));hold on;
% end
% box off;
%计算部分
dt = 1/(xnum*ynum*10);
gama = 1.26;
T0 = 3200;
P0 = 6*1e6;
R = 287.314;
eps = 1;
epsmax = 1e-6;

%入口边界条件
Main = 0.4357;
Tin = T0/(1+(gama-1)/2*Main^2); %静温
Pin = P0/(1+(gama-1)/2*Main^2)^(gama/(gama-1));%静压
ain = (gama*R*T0)^0.5;%声速
uin = Main*ain; %速度
vin = 0;
rouin = (P0-Pin)*gama*Main^2/uin^2/((1+(gama-1)/2*Main^2)^(gama/(gama-1))-1);  %密度
Ein = Pin/(gama-1)/rouin+1/2*uin^2;  %总能
Hin = Ein+Pin/rouin;
%初始化
rou = ones(xnum,ynum)*rouin;
p = ones(xnum,ynum)*Pin;
u = ones(xnum,ynum)*uin;
v = ones(xnum,ynum)*vin;
E = ones(xnum,ynum)*Ein;
H = ones(xnum,ynum)*Hin;
rou_n = rou;
p_n = p;
u_n = u;
v_n = v;
E_n = E;
for step = 0:dt:1
    if eps<epsmax
        break
    end
    rou_n = rou;
    p_n = p;
    u_n = u;
    v_n = v;
    E_n = E;
    H_n = H;
    %内部节点
    for i = 2:1:xnum-1
        for j = 2:1:ynum-1
            dx = (X(i+1,j)-X(i-1,j))/2;
            dy = (Y(i,j+1)-Y(i,j-1))/2;
            %界面通量计算
            rouL = sqrt(rou_n(i-1,j));
            rouR = sqrt(rou_n(i,j)); 
            rouW=rouL*rouR;
            uW=(u_n(i-1,j)*rouL+u_n(i,j)*rouR)/(rouL+rouR);
            vW=(v_n(i-1,j)*rouL+v_n(i,j)*rouR)/(rouL+rouR);
            HW=(H_n(i-1,j)*rouL+H_n(i,j)*rouR)/(rouL+rouR);
            pW=(HW-0.5*(uW^2+vW^2))*rouW*(gama-1)/gama;
            
            rouL = sqrt(rou_n(i,j));
            rouR = sqrt(rou_n(i+1,j)); 
            rouE=rouL*rouR;
            uE=(u_n(i,j)*rouL+u_n(i+1,j)*rouR)/(rouL+rouR);
            vE=(v_n(i,j)*rouL+v_n(i+1,j)*rouR)/(rouL+rouR);
            HE=(H_n(i,j)*rouL+H_n(i+1,j)*rouR)/(rouL+rouR);
            pE=(HE-0.5*(uE^2+vE^2))*rouE*(gama-1)/gama;
                       
            rouL = sqrt(rou_n(i,j-1));
            rouR = sqrt(rou_n(i,j)); 
            rouS=rouL*rouR;
            uS=(u_n(i,j-1)*rouL+u_n(i,j)*rouR)/(rouL+rouR);
            vS=(v_n(i,j-1)*rouL+v_n(i,j)*rouR)/(rouL+rouR);
            HS=(H_n(i,j-1)*rouL+H_n(i,j)*rouR)/(rouL+rouR);
            pS=(HS-0.5*(uS^2+vS^2))*rouS*(gama-1)/gama;
            
            rouL = sqrt(rou_n(i,j));
            rouR = sqrt(rou_n(i,j+1)); 
            rouN=rouL*rouR;
            uN=(u_n(i,j)*rouL+u_n(i,j+1)*rouR)/(rouL+rouR);
            vN=(v_n(i,j)*rouL+v_n(i,j+1)*rouR)/(rouL+rouR);
            HN=(H_n(i,j)*rouL+H_n(i,j+1)*rouR)/(rouL+rouR);
            pN=(HN-0.5*(uN^2+vN^2))*rouN*(gama-1)/gama;
            
            %一阶时间离散
            rou(i,j) = rou_n(i,j)-(rouE*uE-rouW*uW)/dx*dt+(rouN*vN-rouS*vS)/dy*dt;
            u(i,j) = (rou_n(i,j)*u_n(i,j)-(rouE*uE^2+pE-rouW*uW^2-pW)/dx*dt+(rouN*uN*vN-rouS*uS*vS)/dy*dt)/rou(i,j);
            v(i,j) = (rou_n(i,j)*v_n(i,j)-(rouE*uE*vE-rouW*uW*vW)/dx*dt+(rouN*vN^2+pN-rouS*vS^2-pS)/dy*dt)/rou(i,j);
            E(i,j) = (rou_n(i,j)*E_n(i,j)-(rouE*uE*HE-rouW*uW*HW)/dx*dt+(rouN*vN*HN-rouS*vS*HS)/dy*dt)/rou(i,j);
            p(i,j) = (E(i,j)-0.5*(u(i,j)^2+v(i,j)^2))*rou(i,j)*(gama-1);
            H(i,j) = E(i,j)+p(i,j)/rou(i,j);            
            
        end
    end
    %壁面
    for i = 2:1:xnum-1
        dx = X(i+1,1)-X(i,1);
        dy = Y(i,2)-Y(i,1);
        l=(dx^2+dy^2)^0.5;
        nx=-dy/l;
        ny=dx/l;
        u(i,1) = u(i,2)*ny^2-v(i,2)*nx*ny;
        v(i,1) = -u(i,2)*nx*ny+v(i,2)*nx^2;
        a = (p(i,2)/rou(i,2))^0.5;
        rou(i,1) = ((gama-1)^2*(u(i,2)*nx+v(i,2)*ny-2*a/(gama-1))^2*rou(i,2)^gama/4/gama/p(i,2))^(1/(gama-1));
        p(i,1) = p(i,2)*rou(i,1)^gama/rou(i,2)^gama;
        E(i,1) = p(i,1)/(gama-1)/rou(i,1)+0.5*(u(i,1)^2+v(i,1)^2);
        H(i,1) = E(i,1)+p(i,1)/rou(i,1);
        
        dx = X(i+1,ynum)-X(i,ynum);
        dy = Y(i,ynum)-Y(i,ynum-1);
        l=(dx^2+dy^2)^0.5;
        nx=dy/l;
        ny=-dx/l;
        u(i,ynum) = u(i,ynum-1)*ny^2-v(i,ynum-1)*nx*ny;
        v(i,ynum) = -u(i,ynum-1)*nx*ny+v(i,ynum-1)*nx^2;
        a = (p(i,ynum-1)/rou(i,ynum-1))^0.5;
        rou(i,ynum) = ((gama-1)^2*(u(i,ynum-1)*nx+v(i,ynum-1)*ny-2*a/(gama-1))^2*rou(i,ynum-1)^gama/4/gama/p(i,ynum-1))^(1/(gama-1));
        p(i,ynum) = p(i,ynum-1)*rou(i,ynum)^gama/rou(i,ynum-1)^gama;
        E(i,ynum) = p(i,ynum)/(gama-1)/rou(i,ynum)+0.5*(u(i,ynum)^2+v(i,ynum)^2);
        H(i,ynum) = E(i,ynum)+p(i,ynum)/rou(i,ynum);
    end
    
    %出口程序中直接外推
    rou(xnum,:) = rou(xnum-1,:);
    u(xnum,:) = u(xnum-1,:);
    v(xnum,:) = v(xnum-1,:);
    p(xnum,:) = p(xnum-1,:);
    E(xnum,:) = E(xnum-1,:);
    H(xnum,:) = H(xnum-1,:);
  
    eps = max([norm(u-u_n,1)/norm(u,1),norm(v-v_n,1)/norm(v,1),norm(rou-rou_n,1)/norm(rou,1),norm(p-p_n,1)/norm(p,1)])
end
if eps<epsmax
    disp('done!');
    save 'data_no6.mat'
    pcolor(X,Y,u),colorbar;
end






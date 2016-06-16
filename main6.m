clear all;
clc;
%生成网格
[X,Y,xnum,ynum] = getGrid();
%画网格图
%plotGrid(X,Y)
%计算部分
%物性参数
gama = 1.26;
T0 = 3200;
P0 = 6*1e6;
R = 287.314;
Pref = 1e5;
Rouref = 1;
Uref = (Pref/Rouref);
Lref = 5e-3;
eps = 1;
epsmax = 1e-6;
dt = 1/(xnum*ynum*10);
%入口边界条件
Main = 0.4357;
Tin = T0/(1+(gama-1)/2*Main^2); %静温
pin = P0/(1+(gama-1)/2*Main^2)^(gama/(gama-1));%静压
ain = (gama*R*T0)^0.5;%声速
uin = Main*ain; %速度
vin = 0;
rouin = (P0-pin)*gama*Main^2/uin^2/((1+(gama-1)/2*Main^2)^(gama/(gama-1))-1);  %密度
Ein = pin/(gama-1)/rouin+1/2*uin^2;  %总能
Hin = Ein+pin/rouin;
%初始化
rou = ones(xnum,ynum)*rouin/Rouref;
p = ones(xnum,ynum)*pin/Pref;
u = ones(xnum,ynum)*uin/Uref;
v = ones(xnum,ynum)*vin/Uref;
E = ones(xnum,ynum)*Ein;
H = ones(xnum,ynum)*Hin;
rou(xnum1:xnum1+xnum2,:) = rouin*1.4/Rouref;
p(xnum1:xnum1+xnum2,:) = pin*0.7/Pref;
u(xnum1:xnum1+xnum2,:) = uin*1/Uref;
rou(xnum2:xnum,:) = rouin*0.8/Rouref;
p(xnum2:xnum,:) = pin*0.5/Pref;
u(xnum2:xnum,:) = uin*0.5/Uref;
p_n = p;
u_n = u;
v_n = v;
E_n = E;
step = 0;
for time = 0:dt:1
    if eps<epsmax
        break
    end
    step =step+1;
    if step ==4
        break
    end
    rou_n = rou+(ones(xnum,ynum)*1e-4);
    p_n = p+(ones(xnum,ynum)*1e-4);
    u_n = u;
    v_n = v;
    E_n = E;
    H_n = H;
    %壁面
    for i = 2:1:xnum-1
        dx = X(i+1,1)-X(i,1);
        dy = Y(i,2)-Y(i,1);
        l=(dx^2+dy^2)^0.5;
        nx=-dy/l;
        ny=dx/l;
        u(i,1) = u(i,2)*ny^2-v(i,2)*nx*ny;
        v(i,1) = -u(i,2)*nx*ny+v(i,2)*nx^2;
        rou(i,2) = max(rou(i,2),1e-4);
        a = (p(i,2)/rou(i,2))^0.5;
        rou(i,1) = ((gama-1)^2*(u(i,2)*nx+v(i,2)*ny-2*a/(gama-1))^2*rou(i,2)^gama/4/gama/p(i,2))^(1/(gama-1));
        rou(i,1) = max(rou(i,1),1e-4);
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
    
    %出口直接外推
    rou(xnum,:) = rou(xnum-1,:);
    u(xnum,:) = u(xnum-1,:);
    v(xnum,:) = v(xnum-1,:);
    p(xnum,:) = p(xnum-1,:);
    E(xnum,:) = E(xnum-1,:);
    H(xnum,:) = H(xnum-1,:);
    
    %入口给定
    rou(1,:) = rouin;
    u(1,:) = uin;
    v(1,:) = vin;
    p(1,:) = pin;
    E(1,:) = Ein;
    H(1,:) = Hin;
    
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
            rou(i,j) = rou_n(i,j)-(rouE*uE-rouW*uW)/dx*dt-(rouN*vN-rouS*vS)/dy*dt;
            u(i,j) = (rou_n(i,j)*u_n(i,j)-(rouE*uE^2+pE-rouW*uW^2-pW)/dx*dt-(rouN*uN*vN-rouS*uS*vS)/dy*dt)/max(rou(i,j),1e-4);
            v(i,j) = (rou_n(i,j)*v_n(i,j)-(rouE*uE*vE-rouW*uW*vW)/dx*dt-(rouN*vN^2+pN-rouS*vS^2-pS)/dy*dt)/max(rou(i,j),1e-4);
            E(i,j) = (rou_n(i,j)*E_n(i,j)-(rouE*uE*HE-rouW*uW*HW)/dx*dt-(rouN*vN*HN-rouS*vS*HS)/dy*dt)/max(rou(i,j),1e-4);
            p(i,j) = (E(i,j)-0.5*(u(i,j)^2+v(i,j)^2))*rou(i,j)*(gama-1);
            H(i,j) = E(i,j)+p(i,j)/max(rou(i,j),1e-4);            
            
        end
    end
    
  
    eps = max([norm(u-u_n,1)/norm(u,1),norm(v-v_n,1)/norm(v,1),norm(rou-rou_n,1)/norm(rou,1),norm(p-p_n,1)/norm(p,1)])
end
if eps<epsmax
    disp('done!');
    save 'data_no6.mat'
    pcolor(X,Y,u),colorbar;
end






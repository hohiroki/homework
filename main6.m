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
Uref = (Pref/Rouref)^0.5;
Lref = 5e-3;
eps = 1;
epsmax = 1e-6;
dt = 1/(xnum*ynum*10);
%入口无穷远条件
Main = 0.4357;
Tin = T0/(1+(gama-1)/2*Main^2); %静温
pin = P0/(1+(gama-1)/2*Main^2)^(gama/(gama-1));%静压
rouin = pin/R/Tin;  %密度
ain = (gama*pin/rouin)^0.5;%声速
uin = Main*ain; %速度
vin = 0;
Ein = pin/(gama-1)/rouin+0.5*uin^2;  %总能
Hin = Ein+pin/rouin;
%无量纲化入口无穷远条件：
pin_unit = pin/Pref;ain_unit = ain/Uref;uin_unit = ain_unit*Main;vin_unit = 0;rouin_unit = rouin/Rouref;
Ein_unit = pin_unit/(gama-1)/rouin_unit+0.5*uin_unit^2;Hin_unit = Ein_unit+pin_unit/rouin_unit;
S_unit = pin_unit/(rouin_unit)^gama;%无量纲商
%初始化
rou_init = ones(xnum,ynum)*rouin_unit;
p_init = ones(xnum,ynum)*pin_unit;
u_init = ones(xnum,ynum)*0;
v_init = ones(xnum,ynum)*0;
E_init = ones(xnum,ynum)*Ein_unit;
H_init = ones(xnum,ynum)*Hin_unit;
a_init = (gama*p_init./rou_init).^0.5;
%初始化n+1时刻，n时刻计算变量
rou = rou_init;
p = p_init;
u = u_init;
v = v_init;
E = E_init;
H = H_init;
a = (gama*p./rou).^0.5;
rou_n = rou;
u_n = u;
v_n = v;
E_n = E;
H_n = H;
a_n = a;
step = 0;
%单位行向量
ex = linspace(1,1,length(u(1,:)));
ey = linspace(1,1,length(u(:,1)))';
for time = 0:dt:1
    if eps<epsmax
        break
    end
    step =step+1;
    if step ==4
        break
    end
    %入口条件计算
    u(1,:) = (uin_unit-2*ain_unit*ex/(gama-1)+u(2,:)+2*a(2,:)/(gama-1))/2;
    v(1,:) = 0;
    a(1,:) = (u(2,:)+2*a(2,:)/(gama-1)-uin_unit*ey+2*ain_unit*ey/(gama-1))*(gama-1)/4;
    rou(1,:) = ((a(1,:).^2/gama/S_unit)).^(1/(gama-1));
    p(1,:) = rou(1,:).*a(1,:).^2/gama;
    E(1,:) = p(1,:)/(gama-1)./rou(1,:)+0.5*u(1,:).^2;  
    H(1,:) = E(1,:)+p(1,:)./rou(1,:);
    %壁面
    for i = 1:1:xnum-1
        [nx,ny] = getnxny([X(i,1),X(i+1,1)],[Y(i,1),Y(i+1,1)]);
        u(i,1) = u(i,2)*ny^2-v(i,2)*nx*ny;
        v(i,1) = -u(i,2)*nx*ny+v(i,2)*nx^2;
        a_temp = (p(i,2)/rou(i,2))^0.5;
        rou(i,1) = ((gama-1)^2*(u(i,2)*nx+v(i,2)*ny-2*a_temp/(gama-1))^2*rou(i,2)^gama/4/gama/p(i,2))^(1/(gama-1));
        p(i,1) = p(i,2)*rou(i,1)^gama/rou(i,2)^gama;
        E(i,1) = p(i,1)/(gama-1)/rou(i,1)+0.5*u(i,1)^2;  
        H(i,1) = E(i,1)+p(i,1)/rou(i,1);
        
        [nx,ny] = getnxny([X(i,ynum),X(i+1,ynum)],[Y(i,ynum),Y(i+1,ynum)]);
        nx = -1*nx;ny = -1*ny;
        u(i,ynum) = u(i,ynum-1)*ny^2-v(i,ynum-1)*nx*ny;
        v(i,ynum) = -u(i,ynum-1)*nx*ny+v(i,ynum-1)*nx^2;
        a_temp = (p(i,ynum-1)/rou(i,ynum-1))^0.5;
        rou(i,ynum) = ((gama-1)^2*(u(i,ynum-1)*nx+v(i,ynum-1)*ny-2*a_temp/(gama-1))^2*rou(i,ynum-1)^gama/4/gama/p(i,2))^(1/(gama-1));
        p(i,ynum) = p(i,ynum-1)*rou(i,ynum)^gama/rou(i,ynum-1)^gama;
        E(i,ynum) = p(i,ynum)/(gama-1)/rou(i,ynum)+0.5*u(i,ynum)^2;  
        H(i,ynum) = E(i,ynum)+p(i,ynum)/rou(i,ynum);
    end
    %内部节点
    for i = 2:1:xnum-1
        for j = 2:1:ynum-1
            %界面面积
            area = getArea([X(i-1,j),X(i,j-1),X(i+1,j),X(i,j+1)],[Y(i-1,j),Y(i,j-1),Y(i+1,j),Y(i,j+1)]);
            %界面通量计算
            %i+1/2,j处通量
            x_temp_v = [(X(i,j)+X(i+1,j)+X(i+1,j+1)+X(i,j+1))/4,(X(i,j)+X(i,j-1)+X(i+1,j-1)+X(i+1,j))/4];
            y_temp_v = [(Y(i,j)+Y(i+1,j)+Y(i+1,j+1)+Y(i,j+1))/4,(Y(i,j)+Y(i,j-1)+Y(i+1,j-1)+Y(i+1,j))/4];
            dl = ((x_tempv(1)-x_tempv(0))^2+(y_tempv(1)-y_tempv(0))^2)^0.5;
            [nx,ny] = getnxny(x_temp_v,y_temp_v);
            %密度加权量
            rouL = sqrt(rou_n(i,j));
            rouR = sqrt(rou_n(i+1,j)); 
            rouE=rouL*rouR;
            uE=(u_n(i,j)*rouL+u_n(i+1,j)*rouR)/(rouL+rouR);
            vE=(v_n(i,j)*rouL+v_n(i+1,j)*rouR)/(rouL+rouR);
            HE=(H_n(i,j)*rouL+H_n(i+1,j)*rouR)/(rouL+rouR);
            aE = ((gama-1)*(HE-0.5*(uE^2+vE^2)))^0.5;
            pE=(HE-0.5*(uE^2+vE^2))*rouE*(gama-1)/gama;
 
            
            FL=[rou(i,j)*u(i,j),rou(i,j)*u(i,j)^2+p(i,j),rou(i,j)*u(i,j)*v(i,j),rou(i,j)*u(i,j)*H(i,j)]';
            FR=[rou(i+1,j)*u(i+1,j),rou(i+1,j)*u(i+1,j)^2+p(i+1,j),rou(i+1,j)*u(i+1,j)*v(i+1,j),rou(i+1,j)*u(i+1,j)*H(i+1,j)]';
            GL = [rou(i,j)*v(i,j),rou(i,j)*u(i,j)*v(i,j),rou(i,j)*v(i,j)^2+p(i,j),rou(i,j)*v(i,j)*H(i,j)]';
            GR = [rou(i+1,j)*v(i+1,j),rou(i+1,j)*u(i+1,j)*v(i+1,j),rou(i+1,j)*v(i+1,j)^2+p(i+1,j),rou(i+1,j)*v(i+1,j)*H(i+1,j)]';
            F_langta = abs([uE-aE,uE,uE,uE+aE]);
            G_langta = abs([vE-aE,vE,vE,vE+aE]);
            F_afa = [(p(i+1,j)-p(i,j)-rouE*aE*(u(i+1,j)-u(i,j)))/2/aE^2,rou(i+1,j)-rou(i,j)-(p(i+1,j)-p(i,j))/aE^2,rouE*(v(i+1,j)-v(i,j)),(p(i+1,j)-p(i,j)+rouE*aE*(u(i+1,j)-u(i,j)))/2/aE^2];
            F_E = 0.5*(FL+FR)-0.5* 
            
            
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
            v(i,j) = (rou_n(i,j)*v_n(i,j)-(rouE*uE*vE-rouW*uW*vE)/dx*dt-(rouN*vN^2+pN-rouS*vS^2-pS)/dy*dt)/max(rou(i,j),1e-4);
            E(i,j) = (rou_n(i,j)*E_n(i,j)-(rouE*uE*HE-rouW*uW*HW)/dx*dt-(rouN*vN*HN-rouS*vS*HS)/dy*dt)/max(rou(i,j),1e-4);
            p(i,j) = (E(i,j)-0.5*(u(i,j)^2+v(i,j)^2))*rou(i,j)*(gama-1);
            H(i,j) = E(i,j)+p(i,j)/max(rou(i,j),1e-4);            
            
        end
    end
    %出口直接外推
    u(xnum,:) = u(xnum-1,:);v(xnum,:) = v(xnum-1,:);a(xnum,:) = a(xnum-1,:);rou(xnum,:) = rou(xnum-1,:);
    p(xnum,:) = p(xnum-1,:);E(xnum,:) = E(xnum-1,:);H(xnum,:) = H(xnum-1,:);
    %残差
    eps = max([norm(u-u_n,1)/norm(u,1),norm(v-v_n,1)/norm(v,1),norm(rou-rou_n,1)/norm(rou,1),norm(p-p_n,1)/norm(p,1)])
    %存储
    rou_n = rou;
    u_n = u;
    v_n = v;
    E_n = E;
    H_n = H;
    a_n = a;
end
if eps<epsmax
    disp('done!');
    save 'data_no6.mat'
    pcolor(X,Y,u),colorbar;
end






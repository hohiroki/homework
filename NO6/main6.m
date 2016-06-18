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
dt = 1/(xnum*ynum);
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
u_init = ones(xnum,ynum)*uin_unit*0.5;
v_init = ones(xnum,ynum)*uin_unit*0.01;
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
p_n = p;
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
    if step ==20
        break
    end
    %入口条件计算
    u(1,:) = (uin_unit-2*ain_unit*ex/(gama-1)+u(2,:)+2*a(2,:)/(gama-1))/2;
    v(1,:) = 0;
    a(1,:) = (u(2,:)+2*a(2,:)/(gama-1)-uin_unit*ex+2*ain_unit*ex/(gama-1))*(gama-1)/4;
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
        rou(i,ynum) = ((gama-1)^2*(u(i,ynum-1)*nx+v(i,ynum-1)*ny-2*a_temp/(gama-1))^2*rou(i,ynum-1)^gama/4/gama/p(i,ynum-1))^(1/(gama-1));
        p(i,ynum) = p(i,ynum-1)*rou(i,ynum)^gama/rou(i,ynum-1)^gama;
        E(i,ynum) = p(i,ynum)/(gama-1)/rou(i,ynum)+0.5*u(i,ynum)^2;  
        H(i,ynum) = E(i,ynum)+p(i,ynum)/rou(i,ynum);
    end
    %内部节点
    for i = 2:1:xnum-1
        for j = 2:1:ynum-1
            %界面面积
            area = getArea([X(i-1,j),X(i,j-1),X(i+1,j),X(i,j+1)],[Y(i-1,j),Y(i,j-1),Y(i+1,j),Y(i,j+1)])/Lref/Lref;
            %界面角点坐标
            x1=(X(i,j)+X(i+1,j)+X(i+1,j+1)+X(i,j+1))/4;
            y1=(Y(i,j)+Y(i+1,j)+Y(i+1,j+1)+Y(i,j+1))/4;
            x2=(X(i,j)+X(i-1,j)+X(i-1,j+1)+X(i,j+1))/4;
            y2=(Y(i,j)+Y(i-1,j)+Y(i-1,j+1)+Y(i,j+1))/4;
            x3=(X(i,j)+X(i,j-1)+X(i-1,j-1)+X(i-1,j))/4;
            y3=(Y(i,j)+Y(i,j-1)+Y(i-1,j-1)+Y(i-1,j))/4;
            x4=(X(i,j)+X(i,j-1)+X(i+1,j-1)+X(i+1,j))/4;
            y4=(Y(i,j)+Y(i,j-1)+Y(i+1,j-1)+Y(i+1,j))/4;
            %界面通量计算
            %i+1/2,j处通量
            x_temp_v = [x1,x4];
            y_temp_v = [y1,y4];
            dlE = ((x_temp_v(2)-x_temp_v(1))^2+(y_temp_v(2)-y_temp_v(1))^2)^0.5/Lref;
            [nxE,nyE] = getnxny(x_temp_v,y_temp_v);
            %密度加权量
            [rouE,uE,vE,HE,aE,pE]=get_WRou(rou_n,u_n,v_n,H_n,i,j,i+1,j);
            [F_E,G_E] = get_FGrebuilt(rou_n,u_n,v_n,H_n,p_n,rouE,uE,vE,HE,aE,i,j,i+1,j);      
            H_passE=(nxE*F_E+nyE*G_E)*dlE;        
            
            %i-1/2,j处通量
            x_temp_v = [x3,x2];
            y_temp_v = [y3,y2];
            dlW = ((x_temp_v(2)-x_temp_v(1))^2+(y_temp_v(2)-y_temp_v(1))^2)^0.5/Lref;
            [nxW,nyW] = getnxny(x_temp_v,y_temp_v);
            [rouW,uW,vW,HW,aW,pW]=get_WRou(rou_n,u_n,v_n,H_n,i-1,j,i,j);
            [F_W,G_W] = get_FGrebuilt(rou_n,u_n,v_n,H_n,p_n,rouW,uW,vW,HW,aW,i-1,j,i,j);      
            H_passW=(nxW*F_W+nyW*G_W)*dlW;
            
            %i,j-1/2处通量
            x_temp_v = [x4,x3];
            y_temp_v = [y4,y3];
            dlS = ((x_temp_v(2)-x_temp_v(1))^2+(y_temp_v(2)-y_temp_v(1))^2)^0.5/Lref;
            [nxS,nyS] = getnxny(x_temp_v,y_temp_v);
            [rouS,uS,vS,HS,aS,pS]=get_WRou(rou_n,u_n,v_n,H_n,i,j-1,i,j);
            [F_S,G_S] = get_FGrebuilt(rou_n,u_n,v_n,H_n,p_n,rouS,uS,vS,HS,aS,i,j-1,i,j);      
            H_passS=(nxS*F_S+nyS*G_S)*dlS;
            
            %i,j+1/2处通量
            x_temp_v = [x2,x1];
            y_temp_v = [y2,y1];
            dlN = ((x_temp_v(2)-x_temp_v(1))^2+(y_temp_v(2)-y_temp_v(1))^2)^0.5/Lref;
            [nxN,nyN] = getnxny(x_temp_v,y_temp_v);
            [rouN,uN,vN,HN,aN,pN]=get_WRou(rou_n,u_n,v_n,H_n,i,j,i,j+1);
            [F_N,G_N] = get_FGrebuilt(rou_n,u_n,v_n,H_n,p_n,rouN,uN,vN,HN,aN,i,j,i,j+1);      
            H_passN=(nxN*F_N+nyN*G_N)*dlN;
            
            %一阶时间离散
            du_dt = dt*(H_passW+H_passE+H_passN+H_passS)/area;
            rou(i,j) = max(rou_n(i,j)-du_dt(1),1e-4);
            u(i,j) = (u_n(i,j)*rou_n(i,j)-du_dt(2))/du_dt(1);
            v(i,j) = (v_n(i,j)*rou_n(i,j)-du_dt(3))/du_dt(1);
            E(i,j) = (E_n(i,j)*rou_n(i,j)-du_dt(4))/du_dt(1);
            p(i,j) = (E(i,j)-0.5*(u(i,j)^2+v(i,j)^2))*rou(i,j)*(gama-1);
            H(i,j) = E(i,j)+p(i,j)/rou(i,j);            
            
        end
    end
    %出口直接外推
    u(xnum,:) = u(xnum-1,:);v(xnum,:) = v(xnum-1,:);a(xnum,:) = a(xnum-1,:);rou(xnum,:) = rou(xnum-1,:);
    p(xnum,:) = p(xnum-1,:);E(xnum,:) = E(xnum-1,:);H(xnum,:) = H(xnum-1,:);
    %残差
    eps = max([norm(u-u_n,1)/norm(u,1),norm(v-v_n,1)/norm(v,1),norm(rou-rou_n,1)/norm(rou,1),norm(E-E_n,1)/norm(E,1)])
    %存储
    rou_n = rou;
    u_n = u;
    v_n = v;
    E_n = E;
    H_n = H;
    a_n = a;
    p_n = p;
end
if eps<epsmax
    disp('done!');
    save 'data_no6.mat'
    pcolor(X,Y,u),colorbar;
end






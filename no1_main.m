clear all;
clc;
%format long;
%general contral
delta_max = 1e-8;
step_max = 600;

%boudary setting
U = 0.002;
Lx = 0.5;
Ly = 0.1;
mu = 0.001;
rou = 998;
%girid setting
xgrid_num = 100;
ygrid_num = 100;
% dx_averge = Lx/xgrid_num;
% dy_averge = Ly/ygrid_num;
% x_min = dx_averge/10;
% y_min = dy_averge/20;
% Lx_min = Lx*0.2;
% Ly_min = Ly*0.3;
% x_min_grid = xgrid_num*0.4;
% y_min_grid = ygrid_num*0.8;
% 
% dx_normal = (Lx-Lx_min)/(xgrid_num-x_min_grid);
% dy_normal = (Ly-Ly_min)/(ygrid_num-y_min_grid);
% x_ratio = (dx_normal/x_min)^(1/(x_min_grid-1));
% y_ratio = (dy_normal/y_min)^(1/(y_min_grid-1));
% 
% dx_grid = linspace(1,1,xgrid_num)*dx_normal;
% dy_grid = linspace(1,1,ygrid_num)*dy_normal;
% 
% x_grid = linspace(0,0,xgrid_num+1);
% y_grid = linspace(0,0,ygrid_num+1);
% 
% dx_grid(1) = x_min;
% dy_grid(1) = y_min;
% for i =2:1:x_min_grid
%     dx_grid(i) = x_min * x_ratio^(i - 1);
%     x_grid(i) = x_grid(i-1)+dx_grid(i-1);
% end
% for i =x_min_grid+1:1:xgrid_num+1
%     x_grid(i) = x_grid(i-1)+dx_grid(i-1);
% end
% for i =2:1:y_min_grid
%     dy_grid(i) = y_min * y_ratio^(i - 1);
%     y_grid(i) = y_grid(i-1)+dy_grid(i-1);
% end
% for i =y_min_grid+1:1:ygrid_num+1
%     y_grid(i) = y_grid(i-1)+dy_grid(i-1);
% end
%%%%%-------%%%%%%
dx_grid = linspace(1,1,xgrid_num+1)*Lx/xgrid_num;
dy_grid = linspace(1,1,ygrid_num+1)*Ly/ygrid_num;
x_grid = linspace(0,0,xgrid_num+1);
y_grid = linspace(0,0,ygrid_num+1);
x_grid(1)=0;
y_grid(1)=0;
for i=2:1:xgrid_num+1
    x_grid(i) =  x_grid(i-1)+dx_grid(i-1);
end
for i=2:1:ygrid_num+1
    y_grid(i) =  y_grid(i-1)+dy_grid(i-1);
end
%%%%%-------%%%%%%%%
%show grid
% [X,Y] = meshgrid(x_grid,y_grid);
% Z = X.*0;
% surf(X,Y,Z);

%caculation
%initial
u=ones(xgrid_num+1,ygrid_num+1)*U;
v=zeros(xgrid_num+1,ygrid_num+1);
u0 = u;
v0 = v;
delta = 1;
step = 0;

    %add visual point
dy_grid2=[dy_grid,dy_grid(ygrid_num)];
dx_grid2=[dx_grid,dx_grid(ygrid_num)];

%pE&mE
while delta>delta_max&&step<step_max
    step = step+1;
    u0 = u;
    v0 = v;
    %boundary
    u(1,:) = U;
    v(1,:) = 0;
    u(:,1) = 0;
    v(:,1) = 0;
    aE2=u;
    aW2 = u;
    aN2 = u; 
    aS2 = u;
    aP2 = u;
    
    %for upwind
    for i=2:1:xgrid_num
        for j=2:1:ygrid_num
            Fe = rou*(u(i,j)+u(i+1,j))/2;
            Fw = rou*(u(i-1,j)+u(i,j))/2;
            Fn = rou*(v(i,j)+v(i,j+1))/2;
            Fs = rou*(v(i,j-1)+v(i,j))/2;
            
            Dn = mu/dy_grid2(j);
            Ds = mu/dy_grid2(j-1);
            Pn = Fn/Dn;
            Ps = Fs/Ds;
            
            %An = 1+max([-Pn,0]);
            %As = 1+max([-Ps,0]);
            As = max([0,1-0.5*abs(Ps)]);
            Bs = As+Ps;
            An = max([0,1-0.5*abs(Pn)]);
            
            
            dx = (dx_grid2(i)+dx_grid2(i-1))/2.0;
            dy = (dy_grid2(j)+dy_grid2(j-1))/2.0;
            
            
            aE = max([-Fe,0])*dy;
            aW = max([Fw,0])*dy;
            aN = Dn*An*dx;
            %aN = 0; 
            aS = Ds*Bs*dx;
            aP = aE+aW+aN+aS+(Fe-Fw)*dy+(Fn-Fs)*dx;
            u(i,j)=(aE*u(i+1,j)+aW*u(i-1,j)+aN*u(i,j+1)+aS*u(i,j-1))/aP;         
            aE2(i,j)=aE;
            aW2(i,j) = aW;
            aN2(i,j) = aN; 
            aS2(i,j) = aS;
            aP2(i,j) = aP;
        end
    end
    %v方向速度，连续方程给出
    for i=2:1:xgrid_num
        for j=2:1:ygrid_num
            v(i,j+1) = v(i,j-1)-dy*(u(i+1,j)-u(i-1,j))/dx;
        end
    end
    %at x=xgrid_num+1;
        %充分发展
    u(xgrid_num+1,:)=u(xgrid_num,:);
    v(xgrid_num+1,:)=v(xgrid_num,:);
    %at y=ygrid_num+1
        %无速度梯度
    v(:,ygrid_num+1)=v(:,ygrid_num);
    u(:,ygrid_num+1)=u(:,ygrid_num);
        %排挤速度
%     for i = 2:1:xgrid_num+1
%         v(i,ygrid_num+1) = ((1-u(i,:)/U)*dy_grid2'-(1-u(i-1,:)/U)*dy_grid2')/dx_grid(i-1);
%        
%     end
        %连续方程
%     for i = 2:1:xgrid_num+1
%           dx = (dx_grid2(i)+dx_grid2(i-1))/2;
%           dy = (dy_grid2(xgrid_num+1)+dy_grid2(xgrid_num))/2;
%         v(i,ygrid_num+1) =  v(i,ygrid_num)-dy*(u(i,ygrid_num+1)-u(i-1,ygrid_num+1))/dx; 
%     end         
    delta = max(max(abs(u(2:xgrid_num,2:ygrid_num)-u0(2:xgrid_num,2:ygrid_num))))/U+max(max(abs(v(2:xgrid_num,2:ygrid_num)-v0(2:xgrid_num,2:ygrid_num))))/U
end
%
if step==step_max
    disp('step over maxstep');
    disp(delta);
    result = 1;
    return
    
else
    disp('done,delta=');
    disp(double(delta));
    result = [u,v];
    save 'data.mat';
    [X,Y] = meshgrid(x_grid,y_grid);
    contourf(X,Y,u',15);
    shading flat;
end

% function A = get_A_up(P)
% A = 1;
% end
% function A = get_A_exp(P)
% A = abs(P)/(exp(abs(P))-1);
% end
% function A= get_A_mix(P)
% A = max([0,1-0.5*abs(P)]);
% end
% function A = get_A_eq(P)
% A = max([0,(1-0.1*abs(P))^5]);
% end
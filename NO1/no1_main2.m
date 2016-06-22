clear all;
clc;
%format long;
%general contral
delta_max = 1e-5;
step_max = 1800;

%boudary setting
U = 0.002;
Lx = 0.5;
Ly = 0.1;
mu = 1e-3;
rou = 998;
Cmu = 0.09;
Ce1 = 1.45;
Ce2 = 1.9;
ok = 1.0;
oe = 1.3;
%girid setting
xgrid_num = 100;
ygrid_num = 100;
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
%initial
u=ones(xgrid_num+1,ygrid_num+1)*U;
v=zeros(xgrid_num+1,ygrid_num+1);
u0 = u;
v0 = v;
delta = 1;
step = 0;

%%%%%
k = u.*0;
e =k.*0;
mut = k.*0;
u0 = u;
v0 = v;
k0=k;
e0=e;
mut0 = mut;
%%%%%
%add visual point
dy_grid2=[dy_grid,dy_grid(ygrid_num)];
dx_grid2=[dx_grid,dx_grid(ygrid_num)];

%pE&mE
while delta>delta_max&&step<step_max
    step = step+1;
    u0 = u;
    v0 = v;
    k0=k;
    e0=e;
    mut0 = mut;
    
    %boundary
    u(1,:) = U;
    v(1,:) = 0;
    u(:,1) = 0;
    v(:,1) = 0;
    %%%%%%
    k(1,4:10) = U*U*0.002;
    e(1,4:10) = 2*k(1,4:10).^1.5*Cmu^0.75/0.095/0.42;
    u(:,1) = 0;
    v(:,1) = 0;
    k(:,1) = 0;
    e(:,1) = 0; 
    %u
    for i=2:1:xgrid_num
        for j=2:1:ygrid_num
            Fe = rou*(u(i,j)+u(i+1,j))/2;
            Fw = rou*(u(i-1,j)+u(i,j))/2;
            Fn = rou*(v(i,j)+v(i,j+1))/2;
            Fs = rou*(v(i,j-1)+v(i,j))/2;
            
            Dn = (mu+mut(i,j))/dy_grid2(j);
            Ds = (mu+mut(i,j))/dy_grid2(j-1);
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
        end
    end
    %v方向速度，连续方程给出
    u(xgrid_num+1,:)=u(xgrid_num,:);
    u(:,ygrid_num+1)=u(:,ygrid_num);
    for i=2:1:xgrid_num
        for j=2:1:ygrid_num
            v(i,j+1) = v(i,j-1)-dy*(u(i+1,j)-u(i-1,j))/dx;
        end
    end
    %at x=xgrid_num+1;
    v(xgrid_num+1,:)=v(xgrid_num,:);
    %at y=ygrid_num+1
    v(:,ygrid_num+1)=v(:,ygrid_num);

    %%%%%k
    for i=2:1:xgrid_num
        for j=2:1:ygrid_num
            Fe = rou*(u(i,j)+u(i+1,j))/2;
            Fw = rou*(u(i-1,j)+u(i,j))/2;
            Fn = rou*(v(i,j)+v(i,j+1))/2;
            Fs = rou*(v(i,j-1)+v(i,j))/2;
            dx = (dx_grid2(i)+dx_grid2(i-1))/2.0;
            dy = (dy_grid2(j)+dy_grid2(j-1))/2.0;
                        
            aE = max([-Fe,0])*dy;
            aW = max([Fw,0])*dy;
            aN = Dn*An*dx; 
            aS = Ds*Bs*dx;
            aP = aE+aW+aN+aS+(Fe-Fw)*dy+(Fn-Fs)*dx;
            %k PE
            G = rou*mut(i,j)*(2*(u(i,j)-u(i-1,j))^2/dx^2+2*(v(i,j)-v(i-1,j))^2/dy^2+...
            ((u(i,j)-u(i-1,j))/dx+(v(i,j)-v(i-1,j))/dy)^2);
            Sc_k = G-mu*e(i,j)+k0(i,j);
            Sp_k = -1;
            
            Dn = (mu+mut(i,j)/ok)/dy_grid2(j);
            Ds = (mu+mut(i,j)/ok)/dy_grid2(j-1);
            Pn = Fn/Dn;
            Ps = Fs/Ds;
            
            As = max([0,1-0.5*abs(Ps)]);
            Bs = As+Ps;
            An = max([0,1-0.5*abs(Pn)]);
                        
            aE = max([-Fe,0])*dy;
            aW = max([Fw,0])*dy;
            aN = Dn*An*dx; 
            aS = Ds*Bs*dx;
            
            aP_k = aE+aW+aN+aS+(Fe-Fw)*dy+(Fn-Fs)*dx-Sp_k*dx*dy;
            k(i,j) = (aE*k(i+1,j)+aW*k(i-1,j)+aN*k(i,j+1)+aS*k(i,j-1)+Sc_k*dx*dy)/aP_k;
                       
        end
    end
    k(xgrid_num+1,:)=k(xgrid_num,:);
    k(:,ygrid_num+1)=k(:,ygrid_num);
%%%%% e
    for i=2:1:xgrid_num
        for j=2:1:ygrid_num
            Fe = rou*(u(i,j)+u(i+1,j))/2;
            Fw = rou*(u(i-1,j)+u(i,j))/2;
            Fn = rou*(v(i,j)+v(i,j+1))/2;
            Fs = rou*(v(i,j-1)+v(i,j))/2;
            dx = (dx_grid2(i)+dx_grid2(i-1))/2.0;
            dy = (dy_grid2(j)+dy_grid2(j-1))/2.0;
                        
            G = rou*mut(i,j)*(2*(u(i,j)-u(i-1,j))^2/dx^2+2*(v(i,j)-v(i-1,j))^2/dy^2+...
            ((u(i,j)-u(i-1,j))/dx+(v(i,j)-v(i-1,j))/dy)^2);
            Dn = (mu+mut(i,j)/oe)/dy_grid2(j);
            Ds = (mu+mut(i,j)/oe)/dy_grid2(j-1);
            Pn = Fn/Dn;
            Ps = Fs/Ds;
            
            As = max([0,1-0.5*abs(Ps)]);
            Bs = As+Ps;
            An = max([0,1-0.5*abs(Pn)]);
                        
            aE = max([-Fe,0])*dy;
            aW = max([Fw,0])*dy;
            aN = Dn*An*dx; 
            aS = Ds*Bs*dx;
            Sc_e = Ce1*G*e0(i,j)/k(i,j);
            Sp_e = -rou*Ce2*e0(i,j)/k(i,j);
            aP_e = aE+aW+aN+aS+(Fe-Fw)*dy+(Fn-Fs)*dx-Sp_e*dx*dy;
            e(i,j) = (aE*e(i+1,j)+aW*e(i-1,j)+aN*e(i,j+1)+aS*e(i,j-1)+Sc_e*dx*dy)/aP_e;                  
        end
    end
    e(xgrid_num+1,:)=e(xgrid_num,:);
    e(:,ygrid_num+1)=e(:,ygrid_num);
    for i=1:1:xgrid_num+1
        for j=1:1:ygrid_num+1
            if k(i,j)==0||e(i,j)==0
                mut(i,j) =0;
            else
                mut(i,j)= Cmu*k(i,j)^2./e(i,j);
            end
        end
    end
    %%%%%
        
    delta = max(max(abs(u(2:xgrid_num,2:ygrid_num)-u0(2:xgrid_num,2:ygrid_num))))/U+max(max(abs(v(2:xgrid_num,2:ygrid_num)-v0(2:xgrid_num,2:ygrid_num))))/U
end
%
if step==step_max
    disp('step over maxstep');
    disp(delta);
    result = 1;
    
% else
%     disp('done,delta=');
%     disp(double(delta));
%     result = [u,v];
%     save 'data.mat';
%     [X,Y] = meshgrid(x_grid,y_grid);
%     contourf(X,Y,u',15);
%     shading flat;
end


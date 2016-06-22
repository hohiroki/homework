%边界层后处理
function get_fig(x_grid,y_grid,u,v)
[X,Y] = meshgrid(x_grid,y_grid);
    contourf(X,Y,u',15);
    axis equal
    set(gca,'XLim',[0 0.5]);
    set(gca,'YLim',[0 0.1]);
    set(gca,'XTick',[0:0.05:0.5]);
    set(gca,'YTick',[0:0.025:0.1]);
    shading flat;
    hcb=colorbar;
    set(hcb,'ytick',[0:0.0004:0.002]);
    
 contourf(X,Y,v',15);
    axis equal
    set(gca,'XLim',[0 0.5]);
    set(gca,'YLim',[0 0.1]);
    set(gca,'XTick',[0:0.05:0.5]);
    set(gca,'YTick',[0:0.025:0.1]);
    shading flat;
    hcb=colorbar;
    set(hcb,'ytick',[0:max(max(v))/6:max(max(v))]);
    
  y_u99=linspace(0,0,xgrid_num+1);
  
    for i=2:1:xgrid_num+1
      for j=1:1:ygrid_num
          if u(i,j)<=0.99*U &&u(i,j+1)>=0.99*U
              y_u99(i)=y_grid(j);
              break
          end
      end
    end
  y_ref = 4.91*(mu*x_grid/rou/U).^0.5;
  plot(x_grid,y_u99);hold on
  plot(x_grid,y_ref,'Color','g');
  axis equal
  set(gca,'XLim',[0 0.5]);
  set(gca,'YLim',[0 0.1]);
  set(gca,'XTick',[0:0.05:0.5]);
  set(gca,'YTick',[0:0.025:0.1]);
  legend('数值解','布拉修斯解');
end
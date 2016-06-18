function [nx,ny] = getnxny(x,y)
if ~isvector(x)
    error('x is not vector')
end
if length(x)~=2
    error('x length is not 2')
end
if ~isvector(y)
    error('y is not vector')
end
if length(y)~=2
    error('y length is not 2')
end
v = [x(2)-x(1),y(2)-y(1)];
n = null(v);
nx = n(1);
ny = n(2);
if v(1)>0&&v(2)>0
    ny=max(ny,-ny);
    nx=min(nx,-nx);
elseif v(1)<0&&v(2)>0
    ny=min(ny,-ny);
    nx=max(nx,-nx);
elseif v(1)<0&&v(2)<0
    ny=min(ny,-ny);
    nx=max(nx,-nx);
elseif v(1)>0&&v(2)<0
    ny=max(ny,-ny);
    nx=max(nx,-nx);
elseif v(1)==0&&v(2)~=0
    ny=0;
    nx=-v(2)/abs(v(2));
elseif v(2)==0&&v(1)~=0
    ny=v(1)/abs(v(1));
    nx=0;
end
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
end
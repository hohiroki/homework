function area = getArea(x,y)
if ~isvector(x)
    error('x is not vector')
end
if length(x)~=4
    error('x length is not 4')
end
if ~isvector(y)
    error('y is not vector')
end
if length(y)~=4
    error('y length is not 4')
end
area = 0;
for i=1:1:3
    area =area + x(i)*y(i+1)-x(i+1)*y(i);
end
area = 0.5*abs(area);
end
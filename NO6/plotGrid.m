function plotGrid(X,Y)
if size(X)~=size(Y)
    error('X and Y not match')
end
for j=1:1:length(Y(1,:))    
   plot(X(:,j),Y(:,j));hold on; 
end
for i = 1:1:length(X(:,1))
    plot(X(i,:),Y(i,:));hold on;
end
box off;
end
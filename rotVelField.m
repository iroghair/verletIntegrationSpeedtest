x = linspace(0,1,61);
y = x;
[X,Y] = meshgrid(x,y);
clear x,y;

ux = -2*(sin(pi*X))^2.*sin(pi*Y).*cos(pi*Y);
uy =  2*(sin(pi*Y))^2.*sin(pi*X).*cos(pi*X);
% ux = 1000*ones(size(ux));

quiver(X,Y,ux,uy);
axis equal tight
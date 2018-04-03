%% Question 4: divergence theorem getting a graph first 
fx = @(x,y,z) x + y.^1.*z.^1; 
fy = @(x,y,z) -3*x.*y.^1; 
fz = @(x,y,z) x.*z.^2; 

yrange = linspace(-4, 4, 20); 
zrange = linspace(-4, 4, 20); 
[yz, zy] = meshgrid(yrange, zrange); 
plane0 = @(y,z) zeros(size(y)); 
plane2 = @(y,z) 5*ones(size(y)); 

urange = linspace(0, 2*pi, 20); 
vrange = linspace(0, 5, 20); 
[uv, vu] = meshgrid(urange, vrange); 
cx = @(u,v) v; 
cy = @(u,v) 4.*cos(u); 
cz = @(u,v) 4.*sin(u); 


hold on; 
surf(cx(uv,vu), cy(uv,vu), cz(uv,vu)); 

x = cx(uv,vu); y = cy(uv,vu); z = cz(uv,vu); 

quiver3(x, y, z, fx(x,y,z), fy(x,y,z), fz(x,y,z), 3, 'b'); 
quiver3(plane0(yz,zy), yz, zy, fx(plane0(yz,zy),yz,zy), fy(plane0(yz,zy),yz,zy), fz(plane0(yz,zy),yz,zy), 2, 'b'); 
quiver3(plane2(yz,zy), yz, zy, fx(plane2(yz,zy),yz,zy), fy(plane2(yz,zy),yz,zy), fz(plane2(yz,zy),yz,zy), 2, 'b');  

hold off; 
xlabel('x'); 
ylabel('y'); 
zlabel('z'); 
axis([-0.5 5, -5 5, -4 4.5]); 
title('Ahmed Fuad Ali, 400075937');  

%% Question 5: Divergence theorem right side integral 
divF = @(x,y,z) 1 - (3.*x) + (2.*z.*x);
u = @(x,y) -sqrt(9-y.^2);
v = @(x,y) sqrt(9-y.^2); 
A = integral3(divF,0,5,-3,3,u,v)

%% Question 6: left hand side of the divergence theorem 
u2 = @(x,y,z) -sqrt(9-x.^2);
v2 = @(x,y,z) sqrt(9-x.^2);

divF1 = @(y,z) 5 + (y.^1).*(z.^2);
divF2 = @(y,z) 1*(y.^1).*(z.^2);

answer_1 = integral2(divF1,-4,4,u2,v2);
answer_2 = integral2(divF2,-4,4,u2,v2);
answer_3 = A - answer_1 - answer_2;

A2 = [answer_1, answer_2, answer_3]

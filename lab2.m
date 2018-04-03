%% Question 1

% create a contour plot of the given function, 100x100 grid
[x,y] = meshgrid(-1:0.02:1, -1:0.02:1);
f = @(x,y)3*(1-3*x).^2.*exp(-2*x.^2 - (3*y+1).^2)- 10*(3*x/5 - 27*x.^3 - 243*y.^5).*exp(-3*x.^2-5*y.^2)- 1/3*exp(-(3*x+1).^2 - 9*y.^2) + (x.^2+y.^2) - 1;
contourf(x,y,f(x,y),20);
colorbar;
title('Ahmed Fuad Ali, 400075937');
xlabel('x');
ylabel('y');


% now, also graph the vector field given by the gradient vector on the same
% graph, but plot this on a 20x20 grid
[x,y] = meshgrid(-1:2/20:1, -1:2/20:1); 
fx = @(x,y)-(18*(1-3*x)).*exp(-2*x.^2-(3*y+1).^2) ... 
- 6*(1-3*x).^2.*2.*x.*exp(-2*x.^2-(3*y+1).^2) ... 
- (10*(3/5-81*x.^2)).*exp(-3*x.^2-5*y.^2) ... 
+ (20*((3/5)*x-27*x.^3-243*y.^5)).*3.*x.*exp(-3*x.^2-5*y.^2) ... 
- (1/3*(-18*x-6)).*exp(-(3*x+1).^2-9*y.^2)+2*x; 

fy = @(x,y)3*(1-3*x).^2.*(-18*y-6).*exp(-2*x.^2-(3*y+1).^2) ... 
+ 12150*y.^4.*exp(-3*x.^2-5*y.^2) ... 
+ (20*((3/5)*x-27*x.^3-243*y.^5)).*5.*y.*exp(-3*x.^2-5*y.^2) ... 
+ 6*y.*exp(-(3*x+1).^2-9*y.^2)+2*y;

hold on;
quiver(x,y,fx(x,y),fy(x,y),'k')
hold off;

%% Question 2
% final x y and z coords of the hiker
x0 = [-0.18, 0.66];

rprime = @(t,x)[-fx(x(1),x(2));-fy(x(1),x(2))];
tDomain = [0:0.001:5.0]; % tdomain specified by question
[t r] = ode45(rprime,tDomain,x0);
px = r(end,1);
py = r(end,2); % end command allows us to get the final value
pz = f(px,py); % to get final values, just print these


%% Question 3
% plot the trajectory above, thickness =2, black for positive, white for
% negative 

rx = r(:,1); % get all the x values and y values of r
ry = r(:,2);

rpos = [];
rneg = [];
for i = 1:length(rx)
    if f(rx(i),ry(i))>=0
        rpos = [rpos;rx(i),ry(i)];
    else
        rneg = [rneg;rx(i),ry(i)];    
    end
end    

hold on;
plot(rpos(:,1),rpos(:,2),'k','LineWidth',2); % k is black w is white
plot(rneg(:,1),rneg(:,2),'w','LineWidth',2); % line thickness of 2
xlim([-1,1]);
ylim([-1,1]);
hold off;

%% Question 4
% calculate distance bw initial and final pos of hiker
distance = [px,py,pz] - [x0,f(x0(1),x0(2))];  % x0 is initial, px is final
distance = norm(distance); % to output these just print

%% Question 5 
% calculate elevation loss of hiker

elevLoss = f(x0(1), x0(2)) - pz; % to see this just print

%% Question 6
% calculate total distance
dx = 0;
dy = 0;
dz = 0;  % initialize the distances
totaldis = 0;
for i=1:(length(rx)-1) % minus one so we don't hit max
    dx = rx(i+1) - rx(i);
    dy = ry(i+1) - ry(i); % continually sum distances
    dz = f(rx(i+1),ry(i+1)) - f(rx(i),ry(i));
    totaldis = totaldis + norm([dx;dy;dz]);
end
totaldis; % to output these, just print




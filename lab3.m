%% Question 1
% draw a graph of elliptic torus
u = linspace(0, 2*pi, 35);
v = linspace(0, 2*pi, 35);
[u,v] = meshgrid(u,v);
x = (3+2.*cos(v)).*cos(u);
y = (3+2.*cos(v)).*sin(u);
z = 4.*sin(v);
surf(x,y,z);
title('Ahmed Fuad Ali, 400075937')

%% Question 2
% calculate surface area of the elliptic torus given above, and make sure
% to have the one side of the integral

% initialize new space so we can take derivatives
syms u v;

Rx = (6+2.*cos(v)).*cos(u);
Ry = (6+2.*cos(v)).*sin(u);
Rz = 4.*sin(v);

% find the cross product of so we can find normal vector
X = [diff(Rx,u), diff(Ry,u), diff(Rz,u)]; 
Y = [diff(Rx,v), diff(Ry,v), diff(Rz,v)]; 
product = cross(X,Y);

% find magnitude of that vector
mag = sqrt((product(1)^2 + product(2)^2 + product(3)^2)); 

% now get the integral
integrand = matlabFunction(mag);
c = 0;
d = 2*pi;
integral2(integrand,0,2*pi,c,d)

% this should be your answer so display it

%% Question 3
% evaluate right hand integral by first getting it (a) and then integrating
% it over the appropriate region

% first, clear our previous vector spaces and initialize a new space
clear x y 
syms x y 

% now get the new integrands by differentiating
Q = cos(x.^6).*log(y.^3.*exp(x)+1); 
P = exp(x.^2)/y.^3;
Qx = diff(Q,x); 
Py = diff(P,y); 
linearintegrand = Qx-Py

%3b 
issa2 = matlabFunction(linearintegrand);
a = -1; 
b = 1; 
c = 3; 
d = @(x) 4-x.^(2); 
integral2(issa2,a,b,c,d)


%% Question 4
% get the left side of Green's theorem instead this time
% you can check your answer by getting same sum bc c1 + c2 = above sum

f=@(x)exp(x.^2)/3^3; 
c1=quad(f,-1,1);
c1
c2=1.2138-c1; % by definition of integral, c2 has to be this
c2

%% Question 5
% plot the force field and contour function c 

x = -1:0.1:1;
y1 = 4-x.^2; 
y2 = 4 + (x.^2 - x.^2); 

plot(x,y1,'k','LineWidth',2);
hold on;
plot(x,y2,'k','LineWidth',2);
grid on

%first make the vectors we need for the semi circle
for i = 1:21
    F_xcom(i) = (6*x(i)^2 + 5*y1(i))*exp(sqrt(x(i)+1));
    F_ycom(i) = (3*x(i)^2 + x(i)*y1(i))*exp(sqrt(y1(i)));
end

% second, do it for the line 
for i = 1:21
    F_x2com(i) = (6*x(i)^2 + 6*y2(i))*exp(sqrt(x(i)+1));
    F_y2com(i) = (6*x(i)^2 + x(i)*y2(i))*exp(sqrt(y2(i)));
end

% plot the semicircle
quiver(x,y1,F_xcom,F_ycom,'m'); 
hold on;
% plot the line
quiver(x,y2,F_x2com,F_y2com,'g');


xlabel('x');
ylabel('y');
xlim([-1.5 1.5]); 
ylim([2.5 4.5]);
title('Ahmed Fuad Ali, 400075937');

%% Question 6
% find work done by force field from above quesiton
clear x y
syms x y 

P1 = (6*x.^(2)+5*y)*exp(sqrt(x+1)); 
Q1 = (3*x.^(2)+x*y)*exp(sqrt(y)); 
q1x = diff(Q1,x); 
p1y = diff(P1,y); 

% by stokes theorem, the integrand should be
integrand = q1x-p1y; 

% now take the integral with the appropriate bounds
issafunction = matlabFunction(integrand); 
a = -1; 
b = 1; 
c = 3; 
d = @(x) 4-x.^(2); 
integral2(issafunction,a,b,c,d);

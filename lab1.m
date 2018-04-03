%% PROBLEM # 1
% Make a vector function and plot the curve, turn in image
fileID = fopen('lab1.txt');
data = fscanf(fileID,'%u%g%g%g', [4 2000]);
fclose(fileID);

% file opened, so plot it
t = data(1,:);
ft = 4*data(2,:);
gt = 8*data(3,:);
ht = 5*data(4,:);
N = 2000;

plot3(ft,gt,ht,'Linewidth',2)
xlabel('x'), ylabel('y'), zlabel('z');
title('Ahmed Fuad Ali, 400075937');

%% PROBLEM #2
% Plot the coressponding tangent and normal vectors for all i
% use the forward difference quotients to calculate it

df = zeros(1, N-1);
dg = zeros(1, N-1);
dh = zeros(1, N-1);
tanx = zeros(1, N-1);
tany = zeros(1, N-1);
tanz = zeros(1, N-1);
rSpeed = zeros(1, N-1);
dtanx = zeros(1, N-1);
dtany = zeros(1, N-1);
dtanz = zeros(1, N-1);
tanSpeed = zeros(1, N-1);
normx = zeros(1, N-1);
normy = zeros(1, N-1);
normz = zeros(1, N-1);


% to calc tans, need r derivatives and r speed
for i=1:N-1
  df(i) = (ft(i+1) - ft(i))/(t(i+1)-t(i));
  dg(i) = (gt(i+1) - gt(i))/(t(i+1)-t(i));
  dh(i) = (ht(i+1) - ht(i))/(t(i+1)-t(i));
  rSpeed(i) = sqrt((df(i))^2+(dg(i))^2+(dh(i))^2);
  tanx(i) = df(i)/rSpeed(i);
  tany(i) = dg(i)/rSpeed(i);
  tanz(i) = dh(i)/rSpeed(i);
 
end

% Norms need tan derivatives and tanSpeed
for i=1:N-2
  dtanx(i) = (tanx(i+1) - tanx(i))/(t(i+1)-t(i));
  dtany(i) = (tany(i+1) - tany(i))/(t(i+1)-t(i));
  dtanz(i) = (tanz(i+1) - tanz(i))/(t(i+1)-t(i));
end 

% calculate unit normals
for i=1:N-1
  tanSpeed(i) = sqrt((dtanx(i))^2+(dtany(i))^2+(dtanz(i))^2);
  normx(i) = dtanx(i)/tanSpeed(i);
  normy(i) = dtany(i)/tanSpeed(i);
  normz(i) = dtanz(i)/tanSpeed(i);
end 

%% PROBLEM 3
% Plot the tangent vecotrs, unit normals, and binormal vectors along C
% using quiver for t = 50, 250, 450

Bx = zeros(1,N-2);
By = zeros(1,N-2);
Bz = zeros(1,N-2);

% calc the binormal vector
for i=1:N-1
    Bx(i) = tany(i)*normz(i) - normy(i)*tanz(i);
    By(i) = tanz(i)*normx(i) - normz(i)*tanx(i);
    Bz(i) = tanx(i)*normy(i) - normx(i)*tany(i);
end

hold on;
quiver3(ft(1:N-1),gt(1:N-1),ht(1:N-1),tanx(1:N-1),tany(1:N-1),tanz(1:N-1),0.5);
quiver3(ft(1:N-2),gt(1:N-2),ht(1:N-2),normx(1:N-2),normy(1:N-2),normz(1:N-2),0.5);
quiver3(ft(1:N-2),gt(1:N-2),ht(1:N-2),Bx(1:N-2),By(1:N-2),Bz(1:N-2),0.5);
hold off;


%% Problem 4
curvature=sqrt(dtanx.^2+dtany.^2+dtanz.^2)./sqrt(df(1:1999).^2+dg(1:1999).^2+dh(1:1999).^2);
maxval = max(curvature);
index = find(curvature == maxval);



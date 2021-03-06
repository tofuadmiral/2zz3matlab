clc; close all; clear all;
%% Question 1
% plot the fourier transforms in the negative direction for all of them,
% and also make sure to do it for even, odd and identity shift functions. 
func1 = zeros(1,7); 
x = linspace(0,2,201); 

f = sqrt(8+x.^3); 
x2 = linspace(0,-2,201);
f2 = sqrt(8+x2.^3); 
a0 = 2/2*trapz(x,f);

% they asked for 7 curves and my p was 2
for j = 1:7 
    y = f.*cos(2*j*pi*x/2);
    func1(j) = 2/2*trapz(x,y); 
    y1 = f.*sin (2*j*pi*x/2);
    b(j) = 2/2*trapz(x,y1);
end

fh = zeros(7,801);
x3 = linspace(-4,4,801); 

% now we actually need to get the fouriers
for i = 1:7 
for j = 1:i 
    fh(i,:) = fh(i,:) + func1(j)*cos(2*j*pi*x3/2); %p = 3 so x3/3 
end
end

fh = a0/2 + fh; 
 
plot(x,f,'r','LineWidth',3);
hold on;
plot(x+2,f,'b','LineWidth',3);
hold on;
plot(x3,fh,'k','LineWidth',1);


xlim = ([-4 4]);
xlabel('X');
ylabel('Y');
title('Ahmed Fuad Ali, 400075937');
legend('f','fh','f2','f2h','Location','best');



%% Question 2
theerror = fh(7,501);
theoreticalanswer = sqrt(9); 
A = abs(theoreticalanswer-theerror);


%% Question 3
% plot the construction of a discrete fourier transform of the given
% function 

% first we want to set up our environment and initialize some variables
clc; close all; clear all;
x = linspace(0,2*pi,200);
N = length(x);

% this was my function
f = 9*sin(15*pi*x) + 2*cos(8*pi*x);

% now we need to use the stuff from the tutorial 
F = fft(f); 

abs_F = abs(fftshift(F));
alpha = (0:N-1)-(N/2);
plot(alpha,abs_F);

% now plot it
semilogy(alpha,abs_F);
title('Ahmed Fuad Ali, 400075937');
xlabel('x');
ylabel('y');




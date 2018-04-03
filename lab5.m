%% Question 1
%Format will be f(x) = ..., 0 < x < p (my p was 3)
%truncated after 1,2..n terms, for me n = 7

a = zeros(1,7); %7 terms so from 1 to 7 for me

x = linspace(0,2,201); %200 elements for me and from 0 to 3
%Function
f = sqrt(8+x.^3); %The one which is given, does NOT change

%Extension
x2 = linspace(0,-2,201); %This is the extension
f2 = sqrt(8-x2.^3); %ONLY THIS CHANGES FOR OTHER GRAPHS


a0 = 2/2*trapz(x,f); %p = 3 so 2/p

for n = 1:7 %7 curves for me
    y = f.*cos(n*pi*x/2); %my p was 3 so its divided by 3 
    a(n) = 2/2*trapz(x,y); %p = 3
    %Bottom for part C    
end

%Estimate values of fourier series of f

fh = zeros(7,801); %800 elements for fh and f2h
x3 = linspace(-4,4,801); % 0 to 3 is 200, 3 to 6 is 400, so 800
%change -6,6 according to your interval [it's given]


for m = 1:7 %7 curves
for n = 1:m %dont touch this
    fh(m,:) = fh(m,:) + a(n)*cos(n*pi*x3/2); %p = 3 so x3/3 
end
end

fh = a0/2 + fh; 
 
plot(x,f,'r','LineWidth',3);
hold on;
plot(x2,f2,'b','LineWidth',3);
hold on;
plot(x3,fh,'k','LineWidth',1);


xlim = ([-4 4]);  %Change interval according to yours
xlabel('X');
ylabel('Y');
title('Ahmed Fuad Ali, 400075937');
legend('f','fh','f2','f2h','Location','best');

%----------- end
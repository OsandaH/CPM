% 3 point and 5 point formula for numerical differentiation

clear;
clc;

x = 1;
h = 0.2;
exact = cos(x);
for i =1:5;
    h = h/2;

    fp3 = (sin(x+h) - sin(x-h))/(2*h);
    diff1(i) = exact-fp3;
    fp5 = (sin(x-2*h) - 8*sin(x-h) + 8*sin(x+h) - sin(x+2*h))/(12*h);
    diff2(i) = exact - fp5;
end
disp(i);
disp(diff1);
disp(diff2);

%% Numerical Integration 
% Simpsons Rule

clear;
clc

n = 2;
exact = exp(1)-1;

for k = 1:5;
    sum = exp(0);
    n = n*2;
    h(k)= 1/n;
    fac = 2;
    for i = 1:n-1;
        if(fac == 2);
            fac = 4;
        else;
            fac = 2;
        end;
    
        x = i*h(k);
        sum = sum + exp(x)*fac;
    end
    sum = sum + exp(1);
    result = h(k)*sum/3;
    diff3(k) = exact - result;
end;
disp(h)
disp(diff3);

%% Root Finding 

%f(x) = x^2 -5 step size = 0.5 tolerance = 10^-6 initial guess x=1
%Brute Force Method

clear
clc
tolx = 1e-6;
x = 1;
dx = 0.5;

iter = 0;
fold = x*x - 5;
while abs(dx)>tolx
    iter = iter+1;
    x = x+dx;
    fnew = x*x-5;
    roots(iter) = x;
        if(fold * fnew <0)
            x = x-dx;
            dx = dx/2;
        end
end
exact = sqrt(5)
disp(roots);

%%
% Secant Method
%x0=1 x1=1.5

clc
clear

tolx = 1e-6;
x1 = 1;
x2 = 1.5;
dx = x2-x1;
iter = 0;
while abs(dx) > tolx
    iter = iter + 1;
    xsave = x2;
    xf1 = x1*x1-5;
    xf2 = x2*x2-5;
    x2 = x2 - ((xf2*dx)/(xf2-xf1));
    x1= xsave;
    dx = x2-x1;
    roots(iter) = x2;
end
disp(roots);
%%
% Newton Raphson Method

clc
clear

tolx = 1e-6;
x = 1;
dx = 0.5;
iter = 0;
while abs(dx) > tolx
    iter = iter + 1;
    xold = x;
    x = x - ((x*x-5)/(2*x));
    dx = xold - x;
    roots(iter) = x;
end
disp(roots);

%% Initial Value Problem
% Eulers Method

clc;
clear;

h = 0.002;
n = 1/h;
y = 1;
for i = 0:n-1
    x = i*h;
    y = y+h*(-x*y);
    x = x+h;
    xx(i+1) = x;
    y1(i+1) = exp(-0.5*x*x);
    y2(i+1) = y;
end
plot(xx,y1);hold on;
plot(xx,y2,'.r');
hold off;

%% Runge Kutta 2nd Order

clear;
clc
h = 0.02;
n = 1/h;
y = 1;

for i=0:n-1
    x = i*h;
    k1 = h*(-x*y);
    y = y+h*(-(x-(h/2))*(y+(k1/2)));
    %x=x+h; 
    % xx(i+1)=x; 
    % y1(i+1)=exp(-0.5*x*x); 
    % y2(i+1)=y;
end
for j=n:-1:1
    x=j*h; k1=-h*(-x*y); y=y-h*(-(x-(h/2))*(y+(k1/2)));
    x=x-h; xx(j)=x; y1(j)=exp(-0.5*x*x); y2(j)=y;
end

plot(xx,y1);hold on;
plot(xx,y2,'.r');
hold off;

%% In built methods

tspan = 0:0.1:5;
x0 = 0;
odefun = @(t, x) 3 * exp(-x);
[t, x] = ode45(odefun, tspan, x0);

y = -3*exp(-tspan)+3; %analytical 

hold on;
plot(t,x,'.r');
plot(t,y);

%%
clc
clear;
x = 0:0.1:3;
y0 = 1;

myfun = @(x,y) -x*y;

[x,y] = ode45(myfun,x,y0);

y2 = exp(-0.5*x.^2);
plot(x,y,'.r'); hold on;
plot(x,y2,'b');hold off;

%% 2nd order ODE  

function dydt = myfun3(t,y)
    dy1dt = y(2);
    dy2dt = -5*y(2) + 4*y(1) + sin(10*t);
    dydt = [dy1dt; dy2dt];
end


t = 0:0.01:3;
y0 = 0;
z0 = 0;

[t,y] = ode45(@myfun3,t,[y0 z0]);

plot(t,y(:,1));

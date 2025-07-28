clc;
clear all;

f = @(x)  x.^4 - 7.223*x.^3 + 13.447*x.^2 + 0.672*x - 10.223;
df = @(x) 4*x.^3 - 3*7.223*x.^2 + 2*13.447*x + 0.672;

tolerance = 1e-6;
iter = 100;

root = newton_method(f, df, 4, tolerance, iter);
fprintf('Positive Root ≈ %.6f\n', root);

function root = newton_method(f,df,x0,tolerance,iter)
    x = x0;
    for i = 1:iter
        xnew = x - f(x)/df(x);
        if abs(xnew - x)<tolerance
            root = xnew;
            return;
        end
        x = xnew;
    end
    error('Newton-Raphson did not converge.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;

% Define the function
f = @(x) x.^4 - 7.223*x.^3 + 13.447*x.^2 + 0.672*x - 10.223;

% Initial guesses (choose near the expected positive root)
x0 = 4;    % first guess
x1 = 8;    % second guess

tol = 1e-6;
max_iter = 100;

for i = 1:max_iter
    f0 = f(x0);
    f1 = f(x1);
    
    % Secant formula
    x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
    
    % Check convergence
    if abs(x2 - x1) < tol
        fprintf('Positive Root ≈ %.6f found after %d iterations\n', x2, i);
        break;
    end
    
    % Update guesses
    x0 = x1;
    x1 = x2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;

% Parameters
a = 1; b = 8; N = 100;  % number of internal points
h = (b-a)/(N+1);
x = a:h:b;

% Initialize matrix and RHS
A = zeros(N,N);
rhs = zeros(N,1);

for i = 1:N
    xi = x(i+1); % skip x0=1
    alpha = 1/(h^2) - 1/(2*h*xi);
    beta  = -2/(h^2) + 1;
    gamma = 1/(h^2) + 1/(2*h*xi);

    if i > 1
        A(i,i-1) = alpha;
    else
        rhs(i) = rhs(i) - alpha*1;  % y(1) = 1
    end

    A(i,i) = beta;

    if i < N
        A(i,i+1) = gamma;
    else
        rhs(i) = rhs(i) - gamma*0;  % y(8) = 0
    end
end

% Solve
Y = A\rhs;

% Full solution including boundaries
y = [1; Y; 0];

% Display results
disp([x' y])
plot(x,y,'-o'), xlabel('x'), ylabel('y'), title('BVP (i) Solution - Matrix Method')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;

% Shooting parameters
xspan = [0 10];
target = 2; % y'(∞) ≈ y'(10) = 2

% Guess initial y3(0)
y3_guess = 1;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);

% Shooting iteration
for iter = 1:20
    [x, Y] = ode45(@(x,Y) odefun(x,Y), xspan, [0; 0; y3_guess], options);
    yprime_end = Y(end,2);

    if abs(yprime_end - target) < 1e-4
        fprintf('Converged: y3(0) = %.6f\n', y3_guess);
        break;
    else
        % Adjust guess (simple secant approach)
        y3_guess = y3_guess * target / yprime_end;
    end
end

plot(x, Y(:,1)), xlabel('x'), ylabel('y'), title('BVP (ii) Shooting Method')

%% ODE function
function dY = odefun(x,Y)
    y1 = Y(1); y2 = Y(2); y3 = Y(3);
    dY = [y2; y3; -y1*y3];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;

% Time span
t0 = 0; tf = 5; h = 0.1;
t = t0:h:tf;
n = length(t);

% Initial Conditions
y1_euler = zeros(1,n); y2_euler = zeros(1,n);
y1_rk2 = zeros(1,n); y2_rk2 = zeros(1,n);

% Set initial values
y1_euler(1) = 2; y2_euler(1) = 0;
y1_rk2(1) = 2; y2_rk2(1) = 0;

% Euler Method
for i = 1:n-1
    dy1 = y2_euler(i);
    dy2 = -4*y2_euler(i) - 2*y1_euler(i);

    y1_euler(i+1) = y1_euler(i) + h*dy1;
    y2_euler(i+1) = y2_euler(i) + h*dy2;
end

% RK2 Method (Heun’s method)
for i = 1:n-1
    % k1
    k1_y1 = y2_rk2(i);
    k1_y2 = -4*y2_rk2(i) - 2*y1_rk2(i);

    % Predictor
    y1_pred = y1_rk2(i) + h*k1_y1;
    y2_pred = y2_rk2(i) + h*k1_y2;

    % k2
    k2_y1 = y2_pred;
    k2_y2 = -4*y2_pred - 2*y1_pred;

    % Update
    y1_rk2(i+1) = y1_rk2(i) + h*(k1_y1 + k2_y1)/2;
    y2_rk2(i+1) = y2_rk2(i) + h*(k1_y2 + k2_y2)/2;
end

% Analytical solution (for accuracy estimation)
% Characteristic eq: r^2 + 4r + 2 = 0 --> r = -2 ± sqrt(2)
r1 = -2 + sqrt(2); r2 = -2 - sqrt(2);
A = (r1*y1_euler(1) - y2_euler(1)) / (r1 - r2);
B = y1_euler(1) - A;
y_true = A*exp(r1*t) + B*exp(r2*t);

% Error estimation
error_euler = abs(y1_euler - y_true);
error_rk2 = abs(y1_rk2 - y_true);

% Plot
figure;
plot(t, y_true, 'k-', t, y1_euler, 'r--', t, y1_rk2, 'b-.');
legend('Analytical', 'Euler', 'RK2');
xlabel('t'); ylabel('y(t)');
title('Comparison of Euler and RK2 Methods');

figure;
plot(t, error_euler, 'r--', t, error_rk2, 'b-.');
legend('Euler Error', 'RK2 Error');
xlabel('t'); ylabel('Absolute Error');
title('Error Comparison');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

f = @(x, y) 2 .* x .* y + 2 .* x - x.^2 - 2 .* y.^2 + 72;

%setting up the domain
x0 = 0; x1 = 8;
y0 = 0; y1 = 6;
%number of subintervals
nx = 8;
ny = 8;
%size between points
dx = (x1 - x0)/nx;
dy = (y1 - y0)/ny;
%creating the subintervals
domainx = x0:dx:x1;
domainy = y0:dy:y1;
%creating coeffients matrix
coeff_x = 2*ones(1,nx+1);
coeff_y = 2*ones(1,ny+1);
%making odd term as 4
coeff_x(2:2:nx) = 4;
coeff_y(2:2:ny) = 4;
%making first term and last term as 1
coeff_x([1 nx+1])=1;
coeff_y([1 ny+1])=1;

[X,Y] = meshgrid(domainx,domainy);
F=f(X,Y);

sn = (dx * dy / 9) * (coeff_y * F * coeff_x');

Area = (x1 - x0) * (y1 - y0);
T_avg=sn/Area





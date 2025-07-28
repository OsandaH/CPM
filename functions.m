%%%%%% Finding Roots %%%%%%
%% secant method

function root = secant_method(f, x0, x1, tol, max_iter)
    for i = 1:max_iter
        fx0 = f(x0);
        fx1 = f(x1);
        if fx1 - fx0 == 0
            error('Zero denominator in Secant method');
        end
        x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
        if abs(x2 - x1) < tol
            root = x2;
            return
        end
        x0 = x1;
        x1 = x2;
    end
    root = x1;
end

%% newton raphson method

function root = newton_raphson(f, df, x0, tol, max_iter)
    for i = 1:max_iter
        x1 = x0 - f(x0)/df(x0);
        if abs(x1 - x0) < tol
            root = x1;
            return
        end
        x0 = x1;
    end
    root = x0; % return best approximation if convergence not reached
end


% Tolerance and iteration limits
tol = 1e-6;
max_iter = 100;

% Function (i)
f1 = @(x) x.^4 - 8.6*x.^3 - 35.51*x.^2 + 464.4*x - 998.46;
df1 = @(x) 4*x.^3 - 3*8.6*x.^2 - 2*35.51*x + 464.4;

% Newton-Raphson
root1_NR = newton_raphson(f1, df1, 6, tol, max_iter);
% Secant
root1_Sec = secant_method(f1, 5, 7, tol, max_iter);

fprintf('Root for function (i): Newton-Raphson = %.6f, Secant = %.6f\n', root1_NR, root1_Sec);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% approximate the integral of f from a to b %%% boundary value Problem 

%%%%% simpsons rule %%%%%

function I = simpsons_rule(f, a, b, n)
    % Simpson's Rule to approximate the integral of f from a to b
    % n must be an even number

    if mod(n, 2) ~= 0
        error('Number of intervals (n) must be even');
    end

    h = (b - a) / n;         % Step size
    x = a:h:b;               % Divide interval [a, b] into n parts
    y = f(x);                % Evaluate function at all x points

    % Apply Simpson's Rule formula
    I = h/3 * (y(1) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2)) + y(end));
end

% (i) ∫₀^{π/2} [sin(2cos(x)) * sin^2(x)] dx
f1 = @(x) sin(2*cos(x)) .* (sin(x)).^2;                 % Function to integrate
a = 0;                                                  % Start of interval
b = pi;                                                 % End of interval
n = 100;                                                % Number of intervals (even number)

result = simpsons_simple(f, a, b, n);
fprintf('Result = %.6f\n', result);


% (iv) ∫_{x=2}^{3} ∫_{y=x}^{2x^3} (x^2 + y) dy dx
% Outer integral over x, inner integral over y

f4 = @(x) arrayfun(@(xval) ...
    simpsons_rule(@(y) xval^2 + y, xval, 2*xval^3, n), x);  % Inner integral as function of x

I4 = simpsons_rule(f4, 2, 3, n);
fprintf('Integral (iv): %.6f\n', I4);


%%%%%%% Simpsons 2D %%%%%%%%%%

function I = simpsons_2d(f, ax, bx, ay, by, nx, ny)
    % Simpson's Rule for double integrals over [ax,bx] x [ay,by]
    % nx and ny must be even

    if mod(nx,2) ~= 0 || mod(ny,2) ~= 0
        error('nx and ny must be even');
    end

    hx = (bx - ax) / nx;
    hy = (by - ay) / ny;

    x = ax:hx:bx;
    y = ay:hy:by;
    
    % Create 2D grid
    [X, Y] = meshgrid(x, y);
    Z = f(X, Y);

    % Simpson's weights
    wx = ones(1, nx+1);      % nx+1 points
    wx(2:2:end-1) = 4;
    wx(3:2:end-2) = 2;

    wy = ones(1, ny+1);      % ny+1 points
    wy(2:2:end-1) = 4;
    wy(3:2:end-2) = 2;

    % Apply weights to the 2D grid
    W = wy' * wx;

    % Simpson's 2D formula
    I = hx * hy / 9 * sum(sum(W .* Z));
end


% Define the function to integrate
f = @(x, y) x.^2 + y.^2;

% Set integration bounds and subdivisions
ax = 0; bx = 2;     % x from 0 to 2
ay = 0; by = 2;     % y from 0 to 2
nx = 10; ny = 10;   % Must be even

% Call the function
I = simpsons_2d(f, ax, bx, ay, by, nx, ny);

% Display the result
fprintf('Approximate integral value = %.6f\n', I);


%%%%%% Initial Value Problem %%%%%%

%% Common Derivative %%
function dy = ode_system(t, y)
    dy = zeros(2,1);
    dy(1) = y(2);                     % dy1/dt = y2
    dy(2) = -2*y(2) - 4*y(1);     % dy2/dt = -2y2 - 4y1 
end

%% Euler Method %%
function [t, Y] = euler_method(f, y0, t0, tf, h)
    t = t0:h:tf;
    n = length(t);
    Y = zeros(n, length(y0));
    Y(1,:) = y0;
    
    for i = 1:n-1
        Y(i+1,:) = Y(i,:) + h * f(t(i), Y(i,:)')';
    end
end

%% RK 2 %%

function [t, Y] = rk2_method(f, y0, t0, tf, h)
    t = t0:h:tf;
    n = length(t);
    Y = zeros(n, length(y0));
    Y(1,:) = y0;
    
    for i = 1:n-1
        k1 = f(t(i), Y(i,:)')';
        k2 = f(t(i)+h, Y(i,:) + h*k1)';
        Y(i+1,:) = Y(i,:) + h/2 * (k1 + k2);
    end
end


%% RK 4 %%

function [t, Y] = rk4_method(f, y0, t0, tf, h)
    t = t0:h:tf;
    n = length(t);
    Y = zeros(n, length(y0));
    Y(1,:) = y0;
    
    for i = 1:n-1
        k1 = f(t(i), Y(i,:)')';
        k2 = f(t(i)+h/2, Y(i,:) + h/2*k1)';
        k3 = f(t(i)+h/2, Y(i,:) + h/2*k2)';
        k4 = f(t(i)+h, Y(i,:) + h*k3)';
        Y(i+1,:) = Y(i,:) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
    end
end

% Parameters
t0 = 0; tf = 5; h = 0.1;
y0 = [2; 0]; % y(0) = 2, y'(0) = 0

% Solve using different methods
[t_e, Y_e] = euler_method(@ode_system, y0, t0, tf, h);
[t_rk2, Y_rk2] = rk2_method(@ode_system, y0, t0, tf, h);
[t_rk4, Y_rk4] = rk4_method(@ode_system, y0, t0, tf, h);

% Plot
figure;
plot(t_e, Y_e(:,1), 'r--', t_rk2, Y_rk2(:,1), 'b-.', t_rk4, Y_rk4(:,1), 'k', 'LineWidth', 2);
legend('Euler', 'RK2', 'RK4 (Reference)');
xlabel('t'); ylabel('y(t)');
title('Comparison of Euler, RK2, and RK4 Methods');
grid on;

%%%%%%%%%%%%%%% Boundary Vlaue Problem %%%%%%%%%%%%%%%%%%%
%% (d3y/dx3) + y*d(2y/dx2) = 0 y(0) = 0 y'(0) = 0 y(inf) = 2

function bvp_shooting_method()
    % Target boundary
    x0 = 0; xf = 10; h = 0.01;
    x = x0:h:xf;
    
    % Initial conditions: y1 = y, y2 = y', y3 = y''
    y1_0 = 0; y2_0 = 0;

    % Try initial guesses for y3(0)
    g1 = 0; g2 = 5;  % two guesses for y3(0)
    
    % Solve for both guesses
    y1g1 = run_rk3(x, [y1_0; y2_0; g1]);
    y1g2 = run_rk3(x, [y1_0; y2_0; g2]);

    % Secant method to refine guess
    tol = 1e-4;
    while abs(y1g2(end,2) - 2) > tol
        g = g2 - (y1g2(end,2) - 2) * (g2 - g1) / (y1g2(end,2) - y1g1(end,2));
        g1 = g2; g2 = g;
        y1g1 = y1g2;
        y1g2 = run_rk3(x, [y1_0; y2_0; g2]);
    end
    
    % Final solution
    y_final = y1g2(:,1);
    
    % Plot
    plot(x, y_final, 'r', 'LineWidth', 2);
    xlabel('x'); ylabel('y(x)');
    title('Shooting Method Solution');
    grid on;
end

% Runge-Kutta 4th order solver for 3-variable system
function Y = run_rk3(x, y0)
    n = length(x);
    Y = zeros(n, 3);
    Y(1,:) = y0';

    h = x(2) - x(1);

    for i = 1:n-1
        k1 = ode3(x(i), Y(i,:)') * h;
        k2 = ode3(x(i)+h/2, Y(i,:)' + k1/2) * h;
        k3 = ode3(x(i)+h/2, Y(i,:)' + k2/2) * h;
        k4 = ode3(x(i)+h, Y(i,:)' + k3) * h;
        Y(i+1,:) = Y(i,:) + (k1' + 2*k2' + 2*k3' + k4') / 6;
    end
end

function dy = ode3(x, y)
    dy = zeros(3,1);
    dy(1) = y(2);
    dy(2) = y(3);
    dy(3) = -y(1)*y(3);
end


%%%% Interpolation %%%
% Define g(x)
g = @(x) sin(4*pi*x);

% Interval
a = 0; b = 1;

% Define 11 equally spaced sample points (xk)
k = 0:10;
x = a + k * (b - a) / 10;
y = g(x);  % Function values at x

% Interpolation points
xi = linspace(a, b, 1000);  % Fine grid for plotting and error

% Interpolation using 'cubic'
yi_cubic = interp1(x, y, xi, 'cubic');
err_cubic = norm(abs(yi_cubic - g(xi)), 'inf');

% Interpolation using 'spline'
yi_spline = interp1(x, y, xi, 'spline');
err_spline = norm(abs(yi_spline - g(xi)), 'inf');

% Display errors
fprintf('Max error with cubic:  %.6e\n', err_cubic);
fprintf('Max error with spline: %.6e\n', err_spline);

% Plot results
figure;
plot(xi, g(xi), 'k-', 'LineWidth', 1.5); hold on;
plot(xi, yi_cubic, 'r--', 'LineWidth', 1.5);
plot(xi, yi_spline, 'b-.', 'LineWidth', 1.5);
legend('Original g(x)', 'Cubic Interpolation', 'Spline Interpolation');
xlabel('x'); ylabel('y');
title('Interpolation of g(x) = sin(4πx)');
grid on;



%%%%%%% 
% Given data
years = 1950:10:1990;           % 5 columns (time)
service = 10:10:30;             % 3 rows (years of service)

wage = [150.697 199.592 187.625 179.323 195.072;   % service = 10
        250.287 203.212 179.092 322.767 226.505;   % service = 20
        153.706 426.730 249.633 120.281 598.243];  % service = 30

% Interpolation point
year_interp = 1975;
service_interp = 15;

% Interpolate using interp2
wage_interp = interp2(years, service, wage, year_interp, service_interp, 'spline');

% Display result
fprintf('Estimated wage in 1975 for 15 years of service: %.3f\n', wage_interp);


%%%%%%%%%%%%%%% Boundary Vlaue Problem %%%%%%%%%%%%%%%%%%%

function bvp_shooting_method2()
    % Shooting method for: y'' + 2x y' = 0, y(0) = 1, y(xf) = 0
    % We'll solve for xf = 1.5 and xf = 3
    fprintf('Solving for x = 1.5...\n');
    solve_bvp(1.5);

    fprintf('\nSolving for x = 3...\n');
    solve_bvp(3);
end

function solve_bvp(xf)
    x0 = 0; h = 0.01;
    x = x0:h:xf;

    y1_0 = 1;  % y(0)
    g1 = -1; g2 = -5; % Initial guesses for y'(0)

    y1g1 = run_rk4(x, [y1_0; g1]);
    y1g2 = run_rk4(x, [y1_0; g2]);

    tol = 1e-5;
    while abs(y1g2(end,1)) > tol
        g = g2 - y1g2(end,1)*(g2 - g1)/(y1g2(end,1) - y1g1(end,1));
        g1 = g2;
        g2 = g;
        y1g1 = y1g2;
        y1g2 = run_rk4(x, [y1_0; g2]);
    end

    % Final solution
    y_final = y1g2(:,1);

    plot(x, y_final, 'DisplayName', sprintf('x_f = %.1f', xf)); hold on;
    xlabel('x'); ylabel('y(x)');
    title('Shooting Method Solution for BVP');
    legend show;
    grid on;
end

function Y = run_rk4(x, y0)
    n = length(x);
    Y = zeros(n, 2);
    Y(1,:) = y0';
    h = x(2) - x(1);

    for i = 1:n-1
        k1 = ode_rhs(x(i), Y(i,:)') * h;
        k2 = ode_rhs(x(i)+h/2, Y(i,:)' + k1/2) * h;
        k3 = ode_rhs(x(i)+h/2, Y(i,:)' + k2/2) * h;
        k4 = ode_rhs(x(i)+h, Y(i,:)' + k3) * h;
        Y(i+1,:) = Y(i,:) + (k1' + 2*k2' + 2*k3' + k4')/6;
    end
end

function dy = ode_rhs(x, y)
    % y = [y1; y2], dy1/dx = y2, dy2/dx = -2*x*y2
    dy = zeros(2,1);
    dy(1) = y(2);
    dy(2) = -2*x*y(2);
end


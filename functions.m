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


%% approximate the integral of f from a to b

%%% simpsons rule

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


%%%%%%% Simposns 2D %%%%%%%%%%

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

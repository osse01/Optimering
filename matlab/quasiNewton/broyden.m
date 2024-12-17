%% Broyden-family Quasi Newton method
clc
clear
close all

maxiter = 200;
epsilon = 1e-6; % Tolerance
xi = 1; % 0 = DFP and 1 = BFGS

f = @(x) (x(1) - 2)^4 + (x(1) - 2*x(2))^2;

g = @(x) [4*( x(1)-2)^3 + 2*( x(1) - 2*x(2))
           -4*( x(1) - 2*x(2) )];
g2 = @(x) [12*(x(1) - 2)^2 + 2,  -4
                -4                8];

xk = [0 3]'; % Initial guess
Dk = [1 0
      0 1]; % Initial hessian inverse approximation

% For plotting
fvals = [f(xk)];
xvals = [xk];

for i = 1:maxiter
    % Variables
    alpha_k = armijo(xk, f, g, g2);
    xnew = xk - alpha_k*Dk*g(xk);
    gnew = g(xnew);
    
    pk = xnew - xk;
    qk = gnew - g(xk);
    tauk = qk'*Dk*qk;
    vk = pk/(pk'*qk) - Dk*qk/tauk;
    BFGS = tauk*(vk*vk');
    Dk = Dk + pk*pk' / (pk'*qk) - Dk*qk*qk'*Dk / (qk'*Dk*qk) + xi*BFGS;
    
    xk = xnew;
    if ( f(xk) < epsilon)
        break;
    end
    fvals = [fvals f(xk)];
    xvals = [xvals xk];
end

disp("Optimal function value: " + fvals(end))
disp("Optimal x value: ")
disp(xvals(:,end))

% Plot the objective function values
figure;
semilogy(fvals, '-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Iteration');
ylabel('Objective Function Value');
title('Convergence of Objective Function');
grid on;

% Plot the trajectory in R^2
figure;
plot(xvals(1, :), xvals(2, :), '-o', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
scatter(2, 1, 100, 'r', 'filled', 'DisplayName', 'Optimal Point');
xlabel('x_1');
ylabel('x_2');
title('Optimization Path in R^2');
legend('Path', 'Optimal Point');
grid on;
axis equal;
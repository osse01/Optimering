% Setup
clear 
close all
clc

% Variables
x_0 = [0, 3]';
x_star = [2, 1]';
epsilon = 1e-9;
x_curr = x_0;
x_prev = x_0;

% Define functions
% Function to minimize.
f = @(x) (x(1) - 2)^4 + (x(1) - 2*x(2))^2;
% Gradient of function to minimize.
g = @(x) [(4*(x(1)-2)^3 + 2*(x(1) - 2*x(2) )), (-4*(x(1) - 2*x(2)))]';

% Conjugate coefficient
b = @(x_curr, x_prev) ( g(x_curr)'*(g(x_curr) - g(x_prev)) ) / ( norm(g(x_prev))^2 );

d_prev = -g(x_0);

% Storage for plotting
objective_values = [];
trajectory = [];  % Store the trajectory in R^2
iterations = 0;

while norm(g(x_curr)) >= epsilon
    % Store values for plotting
    objective_values = [objective_values, f(x_curr)];
    trajectory = [trajectory, x_curr];
    iterations = iterations + 1;
    
    % Define phi(alpha) based on the current search direction
    phi = @(alpha) f(x_curr + alpha * d_prev);
    
    % Use fminbnd to find the optimal step length alpha
    alpha_opt = fminbnd(phi, 0, 1e6);
    
    % Update point
    x_prev = x_curr;
    x_curr = x_curr + alpha_opt * d_prev;
    
    % Update search direction
    d_curr = -g(x_curr) + b(x_curr, x_prev) * d_prev;
    d_prev = d_curr; 
end
% Store the final values
objective_values = [objective_values, f(x_curr)];
trajectory = [trajectory, x_curr];

% Results
disp("Summary:")
disp("Final function value:")
disp(f(x_curr))
disp("Final point:")
disp(x_curr)

% Plot the objective function values
figure;
semilogy(0:iterations, objective_values, '-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Iteration');
ylabel('Objective Function Value');
title('Convergence of Objective Function');
grid on;

% Plot the trajectory in R^2
figure;
plot(trajectory(1, :), trajectory(2, :), '-o', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
scatter(x_star(1), x_star(2), 100, 'r', 'filled', 'DisplayName', 'Optimal Point');
xlabel('x_1');
ylabel('x_2');
title('Optimization Path in R^2');
legend('Path', 'Optimal Point');
grid on;
axis equal;

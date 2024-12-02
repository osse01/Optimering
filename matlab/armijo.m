% Setup
clear
close all
clc

% Variables
x_0 = [0, 1000]';
x_curr = x_0;
x_star = [2, 1]';
epsilon = 1e-3;
sigma = 0.5; % Must be in [0,1]
run = true;
vals = zeros(1e4, 1);

% Define functions
% Function.
f = @(x) (x(1) - 2)^4 + (x(1) - 2*x(2))^2;
% Gradient of function.
g = @(x) [4*(x(1)-2)^3 + 2*(x(1) - 2*x(2)); -4*(x(1) - 2*x(2))];

vals(1, :) = f(x_curr);

% Initial Newton Step.
x_curr = x_curr - f(x_curr)/norm(g(x_curr));

i = 2;
while (run)
    vals(i, :) = f(x_curr);
    d = -g(x_curr);  % Descent direction
    t = 1e0;         % Initial step size

    while (f(x_curr + t * d) >= f(x_curr) + sigma * g(x_curr)' * d * t)
        t = t * 0.5;  
    end
    x_curr = x_curr + t * d;  % Update
    if (norm(g(x_curr)) <= epsilon)
        run = false;
    end
    i = i + 1;
end
values = vals(1:ceil(length(find(vals))/2), :);
iterations = linspace(1, ceil(length(find(vals))/2), ceil(length(find(vals))/2))';

disp('Optimal x:');
disp(x_curr);

% Plot values against iterations.
disp("Plot!")
% figure
plot(iterations, values, '-o')
xlabel('Iterations')
ylabel('Values')
title('Values vs Iterations')
grid on

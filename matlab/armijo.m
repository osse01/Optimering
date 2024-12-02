% Setup
clear
close all
clc

% Variables
x_curr = [0, 3]';       % Initial guess.
epsilon = 1e-3;         % Error tolerance.
sigma = 0.5;            % Must be in [0,1].
vals = zeros(1e4, 1);   % To Store Values for Plot Later.
iter = 1;               % Iteration Counter.

% Define functions.
f = @(x) (x(1) - 2)^4 + (x(1) - 2*x(2))^2;
g = @(x) [4*(x(1)-2)^3 + 2*(x(1) - 2*x(2)); -4*(x(1) - 2*x(2))];
h = @(x) [12*(x(1)-2)^2 + 2, -4;-4, 8];

while (true)
    vals(iter, :) = f(x_curr);
    % Initial Newton Step.
    t = (norm(g(x_curr))^2)/( g(x_curr)'*(h(x_curr)*g(x_curr) ));

    % Descent Direction.
    d = -g(x_curr);

    while (f(x_curr + t * d) >= f(x_curr) + sigma * g(x_curr)' * d * t)
        t = t * 0.5;  
    end
    % Update x.
    x_curr = x_curr + t * d;  
    if (norm(g(x_curr)) <= epsilon)
        break;
    end
    iter = iter + 1;
end
values = vals(1:ceil(length(find(vals))/2), :);
iterations = linspace(1, ceil(length(find(vals))/2), ceil(length(find(vals))/2))';

disp('Optimal x:');
disp(x_curr);

% Plot values against iterations.
disp("Plot!")
plot(iterations, values, '-o')
xlabel('Iterations')
ylabel('Values')
title('Values vs Iterations')
grid on

% Setup
clear 
close all
clc

% Variables
x_0 = [0, 3]';
x_star = [2, 1]';
epsilon = 1e-16;
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

while norm(g(x_curr)) >= epsilon
    % Define phi(alpha) based on the current search direction
    phi = @(alpha) f(x_curr + alpha * d_prev);
    
    % Use fminbnd to find the optimal step length alpha
    % From matlab page,
    % fminbnd is a function file. The algorithm is based on golden section search
    % and parabolic interpolation. Unless the left endpoint x1 is very close to 
    % the right endpoint x2, fminbnd never evaluates fun at the endpoints, so fun
    % need only be defined for x in the interval x1 < x < x2.
    % https://se.mathworks.com/help/matlab/ref/fminbnd.html
    alpha_opt = fminbnd(phi, 0, 1e6)
    
    % Update point
    x_prev = x_curr;
    x_curr = x_curr + alpha_opt * d_prev;
    
    % Update search direction
    d_curr = -g(x_curr) + b(x_curr, x_prev) * d_prev;
    d_prev = d_curr; 
end

% Results
disp("Summary:")
disp("Final function value:")
disp(f(x_curr))
disp("Final point:")
disp(x_curr)

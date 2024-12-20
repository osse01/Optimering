%% Gauss Newton method for data regression
% Also known as Levenberg-Marquardt
clc
clear
close all

% Data
t = [1,2,4,5,8];
y = [3,4,6,11,20];

g = @(x) x(1)*exp(x(2).*t) - y ;
f = @(x) sum(g(x).^2);
grad =  @(x) [exp(x(2).*t); x(1).*t.*exp(x(2).*t) ];

% Initial values
x_k = [2.5, 0.25]'; % Start point
epsilon = 1e-12;
maxIter = 100;
alpha = 1e0;
C = 100;
delta = C*eye(2);

fvals = [f(x_k)];
xvals = [x_k];
for iter = 1:maxIter
    disp("Iter: " + iter + "-------------------")
    g_k = g(x_k)';
    gradg = grad(x_k);
    G = gradg*gradg';
    deter = det(G)

    % When G near singular => delta more significance
    if abs(deter) < 1e-6 || isnan(deter)
        delta = 1e6*eye(2)
    else
        delta = 1/deter * eye(2) 
    end

    % Update x_k
    alpha = 0.999*alpha; % Damping factor
    x_k = x_k - alpha*(G+delta)\(gradg*g_k)

    % Stopping criteria
    gradf = grad(x_k)*g(x_k)';
    disp("|| grad(f)||: " + norm(gradf));
    if norm(gradf) < epsilon
        break;
    end

    fvals = [fvals, f(x_k)];
    xvals = [xvals, x_k];
end

tplot = 0:0.01:9;
func = @(t) x_k(1)*exp(x_k(2).*t);

figure;
plot(tplot,func(tplot),'LineWidth', 2, 'MarkerSize', 6);
hold on
scatter(t,y,'r', 'filled')
xlabel('t');
ylabel('Regression function');
title('Best fit function for data points')

figure;
plot(fvals, '-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Iteration');
ylabel('Objective Function Value');
title('Convergence of Objective Function');
grid on;

% Plot the trajectory in R^2
figure;
plot(xvals(1, :), xvals(2, :), '-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('x_1');
ylabel('x_2');
title('Optimization Path in R^2');
grid on;
axis equal;
%% DFG -- Davidon Fletcher Powell
format rat
A = [2 0
     0 8];
b = [4 8]';
maxiter = 10;
epsilon = 1e-4;

f = @(x) 1/2*x'*A*x-b'*x;

xk = [0 0]';
gk = A*xk - b;
Dk = [1 0
      0 1];

fvals = [0];
xvals = [xk];
for i = 1:maxiter
    alpha_k = gk'*Dk*gk / ( gk'*Dk'*A*Dk*gk );
    
    xnew = xk - alpha_k*Dk*gk;
    gnew = A*xnew - b;
    
    pk = xnew - xk;
    qk = gnew - gk;
    Dk = Dk + pk*pk' / (pk'*qk) - Dk*qk*qk'*Dk / (qk'*Dk*qk);
    
    xk = xnew;
    gk = gnew;
    if (i>1 && fvals(i-1)-f(xk) < epsilon)
        break;
    end
    fvals = [fvals f(xk)];
    xvals = [xvals xk];
end

% Plot the objective function values
figure;
plot(fvals, '-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Iteration');
ylabel('Objective Function Value');
title('Convergence of Objective Function');
grid on;

% Plot the trajectory in R^2
figure;
plot(xvals(1, :), xvals(2, :), '-o', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
scatter(xvals(1,end), xvals(2,end), 100, 'r', 'filled', 'DisplayName', 'Optimal Point');
xlabel('x_1');
ylabel('x_2');
title('Optimization Path in R^2');
legend('Path', 'Optimal Point');
grid on;
axis equal;
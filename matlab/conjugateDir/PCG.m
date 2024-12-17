%% Preconditioned Conjugate Direction Method
% Define Q, b och H
Q = [1001, 1, 1, 1, 1;
     1, 999, 1, 1, 1;
     1, 1, 101, 1, 1;
     1, 1, 1, 99, 1;
     1, 1, 1, 1, 10];
b = [1; 1; 1; 1; 1];

S = diag([1/1000, 1/1000, 1/100, 1/100, 1/10]);
H = S*S';

% Start point and initalization
x_k = zeros(5, 1); % Start in origin
g_k = Q * x_k - b; % Initial residual
d_k = -H * g_k;
max_iter = 100;
tolerance = 1e-10;

% To store the function values and gradient norms
f_values = [];
gradient_norms = [];

% Optimal sol (x*)
x_star = Q \ b;
f_star = 0.5 * x_star' * Q * x_star - b' * x_star;

% Iterate PCG
for k = 1:max_iter
    gradient_norm = norm(g_k);
    f_k = 0.5 * x_k' * Q * x_k - b' * x_k;

    % Save values
    f_values = [f_values; f_k - f_star];
    gradient_norms = [gradient_norms; gradient_norm];

    % Check convergence
    if gradient_norm < tolerance
        break;
    end

    % Calc steplength
    alpha_k = -(g_k' * d_k) / (d_k' * Q * d_k);

    % Update sol
    x_k = x_k + alpha_k * d_k;

    % Update residual
    g_new = Q*x_k - b;

    % Update beta
    beta_k = (g_new' * H * g_new) / (g_k' * H * g_k);

    % Update search direction
    d_k = -H * g_new + beta_k * d_k;

    % Prepare next iteration
    g_k = g_new;
end

% Ploa results
iterations = 1:length(f_values);

figure;

% Difference in function values
subplot(1, 2, 1);
semilogy(iterations, f_values, '-o');
xlabel('Iterations');
ylabel('f(x_k) - f(x^*)');
title('Convergence of function values (PCG)');
grid on;

% Gradient norm
subplot(1, 2, 2);
semilogy(iterations, gradient_norms, '-o');
xlabel('Iterations');
ylabel('||\nabla f(x_k)||');
title('Convergence of gradient norm (PCG)');
grid on;

hold off
disp('Eigenvalues for Q:');
disp(eig_Q);

disp('Eigenvalues for H*Q:');
disp(eig_HQ);
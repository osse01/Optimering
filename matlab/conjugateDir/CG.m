%% Conjugate gradient method
% Definine Q och b
Q = [1001, 1, 1, 1, 1;
     1, 999, 1, 1, 1;
     1, 1, 101, 1, 1;
     1, 1, 1, 99, 1;
     1, 1, 1, 1, 10];
b = [1; 1; 1; 1; 1];

% Start point och initialization
x_k = zeros(5, 1); % Start in origo
r_k = b - Q * x_k; % Initial residual
p_k = r_k;
max_iter = 100;
tolerance = 1e-10;

% To store the function values and gradient norms
f_values = [];
gradient_norms = [];

% Optimal sol (x*)
x_star = Q \ b;
f_star = 0.5 * x_star' * Q * x_star - b' * x_star;

% Iterate conjugate gradient method
for k = 1:max_iter
    gradient_norm = norm(r_k);
    f_k = 0.5 * x_k' * Q * x_k - b' * x_k;

    % Save values
    f_values = [f_values; f_k - f_star];
    gradient_norms = [gradient_norms; gradient_norm];

    % Check convergence
    if gradient_norm < tolerance
        break;
    end

    % Calc steplength
    alpha_k = (r_k' * r_k) / (p_k' * Q * p_k);

    % Update sol
    x_k = x_k + alpha_k * p_k;

    % Update residual
    r_new = r_k - alpha_k * Q * p_k;

    % Update beta
    beta_k = (r_new' * r_new) / (r_k' * r_k);

    % Update search dir
    p_k = r_new + beta_k * p_k;

    % Prepare for next iteration
    r_k = r_new;
end

% Plot results
iterations = 1:length(f_values);

figure;

% Diff in function values
subplot(1, 2, 1);
semilogy(iterations, f_values, '-o');
xlabel('Iterations');
ylabel('f(x_k) - f(x^*)');
title('Convergence of function values');
grid on;

% Gradient norm
subplot(1, 2, 2);
semilogy(iterations, gradient_norms, '-o');
xlabel('Iterations');
ylabel('||\nabla f(x_k)||');
title('Convergence of gradient norm');
grid on;


function t = armijo(x_curr, f, g, h)
% Variables
epsilon = 1e-9;         % Error tolerance.
sigma = 0.5;            % Must be in [0,1].
maxiter = 50;

for i = 1:maxiter
    % Initial Newton Step.
    t = (norm(g(x_curr))^2)/( g(x_curr)'*(h(x_curr)*g(x_curr) ));

    % Descent Direction.
    d = -g(x_curr);

    while (f(x_curr + t * d) >= f(x_curr) + sigma * g(x_curr)' * d * t)
        t = t * 0.5;  
    end
    gk = g(x_curr);
    % Update x.
    x_curr = x_curr + t * d;  
    if (norm(gk - g(x_curr)) <= epsilon)
        break;
    end
end
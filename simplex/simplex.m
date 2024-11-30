function x = simplex( f, x )
% Simplex alg. that return the new array x based on the function f.

% Max
[f_max,kMax] = max( f( x(1,1:end),x(2,1:end) ) );
x_max = x(1:end,kMax);
% Min
[f_min,kMin] = min( f( x(1,1:end),x(2,1:end) ) );
x_min = x(1:end,kMin);

% Middle point
[n,~] = size(x);
x_hat = 1/n *(sum( [x(1,1:end)', x(2,1:end)'] ) - x_max')';

% Mirror point
x_ref = 2*x_hat - x_max;


f_high = max( f(x(1,1:end),x(2,1:end)) ~= f_max ); % The next highest f.
f_ref = f(x_ref(1), x_ref(2)); % Reference value.

% The diffrent cases:
if ( f_ref < f_min ) % x_ref best!
    x_exp = x_ref + (x_ref - x_hat);
    if ( x_exp < f_ref )
        x_new = x_exp;
    else
        x_new = x_ref;
    end
elseif ( f_high > f_ref )
    x_new = x_ref;
elseif ( f_ref >= f_high )
    if ( f_ref >= f_max )
        x_new = 0.5 * (x_max+x_hat);
    else
        x_new = 0.5 * (x_ref + x_hat);
    end
end
x(1:end,kMax) = x_new;
end
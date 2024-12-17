%% min f(x) = 1/2 x'Qx - b'x using conjugate directions
% Starting with the unit vectors as start directions

Q = [2   2   0;
     2   4  -1;
     0 -1  2 ]
b = [-1 -3 1]';

xi = [[1 0 0]' [0 1 0]' [0 0 1]']; % Start with unit vectors
d =  [[1 0 0]' [0 0 0]' [0 0 0]'];       % First con. dir. as xi(1)

a_21 = -xi(:,2)'*Q*d(:,1) / ( d(:,1)'*Q*d(:,1) )
d(:,2) = xi(:,2) + a_21*d(:,1)

a_31 = -xi(:,3)'*Q*d(:,1) / ( d(:,1)'*Q*d(:,1) )
a_32 = -xi(:,3)'*Q*d(:,2) / ( d(:,2)'*Q*d(:,2) )
d(:,3) = xi(:,3) + a_31*d(:,1) + a_32*d(:,2)

% Test if they are conjugate
%d(:,1)'*Q*d(:,2)
%d(:,1)'*Q*d(:,3)
%d(:,2)'*Q*d(:,3)

% Find optimum
x_1 = [0 0 0]';
g_1 = -b;
alpha_1 = -g_1'*d(:,1)/( d(:,1)'*Q*d(:,1) )
x_2 = x_1 + alpha_1 * d(:,1)

g_2 = Q*x_2 - b;
alpha_2 = -g_2'*d(:,2)/( d(:,2)'*Q*d(:,2) )
x_3 = x_2 + alpha_2 * d(:,2)

g_3 = Q*x_3 - b;
alpha_3 = -g_3'*d(:,3)/( d(:,3)'*Q*d(:,3) )
x_4 = x_3 + alpha_3 * d(:,3)


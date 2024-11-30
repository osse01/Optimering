%% min f(x) = x'Qx - b'x using conjugate directions
% Starting with the unit vectors as start directions

Q = [1   1   1;
     1   2  -1/2;
     0 -1/2  1 ]
b = [-1 -3 1]';

xi = [[1 0 0]' [0 1 0]' [0 0 1]']; % Start with unit vectors
d = [[1 0 0]' [0 0 0]' [0 0 0]'];       % First con. dir. as xi(1)

a_21 = -xi(:,2)'*Q*d(:,1) / ( d(:,1)'*Q*d(:,1) )
d(:,2) = xi(:,2) + a_21*d(:,1)

a_31 = -xi(:,3)'*Q*d(:,1) / ( d(:,1)'*Q*d(:,1) )
a_32 = -xi(:,3)'*Q*d(:,2) / ( d(:,2)'*Q*d(:,2) )
d(:,3) = xi(:,3) + a_31*d(:,1) + a_32*d(:,2)

d(:,1)'*Q*d(:,2)
d(:,1)'*Q*d(:,3)
d(:,2)'*Q*d(:,3)
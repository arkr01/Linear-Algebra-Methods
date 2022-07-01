function [Q , R] = house(A)
    [m , n] = size(A);
    Q = eye(m , m);
    for k = 1 : n
        z = A (k : m, k);
        e1 = [1; zeros(m -k ,1)];
        u = z + sign(z(1)) * norm(z) * e1;
        u = u / norm(u);
        A(k : m, k : n) = A(k : m, k : n) - 2 * u * (u' * A(k : m, k : n));
        r = length(u);
        Q = blkdiag(eye(m - r), eye(r , r) - 2 * u * (u')) * Q;
    end
    Q = Q';
    R = A;
end
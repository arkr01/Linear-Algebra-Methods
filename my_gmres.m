function [x, rel_res, k, flag] = my_gmres(A, b, x0, tol, iter)
    flag = false;
    rel_res = zeros(iter, 1);
    Q = zeros(length(b), iter + 1);
    H = zeros(size(Q));
    
    r0 = b - A * x0;    
    nr = norm(r0);
    Q(:, 1) = r0 / nr;
    nb = norm(b);
    
    for k=1:iter
        [H(1 : k + 1, k), Q] = arnoldiOne(A, Q, k);
        [U, R] = householder(H(1 : k + 1, 1 : k));
        U = U(1 : k + 1, 1 : k);
        R = R(1 : k, 1 : k);
        e1 = [1; zeros(k, 1)];
        y = nr * (R \ ((U.') * e1));
        x = x0 + Q(:, 1 : k)*y;
        rkn = nr*sqrt(1-norm((U.') * e1)^2);
        rel_res(k) = rkn / nb;
        if rel_res(k) <= tol
            flag = true;
            break
        end
    end
    rel_res = rel_res(1 : k);
end
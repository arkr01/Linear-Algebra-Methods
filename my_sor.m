function [x, rel_res, num] = my_sor(A, b, x, max_iter)
    N = 15;
    omega = 2 / (1 + sin(pi / (N + 1)));
    E = tril(A, -1);
    D = diag(diag(A));
    M = (1 / omega) * D + E;
    epsilon = 10^(-8);
    num = 0;
    rel_res = zeros(1, max_iter);
    while true
        num = num + 1;
        r = b - A * x;
        x = x + fwd_sub(M, r);
        rel_res(num) = norm(r) / norm(b);
        if rel_res(num) <= epsilon || num == max_iter
            break;
        end
    end
end
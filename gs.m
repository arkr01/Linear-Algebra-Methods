function [x, rel_res, num] = gs(A, b, x, max_iter)
    dpe = tril(A);
    epsilon = 10^(-8);
    num = 0;
    rel_res = zeros(1, max_iter);
    while true
        num = num + 1;
        r = b - A * x;
        x = x + fwd_sub(dpe, r);
        rel_res(num) = norm(r) / norm(b);
        if rel_res(num) <= epsilon || num == max_iter
            break;
        end
    end
end
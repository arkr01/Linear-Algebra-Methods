function [x, rel_res, k] = restarted_gmres(A, b, x0, tol, m, iter)
    rel_res = zeros(iter, 1);
    [x, new_res, k, flag] = my_gmres(A, b, x0, tol, m);
    rel_res(1 : k) = new_res;
    while ~flag
        [x, new_res, l, flag] = my_gmres(A, b, x, tol, m);
        rel_res(k + 1 : k + l) = new_res;
        k = k + l;
    end
    rel_res = rel_res(1 : k);
end
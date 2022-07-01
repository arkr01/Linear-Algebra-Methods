function [eigvecs, Z, eigvals, explained] = my_pca(x)
    % INPUT x - data as matrix of dimension n_samples x n_features


    mu = mean(x);
%     x = x - mu;  % Unnecessary technically, included as a talking point
    c = cov(x);  % Centres data already --> E(X-EX)=EX-EX=0
    [eigvecs, eigvals] = eig(c);
    eigvals = diag(eigvals);

    % 0 eigenvalues explain none of the variance, these eigenpairs are
    % useless
    %
    % Remove 0 eigenvalues (and their corresponding eigenvectors), by
    % filtering out (basically) zero
    near_zero = 10 ^ (-10);
    indices = find(eigvals < near_zero);
    eigvals(indices) = [];
    eigvecs(:, indices) = [];
    
    % Sort eigenpairs in descending order
    [eigvals, indices] = sort(eigvals, 'descend');
    eigvecs = eigvecs(:, indices);

    % Calculate projections to new space (first k columns --> k dimension reduction)
    Z = ((eigvecs.') * x.').';

    % Calculate percentage of total variance explained by each principal
    % component (eigenvalue explains/measures variance explained by
    % corresponding eigenvector, i.e. the principal component)
    explained = 100 * (cumsum(eigvals) / sum(eigvals));
end
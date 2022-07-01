function [posteriors, likelihoods] = parametric_classifier(train, x)
%parametric_classifier 1-D, 2-class, real-valued Probabilistic Parametric
%   Classfier.
%   parametric_classifier(train, x) returns posterior and likelihood
%   probabilities of input x to classify, given training data (of the
%   form (data, labels))

    % Split training data into data and labels
    data = train(:, 1);
    labels = train(:, 2);
    n = length(labels);

    % Split data into classes + Initialise probabilities
    classed = NaN(n, 2);
    priors = (1 / 2) * ones(2, 1);
    posteriors = zeros(size(priors));
    likelihoods = zeros(size(posteriors));

    for i = 1 : n
        classed(i, labels(i)) = data(i);
    end
    
    % Calculate likelihoods
    for i = 1 : 2
        classed_actual = classed(:, i);
        classed_actual = classed_actual(~isnan(classed_actual));

        max_like = [mean(classed_actual), std(classed_actual, 1)];
%         max_like = mle(classed_actual);
        likelihoods(i) = normpdf(x, max_like(1), max_like(2));
    end
    
    % Calculate posteriors using Bayes Rule
    % Can be simplified to 'posteriors = likelihoods ./ sum(likelihoods)'
    % as equal priors factorise and cancel
    posteriors = (likelihoods .* priors) ./ (sum(likelihoods .* priors));
end
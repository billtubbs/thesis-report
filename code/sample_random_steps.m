function x = sample_random_steps(s, p, sigma)
% sample_random_steps(s, p, sigma)
%
% Generate an array of values for one or more random
% variables x(k,i) as defined in section 2.1 of MacGregor
% et al. (1984).
%
% Arguments:
%   s : array or scalar
%       the size of the output array [nT na] where nT is the
%       number of samples and na is the number of indepedent
%       sequences (na >= 1).  If s is a scalar then it is 
%       assumed that na = 1 and s = nT.
%   p : scalar (if na == 1) or vector (na > 1)
%       p1(i) is the probability of x(k,i) being a non-zero
%       value (0 < p1 < 1).
%   sigma : (optional) scalar (if na == 1) or vector (na > 1)
%       sigma(i) is the standard deviation of the non-zero
%       values of x(k,i). If sigma is not provided, sigma(i)
%       = 1 is used.
%
% Examples:
% >> sample_random_steps(5, 0.2)
% 
% ans =
% 
%          0
%          0
%    -1.3077
%          0
%          0
% 
% >> sample_random_steps([5 2], [0.8; 0.2], [10; 1])
% 
% ans =
% 
%          0    0.0749
%          0    1.1233
%     6.3702   -0.0331
%          0         0
%          0   -0.0977
% 
    if isscalar(s)
        s = [s 1];
    end
    if nargin < 3
        sigma = ones(s(2), 1);
    end
    if isscalar(sigma)
        sigma = sigma*ones(s(2), 1);
    end
    x = zeros(s);
    r = rand(s) < p';
    for i=1:size(x, 2)
        n = sum(r(:, i));
        x(r(:,i), i) = sigma(i) * randn(n, 1);
    end
end

function D1 = seqlhd_alg1(n1, d, m1)
% SEQLHD_ALG1 generates a sequential Latin hypercube design matrix
%   D1 = seqlhd_alg1(n1, m1, d) generates a sequential Latin hypercube design
%   matrix based on the input parameters.
%
%   Input arguments:
%       n1: Positive integer indicating the number of runs.
%       d: Positive integer indicating the number of dimensions.
%       m1: (Optional) Positive integer indicating the distance parameter. Default: ceil(n1/2)
%
%   Output:
%       D1: Matrix of size n1*d.
%
%   Example 1:
%       n1 = 10;
%       d = 3;
%       D1 = seqlhd_alg1(n1, d); % m1 will default to ceil(n1/2)
%   Example 2:
%       n1 = 10;
%       d = 3;
%       m1 = 4;
%       D1 = seqlhd_alg1(n1, d, m1); 
%
%   Notes:
%       - Input parameters n1 and m1 must be positive integers.
%       - Input parameter d represents the number of dimensions and must be a positive integer.
%       - The output matrix D1 has a size of n1*d.
%
%   Author: gongjunpeng21@mails.ucas.ac.cn
%   Date: 2023.07.16

% Check if m1 is provided, otherwise set default value
if nargin < 3
    m1 = ceil(n1/5)+2;
end

% Initialize matrix D1
D1 = zeros(n1, d);

% Generate sequential Latin hypercube design matrix
for oo = 1:d
    D1(1, oo) = (1 - rand) / n1;
    
    for i = 2:n1
        D1(i, oo) = D1(i-1, oo) + 1/n1 - rand*m1/((m1+n1)*n1);
        
        if D1(i, oo) <= (i-1)/n1
            D1(i, oo) = D1(i, oo) + 1/n1;
        end
    end
    
    t = randperm(n1);
    D1(:, oo) = D1(t, oo);
end

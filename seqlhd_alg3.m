function E3 = seqlhd_alg3(E2, n3)
% SEQLHD_ALG3 extends a Latin hypercube design matrix
%   E3 = seqlhd_alg3(E2, n3) extends the input Latin hypercube design matrix E2
%   by adding n3 additional samples to it.
%
%   Input arguments:
%       E2: N2-by-d matrix representing the original Latin hypercube design.
%       n3: Positive integer indicating the number of additional samples to be added.
%
%   Output:
%       E2: (N2+n3)-by-d matrix representing the extended Latin hypercube design.
%
%   Example:
%       E2 = [0.2 0.5; 0.6 0.8];
%       n3 = 3;
%       E3 = seqlhd_alg3(E2, n3);
%
%   Notes:
%       - Input matrix E2 must have dimensions N2-by-d, where N2 is the number of existing samples and d is the number of dimensions.
%       - Input parameter n3 must be a positive integer.
%       - The output matrix E2 has dimensions (N2+n3)-by-d.
%
%   Author: gongjunpeng21@mails.ucas.ac.cn
%   Date: 2023.07.16
%%
[N2, d] = size(E2);
N3 = n3 + N2;
D2_index = ceil(E2 * N3); 
E3 = zeros(N3, d);

% Extend the Latin hypercube design matrix
for i = 1:d
    % Copy the existing samples from E2 to E2 based on D2_index
    for j = 1:N2
        E3(D2_index(j, i), i) = E2(j, i);
    end
    
    % Fill in the remaining empty cells in E2
    for j = 1:N3
        if E3(j, i) == 0
            E3(j, i) = (j - rand) / N3;
        end
    end
    
    % Generate the additional samples for E2
    D2_i = setdiff(E3(:, i), E2(:, i));
    t = randperm(n3);
    D2_i = D2_i(t);
    E3(:, i) = [E2(:, i); D2_i];
end

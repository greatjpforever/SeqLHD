function E2 = seqlhd_alg2(E1, n2, m1)
% SEQLHD_ALG2 extends a Latin hypercube design matrix.
%
%   E2 = seqlhd_alg3(E1, n2, m1) extends the input Latin hypercube design matrix E1
%   by adding n2 additional samples to it.
%
%   This algorithm requires that E1 is generated from seqlhd_alg1.
%   It also requires the corresponding m1 parameter to be provided, 
%   unless the default m1 value was not modified when using seqlhd_alg1.
%
%   Input arguments:
%       E1: n1-by-d matrix representing the original Latin hypercube design.
%       n2: Positive integer indicating the number of additional samples to be added.
%       m1: Positive integer indicating the distance parameter of E1. Default: ceil(n1/5)+2
%
%   Output:
%       E2: (n1+n2)-by-d matrix representing the extended Latin hypercube design.
%
%   Example:
%       E1 = [0.2 0.5; 0.6 0.8];
%       n2 = 3;
%       m1 = 1;
%       E2 = seqlhd_alg2(E1, n2, m1);
%
%   Notes:
%       - Input matrix E1 must be generated from seqlhd_alg1.
%       - Input parameter n2 must be a positive integer.
%       - Input parameter m1 represents the distance parameter of E1. 
%       - If m1 is not provided, it will default to ceil(n1/2).
%       - The output matrix E2 has dimensions (n1+n2)-by-d.
%
%   Author: gongjunpeng21@mails.ucas.ac.cn
%   Date: 2023.07.16

[n1, d] = size(E1);
N2 = n1 + n2;

if nargin < 3
    m1 = ceil(n1/5)+2;
end

E2 = zeros(N2, d);
H = zeros(n2, 1);

% Error test
if n2 < m1
    error('The number of additional samples n2 must be greater than or equal to m1.');
end

for oo = 1:d
    D1_index = ceil(E1(:, oo) * N2); % J_oo of Algorithm 2
    for i = 1:n1
        E2(D1_index(i), oo) = E1(i, oo); % D_1 into E2
    end

    k = 1;
    for i = 1:N2
        if E2(i, oo) == 0 % D2
            coin = rand;
            if coin <= 0.5 
                if i == 1
                    E2(i, oo) = (1-rand) / (2 * N2);        % 0\leq rand <1
                elseif ismember(i - 1, D1_index)            % E2(i-1,oo) is in D1
                    E2(i, oo) = prob_of_x(n1, n2, m1, E2(i-1, oo), i) / (2 * N2) + (i - 1) / N2;
                else % E2(i-1,oo) is not in D1
                    if E2(i - 1, oo) <= (i - 1.5) / N2
                        E2(i, oo) = (1-rand) / (2 * N2) + (i - 1) / N2;
                    else
                        E2(i, oo) = E2(i - 1, oo) + 1 / (2 * N2);
                    end
                end
            else 
                if i == N2
                    E2(i, oo) = 1 - rand / (2 * N2);
                elseif E2(i + 1, oo) ~= 0
                    E2(i, oo) = prob_of_x(n1, n2, m1, E2(i + 1, oo), i) / (2 * N2) + (2 * i - 1) / (2 * N2);
                else
                    E2(i, oo) = i / N2 - rand / (2 * N2);
                end
            end

            H(k) = E2(i, oo);
            k = k + 1;
        end
    end

    % Permutate E2(:,oo) randomly
    tt = randperm(n2);
    E2_diff = H(tt);
    E2(:, oo) = [E1(:, oo); E2_diff];
end
end


%% 
function p = prob_of_x(n1,n2,nf_1,x,j)
%%
% \varphi_{j-1}^{(1)}(x) or \varphi_{j+1}^{(2)}(x)
% proposition 2 and 3
%%
% 2023.02.16
% Email：gongjunpeng21@mails.ucas.ac.cn
%%
% clc
% clear
% close all
% n1 = 3; 
% n2 = 4;   
% nf_1 = 1;  
% x = 0.3188;
% j = 4;
% A_2:ans = 0.1884 
%%
N2 =n1+n2;

if x<=(j-1)/N2 % \varphi_{j-1}^{(1)}(x)  j=\barj+1
    i = ceil((j-2)*n1/N2);
    if j/N2>=(i+1)/n1 %% A_3 of Table 1: 两个1/N2的区间横跨三个1/n1区间
        p = phi_3interval(i,j,n1,n2,nf_1,x)/phi_3interval(i,j,n1,n2,nf_1,(j-1)/N2);
    else  %% 两个1/N2的区间横跨一个或两个1/n1区间
        if i/n1 == (j-2)/N2  || i/n1>=j/N2 %% A_1 of Table 1: 两个1/N2的区间在同一1/n1区间中
            p = (x-(j-2)/N2)*N2;
        elseif (j-2)/N2<i/n1 && i/n1<(j-1)/N2 %% A_2 of Table 1: 两个1/N2的区间横跨两个1/n1区间且左侧区间被分割
            p = phi_leftsep(i,j,n1,n2,nf_1,x)/phi_leftsep(i,j,n1,n2,nf_1,(j-1)/N2);
        else     %% A_4 of Table 1: 两个1/N2的区间横跨两个1/n1区间且右侧区间被分割
            p = phi_rightsep(i,j,n1,n2,nf_1,x)/phi_rightsep(i,j,n1,n2,nf_1,(j-1)/N2);
        end
    end
else     % \varphi_{j+1}^{(2)}(x)  j=\barj-1
    i = ceil((j-1)*n1/N2);
    if (j+1)/N2 >= (i+1)/n1    %% 两个1/N2的区间横跨三个1/n1区间
        p = 1-phi_3interval(ceil((N2-1-j)*n1/N2),N2+1-j,n1,n2,nf_1,1-x)/phi_3interval(ceil((N2-1-j)*n1/N2),N2+1-j,n1,n2,nf_1,(N2-j)/N2);
    else     %% 两个1/N2的区间横跨一个或两个1/n1区间
        if i/n1 == (j-1)/N2 || i/n1 >= (j+1)/N2   %% 两个1/N2的区间在同一1/n1区间中
            p = (x-j/N2)*N2;    
        elseif  i/n1 > (j-1)/N2 && i/n1 <= j/N2  %% 两个1/N2的区间横跨两个1/n1区间且左侧区间被分割
            p = 1-phi_rightsep(ceil((N2-1-j)*n1/N2),N2+1-j,n1,n2,nf_1,1-x)/phi_rightsep(ceil((N2-1-j)*n1/N2),N2+1-j,n1,n2,nf_1,(N2-j)/N2);
        else  %% 两个1/N2的区间横跨两个1/n1区间且右侧区间被分割
            p = 1-phi_leftsep(ceil((N2-1-j)*n1/N2),N2+1-j,n1,n2,nf_1,1-x)/phi_leftsep(ceil((N2-1-j)*n1/N2),N2+1-j,n1,n2,nf_1,(N2-j)/N2);
        end
    end
end
end

%%
function y = phi_rightsep(i,j,n1,n2,nf,x)
N2 = n1+n2;
h1 = nf/(n1*(n1+nf));
low = max(0,h1-(j/N2-i/n1));
knot_original = [(i-1)/n1,(i-1)/n1+h1-low,j/N2-1/(n1+nf)-(h1-low),j/N2-1/(n1+nf),(j-1)/N2];
bottom_original = [h1,low,low,h1,h1];
for oo = 2:5
    if (j-2)/N2 >= knot_original(oo-1) && (j-2)/N2<knot_original(oo)
        knot = [(j-2)/N2,knot_original(oo:5)];
        b = bottom_original(oo-1)+((j-2)/N2-knot_original(oo-1))*(bottom_original(oo)-bottom_original(oo-1))/(knot_original(oo)-knot_original(oo-1));
        bottom = [b,bottom_original(oo:5)];
    end
end
t = length(knot);
y = 0;
for oo = 1:t-1
    if x > knot(oo+1)
        y = y+(knot(oo+1)-knot(oo))/2*(bottom(oo+1)+bottom(oo));
    elseif knot(oo+1)-knot(oo) ~= 0
        y = y+(x-knot(oo))/2*(bottom(oo)+(bottom(oo+1)-bottom(oo))/(knot(oo+1)-knot(oo))*(x-knot(oo))+bottom(oo));
        break;
    end
end
end
%% 
function y = phi_leftsep(i,j,n1,n2,nf_1,x)
N2 = n1+n2;
flag = 1/n1-2/N2;
h1 = nf_1/(n1*(n1+nf_1));
knot = zeros(1,5);
knot(4:5) = [i/n1,(j-1)/N2];
bottom = zeros(1,5);
bottom(4:5) = [h1,h1];
a = max(0,h1+(j-2)/N2-i/n1);
b = min(0,h1+(j-2)/N2-i/n1);
if flag>=h1   % 2/N_2<= 1/(n_1+n_{f_1})
    knot(1:3) = [(j-2)/N2,(j-2)/N2,(j-2)/N2];
    bottom(1:3) = [h1,h1,h1];
elseif flag >= a   % 2/N_2<= 1/n_1 & j/N_2<=(i+1)/n_1-h_1
    knot(1:3) = [(j-2)/N2,(j-2)/N2,(j-2)/N2+h1-flag];
    bottom(1:3) = [flag,flag,h1];
elseif flag >= b    
    if b == 0
        knot(1:3) = [(j-2)/N2,(j-2)/N2,i/n1];
        bottom(1:3) = [flag,flag,flag+i/n1-(j-2)/N2];
    else
        knot(1:3) = [(j-2)/N2,(j-2)/N2-flag,(j-2)/N2-flag+h1];
        bottom(1:3) = [0,0,h1];        
    end
else
    knot(1:3) = [(j-2)/N2,(j-2)/N2-flag,i/n1];
    bottom(1:3) = [0,0,i/n1-(j-2)/N2+flag];            
end
y = 0;
for oo = 1:4
    if x > knot(oo+1)
        y = y+(knot(oo+1)-knot(oo))/2*(bottom(oo+1)+bottom(oo));
    elseif knot(oo+1)-knot(oo)~=0
        y = y+(x-knot(oo))/2*(bottom(oo)+(bottom(oo+1)-bottom(oo))/(knot(oo+1)-knot(oo))*(x-knot(oo))+bottom(oo));
        break;
    end
end
end
%%
function y = phi_3interval(i,j,n1,n2,nf,x)
N2 = n1+n2;
h1 = nf/(n1*(n1+nf));
low = max(0,h1-(j/N2-(i+1)/n1));
knot = [i/n1, i/n1+h1-low, j/N2-1/(n1+nf)-(h1-low), j/N2-1/(n1+nf), (j-1)/N2];
bottom = [h1,low,low,h1,h1];
y = 0;
for oo = 1:4
    if x > knot(oo+1)
        y = y+(knot(oo+1)-knot(oo))/2*(bottom(oo+1)+bottom(oo));
    elseif knot(oo+1)-knot(oo) ~= 0
        y = y+(x-knot(oo))/2*(bottom(oo)+(bottom(oo+1)-bottom(oo))/(knot(oo+1)-knot(oo))*(x-knot(oo))+bottom(oo));
        break;
    end
end
end




    
    















clc
clear all
close all

%% 定义测试函数和参数
test_functions = {@func_log5X, @func_log_sqrt, @func_xsinx, @func_borehole, @func_Matyas, @func_Levy};
function_names = {'log5x', 'log_sqrt', 'xsinx', 'borehole', 'Matyas', 'Levy'};
d_values = [5, 2, 2, 8, 2, 2];
mu_values = [-50, -38.6304931998928, 1438.343187376281, 77.988249653226120, 17.333333333344665, 26.312023687344478];

%% Test 1
% Parameter
N1 = 5;
N2 = 11;
N3 = 25;
iter = 5*10^5;
folder_path = sprintf('data_compare_method_N1=%d_N2=%d_N3=%d_iter=%d', N1, N2, N3, iter);
if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end

% 运行测试函数并保存结果
for i = 1:length(test_functions)
    d = d_values(i);
    mu = mu_values(i);
    function_name = function_names{i};
    
    [Rmse, Var, Bias] = compare_method(test_functions{i}, mu, N1, N2, N3, d, iter);
    save(fullfile(folder_path, [function_name, '.mat']), 'Rmse', 'Var', 'Bias');
end
% RMSE
generate_rmse_table(folder_path, function_names);


%% Test 2
% Parameter
N1 = 20;
N2 = 50;
N3 = 100;
iter = 10^5;

folder_path = sprintf('data_compare_method_N1=%d_N2=%d_N3=%d_iter=%d', N1, N2, N3, iter);
if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end
% 运行测试函数并保存结果
for i = 1:length(test_functions)
    d = d_values(i);
    mu = mu_values(i);
    function_name = function_names{i};
    
    [Rmse, Var, Bias] = compare_method(test_functions{i}, mu, N1, N2, N3, d, iter);
    save(fullfile(folder_path, [function_name, '.mat']), 'Rmse', 'Var', 'Bias');
end
% RMSE
generate_rmse_table(folder_path, function_names);




%% Test 3
% Parameter
N1 = 50;
N2 = 120;
N3 = 250;
iter = 10^3;

folder_path = sprintf('data_compare_method_N1=%d_N2=%d_N3=%d_iter=%d', N1, N2, N3, iter);
if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end
% 运行测试函数并保存结果
for i = 1:length(test_functions)
    d = d_values(i);
    mu = mu_values(i);
    function_name = function_names{i};
    
    [Rmse, Var, Bias] = compare_method(test_functions{i}, mu, N1, N2, N3, d, iter);
    save(fullfile(folder_path, [function_name, '.mat']), 'Rmse', 'Var', 'Bias');
end
% RMSE
generate_rmse_table(folder_path, function_names);






%% def function
function generate_rmse_table(folder_path, function_names)
    % 初始化表格字符串
    table_str_rmse = '';
    table_str_rmse = [table_str_rmse, '\\begin{tabular}{cccccccc}\n'];
    table_str_rmse = [table_str_rmse, ' Function & Stage & IID & DLHD & ILHD & SD & FNLHD & Ours '];
    table_str_rmse = [table_str_rmse, '\\\\ \n'];

    for i = 1:length(function_names)
        % 读取mat文件中的Rmse数据
        file_path = fullfile(folder_path, [function_names{i}, '.mat']);
        load(file_path, 'Rmse');
        Rmse = [Rmse(:, 1), Rmse(:, 3), Rmse(:, 2), Rmse(:, 4:6)];

        % 生成latex表格字符串
        Fi = ['$f_', num2str(i), '(x)$']; % 将数字 i 转换为字符并连接
        table_str_rmse = [table_str_rmse, Fi, ' & 1 & '];

        % 第一行
        for k = 1:size(Rmse, 2)
            table_str_rmse = [table_str_rmse, sprintf('%.3f', Rmse(1, k)), ' & '];
        end
        table_str_rmse = [table_str_rmse, '\\\\ \n & 2 & '];
        % 第二行
        for k = 1:size(Rmse, 2)
            table_str_rmse = [table_str_rmse, sprintf('%.3f', Rmse(2, k)), ' & '];
        end
        % 第三行
        if size(Rmse, 1) == 3
            table_str_rmse = [table_str_rmse, '\\\\ \n & 3 & '];
            [min_values_3, min_indices_3] = mink(nonzeros(Rmse(3, :)), 1);
            idx_3 = find(Rmse(3, :) <= min_values_3 * 1);
            for k = 1:size(Rmse, 2)
                table_str_rmse = [table_str_rmse, sprintf('%.3f', Rmse(3, k)), ' & '];
            end
        end
        table_str_rmse = [table_str_rmse, '\\\\\n'];
    end
    table_str_rmse = [table_str_rmse, '\\end{tabular}\n'];
    table_str_rmse = strrep(table_str_rmse, '& \\\\', '\\\\');

    % 保存表格到文件
    file_path = fullfile(folder_path, 'Rmse_table.txt');
    fid = fopen(file_path, 'w');
    fprintf(fid, table_str_rmse);
    fclose(fid);
end

%% f1
function y = func_log5X(x)
y = 10*log(x(:,1).*x(:,2).*x(:,3).*x(:,4).*x(:,5));
end
%% f2
function f=func_log_sqrt(x)
f=100*log(1./sqrt(x(:,1)+1./sqrt(x(:,2))));
end
%% f3
function y = func_xsinx(x)
y = 10*(30+x(:,1).*sin(x(:,1))).*(4+exp(-x(:,2).^2));
end
%% f4
function f=func_borehole(x)
a=[0.05 100 63070 990 63.1 700 1120 9855];
b=[0.15 50000 115600 1110 116 820 1680 12045];
[n,~] = size(x);
for i = 1:n
    x(i,:) = x(i,:).*(b-a)+a;
end
f=(2*pi*x(:,3).*(x(:,4)-x(:,6)))./(log10(x(:,2)./x(:,1)).*(1+x(:,3)./x(:,5)+2*x(:,3).*x(:,7)./(x(:,1).^2.*x(:,8).*log10(x(:,2)./x(:,1)))));
end
%% f5
function y = func_Matyas(x)
x = (x-0.5)*20;
y = 0.26*(x(:,1).^2+x(:,2).^2)-0.48*x(:,1).*x(:,2);
end
%% f6
function y = func_Levy(x)
x = 1+((x-0.5)*20-1)/4;
y =sin(pi*x(:,1)).^2+(x(:,1)-1).^2.*(1+10*sin(pi*x(:,2)).^2)+(x(:,2)-1).^2.*(1+10*sin(2*pi*x(:,2)).^2);
end


%% Function of designs_construct
%% SD
function D = design_SD(n1,n2,n3,m)
%% 生成第一层实验
N1 = n1;
D1 = zeros(N1,m);
for i = 1:m
    D1(1:N1,i) = (randperm(N1)'-rand(N1,1))/N1;
end
%% 生成第二层实验
D2 = zeros(n2-n1,m);
N2 = N1+ lcm(N1,n2-n1);  % n_1+L(n_1,n_A)
l2 = 1/N2;
z = [0:N2-1]/N2;
for j = 1:m
    dj = floor(D1(:,j)*N2)/N2;
    Bj = setdiff(z,dj);
    Bj = sort(Bj); %% lcm(N1,n2-n1)
    D2(:,j) = (randperm(n2-n1)'-rand(n2-n1,1))/(n2-n1);
    for i=1:n2-n1
        k = floor((N2-n1)*D2(i,j))+1;
        r = Bj(k);
        D2(i,j)= l2*rand+r;
    end
    D2(:,j) = D2(randperm(n2-n1),j);
end
if n3==n2
    D=[D1;D2];
else
    %% 生成第三层实验
    if N2 - n2 == 0
        N3 = N2 + lcm(N2, n3-n2);
    else
        N3 = N2 * lcm(N2-n2, n3-n2);
    end
    l3 = 1/N3;
    D3 = zeros(n3-n2,m);
    D_cup = [D1;D2];
    %%
    for j = 1:m
        z = [0:N3-1]/N3;
        dj = floor(D_cup(:,j)*N3)/N3;
        Bj = setdiff(z,dj);
        Bj = sort(Bj);
        D3(:,j) = (randperm(n3-n2)'-rand(n3-n2,1))/(n3-n2);
        for i=1:n3-n2
            k = floor((N3-n2)*D3(i,j))+1;
            r = Bj(k);
            D3(i,j)= l3*rand+r;
        end
        D3(:,j) = D3(randperm(n3-n2),j);
    end
    % 拼接三层实验
    D = [D1;D2;D3];
end
end

%% FNLHD
function D = design_FNLHD(n1,n2,n3,m)
s = lcm(lcm(n1,n2),n3);
t1 = s/n1;
t2 = s/n2;
t3 = s/n3;

D1 = zeros(n1,m);
D2 = zeros(n2-n1,m);
D3 = zeros(n3-n2,m);
D = zeros(n3,m);

for j=1:m
    perm = randperm(n1);
    oo = -1;
    while oo < 0
        % step 2
        for i = 1:n1
            % first layer, draw permutation from Z_n1;
            D1(i,j) = randsample((perm(i)-1)*t1+1:perm(i)*t1,1);
        end
        % step 3
        for l=1:n2
            Bl = ((l-1)*t2+1):l*t2;
            num_set = numel(setdiff(Bl,D1(:,j)));
            if num_set< t2-1
                oo = -1;
                break
            else
                oo = 1;
            end
        end
    end
    % repeat step 1-3 to generate the second layer using the remaining points
    % step 1
    set_2 = setdiff(1:n2,ceil(D1(:,j)/t2));
    perm = randperm(numel(set_2));
    set_2perm = set_2(perm);
    oo = -1;
    while oo < 0
        % step 2
        for i = 1:n2-n1 % step 5
            D2(i,j) = randsample((set_2perm(i)-1)*t2+1:set_2perm(i)*t2,1);
        end
        C = [D1(:,j);D2(:,j)];
        % step 3
        for l=1:n3 %  check for C
            Bl = ((l-1)*t3+1):l*t3;
            num_set = numel(setdiff(Bl,C));
            if num_set< t3-1
                oo = -1;
                break
            else
                oo = 1;
            end
        end
    end
    % generate the third layer using the remaining points
    % step 1
    set_3 = setdiff(1:n3,ceil(C/t3));
    perm = randperm(numel(set_3));
    set_3perm = set_3(perm);
    % step 2
    for i = 1:n3-n2
        if t3 == 1
            D3(i,j) = set_3perm(i);
        else
            D3(i,j) = randsample((set_3perm(i)-1)*t3+1:set_3perm(i)*t3,1);
        end
    end
    D(:,j) = ([D1(:,j);D2(:,j);D3(:,j)]-rand(n3,1))/s;
end
end

%% Seq LHD
function E3 = N3_generation(n1,n2,n3,m1,m2,d)
N2 = n2+n1;
N3 = N2+n3;
%%
E2 = N2_controlled_begin(n1,n2,m1,m2,d);
D2_index = ceil(E2*N3);
E3 = zeros(N3,d);
for i = 1:d
    for j = 1:N2
        E3(D2_index(j,i),i) = E2(j,i);
    end
    for j = 1:N3
        if E3(j,i) == 0
            E3(j,i) = (j-rand)/N3;
        end
    end
    D3_i = setdiff(E3(:,i),E2(:,i));
    t = randperm(N3-N2);
    D3_i = D3_i(t);
    E3(:,i) = [E2(:,i);D3_i];
end
end

function E2 = N2_controlled_begin(n1,n2,m1,m2,d)
N2 = n1+n2;
D1 = N1_controlled_begin(n1,m1,d);
E2 = zeros(N2,d); % E2
h2 = m2/((2*N2)*(2*N2+m2));
H = zeros(n2,1);
for oo = 1:d
    D1_index = ceil(D1(:,oo)*N2); %% J_oo of Algorithm 2
    for i = 1:n1
        E2(D1_index(i),oo) = D1(i,oo); %% D_1 into E2
    end
    k = 1;
    for i = 1:N2
        if E2(i,oo) == 0  % D2
            coin = rand;
            if coin<=0.5       % coin<= 0.5
                if i == 1
                    E2(i,oo) = rand/(2*N2);
                elseif ismember(i-1,D1_index) %% E2(i-1,oo) is in D1
                    E2(i,oo) = prob_of_x(n1,n2,m1,E2(i-1,oo),i)/(2*N2)+(i-1)/N2-rand*h2;
                else  %% E2(i-1,oo) is not in D1
                    if E2(i-1,oo)<(i-1.5)/N2
                        E2(i,oo) = rand/(2*N2)+(i-1)/N2;
                    else
                        E2(i,oo) = E2(i-1,oo)+1/(2*N2)-rand*h2;
                    end
                end
                % make E2(i,oo) fall into ((i-1)/N2,i/N2]
                if E2(i,oo)<=(i-1)/N2
                    E2(i,oo) = E2(i,oo)+1/(2*N2);
                end
            else    % coin>0.5
                if i == N2
                    E2(i,oo) = 1-rand/(2*N2);
                elseif E2(i+1,oo) ~= 0
                    E2(i,oo) = prob_of_x(n1,n2,m1,E2(i+1,oo),i)/(2*N2)+(2*i-1)/(2*N2)+rand*h2;
                else
                    E2(i,oo) = i/N2-rand/(2*N2);
                end
                % make E2(i,oo) fall into ((i-1)/N2,i/N2]
                if E2(i,oo)>i/N2
                    E2(i,oo) = E2(i,oo)-1/(2*N2);
                end
            end
            H(k) = E2(i,oo);
            k = k+1;
        end
    end
    %% permutate D1(:,oo) randomly
    t = randperm(n1);
    D1(:,oo) = D1(t,oo);
    %% permutate E2(:,oo) randomly
    tt = randperm(n2);
    E2_diff = H(tt);
    E2(:,oo) = [D1(:,oo);E2_diff];
end
end

function E1 = N1_controlled_begin(n1,m1,d)
E1 = zeros(n1,d);
for oo = 1:d
    E1(1,oo) = (1-rand)/n1;
    for i = 2:n1
        E1(i,oo) = E1(i-1,oo)+1/n1-rand*m1/((m1+n1)*n1);
        if E1(i,oo)<=(i-1)/n1
            E1(i,oo) = E1(i,oo)+1/n1;
        end
    end
    t = randperm(n1);
    E1(:,oo) = E1(t,oo);
end
end

function p = prob_of_x(n1,n2,m,x,j)
N2 =n1+n2;
if x<=(j-1)/N2 % \varphi_{j-1}^{(1)}(x);
    i = ceil((j-2)*n1/N2);
    if j/N2>=(i+1)/n1       %% 两个1/N2的区间横跨三个1/n1区间
        p = phi_3interval(i,j,n1,n2,m,x)/phi_3interval(i,j,n1,n2,m,(j-1)/N2);
    else  %% 两个1/N2的区间横跨一个或两个1/n1区间
        if i/n1 == (j-2)/N2  || i/n1>=j/N2 %% 两个1/N2的区间在同一1/n1区间中
            p = (x-(j-2)/N2)*N2;
        elseif (j-2)/N2<i/n1 && i/n1<(j-1)/N2 %% 两个1/N2的区间横跨两个1/n1区间且左侧区间被分割
            p = phi_leftsep(i,j,n1,n2,m,x)/phi_leftsep(i,j,n1,n2,m,(j-1)/N2);
        else     %% 两个1/N2的区间横跨两个1/n1区间且右侧区间被分割
            p = phi_rightsep(i,j,n1,n2,m,x)/phi_rightsep(i,j,n1,n2,m,(j-1)/N2);
        end
    end
else     % \varphi_{j+1}^{(2)}(x);
    i = ceil((j-1)*n1/N2);
    if (j+1)/N2 >= (i+1)/n1    %% 两个1/N2的区间横跨三个1/n1区间
        p = 1-phi_3interval(ceil((N2-1-j)*n1/N2),N2+1-j,n1,n2,m,1-x)/phi_3interval(ceil((N2-1-j)*n1/N2),N2+1-j,n1,n2,m,(N2-j)/N2);
    else     %% 两个1/N2的区间横跨一个或两个1/n1区间
        if i/n1 == (j-1)/N2 || i/n1 >= (j+1)/N2   %% 两个1/N2的区间在同一1/n1区间中
            p = (x-j/N2)*N2;
        elseif  i/n1 > (j-1)/N2 && i/n1 <= j/N2  %% 两个1/N2的区间横跨两个1/n1区间且左侧区间被分割
            p = 1-phi_rightsep(ceil((N2-1-j)*n1/N2),N2+1-j,n1,n2,m,1-x)/phi_rightsep(ceil((N2-1-j)*n1/N2),N2+1-j,n1,n2,m,(N2-j)/N2);
        else  %% 两个1/N2的区间横跨两个1/n1区间且右侧区间被分割
            p = 1-phi_leftsep(ceil((N2-1-j)*n1/N2),N2+1-j,n1,n2,m,1-x)/phi_leftsep(ceil((N2-1-j)*n1/N2),N2+1-j,n1,n2,m,(N2-j)/N2);
        end
    end
end
end

function y = phi_rightsep(i,j,n1,n2,m,x)
N2 = n1+n2;
h1 = m/(n1*(n1+m));
low = max(0,h1-(j/N2-i/n1));
knot_original = [(i-1)/n1,(i-1)/n1+h1-low,j/N2-1/(n1+m)-(h1-low),j/N2-1/(n1+m),(j-1)/N2];
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

function y = phi_leftsep(i,j,n1,n2,m,x)
N2 = n1+n2;
flag = 1/n1-2/N2;
h1 = m/(n1*(n1+m));
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

function y = phi_3interval(i,j,n1,n2,m,x)
N2 = n1+n2;
h1 = m/(n1*(n1+m));
low = max(0,h1-(j/N2-(i+1)/n1));
knot = [i/n1, i/n1+h1-low, j/N2-1/(n1+m)-(h1-low), j/N2-1/(n1+m), (j-1)/N2];
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

%% compare_method
function [Rmse,Var,Bias] = compare_method(func,mu,N1,N2,N3,d,iter)
%% IID
mean_iid = zeros(3,iter);
for oo = 1:iter
    D = zeros(N3,d);
    D(1:N1,:) = rand(N1,d);
    D(N1+1:N2,:) = rand(N2-N1,d);
    D(N2+1:N3,:) = rand(N3-N2,d);
    y = func(D);
    mean1 = mean(y(1:N1));
    mean2 = mean(y(1:N2));
    mean3 = mean(y);
    mean_iid(:,oo) = [mean1;mean2;mean3];
end
RMSE_iid = zeros(3,1);
sd_iid = zeros(3,1);
bias_iid = zeros(3,1);
for i = 1:3
    RMSE_iid(i) = sqrt((mean_iid(i,:)-mu)*(mean_iid(i,:)-mu)'/iter);
    sd_iid(i) = sqrt(var(mean_iid(i,:)));
    bias_iid(i) = mean(mean_iid(i,:)-mu);
end

%% AddLHD i阶段为独立的n_i次LHD
mean_Add = zeros(3,iter);
for oo = 1:iter
    D = zeros(N3,d);
    for i = 1:d
        D(1:N1,i) = (randperm(N1)'-rand(N1,1))/N1;
        D(N1+1:N2,i) = (randperm(N2-N1)'-rand(N2-N1,1))/(N2-N1);
        D(N2+1:N3,i) = (randperm(N3-N2)'-rand(N3-N2,1))/(N3-N2);
    end
    y = func(D);
    mean1 = mean(y(1:N1));
    mean2 = mean(y(1:N2));
    mean3 = mean(y);
    mean_Add(:,oo) = [mean1;mean2;mean3];
end
RMSE_Add = zeros(3,1);
sd_Add = zeros(3,1);
bias_Add = zeros(3,1);
for i = 1:3
    RMSE_Add(i) = sqrt((mean_Add(i,:)-mu)*(mean_Add(i,:)-mu)'/iter);
    sd_Add(i) = sqrt(var(mean_Add(i,:)));
    bias_Add(i) = mean(mean_Add(i,:)-mu);
end
%% diviLHDLHD
mean_diviLHD = zeros(3,iter);
for oo = 1:iter
    D = zeros(N3,d);
    for i = 1:d
        D(1:N3,i) = (randperm(N3)'-rand(N3,1))/N3;
    end
    y = func(D);
    mean1 = mean(y(1:N1));
    mean2 = mean(y(1:N2));
    mean3 = mean(y);
    mean_diviLHD(:,oo) = [mean1;mean2;mean3];
end
RMSE_diviLHD = zeros(3,1);
sd_diviLHD = zeros(3,1);
bias_diviLHD = zeros(3,1);
for i = 1:3
    RMSE_diviLHD(i) = sqrt((mean_diviLHD(i,:)-mu)*(mean_diviLHD(i,:)-mu)'/iter);
    sd_diviLHD(i) = sqrt(var(mean_diviLHD(i,:)));
    bias_diviLHD(i) = mean(mean_diviLHD(i,:)-mu);
end

%% SD
mean_sd = zeros(3,iter);
for oo = 1:iter
    D = design_SD(N1,N2,N3,d); % SD by kong
    y = func(D);
    mean1 = mean(y(1:N1));
    mean2 = mean(y(1:N2));
    mean3 = mean(y);
    mean_sd(:,oo) = [mean1;mean2;mean3];
end
RMSE_sd = zeros(3,1);
sd_sd = zeros(3,1);
bias_sd = zeros(3,1);
for i = 1:3
    RMSE_sd(i) = sqrt((mean_sd(i,:)-mu)*(mean_sd(i,:)-mu)'/iter);
    sd_sd(i) = sqrt(var(mean_sd(i,:)));
    bias_sd(i) = mean(mean_sd(i,:)-mu);
end
%% NLH_chen
mean_nl = zeros(3,iter);
for oo = 1:iter
    D = design_FNLHD(N1,N2,N3,d);  % chen的方法：flexible NLHD
    y = func(D);
    mean1 = mean(y(1:N1));
    mean2 = mean(y(1:N2));
    mean3 = mean(y);
    mean_nl(:,oo) = [mean1;mean2;mean3];
end
RMSE_nl = zeros(3,1);
sd_nl = zeros(3,1);
bias_nl = zeros(3,1);
for i = 1:3
    RMSE_nl(i) = sqrt((mean_nl(i,:)-mu)*(mean_nl(i,:)-mu)'/iter);
    sd_nl(i) = sqrt(var(mean_nl(i,:)));
    bias_nl(i) = mean(mean_nl(i,:)-mu);
end
%% Sequential LHD
n1= N1;
n2 = N2-N1;
n3 = N3-N2;
m1 = ceil(n1/5)+2;
m2 = 0;
mean_seq = zeros(3,iter);
for oo = 1:iter
    D = N3_generation(n1,n2,n3,m1,m2,d);
    y = func(D);
    mean1 = mean(y(1:N1));
    mean2 = mean(y(1:N2));
    mean3 = mean(y);
    mean_seq(:,oo) = [mean1;mean2;mean3];
end

RMSE_seq = zeros(3,1);
sd_seq = zeros(3,1);
bias_seq = zeros(3,1);
for i = 1:3
    RMSE_seq(i) = sqrt((mean_seq(i,:)-mu)*(mean_seq(i,:)-mu)'/iter);
    sd_seq(i) = sqrt(var(mean_seq(i,:)));
    bias_seq(i) = mean(mean_seq(i,:)-mu);
end

Rmse = [RMSE_iid,RMSE_Add,RMSE_diviLHD,RMSE_sd,RMSE_nl,RMSE_seq];
Var = [sd_iid,sd_Add,sd_diviLHD,sd_sd,sd_nl,sd_seq];
Bias = [bias_iid, bias_Add,bias_diviLHD,bias_sd,bias_nl,bias_seq];
end



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










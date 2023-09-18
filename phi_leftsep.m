function y = phi_leftsep(i,j,n1,n2,nf_1,x)
%%
% clc
% clear all
% n1 = 3; 
% n2 = 4;   
% nf_1 = 1;  
% x = 3/7;
% j = 4;
% i=1;
% ans = 0.0113
%%
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






    
    
function y = phi_3interval(i,j,n1,n2,nf,x)
%%
% 2018.07.23
% Email: xujin_nudt@163.com
%%
% 两个1/N2区间横跨三个1/n1区间
%%
%%
% n1 = 4;
% n2 = 2;
% nf = 1;
% x = 0.3257;
% j = 3;
% i = 1;
%%
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
